// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.archaeopteryx;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GraphicsEnvironment;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRootPane;
import javax.swing.KeyStroke;
import javax.swing.UIManager;

import org.forester.archaeopteryx.NodeDataDraft.Problem;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The per-node window: a {@link NodeDataForm} plus the chrome around it. {@link Mode#VIEW} ("Show node data") is
 * read-only with a single Close button; {@link Mode#EDIT} ("Edit node data") adds a status line (the first
 * validation problem, or "Unsaved changes") and a <b>Write to Tree</b> button, which is the default button and the
 * only thing that ever touches the tree. Closing with unsaved edits asks: Write and Close / Discard / Cancel. The
 * window is non-modal so several can stay open; the tree panel tracks them and closes them all when an undo or
 * redo swaps the tree underneath (their node would belong to the replaced tree).
 */
final class NodeFrame extends JFrame {

    private static final long  serialVersionUID = -6943510233968557246L;
    private final TreePanel    _tree_panel;
    /** The slot this window was opened into -- only used to cascade the window position. */
    private final int          _index;
    private final NodeDataForm _form;
    private final JLabel       _status;
    private final JButton      _write;
    private final String       _base_title;

    NodeFrame( final PhylogenyNode n, final TreePanel tp, final int index, final NodeDataForm.Mode mode ) {
        _tree_panel = tp;
        _index = index;
        _form = new NodeDataForm( n, tp, mode );
        _base_title = ( _form.isEditable() ? "Edit Node: " : "Node: " ) + NodeDataForm.nodeLabel( n );
        setTitle( _base_title );
        setDefaultCloseOperation( DO_NOTHING_ON_CLOSE );
        getContentPane().setLayout( new BorderLayout() );
        getContentPane().add( _form, BorderLayout.CENTER );
        // -- footer: status (left) + buttons (right) --
        final JPanel footer = new JPanel( new BorderLayout( 12, 0 ) );
        footer.setBorder( BorderFactory.createCompoundBorder( BorderFactory
                .createMatteBorder( 1, 0, 0, 0, borderColor() ), BorderFactory.createEmptyBorder( 8, 14, 10, 14 ) ) );
        _status = new JLabel( " " );
        _status.setForeground( mutedColor() );
        footer.add( _status, BorderLayout.CENTER );
        final JPanel buttons = new JPanel( new FlowLayout( FlowLayout.RIGHT, 8, 0 ) );
        final JButton close = new JButton( "Close" );
        close.addActionListener( e -> requestClose() );
        if ( _form.isEditable() ) {
            _write = new JButton( "Write to Tree" );
            _write.setToolTipText( "Validate and write every change in this window to the tree (one undo step)" );
            _write.addActionListener( e -> writeNow() );
            buttons.add( _write );
            getRootPane().setDefaultButton( _write );
        }
        else {
            _write = null;
        }
        buttons.add( close );
        footer.add( buttons, BorderLayout.EAST );
        getContentPane().add( footer, BorderLayout.SOUTH );
        // -- keys: Esc / Cmd-W close, Cmd-Enter writes --
        final int menu_mask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
        bind( KeyStroke.getKeyStroke( KeyEvent.VK_ESCAPE, 0 ), "close" );
        bind( KeyStroke.getKeyStroke( KeyEvent.VK_W, menu_mask ), "close" );
        getRootPane().getActionMap().put( "close", action( this::requestClose ) );
        if ( _form.isEditable() ) {
            bind( KeyStroke.getKeyStroke( KeyEvent.VK_ENTER, menu_mask ), "write" );
            getRootPane().getActionMap().put( "write", action( this::writeNow ) );
        }
        _form.addChangeListener( this::refreshState );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                requestClose();
            }
        } );
        refreshState();
        sizeAndPlace();
        setVisible( true );
    }

    private void bind( final KeyStroke ks, final String action ) {
        getRootPane().getInputMap( JComponent.WHEN_IN_FOCUSED_WINDOW ).put( ks, action );
    }

    private static AbstractAction action( final Runnable r ) {
        return new AbstractAction() {

            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed( final java.awt.event.ActionEvent e ) {
                r.run();
            }
        };
    }

    NodeDataForm getForm() {
        return _form;
    }

    boolean isEditable() {
        return _form.isEditable();
    }

    boolean isDirty() {
        return _form.isDirty();
    }

    /** The Write button's action (also Cmd-Enter and the default button): write if valid, else show why not. */
    boolean writeNow() {
        final boolean ok = _form.write();
        refreshState();
        if ( ok ) {
            _status.setText( "Written to the tree." );
            _status.setForeground( mutedColor() );
        }
        return ok;
    }

    /** Close, asking first when there are unsaved edits (Write and Close / Discard Changes / Cancel). */
    void requestClose() {
        if ( _form.isDirty() ) { // never true in VIEW mode
            final Object[] options = { "Write and Close", "Discard Changes", "Cancel" };
            final int r = JOptionPane.showOptionDialog( this,
                                                        "This node has changes that have not been written to the tree.",
                                                        "Unsaved Changes",
                                                        JOptionPane.YES_NO_CANCEL_OPTION,
                                                        JOptionPane.WARNING_MESSAGE,
                                                        null,
                                                        options,
                                                        options[ 0 ] );
            if ( r == 0 ) {
                if ( !writeNow() ) {
                    return; // a validation problem -- it is now shown in the status line; stay open
                }
            }
            else if ( r != 1 ) {
                return; // Cancel (or the dialog was dismissed)
            }
        }
        close();
    }

    /** Closes unconditionally: releases the tree panel's slot and disposes the window. */
    void close() {
        if ( _tree_panel != null ) {
            _tree_panel.removeEditNodeFrame( this );
        }
        dispose();
    }

    /** For tests: the status line's current text. */
    String statusTextForTest() {
        return _status.getText();
    }

    /** For tests: the Write button (null in VIEW mode). */
    JButton writeButtonForTest() {
        return _write;
    }

    private void refreshState() {
        final boolean dirty = _form.isDirty();
        if ( dirty != isDocumentModifiedMark() ) { // title + OS mark go to the native peer: only on a real flip
            setTitle( ( dirty ? "• " : "" ) + _base_title );
            getRootPane().putClientProperty( "Window.documentModified", dirty );
        }
        if ( !_form.isEditable() ) {
            return;
        }
        final List<Problem> problems = _form.problems();
        final String text;
        final Color color;
        if ( !problems.isEmpty() ) {
            text = problems.get( 0 ).message;
            color = errorColor();
        }
        else if ( dirty ) {
            text = "Unsaved changes";
            color = mutedColor();
        }
        else {
            text = " ";
            color = mutedColor();
        }
        if ( !text.equals( _status.getText() ) ) {
            _status.setText( text );
            _status.setForeground( color );
        }
        _write.setEnabled( dirty && problems.isEmpty() );
    }

    /** Pack, then keep the window inside the usable screen and cascade it off the tree panel by its slot index. */
    private void sizeAndPlace() {
        pack();
        final int em = getFont() != null ? getFont().getSize() : 13;
        final Rectangle usable = GraphicsEnvironment.getLocalGraphicsEnvironment().getMaximumWindowBounds();
        final Dimension pref = getSize();
        final int w = Math.min( Math.max( pref.width, em * 44 ), (int) ( usable.width * 0.9 ) );
        final int h = Math.min( Math.max( pref.height, em * 24 ), (int) ( usable.height * 0.9 ) );
        setSize( w, h );
        setMinimumSize( new Dimension( Math.min( em * 30, w ), Math.min( em * 16, h ) ) );
        if ( ( _tree_panel != null ) && _tree_panel.isShowing() ) {
            setLocationRelativeTo( _tree_panel );
        }
        else {
            setLocationRelativeTo( null );
        }
        final Point p = getLocation();
        p.translate( _index * 24, _index * 24 );
        p.x = Math.max( usable.x, Math.min( p.x, usable.x + usable.width - w ) );
        p.y = Math.max( usable.y, Math.min( p.y, usable.y + usable.height - h ) );
        setLocation( p );
    }

    private static Color mutedColor() {
        final Color c = UIManager.getColor( "Label.disabledForeground" );
        return ( c != null ) ? c : Color.GRAY;
    }

    private static Color borderColor() {
        final Color c = UIManager.getColor( "Component.borderColor" );
        return ( c != null ) ? c : Color.LIGHT_GRAY;
    }

    /** Error TEXT colour (FlatLaf's action red reads on both themes; the border red is too dim on dark). */
    private static Color errorColor() {
        Color c = UIManager.getColor( "Actions.Red" );
        if ( c == null ) {
            c = UIManager.getColor( "Component.error.focusedBorderColor" );
        }
        return ( c != null ) ? c : new Color( 0xD0342C );
    }

    /** Whether the frame's root pane shows the OS "document modified" mark (macOS close-button dot). */
    boolean isDocumentModifiedMark() {
        final JRootPane rp = getRootPane();
        return Boolean.TRUE.equals( rp.getClientProperty( "Window.documentModified" ) );
    }
}
