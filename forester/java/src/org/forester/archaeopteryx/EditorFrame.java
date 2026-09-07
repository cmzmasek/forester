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

import static org.forester.archaeopteryx.FormWidgets.borderColor;
import static org.forester.archaeopteryx.FormWidgets.errorColor;
import static org.forester.archaeopteryx.FormWidgets.mutedColor;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
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

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRootPane;
import javax.swing.KeyStroke;

import org.forester.archaeopteryx.NodeDataDraft.Problem;

/**
 * The window chrome shared by the editor windows (node data, tree properties): the form in the middle, a footer
 * with a status line (the first validation problem in red, or "Unsaved changes") and the buttons -- <b>Write to
 * Tree</b> (the default button, and the only thing that ever touches the tree) plus Close. Esc and Cmd-W close,
 * Cmd-Enter writes. Closing with unsaved edits asks: Write and Close / Discard Changes / Cancel. The title gets a
 * leading bullet and the OS "document modified" mark while there are unsaved edits. A read-only form gets only
 * the Close button.
 */
abstract class EditorFrame extends JFrame {

    /** What the chrome needs from the form it wraps. */
    interface Form {

        JComponent component();

        boolean isEditable();

        boolean isDirty();

        List<Problem> problems();

        /** Writes if valid; returns false (and touches nothing) on a validation problem. */
        boolean write();

        void addChangeListener( Runnable r );

        /** A reason the form cannot write right now (shown in the status line, Write disabled), or null. */
        default String notice() {
            return null;
        }
    }

    private static final long serialVersionUID = 1L;
    private final Form        _form;
    private final JLabel      _status;
    private final JButton     _write;
    private final String      _unsaved_message;
    private String            _base_title;

    /**
     * @param form            the form to wrap
     * @param base_title      the window title (a bullet is prefixed while dirty)
     * @param unsaved_message the first line of the close confirmation ("This node has changes ...")
     */
    EditorFrame( final Form form, final String base_title, final String unsaved_message ) {
        _form = form;
        _base_title = base_title;
        _unsaved_message = unsaved_message;
        setTitle( base_title );
        setDefaultCloseOperation( DO_NOTHING_ON_CLOSE );
        getContentPane().setLayout( new BorderLayout() );
        getContentPane().add( form.component(), BorderLayout.CENTER );
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
        if ( form.isEditable() ) {
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
        getRootPane().getActionMap().put( "close", FormWidgets.action( this::requestClose ) );
        if ( form.isEditable() ) {
            bind( KeyStroke.getKeyStroke( KeyEvent.VK_ENTER, menu_mask ), "write" );
            getRootPane().getActionMap().put( "write", FormWidgets.action( this::writeNow ) );
        }
        form.addChangeListener( this::refreshState );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                requestClose();
            }
        } );
        refreshState();
    }

    /** Maps {@code ks} to the named root-pane action while the window is focused. */
    protected void bind( final KeyStroke ks, final String action ) {
        getRootPane().getInputMap( JComponent.WHEN_IN_FOCUSED_WINDOW ).put( ks, action );
    }

    /** Called once, from {@link #close()}, after the window is done: release whatever tracks it. */
    protected abstract void onClosed();

    boolean isEditable() {
        return _form.isEditable();
    }

    boolean isDirty() {
        return _form.isDirty();
    }

    /** Renames the window (the dirty bullet is re-applied). */
    void setBaseTitle( final String base_title ) {
        _base_title = base_title;
        setTitle( ( isDocumentModifiedMark() ? "• " : "" ) + _base_title );
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
        if ( _form.isDirty() ) { // never true for a read-only form
            final Object[] options = { "Write and Close", "Discard Changes", "Cancel" };
            final int r = JOptionPane.showOptionDialog( this,
                                                        _unsaved_message,
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

    /** Closes unconditionally: tells the owner ({@link #onClosed()}) and disposes the window. */
    void close() {
        onClosed();
        dispose();
    }

    /** For tests: the status line's current text. */
    String statusTextForTest() {
        return _status.getText();
    }

    /** For tests: the Write button (null for a read-only form). */
    JButton writeButtonForTest() {
        return _write;
    }

    /** Re-derives the title mark, the status line and the Write button from the form's state. */
    protected void refreshState() {
        final boolean dirty = _form.isDirty();
        if ( dirty != isDocumentModifiedMark() ) { // title + OS mark go to the native peer: only on a real flip
            setTitle( ( dirty ? "• " : "" ) + _base_title );
            getRootPane().putClientProperty( "Window.documentModified", dirty );
        }
        final String notice = _form.notice();
        if ( notice != null ) { // the form cannot write at all right now (e.g. its node is gone) -- any mode
            if ( !notice.equals( _status.getText() ) ) {
                _status.setText( notice );
                _status.setForeground( errorColor() );
            }
            if ( _write != null ) {
                _write.setEnabled( false );
            }
            return;
        }
        if ( !_form.isEditable() ) {
            if ( !" ".equals( _status.getText() ) ) { // a notice that has since been lifted
                _status.setText( " " );
                _status.setForeground( mutedColor() );
            }
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

    /**
     * Pack, then keep the window inside the usable screen, place it over {@code anchor} (or the screen centre) and
     * cascade it by {@code cascade_index} steps. {@code min_em_w}/{@code min_em_h} are the minimum size in ems.
     */
    protected void sizeAndPlace( final Component anchor, final int cascade_index, final int min_em_w,
                                 final int min_em_h ) {
        pack();
        final int em = getFont() != null ? getFont().getSize() : 13;
        final Rectangle usable = GraphicsEnvironment.getLocalGraphicsEnvironment().getMaximumWindowBounds();
        final Dimension pref = getSize();
        final int w = Math.min( Math.max( pref.width, em * min_em_w ), (int) ( usable.width * 0.9 ) );
        final int h = Math.min( Math.max( pref.height, em * min_em_h ), (int) ( usable.height * 0.9 ) );
        setSize( w, h );
        setMinimumSize( new Dimension( Math.min( em * 30, w ), Math.min( em * 16, h ) ) );
        if ( ( anchor != null ) && anchor.isShowing() ) {
            setLocationRelativeTo( anchor );
        }
        else {
            setLocationRelativeTo( null );
        }
        final Point p = getLocation();
        p.translate( cascade_index * 24, cascade_index * 24 );
        p.x = Math.max( usable.x, Math.min( p.x, usable.x + usable.width - w ) );
        p.y = Math.max( usable.y, Math.min( p.y, usable.y + usable.height - h ) );
        setLocation( p );
    }

    /** Whether the frame's root pane shows the OS "document modified" mark (macOS close-button dot). */
    boolean isDocumentModifiedMark() {
        final JRootPane rp = getRootPane();
        return Boolean.TRUE.equals( rp.getClientProperty( "Window.documentModified" ) );
    }
}
