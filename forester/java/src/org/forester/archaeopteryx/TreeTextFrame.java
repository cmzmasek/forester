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

import static org.forester.archaeopteryx.FormWidgets.accentColor;
import static org.forester.archaeopteryx.FormWidgets.borderColor;
import static org.forester.archaeopteryx.FormWidgets.monoFont;
import static org.forester.archaeopteryx.FormWidgets.mutedColor;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.Toolkit;
import java.awt.datatransfer.StringSelection;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.KeyStroke;
import javax.swing.Timer;
import javax.swing.text.BadLocationException;
import javax.swing.text.DefaultCaret;
import javax.swing.text.DefaultHighlighter;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;

import org.forester.archaeopteryx.FormWidgets.Header;
import org.forester.archaeopteryx.TreeText.Format;
import org.forester.archaeopteryx.TreeText.Span;
import org.forester.phylogeny.Phylogeny;

/**
 * The "Tree as Text" window (View &gt; as phyloXML / as Newick / as Nexus): one window per tree panel with a
 * format switcher, monospace text with the markup muted so the names stand out, a find field (Cmd-F, Enter /
 * Shift-Enter step through the hits), a wrap toggle, Copy and Save As. The text of a format is generated when
 * its tab is first shown and again after the tree changed (the panel marks the window stale; regeneration is
 * coalesced and happens when the window is looked at).
 */
final class TreeTextFrame extends JFrame {

    private static final long          serialVersionUID = 1L;
    /** How long after the last change the window waits before regenerating the visible text. */
    static final int                   REFRESH_DELAY_MS = 300;
    /** Hits beyond this are not highlighted (a single-letter search in a huge file). */
    static final int                   MAX_HITS         = 5000;
    private final TreePanel            _tree_panel;
    private final Header               _header;
    private final JTabbedPane          _tabs;
    private final Map<Format, JTextPane> _panes         = new EnumMap<>( Format.class );
    private final Map<Format, String>  _texts           = new EnumMap<>( Format.class );
    private final Map<Format, Boolean> _wrap            = new EnumMap<>( Format.class );
    private final JTextField           _find;
    private final JLabel               _find_count;
    private final JCheckBox            _wrap_box;
    private final JLabel               _status;
    private final Timer                _refresh;
    private final List<Integer>        _hits            = new ArrayList<>();
    private int                        _hit             = -1;
    private boolean                    _stale           = false;

    TreeTextFrame( final TreePanel tp, final Format initial ) {
        _tree_panel = tp;
        setDefaultCloseOperation( DO_NOTHING_ON_CLOSE );
        getContentPane().setLayout( new BorderLayout() );
        _header = new Header( titleText(), "" );
        getContentPane().add( _header, BorderLayout.NORTH );
        _tabs = new JTabbedPane();
        for( final Format f : Format.values() ) {
            final JTextPane pane = newPane();
            _panes.put( f, pane );
            _wrap.put( f, f.wrapsByDefault() );
            final JScrollPane sp = new JScrollPane( pane );
            sp.setBorder( BorderFactory.createMatteBorder( 1, 0, 1, 0, borderColor() ) );
            sp.getVerticalScrollBar().setUnitIncrement( 16 );
            _tabs.addTab( f.label, sp );
        }
        _tabs.setSelectedIndex( initial.ordinal() );
        _tabs.addChangeListener( e -> showFormat( format() ) );
        getContentPane().add( _tabs, BorderLayout.CENTER );
        // -- footer: find + wrap (left), status (centre), buttons (right) --
        final JPanel footer = new JPanel( new BorderLayout( 12, 0 ) );
        footer.setBorder( BorderFactory.createEmptyBorder( 8, 14, 10, 14 ) );
        final JPanel left = new JPanel( new FlowLayout( FlowLayout.LEFT, 8, 0 ) );
        _find = new JTextField( 16 );
        _find.putClientProperty( "JTextField.placeholderText", "Find" );
        _find.putClientProperty( "JTextField.showClearButton", true );
        _find.putClientProperty( "JTextField.selectAllOnFocusPolicy", "never" );
        _find.getDocument().addDocumentListener( FormWidgets.onChange( this::search ) );
        _find.addActionListener( e -> step( ( e.getModifiers() & java.awt.event.ActionEvent.SHIFT_MASK ) != 0 ? -1 : 1 ) );
        left.add( _find );
        _find_count = new JLabel( " " );
        _find_count.setForeground( mutedColor() );
        left.add( _find_count );
        _wrap_box = new JCheckBox( "Wrap lines", initial.wrapsByDefault() );
        _wrap_box.setFocusable( false );
        _wrap_box.addActionListener( e -> setWrap( _wrap_box.isSelected() ) );
        left.add( _wrap_box );
        footer.add( left, BorderLayout.WEST );
        _status = new JLabel( " " );
        _status.setForeground( mutedColor() );
        footer.add( _status, BorderLayout.CENTER );
        final JPanel buttons = new JPanel( new FlowLayout( FlowLayout.RIGHT, 8, 0 ) );
        final JButton copy = new JButton( "Copy" );
        copy.setToolTipText( "Copy the whole text to the clipboard (Cmd-C copies a selection)" );
        copy.addActionListener( e -> copyAll() );
        buttons.add( copy );
        final JButton save = new JButton( "Save As…" );
        save.setToolTipText( "Write this text to a file" );
        save.addActionListener( e -> saveAs() );
        buttons.add( save );
        final JButton close = new JButton( "Close" );
        close.addActionListener( e -> close() );
        buttons.add( close );
        footer.add( buttons, BorderLayout.EAST );
        getContentPane().add( footer, BorderLayout.SOUTH );
        // -- keys --
        final int menu_mask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
        bind( KeyStroke.getKeyStroke( KeyEvent.VK_ESCAPE, 0 ), "close", this::close );
        bind( KeyStroke.getKeyStroke( KeyEvent.VK_W, menu_mask ), "close", this::close );
        bind( KeyStroke.getKeyStroke( KeyEvent.VK_F, menu_mask ), "find", () -> {
            _find.requestFocusInWindow();
            _find.selectAll();
        } );
        bind( KeyStroke.getKeyStroke( KeyEvent.VK_G, menu_mask ), "next", () -> step( 1 ) );
        bind( KeyStroke.getKeyStroke( KeyEvent.VK_G, menu_mask | KeyEvent.SHIFT_DOWN_MASK ), "prev", () -> step( -1 ) );
        _refresh = new Timer( REFRESH_DELAY_MS, e -> refreshNow() );
        _refresh.setRepeats( false );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                close();
            }

            @Override
            public void windowActivated( final WindowEvent e ) {
                if ( _stale ) {
                    _refresh.stop();
                    refreshNow();
                }
            }
        } );
        showFormat( initial );
        pack();
        final int em = getFont() != null ? getFont().getSize() : 13;
        setSize( em * 56, em * 38 );
        setMinimumSize( new java.awt.Dimension( em * 30, em * 18 ) );
        if ( ( tp != null ) && tp.isShowing() ) {
            setLocationRelativeTo( tp );
        }
        else {
            setLocationRelativeTo( null );
        }
        setVisible( true );
    }

    // ------------------------------------------------------------------ API (TreePanel + tests)
    Format format() {
        return Format.values()[ _tabs.getSelectedIndex() ];
    }

    /** Switches to {@code f} (generating its text if needed) and brings the window forward. */
    void show( final Format f ) {
        _tabs.setSelectedIndex( f.ordinal() );
        toFront();
        requestFocus();
    }

    /** The tree changed (or was replaced): regenerate soon (now, if the window is active). Cheap to call often. */
    void markStale() {
        _stale = true;
        _texts.clear();
        _refresh.restart();
    }

    /** Regenerates the visible format's text right away (the other formats regenerate when shown). */
    void refreshNow() {
        _refresh.stop();
        _stale = false;
        _texts.clear();
        _header.setTitle( titleText() );
        showFormat( format() );
    }

    /** The text currently shown. */
    String textForTest() {
        return _panes.get( format() ).getText();
    }

    String subtitleForTest() {
        return _header.getSubtitle();
    }

    boolean isStaleForTest() {
        return _stale;
    }

    int hitCountForTest() {
        return _hits.size();
    }

    int currentHitForTest() {
        return _hit;
    }

    JTextField findFieldForTest() {
        return _find;
    }

    void setFindTextForTest( final String s ) {
        _find.setText( s );
    }

    void stepForTest( final int dir ) {
        step( dir );
    }

    boolean isWrapForTest() {
        return _wrap.get( format() );
    }

    void setWrapForTest( final boolean wrap ) {
        _wrap_box.setSelected( wrap );
        setWrap( wrap );
    }

    /** The tint colour at {@code offset} of the shown text (Label.foreground for plain data). */
    Color colorAtForTest( final int offset ) {
        final StyledDocument doc = _panes.get( format() ).getStyledDocument();
        final Color c = StyleConstants.getForeground( doc.getCharacterElement( offset ).getAttributes() );
        return c;
    }

    void close() {
        _refresh.stop();
        if ( _tree_panel != null ) {
            _tree_panel.treeTextFrameClosed( this );
        }
        dispose();
    }

    // ------------------------------------------------------------------ text
    private String titleText() {
        final Phylogeny phy = ( _tree_panel != null ) ? _tree_panel.getPhylogeny() : null;
        final String name = ( phy != null ) ? phy.getName() : null;
        return ( ( name == null ) || name.trim().isEmpty() ) ? TreePropertiesForm.UNTITLED : name.trim();
    }

    private void showFormat( final Format f ) {
        final JTextPane pane = _panes.get( f );
        String text = _texts.get( f );
        if ( text == null ) {
            final Phylogeny phy = ( _tree_panel != null ) ? _tree_panel.getPhylogeny() : null;
            final org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE style = ( ( _tree_panel != null )
                    && ( _tree_panel.getOptions() != null ) ) ? _tree_panel.getOptions().getNhConversionSupportValueStyle()
                            : org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE;
            text = TreeText.render( phy, f, style );
            _texts.put( f, text );
            fill( pane, text, f );
        }
        setTitle( "Tree as " + f.label + ": " + titleText() );
        _header.setSubtitle( f.label + " · " + TreeText.sizeText( text ) );
        _wrap_box.setSelected( _wrap.get( f ) );
        _status.setText( " " );
        search();
    }

    /** Puts {@code text} into {@code pane}, tinted by {@link TreeText#spans}. */
    private void fill( final JTextPane pane, final String text, final Format f ) {
        final StyledDocument doc = pane.getStyledDocument();
        try {
            doc.remove( 0, doc.getLength() );
            doc.insertString( 0, text, plainStyle() );
            final List<Span> spans = TreeText.spans( text, f );
            if ( !spans.isEmpty() ) {
                final SimpleAttributeSet markup = new SimpleAttributeSet();
                StyleConstants.setForeground( markup, mutedColor() );
                final SimpleAttributeSet comment = new SimpleAttributeSet( markup );
                StyleConstants.setItalic( comment, true );
                final SimpleAttributeSet keyword = new SimpleAttributeSet();
                StyleConstants.setForeground( keyword, accentColor() );
                StyleConstants.setBold( keyword, true );
                for( final Span s : spans ) {
                    final SimpleAttributeSet a;
                    switch ( s.kind ) {
                        case COMMENT:
                            a = comment;
                            break;
                        case KEYWORD:
                            a = keyword;
                            break;
                        default:
                            a = markup;
                    }
                    doc.setCharacterAttributes( s.start, s.end - s.start, a, false );
                }
            }
        }
        catch ( final BadLocationException e ) {
            pane.setText( text );
        }
        pane.setCaretPosition( 0 );
    }

    private static SimpleAttributeSet plainStyle() {
        final SimpleAttributeSet plain = new SimpleAttributeSet();
        final Color fg = javax.swing.UIManager.getColor( "TextArea.foreground" );
        StyleConstants.setForeground( plain, ( fg != null ) ? fg : Color.BLACK );
        return plain;
    }

    private JTextPane newPane() {
        final JTextPane pane = new JTextPane() {

            private static final long serialVersionUID = 1L;

            @Override
            public boolean getScrollableTracksViewportWidth() {
                final Boolean wrap = _wrap.get( paneFormat( this ) ); // null while the window is being built
                if ( ( wrap == null ) || wrap ) {
                    return true;
                }
                final Component parent = getParent();
                return ( parent == null ) || ( getUI().getPreferredSize( this ).width <= parent.getWidth() );
            }
        };
        pane.setEditable( false );
        pane.setFont( monoFont( pane.getFont() ) );
        pane.setBorder( BorderFactory.createEmptyBorder( 8, 10, 8, 10 ) );
        if ( pane.getCaret() instanceof DefaultCaret ) {
            ( (DefaultCaret) pane.getCaret() ).setUpdatePolicy( DefaultCaret.NEVER_UPDATE );
        }
        return pane;
    }

    private Format paneFormat( final JTextPane pane ) {
        for( final Map.Entry<Format, JTextPane> e : _panes.entrySet() ) {
            if ( e.getValue() == pane ) {
                return e.getKey();
            }
        }
        return Format.PHYLOXML;
    }

    private void setWrap( final boolean wrap ) {
        _wrap.put( format(), wrap );
        final JTextPane pane = _panes.get( format() );
        pane.revalidate();
        pane.repaint();
    }

    // ------------------------------------------------------------------ find
    /** Recomputes the hits of the find text in the shown text and highlights them (from scratch). */
    private void search() {
        final JTextPane pane = _panes.get( format() );
        pane.getHighlighter().removeAllHighlights();
        _hits.clear();
        _hit = -1;
        final String needle = _find.getText();
        if ( needle.isEmpty() ) {
            _find_count.setText( " " );
            return;
        }
        final String hay = pane.getText().toLowerCase( Locale.ROOT );
        final String n = needle.toLowerCase( Locale.ROOT );
        int from = 0;
        int total = 0;
        final Color c = accentColor();
        final DefaultHighlighter.DefaultHighlightPainter painter = new DefaultHighlighter.DefaultHighlightPainter(
                new Color( c.getRed(), c.getGreen(), c.getBlue(), 70 ) );
        while( true ) {
            final int i = hay.indexOf( n, from );
            if ( i < 0 ) {
                break;
            }
            ++total;
            if ( _hits.size() < MAX_HITS ) {
                _hits.add( i );
                try {
                    pane.getHighlighter().addHighlight( i, i + n.length(), painter );
                }
                catch ( final BadLocationException e ) {
                    // a stale offset (the document changed underneath) -- skip this hit
                }
            }
            from = i + Math.max( 1, n.length() );
        }
        if ( total == 0 ) {
            _find_count.setText( "no matches" );
            _find_count.setForeground( FormWidgets.errorColor() );
        }
        else {
            _find_count.setForeground( mutedColor() );
            step( 1 ); // select the first hit
        }
    }

    /** Selects the next ({@code dir} &gt; 0) or previous hit, wrapping around. */
    private void step( final int dir ) {
        if ( _hits.isEmpty() ) {
            return;
        }
        _hit = ( ( _hit + dir ) % _hits.size() + _hits.size() ) % _hits.size();
        final int at = _hits.get( _hit );
        final JTextPane pane = _panes.get( format() );
        final int len = _find.getText().length();
        try {
            pane.scrollRectToVisible( pane.modelToView2D( at ).getBounds() );
        }
        catch ( final BadLocationException e ) {
            // ignore: the selection below still moves
        }
        pane.setCaretPosition( at );
        pane.moveCaretPosition( at + len );
        _find_count.setText( ( _hit + 1 ) + " of " + _hits.size() );
    }

    // ------------------------------------------------------------------ copy / save
    private void copyAll() {
        try {
            Toolkit.getDefaultToolkit().getSystemClipboard().setContents( new StringSelection( textForTest() ), null );
            _status.setText( "Copied to the clipboard." );
        }
        catch ( final IllegalStateException e ) {
            _status.setText( "The clipboard is not available right now." );
        }
    }

    private void saveAs() {
        final Phylogeny phy = ( _tree_panel != null ) ? _tree_panel.getPhylogeny() : null;
        final File tree_file = ( _tree_panel != null ) ? _tree_panel.getTreeFile() : null;
        final MainFrame mf = ( ( _tree_panel != null ) && ( _tree_panel.getMainPanel() != null ) )
                ? _tree_panel.getMainPanel().getMainFrame() : null;
        final JFileChooser fc = new JFileChooser();
        fc.setDialogTitle( "Save " + format().label + " text as" );
        if ( mf != null ) {
            fc.setCurrentDirectory( mf.getCurrentDir( DirectoryPreferences.Category.SAVE ) );
        }
        else if ( tree_file != null ) {
            fc.setCurrentDirectory( tree_file.getParentFile() );
        }
        fc.setSelectedFile( new File( fc.getCurrentDirectory(),
                                      TreeText.suggestedFileName( phy, tree_file, format() ) ) );
        if ( fc.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) {
            return;
        }
        final File out = fc.getSelectedFile();
        if ( out.exists() ) {
            final int r = JOptionPane.showConfirmDialog( this, "Overwrite \"" + out.getName() + "\"?", "File exists",
                                                         JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE );
            if ( r != JOptionPane.OK_OPTION ) {
                return;
            }
        }
        try {
            Files.write( out.toPath(), textForTest().getBytes( StandardCharsets.UTF_8 ) );
            _status.setText( "Saved to " + out.getName() + "." );
            if ( mf != null ) {
                mf.setCurrentDir( DirectoryPreferences.Category.SAVE, out.getParentFile() );
            }
        }
        catch ( final IOException e ) {
            JOptionPane.showMessageDialog( this, "Could not write \"" + out + "\":\n" + e.getMessage(), "Save failed",
                                           JOptionPane.ERROR_MESSAGE );
        }
    }

    private void bind( final KeyStroke ks, final String name, final Runnable r ) {
        getRootPane().getInputMap( JComponent.WHEN_IN_FOCUSED_WINDOW ).put( ks, name );
        getRootPane().getActionMap().put( name, FormWidgets.action( r ) );
    }
}
