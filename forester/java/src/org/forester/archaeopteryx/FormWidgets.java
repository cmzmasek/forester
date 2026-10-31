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
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Rectangle;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.DefaultCaret;
import javax.swing.text.JTextComponent;

/**
 * The shared building blocks of the data windows (node data, tree properties, tree as text): theme colours, the
 * header (title + muted subtitle), collapsible {@link Section}s, the label/value {@link Grid}, read-only values
 * that still read as text, editable fields with placeholders, and the width-tracking scrolling page. One look for
 * all of them, defined once.
 */
final class FormWidgets {

    /** The header title's size relative to the body font: a heading, not a banner. */
    static final float TITLE_SCALE      = 1.15f;
    /** How far a section's body is indented under its heading. */
    static final int   SECTION_INDENT   = 18;
    /** The preferred width a read-only value reports (it stretches to the row anyway). */
    static final int   VIEW_VALUE_WIDTH = 120;

    private FormWidgets() {
    }

    // ------------------------------------------------------------------ colours + fonts
    static Color mutedColor() {
        final Color c = UIManager.getColor( "Label.disabledForeground" );
        return ( c != null ) ? c : Color.GRAY;
    }

    static Color borderColor() {
        final Color c = UIManager.getColor( "Component.borderColor" );
        return ( c != null ) ? c : Color.LIGHT_GRAY;
    }

    static Color accentColor() {
        final Color c = UIManager.getColor( "Component.accentColor" );
        return ( c != null ) ? c : new Color( 0x2675BF );
    }

    /** Error TEXT colour (FlatLaf's action red reads on both themes; the border red is too dim on dark). */
    static Color errorColor() {
        Color c = UIManager.getColor( "Actions.Red" );
        if ( c == null ) {
            c = UIManager.getColor( "Component.error.focusedBorderColor" );
        }
        return ( c != null ) ? c : new Color( 0xD0342C );
    }

    static Font monoFont( final Font base ) {
        return new Font( Font.MONOSPACED, Font.PLAIN, base.getSize() );
    }

    // ------------------------------------------------------------------ small helpers
    /** A document listener that runs {@code r} on every kind of change. */
    static DocumentListener onChange( final Runnable r ) {
        return new DocumentListener() {

            @Override
            public void insertUpdate( final DocumentEvent e ) {
                r.run();
            }

            @Override
            public void removeUpdate( final DocumentEvent e ) {
                r.run();
            }

            @Override
            public void changedUpdate( final DocumentEvent e ) {
                r.run();
            }
        };
    }

    static AbstractAction action( final Runnable r ) {
        return new AbstractAction() {

            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed( final java.awt.event.ActionEvent e ) {
                r.run();
            }
        };
    }

    // ------------------------------------------------------------------ header
    /** The window header: a bold title only slightly larger than the body text, over a muted one-line subtitle. */
    static final class Header extends JPanel {

        private static final long serialVersionUID = 1L;
        private final JLabel      _title;
        private final JLabel      _subtitle;

        Header( final String title, final String subtitle ) {
            super( new BorderLayout( 0, 2 ) );
            setBorder( BorderFactory.createEmptyBorder( 12, 14, 8, 14 ) );
            _title = new JLabel( title );
            _title.setFont( _title.getFont().deriveFont( Font.BOLD, _title.getFont().getSize2D() * TITLE_SCALE ) );
            _subtitle = new JLabel( subtitle );
            _subtitle.setForeground( mutedColor() );
            add( _title, BorderLayout.NORTH );
            add( _subtitle, BorderLayout.CENTER );
        }

        void setTitle( final String title ) {
            _title.setText( title );
        }

        void setSubtitle( final String subtitle ) {
            _subtitle.setText( subtitle );
        }

        String getTitle() {
            return _title.getText();
        }

        String getSubtitle() {
            return _subtitle.getText();
        }
    }

    // ------------------------------------------------------------------ fields + values
    /** An editable single-line field: placeholder hint, and a click puts the caret where you click (no select-all). */
    static JTextField editField( final String value, final String placeholder ) {
        final JTextField tf = new JTextField( value, 10 ); // a column count: the TEXT must not size the grid
        if ( placeholder != null ) {
            tf.putClientProperty( "JTextField.placeholderText", placeholder );
        }
        tf.putClientProperty( "JTextField.selectAllOnFocusPolicy", "never" );
        return tf;
    }

    /** A read-only value that still looks like text (selectable, copyable), not like a disabled field. */
    static JTextComponent viewValue( final String value, final boolean multiline ) {
        final JTextComponent tc;
        if ( multiline ) {
            final JTextArea ta = new JTextArea( value ) {

                private static final long serialVersionUID = 1L;

                @Override
                public Dimension getPreferredSize() {
                    // the HEIGHT follows the wrapped text at the current width; the WIDTH must never be the text's
                    // (a wrapped area reports its last laid-out width, which would ratchet the grid wider than the
                    // window -- and GridBagLayout answers an over-wide grid by collapsing every row to its minimum)
                    final Dimension d = super.getPreferredSize();
                    return new Dimension( VIEW_VALUE_WIDTH, d.height );
                }
            };
            ta.setLineWrap( true );
            ta.setWrapStyleWord( true );
            tc = ta;
        }
        else {
            tc = new JTextField( value ) {

                private static final long serialVersionUID = 1L;

                @Override
                public Dimension getPreferredSize() {
                    final Dimension d = super.getPreferredSize();
                    return new Dimension( Math.min( d.width, VIEW_VALUE_WIDTH ), d.height );
                }
            };
        }
        tc.setFocusable( false ); // read-only: no caret, no focus ring; the value is display, not input
        tc.setEditable( false );
        tc.setOpaque( false );
        tc.setBorder( BorderFactory.createEmptyBorder( 2, 0, 2, 0 ) );
        tc.setForeground( UIManager.getColor( "Label.foreground" ) );
        if ( tc.getCaret() instanceof DefaultCaret ) {
            // a read-only value must never scroll the page to itself (a long mol seq would open the window scrolled)
            ( (DefaultCaret) tc.getCaret() ).setUpdatePolicy( DefaultCaret.NEVER_UPDATE );
        }
        tc.setCaretPosition( 0 );
        return tc;
    }

    /** The client property a link value carries: the address it opens (read by tests and by nothing else). */
    static final String LINK_URL_PROPERTY = "aptx.link.url";

    /** The theme's link colour (FlatLaf's {@code Component.linkColor}, readable on light and dark). */
    static Color linkColor() {
        final Color c = UIManager.getColor( "Component.linkColor" );
        return ( c != null ) ? c : accentColor();
    }

    /** A plain single left click: no other button, no modifier, not the second of a double click. */
    static boolean isPlainSingleLeftClick( final MouseEvent e ) {
        return javax.swing.SwingUtilities.isLeftMouseButton( e ) && ( e.getClickCount() == 1 )
                && ( ( e.getModifiersEx() & ( java.awt.event.InputEvent.SHIFT_DOWN_MASK
                        | java.awt.event.InputEvent.CTRL_DOWN_MASK | java.awt.event.InputEvent.ALT_DOWN_MASK
                        | java.awt.event.InputEvent.META_DOWN_MASK ) ) == 0 );
    }

    /**
     * A read-only value ({@link #viewValue}) that opens {@code url} on a plain left click ON ITS TEXT: the link
     * colour, an underline and a hand cursor while the pointer is over the text, and the address as its tooltip.
     * The component itself is stretched across its row, so the blank space beside or below the text is not the
     * link ({@link #isOverText}). The colour is re-applied after every look-and-feel update (a theme switch
     * re-installs the theme's text colour). With a null {@code url} it is the plain value, so a caller need not
     * branch.
     */
    static JTextComponent linkValue( final String text, final String url, final boolean multiline,
                                     final java.util.function.Consumer<String> open ) {
        final JTextComponent tc = viewValue( text, multiline );
        if ( url == null ) {
            return tc;
        }
        tc.setForeground( plainColor( linkColor() ) );
        tc.addPropertyChangeListener( "UI", e -> tc.setForeground( plainColor( linkColor() ) ) );
        tc.setToolTipText( url );
        tc.putClientProperty( LINK_URL_PROPERTY, url );
        final Cursor base = tc.getCursor();
        final MouseAdapter m = new MouseAdapter() {

            private Font _plain; // non-null while the pointer is over the text

            @Override
            public void mouseClicked( final MouseEvent e ) {
                if ( isPlainSingleLeftClick( e ) && isOverText( tc, e.getPoint() ) ) {
                    open.accept( url );
                }
            }

            @Override
            public void mouseMoved( final MouseEvent e ) {
                hover( isOverText( tc, e.getPoint() ) );
            }

            @Override
            public void mouseExited( final MouseEvent e ) {
                hover( false );
            }

            private void hover( final boolean over ) {
                if ( over == ( _plain != null ) ) {
                    return;
                }
                if ( over ) {
                    _plain = tc.getFont();
                    final java.util.Map<java.awt.font.TextAttribute, Object> a = new java.util.HashMap<>();
                    a.put( java.awt.font.TextAttribute.UNDERLINE, java.awt.font.TextAttribute.UNDERLINE_ON );
                    tc.setFont( _plain.deriveFont( a ) );
                    tc.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
                }
                else {
                    tc.setFont( _plain );
                    tc.setCursor( base );
                    _plain = null;
                }
            }
        };
        tc.addMouseListener( m );
        tc.addMouseMotionListener( m );
        return tc;
    }

    /** {@code c} as a plain colour: a look-and-feel update replaces a {@code UIResource} colour, never this. */
    private static Color plainColor( final Color c ) {
        return new Color( c.getRGB(), true );
    }

    /**
     * Whether {@code p} lies on a character of {@code tc}'s text -- not on the blank stretch beside the text, below
     * its last line, or in the border. False for a component that has not been laid out.
     */
    static boolean isOverText( final JTextComponent tc, final java.awt.Point p ) {
        final int n = tc.getDocument().getLength();
        if ( n == 0 ) {
            return false;
        }
        try {
            // the caret position nearest the pointer: the character under it is the one before or the one after
            final int pos = tc.viewToModel2D( p );
            for( int i = Math.max( 0, pos - 1 ); i <= Math.min( n - 1, pos ); ++i ) {
                final java.awt.geom.Rectangle2D box = characterBox( tc, i );
                if ( ( box != null ) && box.contains( p.getX(), p.getY() ) ) {
                    return true;
                }
            }
        }
        catch ( final javax.swing.text.BadLocationException e ) {
            return false;
        }
        return false;
    }

    /** The box character {@code i} of {@code tc} occupies: from its caret position to the next one on its line
     *  (at a line's end, its advance); null for a line break or before layout. */
    static java.awt.geom.Rectangle2D characterBox( final JTextComponent tc, final int i )
            throws javax.swing.text.BadLocationException {
        final String ch = tc.getDocument().getText( i, 1 );
        final java.awt.geom.Rectangle2D a = tc.modelToView2D( i );
        if ( ( a == null ) || ch.equals( "\n" ) ) {
            return null;
        }
        final java.awt.geom.Rectangle2D b = tc.modelToView2D( i + 1 );
        final double x0;
        final double x1;
        if ( ( b != null ) && ( Math.abs( b.getY() - a.getY() ) < 0.5 ) && ( b.getX() != a.getX() ) ) {
            x0 = Math.min( a.getX(), b.getX() );
            x1 = Math.max( a.getX(), b.getX() );
        }
        else { // the last character of a (wrapped) line
            x0 = a.getX();
            x1 = a.getX() + tc.getFontMetrics( tc.getFont() ).stringWidth( ch );
        }
        return new java.awt.geom.Rectangle2D.Double( x0, a.getY(), x1 - x0, a.getHeight() );
    }

    /**
     * Several values, one per line, each a link where {@code urls} has an address for it (null: plain text). Each
     * line wraps on its own, so a long address stays readable.
     */
    static JPanel linkList( final java.util.List<String> texts, final java.util.List<String> urls,
                            final java.util.function.Consumer<String> open ) {
        final JPanel p = new JPanel();
        p.setOpaque( false );
        p.setLayout( new BoxLayout( p, BoxLayout.Y_AXIS ) );
        for( int i = 0; i < texts.size(); ++i ) {
            final JTextComponent line = linkValue( texts.get( i ), urls.get( i ), true, open );
            line.setAlignmentX( Component.LEFT_ALIGNMENT );
            p.add( line );
        }
        return p;
    }

    /** A flat, left-aligned text button for in-page actions ("+ Add ..."). */
    static JButton linkButton( final String text, final Runnable action ) {
        final JButton b = new JButton( text );
        b.putClientProperty( "JButton.buttonType", "toolBarButton" );
        b.setHorizontalAlignment( SwingConstants.LEFT );
        b.setFocusable( false );
        b.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
        b.addActionListener( e -> action.run() );
        return b;
    }

    /** The small "×" button that removes a row or card. */
    static JButton removeButton( final String tooltip, final Runnable action ) {
        final JButton b = new JButton( "×" );
        b.putClientProperty( "JButton.buttonType", "toolBarButton" );
        b.setToolTipText( tooltip );
        b.setFocusable( false );
        b.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
        b.addActionListener( e -> action.run() );
        return b;
    }

    /** The small scrolling box around a multi-line field: its height is fixed by the area's row count (so a long
     *  value scrolls inside it), it hands the mouse wheel to the page when it cannot scroll itself, and it keeps
     *  that height even if the surrounding grid falls back to minimum sizes. */
    static JScrollPane areaScrollPane( final JTextArea ta, final boolean wrap ) {
        final JScrollPane sp = new JScrollPane( ta );
        sp.setVerticalScrollBarPolicy( JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED );
        sp.setHorizontalScrollBarPolicy( wrap ? JScrollPane.HORIZONTAL_SCROLLBAR_NEVER
                : JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED );
        final Dimension pref = sp.getPreferredSize();
        sp.setMinimumSize( new Dimension( 60, pref.height ) );
        forwardWheelWhenIdle( sp );
        return sp;
    }

    /** A nested scroll pane hands the wheel to its parent when it cannot scroll further itself. */
    static void forwardWheelWhenIdle( final JScrollPane inner ) {
        inner.addMouseWheelListener( e -> {
            final JScrollBar bar = inner.getVerticalScrollBar();
            final boolean can_scroll = bar.isVisible()
                    && ( ( e.getWheelRotation() < 0 ) ? ( bar.getValue() > bar.getMinimum() )
                            : ( bar.getValue() + bar.getVisibleAmount() < bar.getMaximum() ) );
            if ( !can_scroll ) {
                final Component parent = inner.getParent();
                if ( parent != null ) {
                    final MouseWheelEvent copy = new MouseWheelEvent( parent, e.getID(), e.getWhen(),
                            e.getModifiersEx(), e.getX(), e.getY(), e.getClickCount(), e.isPopupTrigger(),
                            e.getScrollType(), e.getScrollAmount(), e.getWheelRotation() );
                    parent.dispatchEvent( copy );
                    e.consume();
                }
            }
        } );
    }

    // ------------------------------------------------------------------ page
    /** A fresh, empty page: sections stack vertically, and it tracks the viewport width (nothing scrolls sideways). */
    static JPanel newPage() {
        final JPanel page = new ScrollablePage();
        page.setLayout( new BoxLayout( page, BoxLayout.Y_AXIS ) );
        page.setBorder( BorderFactory.createEmptyBorder( 4, 14, 10, 14 ) );
        return page;
    }

    /** The borderless scroll pane a page lives in. */
    static JScrollPane pageScroller( final JPanel page ) {
        final JScrollPane scroll = new JScrollPane( page );
        scroll.setBorder( BorderFactory.createEmptyBorder() );
        scroll.setHorizontalScrollBarPolicy( JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
        scroll.getVerticalScrollBar().setUnitIncrement( 16 );
        return scroll;
    }

    /** The scrolling page: tracks the viewport WIDTH (so nothing scrolls sideways) but not the height. */
    static final class ScrollablePage extends JPanel implements Scrollable {

        private static final long serialVersionUID = 1L;

        @Override
        public Dimension getPreferredScrollableViewportSize() {
            return getPreferredSize();
        }

        @Override
        public int getScrollableUnitIncrement( final Rectangle r, final int o, final int d ) {
            return 16;
        }

        @Override
        public int getScrollableBlockIncrement( final Rectangle r, final int o, final int d ) {
            return Math.max( 16, r.height - 16 );
        }

        @Override
        public boolean getScrollableTracksViewportWidth() {
            return true;
        }

        @Override
        public boolean getScrollableTracksViewportHeight() {
            return false;
        }
    }

    // ------------------------------------------------------------------ collapsible section
    /** A collapsible section: a bold heading (with a disclosure icon and an optional muted detail) over a body. */
    static final class Section extends JPanel {

        private static final long serialVersionUID = 1L;
        private final JPanel      _body_wrap;
        private final JButton     _header;
        private final JLabel      _detail;
        private boolean           _expanded;

        Section( final String title, final String detail, final JComponent body, final boolean expanded ) {
            super( new BorderLayout() );
            setOpaque( false );
            setAlignmentX( LEFT_ALIGNMENT );
            _expanded = expanded;
            _header = new JButton( title );
            _header.putClientProperty( "JButton.buttonType", "toolBarButton" );
            _header.setHorizontalAlignment( SwingConstants.LEFT );
            _header.setFont( _header.getFont().deriveFont( Font.BOLD ) );
            _header.setFocusable( false );
            _header.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
            _header.addActionListener( e -> toggle() );
            _detail = new JLabel( ( detail == null ) ? "" : detail );
            _detail.setForeground( mutedColor() );
            final JPanel head = new JPanel( new BorderLayout( 6, 0 ) );
            head.setOpaque( false );
            head.add( _header, BorderLayout.WEST );
            head.add( _detail, BorderLayout.CENTER );
            head.setBorder( BorderFactory.createCompoundBorder( BorderFactory
                    .createMatteBorder( 1, 0, 0, 0, borderColor() ), BorderFactory.createEmptyBorder( 4, 0, 4, 0 ) ) );
            head.addMouseListener( new MouseAdapter() {

                @Override
                public void mouseClicked( final MouseEvent e ) {
                    toggle();
                }
            } );
            head.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
            add( head, BorderLayout.NORTH );
            _body_wrap = new JPanel( new BorderLayout() );
            _body_wrap.setOpaque( false );
            _body_wrap.setBorder( BorderFactory.createEmptyBorder( 2, SECTION_INDENT, 10, 0 ) );
            _body_wrap.add( body, BorderLayout.CENTER );
            add( _body_wrap, BorderLayout.CENTER );
            _body_wrap.setVisible( expanded );
            updateIcon();
        }

        String getTitle() {
            return _header.getText().replaceFirst( "^[▾▸] ", "" );
        }

        boolean isExpanded() {
            return _expanded;
        }

        void toggle() {
            setExpanded( !_expanded );
        }

        void setExpanded( final boolean expanded ) {
            _expanded = expanded;
            _body_wrap.setVisible( _expanded );
            updateIcon();
            revalidate(); // the enclosing scroll pane is the validate root: the page re-lays out
            repaint();
        }

        void setDetail( final String detail ) {
            _detail.setText( ( detail == null ) ? "" : detail );
        }

        /** Swaps the section body (a refresh that keeps the heading and the expanded state). */
        void setBody( final JComponent body ) {
            _body_wrap.removeAll();
            _body_wrap.add( body, BorderLayout.CENTER );
            revalidate();
            repaint();
        }

        private void updateIcon() {
            final Icon icon = UIManager.getIcon( _expanded ? "Tree.expandedIcon" : "Tree.collapsedIcon" );
            _header.setIcon( icon );
            if ( icon == null ) {
                _header.setText( ( _expanded ? "▾ " : "▸ " ) + getTitle() );
            }
        }

        @Override
        public Dimension getMaximumSize() {
            final Dimension d = super.getPreferredSize();
            return new Dimension( Integer.MAX_VALUE, d.height );
        }
    }

    // ------------------------------------------------------------------ label/value grid
    /** The label/field grid used inside every section: a fixed-width label column, then the value(s). */
    static final class Grid extends JPanel {

        private static final long serialVersionUID = 1L;
        private final int         _label_width;
        private int               _row             = 0;

        /** @param label_width the label column's width (measure the longest label once, so rows line up) */
        Grid( final int label_width ) {
            super( new GridBagLayout() );
            setOpaque( false );
            _label_width = label_width;
        }

        int rowCount() {
            return _row;
        }

        void row( final String label, final JComponent field, final boolean top_aligned ) {
            row( makeLabel( label ), field, top_aligned );
        }

        void row( final JComponent label, final JComponent field, final boolean top_aligned ) {
            final GridBagConstraints c = base( top_aligned );
            c.gridx = 0;
            c.weightx = 0;
            add( label, c );
            c.gridx = 1;
            c.gridwidth = 3;
            c.weightx = 1;
            c.fill = GridBagConstraints.HORIZONTAL;
            c.insets = new Insets( 3, 0, 3, 0 );
            add( field, c );
            _row++;
        }

        void row( final String label, final JComponent field, final String label2, final JComponent field2 ) {
            final GridBagConstraints c = base( false );
            c.gridx = 0;
            c.weightx = 0;
            add( makeLabel( label ), c );
            c.gridx = 1;
            c.weightx = 0.6;
            c.fill = GridBagConstraints.HORIZONTAL;
            c.insets = new Insets( 3, 0, 3, 0 );
            add( field, c );
            c.gridx = 2;
            c.weightx = 0;
            c.fill = GridBagConstraints.NONE;
            c.insets = new Insets( 3, 14, 3, 8 );
            final JLabel l2 = new JLabel( label2 );
            add( l2, c );
            c.gridx = 3;
            c.weightx = 0.4;
            c.fill = GridBagConstraints.HORIZONTAL;
            c.insets = new Insets( 3, 0, 3, 0 );
            add( field2, c );
            _row++;
        }

        /** A full-width row without a label (a histogram, a note). */
        void span( final JComponent c_full ) {
            final GridBagConstraints c = base( false );
            c.gridx = 0;
            c.gridwidth = 4;
            c.weightx = 1;
            c.fill = GridBagConstraints.HORIZONTAL;
            c.insets = new Insets( 3, 0, 3, 0 );
            add( c_full, c );
            _row++;
        }

        private GridBagConstraints base( final boolean top_aligned ) {
            final GridBagConstraints c = new GridBagConstraints();
            c.gridy = _row;
            c.insets = new Insets( 3, 0, 3, 8 );
            c.anchor = top_aligned ? GridBagConstraints.NORTHWEST : GridBagConstraints.WEST;
            return c;
        }

        private JLabel makeLabel( final String text ) {
            final JLabel l = new JLabel( text );
            final Dimension d = l.getPreferredSize();
            l.setPreferredSize( new Dimension( Math.max( _label_width, d.width ), d.height ) );
            return l;
        }
    }
}
