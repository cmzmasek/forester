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

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.GraphicsConfiguration;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.Window;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.function.BooleanSupplier;
import java.util.function.Supplier;

import javax.swing.BorderFactory;
import javax.swing.DefaultListCellRenderer;
import javax.swing.DefaultListModel;
import javax.swing.JList;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JWindow;
import javax.swing.ListSelectionModel;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;

/**
 * The search value box's type-ahead: the distinct values of the chosen field, filtered as you type, in a list under
 * the box. A PORT of the Archaeopteryx.js suggestion list (archaeopteryx.js {@code openSuggestions} and friends,
 * commit b7e1bc8), so the two viewers suggest the same values the same way -- Christian, 2026-09-13: search is one
 * of the programs' strong points, and JS and desktop must behave exactly alike.
 * <ul>
 * <li>Only the TERM being typed is completed: {@code ','} is OR and {@code '+'} is AND in the plain-text modes, so the
 * term is whatever follows the last of those, and a pick replaces just that term.</li>
 * <li>Values are matched the way the box's MODE will match them -- a prefix for "starts with", a suffix for "ends
 * with", a substring otherwise -- honouring the box's "Match case" checkbox (Christian, 2026-09-13: a JOINT rule, both
 * viewers respect it).</li>
 * <li>At most {@link #MAX_ROWS} rows, in the values' sorted order, the matched part in bold accent; when there are more,
 * a final "N more — keep typing" row.</li>
 * <li>Nothing is shown when nothing matches, or when the one match is exactly what is typed.</li>
 * <li>Arrow keys open the list and move along it (wrapping), Enter picks, Escape closes, Tab leaves the box (which
 * closes the list without picking).</li>
 * </ul>
 * The value set is recomputed once per popup SESSION (on focus-gain / first keystroke) and filtered in memory per
 * keystroke, so it cannot go stale against a tree edit. The popup window is non-focusable, so showing it and clicking a
 * row never steal focus from the value field. Off, as before, for numeric fields, Any Text, the molecular sequence and
 * regex mode (the value supplier is empty then).
 */
final class SearchValueAutocomplete {

    /** Rows shown before the "N more" row -- the JS {@code SUGGEST_MAX_ROWS}. */
    static final int MAX_ROWS = 10;
    /** The widest the list grows (px) -- the JS {@code max-width:320px}; it is never narrower than the box. */
    static final int MAX_WIDTH = 320;

    /** What the list shows: up to {@link #MAX_ROWS} matching values, and how many more matched. */
    static final class Model {

        final List<String> _rows;
        final int          _more;

        Model( final List<String> rows, final int more ) {
            _rows = rows;
            _more = more;
        }

        boolean isEmpty() {
            return _rows.isEmpty();
        }
    }

    private final JTextField           _field;
    private final Supplier<List<String>> _values; // recompute-on-open; empty => suggestions are off right now
    private final Supplier<SearchMode> _mode;     // how the box will match, so the list filters the same way
    private final BooleanSupplier      _case_sensitive; // the "Match case" checkbox: suggestions honour it too
    private final Runnable             _on_accept;

    private JWindow                    _popup;
    private DefaultListModel<String>   _list_model;
    private JList<String>              _list;
    private List<String>               _session_values; // computed once per session; null == no active session
    private Model                      _model = new Model( Collections.<String>emptyList(), 0 );
    private String                     _term  = "";     // the term the shown rows were matched against
    private boolean                    _adjusting;      // guard: our own setText must not re-trigger us

    SearchValueAutocomplete( final JTextField field, final Supplier<List<String>> values,
                             final Supplier<SearchMode> mode, final BooleanSupplier case_sensitive,
                             final Runnable on_accept ) {
        _field = field;
        _values = values;
        _mode = mode;
        _case_sensitive = case_sensitive;
        _on_accept = on_accept;
        install();
    }

    private void install() {
        _field.addFocusListener( new FocusAdapter() {

            @Override
            public void focusGained( final FocusEvent e ) {
                // the JS list opens on focus only when the box already holds text; an empty box waits for a keystroke
                startSession( !_field.getText().isEmpty() );
            }

            @Override
            public void focusLost( final FocusEvent e ) {
                endSession();
            }
        } );
        _field.getDocument().addDocumentListener( new javax.swing.event.DocumentListener() {

            @Override
            public void insertUpdate( final javax.swing.event.DocumentEvent e ) {
                onDocChange();
            }

            @Override
            public void removeUpdate( final javax.swing.event.DocumentEvent e ) {
                onDocChange();
            }

            @Override
            public void changedUpdate( final javax.swing.event.DocumentEvent e ) {
                onDocChange();
            }
        } );
        _field.addKeyListener( new KeyAdapter() {

            @Override
            public void keyPressed( final KeyEvent e ) {
                onKeyPressed( e );
            }
        } );
    }

    // ---- the pure rules (JS parity; each is pinned by SearchAutocompleteTest) ---------------------------------

    /** Where the term being typed starts: just after the last {@code ','} (OR) or {@code '+'} (AND), else 0. */
    static int termStart( final String text ) {
        if ( text == null ) {
            return 0;
        }
        return Math.max( text.lastIndexOf( ',' ), text.lastIndexOf( '+' ) ) + 1;
    }

    /** The term being typed: the text after {@link #termStart}, trimmed. */
    static String currentTerm( final String text ) {
        return ( text == null ) ? "" : text.substring( termStart( text ) ).trim();
    }

    /** {@code s} as compared: itself when matching case, else lower-cased. */
    private static String fold( final String s, final boolean case_sensitive ) {
        return case_sensitive ? s : s.toLowerCase( Locale.ROOT );
    }

    /** Whether {@code value} admits {@code term} the way {@code mode} will match it: a prefix for STARTS_WITH, a suffix
     *  for ENDS_WITH, a substring otherwise -- ignoring case unless {@code case_sensitive}. An empty term admits
     *  everything. */
    static boolean matches( final String value, final String term, final SearchMode mode, final boolean case_sensitive ) {
        if ( term.isEmpty() ) {
            return true;
        }
        final String lv = fold( value, case_sensitive );
        final String needle = fold( term, case_sensitive );
        if ( mode == SearchMode.STARTS_WITH ) {
            return lv.startsWith( needle );
        }
        if ( mode == SearchMode.ENDS_WITH ) {
            return lv.endsWith( needle );
        }
        return lv.contains( needle );
    }

    /** The values admitting {@code term}, in the order given (the values' sorted order -- no re-ranking). */
    static List<String> matches( final List<String> values, final String term, final SearchMode mode,
                                 final boolean case_sensitive ) {
        final List<String> out = new ArrayList<>();
        for ( final String v : values ) {
            if ( matches( v, term, mode, case_sensitive ) ) {
                out.add( v );
            }
        }
        return out;
    }

    /**
     * What the list shows for the box text: the first {@link #MAX_ROWS} matches of the term being typed, and the count
     * of the rest. Empty -- nothing to offer -- when nothing matches, or when the one match is exactly what is typed
     * (compared the way the box matches: case ignored unless {@code case_sensitive}).
     */
    static Model model( final List<String> values, final String text, final SearchMode mode,
                        final boolean case_sensitive ) {
        final String term = currentTerm( text );
        final List<String> all = matches( values, term, mode, case_sensitive );
        if ( all.isEmpty()
                || ( ( all.size() == 1 ) && fold( all.get( 0 ), case_sensitive ).equals( fold( term, case_sensitive ) ) ) ) {
            return new Model( Collections.<String>emptyList(), 0 );
        }
        final List<String> rows = new ArrayList<>( all.subList( 0, Math.min( MAX_ROWS, all.size() ) ) );
        return new Model( rows, all.size() - rows.size() );
    }

    /** The trailing row when more matched than the list shows: {@code "N more — keep typing"}. */
    static String moreText( final int more ) {
        return more + " more — keep typing";
    }

    /** The part of {@code value} the typed term matched -- its FIRST occurrence, compared the way the box matches --
     *  as {@code {at, length}}, or {@code null} when the term is empty or does not occur. (The JS list bolds the
     *  first occurrence whatever the mode, even for "ends with".) */
    static int[] matchSpan( final String value, final String term, final boolean case_sensitive ) {
        if ( term.isEmpty() ) {
            return null;
        }
        final int at = fold( value, case_sensitive ).indexOf( fold( term, case_sensitive ) );
        return ( at < 0 ) ? null : new int[] { at, term.length() };
    }

    /** The box text after picking {@code value}: everything up to and including the last separator, then the value --
     *  so {@code "cat,do"} picking {@code "dog"} gives {@code "cat,dog"} and the OR stays intact. */
    static String pickInto( final String text, final String value ) {
        final String t = ( text == null ) ? "" : text;
        return t.substring( 0, termStart( t ) ) + value;
    }

    /** A row as HTML: the value, escaped, with the matched part in bold {@code accent} (a CSS colour). */
    static String rowHtml( final String value, final String term, final String accent, final boolean case_sensitive ) {
        final int[] span = matchSpan( value, term, case_sensitive );
        if ( span == null ) {
            return "<html>" + escape( value );
        }
        return "<html>" + escape( value.substring( 0, span[ 0 ] ) ) + "<b><font color=\"" + accent + "\">"
                + escape( value.substring( span[ 0 ], span[ 0 ] + span[ 1 ] ) ) + "</font></b>"
                + escape( value.substring( span[ 0 ] + span[ 1 ] ) );
    }

    static String escape( final String s ) {
        return s.replace( "&", "&amp;" ).replace( "<", "&lt;" ).replace( ">", "&gt;" ).replace( "\"", "&quot;" );
    }

    static String hex( final Color c ) {
        return String.format( "#%02x%02x%02x", c.getRed(), c.getGreen(), c.getBlue() );
    }

    // ---- session / filtering -----------------------------------------------------------------------------------

    /** Begins a popup session: compute the field's distinct values once; show the list when {@code show}. */
    private void startSession( final boolean show ) {
        final List<String> vals = ( _values == null ) ? null : _values.get();
        if ( ( vals == null ) || vals.isEmpty() ) {
            _session_values = null; // nothing to suggest for this field/mode
            return;
        }
        _session_values = vals;
        if ( show ) {
            renderFiltered();
        }
    }

    /** Ends the session (focus left the field, or the field/mode changed): drop the cached values and hide. */
    void endSession() {
        _session_values = null;
        hideWindow();
    }

    private void onDocChange() {
        if ( _adjusting || !_field.isFocusOwner() ) {
            return; // ignore our own setText and programmatic edits while unfocused (e.g. Reset)
        }
        if ( _session_values == null ) {
            startSession( true );
        }
        else {
            renderFiltered();
        }
    }

    private SearchMode currentMode() {
        final SearchMode m = ( _mode == null ) ? null : _mode.get();
        return ( m == null ) ? SearchMode.CONTAINS : m;
    }

    private boolean caseSensitive() {
        return ( _case_sensitive != null ) && _case_sensitive.getAsBoolean();
    }

    /** Filters the cached session values by the term being typed and (re)shows the popup. */
    private void renderFiltered() {
        if ( _session_values == null ) {
            return;
        }
        _term = currentTerm( _field.getText() );
        _model = model( _session_values, _field.getText(), currentMode(), caseSensitive() );
        if ( _model.isEmpty() ) {
            hideWindow();
            return;
        }
        ensurePopup();
        _list_model.clear();
        for ( final String s : _model._rows ) {
            _list_model.addElement( s );
        }
        if ( _model._more > 0 ) {
            _list_model.addElement( moreText( _model._more ) );
        }
        _list.clearSelection(); // no pre-selection -> ENTER falls through to a plain search until the user arrows
        positionAndShow();
    }

    private boolean isShowing() {
        return ( _popup != null ) && _popup.isVisible();
    }

    /** The index of the non-selectable "N more" row, or -1. */
    private int moreIndex() {
        return ( _model._more > 0 ) ? _model._rows.size() : -1;
    }

    // ---- key / mouse handling ----------------------------------------------------------------------------------

    private void onKeyPressed( final KeyEvent e ) {
        switch ( e.getKeyCode() ) {
            case KeyEvent.VK_DOWN:
            case KeyEvent.VK_UP:
                // the arrows open the list when it is closed, then move along it
                if ( !isShowing() ) {
                    if ( _session_values == null ) {
                        startSession( false );
                    }
                    renderFiltered();
                    if ( !isShowing() ) {
                        return;
                    }
                }
                moveSelection( ( e.getKeyCode() == KeyEvent.VK_DOWN ) ? 1 : -1 );
                e.consume();
                break;
            case KeyEvent.VK_ESCAPE:
                if ( isShowing() ) {
                    hideWindow(); // closes the list only (consumed, so the view does not reset); the session stays
                    e.consume();
                }
                break;
            case KeyEvent.VK_ENTER: {
                if ( !isShowing() ) {
                    return;
                }
                final String v = selectedValue();
                if ( v != null ) {
                    accept( v, false ); // the value box's own keyReleased listener runs the search on ENTER
                    e.consume();
                }
                else {
                    hideWindow(); // nothing highlighted: let ENTER run a plain search, just close the list
                }
                break;
            }
            // Tab is a focus-traversal key: Swing moves the focus before any key listener sees it, and the list closes
            // with the focus (focusLost -> endSession) -- it never picks (JS parity)
            default:
                break;
        }
    }

    private void moveSelection( final int dir ) {
        final int n = _model._rows.size(); // the "N more" row is never selectable
        if ( n == 0 ) {
            return;
        }
        final int sel = _list.getSelectedIndex();
        final int next;
        if ( dir > 0 ) {
            next = ( sel < 0 ) ? 0 : ( sel + 1 ) % n;
        }
        else {
            next = ( sel < 0 ) ? n - 1 : ( sel - 1 + n ) % n;
        }
        _list.setSelectedIndex( next );
        _list.ensureIndexIsVisible( next );
    }

    private String selectedValue() {
        final int i = _list.getSelectedIndex();
        if ( ( i < 0 ) || ( i >= _model._rows.size() ) ) {
            return null;
        }
        return _model._rows.get( i );
    }

    /** Replaces the term being typed with {@code value} (an exact existing value), closes the popup, and -- when
     *  {@code fire_search} -- runs the search (the click path; the ENTER path lets the box's own keyReleased fire the
     *  search instead, so it is not double-run). */
    private void accept( final String value, final boolean fire_search ) {
        final String text = pickInto( _field.getText(), value );
        _adjusting = true;
        try {
            _field.setText( text );
            _field.setCaretPosition( text.length() );
        }
        finally {
            _adjusting = false;
        }
        hideWindow();
        if ( fire_search && ( _on_accept != null ) ) {
            _on_accept.run();
        }
    }

    // ---- popup window ------------------------------------------------------------------------------------------

    private void ensurePopup() {
        if ( _popup != null ) {
            return;
        }
        final Window owner = SwingUtilities.getWindowAncestor( _field );
        _popup = new JWindow( owner );
        _popup.setFocusableWindowState( false ); // never steal focus from the value field
        _list_model = new DefaultListModel<>();
        _list = new JList<>( _list_model );
        _list.setFont( _field.getFont() );
        _list.setSelectionMode( ListSelectionModel.SINGLE_SELECTION );
        // the highlighted row is accent with white ink, like the JS list's .aptx-active row. The popup window is
        // non-focusable, so the look-and-feel would paint its INACTIVE (grey) selection; FlatLaf paints the selection
        // itself, over the renderer, so it is told the colours for both states (its per-component style property)
        final String accent = hex( FormWidgets.accentColor() );
        _list.putClientProperty( "FlatLaf.style", "selectionBackground: " + accent + "; selectionForeground: #ffffff; "
                + "selectionInactiveBackground: " + accent + "; selectionInactiveForeground: #ffffff" );
        _list.setCellRenderer( new RowRenderer() );
        _list.addMouseListener( new MouseAdapter() {

            @Override
            public void mousePressed( final MouseEvent e ) {
                final int i = _list.locationToIndex( e.getPoint() );
                if ( ( i < 0 ) || ( i >= _model._rows.size() ) || !_list.getCellBounds( i, i ).contains( e.getPoint() ) ) {
                    return;
                }
                accept( _model._rows.get( i ), true );
            }
        } );
        _list.addMouseMotionListener( new java.awt.event.MouseMotionAdapter() {

            @Override
            public void mouseMoved( final MouseEvent e ) {
                final int i = _list.locationToIndex( e.getPoint() );
                if ( ( i >= 0 ) && ( i < _model._rows.size() ) && _list.getCellBounds( i, i ).contains( e.getPoint() ) ) {
                    _list.setSelectedIndex( i );
                }
            }
        } );
        final JScrollPane sp = new JScrollPane( _list, ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER,
                                                ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER );
        sp.setBorder( BorderFactory.createLineBorder( Color.GRAY ) );
        _popup.getContentPane().add( sp );
    }

    /** Under the box and at least as wide as it (at most {@link #MAX_WIDTH}); above it when the screen's bottom is
     *  too close. */
    private void positionAndShow() {
        if ( !_field.isShowing() ) {
            return; // headless / not yet realized -- the in-memory model is up to date, just don't show a window
        }
        _list.setVisibleRowCount( _list_model.size() );
        _popup.pack();
        final int w = Math.max( _field.getWidth(), Math.min( _popup.getWidth(), MAX_WIDTH ) );
        final int h = _popup.getHeight();
        _popup.setSize( w, h );
        final Point p = _field.getLocationOnScreen();
        int y = p.y + _field.getHeight() + 3;
        final GraphicsConfiguration gc = _field.getGraphicsConfiguration();
        if ( gc != null ) {
            final Rectangle screen = gc.getBounds();
            final Insets in = Toolkit.getDefaultToolkit().getScreenInsets( gc );
            final int bottom = ( screen.y + screen.height ) - in.bottom - 8;
            final int top = screen.y + in.top + 8;
            if ( ( ( y + h ) > bottom ) && ( ( p.y - h - 3 ) > top ) ) {
                y = p.y - h - 3;
            }
        }
        _popup.setLocation( p.x, y );
        _popup.setVisible( true );
    }

    private void hideWindow() {
        if ( _popup != null ) {
            _popup.setVisible( false );
        }
    }

    /** A value row with the matched part in bold accent (white on the highlighted row, like the JS); the trailing
     *  "N more" row grey and small, never highlighted. */
    private final class RowRenderer extends DefaultListCellRenderer {

        @Override
        public Component getListCellRendererComponent( final JList<?> list, final Object value, final int index,
                                                       final boolean selected, final boolean focus ) {
            final boolean more = ( index == moreIndex() );
            final Component c = super.getListCellRendererComponent( list, value, index, selected && !more, false );
            if ( more ) {
                c.setForeground( Color.GRAY );
                c.setFont( c.getFont().deriveFont( Font.ITALIC, Math.max( 9f, c.getFont().getSize2D() - 1 ) ) );
                return c;
            }
            // the highlighted row is accent-on-white-ink like the JS list (.aptx-active) -- set here explicitly, because
            // the popup window never has focus and the look-and-feel would otherwise paint its INACTIVE grey selection
            final Color accent = FormWidgets.accentColor();
            if ( selected ) {
                c.setBackground( accent );
                c.setForeground( Color.WHITE );
            }
            setText( rowHtml( String.valueOf( value ), _term, hex( selected ? Color.WHITE : accent ), caseSensitive() ) );
            return c;
        }
    }

    // ---- test hooks (headless-drivable; no window required) -----------------------------------------------------

    /** Recomputes the session values and returns exactly what the popup WOULD show for the current field text.
     *  Does not open a window, so it runs headless. */
    Model modelForTest() {
        final List<String> vals = ( _values == null ) ? null : _values.get();
        _session_values = ( ( vals == null ) || vals.isEmpty() ) ? null : vals;
        if ( _session_values == null ) {
            return new Model( Collections.<String>emptyList(), 0 );
        }
        return model( _session_values, _field.getText(), currentMode(), caseSensitive() );
    }

    /** Simulates picking {@code value} from the popup (replaces the current term + runs the search). */
    void acceptValueForTest( final String value ) {
        accept( value, true );
    }

    /** Whether the list window is showing (needs a realized field). */
    boolean isShowingForTest() {
        return isShowing();
    }

    /** The highlighted row, or -1. */
    int selectedIndexForTest() {
        return ( _list == null ) ? -1 : _list.getSelectedIndex();
    }

    /** The rows the open list holds, the "N more" row included. */
    List<String> shownRowsForTest() {
        final List<String> out = new ArrayList<>();
        if ( _list_model != null ) {
            for ( int i = 0; i < _list_model.size(); ++i ) {
                out.add( _list_model.get( i ) );
            }
        }
        return out;
    }
}
