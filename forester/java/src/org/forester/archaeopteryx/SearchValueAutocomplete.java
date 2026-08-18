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
import java.awt.Point;
import java.awt.Window;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
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
 * A type-to-filter value-suggestion popup for one search value field. As the user focuses / types in the value box,
 * it offers the DISTINCT EXISTING values of the selected field in the current tree (from {@code valuesSupplier}),
 * filtered by a case-honouring substring match with prefix matches ranked first; picking one fills the box (an exact
 * value) and runs the search. Cardinality-agnostic: a tiny value set (e.g. node type = leaf/internal/root) is a browse
 * list, a large one (hundreds of hosts/codes) filters as you type and is capped, so it never renders a giant menu.
 *
 * <p>Design decisions (see the search-autocomplete notes): the value set is recomputed once per popup SESSION (on
 * focus-gain / first keystroke) and filtered in memory per keystroke -- there is NO cross-session cache, so it cannot
 * go stale against a tree/data edit (the recompute-on-open IS the invalidation). The popup window is non-focusable, so
 * showing it and clicking a row never steal focus from the value field. Accepting a value leaves the match MODE
 * unchanged (picking a full value in "contains" is still correct; it is an exact hit in exact/whole-word mode).
 */
final class SearchValueAutocomplete {

    static final int    VISIBLE_ROWS = 10;   // rows shown before the list scrolls
    static final int    MATCH_CAP    = 100;  // never build/show more than this many matches (then the hint row)
    static final String HINT_TEXT    = "… keep typing to narrow"; // non-selectable last row when truncated

    private final JTextField                _field;
    private final Supplier<List<String>>    _values;         // recompute-on-open; empty => autocomplete N/A right now
    private final BooleanSupplier           _case_sensitive;
    private final Runnable                  _on_accept;

    private JWindow                         _popup;
    private DefaultListModel<String>        _model;
    private JList<String>                   _list;
    private List<String>                    _session_values; // computed once per session; null == no active session
    private boolean                         _adjusting;      // guard: our own setText must not re-trigger us
    private int                             _hint_index = -1;

    SearchValueAutocomplete( final JTextField field, final Supplier<List<String>> values,
                             final BooleanSupplier case_sensitive, final Runnable on_accept ) {
        _field = field;
        _values = values;
        _case_sensitive = case_sensitive;
        _on_accept = on_accept;
        install();
    }

    private void install() {
        _field.addFocusListener( new FocusAdapter() {

            @Override
            public void focusGained( final FocusEvent e ) {
                startSession();
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

    // ---- session / filtering -----------------------------------------------------------------------------------

    /** Begins a popup session: compute the field's distinct values once, then show the (browse) list. */
    private void startSession() {
        final List<String> vals = ( _values == null ) ? null : _values.get();
        if ( ( vals == null ) || vals.isEmpty() ) {
            _session_values = null; // nothing to suggest for this field/mode
            return;
        }
        _session_values = vals;
        renderFiltered();
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
            startSession();
        }
        else {
            renderFiltered();
        }
    }

    /** Filters the cached session values by the current text and (re)shows the popup. */
    private void renderFiltered() {
        if ( _session_values == null ) {
            return;
        }
        final boolean cs = ( _case_sensitive != null ) && _case_sensitive.getAsBoolean();
        final List<String> shown = displayModel( _session_values, _field.getText(), cs, MATCH_CAP );
        ensurePopup();
        _model.clear();
        for ( final String s : shown ) {
            _model.addElement( s );
        }
        _hint_index = ( !shown.isEmpty() && HINT_TEXT.equals( shown.get( shown.size() - 1 ) ) )
                ? shown.size() - 1 : -1;
        if ( _model.isEmpty() ) {
            hideWindow();
            return;
        }
        _list.clearSelection(); // no pre-selection -> ENTER falls through to a plain search until the user arrows
        positionAndShow();
    }

    /** The full ordered match list: prefix matches first (then other substring matches), alphabetical within each
     *  group (relies on {@code all} being sorted); empty query returns {@code all} (the browse list). Pure/testable. */
    static List<String> filterValues( final List<String> all, final String typed, final boolean case_sensitive ) {
        final String q = ( typed == null ) ? "" : typed.trim();
        if ( q.isEmpty() ) {
            return all;
        }
        final String needle = case_sensitive ? q : q.toLowerCase( Locale.ROOT );
        final List<String> prefix = new ArrayList<>();
        final List<String> other = new ArrayList<>();
        for ( final String v : all ) {
            final String hay = case_sensitive ? v : v.toLowerCase( Locale.ROOT );
            final int idx = hay.indexOf( needle );
            if ( idx == 0 ) {
                prefix.add( v );
            }
            else if ( idx > 0 ) {
                other.add( v );
            }
        }
        prefix.addAll( other );
        return prefix;
    }

    /** The list the popup actually shows: the filtered matches capped at {@code cap}, with a non-selectable
     *  {@link #HINT_TEXT} row appended when there are more matches than the cap. Pure/testable. */
    static List<String> displayModel( final List<String> all, final String typed, final boolean case_sensitive,
                                      final int cap ) {
        final List<String> matches = filterValues( all, typed, case_sensitive );
        if ( matches.size() > cap ) {
            final List<String> out = new ArrayList<>( matches.subList( 0, cap ) );
            out.add( HINT_TEXT );
            return out;
        }
        return new ArrayList<>( matches );
    }

    // ---- key / mouse handling ----------------------------------------------------------------------------------

    private void onKeyPressed( final KeyEvent e ) {
        if ( ( _popup == null ) || !_popup.isVisible() ) {
            return; // when hidden, let every key behave normally (a keystroke re-opens via the document listener)
        }
        switch ( e.getKeyCode() ) {
            case KeyEvent.VK_DOWN:
                moveSelection( 1 );
                e.consume();
                break;
            case KeyEvent.VK_UP:
                moveSelection( -1 );
                e.consume();
                break;
            case KeyEvent.VK_ESCAPE:
                hideWindow(); // keep the session so the next keystroke re-shows; caret/text untouched
                e.consume();
                break;
            case KeyEvent.VK_ENTER: {
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
            case KeyEvent.VK_TAB: {
                final String v = selectedValue();
                if ( v != null ) {
                    accept( v, true ); // Tab is not consumed -> focus moves on after the value is filled
                }
                break;
            }
            default:
                break;
        }
    }

    private void moveSelection( final int dir ) {
        final int n = _model.size();
        final int max_sel = ( _hint_index >= 0 ) ? _hint_index - 1 : n - 1; // never land on the hint row
        if ( max_sel < 0 ) {
            return;
        }
        int sel = _list.getSelectedIndex();
        if ( dir > 0 ) {
            sel = ( sel < 0 ) ? 0 : sel + 1;
            if ( sel > max_sel ) {
                sel = 0;
            }
        }
        else {
            sel = ( sel < 0 ) ? max_sel : sel - 1;
            if ( sel < 0 ) {
                sel = max_sel;
            }
        }
        _list.setSelectedIndex( sel );
        _list.ensureIndexIsVisible( sel );
    }

    private String selectedValue() {
        final int i = _list.getSelectedIndex();
        if ( ( i < 0 ) || ( i == _hint_index ) ) {
            return null;
        }
        return _model.getElementAt( i );
    }

    /** Fills the value box with {@code value} (an exact existing value), closes the popup, and -- when
     *  {@code fire_search} -- runs the search (used by the click / Tab paths; the ENTER path lets the box's own
     *  keyReleased fire the search instead, so it is not double-run). */
    private void accept( final String value, final boolean fire_search ) {
        _adjusting = true;
        try {
            _field.setText( value );
            _field.setCaretPosition( value.length() );
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
        _model = new DefaultListModel<>();
        _list = new JList<>( _model );
        _list.setFont( _field.getFont() );
        _list.setSelectionMode( ListSelectionModel.SINGLE_SELECTION );
        _list.setCellRenderer( new HintRenderer() );
        _list.addMouseListener( new MouseAdapter() {

            @Override
            public void mousePressed( final MouseEvent e ) {
                final int i = _list.locationToIndex( e.getPoint() );
                if ( ( i < 0 ) || ( i == _hint_index ) || !_list.getCellBounds( i, i ).contains( e.getPoint() ) ) {
                    return;
                }
                accept( _model.getElementAt( i ), true );
            }
        } );
        _list.addMouseMotionListener( new java.awt.event.MouseMotionAdapter() {

            @Override
            public void mouseMoved( final MouseEvent e ) {
                final int i = _list.locationToIndex( e.getPoint() );
                if ( ( i >= 0 ) && ( i != _hint_index ) && _list.getCellBounds( i, i ).contains( e.getPoint() ) ) {
                    _list.setSelectedIndex( i );
                }
            }
        } );
        final JScrollPane sp = new JScrollPane( _list, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
                                                ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER );
        sp.setBorder( BorderFactory.createLineBorder( Color.GRAY ) );
        _popup.getContentPane().add( sp );
    }

    private void positionAndShow() {
        if ( !_field.isShowing() ) {
            return; // headless / not yet realized -- the in-memory model is up to date, just don't show a window
        }
        _list.setVisibleRowCount( Math.min( VISIBLE_ROWS, _model.size() ) );
        _popup.pack();
        final int w = Math.max( _field.getWidth(), _popup.getWidth() );
        _popup.setSize( w, _popup.getHeight() );
        final Point p = _field.getLocationOnScreen();
        _popup.setLocation( p.x, p.y + _field.getHeight() );
        _popup.setVisible( true );
    }

    private void hideWindow() {
        if ( _popup != null ) {
            _popup.setVisible( false );
        }
    }

    /** Greys + italicises the non-selectable {@link #HINT_TEXT} row; every other row renders normally. */
    private final class HintRenderer extends DefaultListCellRenderer {

        @Override
        public Component getListCellRendererComponent( final JList<?> list, final Object value, final int index,
                                                       final boolean selected, final boolean focus ) {
            final boolean hint = ( index == _hint_index );
            final Component c = super.getListCellRendererComponent( list, value, index, selected && !hint, false );
            if ( hint ) {
                c.setForeground( Color.GRAY );
                c.setFont( c.getFont().deriveFont( Font.ITALIC ) );
            }
            return c;
        }
    }

    // ---- test hooks (headless-drivable; no window required) -----------------------------------------------------

    /** Recomputes the session values and returns exactly what the popup WOULD show for the current field text
     *  (filtered, capped, hint row included). Does not open a window, so it runs headless. */
    List<String> displayModelForTest() {
        final List<String> vals = ( _values == null ) ? null : _values.get();
        _session_values = ( ( vals == null ) || vals.isEmpty() ) ? null : vals;
        if ( _session_values == null ) {
            return new ArrayList<>();
        }
        final boolean cs = ( _case_sensitive != null ) && _case_sensitive.getAsBoolean();
        return displayModel( _session_values, _field.getText(), cs, MATCH_CAP );
    }

    /** Simulates picking {@code value} from the popup (fills the field + runs the search). */
    void acceptValueForTest( final String value ) {
        accept( value, true );
    }

    /** Whether the suggestion window is currently visible (exercises the realized show path). */
    boolean isPopupShowingForTest() {
        return ( _popup != null ) && _popup.isVisible();
    }

    /** The rows currently shown in the popup (hint row included), for a realized-frame check. */
    List<String> currentModelForTest() {
        final List<String> out = new ArrayList<>();
        if ( _model != null ) {
            for ( int i = 0; i < _model.size(); i++ ) {
                out.add( _model.getElementAt( i ) );
            }
        }
        return out;
    }
}
