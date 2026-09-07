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

import java.awt.GraphicsEnvironment;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.text.JTextComponent;

import org.forester.archaeopteryx.FormWidgets.Grid;
import org.forester.archaeopteryx.FormWidgets.Header;
import org.forester.archaeopteryx.FormWidgets.Section;

/**
 * Tests for the shared {@link FormWidgets} (needs a display): the header's texts, a section's title / toggle /
 * body swap, the grid's row counting, the read-only value's width cap and non-focusability, the edit field's
 * placeholder, the page scroller's width tracking, and the theme colours never being null.
 */
public final class FormWidgetsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FormWidgets: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> run( ok ) );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void run( final boolean[] ok ) {
        // colours and fonts
        check( ok, "colours", ( FormWidgets.mutedColor() != null ) && ( FormWidgets.borderColor() != null )
                && ( FormWidgets.accentColor() != null ) && ( FormWidgets.errorColor() != null ) );
        check( ok, "mono font", FormWidgets.monoFont( new JLabel().getFont() ).getFamily().toLowerCase()
                .contains( "mono" ) || FormWidgets.monoFont( new JLabel().getFont() ).getName().equals( "Monospaced" ) );
        // header
        final Header h = new Header( "Title", "sub" );
        check( ok, "header texts", "Title".equals( h.getTitle() ) && "sub".equals( h.getSubtitle() ) );
        h.setTitle( "T2" );
        h.setSubtitle( "s2" );
        check( ok, "header set", "T2".equals( h.getTitle() ) && "s2".equals( h.getSubtitle() ) );
        // section
        final JPanel body = new JPanel();
        final Section s = new Section( "Sec", "detail", body, false );
        check( ok, "section title", "Sec".equals( s.getTitle() ) );
        check( ok, "collapsed", !s.isExpanded() && !body.isShowing() );
        s.toggle();
        check( ok, "expanded", s.isExpanded() );
        s.setExpanded( false );
        check( ok, "set collapsed", !s.isExpanded() );
        final JPanel body2 = new JPanel();
        s.setBody( body2 );
        check( ok, "body swapped", ( body2.getParent() != null ) && ( body.getParent() == null ) );
        check( ok, "max height is preferred", s.getMaximumSize().height == s.getPreferredSize().height );
        s.setDetail( null );
        // grid
        final Grid g = new Grid( 80 );
        g.row( "a", new JTextField(), false );
        g.row( "b", new JTextField(), "c", new JTextField() );
        g.span( new JLabel( "wide" ) );
        check( ok, "grid rows", g.rowCount() == 3 );
        check( ok, "grid components: 2 + 4 + 1", g.getComponentCount() == 7 );
        check( ok, "label column width", g.getComponent( 0 ).getPreferredSize().width >= 80 );
        // read-only values
        final JTextComponent one = FormWidgets.viewValue( "x".repeat( 400 ), false );
        check( ok, "single-line value width capped", one.getPreferredSize().width <= FormWidgets.VIEW_VALUE_WIDTH );
        check( ok, "not focusable, not editable", !one.isFocusable() && !one.isEditable() );
        final JTextComponent multi = FormWidgets.viewValue( "line\nline", true );
        check( ok, "multi-line value width fixed", multi.getPreferredSize().width == FormWidgets.VIEW_VALUE_WIDTH );
        check( ok, "multi-line is an area", multi instanceof JTextArea );
        // edit field
        final JTextField ef = FormWidgets.editField( "v", "hint" );
        check( ok, "placeholder", "hint".equals( ef.getClientProperty( "JTextField.placeholderText" ) ) );
        check( ok, "no select-all on focus", "never".equals( ef.getClientProperty( "JTextField.selectAllOnFocusPolicy" ) ) );
        // page + scroller
        final JPanel page = FormWidgets.newPage();
        page.add( s );
        final JScrollPane sp = FormWidgets.pageScroller( page );
        check( ok, "tracks width", ( (javax.swing.Scrollable) page ).getScrollableTracksViewportWidth()
                && !( (javax.swing.Scrollable) page ).getScrollableTracksViewportHeight() );
        check( ok, "no horizontal bar", sp.getHorizontalScrollBarPolicy() == JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
        final JFrame f = new JFrame();
        f.getContentPane().add( sp );
        f.setSize( 300, 200 );
        f.validate();
        check( ok, "page as wide as the viewport", page.getWidth() == sp.getViewport().getWidth() );
        f.dispose();
        // buttons
        final boolean[] ran = { false };
        FormWidgets.linkButton( "+ Add", () -> ran[ 0 ] = true ).doClick();
        check( ok, "link button runs", ran[ 0 ] );
        ran[ 0 ] = false;
        FormWidgets.removeButton( "tip", () -> ran[ 0 ] = true ).doClick();
        check( ok, "remove button runs", ran[ 0 ] );
        final int[] n = { 0 };
        final JTextField t = new JTextField();
        t.getDocument().addDocumentListener( FormWidgets.onChange( () -> n[ 0 ]++ ) );
        t.setText( "a" );
        t.setText( "" );
        check( ok, "onChange fires on insert and remove", n[ 0 ] == 2 );
    }

    private static void check( final boolean[] ok, final String what, final boolean condition ) {
        if ( !condition ) {
            System.out.println( "  [FormWidgetsTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private FormWidgetsTest() {
    }
}
