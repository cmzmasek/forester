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

import java.awt.Component;
import java.awt.GraphicsEnvironment;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Guards that the Tools "Collapse Branches" submenu uses the same font as the regular menu items -- the
 * look-and-feel's menu-item font, the one GUI font since 2026-10-27 ({@code MainFrame.applyGuiFont}); a font
 * set on the submenu by hand would make it look more prominent than its siblings. Needs FlatLaf + a display, so it is a no-op when headless and runs
 * standalone (the established pattern for the GUI tests here).
 */
public final class CollapseMenuFontTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CollapseMenuFont: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tinyTree() }, conf, "cmf" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final JMenu collapse = findSubmenu( mf[ 0 ]._tools_menu, "Collapse Branches" );
                    if ( collapse == null ) {
                        ok[ 0 ] = false;
                        System.out.println( "  'Collapse Branches' submenu not found in Tools" );
                    }
                    else if ( !javax.swing.UIManager.getFont( "MenuItem.font" ).equals( collapse.getFont() ) ) {
                        ok[ 0 ] = false;
                        System.out.println( "  submenu font " + collapse.getFont() + " != menu-item font "
                                + javax.swing.UIManager.getFont( "MenuItem.font" ) );
                    }
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    private static JMenu findSubmenu( final JMenu menu, final String text ) {
        if ( menu != null ) {
            for( final Component c : menu.getMenuComponents() ) {
                if ( ( c instanceof JMenu ) && text.equals( ( (JMenu) c ).getText() ) ) {
                    return (JMenu) c;
                }
            }
        }
        return null;
    }

    private static Phylogeny tinyTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( new PhylogenyNode() );
        root.addAsChild( new PhylogenyNode() );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private CollapseMenuFontTest() {
    }
}
