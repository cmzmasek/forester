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

import java.awt.image.BufferedImage;
import java.io.InputStream;

import javax.imageio.ImageIO;
import javax.swing.Icon;

/**
 * Guards the Help/About changes: the About text keeps the program-info + reference content but no longer has the
 * "For more information & download" or "Documentation" sections (nor their URLs); the About logo is bundled in the
 * jar and loads; and the Help-menu links point at the GitHub README / the archaeopteryx-js repo. Headless.
 */
public final class AboutDialogTest {

    private static final String LOGO    = "/resources/images/archaeopteryx-logo.png";
    private static final String LOGO_2X = "/resources/images/archaeopteryx-logo@2x.png";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "AboutDialog: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return aboutTextTrimmed() && logoBundled() && helpUrls();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected: " + e );
        }
    }

    /** The About text still has the essentials, but the two removed sections (and their URLs) are gone. */
    private static boolean aboutTextTrimmed() {
        final String t = MainFrame.buildAboutText();
        if ( !t.contains( "Archaeopteryx" ) || !t.contains( AptxConstants.VERSION ) ) {
            return fail( "About text should contain the program name and version" );
        }
        if ( !t.contains( "References:" ) || !t.contains( AptxConstants.PHYLOXML_REFERENCE_SHORT ) ) {
            return fail( "About text should keep the References section" );
        }
        if ( !t.contains( AptxConstants.AUTHOR_EMAIL ) ) {
            return fail( "About text should keep the Comments email" );
        }
        if ( t.contains( "For more information" ) ) {
            return fail( "About text should no longer have the 'For more information & download' section" );
        }
        if ( t.contains( "Documentation:" ) ) {
            return fail( "About text should no longer have the 'Documentation' section" );
        }
        if ( t.contains( AptxConstants.APTX_WEB_SITE ) || t.contains( AptxConstants.APTX_DOC_SITE ) ) {
            return fail( "About text should no longer list the website / documentation URLs (they live in the Help menu)" );
        }
        return true;
    }

    /** Both the base and @2x logos are on the classpath (build copied them into the jar) and load as real images,
     *  with the @2x variant exactly twice the base (so aboutLogo's HiDPI image is a valid pair). */
    private static boolean logoBundled() throws Exception {
        final Icon icon = MainFrame.aboutLogo();
        if ( ( icon == null ) || ( icon.getIconWidth() <= 0 ) || ( icon.getIconHeight() <= 0 ) ) {
            return fail( "aboutLogo() must load the bundled logo, got " + icon );
        }
        final BufferedImage base = readBundled( LOGO );
        final BufferedImage x2 = readBundled( LOGO_2X );
        if ( base == null ) {
            return fail( "base logo not bundled at " + LOGO + " (check build.xml copy_resources)" );
        }
        if ( x2 == null ) {
            return fail( "@2x logo not bundled at " + LOGO_2X );
        }
        if ( ( base.getWidth() < 32 ) || ( x2.getWidth() != ( 2 * base.getWidth() ) )
                || ( x2.getHeight() != ( 2 * base.getHeight() ) ) ) {
            return fail( "@2x logo must be exactly twice the base (" + base.getWidth() + "x" + base.getHeight()
                    + " vs " + x2.getWidth() + "x" + x2.getHeight() + ")" );
        }
        return true;
    }

    private static BufferedImage readBundled( final String resource ) throws Exception {
        try ( final InputStream in = MainFrame.class.getResourceAsStream( resource ) ) {
            return ( in == null ) ? null : ImageIO.read( in );
        }
    }

    /** The three Help-menu links point where the change specifies. */
    private static boolean helpUrls() {
        if ( !AptxConstants.APTX_WEB_SITE.contains( "github.com/cmzmasek/forester" ) ) {
            return fail( "'Archaeopteryx Home' should point at the GitHub README, got " + AptxConstants.APTX_WEB_SITE );
        }
        if ( !AptxConstants.APTX_DOC_SITE.contains( "github.com/cmzmasek/forester" ) ) {
            return fail( "'Documentation' should point at the GitHub README, got " + AptxConstants.APTX_DOC_SITE );
        }
        if ( !AptxConstants.APTX_JS_WEB_SITE.contains( "github.com/cmzmasek/archaeopteryx-js" ) ) {
            return fail( "'Archaeopteryx.js' should point at its repository, got " + AptxConstants.APTX_JS_WEB_SITE );
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [AboutDialogTest] " + msg );
        return false;
    }

    private AboutDialogTest() {
    }
}
