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
 * Guards the Help/About changes: the About box is now HTML with clickable links to the website, the source
 * repository, and the phyloXML paper; the "phyloXML location" line is gone; the About logo is bundled in the jar
 * and loads; and the Help-menu "Documentation"/"Archaeopteryx Home" links point at the website. Headless.
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
            return aboutHtmlContent() && logoBundled() && helpUrls();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected: " + e );
        }
    }

    /** The About HTML keeps the essentials, gains clickable website/GitHub/phyloXML-paper links, and drops the
     *  "phyloXML location" line. */
    private static boolean aboutHtmlContent() {
        final String t = MainFrame.buildAboutHtml();
        if ( !t.contains( "Archaeopteryx" ) || !t.contains( AptxConstants.VERSION ) ) {
            return fail( "About should contain the program name and version" );
        }
        if ( !t.contains( "References:" ) || !t.contains( AptxConstants.PHYLOXML_REFERENCE_SHORT ) ) {
            return fail( "About should keep the References section" );
        }
        if ( !t.contains( AptxConstants.AUTHOR_EMAIL ) ) {
            return fail( "About should keep the Comments email" );
        }
        // the website, source repository, and phyloXML paper are now CLICKABLE links (href), not bare text
        if ( !t.contains( "href=\"" + AptxConstants.APTX_WEB_SITE + "\"" ) ) {
            return fail( "About should link to the website, got: " + t );
        }
        if ( !t.contains( "href=\"" + AptxConstants.APTX_GITHUB + "\"" ) ) {
            return fail( "About should link to the source repository (GitHub), got: " + t );
        }
        if ( !t.contains( "href=\"" + AptxConstants.PHYLOXML_REFERENCE_URL + "\"" ) ) {
            return fail( "About should make the phyloXML reference a clickable link to the paper, got: " + t );
        }
        if ( t.contains( "phyloXML location" ) ) {
            return fail( "About should no longer list the phyloXML location" );
        }
        if ( t.contains( "For more information" ) ) {
            return fail( "About should no longer have the 'For more information & download' section" );
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

    /** The Help-menu links point where the change specifies: "Documentation" and "Archaeopteryx Home" at the
     *  website, "Archaeopteryx.js" at its repo, and the About "Source (GitHub)" link at the forester repo. */
    private static boolean helpUrls() {
        if ( !AptxConstants.APTX_WEB_SITE.contains( "cmzmasek.github.io/archaeopteryx" ) ) {
            return fail( "'Archaeopteryx Home' should point at the website, got " + AptxConstants.APTX_WEB_SITE );
        }
        if ( !AptxConstants.APTX_DOC_SITE.contains( "cmzmasek.github.io/archaeopteryx" ) ) {
            return fail( "'Documentation' should point at the website, got " + AptxConstants.APTX_DOC_SITE );
        }
        if ( !AptxConstants.APTX_JS_WEB_SITE.contains( "github.com/cmzmasek/archaeopteryx-js" ) ) {
            return fail( "'Archaeopteryx.js' should point at its repository, got " + AptxConstants.APTX_JS_WEB_SITE );
        }
        if ( !AptxConstants.APTX_GITHUB.contains( "github.com/cmzmasek/forester" ) ) {
            return fail( "the About GitHub link should point at the forester repo, got " + AptxConstants.APTX_GITHUB );
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
