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

import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.image.BufferedImage;

/**
 * Unit test for {@link ImageSelection}, the clipboard {@link java.awt.datatransfer.Transferable} used by "Copy
 * Image to Clipboard": it must offer exactly the image flavor, hand back the SAME image for it, and throw (not
 * return null) for any other flavor. Pure logic -- runs headless, no display needed.
 */
public final class ImageSelectionTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ImageSelection: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            final BufferedImage img = new BufferedImage( 3, 2, BufferedImage.TYPE_INT_RGB );
            final ImageSelection sel = new ImageSelection( img );
            // exactly one flavor, and it is the image flavor
            final DataFlavor[] flavors = sel.getTransferDataFlavors();
            if ( ( flavors == null ) || ( flavors.length != 1 ) || !DataFlavor.imageFlavor.equals( flavors[ 0 ] ) ) {
                System.out.println( "  [ImageSelectionTest] should offer exactly the image flavor" );
                return false;
            }
            if ( !sel.isDataFlavorSupported( DataFlavor.imageFlavor ) ) {
                System.out.println( "  [ImageSelectionTest] the image flavor must be supported" );
                return false;
            }
            if ( sel.isDataFlavorSupported( DataFlavor.stringFlavor ) ) {
                System.out.println( "  [ImageSelectionTest] no flavor other than image should be supported" );
                return false;
            }
            // the image flavor returns the SAME image instance
            if ( sel.getTransferData( DataFlavor.imageFlavor ) != img ) {
                System.out.println( "  [ImageSelectionTest] the image flavor must return the wrapped image" );
                return false;
            }
            // an unsupported flavor must throw, not return null
            try {
                sel.getTransferData( DataFlavor.stringFlavor );
                System.out.println( "  [ImageSelectionTest] an unsupported flavor must throw" );
                return false;
            }
            catch ( final UnsupportedFlavorException expected ) {
                // correct
            }
            return true;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private ImageSelectionTest() {
    }
}
