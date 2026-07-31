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

import java.awt.Image;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;

/**
 * A {@link Transferable} that carries a single {@link Image} for the system clipboard, offered as
 * {@link DataFlavor#imageFlavor} only. Used by the "Copy Image to Clipboard" action so a rendered tree
 * can be pasted into a document or slide. (The JDK ships no public image Transferable; this is the
 * standard tiny implementation.)
 */
final class ImageSelection implements Transferable {

    private final Image _image;

    ImageSelection( final Image image ) {
        _image = image;
    }

    @Override
    public DataFlavor[] getTransferDataFlavors() {
        return new DataFlavor[] { DataFlavor.imageFlavor };
    }

    @Override
    public boolean isDataFlavorSupported( final DataFlavor flavor ) {
        return DataFlavor.imageFlavor.equals( flavor );
    }

    @Override
    public Object getTransferData( final DataFlavor flavor ) throws UnsupportedFlavorException {
        if ( !DataFlavor.imageFlavor.equals( flavor ) ) {
            throw new UnsupportedFlavorException( flavor );
        }
        return _image;
    }
}
