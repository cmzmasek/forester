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

package org.forester.applications;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.archaeopteryx.Configuration;
import org.forester.archaeopteryx.Options;
import org.forester.archaeopteryx.TreeColorSet;

public class phylo2graphics {

    public static void main( final String[] args ) {
        try {
            final Configuration config = new Configuration();
            // Could also read a configuration file with:
            // Configuration config = new Configuration("my_configuration_file.txt", false, false, false);
            config.putDisplayColors( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
            config.putDisplayColors( TreeColorSet.BRANCH, new Color( 0, 0, 0 ) );
            config.putDisplayColors( TreeColorSet.TAXONOMY, new Color( 0, 0, 0 ) );
            config.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            AptxUtil.writePhylogenyToGraphicsFile( new File( "/home/czmasek/tol_117_TEST.xml" ),
                                                   new File( "/home/czmasek/tol_117_TEST_.png" ),
                                                   1000,
                                                   1000,
                                                   GraphicsExportType.PNG,
                                                   config );
            // If the tree 'phy' already exists, can also use this:
            //AptxUtil.writePhylogenyToGraphicsFile( phy,
            //                                       new File( "/home/czmasek/tol_117_TEST_.png" ),
            //                                      1000,
            //                                      1000,
            //                                      GraphicsExportType.PNG,
            //                                      config );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }
}
