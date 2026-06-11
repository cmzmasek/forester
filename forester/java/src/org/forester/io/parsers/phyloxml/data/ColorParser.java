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

package org.forester.io.parsers.phyloxml.data;

import java.awt.Color;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.PhylogenyData;

public class ColorParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new ColorParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private ColorParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhyloXmlDataFormatException {
        int red = 0;
        int green = 0;
        int blue = 0;
        int alpha = -1;
        for( int j = 0; j < element.getNumberOfChildElements(); ++j ) {
            final XmlElement c = element.getChildElement( j );
            if ( c.getQualifiedName().equals( PhyloXmlMapping.COLOR_RED ) ) {
                red = c.getValueAsInt();
            }
            else if ( c.getQualifiedName().equals( PhyloXmlMapping.COLOR_GREEN ) ) {
                green = c.getValueAsInt();
            }
            else if ( c.getQualifiedName().equals( PhyloXmlMapping.COLOR_BLUE ) ) {
                blue = c.getValueAsInt();
            }
            else if ( c.getQualifiedName().equals( PhyloXmlMapping.COLOR_ALPHA ) ) {
                alpha = c.getValueAsInt();
            }
        }
        final BranchColor color = new BranchColor();
        if ( alpha < 0 ) { 
            color.setValue( new Color( red, green, blue ) );
        }
        else { 
            color.setValue( new Color( red, green, blue, alpha ) );
        }
        return color;
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
