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

package org.forester.phylogeny.data;

import java.awt.Color;
import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class BranchColor implements PhylogenyData {

    private Color _color;

    public BranchColor() {
        _color = null;
    }

    public BranchColor( final Color color ) {
        _color = color;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue().toString() );
    }

    @Override
    public StringBuffer asText() {
        return new StringBuffer( getValue().toString() );
    }

    @Override
    /**
     * Not a deep copy.
     *
     */
    public PhylogenyData copy() {
        final BranchColor bc = new BranchColor();
        bc.setValue( getValue() );
        return bc;
    }

    public Color getValue() {
        return _color;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        return getValue().equals( ( ( BranchColor ) data ).getValue() );
    }

    public void setValue( final Color color ) {
        _color = color;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.COLOR );
        PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.COLOR_RED, getValue().getRed() + "", indentation );
        PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.COLOR_GREEN, getValue().getGreen() + "", indentation );
        PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.COLOR_BLUE, getValue().getBlue() + "", indentation );
        if ( getValue().getAlpha() != 255 ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.COLOR_ALPHA, getValue().getAlpha() + "", indentation );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.COLOR );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
