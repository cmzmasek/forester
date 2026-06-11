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

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class BranchWidth implements PhylogenyData {

    public final static double BRANCH_WIDTH_DEFAULT_VALUE = 1.0;
    private final double       _value;

    public BranchWidth() {
        _value = BRANCH_WIDTH_DEFAULT_VALUE;
    }

    public BranchWidth( final double value ) {
        _value = value;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue() + "" );
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        return new BranchWidth( getValue() );
    }

    public double getValue() {
        return _value;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        return getValue() == ( ( BranchWidth ) data ).getValue();
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer w, final int level, final String indentation ) throws IOException {
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( indentation );
        PhylogenyDataUtil.appendElement( w, PhyloXmlMapping.WIDTH, getValue() + "" );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
