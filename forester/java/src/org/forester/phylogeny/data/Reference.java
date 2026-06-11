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

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.util.ForesterUtil;

public class Reference implements PhylogenyData {

    String _desc;
    String _doi;

    public Reference( final String desc ) {
        _desc = desc;
        _doi = "";
    }

    public Reference( final String desc, final String doi ) {
        _desc = desc;
        _doi = doi;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getDescription() );
    }

    @Override
    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( !ForesterUtil.isEmpty( getDoi() ) ) {
            sb.append( "[doi:" );
            sb.append( getDoi() );
            sb.append( "] " );
        }
        sb.append( getDescription() );
        return sb;
    }

    @Override
    public PhylogenyData copy() {
        return new Reference( getDescription(), getDoi() );
    }

    public String getDoi() {
        return _doi;
    }

    public String getDescription() {
        return _desc;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        if ( ( data == null ) || ( getDescription() == null ) ) {
            return false;
        }
        return ( ( Reference ) data ).getDescription().equals( getDescription() )
                && ( ( Reference ) data ).getDoi().equals( getDoi() );
    }

    public void setDoi( final String doi ) throws PhyloXmlDataFormatException {
        if ( !ForesterUtil.isEmpty( doi ) && !PhyloXmlUtil.LIT_REF_DOI_PATTERN.matcher( doi ).matches() ) {
            throw new PhyloXmlDataFormatException( "illegal doi: [" + doi + "]" );
        }
        _doi = doi;
    }

    public void setValue( final String value ) {
        _desc = value;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.REFERENCE, PhyloXmlMapping.REFERENCE_DOI_ATTR, getDoi() );
        if ( !ForesterUtil.isEmpty( getDescription() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.REFERENCE_DESC, getDescription(), indentation );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.REFERENCE );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}