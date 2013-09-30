// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: www.phylosoft.org

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public final class Accession implements PhylogenyData, Comparable<Accession> {

    final private String _comment;
    final private String _source;
    final private String _source_value;
    final private String _value;

    public Accession( final String value, final String source ) {
        _value = value;
        _source = source;
        _comment = "";
        if ( source != null ) {
            _source_value = source + value;
        }
        else {
            _source_value = value;
        }
    }

    public Accession( final String value, final String source, final String comment ) {
        _value = value;
        _source = source;
        _comment = comment;
        if ( source != null ) {
            _source_value = source + value;
        }
        else {
            _source_value = value;
        }
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue() );
    }

    @Override
    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( !ForesterUtil.isEmpty( getSource() ) ) {
            sb.append( getSource() );
            sb.append( ": " );
        }
        sb.append( getValue() );
        if ( !ForesterUtil.isEmpty( getComment() ) ) {
            sb.append( " (" );
            sb.append( getComment() );
            sb.append( ")" );
        }
        return sb;
    }

    @Override
    public int compareTo( final Accession o ) {
        if ( equals( o ) ) {
            return 0;
        }
        return _source_value.compareTo( o._source_value );
    }

    @Override
    public PhylogenyData copy() {
        return new Accession( getValue(), getSource() );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            return false;
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return isEqual( ( Accession ) o );
        }
    }

    public String getComment() {
        return _comment;
    }

    public String getSource() {
        return _source;
    }

    public String getValue() {
        return _value;
    }

    @Override
    public int hashCode() {
        return _source_value.hashCode();
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        if ( this == data ) {
            return true;
        }
        if ( ( data == null ) || ( getValue() == null ) ) {
            return false;
        }
        final Accession a = ( Accession ) data;
        if ( ( getSource() != null ) && ( a.getSource() != null ) ) {
            return ( a.getValue().equals( getValue() ) && a.getSource().equals( getSource() ) );
        }
        return ( a.getValue().equals( getValue() ) );
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        sb.append( ":" );
        sb.append( NHXtags.SEQUENCE_ACCESSION );
        sb.append( ForesterUtil.replaceIllegalNhxCharacters( getValue() ) );
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( ForesterUtil.isEmpty( getSource() ) ) {
            if ( ForesterUtil.isEmpty( getComment() ) ) {
                PhylogenyDataUtil.appendElement( writer,
                                                 PhyloXmlMapping.ACCESSION,
                                                 getValue(),
                                                 PhyloXmlMapping.ACCESSION_SOURCE_ATTR,
                                                 "unknown",
                                                 indentation );
            }
            else {
                PhylogenyDataUtil.appendElement( writer,
                                                 PhyloXmlMapping.ACCESSION,
                                                 getValue(),
                                                 PhyloXmlMapping.ACCESSION_SOURCE_ATTR,
                                                 "unknown",
                                                 PhyloXmlMapping.ACCESSION_COMMENT_ATTR,
                                                 getComment(),
                                                 indentation );
            }
        }
        else {
            if ( ForesterUtil.isEmpty( getComment() ) ) {
                PhylogenyDataUtil.appendElement( writer,
                                                 PhyloXmlMapping.ACCESSION,
                                                 getValue(),
                                                 PhyloXmlMapping.ACCESSION_SOURCE_ATTR,
                                                 getSource(),
                                                 indentation );
            }
            else {
                PhylogenyDataUtil.appendElement( writer,
                                                 PhyloXmlMapping.ACCESSION,
                                                 getValue(),
                                                 PhyloXmlMapping.ACCESSION_SOURCE_ATTR,
                                                 getSource(),
                                                 PhyloXmlMapping.ACCESSION_COMMENT_ATTR,
                                                 getComment(),
                                                 indentation );
            }
        }
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
