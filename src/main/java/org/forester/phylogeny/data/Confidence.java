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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.util.ForesterUtil;

public class Confidence implements PhylogenyData, Comparable<Confidence> {

    public final static int          CONFIDENCE_DEFAULT_VALUE = -9999;
    private double                   _value;
    private double                   _sd;
    private String                   _type;
    public final static NumberFormat FORMATTER;
    static {
        final DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator( '.' );
        FORMATTER = new DecimalFormat( "#.#########", dfs );
        FORMATTER.setMaximumFractionDigits( PhyloXmlUtil.ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT );
    }

    public Confidence() {
        init();
    }

    public Confidence( final double value, final String type ) {
        setValue( value );
        setType( type );
        setStandardDeviation( CONFIDENCE_DEFAULT_VALUE );
    }

    public Confidence( final double value, final String type, final double sd ) {
        setValue( value );
        setType( type );
        setStandardDeviation( sd );
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer().append( ForesterUtil.FORMATTER_6.format( getValue() ) );
    }

    @Override
    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( !ForesterUtil.isEmpty( getType() ) ) {
            sb.append( "[" );
            sb.append( getType() );
            sb.append( "] " );
        }
        sb.append( ForesterUtil.FORMATTER_6.format( getValue() ) );
        if ( getStandardDeviation() != CONFIDENCE_DEFAULT_VALUE ) {
            sb.append( " (sd=" );
            sb.append( getStandardDeviation() );
            sb.append( ")" );
        }
        return sb;
    }

    @Override
    public int compareTo( final Confidence confidence ) {
        if ( this == confidence ) {
            return 0;
        }
        return getType().compareToIgnoreCase( confidence.getType() );
    }

    @Override
    public PhylogenyData copy() {
        return new Confidence( getValue(), getType(), getStandardDeviation() );
    }

    public String getType() {
        return _type;
    }

    public double getValue() {
        return _value;
    }

    public double getStandardDeviation() {
        return _sd;
    }

    public void init() {
        setValue( CONFIDENCE_DEFAULT_VALUE );
        setType( "" );
        setStandardDeviation( CONFIDENCE_DEFAULT_VALUE );
    }

    @Override
    public boolean isEqual( final PhylogenyData confidence ) {
        if ( confidence == null ) {
            return false;
        }
        if ( !( confidence instanceof Confidence ) ) {
            return false;
        }
        final Confidence s = ( Confidence ) confidence;
        if ( s.getValue() != getValue() ) {
            return false;
        }
        if ( !s.getType().equals( getType() ) ) {
            return false;
        }
        return true;
    }

    public void setType( final String type ) {
        _type = type;
    }

    public void setValue( final double value ) {
        _value = value;
    }

    public void setStandardDeviation( final double sd ) {
        _sd = sd;
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        sb.append( NHXtags.SUPPORT );
        sb.append( FORMATTER.format( ForesterUtil.round( getValue(),
                                                         PhyloXmlUtil.ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT ) ) );
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( getValue() == CONFIDENCE_DEFAULT_VALUE ) {
            return;
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        if ( getStandardDeviation() != CONFIDENCE_DEFAULT_VALUE ) {
            PhylogenyDataUtil
            .appendElement( writer,
                            PhyloXmlMapping.CONFIDENCE,
                            FORMATTER.format( ForesterUtil
                                              .round( getValue(), PhyloXmlUtil.ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT ) ),
                                              PhyloXmlMapping.CONFIDENCE_TYPE_ATTR,
                                              ForesterUtil.isEmpty( getType() ) ? "unknown" : getType(),
                                                      PhyloXmlMapping.CONFIDENCE_SD_ATTR,
                                                      String.valueOf( ForesterUtil
                                                                      .round( getStandardDeviation(),
                                                                              PhyloXmlUtil.ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT ) ) );
        }
        else {
            PhylogenyDataUtil
            .appendElement( writer,
                            PhyloXmlMapping.CONFIDENCE,
                            FORMATTER.format( ForesterUtil
                                              .round( getValue(), PhyloXmlUtil.ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT ) ),
                                              PhyloXmlMapping.CONFIDENCE_TYPE_ATTR,
                                              ForesterUtil.isEmpty( getType() ) ? "unknown" : getType() );
        }
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
