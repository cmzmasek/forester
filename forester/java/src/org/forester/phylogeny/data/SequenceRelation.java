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
import java.util.LinkedHashMap;
import java.util.Map;

public class SequenceRelation implements PhylogenyData {

    //public final static Map<String, SEQUENCE_RELATION_TYPE> typesToNames = new LinkedHashMap<String, SEQUENCE_RELATION_TYPE>();
    public final static Map<SEQUENCE_RELATION_TYPE, String> typesToNames                                = new LinkedHashMap<SEQUENCE_RELATION_TYPE, String>();
    public final static String                              SEQUENCE_RELATION_TYPE_ORTHOLOGY            = "orthology";
    public final static String                              SEQUENCE_RELATION_TYPE_ONE_TO_ONE_ORTHOLOGY = "one_to_one_orthology";
    public final static String                              SEQUENCE_RELATION_TYPE_SUPER_ORTHOLOGY      = "super_orthology";
    public final static String                              SEQUENCE_RELATION_TYPE_PARALOGY             = "paralogy";
    public final static String                              SEQUENCE_RELATION_TYPE_ULTRA_PARALOGY       = "ultra_paralogy";
    public final static String                              SEQUENCE_RELATION_TYPE_XENOLOGY             = "xenology";
    public final static String                              SEQUENCE_RELATION_TYPE_UNKNOWN              = "unknown";
    public final static String                              SEQUENCE_RELATION_TYPE_OTHER                = "other";
    private Sequence                                        ref0;
    private Sequence                                        ref1;
    private SEQUENCE_RELATION_TYPE                          type;
    private Double                                          distance;
    private Confidence                                      confidence;
    static {
        typesToNames.put( SEQUENCE_RELATION_TYPE.orthology, SEQUENCE_RELATION_TYPE_ORTHOLOGY );
        typesToNames.put( SEQUENCE_RELATION_TYPE.one_to_one_orthology, SEQUENCE_RELATION_TYPE_ONE_TO_ONE_ORTHOLOGY );
        typesToNames.put( SEQUENCE_RELATION_TYPE.super_orthology, SEQUENCE_RELATION_TYPE_SUPER_ORTHOLOGY );
        typesToNames.put( SEQUENCE_RELATION_TYPE.paralogy, SEQUENCE_RELATION_TYPE_PARALOGY );
        typesToNames.put( SEQUENCE_RELATION_TYPE.ultra_paralogy, SEQUENCE_RELATION_TYPE_ULTRA_PARALOGY );
        typesToNames.put( SEQUENCE_RELATION_TYPE.xenology, SEQUENCE_RELATION_TYPE_XENOLOGY );
        typesToNames.put( SEQUENCE_RELATION_TYPE.unknown, SEQUENCE_RELATION_TYPE_UNKNOWN );
        typesToNames.put( SEQUENCE_RELATION_TYPE.other, SEQUENCE_RELATION_TYPE_OTHER );
    }

    @Override
    public StringBuffer asSimpleText() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public StringBuffer asText() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public PhylogenyData copy() {
        // TODO Auto-generated method stub
        return null;
    }

    public Confidence getConfidence() {
        return confidence;
    }

    public Double getDistance() {
        return distance;
    }

    public Sequence getRef0() {
        return ref0;
    }

    public Sequence getRef1() {
        return ref1;
    }

    public SEQUENCE_RELATION_TYPE getType() {
        return type;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        // TODO Auto-generated method stub
        return false;
    }

    public void setConfidence( final Confidence confidence ) {
        this.confidence = confidence;
    }

    public void setDistance( final Double distance ) {
        this.distance = distance;
    }

    public void setRef0( final Sequence ref0 ) {
        this.ref0 = ref0;
    }

    public void setRef1( final Sequence ref1 ) {
        this.ref1 = ref1;
    }

    public void setType( final SEQUENCE_RELATION_TYPE type ) {
        this.type = type;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        // TODO Auto-generated method stub
    }

    public static String getPrintableNameByType( final SEQUENCE_RELATION_TYPE type ) {
        String s = typesToNames.get( type );
        if ( s != null ) {
            s = s.replace( '_', ' ' );
            if ( ( s.length() > 15 ) && s.toLowerCase().endsWith( "ology" ) ) {
                s = s.substring( 0, s.length() - 5 ) + ".";
            }
        }
        return s;
    }

    public static enum SEQUENCE_RELATION_TYPE {
        orthology, one_to_one_orthology, super_orthology, paralogy, ultra_paralogy, xenology, unknown, other;
    }
}
