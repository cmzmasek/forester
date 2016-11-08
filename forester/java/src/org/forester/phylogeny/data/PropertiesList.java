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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.forester.util.ForesterUtil;

public class PropertiesList implements PhylogenyData {

    private final List<Property>       _properties;
    private final Comparator<Property> comp = new Comparator<Property>() {

                                                @Override
                                                public int compare( final Property p1, final Property p2 ) {
                                                    return p2.getRef().compareTo( p1.getRef() );
                                                }
                                            };

    public PropertiesList() {
        _properties = new ArrayList<Property>();
    }

    public int size() {
        return _properties.size();
    }

    public void addProperty( final Property property ) throws IllegalArgumentException {
        _properties.add( property );
        Collections.sort( _properties, comp );
    }

    @Override
    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        boolean first = true;
        for( final Property p : getProperties() ) {
            if ( first ) {
                first = false;
            }
            else {
                sb.append( "\n" );
            }
            sb.append( p.asText() );
        }
        return sb;
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        final PropertiesList new_one = new PropertiesList();
        for( final Property r : getProperties() ) {
            new_one.addProperty( r );
        }
        return new_one;
    }

    public List<Property> getProperties() {
        return _properties;
    }

    public List<Property> getPropertiesWithGivenReferencePrefix( final String ref_prefix )
            throws IllegalArgumentException {
        if ( ForesterUtil.isEmpty( ref_prefix ) ) {
            throw new IllegalArgumentException( "reference prefix is null or empty" );
        }
        final String my_ref_prefix = new String( ref_prefix.trim() );
        final List<Property> props = new ArrayList<Property>();
        for( final Property p : getProperties() ) {
            if ( p.getRef().startsWith( my_ref_prefix ) ) {
                props.add( p );
            }
        }
        return props;
    }

    public List<Property> getProperties( final String ref ) throws IllegalArgumentException {
        final List<Property> props = new ArrayList<Property>();
        for( final Property p : getProperties() ) {
            if ( p.getRef().equals( ref ) ) {
                props.add( p );
            }
        }
        return props;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( getProperties() != null ) {
            for( final Property p : getProperties() ) {
                p.toPhyloXML( writer, level, indentation );
            }
        }
    }

    @Override
    public String toString() {
        return asSimpleText().toString();
    }
}
