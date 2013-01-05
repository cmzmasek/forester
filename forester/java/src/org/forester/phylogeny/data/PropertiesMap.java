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
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.util.ForesterUtil;

public class PropertiesMap implements PhylogenyData {

    private final SortedMap<String, Property> _properties;

    public PropertiesMap() {
        _properties = new TreeMap<String, Property>();
    }

    public int size() {
        return _properties.size();
    }

    public void addProperty( final Property property ) throws IllegalArgumentException {
        if ( getProperties().containsKey( property.getRef() ) ) {
            throw new IllegalArgumentException( "ref [" + property.getRef() + "] is already present" );
        }
        getProperties().put( property.getRef(), property );
    }

    @Override
    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        boolean first = true;
        for( final String ref : getPropertyRefs() ) {
            if ( first ) {
                first = false;
            }
            else {
                sb.append( " " );
            }
            sb.append( getProperty( ref ).asText() );
        }
        return sb;
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        final PropertiesMap new_one = new PropertiesMap();
        for( final String r : getProperties().keySet() ) {
            new_one.addProperty( getProperties().get( r ) );
        }
        return new_one;
    }

    public SortedMap<String, Property> getProperties() {
        return _properties;
    }

    public Property[] getPropertiesArray() {
        final Property[] a = new Property[ getProperties().size() ];
        int i = 0;
        for( final String ref : getProperties().keySet() ) {
            a[ i++ ] = getProperties().get( ref );
        }
        return a;
    }

    public List<Property> getPropertiesWithGivenReferencePrefix( final String ref_prefix )
            throws IllegalArgumentException {
        if ( ForesterUtil.isEmpty( ref_prefix ) ) {
            throw new IllegalArgumentException( "reference prefix is null or empty" );
        }
        final String my_ref_prefix = new String( ref_prefix.trim() );
        final List<Property> props = new ArrayList<Property>();
        for( final String ref : getProperties().keySet() ) {
            if ( ref.startsWith( my_ref_prefix ) ) {
                props.add( getProperty( ref ) );
            }
        }
        return props;
    }

    public Property getProperty( final String ref ) throws IllegalArgumentException {
        if ( getProperties().containsKey( ref ) ) {
            return getProperties().get( ref );
        }
        else {
            throw new IllegalArgumentException( "reference [" + ref + "] is not present" );
        }
    }

    /**
     * Returns all property refs of this PhylogenyNode as String array.
     */
    public String[] getPropertyRefs() {
        if ( getProperties() == null ) {
            return new String[ 0 ];
        }
        final Property[] properties = getPropertiesArray();
        final String[] refs = new String[ properties.length ];
        for( int i = 0; i < properties.length; ++i ) {
            refs[ i ] = properties[ i ].getRef();
        }
        return refs;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public boolean refExists( final String ref ) {
        if ( getProperties() != null ) {
            for( final String r : getProperties().keySet() ) {
                if ( r.equalsIgnoreCase( ref ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    public Property removeProperty( final String ref ) throws IllegalArgumentException {
        if ( getProperties().containsKey( ref ) ) {
            return getProperties().remove( ref );
        }
        else {
            throw new IllegalArgumentException( "reference [" + ref + "] is not present" );
        }
    }

    public List<String> removePropertiesWithGivenReferencePrefix( final String ref_prefix )
            throws IllegalArgumentException {
        if ( ForesterUtil.isEmpty( ref_prefix ) ) {
            throw new IllegalArgumentException( "reference prefix is null or empty" );
        }
        final String my_ref_prefix = new String( ref_prefix.trim() );
        final List<String> to_remove = new ArrayList<String>();
        for( final String ref : getProperties().keySet() ) {
            if ( ref.startsWith( my_ref_prefix ) ) {
                to_remove.add( ref );
            }
        }
        for( final String ref : to_remove ) {
            getProperties().remove( ref );
        }
        return to_remove;
    }

    @Override
    public StringBuffer toNHX() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( getProperties() != null ) {
            for( final String ref : getProperties().keySet() ) {
                getProperties().get( ref ).toPhyloXML( writer, level, indentation );
            }
        }
    }

    @Override
    public String toString() {
        return asSimpleText().toString();
    }
}
