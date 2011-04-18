// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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
// WWW: www.phylosoft.org/forester

package org.forester.ws.uniprot;

import java.util.List;

public final class UniProtTaxonomyEntry {

    private String _id;
    private String _ac;
    private String _rec_name;
    private String _os_scientific_name;
    private String _os_common_name;
    private String _tax_id;

    private UniProtTaxonomyEntry() {
    }

    public static UniProtTaxonomyEntry createInstanceFromPlainText( final List<String> lines ) {
        final UniProtTaxonomyEntry e = new UniProtTaxonomyEntry();
        for( final String line : lines ) {
            if ( line.startsWith( "ID" ) ) {
                e.setId( line.split( "\\s+" )[ 1 ] );
            }
            else if ( line.startsWith( "AC" ) ) {
                e.setAc( extract( line, "AC", ";" ) );
            }
            else if ( line.startsWith( "DE" ) ) {
                if ( ( line.indexOf( "RecName:" ) > 0 ) && ( line.indexOf( "Full=" ) > 0 ) ) {
                    e.setRecName( extract( line, "Full=", ";" ) );
                }
                if ( ( line.indexOf( "RecName:" ) > 0 ) && ( line.indexOf( "Full=" ) > 0 ) ) {
                    e.setRecName( extract( line, "Full=", ";" ) );
                }
            }
            else if ( line.startsWith( "OS" ) ) {
                if ( line.indexOf( "(" ) > 0 ) {
                    e.setOsScientificName( extract( line, "OS", "(" ) );
                }
                else {
                    e.setOsScientificName( extract( line, "OS", "." ) );
                }
            }
            else if ( line.startsWith( "OX" ) ) {
                e.setTaxId( extract( line, "OX", ";" ) );
            }
        }
        return e;
    }

    private static String extract( final String target, final String a, final String b ) {
        final int i_a = target.indexOf( a );
        final int i_b = target.indexOf( b );
        if ( ( i_a < 0 ) || ( i_b < i_a ) ) {
            throw new IllegalArgumentException( "attempt to extract from [" + target + "] between [" + a + "] and ["
                    + b + "]" );
        }
        return target.substring( i_a + a.length() + 1, i_b - 1 ).trim();
    }

    public String getId() {
        return _id;
    }

    private void setId( final String id ) {
        _id = id;
    }

    public String getAc() {
        return _ac;
    }

    private void setAc( final String ac ) {
        _ac = ac;
    }

    public String getRecName() {
        return _rec_name;
    }

    private void setRecName( final String rec_name ) {
        _rec_name = rec_name;
    }

    public String getOsScientificName() {
        return _os_scientific_name;
    }

    private void setOsScientificName( final String os_scientific_name ) {
        _os_scientific_name = os_scientific_name;
    }

    public String getOsCommonName() {
        return _os_common_name;
    }

    private void setOsCommonName( final String os_common_name ) {
        _os_common_name = os_common_name;
    }

    public String getTaxId() {
        return _tax_id;
    }

    private void setTaxId( final String tax_id ) {
        _tax_id = tax_id;
    }
}
