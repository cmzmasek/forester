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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.ws.wabi;

import java.io.IOException;

import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterUtil;
import org.forester.ws.wabi.TxSearch.TAX_NAME_CLASS;
import org.forester.ws.wabi.TxSearch.TAX_RANK;

public final class WabiTools {

    private static String getATxName( final Taxonomy tax ) throws IOException {
        String name = null;
        if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            name = tax.getScientificName();
        }
        else if ( !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
            name = tax.getCommonName();
        }
        if ( ForesterUtil.isEmpty( name ) ) {
            String id_value = null;
            if ( PhylogenyMethods.isTaxonomyHasIdentifierOfGivenProvider( tax, new String[] { "uniprot", "ncbi" } ) ) {
                id_value = tax.getIdentifier().getValue();
            }
            if ( !ForesterUtil.isEmpty( id_value ) ) {
                name = TxSearch.getTxName( id_value );
            }
        }
        return name;
    }

    public static String[] obtainLineage( final Taxonomy tax ) throws IOException {
        final String name = getATxName( tax );
        String result = null;
        if ( !ForesterUtil.isEmpty( name ) ) {
            result = TxSearch.searchParam( name, TAX_NAME_CLASS.ALL, TAX_RANK.ALL, 2, true );
        }
        if ( !ForesterUtil.isEmpty( result ) ) {
            final String[] lin = TxSearch.getLineage( result );
            if ( lin != null ) {
                final String[] lin_plus_self = new String[ lin.length + 1 ];
                for( int i = 0; i < lin.length; ++i ) {
                    lin_plus_self[ i ] = lin[ i ];
                }
                lin_plus_self[ lin.length ] = name;
                return lin_plus_self;
            }
        }
        return null;
    }

    public static String obtainRank( final Taxonomy tax ) throws IOException {
        final String result = searchParam( tax );
        if ( !ForesterUtil.isEmpty( result ) ) {
            return TxSearch.getTaxonomicRank( result );
        }
        return null;
    }

    private static String searchParam( final Taxonomy tax ) throws IOException {
        final String name = getATxName( tax );
        String result = null;
        if ( !ForesterUtil.isEmpty( name ) ) {
            result = TxSearch.searchParam( name, TAX_NAME_CLASS.ALL, TAX_RANK.ALL, 2, true );
        }
        return result;
    }
}
