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

package org.forester.ws.seqdb;

import java.util.ArrayList;
import java.util.List;

import org.forester.util.ForesterUtil;

public final class UniProtTaxonomy {

    public static final String ARCHAEA            = "Archaea";
    public static final String BACTERIA           = "Bacteria";
    public static final String EUKARYOTA          = "Eukaryota";
    private final List<String> _lineage;
    private final String       _code;
    private final String       _scientific_name;
    private final String       _common_name;
    private final String       _synonym;
    private final String       _rank;
    private final String       _id;
    public final static String CELLULAR_ORGANISMS = "cellular organisms";
    public final static String VIRUSES            = "Viruses";
    public static final String X                  = "x";

    public UniProtTaxonomy( final String line ) {
        final String[] items = line.split( "\t" );
        if ( items.length < 5 ) {
            throw new IllegalArgumentException( "cannot parse uniprot taxonomy from: " + line );
        }
        _id = items[ 0 ].trim();
        _code = items[ 1 ].trim();
        _scientific_name = items[ 2 ].trim();
        _common_name = items[ 3 ].trim();
        _synonym = items[ 4 ].trim();
        if ( items.length > 6 ) {
            _rank = items[ 7 ].trim();
        }
        else {
            _rank = "";
        }
        String[] lin = null;
        if ( items.length > 8 ) {
            lin = items[ 8 ].split( "; " );
        }
        _lineage = new ArrayList<String>();
        if ( ( lin != null ) && ( lin.length > 0 ) ) {
            final List<String> temp = new ArrayList<String>();
            for( final String t : lin ) {
                if ( !ForesterUtil.isEmpty( t ) ) {
                    temp.add( t.trim() );
                }
            }
            for( int i = 0; i < temp.size(); ++i ) {
                if ( ( i == 0 )
                        && ( temp.get( i ).equalsIgnoreCase( EUKARYOTA ) || temp.get( i ).equalsIgnoreCase( BACTERIA ) || temp
                                .get( i ).equalsIgnoreCase( ARCHAEA ) ) ) {
                    _lineage.add( CELLULAR_ORGANISMS );
                }
                _lineage.add( temp.get( i ) );
            }
        }
        if ( _lineage.isEmpty()
                && ( _scientific_name.equalsIgnoreCase( EUKARYOTA ) || _scientific_name.equalsIgnoreCase( BACTERIA ) || _scientific_name
                        .equalsIgnoreCase( ARCHAEA ) ) ) {
            _lineage.add( CELLULAR_ORGANISMS );
        }
        _lineage.add( _scientific_name );
        if ( _lineage.isEmpty() ) {
            throw new IllegalArgumentException( "lineage in a UniProt taxonomy can not be empty\n: " + line );
        }
    }

    public UniProtTaxonomy( final List<String> lineage,
                            final String code,
                            final String common_name,
                            final String scientific_name,
                            final String synonym,
                            final String rank,
                            final String id ) {
        _lineage = lineage;
        _code = code;
        _scientific_name = scientific_name;
        _common_name = common_name;
        _synonym = synonym;
        _rank = rank;
        _id = id;
        if ( ( ( _lineage != null ) && _lineage.isEmpty() )
                || ( ( !ForesterUtil.isEmpty( _lineage ) ) && !_lineage.get( _lineage.size() - 1 )
                        .equalsIgnoreCase( _scientific_name ) ) ) {
            _lineage.add( _scientific_name );
        }
    }

    /**
     * Creates deep copy for all fields, except lineage.
     *
     * @return
     */
    public UniProtTaxonomy copy() {
        return new UniProtTaxonomy( getLineage(),
                                    getCode() != null ? new String( getCode() ) : null,
                                            getCommonName() != null ? new String( getCommonName() ) : null,
                                                    getScientificName() != null ? new String( getScientificName() ) : null,
                                                            getSynonym() != null ? new String( getSynonym() ) : null,
                                                                    getRank() != null ? new String( getRank() ) : null,
                                                                            getId() != null ? new String( getId() ) : null );
    }

    public String getCode() {
        return _code;
    }

    public String getCommonName() {
        return _common_name;
    }

    public String getId() {
        return _id;
    }

    public List<String> getLineage() {
        return _lineage;
    }

    public String getRank() {
        return _rank;
    }

    public String getScientificName() {
        return _scientific_name;
    }

    public String getSynonym() {
        return _synonym;
    }

    public final static UniProtTaxonomy createSpecialFromScientificName( final String sn ) {
        final List<String> lineage = new ArrayList<String>();
        final String code = "";
        final String common_name = "";
        String scientific_name = "";
        final String synonym = "";
        String rank = "";
        String id = "";
        if ( sn.equalsIgnoreCase( BACTERIA ) ) {
            scientific_name = BACTERIA;
            lineage.add( "cellular organisms" );
            rank = "superkingdom";
            id = "2";
        }
        else if ( sn.equalsIgnoreCase( ARCHAEA ) ) {
            scientific_name = ARCHAEA;
            lineage.add( "cellular organisms" );
            rank = "superkingdom";
            id = "2157";
        }
        else if ( sn.equalsIgnoreCase( EUKARYOTA ) ) {
            scientific_name = EUKARYOTA;
            lineage.add( "cellular organisms" );
            rank = "superkingdom";
            id = "2759";
        }
        else if ( sn.equalsIgnoreCase( VIRUSES ) ) {
            scientific_name = VIRUSES;
            rank = "superkingdom";
            id = "10239";
        }
        else if ( sn.equalsIgnoreCase( X ) ) {
            scientific_name = X;
        }
        else {
            throw new IllegalArgumentException( "illegal attempt to make UniProt taxonomy for :" + sn );
        }
        return new UniProtTaxonomy( lineage, code, common_name, scientific_name, synonym, rank, id );
    }
}
