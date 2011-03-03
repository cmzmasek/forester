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

package org.forester.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.ws.uniprot.UniProtTaxonomy;
import org.forester.ws.uniprot.UniProtWsTools;

public final class AncestralTaxonomyInference {

    private static final int                              MAX_CACHE_SIZE           = 100000;
    private static final int                              MAX_TAXONOMIES_TO_RETURN = 100;
    private static final HashMap<String, UniProtTaxonomy> _sn_up_cache_map         = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _code_up_cache_map       = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _cn_up_cache_map         = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _id_up_cache_map         = new HashMap<String, UniProtTaxonomy>();

    synchronized private static void clearCachesIfTooLarge() {
        if ( getSnTaxCacheMap().size() > MAX_CACHE_SIZE ) {
            getSnTaxCacheMap().clear();
        }
        if ( getCnTaxCacheMap().size() > MAX_CACHE_SIZE ) {
            getCnTaxCacheMap().clear();
        }
        if ( getCodeTaxCacheMap().size() > MAX_CACHE_SIZE ) {
            getCodeTaxCacheMap().clear();
        }
        if ( getIdTaxCacheMap().size() > MAX_CACHE_SIZE ) {
            getIdTaxCacheMap().clear();
        }
    }

    synchronized private static HashMap<String, UniProtTaxonomy> getCnTaxCacheMap() {
        return _cn_up_cache_map;
    }

    synchronized private static HashMap<String, UniProtTaxonomy> getCodeTaxCacheMap() {
        return _code_up_cache_map;
    }

    synchronized private static HashMap<String, UniProtTaxonomy> getIdTaxCacheMap() {
        return _id_up_cache_map;
    }

    synchronized private static HashMap<String, UniProtTaxonomy> getSnTaxCacheMap() {
        return _sn_up_cache_map;
    }

    synchronized private static UniProtTaxonomy getTaxonomies( final HashMap<String, UniProtTaxonomy> cache,
                                                               final String query,
                                                               final QUERY_TYPE qt ) throws IOException {
        if ( cache.containsKey( query ) ) {
            return cache.get( query ).copy();
        }
        else {
            List<UniProtTaxonomy> up_taxonomies = null;
            switch ( qt ) {
                case ID:
                    up_taxonomies = getTaxonomiesFromId( query );
                    break;
                case CODE:
                    up_taxonomies = getTaxonomiesFromTaxonomyCode( query );
                    break;
                case SN:
                    up_taxonomies = getTaxonomiesFromScientificName( query );
                    break;
                case CN:
                    up_taxonomies = getTaxonomiesFromCommonName( query );
                    break;
                default:
                    throw new RuntimeException();
            }
            if ( ( up_taxonomies != null ) && ( up_taxonomies.size() == 1 ) ) {
                final UniProtTaxonomy up_tax = up_taxonomies.get( 0 );
                if ( !ForesterUtil.isEmpty( up_tax.getScientificName() ) ) {
                    getSnTaxCacheMap().put( up_tax.getScientificName(), up_tax );
                }
                if ( !ForesterUtil.isEmpty( up_tax.getCode() ) ) {
                    getCodeTaxCacheMap().put( up_tax.getCode(), up_tax );
                }
                if ( !ForesterUtil.isEmpty( up_tax.getCommonName() ) ) {
                    getCnTaxCacheMap().put( up_tax.getCommonName(), up_tax );
                }
                if ( !ForesterUtil.isEmpty( up_tax.getId() ) ) {
                    getIdTaxCacheMap().put( up_tax.getId(), up_tax );
                }
                return up_tax;
            }
            else {
                return null;
            }
        }
    }

    synchronized private static List<UniProtTaxonomy> getTaxonomiesFromCommonName( final String query )
            throws IOException {
        return UniProtWsTools.getTaxonomiesFromCommonNameStrict( query, MAX_TAXONOMIES_TO_RETURN );
    }

    synchronized private static List<UniProtTaxonomy> getTaxonomiesFromId( final String query ) throws IOException {
        return UniProtWsTools.getTaxonomiesFromId( query, MAX_TAXONOMIES_TO_RETURN );
    }

    synchronized private static List<UniProtTaxonomy> getTaxonomiesFromScientificName( final String query )
            throws IOException {
        return UniProtWsTools.getTaxonomiesFromScientificNameStrict( query, MAX_TAXONOMIES_TO_RETURN );
    }

    synchronized private static List<UniProtTaxonomy> getTaxonomiesFromTaxonomyCode( final String query )
            throws IOException {
        return UniProtWsTools.getTaxonomiesFromTaxonomyCode( query, MAX_TAXONOMIES_TO_RETURN );
    }

    synchronized public static SortedSet<String> inferTaxonomyFromDescendents( final Phylogeny phy ) throws IOException {
        clearCachesIfTooLarge();
        final SortedSet<String> not_found = new TreeSet<String>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            // final QUERY_TYPE qt = null;
            // Taxonomy tax = null;
            // if ( node.getNodeData().isHasTaxonomy() ) {
            // tax = node.getNodeData().getTaxonomy();
            // }
            // UniProtTaxonomy up_tax = null;
            // if ( ( tax != null )
            // && ( isHasAppropriateId( tax ) || !ForesterUtil.isEmpty(
            // tax.getScientificName() )
            // || !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ||
            // !ForesterUtil.isEmpty( tax
            // .getCommonName() ) ) ) {
            // final String query = null;
            // up_tax = obtainUniProtTaxonomy( tax, query, qt );
            // if ( up_tax == null ) {
            // not_found.add( query );
            // }
            // else {
            // updateTaxonomy( qt, node, tax, up_tax );
            // }
            // }
            if ( !node.isExternal() ) {
                inferTaxonomyFromDescendents( node, not_found );
            }
        }
        return not_found;
    }

    synchronized private static void inferTaxonomyFromDescendents( final PhylogenyNode n,
                                                                   final SortedSet<String> not_found )
            throws IOException {
        if ( n.isExternal() ) {
            throw new IllegalArgumentException( "attempt to infer taxonomy from descendants of external node" );
        }
        n.getNodeData().setTaxonomy( null );
        final List<PhylogenyNode> descs = n.getDescendants();
        final List<String[]> lineages = new ArrayList<String[]>();
        int shortest_lin_length = Integer.MAX_VALUE;
        for( final PhylogenyNode desc : descs ) {
            if ( desc.getNodeData().isHasTaxonomy()
                    && ( isHasAppropriateId( desc.getNodeData().getTaxonomy() )
                            || !ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getScientificName() )
                            || !ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getTaxonomyCode() ) || !ForesterUtil
                            .isEmpty( desc.getNodeData().getTaxonomy().getCommonName() ) ) ) {
                final QUERY_TYPE qt = null;
                final String query = null;
                final UniProtTaxonomy up_tax = obtainUniProtTaxonomy( desc.getNodeData().getTaxonomy(), query, qt );
                String[] lineage = null;
                if ( up_tax != null ) {
                    lineage = obtainLineagePlusOwnScientificName( up_tax );
                }
                if ( ( lineage == null ) || ( lineage.length < 1 ) ) {
                    not_found.add( desc.getNodeData().getTaxonomy().asText().toString() );
                    return;
                }
                if ( lineage.length < shortest_lin_length ) {
                    shortest_lin_length = lineage.length;
                }
                lineages.add( lineage );
            }
            else {
                String msg = "Node(s) with no or inappropriate taxonomic information found";
                if ( !ForesterUtil.isEmpty( desc.getName() ) ) {
                    msg = "Node " + desc.getName() + " has no or inappropriate taxonomic information";
                }
                throw new IllegalArgumentException( msg );
            }
        }
        String last_common_lineage = null;
        if ( shortest_lin_length > 0 ) {
            I: for( int i = 0; i < shortest_lin_length; ++i ) {
                final String lineage_0 = lineages.get( 0 )[ i ];
                for( int j = 1; j < lineages.size(); ++j ) {
                    if ( !lineage_0.equals( lineages.get( j )[ i ] ) ) {
                        break I;
                    }
                }
                last_common_lineage = lineage_0;
            }
        }
        if ( last_common_lineage == null ) {
            return;
        }
        // if ( !n.getNodeData().isHasTaxonomy() ) {
        // n.getNodeData().setTaxonomy( new Taxonomy() );
        // }
        final Taxonomy tax = new Taxonomy();
        n.getNodeData().setTaxonomy( tax );
        tax.setScientificName( last_common_lineage );
        final UniProtTaxonomy up_tax = obtainUniProtTaxonomyFromSn( last_common_lineage );
        if ( up_tax != null ) {
            if ( !ForesterUtil.isEmpty( up_tax.getRank() ) ) {
                try {
                    tax.setRank( up_tax.getRank().toLowerCase() );
                }
                catch ( final PhyloXmlDataFormatException ex ) {
                    tax.setRank( "" );
                }
            }
            if ( !ForesterUtil.isEmpty( up_tax.getId() ) ) {
                tax.setIdentifier( new Identifier( up_tax.getId(), "uniprot" ) );
            }
            if ( !ForesterUtil.isEmpty( up_tax.getCommonName() ) ) {
                tax.setCommonName( up_tax.getCommonName() );
            }
            if ( !ForesterUtil.isEmpty( up_tax.getSynonym() ) && !tax.getSynonyms().contains( up_tax.getSynonym() ) ) {
                tax.getSynonyms().add( up_tax.getSynonym() );
            }
        }
        for( final PhylogenyNode desc : descs ) {
            if ( !desc.isExternal() && desc.getNodeData().isHasTaxonomy()
                    && desc.getNodeData().getTaxonomy().isEqual( tax ) ) {
                desc.getNodeData().setTaxonomy( null );
            }
        }
    }

    synchronized private static boolean isHasAppropriateId( final Taxonomy tax ) {
        return ( ( tax.getIdentifier() != null ) && ( !ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) && ( tax
                .getIdentifier().getProvider().equalsIgnoreCase( "ncbi" )
                || tax.getIdentifier().getProvider().equalsIgnoreCase( "uniprot" ) || tax.getIdentifier().getProvider()
                .equalsIgnoreCase( "uniprotkb" ) ) ) );
    }

    synchronized public static SortedSet<String> obtainDetailedTaxonomicInformation( final Phylogeny phy )
            throws IOException {
        clearCachesIfTooLarge();
        final SortedSet<String> not_found = new TreeSet<String>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            final QUERY_TYPE qt = null;
            Taxonomy tax = null;
            if ( node.getNodeData().isHasTaxonomy() ) {
                tax = node.getNodeData().getTaxonomy();
            }
            else if ( node.isExternal() ) {
                if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                    not_found.add( node.getName() );
                }
                else {
                    not_found.add( node.toString() );
                }
            }
            UniProtTaxonomy up_tax = null;
            if ( ( tax != null )
                    && ( isHasAppropriateId( tax ) || !ForesterUtil.isEmpty( tax.getScientificName() )
                            || !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) || !ForesterUtil.isEmpty( tax
                            .getCommonName() ) ) ) {
                up_tax = obtainUniProtTaxonomy( tax, null, qt );
                if ( up_tax != null ) {
                    updateTaxonomy( qt, node, tax, up_tax );
                }
                else {
                    not_found.add( tax.toString() );
                }
            }
        }
        return not_found;
    }

    synchronized private static String[] obtainLineagePlusOwnScientificName( final UniProtTaxonomy up_tax ) {
        final String[] lineage = up_tax.getLineage();
        final String[] lin_plus_self = new String[ lineage.length + 1 ];
        for( int i = 0; i < lineage.length; ++i ) {
            lin_plus_self[ i ] = lineage[ i ];
        }
        lin_plus_self[ lineage.length ] = up_tax.getScientificName();
        return lin_plus_self;
    }

    synchronized private static UniProtTaxonomy obtainUniProtTaxonomy( final Taxonomy tax, String query, QUERY_TYPE qt )
            throws IOException {
        if ( isHasAppropriateId( tax ) ) {
            query = tax.getIdentifier().getValue();
            qt = QUERY_TYPE.ID;
            return getTaxonomies( getIdTaxCacheMap(), query, qt );
        }
        else if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            query = tax.getScientificName();
            qt = QUERY_TYPE.SN;
            return getTaxonomies( getSnTaxCacheMap(), query, qt );
        }
        else if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
            query = tax.getTaxonomyCode();
            qt = QUERY_TYPE.CODE;
            return getTaxonomies( getCodeTaxCacheMap(), query, qt );
        }
        else {
            query = tax.getCommonName();
            qt = QUERY_TYPE.CN;
            return getTaxonomies( getCnTaxCacheMap(), query, qt );
        }
    }

    synchronized private static UniProtTaxonomy obtainUniProtTaxonomyFromSn( final String sn ) throws IOException {
        UniProtTaxonomy up_tax = null;
        if ( getSnTaxCacheMap().containsKey( sn ) ) {
            up_tax = getSnTaxCacheMap().get( sn ).copy();
        }
        else {
            final List<UniProtTaxonomy> up_taxonomies = getTaxonomiesFromScientificName( sn );
            if ( ( up_taxonomies != null ) && ( up_taxonomies.size() == 1 ) ) {
                up_tax = up_taxonomies.get( 0 );
                getSnTaxCacheMap().put( sn, up_tax );
                if ( !ForesterUtil.isEmpty( up_tax.getCode() ) ) {
                    getCodeTaxCacheMap().put( up_tax.getCode(), up_tax );
                }
                if ( !ForesterUtil.isEmpty( up_tax.getCommonName() ) ) {
                    getCnTaxCacheMap().put( up_tax.getCommonName(), up_tax );
                }
                if ( !ForesterUtil.isEmpty( up_tax.getId() ) ) {
                    getIdTaxCacheMap().put( up_tax.getId(), up_tax );
                }
            }
        }
        return up_tax;
    }

    synchronized private static void updateTaxonomy( final QUERY_TYPE qt,
                                                     final PhylogenyNode node,
                                                     final Taxonomy tax,
                                                     final UniProtTaxonomy up_tax ) {
        if ( ( qt != QUERY_TYPE.SN ) && !ForesterUtil.isEmpty( up_tax.getScientificName() )
                && ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            tax.setScientificName( up_tax.getScientificName() );
        }
        if ( node.isExternal()
                && ( ( qt != QUERY_TYPE.CODE ) && !ForesterUtil.isEmpty( up_tax.getCode() ) && ForesterUtil
                        .isEmpty( tax.getTaxonomyCode() ) ) ) {
            tax.setTaxonomyCode( up_tax.getCode() );
        }
        if ( ( qt != QUERY_TYPE.CN ) && !ForesterUtil.isEmpty( up_tax.getCommonName() )
                && ForesterUtil.isEmpty( tax.getCommonName() ) ) {
            tax.setCommonName( up_tax.getCommonName() );
        }
        if ( !ForesterUtil.isEmpty( up_tax.getSynonym() ) && !tax.getSynonyms().contains( up_tax.getSynonym() ) ) {
            tax.getSynonyms().add( up_tax.getSynonym() );
        }
        if ( !ForesterUtil.isEmpty( up_tax.getRank() ) && ForesterUtil.isEmpty( tax.getRank() ) ) {
            try {
                tax.setRank( up_tax.getRank().toLowerCase() );
            }
            catch ( final PhyloXmlDataFormatException ex ) {
                tax.setRank( "" );
            }
        }
        if ( ( qt != QUERY_TYPE.ID ) && !ForesterUtil.isEmpty( up_tax.getId() ) && ( tax.getIdentifier() == null ) ) {
            tax.setIdentifier( new Identifier( up_tax.getId(), "uniprot" ) );
        }
    }

    private enum QUERY_TYPE {
        CODE, SN, CN, ID;
    }
}
