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
    private static final int                              MAX_TAXONOMIES_TO_RETURN = 10;
    private static final HashMap<String, UniProtTaxonomy> _sn_up_cache_map         = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _lineage_up_cache_map    = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _code_up_cache_map       = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _cn_up_cache_map         = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _id_up_cache_map         = new HashMap<String, UniProtTaxonomy>();

    synchronized private static void clearCachesIfTooLarge() {
        if ( getSnTaxCacheMap().size() > MAX_CACHE_SIZE ) {
            getSnTaxCacheMap().clear();
        }
        if ( getLineageTaxCacheMap().size() > MAX_CACHE_SIZE ) {
            getLineageTaxCacheMap().clear();
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

    synchronized private static HashMap<String, UniProtTaxonomy> getLineageTaxCacheMap() {
        return _lineage_up_cache_map;
    }

    synchronized private static UniProtTaxonomy getTaxonomies( final HashMap<String, UniProtTaxonomy> cache,
                                                               final Object query,
                                                               final QUERY_TYPE qt ) throws IOException,
            AncestralTaxonomyInferenceException {
        if ( cache.containsKey( query ) ) {
            return cache.get( query ).copy();
        }
        else {
            List<UniProtTaxonomy> up_taxonomies = null;
            switch ( qt ) {
                case ID:
                    up_taxonomies = getTaxonomiesFromId( ( String ) query );
                    break;
                case CODE:
                    up_taxonomies = getTaxonomiesFromTaxonomyCode( ( String ) query );
                    break;
                case SN:
                    up_taxonomies = getTaxonomiesFromScientificName( ( String ) query );
                    break;
                case CN:
                    up_taxonomies = getTaxonomiesFromCommonName( ( String ) query );
                    break;
                case LIN:
                    return obtainUniProtTaxonomyFromLineage( ( List<String> ) query );
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

    synchronized public static void inferTaxonomyFromDescendents( final Phylogeny phy ) throws IOException,
            AncestralTaxonomyInferenceException {
        clearCachesIfTooLarge();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( !node.isExternal() ) {
                inferTaxonomyFromDescendents( node );
            }
        }
    }

    synchronized private static void inferTaxonomyFromDescendents( final PhylogenyNode n ) throws IOException,
            AncestralTaxonomyInferenceException {
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
                            || !ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getLineage() )
                            || !ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getTaxonomyCode() ) || !ForesterUtil
                            .isEmpty( desc.getNodeData().getTaxonomy().getCommonName() ) ) ) {
                final UniProtTaxonomy up_tax = obtainUniProtTaxonomy( desc.getNodeData().getTaxonomy(), null, null );
                if ( ( up_tax == null ) && ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getLineage() ) ) {
                    String desc_str = "";
                    if ( !ForesterUtil.isEmpty( desc.getName() ) ) {
                        desc_str = "\"" + desc.getName() + "\"";
                    }
                    else {
                        desc_str = "[" + desc.getId() + "]";
                    }
                    System.out.println( desc.getNodeData().getTaxonomy().toString() );
                    System.out.println( ForesterUtil.stringListToString( desc.getNodeData().getTaxonomy().getLineage(),
                                                                         "  >  " ) );
                    throw new AncestralTaxonomyInferenceException( "a taxonomy for node " + desc_str
                            + " could not be established from the database" );
                }
                String[] lineage = ForesterUtil.stringListToArray( desc.getNodeData().getTaxonomy().getLineage() );
                if ( ( lineage == null ) || ( lineage.length < 1 ) ) {
                    lineage = ForesterUtil.stringListToArray( up_tax.getLineage() );
                }
                if ( ( lineage == null ) || ( lineage.length < 1 ) ) {
                    throw new AncestralTaxonomyInferenceException( "a taxonomic lineage for node \""
                            + desc.getNodeData().getTaxonomy().toString() + "\" could not be established" );
                }
                if ( lineage.length < shortest_lin_length ) {
                    shortest_lin_length = lineage.length;
                }
                lineages.add( lineage );
            }
            else {
                String node = "";
                if ( !ForesterUtil.isEmpty( desc.getName() ) ) {
                    node = "\"" + desc.getName() + "\"";
                }
                else {
                    node = "[" + desc.getId() + "]";
                }
                //   final List<PhylogenyNode> e = desc.getAllExternalDescendants();
                //TODO remove me!
                //                System.out.println();
                //                int x = 0;
                //                for( final PhylogenyNode object : e ) {
                //                    System.out.println( x + ":" );
                //                    System.out.println( object.getName() + "  " );
                //                    x++;
                //                }
                //                System.out.println();
                //
                throw new AncestralTaxonomyInferenceException( "node " + node
                        + " has no or inappropriate taxonomic information" );
            }
        }
        final List<String> last_common_lineage = new ArrayList<String>();
        String last_common = null;
        if ( shortest_lin_length > 0 ) {
            I: for( int i = 0; i < shortest_lin_length; ++i ) {
                final String lineage_0 = lineages.get( 0 )[ i ];
                for( int j = 1; j < lineages.size(); ++j ) {
                    if ( !lineage_0.equals( lineages.get( j )[ i ] ) ) {
                        break I;
                    }
                }
                last_common_lineage.add( lineage_0 );
                last_common = lineage_0;
            }
        }
        if ( last_common_lineage.isEmpty() ) {
            boolean saw_viruses = false;
            boolean saw_cellular_organism = false;
            for( final String[] lineage : lineages ) {
                if ( lineage.length > 0 ) {
                    if ( lineage[ 0 ].equalsIgnoreCase( UniProtTaxonomy.VIRUSES ) ) {
                        saw_viruses = true;
                    }
                    else if ( lineage[ 0 ].equalsIgnoreCase( UniProtTaxonomy.CELLULAR_ORGANISMS ) ) {
                        saw_cellular_organism = true;
                    }
                    if ( saw_cellular_organism && saw_viruses ) {
                        break;
                    }
                }
            }
            if ( saw_cellular_organism && saw_viruses ) {
                last_common_lineage.add( UniProtTaxonomy.CELLULAR_ORGANISMS );
                last_common = UniProtTaxonomy.CELLULAR_ORGANISMS;
            }
            else {
                String msg = "no common lineage for:\n";
                int counter = 0;
                for( final String[] strings : lineages ) {
                    msg += counter + ": ";
                    ++counter;
                    for( final String string : strings ) {
                        msg += string + " ";
                    }
                    msg += "\n";
                }
                throw new AncestralTaxonomyInferenceException( msg );
            }
        }
        final Taxonomy tax = new Taxonomy();
        n.getNodeData().setTaxonomy( tax );
        tax.setScientificName( last_common );
        final UniProtTaxonomy up_tax = obtainUniProtTaxonomyFromLineage( last_common_lineage );
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
            if ( up_tax.getLineage() != null ) {
                tax.setLineage( new ArrayList<String>() );
                for( final String lin : up_tax.getLineage() ) {
                    if ( !ForesterUtil.isEmpty( lin ) ) {
                        tax.getLineage().add( lin );
                    }
                }
            }
        }
        if ( ForesterUtil.isEmpty( tax.getLineage() ) ) {
            tax.setLineage( new ArrayList<String>() );
            for( final String lin : last_common_lineage ) {
                if ( !ForesterUtil.isEmpty( lin ) ) {
                    tax.getLineage().add( lin );
                }
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

    synchronized public static SortedSet<String> obtainDetailedTaxonomicInformation( final Phylogeny phy,
                                                                                     final boolean delete )
            throws IOException, AncestralTaxonomyInferenceException {
        clearCachesIfTooLarge();
        final SortedSet<String> not_found = new TreeSet<String>();
        List<PhylogenyNode> not_found_external_nodes = null;
        if ( delete ) {
            not_found_external_nodes = new ArrayList<PhylogenyNode>();
        }
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
                if ( delete ) {
                    not_found_external_nodes.add( node );
                }
            }
            UniProtTaxonomy uniprot_tax = null;
            if ( ( tax != null )
                    && ( isHasAppropriateId( tax ) || !ForesterUtil.isEmpty( tax.getScientificName() )
                            || !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) || !ForesterUtil.isEmpty( tax
                            .getCommonName() ) ) ) {
                uniprot_tax = obtainUniProtTaxonomy( tax, null, qt );
                if ( uniprot_tax != null ) {
                    updateTaxonomy( qt, node, tax, uniprot_tax );
                }
                else {
                    not_found.add( tax.toString() );
                    if ( delete && node.isExternal() ) {
                        not_found_external_nodes.add( node );
                    }
                }
            }
        }
        if ( delete ) {
            for( final PhylogenyNode node : not_found_external_nodes ) {
                phy.deleteSubtree( node, true );
            }
            phy.externalNodesHaveChanged();
            phy.hashIDs();
            phy.recalculateNumberOfExternalDescendants( true );
        }
        return not_found;
    }

    synchronized public static UniProtTaxonomy obtainUniProtTaxonomy( final Taxonomy tax, Object query, QUERY_TYPE qt )
            throws IOException, AncestralTaxonomyInferenceException {
        if ( isHasAppropriateId( tax ) ) {
            query = tax.getIdentifier().getValue();
            qt = QUERY_TYPE.ID;
            return getTaxonomies( getIdTaxCacheMap(), query, qt );
        }
        else if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            if ( !ForesterUtil.isEmpty( tax.getLineage() ) ) {
                query = tax.getLineage();
                qt = QUERY_TYPE.LIN;
                return getTaxonomies( getLineageTaxCacheMap(), query, qt );
            }
            else {
                query = tax.getScientificName();
                qt = QUERY_TYPE.SN;
                return getTaxonomies( getSnTaxCacheMap(), query, qt );
            }
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

    synchronized private static UniProtTaxonomy obtainUniProtTaxonomyFromLineage( final List<String> lineage )
            throws AncestralTaxonomyInferenceException, IOException {
        final String lineage_str = ForesterUtil.stringListToString( lineage, ">" );
        UniProtTaxonomy up_tax = null;
        if ( getLineageTaxCacheMap().containsKey( lineage_str ) ) {
            up_tax = getLineageTaxCacheMap().get( lineage_str ).copy();
        }
        else {
            final List<UniProtTaxonomy> up_taxonomies = getTaxonomiesFromScientificName( lineage
                    .get( lineage.size() - 1 ) );
            if ( ( up_taxonomies != null ) && ( up_taxonomies.size() > 0 ) ) {
                for( final UniProtTaxonomy up_taxonomy : up_taxonomies ) {
                    boolean match = true;
                    I: for( int i = 0; i < lineage.size(); ++i ) {
                        if ( !lineage.get( i ).equalsIgnoreCase( up_taxonomy.getLineage().get( i ) ) ) {
                            match = false;
                            break I;
                        }
                    }
                    if ( match ) {
                        if ( up_tax != null ) {
                            throw new AncestralTaxonomyInferenceException( "lineage \""
                                    + ForesterUtil.stringListToString( lineage, " > " ) + "\" is not unique" );
                        }
                        up_tax = up_taxonomy;
                    }
                }
                if ( up_tax == null ) {
                    throw new AncestralTaxonomyInferenceException( "lineage \""
                            + ForesterUtil.stringListToString( lineage, " > " ) + "\" not found" );
                }
                getLineageTaxCacheMap().put( lineage_str, up_tax );
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
        if ( node.isExternal() && ( qt != QUERY_TYPE.CODE ) && !ForesterUtil.isEmpty( up_tax.getCode() )
                && ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
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
        if ( ( qt != QUERY_TYPE.ID ) && !ForesterUtil.isEmpty( up_tax.getId() )
                && ( ( tax.getIdentifier() == null ) || ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) ) {
            tax.setIdentifier( new Identifier( up_tax.getId(), "uniprot" ) );
        }
        if ( up_tax.getLineage() != null ) {
            tax.setLineage( new ArrayList<String>() );
            for( final String lin : up_tax.getLineage() ) {
                if ( !ForesterUtil.isEmpty( lin ) ) {
                    tax.getLineage().add( lin );
                }
            }
        }
    }

    private enum QUERY_TYPE {
        CODE, SN, CN, ID, LIN;
    }
}
