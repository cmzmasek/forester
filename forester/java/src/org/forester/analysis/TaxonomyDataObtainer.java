// $Id:
//
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
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.swing.JOptionPane;

import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.archaeopteryx.TreePanel;
import org.forester.archaeopteryx.tools.RunnableProcess;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.ws.uniprot.UniProtTaxonomy;
import org.forester.ws.uniprot.UniProtWsTools;

public class TaxonomyDataObtainer extends RunnableProcess {

    private static final int                              MAX_CACHE_SIZE           = 100000;
    
    private static final HashMap<String, UniProtTaxonomy> _sn_up_cache_map         = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _lineage_up_cache_map    = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _code_up_cache_map       = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _cn_up_cache_map         = new HashMap<String, UniProtTaxonomy>();
    private static final HashMap<String, UniProtTaxonomy> _id_up_cache_map         = new HashMap<String, UniProtTaxonomy>();

    
    private final Phylogeny            _phy;
    private final MainFrameApplication _mf;
    private final TreePanel            _treepanel;
    private final boolean              _delete;

    public TaxonomyDataObtainer( final MainFrameApplication mf,
                                 final TreePanel treepanel,
                                 final Phylogeny phy,
                                 final boolean delete ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
        _delete = delete;
    }

    public TaxonomyDataObtainer( final MainFrameApplication mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
        _delete = false;
    }

    private String getBaseUrl() {
        return UniProtWsTools.BASE_URL;
    }

    private void execute() {
        start( _mf, "taxonomy data" );
        SortedSet<String> not_found = null;
        try {
            not_found = obtainDetailedTaxonomicInformation( _phy, _delete );
        }
        catch ( final UnknownHostException e ) {
            JOptionPane.showMessageDialog( _mf,
                                           "Could not connect to \"" + getBaseUrl() + "\"",
                                           "Network error during taxonomic information gathering",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            JOptionPane.showMessageDialog( _mf,
                                           e.toString(),
                                           "Failed to obtain taxonomic information",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final AncestralTaxonomyInferenceException e ) {
            e.printStackTrace();
            JOptionPane.showMessageDialog( _mf,
                                           e.toString(),
                                           "Failed to obtain taxonomic information",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        finally {
            end( _mf );
        }
        if ( ( _phy == null ) || _phy.isEmpty() ) {
            try {
                JOptionPane.showMessageDialog( _mf,
                                               "None of the external node taxonomies could be resolved",
                                               "Taxonomy Tool Failed",
                                               JOptionPane.WARNING_MESSAGE );
            }
            catch ( final Exception e ) {
                // Not important if this fails, do nothing. 
            }
            return;
        }
        _treepanel.setTree( _phy );
        _mf.showWhole();
        _treepanel.setEdited( true );
        if ( ( not_found != null ) && ( not_found.size() > 0 ) ) {
            int max = not_found.size();
            boolean more = false;
            if ( max > 20 ) {
                more = true;
                max = 20;
            }
            final StringBuffer sb = new StringBuffer();
            sb.append( "Not all taxonomies could be resolved.\n" );
            if ( not_found.size() == 1 ) {
                if ( _delete ) {
                    sb.append( "The following taxonomy was not found and deleted (if external):\n" );
                }
                else {
                    sb.append( "The following taxonomy was not found:\n" );
                }
            }
            else {
                if ( _delete ) {
                    sb.append( "The following taxonomies were not found and deleted (if external) (total: "
                            + not_found.size() + "):\n" );
                }
                else {
                    sb.append( "The following taxonomies were not found (total: " + not_found.size() + "):\n" );
                }
            }
            int i = 0;
            for( final String string : not_found ) {
                if ( i > 19 ) {
                    break;
                }
                sb.append( string );
                sb.append( "\n" );
                ++i;
            }
            if ( more ) {
                sb.append( "..." );
            }
            try {
                JOptionPane.showMessageDialog( _mf,
                                               sb.toString(),
                                               "Taxonomy Tool Completed",
                                               JOptionPane.WARNING_MESSAGE );
            }
            catch ( final Exception e ) {
                // Not important if this fails, do nothing. 
            }
        }
        else {
            try {
                JOptionPane.showMessageDialog( _mf,
                                               "Taxonomy tool successfully completed",
                                               "Taxonomy Tool Completed",
                                               JOptionPane.INFORMATION_MESSAGE );
            }
            catch ( final Exception e ) {
                // Not important if this fails, do nothing.
            }
        }
    }

    
    synchronized static void clearCachesIfTooLarge() {
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
    
    synchronized private static SortedSet<String> obtainDetailedTaxonomicInformation( final Phylogeny phy,
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
    
    


    public static UniProtTaxonomy obtainUniProtTaxonomy( final Taxonomy tax, Object query, QUERY_TYPE qt )
            throws IOException, AncestralTaxonomyInferenceException {
        if ( TaxonomyDataObtainer.isHasAppropriateId( tax ) ) {
            query = tax.getIdentifier().getValue();
            qt = QUERY_TYPE.ID;
            return getTaxonomies( TaxonomyDataObtainer.getIdTaxCacheMap(), query, qt );
        }
        else if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            if ( !ForesterUtil.isEmpty( tax.getLineage() ) ) {
                query = tax.getLineage();
                qt = QUERY_TYPE.LIN;
                return getTaxonomies( TaxonomyDataObtainer.getLineageTaxCacheMap(), query, qt );
            }
            else {
                query = tax.getScientificName();
                qt = QUERY_TYPE.SN;
                return getTaxonomies( TaxonomyDataObtainer.getSnTaxCacheMap(), query, qt );
            }
        }
        else if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
            query = tax.getTaxonomyCode();
            qt = QUERY_TYPE.CODE;
            return getTaxonomies( TaxonomyDataObtainer.getCodeTaxCacheMap(), query, qt );
        }
        else {
            query = tax.getCommonName();
            qt = QUERY_TYPE.CN;
            return getTaxonomies( TaxonomyDataObtainer.getCnTaxCacheMap(), query, qt );
        }
    }

    
    
    
    static boolean isHasAppropriateId( final Taxonomy tax ) {
        return ( ( tax.getIdentifier() != null ) && ( !ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) && ( tax
                .getIdentifier().getProvider().equalsIgnoreCase( "ncbi" )
                || tax.getIdentifier().getProvider().equalsIgnoreCase( "uniprot" ) || tax.getIdentifier().getProvider()
                .equalsIgnoreCase( "uniprotkb" ) ) ) );
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


    synchronized  static HashMap<String, UniProtTaxonomy> getCnTaxCacheMap() {
        return _cn_up_cache_map;
    }

    synchronized  static HashMap<String, UniProtTaxonomy> getCodeTaxCacheMap() {
        return _code_up_cache_map;
    }

    synchronized  static HashMap<String, UniProtTaxonomy> getIdTaxCacheMap() {
        return _id_up_cache_map;
    }

    synchronized  static HashMap<String, UniProtTaxonomy> getSnTaxCacheMap() {
        return _sn_up_cache_map;
    }

    synchronized  static HashMap<String, UniProtTaxonomy> getLineageTaxCacheMap() {
        return _lineage_up_cache_map;
    }

    
    enum QUERY_TYPE {
        CODE, SN, CN, ID, LIN;
    }
    
    @Override
    public void run() {
        execute();
    }
}