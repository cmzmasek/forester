// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.rio;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.forester.datastructures.IntMatrix;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sdi.GSDI;
import org.forester.sdi.GSDIR;
import org.forester.sdi.SDIException;
import org.forester.sdi.SDIR;
import org.forester.sdi.SDIutil.ALGORITHM;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class RIO {

    public enum REROOTING {
        NONE, BY_ALGORITHM, MIDPOINT, OUTGROUP;
    }
    private Phylogeny[]                _analyzed_gene_trees;
    private List<PhylogenyNode>        _removed_gene_tree_nodes;
    private int                        _ext_nodes;
    private TaxonomyComparisonBase     _gsdir_tax_comp_base;
    private StringBuilder              _log;
    private BasicDescriptiveStatistics _duplications_stats;
    private boolean                    _produce_log;
    private boolean                    _verbose;
    private REROOTING                  _rerooting;

    public RIO( final File gene_trees_file,
                final Phylogeny species_tree,
                final ALGORITHM algorithm,
                final REROOTING rerooting,
                final String outgroup,
                final boolean produce_log,
                final boolean verbose ) throws IOException, SDIException, RIOException {
        checkReRooting( rerooting, outgroup );
        init( produce_log, verbose, rerooting );
        inferOrthologs( gene_trees_file, species_tree, algorithm, outgroup );
    }

    public RIO( final Phylogeny[] gene_trees,
                final Phylogeny species_tree,
                final ALGORITHM algorithm,
                final REROOTING rerooting,
                final String outgroup,
                final boolean produce_log,
                final boolean verbose ) throws IOException, SDIException, RIOException {
        checkReRooting( rerooting, outgroup );
        init( produce_log, verbose, rerooting );
        inferOrthologs( gene_trees, species_tree, algorithm, outgroup );
    }

    public final Phylogeny[] getAnalyzedGeneTrees() {
        return _analyzed_gene_trees;
    }

    /**
     * Returns the numbers of number of ext nodes in gene trees analyzed (after
     * stripping).
     * 
     * @return number of ext nodes in gene trees analyzed (after stripping)
     */
    public final int getExtNodesOfAnalyzedGeneTrees() {
        return _ext_nodes;
    }

    public final TaxonomyComparisonBase getGSDIRtaxCompBase() {
        return _gsdir_tax_comp_base;
    }

    public final StringBuilder getLog() {
        return _log;
    }

    public final List<PhylogenyNode> getRemovedGeneTreeNodes() {
        return _removed_gene_tree_nodes;
    }

    private final void inferOrthologs( final File gene_trees_file,
                                       final Phylogeny species_tree,
                                       final ALGORITHM algorithm,
                                       final String outgroup ) throws SDIException, RIOException,
            FileNotFoundException, IOException {
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final PhylogenyParser p = ParserUtils.createParserDependingOnFileType( gene_trees_file, true );
        if ( p instanceof NHXParser ) {
            final NHXParser nhx = ( NHXParser ) p;
            nhx.setReplaceUnderscores( false );
            nhx.setIgnoreQuotes( true );
            nhx.setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.YES );
        }
        final Phylogeny[] gene_trees = factory.create( gene_trees_file, p );
        inferOrthologs( gene_trees, species_tree, algorithm, outgroup );
    }

    private final void inferOrthologs( final Phylogeny[] gene_trees,
                                       final Phylogeny species_tree,
                                       final ALGORITHM algorithm,
                                       final String outgroup ) throws SDIException, RIOException,
            FileNotFoundException, IOException {
        if ( algorithm == ALGORITHM.SDIR ) {
            // Removes from species_tree all species not found in gene_tree.
            PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( gene_trees[ 0 ], species_tree );
            if ( species_tree.isEmpty() ) {
                throw new RIOException( "failed to establish species based mapping between gene and species trees" );
            }
        }
        if ( _produce_log ) {
            _log = new StringBuilder();
            if ( _rerooting == REROOTING.BY_ALGORITHM ) {
                writeLogSubHeader();
            }
        }
        _analyzed_gene_trees = new Phylogeny[ gene_trees.length ];
        int gene_tree_ext_nodes = 0;
        if ( _verbose ) {
            System.out.println();
        }
        for( int i = 0; i < gene_trees.length; ++i ) {
            final Phylogeny gt = gene_trees[ i ];
            if ( _verbose ) {
                ForesterUtil.updateProgress( ( double ) i / gene_trees.length );
            }
            if ( i == 0 ) {
                gene_tree_ext_nodes = gt.getNumberOfExternalNodes();
            }
            else if ( gene_tree_ext_nodes != gt.getNumberOfExternalNodes() ) {
                throw new RIOException( "gene tree #" + ( i + 1 ) + " has a different number of external nodes ("
                        + gt.getNumberOfExternalNodes() + ") than the preceding gene trees (" + gene_tree_ext_nodes
                        + ")" );
            }
            if ( algorithm == ALGORITHM.SDIR ) {
                // Removes from gene_tree all species not found in species_tree.
                PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gt );
                if ( gt.isEmpty() ) {
                    throw new RIOException( "failed to establish species based mapping between gene and species trees" );
                }
            }
            _analyzed_gene_trees[ i ] = performOrthologInference( gt, species_tree, algorithm, outgroup, i );
        }
        if ( _verbose ) {
            System.out.println();
            System.out.println();
        }
    }

    private final void init( final boolean produce_log, final boolean verbose, final REROOTING rerooting ) {
        _produce_log = produce_log;
        _verbose = verbose;
        _rerooting = rerooting;
        _ext_nodes = -1;
        _log = null;
        _gsdir_tax_comp_base = null;
        _analyzed_gene_trees = null;
        _removed_gene_tree_nodes = null;
        _duplications_stats = new BasicDescriptiveStatistics();
    }

    private final Phylogeny performOrthologInference( final Phylogeny gene_tree,
                                                      final Phylogeny species_tree,
                                                      final ALGORITHM algorithm,
                                                      final String outgroup,
                                                      final int i ) throws SDIException, RIOException {
        final Phylogeny assigned_tree;
        switch ( algorithm ) {
            case SDIR: {
                assigned_tree = performOrthologInferenceBySDI( gene_tree, species_tree );
                break;
            }
            case GSDIR: {
                assigned_tree = performOrthologInferenceByGSDI( gene_tree, species_tree, outgroup, i );
                break;
            }
            default: {
                throw new IllegalArgumentException( "illegal algorithm: " + algorithm );
            }
        }
        if ( i == 0 ) {
            _ext_nodes = assigned_tree.getNumberOfExternalNodes();
        }
        else if ( _ext_nodes != assigned_tree.getNumberOfExternalNodes() ) {
            throw new RIOException( "after stripping gene tree #" + ( i + 1 )
                    + " has a different number of external nodes (" + assigned_tree.getNumberOfExternalNodes()
                    + ") than the preceding gene trees (" + _ext_nodes + ")" );
        }
        return assigned_tree;
    }

    private final Phylogeny performOrthologInferenceByGSDI( final Phylogeny gene_tree,
                                                            final Phylogeny species_tree,
                                                            final String outgroup,
                                                            final int i ) throws SDIException, RIOException {
        final Phylogeny assigned_tree;
        if ( _rerooting == REROOTING.BY_ALGORITHM ) {
            final GSDIR gsdir = new GSDIR( gene_tree, species_tree, true, i == 0 );
            final List<Phylogeny> assigned_trees = gsdir.getMinDuplicationsSumGeneTrees();
            if ( i == 0 ) {
                _removed_gene_tree_nodes = gsdir.getStrippedExternalGeneTreeNodes();
                for( final PhylogenyNode r : _removed_gene_tree_nodes ) {
                    if ( !r.getNodeData().isHasTaxonomy() ) {
                        throw new RIOException( "node with no (appropriate) taxonomic information found in gene tree #1: "
                                + r.toString() );
                    }
                }
            }
            final List<Integer> shortests = GSDIR.getIndexesOfShortestTree( assigned_trees );
            assigned_tree = assigned_trees.get( shortests.get( 0 ) );
            if ( _produce_log ) {
                writeStatsToLog( i, gsdir, shortests );
            }
            if ( i == 0 ) {
                _gsdir_tax_comp_base = gsdir.getTaxCompBase();
            }
            _duplications_stats.addValue( gsdir.getMinDuplicationsSum() );
        }
        else {
            if ( _rerooting == REROOTING.MIDPOINT ) {
                PhylogenyMethods.midpointRoot( gene_tree );
            }
            else if ( _rerooting == REROOTING.OUTGROUP ) {
                PhylogenyNode n;
                try {
                    n = gene_tree.getNode( outgroup );
                }
                catch ( IllegalArgumentException e ) {
                    throw new RIOException( "failed to perform re-rooting by outgroup: " + e.getLocalizedMessage() );
                }
                gene_tree.reRoot( n );
            }
            final GSDI gsdi = new GSDI( gene_tree, species_tree, true, true, true );
            _removed_gene_tree_nodes = gsdi.getStrippedExternalGeneTreeNodes();
            for( final PhylogenyNode r : _removed_gene_tree_nodes ) {
                if ( !r.getNodeData().isHasTaxonomy() ) {
                    throw new RIOException( "node with no (appropriate) taxonomic information found in gene tree #1: "
                            + r.toString() );
                }
            }
            assigned_tree = gene_tree;
            if ( i == 0 ) {
                _gsdir_tax_comp_base = gsdi.getTaxCompBase();
            }
            _duplications_stats.addValue( gsdi.getDuplicationsSum() );
        }
        return assigned_tree;
    }

    private final Phylogeny performOrthologInferenceBySDI( final Phylogeny gene_tree, final Phylogeny species_tree )
            throws SDIException {
        final SDIR sdir = new SDIR();
        return sdir.infer( gene_tree, species_tree, false, true, true, true, 1 )[ 0 ];
    }

    private void writeLogSubHeader() {
        _log.append( "#" );
        _log.append( "\t" );
        _log.append( "with minimal number of duplications" );
        _log.append( "/" );
        _log.append( "root placements" );
        _log.append( "\t[" );
        _log.append( "min" );
        _log.append( "-" );
        _log.append( "max" );
        _log.append( "]\t<" );
        _log.append( "shortest" );
        _log.append( ">" );
        _log.append( ForesterUtil.LINE_SEPARATOR );
    }

    private final void writeStatsToLog( final int i, final GSDIR gsdir, final List<Integer> shortests ) {
        final BasicDescriptiveStatistics stats = gsdir.getDuplicationsSumStats();
        _log.append( i );
        _log.append( "\t" );
        _log.append( gsdir.getMinDuplicationsSumGeneTrees().size() );
        _log.append( "/" );
        _log.append( stats.getN() );
        _log.append( "\t[" );
        _log.append( ( int ) stats.getMin() );
        _log.append( "-" );
        _log.append( ( int ) stats.getMax() );
        _log.append( "]\t<" );
        _log.append( shortests.size() );
        _log.append( ">" );
        _log.append( ForesterUtil.LINE_SEPARATOR );
    }

    public final static IntMatrix calculateOrthologTable( final Phylogeny[] analyzed_gene_trees, final boolean sort )
            throws RIOException {
        final List<String> labels = new ArrayList<String>();
        final Set<String> labels_set = new HashSet<String>();
        String label;
        for( final PhylogenyNode n : analyzed_gene_trees[ 0 ].getExternalNodes() ) {
            if ( n.getNodeData().isHasSequence() && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getName() ) ) {
                label = n.getNodeData().getSequence().getName();
            }
            else if ( n.getNodeData().isHasSequence()
                    && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getSymbol() ) ) {
                label = n.getNodeData().getSequence().getSymbol();
            }
            else if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                label = n.getName();
            }
            else {
                throw new RIOException( "node " + n + " has no appropriate label" );
            }
            if ( labels_set.contains( label ) ) {
                throw new RIOException( "label " + label + " is not unique" );
            }
            labels_set.add( label );
            labels.add( label );
        }
        if ( sort ) {
            Collections.sort( labels );
        }
        final IntMatrix m = new IntMatrix( labels );
        int counter = 0;
        for( final Phylogeny gt : analyzed_gene_trees ) {
            counter++;
            PhylogenyMethods.preOrderReId( gt );
            final HashMap<String, PhylogenyNode> map = PhylogenyMethods.createNameToExtNodeMap( gt );
            for( int x = 0; x < m.size(); ++x ) {
                final String mx = m.getLabel( x );
                final PhylogenyNode nx = map.get( mx );
                if ( nx == null ) {
                    throw new RIOException( "node \"" + mx + "\" not present in gene tree #" + counter );
                }
                String my;
                PhylogenyNode ny;
                for( int y = 0; y < m.size(); ++y ) {
                    my = m.getLabel( y );
                    ny = map.get( my );
                    if ( ny == null ) {
                        throw new RIOException( "node \"" + my + "\" not present in gene tree #" + counter );
                    }
                    if ( !PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( nx, ny ).isDuplication() ) {
                        m.inreaseByOne( x, y );
                    }
                }
            }
        }
        return m;
    }

    private final static void checkReRooting( final REROOTING rerooting, final String outgroup ) {
        if ( !ForesterUtil.isEmpty( outgroup ) && rerooting != REROOTING.OUTGROUP ) {
            throw new IllegalArgumentException( "can only use outgroup when re-rooting in this manner" );
        }
    }

    public BasicDescriptiveStatistics getDuplicationsStatistics() {
        return _duplications_stats;
    }
}
