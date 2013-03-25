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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.rio;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.datastructures.IntMatrix;
import org.forester.io.parsers.IteratingPhylogenyParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sdi.GSDI;
import org.forester.sdi.GSDIR;
import org.forester.sdi.SDIException;
import org.forester.sdi.SDIR;
import org.forester.sdi.SDIutil;
import org.forester.sdi.SDIutil.ALGORITHM;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class RIO {

    public static final int                  DEFAULT_RANGE = -1;
    private static final int                 END_OF_GT     = Integer.MAX_VALUE;
    private static IntMatrix                 _m;
    private Phylogeny[]                      _analyzed_gene_trees;
    private List<PhylogenyNode>              _removed_gene_tree_nodes;
    private int                              _ext_nodes;
    private int                              _int_nodes;
    private TaxonomyComparisonBase           _gsdir_tax_comp_base;
    private final StringBuilder              _log;
    private final BasicDescriptiveStatistics _duplications_stats;
    private final boolean                    _produce_log;
    private final boolean                    _verbose;
    private final REROOTING                  _rerooting;
    private final Phylogeny                  _species_tree;
    private Phylogeny                        _min_dub_gene_tree;

    private RIO( final IteratingPhylogenyParser p,
                 final Phylogeny species_tree,
                 final ALGORITHM algorithm,
                 final REROOTING rerooting,
                 final String outgroup,
                 int first,
                 int last,
                 final boolean produce_log,
                 final boolean verbose,
                 final boolean transfer_taxonomy ) throws IOException, SDIException, RIOException {
        if ( ( last == DEFAULT_RANGE ) && ( first >= 0 ) ) {
            last = END_OF_GT;
        }
        else if ( ( first == DEFAULT_RANGE ) && ( last >= 0 ) ) {
            first = 0;
        }
        removeSingleDescendentsNodes( species_tree, verbose );
        p.reset();
        checkPreconditions( p, species_tree, rerooting, outgroup, first, last );
        _produce_log = produce_log;
        _verbose = verbose;
        _rerooting = rerooting;
        _ext_nodes = -1;
        _int_nodes = -1;
        _log = new StringBuilder();
        _gsdir_tax_comp_base = null;
        _analyzed_gene_trees = null;
        _removed_gene_tree_nodes = null;
        _duplications_stats = new BasicDescriptiveStatistics();
        p.reset();
        inferOrthologs( p, species_tree, algorithm, outgroup, first, last, transfer_taxonomy );
        _species_tree = species_tree;
    }

    private RIO( final Phylogeny[] gene_trees,
                 final Phylogeny species_tree,
                 final ALGORITHM algorithm,
                 final REROOTING rerooting,
                 final String outgroup,
                 int first,
                 int last,
                 final boolean produce_log,
                 final boolean verbose,
                 final boolean transfer_taxonomy ) throws IOException, SDIException, RIOException {
        if ( ( last == DEFAULT_RANGE ) && ( first >= 0 ) ) {
            last = gene_trees.length - 1;
        }
        else if ( ( first == DEFAULT_RANGE ) && ( last >= 0 ) ) {
            first = 0;
        }
        removeSingleDescendentsNodes( species_tree, verbose );
        checkPreconditions( gene_trees, species_tree, rerooting, outgroup, first, last );
        _produce_log = produce_log;
        _verbose = verbose;
        _rerooting = rerooting;
        _ext_nodes = -1;
        _int_nodes = -1;
        _log = new StringBuilder();
        _gsdir_tax_comp_base = null;
        _analyzed_gene_trees = null;
        _removed_gene_tree_nodes = null;
        _duplications_stats = new BasicDescriptiveStatistics();
        inferOrthologs( gene_trees, species_tree, algorithm, outgroup, first, last, transfer_taxonomy );
        _species_tree = species_tree;
    }

    public final Phylogeny[] getAnalyzedGeneTrees() {
        return _analyzed_gene_trees;
    }

    public final BasicDescriptiveStatistics getDuplicationsStatistics() {
        return _duplications_stats;
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

    /**
     * Returns the numbers of number of int nodes in gene trees analyzed (after
     * stripping).
     * 
     * @return number of int nodes in gene trees analyzed (after stripping)
     */
    public final int getIntNodesOfAnalyzedGeneTrees() {
        return _int_nodes;
    }

    public final StringBuilder getLog() {
        return _log;
    }

    final public Phylogeny getMinDuplicationsGeneTree() {
        return _min_dub_gene_tree;
    }

    public final IntMatrix getOrthologTable() {
        return _m;
    }

    public final List<PhylogenyNode> getRemovedGeneTreeNodes() {
        return _removed_gene_tree_nodes;
    }

    public final Phylogeny getSpeciesTree() {
        return _species_tree;
    }

    private final void inferOrthologs( final IteratingPhylogenyParser parser,
                                       final Phylogeny species_tree,
                                       final ALGORITHM algorithm,
                                       final String outgroup,
                                       int first,
                                       final int last,
                                       final boolean transfer_taxonomy ) throws SDIException, RIOException,
            FileNotFoundException, IOException {
        if ( !parser.hasNext() ) {
            throw new RIOException( "no gene trees to analyze" );
        }
        if ( log() ) {
            preLog( -1, species_tree, algorithm, outgroup );
        }
        if ( _verbose ) {
            System.out.println();
        }
        final DecimalFormat pf = new java.text.DecimalFormat( "000" );
        int gene_tree_ext_nodes = 0;
        int i = 0;
        int counter = 0;
        final boolean no_range = ( first < 0 ) || ( last < first );
        while ( parser.hasNext() ) {
            final Phylogeny gt = parser.next();
            if ( no_range || ( ( i >= first ) && ( i <= last ) ) ) {
                if ( gt.isEmpty() ) {
                    throw new RIOException( "gene tree #" + i + " is empty" );
                }
                if ( gt.getNumberOfExternalNodes() == 1 ) {
                    throw new RIOException( "gene tree #" + i + " has only one external node" );
                }
                if ( _verbose ) {
                    ForesterUtil.updateProgress( i, pf );
                }
                if ( counter == 0 ) {
                    if ( algorithm == ALGORITHM.SDIR ) {
                        // Removes from species_tree all species not found in gene_tree.
                        PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( gt, species_tree );
                        if ( species_tree.isEmpty() ) {
                            throw new RIOException( "failed to establish species based mapping between gene and species trees" );
                        }
                    }
                    gene_tree_ext_nodes = gt.getNumberOfExternalNodes();
                }
                else if ( gene_tree_ext_nodes != gt.getNumberOfExternalNodes() ) {
                    throw new RIOException( "gene tree #" + i + " has a different number of external nodes ("
                            + gt.getNumberOfExternalNodes() + ") than the preceding gene tree(s) ("
                            + gene_tree_ext_nodes + ")" );
                }
                if ( algorithm == ALGORITHM.SDIR ) {
                    // Removes from gene_tree all species not found in species_tree.
                    PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gt );
                    if ( gt.isEmpty() ) {
                        throw new RIOException( "failed to establish species based mapping between gene and species trees" );
                    }
                }
                final Phylogeny analyzed_gt = performOrthologInference( gt,
                                                                        species_tree,
                                                                        algorithm,
                                                                        outgroup,
                                                                        counter,
                                                                        transfer_taxonomy );
                RIO.calculateOrthologTable( analyzed_gt, true, counter );
                ++counter;
            }
            ++i;
        }
        if ( ( first >= 0 ) && ( counter == 0 ) && ( i > 0 ) ) {
            throw new RIOException( "attempt to analyze first gene tree #" + first + " in a set of " + i );
        }
        if ( no_range ) {
            first = 0;
        }
        if ( log() ) {
            postLog( species_tree, first, first + counter - 1 );
        }
        if ( _verbose ) {
            System.out.println();
            System.out.println();
        }
    }

    private final void inferOrthologs( final Phylogeny[] gene_trees,
                                       final Phylogeny species_tree,
                                       final ALGORITHM algorithm,
                                       final String outgroup,
                                       final int first,
                                       final int last,
                                       final boolean transfer_taxonomy ) throws SDIException, RIOException,
            FileNotFoundException, IOException {
        if ( algorithm == ALGORITHM.SDIR ) {
            // Removes from species_tree all species not found in gene_tree.
            PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( gene_trees[ 0 ], species_tree );
            if ( species_tree.isEmpty() ) {
                throw new RIOException( "failed to establish species based mapping between gene and species trees" );
            }
        }
        final Phylogeny[] my_gene_trees;
        if ( ( first >= 0 ) && ( last >= first ) && ( last < gene_trees.length ) ) {
            my_gene_trees = new Phylogeny[ ( 1 + last ) - first ];
            int c = 0;
            for( int i = first; i <= last; ++i ) {
                my_gene_trees[ c++ ] = gene_trees[ i ];
            }
        }
        else {
            my_gene_trees = gene_trees;
        }
        if ( log() ) {
            preLog( gene_trees.length, species_tree, algorithm, outgroup );
        }
        if ( _verbose && ( my_gene_trees.length >= 4 ) ) {
            System.out.println();
        }
        _analyzed_gene_trees = new Phylogeny[ my_gene_trees.length ];
        int gene_tree_ext_nodes = 0;
        for( int i = 0; i < my_gene_trees.length; ++i ) {
            final Phylogeny gt = my_gene_trees[ i ];
            if ( gt.isEmpty() ) {
                throw new RIOException( "gene tree #" + i + " is empty" );
            }
            if ( gt.getNumberOfExternalNodes() == 1 ) {
                throw new RIOException( "gene tree #" + i + " has only one external node" );
            }
            if ( _verbose && ( my_gene_trees.length > 4 ) ) {
                ForesterUtil.updateProgress( ( ( double ) i ) / my_gene_trees.length );
            }
            if ( i == 0 ) {
                gene_tree_ext_nodes = gt.getNumberOfExternalNodes();
            }
            else if ( gene_tree_ext_nodes != gt.getNumberOfExternalNodes() ) {
                throw new RIOException( "gene tree #" + i + " has a different number of external nodes ("
                        + gt.getNumberOfExternalNodes() + ") than the preceding gene tree(s) (" + gene_tree_ext_nodes
                        + ")" );
            }
            if ( algorithm == ALGORITHM.SDIR ) {
                // Removes from gene_tree all species not found in species_tree.
                PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gt );
                if ( gt.isEmpty() ) {
                    throw new RIOException( "failed to establish species based mapping between gene and species trees" );
                }
            }
            _analyzed_gene_trees[ i ] = performOrthologInference( gt,
                                                                  species_tree,
                                                                  algorithm,
                                                                  outgroup,
                                                                  i,
                                                                  transfer_taxonomy );
        }
        if ( log() ) {
            postLog( species_tree, first, last );
        }
        if ( _verbose && ( my_gene_trees.length > 4 ) ) {
            System.out.println();
            System.out.println();
        }
    }

    private final boolean log() {
        return _produce_log;
    }

    private final void log( final String s ) {
        _log.append( s );
        _log.append( ForesterUtil.LINE_SEPARATOR );
    }

    private final void logRemovedGeneTreeNodes() {
        log( "Species stripped from gene trees:" );
        final SortedSet<String> rn = new TreeSet<String>();
        for( final PhylogenyNode n : getRemovedGeneTreeNodes() ) {
            final Taxonomy t = n.getNodeData().getTaxonomy();
            switch ( getGSDIRtaxCompBase() ) {
                case CODE: {
                    rn.add( t.getTaxonomyCode() );
                    break;
                }
                case ID: {
                    rn.add( t.getIdentifier().toString() );
                    break;
                }
                case SCIENTIFIC_NAME: {
                    rn.add( t.getScientificName() );
                    break;
                }
            }
        }
        for( final String s : rn ) {
            log( s );
        }
        log( "" );
    }

    private final Phylogeny performOrthologInference( final Phylogeny gene_tree,
                                                      final Phylogeny species_tree,
                                                      final ALGORITHM algorithm,
                                                      final String outgroup,
                                                      final int i,
                                                      final boolean transfer_taxonomy ) throws SDIException,
            RIOException {
        final Phylogeny assigned_tree;
        switch ( algorithm ) {
            case SDIR: {
                assigned_tree = performOrthologInferenceBySDI( gene_tree, species_tree );
                break;
            }
            case GSDIR: {
                assigned_tree = performOrthologInferenceByGSDI( gene_tree, species_tree, outgroup, i, transfer_taxonomy );
                break;
            }
            default: {
                throw new IllegalArgumentException( "illegal algorithm: " + algorithm );
            }
        }
        if ( i == 0 ) {
            _ext_nodes = assigned_tree.getNumberOfExternalNodes();
            _int_nodes = assigned_tree.getNumberOfInternalNodes();
        }
        else if ( _ext_nodes != assigned_tree.getNumberOfExternalNodes() ) {
            throw new RIOException( "after stripping gene tree #" + i + " has a different number of external nodes ("
                    + assigned_tree.getNumberOfExternalNodes() + ") than the preceding gene tree(s) (" + _ext_nodes
                    + ")" );
        }
        return assigned_tree;
    }

    private final Phylogeny performOrthologInferenceByGSDI( final Phylogeny gene_tree,
                                                            final Phylogeny species_tree,
                                                            final String outgroup,
                                                            final int i,
                                                            final boolean transfer_taxonomy ) throws SDIException,
            RIOException {
        final Phylogeny assigned_tree;
        final int dups;
        if ( _rerooting == REROOTING.BY_ALGORITHM ) {
            final GSDIR gsdir = new GSDIR( gene_tree, species_tree, true, i == 0, transfer_taxonomy );
            assigned_tree = gsdir.getMinDuplicationsSumGeneTree();
            if ( i == 0 ) {
                _removed_gene_tree_nodes = gsdir.getStrippedExternalGeneTreeNodes();
                for( final PhylogenyNode r : _removed_gene_tree_nodes ) {
                    if ( !r.getNodeData().isHasTaxonomy() ) {
                        throw new RIOException( "node with no (appropriate) taxonomic information found in gene tree #"
                                + i + ": " + r.toString() );
                    }
                }
            }
            if ( i == 0 ) {
                _gsdir_tax_comp_base = gsdir.getTaxCompBase();
            }
            dups = gsdir.getMinDuplicationsSum();
        }
        else {
            if ( _rerooting == REROOTING.MIDPOINT ) {
                PhylogenyMethods.midpointRoot( gene_tree );
            }
            else if ( _rerooting == REROOTING.OUTGROUP ) {
                final PhylogenyNode n = gene_tree.getNode( outgroup );
                gene_tree.reRoot( n );
            }
            final GSDI gsdi = new GSDI( gene_tree, species_tree, true, true, true, transfer_taxonomy );
            _removed_gene_tree_nodes = gsdi.getStrippedExternalGeneTreeNodes();
            for( final PhylogenyNode r : _removed_gene_tree_nodes ) {
                if ( !r.getNodeData().isHasTaxonomy() ) {
                    throw new RIOException( "node with no (appropriate) taxonomic information found in gene tree #" + i
                            + ": " + r.toString() );
                }
            }
            assigned_tree = gene_tree;
            if ( i == 0 ) {
                _gsdir_tax_comp_base = gsdi.getTaxCompBase();
            }
            dups = gsdi.getDuplicationsSum();
        }
        if ( ( i == 0 ) || ( dups < _duplications_stats.getMin() ) ) {
            _min_dub_gene_tree = assigned_tree;
        }
        _duplications_stats.addValue( dups );
        return assigned_tree;
    }

    private final Phylogeny performOrthologInferenceBySDI( final Phylogeny gene_tree, final Phylogeny species_tree )
            throws SDIException {
        final SDIR sdir = new SDIR();
        return sdir.infer( gene_tree, species_tree, false, true, true, true, 1 )[ 0 ];
    }

    private final void postLog( final Phylogeny species_tree, final int first, final int last ) {
        log( "" );
        if ( ( getRemovedGeneTreeNodes() != null ) && ( getRemovedGeneTreeNodes().size() > 0 ) ) {
            logRemovedGeneTreeNodes();
        }
        log( "Species tree external nodes (after stripping)   : " + species_tree.getNumberOfExternalNodes() );
        log( "Species tree polytomies (after stripping)       : "
                + PhylogenyMethods.countNumberOfPolytomies( species_tree ) );
        log( "Taxonomy linking based on                       : " + getGSDIRtaxCompBase() );
        final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.#" );
        if ( ( first >= 0 ) && ( last >= 0 ) ) {
            log( "Gene trees analyzed range                       : " + first + "-" + last );
        }
        log( "Gene trees analyzed                             : " + _duplications_stats.getN() );
        log( "Mean number of duplications                     : " + df.format( _duplications_stats.arithmeticMean() )
                + " (sd: " + df.format( _duplications_stats.sampleStandardDeviation() ) + ")" + " ("
                + df.format( ( 100.0 * _duplications_stats.arithmeticMean() ) / getIntNodesOfAnalyzedGeneTrees() )
                + "%)" );
        if ( _duplications_stats.getN() > 3 ) {
            log( "Median number of duplications                   : " + df.format( _duplications_stats.median() )
                    + " (" + df.format( ( 100.0 * _duplications_stats.median() ) / getIntNodesOfAnalyzedGeneTrees() )
                    + "%)" );
        }
        log( "Minimum duplications                            : " + ( int ) _duplications_stats.getMin() + " ("
                + df.format( ( 100.0 * _duplications_stats.getMin() ) / getIntNodesOfAnalyzedGeneTrees() ) + "%)" );
        log( "Maximum duplications                            : " + ( int ) _duplications_stats.getMax() + " ("
                + df.format( ( 100.0 * _duplications_stats.getMax() ) / getIntNodesOfAnalyzedGeneTrees() ) + "%)" );
        log( "Gene tree internal nodes                        : " + getIntNodesOfAnalyzedGeneTrees() );
        log( "Gene tree external nodes                        : " + getExtNodesOfAnalyzedGeneTrees() );
    }

    private final void preLog( final int gene_trees,
                               final Phylogeny species_tree,
                               final ALGORITHM algorithm,
                               final String outgroup ) {
        if ( gene_trees > 0 ) {
            log( "Number of gene trees (total)                    : " + gene_trees );
        }
        log( "Algorithm                                       : " + algorithm );
        log( "Species tree external nodes (prior to stripping): " + species_tree.getNumberOfExternalNodes() );
        log( "Species tree polytomies (prior to stripping)    : "
                + PhylogenyMethods.countNumberOfPolytomies( species_tree ) );
        String rs = "";
        switch ( _rerooting ) {
            case BY_ALGORITHM: {
                rs = "minimizing duplications";
                break;
            }
            case MIDPOINT: {
                rs = "midpoint";
                break;
            }
            case OUTGROUP: {
                rs = "outgroup: " + outgroup;
                break;
            }
            case NONE: {
                rs = "none";
                break;
            }
        }
        log( "Re-rooting                                      : " + rs );
    }

    public final static IntMatrix calculateOrthologTable( final Phylogeny[] analyzed_gene_trees, final boolean sort )
            throws RIOException {
        final List<String> labels = new ArrayList<String>();
        final Set<String> labels_set = new HashSet<String>();
        for( final PhylogenyNode n : analyzed_gene_trees[ 0 ].getExternalNodes() ) {
            final String label = obtainLabel( labels_set, n );
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
            updateCounts( m, counter, gt );
        }
        return m;
    }

    public final static RIO executeAnalysis( final File gene_trees_file,
                                             final File species_tree_file,
                                             final ALGORITHM algorithm,
                                             final REROOTING rerooting,
                                             final String outgroup,
                                             final int first,
                                             final int last,
                                             final boolean produce_log,
                                             final boolean verbose,
                                             final boolean transfer_taxonomy ) throws IOException, SDIException,
            RIOException {
        final Phylogeny[] gene_trees = parseGeneTrees( gene_trees_file );
        if ( gene_trees.length < 1 ) {
            throw new RIOException( "\"" + gene_trees_file + "\" is devoid of appropriate gene trees" );
        }
        final Phylogeny species_tree = SDIutil.parseSpeciesTree( gene_trees[ 0 ],
                                                                 species_tree_file,
                                                                 false,
                                                                 true,
                                                                 TAXONOMY_EXTRACTION.NO );
        return new RIO( gene_trees,
                        species_tree,
                        algorithm,
                        rerooting,
                        outgroup,
                        first,
                        last,
                        produce_log,
                        verbose,
                        transfer_taxonomy );
    }

    public final static RIO executeAnalysis( final File gene_trees_file,
                                             final Phylogeny species_tree,
                                             final ALGORITHM algorithm,
                                             final REROOTING rerooting,
                                             final String outgroup,
                                             final boolean produce_log,
                                             final boolean verbose,
                                             final boolean transfer_taxonomy ) throws IOException, SDIException,
            RIOException {
        return new RIO( parseGeneTrees( gene_trees_file ),
                        species_tree,
                        algorithm,
                        rerooting,
                        outgroup,
                        DEFAULT_RANGE,
                        DEFAULT_RANGE,
                        produce_log,
                        verbose,
                        transfer_taxonomy );
    }

    public final static RIO executeAnalysis( final File gene_trees_file,
                                             final Phylogeny species_tree,
                                             final ALGORITHM algorithm,
                                             final REROOTING rerooting,
                                             final String outgroup,
                                             final int first,
                                             final int last,
                                             final boolean produce_log,
                                             final boolean verbose,
                                             final boolean transfer_taxonomy ) throws IOException, SDIException,
            RIOException {
        return new RIO( parseGeneTrees( gene_trees_file ),
                        species_tree,
                        algorithm,
                        rerooting,
                        outgroup,
                        first,
                        last,
                        produce_log,
                        verbose,
                        transfer_taxonomy );
    }

    public final static RIO executeAnalysis( final IteratingPhylogenyParser p,
                                             final File species_tree_file,
                                             final ALGORITHM algorithm,
                                             final REROOTING rerooting,
                                             final String outgroup,
                                             final int first,
                                             final int last,
                                             final boolean produce_log,
                                             final boolean verbose,
                                             final boolean transfer_taxonomy ) throws IOException, SDIException,
            RIOException {
        final Phylogeny g0 = p.next();
        if ( ( g0 == null ) || g0.isEmpty() || ( g0.getNumberOfExternalNodes() < 2 ) ) {
            throw new RIOException( "input file does not seem to contain any gene trees" );
        }
        final Phylogeny species_tree = SDIutil.parseSpeciesTree( g0,
                                                                 species_tree_file,
                                                                 false,
                                                                 true,
                                                                 TAXONOMY_EXTRACTION.NO );
        p.reset();
        return new RIO( p,
                        species_tree,
                        algorithm,
                        rerooting,
                        outgroup,
                        first,
                        last,
                        produce_log,
                        verbose,
                        transfer_taxonomy );
    }

    public final static RIO executeAnalysis( final IteratingPhylogenyParser p,
                                             final Phylogeny species_tree,
                                             final ALGORITHM algorithm,
                                             final REROOTING rerooting,
                                             final String outgroup,
                                             final boolean produce_log,
                                             final boolean verbose,
                                             final boolean transfer_taxonomy ) throws IOException, SDIException,
            RIOException {
        return new RIO( p,
                        species_tree,
                        algorithm,
                        rerooting,
                        outgroup,
                        DEFAULT_RANGE,
                        DEFAULT_RANGE,
                        produce_log,
                        verbose,
                        transfer_taxonomy );
    }

    public final static RIO executeAnalysis( final IteratingPhylogenyParser p,
                                             final Phylogeny species_tree,
                                             final ALGORITHM algorithm,
                                             final REROOTING rerooting,
                                             final String outgroup,
                                             final int first,
                                             final int last,
                                             final boolean produce_log,
                                             final boolean verbose,
                                             final boolean transfer_taxonomy ) throws IOException, SDIException,
            RIOException {
        return new RIO( p,
                        species_tree,
                        algorithm,
                        rerooting,
                        outgroup,
                        first,
                        last,
                        produce_log,
                        verbose,
                        transfer_taxonomy );
    }

    public final static RIO executeAnalysis( final Phylogeny[] gene_trees, final Phylogeny species_tree )
            throws IOException, SDIException, RIOException {
        return new RIO( gene_trees,
                        species_tree,
                        ALGORITHM.GSDIR,
                        REROOTING.BY_ALGORITHM,
                        null,
                        DEFAULT_RANGE,
                        DEFAULT_RANGE,
                        false,
                        false,
                        false );
    }

    public final static RIO executeAnalysis( final Phylogeny[] gene_trees,
                                             final Phylogeny species_tree,
                                             final ALGORITHM algorithm,
                                             final REROOTING rerooting,
                                             final String outgroup,
                                             final boolean produce_log,
                                             final boolean verbose,
                                             final boolean transfer_taxonomy ) throws IOException, SDIException,
            RIOException {
        return new RIO( gene_trees,
                        species_tree,
                        algorithm,
                        rerooting,
                        outgroup,
                        DEFAULT_RANGE,
                        DEFAULT_RANGE,
                        produce_log,
                        verbose,
                        transfer_taxonomy );
    }

    public final static RIO executeAnalysis( final Phylogeny[] gene_trees,
                                             final Phylogeny species_tree,
                                             final ALGORITHM algorithm,
                                             final REROOTING rerooting,
                                             final String outgroup,
                                             final int first,
                                             final int last,
                                             final boolean produce_log,
                                             final boolean verbose,
                                             final boolean transfer_taxonomy ) throws IOException, SDIException,
            RIOException {
        return new RIO( gene_trees,
                        species_tree,
                        algorithm,
                        rerooting,
                        outgroup,
                        first,
                        last,
                        produce_log,
                        verbose,
                        transfer_taxonomy );
    }

    private final static void calculateOrthologTable( final Phylogeny g, final boolean sort, final int counter )
            throws RIOException {
        if ( counter == 0 ) {
            final List<String> labels = new ArrayList<String>();
            final Set<String> labels_set = new HashSet<String>();
            for( final PhylogenyNode n : g.getExternalNodes() ) {
                final String label = obtainLabel( labels_set, n );
                labels_set.add( label );
                labels.add( label );
            }
            if ( sort ) {
                Collections.sort( labels );
            }
            _m = new IntMatrix( labels );
        }
        updateCounts( _m, counter, g );
    }

    private final static void checkPreconditions( final IteratingPhylogenyParser p,
                                                  final Phylogeny species_tree,
                                                  final REROOTING rerooting,
                                                  final String outgroup,
                                                  final int first,
                                                  final int last ) throws RIOException, IOException {
        final Phylogeny g0 = p.next();
        if ( ( g0 == null ) || g0.isEmpty() ) {
            throw new RIOException( "input file does not seem to contain any gene trees" );
        }
        if ( g0.getNumberOfExternalNodes() < 2 ) {
            throw new RIOException( "input file does not seem to contain any useable gene trees" );
        }
        if ( !species_tree.isRooted() ) {
            throw new RIOException( "species tree is not rooted" );
        }
        if ( !( ( last == DEFAULT_RANGE ) && ( first == DEFAULT_RANGE ) )
                && ( ( last < first ) || ( last < 0 ) || ( first < 0 ) ) ) {
            throw new RIOException( "attempt to set range (0-based) of gene to analyze to: from " + first + " to "
                    + last );
        }
        if ( ( rerooting == REROOTING.OUTGROUP ) && ForesterUtil.isEmpty( outgroup ) ) {
            throw new RIOException( "outgroup not set for midpoint rooting" );
        }
        if ( ( rerooting != REROOTING.OUTGROUP ) && !ForesterUtil.isEmpty( outgroup ) ) {
            throw new RIOException( "outgroup only used for midpoint rooting" );
        }
        if ( ( rerooting == REROOTING.MIDPOINT ) && ( PhylogenyMethods.calculateMaxDistanceToRoot( g0 ) <= 0 ) ) {
            throw new RIOException( "attempt to use midpoint rooting on gene trees which seem to have no (positive) branch lengths (cladograms)" );
        }
        if ( rerooting == REROOTING.OUTGROUP ) {
            try {
                g0.getNode( outgroup );
            }
            catch ( final IllegalArgumentException e ) {
                throw new RIOException( "cannot perform re-rooting by outgroup: " + e.getLocalizedMessage() );
            }
        }
    }

    private final static void checkPreconditions( final Phylogeny[] gene_trees,
                                                  final Phylogeny species_tree,
                                                  final REROOTING rerooting,
                                                  final String outgroup,
                                                  final int first,
                                                  final int last ) throws RIOException {
        if ( !species_tree.isRooted() ) {
            throw new RIOException( "species tree is not rooted" );
        }
        if ( !( ( last == DEFAULT_RANGE ) && ( first == DEFAULT_RANGE ) )
                && ( ( last < first ) || ( last >= gene_trees.length ) || ( last < 0 ) || ( first < 0 ) ) ) {
            throw new RIOException( "attempt to set range (0-based) of gene to analyze to: from " + first + " to "
                    + last + " (out of " + gene_trees.length + ")" );
        }
        if ( ( rerooting == REROOTING.OUTGROUP ) && ForesterUtil.isEmpty( outgroup ) ) {
            throw new RIOException( "outgroup not set for midpoint rooting" );
        }
        if ( ( rerooting != REROOTING.OUTGROUP ) && !ForesterUtil.isEmpty( outgroup ) ) {
            throw new RIOException( "outgroup only used for midpoint rooting" );
        }
        if ( ( rerooting == REROOTING.MIDPOINT )
                && ( PhylogenyMethods.calculateMaxDistanceToRoot( gene_trees[ 0 ] ) <= 0 ) ) {
            throw new RIOException( "attempt to use midpoint rooting on gene trees which seem to have no (positive) branch lengths (cladograms)" );
        }
        if ( rerooting == REROOTING.OUTGROUP ) {
            try {
                gene_trees[ 0 ].getNode( outgroup );
            }
            catch ( final IllegalArgumentException e ) {
                throw new RIOException( "cannot perform re-rooting by outgroup: " + e.getLocalizedMessage() );
            }
        }
    }

    private final static String obtainLabel( final Set<String> labels_set, final PhylogenyNode n ) throws RIOException {
        String label;
        if ( n.getNodeData().isHasSequence() && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getName() ) ) {
            label = n.getNodeData().getSequence().getName();
        }
        else if ( n.getNodeData().isHasSequence() && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getSymbol() ) ) {
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
        return label;
    }

    private final static Phylogeny[] parseGeneTrees( final File gene_trees_file ) throws FileNotFoundException,
            IOException {
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final PhylogenyParser p = ParserUtils.createParserDependingOnFileType( gene_trees_file, true );
        if ( p instanceof NHXParser ) {
            final NHXParser nhx = ( NHXParser ) p;
            nhx.setReplaceUnderscores( false );
            nhx.setIgnoreQuotes( true );
            nhx.setTaxonomyExtraction( TAXONOMY_EXTRACTION.AGRESSIVE );
        }
        else if ( p instanceof NexusPhylogeniesParser ) {
            final NexusPhylogeniesParser nex = ( NexusPhylogeniesParser ) p;
            nex.setReplaceUnderscores( false );
            nex.setIgnoreQuotes( true );
            nex.setTaxonomyExtraction( TAXONOMY_EXTRACTION.AGRESSIVE );
        }
        return factory.create( gene_trees_file, p );
    }

    private final static void removeSingleDescendentsNodes( final Phylogeny species_tree, final boolean verbose ) {
        final int o = PhylogenyMethods.countNumberOfOneDescendantNodes( species_tree );
        if ( o > 0 ) {
            if ( verbose ) {
                System.out.println( "warning: species tree has " + o
                        + " internal nodes with only one descendent which are therefore going to be removed" );
            }
            PhylogenyMethods.deleteInternalNodesWithOnlyOneDescendent( species_tree );
        }
    }

    private final static void updateCounts( final IntMatrix m, final int counter, final Phylogeny g )
            throws RIOException {
        PhylogenyMethods.preOrderReId( g );
        final HashMap<String, PhylogenyNode> map = PhylogenyMethods.createNameToExtNodeMap( g );
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

    public enum REROOTING {
        NONE, BY_ALGORITHM, MIDPOINT, OUTGROUP;
    }
}
