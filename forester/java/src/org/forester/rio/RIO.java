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
import org.forester.sdi.GSDIR;
import org.forester.sdi.SDI.ALGORITHM;
import org.forester.sdi.SDI.TaxonomyComparisonBase;
import org.forester.sdi.SDIException;
import org.forester.sdi.SDIR;
import org.forester.util.ForesterUtil;

public final class RIO {

    private final static boolean   ROOT_BY_MINIMIZING_SUM_OF_DUPS = true;
    private final static boolean   ROOT_BY_MINIMIZING_TREE_HEIGHT = true;
    private Phylogeny[]            _analyzed_gene_trees;
    private List<PhylogenyNode>    _removed_gene_tree_nodes;
    private int                    _samples;
    private int                    _ext_nodes;
    private TaxonomyComparisonBase _gsdir_tax_comp_base;

    public RIO( final File gene_trees_file, final Phylogeny species_tree, final ALGORITHM algorithm )
            throws IOException, SDIException, RIOException {
        init();
        inferOrthologs( gene_trees_file, species_tree, algorithm );
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

    public final int getNumberOfSamples() {
        return _samples;
    }

    public final List<PhylogenyNode> getRemovedGeneTreeNodes() {
        return _removed_gene_tree_nodes;
    }

    private final void inferOrthologs( final File gene_trees_file,
                                       final Phylogeny species_tree,
                                       final ALGORITHM algorithm ) throws SDIException, RIOException,
            FileNotFoundException, IOException {
        // Read in first tree to get its sequence names
        // and strip species_tree.
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final PhylogenyParser p = ParserUtils.createParserDependingOnFileType( gene_trees_file, true );
        if ( p instanceof NHXParser ) {
            final NHXParser nhx = ( NHXParser ) p;
            nhx.setReplaceUnderscores( false );
            nhx.setIgnoreQuotes( true );
            nhx.setTaxonomyExtraction( NHXParser.TAXONOMY_EXTRACTION.YES );
        }
        final Phylogeny[] gene_trees = factory.create( gene_trees_file, p );
        // Removes from species_tree all species not found in gene_tree.
        final List<PhylogenyNode> _removed_species_tree_ext_nodes = PhylogenyMethods
                .taxonomyBasedDeletionOfExternalNodes( gene_trees[ 0 ], species_tree );
        if ( species_tree.isEmpty() ) {
            throw new RIOException( "failed to establish species based mapping between gene and species trees" );
        }
        _analyzed_gene_trees = new Phylogeny[ gene_trees.length ];
        int i = 0;
        int gene_tree_ext_nodes = 0;
        for( final Phylogeny gt : gene_trees ) {
            if ( algorithm == ALGORITHM.SDIR ) {
                // Removes from gene_tree all species not found in species_tree.
                PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gt );
                if ( gt.isEmpty() ) {
                    throw new RIOException( "failed to establish species based mapping between gene and species trees" );
                }
                if ( i == 0 ) {
                    gene_tree_ext_nodes = gt.getNumberOfExternalNodes();
                }
                else if ( gene_tree_ext_nodes != gt.getNumberOfExternalNodes() ) {
                    throw new RIOException( "(cleaned up) gene tree #" + ( i + 1 )
                            + " has a different number of external nodes (" + gt.getNumberOfExternalNodes()
                            + ") than those gene trees preceding it (" + gene_tree_ext_nodes + ")" );
                }
            }
            _analyzed_gene_trees[ i ] = performOrthologInference( gt, species_tree, algorithm, i );
            ++i;
        }
        setNumberOfSamples( gene_trees.length );
    }

    private final void init() {
        _samples = 1;
        _ext_nodes = 0;
    }

    private final Phylogeny performOrthologInference( final Phylogeny gene_tree,
                                                      final Phylogeny species_tree,
                                                      final ALGORITHM algorithm,
                                                      final int i ) throws SDIException, RIOException {
        final Phylogeny assigned_tree;
        switch ( algorithm ) {
            case SDIR: {
                final SDIR sdir = new SDIR();
                assigned_tree = sdir.infer( gene_tree,
                                            species_tree,
                                            false,
                                            RIO.ROOT_BY_MINIMIZING_SUM_OF_DUPS,
                                            RIO.ROOT_BY_MINIMIZING_TREE_HEIGHT,
                                            true,
                                            1 )[ 0 ];
                break;
            }
            case GSDIR: {
                final GSDIR gsdir = new GSDIR( gene_tree, species_tree, true, i == 0 );
                assigned_tree = gsdir.getMinDuplicationsSumGeneTrees().get( 0 );
                _gsdir_tax_comp_base = gsdir.getTaxCompBase();
                break;
            }
            default: {
                throw new IllegalArgumentException( "illegal algorithm: " + algorithm );
            }
        }
        setExtNodesOfAnalyzedGeneTrees( assigned_tree.getNumberOfExternalNodes() );
        return assigned_tree;
    }

    private final void setExtNodesOfAnalyzedGeneTrees( final int i ) {
        _ext_nodes = i;
    }

    private final void setNumberOfSamples( int i ) {
        if ( i < 1 ) {
            i = 1;
        }
        _samples = i;
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
                throw new IllegalArgumentException( "node " + n + " has no appropriate label" );
            }
            if ( labels_set.contains( label ) ) {
                throw new IllegalArgumentException( "label " + label + " is not unique" );
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
}
