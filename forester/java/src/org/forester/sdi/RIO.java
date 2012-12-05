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

package org.forester.sdi;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public final class RIO {

    private final static boolean                      ROOT_BY_MINIMIZING_SUM_OF_DUPS = true;
    private final static boolean                      ROOT_BY_MINIMIZING_TREE_HEIGHT = true;
    private Phylogeny[]                               _analyzed_gene_trees;
    private HashMap<String, HashMap<String, Integer>> _o_maps;
    private HashMap<String, HashMap<String, Integer>> _so_maps;
    private HashMap<String, HashMap<String, Integer>> _up_maps;
    private List<String>                              _seq_names;
    private int                                       _samples;
    private int                                       _ext_nodes;

    /**
     * Default constructor.
     * @throws SDIException 
     * @throws IOException 
     * @throws RioException 
     */
    public RIO( final File gene_trees_file, final Phylogeny species_tree, final String query ) throws IOException,
            SDIException, RioException {
        if ( ForesterUtil.isEmpty( query ) ) {
            throw new IllegalArgumentException( "query is empty" );
        }
        init();
        inferOrthologs( gene_trees_file, species_tree, query );
    }

    public RIO( final File gene_trees_file, final Phylogeny species_tree ) throws IOException, SDIException,
            RioException {
        init();
        inferOrthologs( gene_trees_file, species_tree, null );
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

    /**
     * Returns a HashMap containing the inferred "ultra paralogs" of the
     * external gene tree node with the sequence name seq_name. Sequence names
     * are the keys (String), numbers of observations are the values (Int).
     * "ultra paralogs" are to be inferred by method "inferOrthologs". Throws an
     * exception if seq_name is not found. 
     * 
     * @param seq_name
     *            sequence name of a external node of the gene trees
     * @return HashMap containing the inferred ultra paralogs
     *         (name(String)->value(Int))
     */
    public final HashMap<String, Integer> getInferredUltraParalogs( final String seq_name ) {
        if ( _up_maps == null ) {
            return null;
        }
        return _up_maps.get( seq_name );
    }

    public final int getNumberOfSamples() {
        return _samples;
    }

    /**
     * Returns a String containg the names of orthologs of the PhylogenyNode
     * with seq name query_name. The String also contains how many times a
     * particular ortholog has been observed.
     * <p>
     * <ul>
     * The output order is (per line): Name, Ortholog, Subtree neighbor, Super
     * ortholog, Distance
     * </ul>
     * <p>
     * The sort priority of this is determined by sort in the following manner:
     * <ul>
     * <li>0 : Ortholog
     * <li>1 : Ortholog, Super ortholog
     * <li>2 : Super ortholog, Ortholog
     * </ul>
     * <p>
     * Returns "-" if no putative orthologs have been found (given
     * threshold_orthologs).
     * <p>
     * Orthologs are to be inferred by method "inferOrthologs".
     * <p>
     * (Last modified: 05/08/01)
     * 
     * @param query_name
     *            sequence name of a external node of the gene trees
     * @param sort
     *            order and sort priority
     * @param threshold_orthologs
     *            the minimal number of observations for a a sequence to be
     *            reported as orthologous, in percents (0.0-100.0%)
     * @param threshold_subtreeneighborings
     *            the minimal number of observations for a a sequence to be
     *            reported as orthologous, in percents (0.0-100.0%)
     * @return String containing the inferred orthologs, String containing "-"
     *         if no orthologs have been found null in case of error
     */
    public final StringBuffer inferredOrthologsToString( final String query_name, int sort, double threshold_orthologs ) {
        HashMap<String, Integer> o_hashmap = null;
        HashMap<String, Integer> s_hashmap = null;
        String name = "";
        double o = 0.0; // Orthologs.
        double s = 0.0; // Super orthologs.
        double value1 = 0.0;
        double value2 = 0.0;
        final ArrayList<ResultLine> nv = new ArrayList<ResultLine>();
        if ( ( _o_maps == null ) || ( _so_maps == null ) ) {
            throw new RuntimeException( "orthologs have not been calculated (successfully)" );
        }
        if ( ( sort < 0 ) || ( sort > 2 ) ) {
            sort = 1;
        }
        if ( threshold_orthologs < 0.0 ) {
            threshold_orthologs = 0.0;
        }
        else if ( threshold_orthologs > 100.0 ) {
            threshold_orthologs = 100.0;
        }
        o_hashmap = getInferredOrthologs( query_name );
        s_hashmap = getInferredSuperOrthologs( query_name );
        if ( ( o_hashmap == null ) || ( s_hashmap == null ) ) {
            throw new RuntimeException( "Orthologs for " + query_name + " were not established" );
        }
        final StringBuffer orthologs = new StringBuffer();
        if ( _seq_names.size() > 0 ) {
            I: for( int i = 0; i < _seq_names.size(); ++i ) {
                name = _seq_names.get( i );
                if ( name.equals( query_name ) ) {
                    continue I;
                }
                o = getBootstrapValueFromHash( o_hashmap, name );
                if ( o < threshold_orthologs ) {
                    continue I;
                }
                s = getBootstrapValueFromHash( s_hashmap, name );
                switch ( sort ) {
                    case 0:
                        nv.add( new ResultLine( name, o, 5 ) );
                        break;
                    case 1:
                        nv.add( new ResultLine( name, o, s, 5 ) );
                        break;
                    case 2:
                        nv.add( new ResultLine( name, s, o, 5 ) );
                        break;
                    default:
                        nv.add( new ResultLine( name, o, 5 ) );
                }
            } // End of I for loop.
            if ( ( nv != null ) && ( nv.size() > 0 ) ) {
                orthologs.append( "seq name\t\tortho\ts-ortho" + ForesterUtil.LINE_SEPARATOR );
                final ResultLine[] nv_array = new ResultLine[ nv.size() ];
                for( int j = 0; j < nv.size(); ++j ) {
                    nv_array[ j ] = nv.get( j );
                }
                Arrays.sort( nv_array );
                for( final ResultLine element : nv_array ) {
                    name = element.getKey();
                    value1 = element.getValue1();
                    value2 = element.getValue2();
                    orthologs.append( addNameAndValues( name, value1, value2, sort ) );
                }
            }
        }
        // No orthologs found.
        if ( ( orthologs == null ) || ( orthologs.length() < 1 ) ) {
            orthologs.append( "-" );
        }
        return orthologs;
    }

    /**
     * Returns a String containg the names of orthologs of the PhylogenyNode
     * with seq name query_name. The String also contains how many times a
     * particular ortholog has been observed. Returns "-" if no putative
     * orthologs have been found (given threshold_orthologs).
     * <p>
     * Orthologs are to be inferred by method "inferOrthologs".
     * 
     * @param query_name
     *            sequence name of a external node of the gene trees
     * @param return_dists
     * @param threshold_ultra_paralogs
     *            between 1 and 100
     * @return String containing the inferred orthologs, String containing "-"
     *         if no orthologs have been found null in case of error
     */
    public final String inferredUltraParalogsToString( final String query_name, double threshold_ultra_paralogs ) {
        HashMap<String, Integer> sp_hashmap = null;
        String name = "", ultra_paralogs = "";
        int sort = 0;
        double sp = 0.0;
        double value1 = 0.0;
        double value2 = 0.0;
        final List<ResultLine> nv = new ArrayList<ResultLine>();
        if ( threshold_ultra_paralogs < 1.0 ) {
            threshold_ultra_paralogs = 1.0;
        }
        else if ( threshold_ultra_paralogs > 100.0 ) {
            threshold_ultra_paralogs = 100.0;
        }
        if ( _up_maps == null ) {
            throw new RuntimeException( "Ultra paralogs have not been calculated (successfully)." );
        }
        sp_hashmap = getInferredUltraParalogs( query_name );
        if ( sp_hashmap == null ) {
            throw new RuntimeException( "Ultra paralogs for " + query_name + " were not established" );
        }
        if ( _seq_names.size() > 0 ) {
            I: for( int i = 0; i < _seq_names.size(); ++i ) {
                name = _seq_names.get( i );
                if ( name.equals( query_name ) ) {
                    continue I;
                }
                sp = getBootstrapValueFromHash( sp_hashmap, name );
                if ( sp < threshold_ultra_paralogs ) {
                    continue I;
                }
                nv.add( new ResultLine( name, sp, 5 ) );
            } // End of I for loop.
            if ( ( nv != null ) && ( nv.size() > 0 ) ) {
                final ResultLine[] nv_array = new ResultLine[ nv.size() ];
                for( int j = 0; j < nv.size(); ++j ) {
                    nv_array[ j ] = nv.get( j );
                }
                Arrays.sort( nv_array );
                sort = 90;
                for( final ResultLine element : nv_array ) {
                    name = element.getKey();
                    value1 = element.getValue1();
                    value2 = element.getValue2();
                    ultra_paralogs += addNameAndValues( name, value1, value2, sort );
                }
            }
        }
        // No ultra paralogs found.
        if ( ( ultra_paralogs == null ) || ( ultra_paralogs.length() < 1 ) ) {
            ultra_paralogs = "-";
        }
        return ultra_paralogs;
    }

    // Helper method for inferredOrthologsToString.
    // inferredOrthologsToArrayList,
    // and inferredUltraParalogsToString.
    private final double getBootstrapValueFromHash( final HashMap<String, Integer> h, final String name ) {
        if ( !h.containsKey( name ) ) {
            return 0.0;
        }
        final int i = h.get( name );
        return ( ( i * 100.0 ) / getNumberOfSamples() );
    }

    /**
     * Returns a HashMap containing the inferred orthologs of the external gene
     * tree node with the sequence name seq_name. Sequence names are the keys
     * (String), numbers of observations are the values (Int). Orthologs are to
     * be inferred by method "inferOrthologs". Throws an exception if seq_name
     * is not found.
     * 
     * @param seq_name
     *            sequence name of a external node of the gene trees
     * @return HashMap containing the inferred orthologs
     *         (name(String)->value(Int))
     */
    private final HashMap<String, Integer> getInferredOrthologs( final String seq_name ) {
        if ( _o_maps == null ) {
            return null;
        }
        return _o_maps.get( seq_name );
    }

    /**
     * Returns a HashMap containing the inferred "super orthologs" of the
     * external gene tree node with the sequence name seq_name. Sequence names
     * are the keys (String), numbers of observations are the values (Int).
     * Super orthologs are to be inferred by method "inferOrthologs". Throws an
     * exception if seq_name is not found.
     * 
     * @param seq_name
     *            sequence name of a external node of the gene trees
     * @return HashMap containing the inferred super orthologs
     *         (name(String)->value(Int))
     */
    private final HashMap<String, Integer> getInferredSuperOrthologs( final String seq_name ) {
        if ( _so_maps == null ) {
            return null;
        }
        return _so_maps.get( seq_name );
    }

    /**
     * Infers the orthologs (as well the "super orthologs", the "subtree
     * neighbors", and the "ultra paralogs") for each external node of the gene
     * Trees in multiple tree File gene_trees_file (=output of PHYLIP NEIGHBOR,
     * for example). Tallies how many times each sequence is (super-)
     * orthologous towards the query. Tallies how many times each sequence is
     * ultra paralogous towards the query. Tallies how many times each sequence
     * is a subtree neighbor of the query. Gene duplications are inferred using
     * SDI. Modifies its argument species_tree. Is a little faster than
     * "inferOrthologs(File,Phylogeny)" since orthologs are only inferred for
     * query.
     * <p>
     * To obtain the results use the methods listed below.
     * 
     * @param gene_trees_file
     *            a File containing gene Trees in NH format, which is the result
     *            of performing a bootstrap analysis in PHYLIP
     * @param species_tree
     *            a species Phylogeny, which has species names in its species
     *            fields
     * @param query
     *            the sequence name of the squence whose orthologs are to be
     *            inferred
     * @throws SDIException 
     * @throws RioException 
     * @throws IOException 
     * @throws FileNotFoundException 
     */
    private final void inferOrthologs( final File gene_trees_file, final Phylogeny species_tree, final String query )
            throws SDIException, RioException, FileNotFoundException, IOException {
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
        PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( gene_trees[ 0 ], species_tree );
        if ( species_tree.isEmpty() ) {
            throw new RioException( "failed to establish species based mapping between gene and species trees" );
        }
        if ( !ForesterUtil.isEmpty( query ) ) {
            PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gene_trees[ 0 ] );
            if ( gene_trees[ 0 ].isEmpty() ) {
                throw new RioException( "failed to establish species based mapping between gene and species trees" );
            }
            _seq_names = getAllExternalSequenceNames( gene_trees[ 0 ] );
            if ( ( _seq_names == null ) || ( _seq_names.size() < 1 ) ) {
                throw new RioException( "could not get sequence names" );
            }
            _o_maps = new HashMap<String, HashMap<String, Integer>>();
            _so_maps = new HashMap<String, HashMap<String, Integer>>();
            _up_maps = new HashMap<String, HashMap<String, Integer>>();
            _o_maps.put( query, new HashMap<String, Integer>( _seq_names.size() ) );
            _so_maps.put( query, new HashMap<String, Integer>( _seq_names.size() ) );
            _up_maps.put( query, new HashMap<String, Integer>( _seq_names.size() ) );
        }
        _analyzed_gene_trees = new Phylogeny[ gene_trees.length ];
        int c = 0;
        for( final Phylogeny gt : gene_trees ) {
            // Removes from gene_tree all species not found in species_tree.
            PhylogenyMethods.taxonomyBasedDeletionOfExternalNodes( species_tree, gt );
            if ( gt.isEmpty() ) {
                throw new RioException( "failed to establish species based mapping between gene and species trees" );
            }
            _analyzed_gene_trees[ c++ ] = inferOrthologsHelper( gt, species_tree, query );
        }
        setNumberOfSamples( gene_trees.length );
    }

    // Helper method which performs the actual ortholog inference for
    // the external node with seqname query.
    private final Phylogeny inferOrthologsHelper( final Phylogeny gene_tree,
                                                  final Phylogeny species_tree,
                                                  final String query ) throws SDIException, RioException {
        final SDIR sdiunrooted = new SDIR();
        final Phylogeny assigned_tree = sdiunrooted.infer( gene_tree,
                                                           species_tree,
                                                           false,
                                                           RIO.ROOT_BY_MINIMIZING_SUM_OF_DUPS,
                                                           RIO.ROOT_BY_MINIMIZING_TREE_HEIGHT,
                                                           true,
                                                           1 )[ 0 ];
        setExtNodesOfAnalyzedGeneTrees( assigned_tree.getNumberOfExternalNodes() );
        if ( !ForesterUtil.isEmpty( query ) ) {
            final List<PhylogenyNode> nodes = getNodesViaSequenceName( assigned_tree, query );
            if ( nodes.size() > 1 ) {
                throw new RioException( "node named [" + query + "] not unique" );
            }
            else if ( nodes.isEmpty() ) {
                throw new RioException( "no node containing a sequence named [" + query + "] found" );
            }
            final PhylogenyNode query_node = nodes.get( 0 );
            updateCounts( _o_maps, query, PhylogenyMethods.getOrthologousNodes( assigned_tree, query_node ) );
            updateCounts( _so_maps, query, PhylogenyMethods.getSuperOrthologousNodes( query_node ) );
            updateCounts( _up_maps, query, PhylogenyMethods.getUltraParalogousNodes( query_node ) );
        }
        return assigned_tree;
    }

    private final void init() {
        _o_maps = null;
        _so_maps = null;
        _up_maps = null;
        _seq_names = null;
        _samples = 1;
        _ext_nodes = 0;
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

    // Helper for doInferOrthologs( Phylogeny, Phylogeny, String )
    // and doInferOrthologs( Phylogeny, Phylogeny ).
    private final void updateCounts( final HashMap<String, HashMap<String, Integer>> counter_map,
                                     final String query_seq_name,
                                     final List<PhylogenyNode> nodes ) {
        final HashMap<String, Integer> hash_map = counter_map.get( query_seq_name );
        if ( hash_map == null ) {
            throw new RuntimeException( "unexpected error in updateCounts" );
        }
        for( int j = 0; j < nodes.size(); ++j ) {
            String seq_name;
            if ( ( nodes.get( j ) ).getNodeData().isHasSequence()
                    && !ForesterUtil.isEmpty( ( nodes.get( j ) ).getNodeData().getSequence().getName() ) ) {
                seq_name = ( nodes.get( j ) ).getNodeData().getSequence().getName();
            }
            else {
                seq_name = ( nodes.get( j ) ).getName();
            }
            if ( hash_map.containsKey( seq_name ) ) {
                hash_map.put( seq_name, hash_map.get( seq_name ) + 1 );
            }
            else {
                hash_map.put( seq_name, 1 );
            }
        }
    }

    public final static IntMatrix calculateOrthologTable( final Phylogeny[] analyzed_gene_trees ) throws RioException {
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
                    throw new RioException( "node \"" + mx + "\" not present in gene tree #" + counter );
                }
                for( int y = 0; y < m.size(); ++y ) {
                    final String my = m.getLabel( y );
                    final PhylogenyNode ny = map.get( my );
                    if ( ny == null ) {
                        throw new RioException( "node \"" + my + "\" not present in gene tree #" + counter );
                    }
                    if ( !PhylogenyMethods.calculateLCAonTreeWithIdsInPreOrder( nx, ny ).isDuplication() ) {
                        m.inreaseByOne( x, y );
                    }
                }
            }
        }
        return m;
    }

    /**
     * Returns the order in which ortholog (o), "super ortholog" (s) and
     * distance (d) are returned and sorted (priority of sort always goes from
     * left to right), given sort. For the meaning of sort
     * 
     * @see #inferredOrthologsToString(String,int,double,double)
     *      
     * @param sort
     *            determines order and sort priority
     * @return String indicating the order
     */
    public final static String getOrder( final int sort ) {
        String order = "";
        switch ( sort ) {
            case 0:
                order = "orthologies";
                break;
            case 1:
                order = "orthologies > super orthologies";
                break;
            case 2:
                order = "super orthologies > orthologies";
                break;
            default:
                order = "orthologies";
                break;
        }
        return order;
    }

    public final static StringBuffer getOrderHelp() {
        final StringBuffer sb = new StringBuffer();
        sb.append( "  0: orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  1: orthologies > super orthologies" + ForesterUtil.LINE_SEPARATOR );
        sb.append( "  2: super orthologies > orthologies" + ForesterUtil.LINE_SEPARATOR );
        return sb;
    }

    // Helper method for inferredOrthologsToString
    // and inferredUltraParalogsToString.
    private final static String addNameAndValues( final String name,
                                                  final double value1,
                                                  final double value2,
                                                  final int sort ) {
        final java.text.DecimalFormat df = new java.text.DecimalFormat( "0.#####" );
        df.setDecimalSeparatorAlwaysShown( false );
        String line = "";
        if ( name.length() < 8 ) {
            line += ( name + "\t\t\t" );
        }
        else if ( name.length() < 16 ) {
            line += ( name + "\t\t" );
        }
        else {
            line += ( name + "\t" );
        }
        switch ( sort ) {
            case 0:
                line += addToLine( value1, df );
                line += "-\t";
                break;
            case 1:
                line += addToLine( value1, df );
                line += addToLine( value2, df );
                break;
            case 2:
                line += addToLine( value2, df );
                line += addToLine( value1, df );
                break;
            case 90:
                line += addToLine( value1, df );
                line += "-\t";
                break;
            case 91:
                line += addToLine( value1, df );
                line += addToLine( value2, df );
                break;
        }
        line += ForesterUtil.LINE_SEPARATOR;
        return line;
    }

    // Helper for addNameAndValues.
    private final static String addToLine( final double value, final java.text.DecimalFormat df ) {
        String s = "";
        if ( value != ResultLine.DEFAULT ) {
            s = df.format( value ) + "\t";
        }
        else {
            s = "-\t";
        }
        return s;
    }

    private final static List<String> getAllExternalSequenceNames( final Phylogeny phy ) throws RioException {
        final List<String> names = new ArrayList<String>();
        for( final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeData().isHasSequence() && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getName() ) ) {
                names.add( n.getNodeData().getSequence().getName() );
            }
            else if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                names.add( n.getName() );
            }
            else {
                throw new RioException( "node has no (sequence) name: " + n );
            }
        }
        return names;
    }

    private final static List<PhylogenyNode> getNodesViaSequenceName( final Phylogeny phy, final String seq_name ) {
        final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeData().isHasSequence() && n.getNodeData().getSequence().getName().equals( seq_name ) ) {
                nodes.add( n );
            }
            if ( !n.getNodeData().isHasSequence() && n.getName().equals( seq_name ) ) {
                nodes.add( n );
            }
        }
        return nodes;
    }

    private final class ResultLine implements Comparable<ResultLine> {

        public static final int DEFAULT = -999;
        private final String    _key;
        private final double    _value1;
        private final double    _value2;
        private int[]           _p;

        ResultLine( final String name, final double value1, final double value2, final int c ) {
            setSigns();
            _key = name;
            _value1 = value1;
            _value2 = value2;
            if ( ( c >= 0 ) && ( c <= 2 ) ) {
                _p[ c ] = -1;
            }
        }

        ResultLine( final String name, final double value1, final int c ) {
            setSigns();
            _key = name;
            _value1 = value1;
            _value2 = ResultLine.DEFAULT;
            if ( c == 0 ) {
                _p[ 0 ] = -1;
            }
        }

        @Override
        public int compareTo( final ResultLine n ) {
            if ( ( getValue1() != ResultLine.DEFAULT ) && ( n.getValue1() != ResultLine.DEFAULT ) ) {
                if ( getValue1() < n.getValue1() ) {
                    return _p[ 0 ];
                }
                if ( getValue1() > n.getValue1() ) {
                    return ( -_p[ 0 ] );
                }
            }
            if ( ( getValue2() != ResultLine.DEFAULT ) && ( n.getValue2() != ResultLine.DEFAULT ) ) {
                if ( getValue2() < n.getValue2() ) {
                    return _p[ 1 ];
                }
                if ( getValue2() > n.getValue2() ) {
                    return ( -_p[ 1 ] );
                }
            }
            return ( getKey().compareTo( n.getKey() ) );
        }

        String getKey() {
            return _key;
        }

        double getValue1() {
            return _value1;
        }

        double getValue2() {
            return _value2;
        }

        private void setSigns() {
            _p = new int[ 2 ];
            _p[ 0 ] = _p[ 1 ] = +1;
        }
    } // ResultLine
}
