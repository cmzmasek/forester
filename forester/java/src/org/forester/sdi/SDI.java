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

import java.util.HashMap;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;
import org.forester.util.ForesterUtil;

/*
 * Implements our algorithm for speciation - duplication inference (SDI). <p>
 * Reference: </p> <ul> <li>Zmasek, C.M. and Eddy, S.R. (2001) "A simple
 * algorithm to infer gene duplication and speciation events on a gene tree".
 * Bioinformatics, in press. </ul> <p> The initialization is accomplished by:
 * </p> <ul> <li>method "linkExtNodesOfG()" of class SDI: setting the links for
 * the external nodes of the gene tree <li>"preorderReID(int)" from class
 * Phylogeny: numbering of nodes of the species tree in preorder <li>the
 * optional stripping of the species tree is accomplished by method
 * "stripTree(Phylogeny,Phylogeny)" of class Phylogeny </ul> <p> The recursion
 * part is accomplished by this class' method
 * "geneTreePostOrderTraversal(PhylogenyNode)". <p> Requires JDK 1.2 or greater.
 * 
 * @see SDI#linkNodesOfG()
 * 
 * @see Phylogeny#preorderReID(int)
 * 
 * @see
 * PhylogenyMethods#taxonomyBasedDeletionOfExternalNodes(Phylogeny,Phylogeny)
 * 
 * @see #geneTreePostOrderTraversal(PhylogenyNode)
 * 
 * @author Christian M. Zmasek
 * 
 * @version 1.102 -- last modified: 10/02/01
 */
public class SDI {

    final Phylogeny _gene_tree;
    final Phylogeny _species_tree;
    int             _duplications_sum; // Sum of duplications.
    int             _mapping_cost;    // Mapping cost "L".

    /**
     * Constructor which sets the gene tree and the species tree to be compared.
     * species_tree is the species tree to which the gene tree gene_tree will be
     * compared to - with method "infer(boolean)". Both Trees must be completely
     * binary and rooted. The actual inference is accomplished with method
     * "infer(boolean)". The mapping cost L can then be calculated with method
     * "computeMappingCost()".
     * <p>
     * (Last modified: 01/11/01)
     * 
     * @see #infer(boolean)
     * @see SDI#computeMappingCostL()
     * @param gene_tree
     *            reference to a rooted binary gene Phylogeny to which assign
     *            duplication vs speciation, must have species names in the
     *            species name fields for all external nodes
     * @param species_tree
     *            reference to a rooted binary species Phylogeny which might get
     *            stripped in the process, must have species names in the
     *            species name fields for all external nodes
     * @throws SDIException 
     */
    public SDI( final Phylogeny gene_tree, final Phylogeny species_tree ) throws SDIException {
        if ( species_tree.isEmpty() || gene_tree.isEmpty() ) {
            throw new IllegalArgumentException( "attempt to infer duplications using empty tree(s)" );
        }
        if ( !species_tree.isRooted() ) {
            throw new IllegalArgumentException( "attempt to infer duplications on unrooted species tree" );
        }
        _gene_tree = gene_tree;
        _species_tree = species_tree;
        _mapping_cost = -1;
        _duplications_sum = 0;
        PhylogenyMethods.preOrderReId( getSpeciesTree() );
        linkNodesOfG();
        geneTreePostOrderTraversal( getGeneTree().getRoot() );
    }

    /**
     * Computes the cost of mapping the gene tree gene_tree onto the species
     * tree species_tree. Before this method can be called, the mapping has to
     * be calculated with method "infer(boolean)".
     * <p>
     * Reference. Zhang, L. (1997) On a Mirkin-Muchnik-Smith Conjecture for
     * Comparing Molecular Phylogenies. Journal of Computational Biology 4
     * 177-187.
     * 
     * @return the mapping cost "L"
     */
    public int computeMappingCostL() {
        _species_tree.levelOrderReID();
        _mapping_cost = 0;
        computeMappingCostHelper( _gene_tree.getRoot() );
        return _mapping_cost;
    }

    /**
     * Returns the number of duplications.
     * 
     * @return number of duplications
     */
    public int getDuplicationsSum() {
        return _duplications_sum;
    }

    /**
     * Returns the gene tree.
     * 
     * @return gene tree
     */
    public Phylogeny getGeneTree() {
        return _gene_tree;
    }

    /**
     * Returns the species tree.
     * 
     * @return species tree
     */
    public Phylogeny getSpeciesTree() {
        return _species_tree;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( getClass() );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        sb.append( "Duplications sum                   : " + getDuplicationsSum() );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        sb.append( "mapping cost L                     : " + computeMappingCostL() );
        return sb.toString();
    }

    /**
     * Traverses the subtree of PhylogenyNode g in postorder, calculating the
     * mapping function M, and determines which nodes represent speciation
     * events and which ones duplication events.
     * <p>
     * Preconditions: Mapping M for external nodes must have been calculated and
     * the species tree must be labelled in preorder.
     * <p>
     * (Last modified: 01/11/01)
     * 
     * @param g
     *            starting node of a gene tree - normally the root
     */
    void geneTreePostOrderTraversal( final PhylogenyNode g ) {
        PhylogenyNode a, b;
        if ( !g.isExternal() ) {
            geneTreePostOrderTraversal( g.getChildNode( 0 ) );
            geneTreePostOrderTraversal( g.getChildNode( 1 ) );
            a = g.getChildNode( 0 ).getLink();
            b = g.getChildNode( 1 ).getLink();
            while ( a != b ) {
                if ( a.getId() > b.getId() ) {
                    a = a.getParent();
                }
                else {
                    b = b.getParent();
                }
            }
            g.setLink( a );
            // Determines whether dup. or spec.
            Event event = null;
            if ( ( a == g.getChildNode( 0 ).getLink() ) || ( a == g.getChildNode( 1 ).getLink() ) ) {
                event = Event.createSingleDuplicationEvent();
                ++_duplications_sum;
            }
            else {
                event = Event.createSingleSpeciationEvent();
            }
            g.getNodeData().setEvent( event );
        }
    } // geneTreePostOrderTraversal( PhylogenyNode )

    /**
     * Calculates the mapping function for the external nodes of the gene tree:
     * links (sets the field "link" of PhylogenyNode) each external
     * PhylogenyNode of gene_tree to the external PhylogenyNode of species_tree
     * which has the same species name.
     * @throws SDIException 
     */
    final void linkNodesOfG() throws SDIException {
        final Map<String, PhylogenyNode> speciestree_ext_nodes = new HashMap<String, PhylogenyNode>();
        final TaxonomyComparisonBase tax_comp_base = determineTaxonomyComparisonBase();
        // Put references to all external nodes of the species tree into a map.
        // Stringyfied taxonomy is the key, node is the value.
        for( final PhylogenyNodeIterator iter = _species_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode s = iter.next();
            final String tax_str = SDIutil.taxonomyToString( s, tax_comp_base );
            if ( speciestree_ext_nodes.containsKey( tax_str ) ) {
                throw new IllegalArgumentException( "taxonomy [" + s.getNodeData().getTaxonomy()
                        + "] is not unique in species phylogeny" );
            }
            speciestree_ext_nodes.put( tax_str, s );
        }
        // Retrieve the reference to the node with a matching stringyfied taxonomy.
        for( final PhylogenyNodeIterator iter = _gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode g = iter.next();
            final String tax_str = SDIutil.taxonomyToString( g, tax_comp_base );
            final PhylogenyNode s = speciestree_ext_nodes.get( tax_str );
            if ( s == null ) {
                throw new IllegalArgumentException( "taxonomy [" + g.getNodeData().getTaxonomy()
                        + "] not present in species tree" );
            }
            g.setLink( s );
        }
    }

    /**
     * Updates the mapping function M after the root of the gene tree has been
     * moved by one branch. It calculates M for the root of the gene tree and
     * one of its two children.
     * <p>
     * To be used ONLY by method "SDIunrooted.fastInfer(Phylogeny,Phylogeny)".
     * <p>
     * (Last modfied: 10/02/01)
     * 
     * @param prev_root_was_dup
     *            true if the previous root was a duplication, false otherwise
     * @param prev_root_c1
     *            child 1 of the previous root
     * @param prev_root_c2
     *            child 2 of the previous root
     * @return number of duplications which have been assigned in gene tree
     */
    int updateM( final boolean prev_root_was_dup, final PhylogenyNode prev_root_c1, final PhylogenyNode prev_root_c2 ) {
        final PhylogenyNode root = getGeneTree().getRoot();
        if ( ( root.getChildNode1() == prev_root_c1 ) || ( root.getChildNode2() == prev_root_c1 ) ) {
            calculateMforNode( prev_root_c1 );
        }
        else {
            calculateMforNode( prev_root_c2 );
        }
        Event event = null;
        if ( prev_root_was_dup ) {
            event = Event.createSingleDuplicationEvent();
        }
        else {
            event = Event.createSingleSpeciationEvent();
        }
        root.getNodeData().setEvent( event );
        calculateMforNode( root );
        return getDuplicationsSum();
    } // updateM( boolean, PhylogenyNode, PhylogenyNode )

    // Helper method for updateM( boolean, PhylogenyNode, PhylogenyNode )
    // Calculates M for PhylogenyNode n, given that M for the two children
    // of n has been calculated.
    // (Last modified: 10/02/01)
    private void calculateMforNode( final PhylogenyNode n ) {
        if ( !n.isExternal() ) {
            final boolean was_duplication = n.isDuplication();
            PhylogenyNode a = n.getChildNode1().getLink();
            PhylogenyNode b = n.getChildNode2().getLink();
            while ( a != b ) {
                if ( a.getId() > b.getId() ) {
                    a = a.getParent();
                }
                else {
                    b = b.getParent();
                }
            }
            n.setLink( a );
            Event event = null;
            if ( ( a == n.getChildNode1().getLink() ) || ( a == n.getChildNode2().getLink() ) ) {
                event = Event.createSingleDuplicationEvent();
                if ( !was_duplication ) {
                    ++_duplications_sum;
                }
            }
            else {
                event = Event.createSingleSpeciationEvent();
                if ( was_duplication ) {
                    --_duplications_sum;
                }
            }
            n.getNodeData().setEvent( event );
        }
    } // calculateMforNode( PhylogenyNode )

    // Helper method for "computeMappingCost()".
    private void computeMappingCostHelper( final PhylogenyNode g ) {
        if ( !g.isExternal() ) {
            computeMappingCostHelper( g.getChildNode1() );
            computeMappingCostHelper( g.getChildNode2() );
            if ( ( g.getLink() != g.getChildNode1().getLink() ) && ( g.getLink() != g.getChildNode2().getLink() ) ) {
                _mapping_cost += ( ( g.getChildNode1().getLink().getId() + g.getChildNode2().getLink().getId() )
                        - ( 2 * g.getLink().getId() ) - 2 );
            }
            else if ( ( g.getLink() != g.getChildNode1().getLink() ) && ( g.getLink() == g.getChildNode2().getLink() ) ) {
                _mapping_cost += ( ( g.getChildNode1().getLink().getId() - g.getLink().getId() ) + 1 );
            }
            else if ( ( g.getLink() == g.getChildNode1().getLink() ) && ( g.getLink() != g.getChildNode2().getLink() ) ) {
                _mapping_cost += ( ( g.getChildNode2().getLink().getId() - g.getLink().getId() ) + 1 );
            }
            else {
                _mapping_cost++;
            }
        }
    }

    private TaxonomyComparisonBase determineTaxonomyComparisonBase() {
        TaxonomyComparisonBase base = null;
        boolean all_have_id = true;
        boolean all_have_code = true;
        boolean all_have_sn = true;
        for( final PhylogenyNodeIterator iter = _species_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                final Taxonomy tax = n.getNodeData().getTaxonomy();
                if ( ( tax.getIdentifier() == null ) || ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) {
                    all_have_id = false;
                }
                if ( ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                    all_have_code = false;
                }
                if ( ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                    all_have_sn = false;
                }
            }
            else {
                throw new IllegalArgumentException( "species tree node [" + n + "] has no taxonomic data" );
            }
        }
        for( final PhylogenyNodeIterator iter = _gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                final Taxonomy tax = n.getNodeData().getTaxonomy();
                if ( ( tax.getIdentifier() == null ) || ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) {
                    all_have_id = false;
                }
                if ( ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                    all_have_code = false;
                }
                if ( ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                    all_have_sn = false;
                }
            }
            else {
                throw new IllegalArgumentException( "gene tree node [" + n + "] has no taxonomic data" );
            }
        }
        if ( all_have_id ) {
            base = TaxonomyComparisonBase.ID;
        }
        else if ( all_have_code ) {
            base = TaxonomyComparisonBase.CODE;
        }
        else if ( all_have_sn ) {
            base = TaxonomyComparisonBase.SCIENTIFIC_NAME;
        }
        else {
            throw new IllegalArgumentException( "gene tree and species tree have incomparable taxonomies" );
        }
        return base;
    }

    /**
     * Calculates the mapping function for the external nodes of the gene tree:
     * links (sets the field "link" of PhylogenyNode) each external by taxonomy
     * identifier
     * PhylogenyNode of gene_tree to the external PhylogenyNode of species_tree
     * which has the same species name.
     * Olivier CHABROL : olivier.chabrol@univ-provence.fr
     */
    private final void linkNodesOfGByTaxonomyIdentifier() {
        final HashMap<String, PhylogenyNode> speciestree_ext_nodes = new HashMap<String, PhylogenyNode>();
        if ( _species_tree.getFirstExternalNode().isRoot() ) {
            speciestree_ext_nodes.put( _species_tree.getFirstExternalNode().getNodeData().getTaxonomy().getIdentifier()
                    .getValue(), _species_tree.getFirstExternalNode() );
        }
        else {
            for( final PhylogenyNodeIterator iter = _species_tree.iteratorExternalForward(); iter.hasNext(); ) {
                final PhylogenyNode s = iter.next();
                speciestree_ext_nodes.put( s.getNodeData().getTaxonomy().getIdentifier().getValue(), s );
            }
        }
        for( final PhylogenyNodeIterator iter = _gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode g = iter.next();
            final PhylogenyNode s = speciestree_ext_nodes
                    .get( g.getNodeData().getTaxonomy().getIdentifier().getValue() );
            if ( s == null ) {
                String message = "species [" + g.getNodeData().getTaxonomy().getIdentifier().getValue();
                message += "] not present in species tree";
                throw new IllegalArgumentException( message );
            }
            g.setLink( s );
        }
    }
} // End of class SDIse.
