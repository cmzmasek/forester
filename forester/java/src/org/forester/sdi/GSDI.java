// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/*
 * Implements our algorithm for speciation - duplication inference (SDI). <p>
 * The initialization is accomplished by: </p> <ul> <li>method
 * "linkExtNodesOfG()" of class SDI: setting the links for the external nodes of
 * the gene tree <li>"preorderReID(int)" from class Phylogeny: numbering of
 * nodes of the species tree in preorder <li>the optional stripping of the
 * species tree is accomplished by method "stripTree(Phylogeny,Phylogeny)" of
 * class Phylogeny </ul> <p> The recursion part is accomplished by this class'
 * method "geneTreePostOrderTraversal(PhylogenyNode)". <p> Requires JDK 1.5 or
 * greater.
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
 */
public class GSDI extends SDI {

    private final boolean             _most_parsimonious_duplication_model;
    private final boolean             _strip_gene_tree;
    private final boolean             _strip_species_tree;
    protected int                       _speciation_or_duplication_events_sum;
    protected int                       _speciations_sum;
    private final List<PhylogenyNode> _stripped_gene_tree_nodes;
    private final List<PhylogenyNode> _stripped_species_tree_nodes;
    private final Set<PhylogenyNode>  _mapped_species_tree_nodes;
    private TaxonomyComparisonBase    _tax_comp_base;
    private final SortedSet<String>   _scientific_names_mapped_to_reduced_specificity;

    public GSDI( final Phylogeny gene_tree,
                 final Phylogeny species_tree,
                 final boolean most_parsimonious_duplication_model,
                 final boolean strip_gene_tree,
                 final boolean strip_species_tree ) throws SDIException {
        super( gene_tree, species_tree );
        _speciation_or_duplication_events_sum = 0;
        _speciations_sum = 0;
        _most_parsimonious_duplication_model = most_parsimonious_duplication_model;
        _duplications_sum = 0;
        _strip_gene_tree = strip_gene_tree;
        _strip_species_tree = strip_species_tree;
        _stripped_gene_tree_nodes = new ArrayList<PhylogenyNode>();
        _stripped_species_tree_nodes = new ArrayList<PhylogenyNode>();
        _mapped_species_tree_nodes = new HashSet<PhylogenyNode>();
        _scientific_names_mapped_to_reduced_specificity = new TreeSet<String>();
        linkNodesOfG();
        PhylogenyMethods.preOrderReId( getSpeciesTree() );
        geneTreePostOrderTraversal();
    }

    GSDI( final Phylogeny gene_tree, final Phylogeny species_tree, final boolean most_parsimonious_duplication_model )
            throws SDIException {
        this( gene_tree, species_tree, most_parsimonious_duplication_model, false, false );
    }
    
    public GSDI( final Phylogeny gene_tree,
                 final Phylogeny species_tree,
                 final boolean most_parsimonious_duplication_model,
                 final boolean strip_gene_tree,
                 final boolean strip_species_tree,
                 int x ) throws SDIException {
        super( gene_tree, species_tree );
        _speciation_or_duplication_events_sum = 0;
        _speciations_sum = 0;
        _most_parsimonious_duplication_model = most_parsimonious_duplication_model;
        _duplications_sum = 0;
        _strip_gene_tree = strip_gene_tree;
        _strip_species_tree = strip_species_tree;
        _stripped_gene_tree_nodes = new ArrayList<PhylogenyNode>();
        _stripped_species_tree_nodes = new ArrayList<PhylogenyNode>();
        _mapped_species_tree_nodes = new HashSet<PhylogenyNode>();
        _scientific_names_mapped_to_reduced_specificity = new TreeSet<String>();
       
    }

    // s is the node on the species tree g maps to.
    private final void determineEvent( final PhylogenyNode s, final PhylogenyNode g ) {
        boolean oyako = false;
        if ( ( g.getChildNode1().getLink() == s ) || ( g.getChildNode2().getLink() == s ) ) {
            oyako = true;
        }
        if ( g.getLink().getNumberOfDescendants() == 2 ) {
            if ( oyako ) {
                g.getNodeData().setEvent( createDuplicationEvent() );
            }
            else {
                g.getNodeData().setEvent( createSpeciationEvent() );
            }
        }
        else {
            if ( oyako ) {
                final Set<PhylogenyNode> set = new HashSet<PhylogenyNode>();
                for( PhylogenyNode n : g.getChildNode1().getAllExternalDescendants() ) {
                    n = n.getLink();
                    while ( n.getParent() != s ) {
                        n = n.getParent();
                        if ( n.isRoot() ) {
                            break;
                        }
                    }
                    set.add( n );
                }
                boolean multiple = false;
                for( PhylogenyNode n : g.getChildNode2().getAllExternalDescendants() ) {
                    n = n.getLink();
                    while ( n.getParent() != s ) {
                        n = n.getParent();
                        if ( n.isRoot() ) {
                            break;
                        }
                    }
                    if ( set.contains( n ) ) {
                        multiple = true;
                        break;
                    }
                }
                if ( multiple ) {
                    g.getNodeData().setEvent( createDuplicationEvent() );
                }
                else {
                    if ( _most_parsimonious_duplication_model ) {
                        g.getNodeData().setEvent( createSpeciationEvent() );
                    }
                    else {
                        g.getNodeData().setEvent( createSingleSpeciationOrDuplicationEvent() );
                    }
                }
            }
            else {
                g.getNodeData().setEvent( createSpeciationEvent() );
            }
        }
    }

    /**
     * Traverses the subtree of PhylogenyNode g in postorder, calculating the
     * mapping function M, and determines which nodes represent speciation
     * events and which ones duplication events.
     * <p>
     * Preconditions: Mapping M for external nodes must have been calculated and
     * the species tree must be labeled in preorder.
     * <p>
     * 
     */
    final void geneTreePostOrderTraversal() {
        for( final PhylogenyNodeIterator it = getGeneTree().iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode g = it.next();
            if ( g.isInternal() ) {
                PhylogenyNode s1 = g.getChildNode1().getLink();
                PhylogenyNode s2 = g.getChildNode2().getLink();
                while ( s1 != s2 ) {
                    if ( s1.getId() > s2.getId() ) {
                        s1 = s1.getParent();
                    }
                    else {
                        s2 = s2.getParent();
                    }
                }
                g.setLink( s1 );
                determineEvent( s1, g );
            }
        }
    }

    private final Event createDuplicationEvent() {
        final Event event = Event.createSingleDuplicationEvent();
        ++_duplications_sum;
        return event;
    }

    private final Event createSingleSpeciationOrDuplicationEvent() {
        final Event event = Event.createSingleSpeciationOrDuplicationEvent();
        ++_speciation_or_duplication_events_sum;
        return event;
    }

    private final Event createSpeciationEvent() {
        final Event event = Event.createSingleSpeciationEvent();
        ++_speciations_sum;
        return event;
    }

    public final int getSpeciationOrDuplicationEventsSum() {
        return _speciation_or_duplication_events_sum;
    }

    public final int getSpeciationsSum() {
        return _speciations_sum;
    }

    /**
     * This allows for linking of internal nodes of the species tree (as opposed
     * to just external nodes, as in the method it overrides.
     * @throws SDIException 
     * 
     */
    @Override
    final void linkNodesOfG() throws SDIException {
        final Map<String, PhylogenyNode> species_to_node_map = new HashMap<String, PhylogenyNode>();
        final List<PhylogenyNode> species_tree_ext_nodes = new ArrayList<PhylogenyNode>();
        _tax_comp_base = determineTaxonomyComparisonBase( _gene_tree );
        // Stringyfied taxonomy is the key, node is the value.
        for( final PhylogenyNodeIterator iter = _species_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode s = iter.next();
            species_tree_ext_nodes.add( s );
            if ( s.getNodeData().isHasTaxonomy() ) {
                final String tax_str = taxonomyToString( s, _tax_comp_base );
                if ( !ForesterUtil.isEmpty( tax_str ) ) {
                    if ( species_to_node_map.containsKey( tax_str ) ) {
                        throw new SDIException( "taxonomy \"" + s + "\" is not unique in species tree" );
                    }
                    species_to_node_map.put( tax_str, s );
                }
            }
        }
        // Retrieve the reference to the node with a matching stringyfied taxonomy.
        for( final PhylogenyNodeIterator iter = _gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode g = iter.next();
            if ( !g.getNodeData().isHasTaxonomy() ) {
                if ( _strip_gene_tree ) {
                    _stripped_gene_tree_nodes.add( g );
                }
                else {
                    throw new SDIException( "gene tree node \"" + g + "\" has no taxonomic data" );
                }
            }
            else {
                final String tax_str = taxonomyToString( g, _tax_comp_base );
                if ( ForesterUtil.isEmpty( tax_str ) ) {
                    if ( _strip_gene_tree ) {
                        _stripped_gene_tree_nodes.add( g );
                    }
                    else {
                        throw new SDIException( "gene tree node \"" + g + "\" has no appropriate taxonomic data" );
                    }
                }
                else {
                    PhylogenyNode s = species_to_node_map.get( tax_str );
                    if ( ( _tax_comp_base == TaxonomyComparisonBase.SCIENTIFIC_NAME ) && ( s == null )
                            && ( ForesterUtil.countChars( tax_str, ' ' ) > 1 ) ) {
                        s = tryMapByRemovingOverlySpecificData( species_to_node_map, tax_str );
                    }
                    if ( s == null ) {
                        if ( _strip_gene_tree ) {
                            _stripped_gene_tree_nodes.add( g );
                        }
                        else {
                            throw new SDIException( "taxonomy \"" + g.getNodeData().getTaxonomy()
                                    + "\" not present in species tree" );
                        }
                    }
                    else {
                        g.setLink( s );
                        _mapped_species_tree_nodes.add( s );
                    }
                }
            }
        } // for loop
        if ( _strip_gene_tree ) {
            stripGeneTree();
            if ( getGeneTree().isEmpty() || ( getGeneTree().getNumberOfExternalNodes() < 2 ) ) {
                throw new SDIException( "species could not be mapped between gene tree and species tree" );
            }
        }
        if ( _strip_species_tree ) {
            stripSpeciesTree( species_tree_ext_nodes );
        }
    }

    private final PhylogenyNode tryMapByRemovingOverlySpecificData( final Map<String, PhylogenyNode> species_to_node_map,
                                                                    final String tax_str ) {
        PhylogenyNode s = tryMapByRemovingOverlySpecificData( species_to_node_map, tax_str, " (" );
        if ( s == null ) {
            if ( ForesterUtil.countChars( tax_str, ' ' ) == 2 ) {
                final String new_tax_str = tax_str.substring( 0, tax_str.lastIndexOf( ' ' ) ).trim();
                s = species_to_node_map.get( new_tax_str );
                if ( s != null ) {
                    addScientificNamesMappedToReducedSpecificity( tax_str, new_tax_str );
                }
            }
        }
        if ( s == null ) {
            for( final String t : new String[] { " subspecies ", " strain ", " variety ", " varietas ", " subvariety ",
                    " form ", " subform ", " cultivar ", " section ", " subsection " } ) {
                s = tryMapByRemovingOverlySpecificData( species_to_node_map, tax_str, t );
                if ( s != null ) {
                    break;
                }
            }
        }
        return s;
    }

    private final PhylogenyNode tryMapByRemovingOverlySpecificData( final Map<String, PhylogenyNode> species_to_node_map,
                                                                    final String tax_str,
                                                                    final String term ) {
        final int i = tax_str.indexOf( term );
        if ( i > 4 ) {
            final String new_tax_str = tax_str.substring( 0, i ).trim();
            final PhylogenyNode s = species_to_node_map.get( new_tax_str );
            if ( s != null ) {
                addScientificNamesMappedToReducedSpecificity( tax_str, new_tax_str );
            }
            return s;
        }
        return null;
    }

    private final void addScientificNamesMappedToReducedSpecificity( final String s1, final String s2 ) {
        _scientific_names_mapped_to_reduced_specificity.add( s1 + " -> " + s2 );
    }

    public final SortedSet<String> getReMappedScientificNamesFromGeneTree() {
        return _scientific_names_mapped_to_reduced_specificity;
    }

    public TaxonomyComparisonBase getTaxCompBase() {
        return _tax_comp_base;
    }

    private void stripSpeciesTree( final List<PhylogenyNode> species_tree_ext_nodes ) {
        for( final PhylogenyNode s : species_tree_ext_nodes ) {
            if ( !_mapped_species_tree_nodes.contains( s ) ) {
                _species_tree.deleteSubtree( s, true );
                _stripped_species_tree_nodes.add( s );
            }
        }
        _species_tree.clearHashIdToNodeMap();
        _species_tree.externalNodesHaveChanged();
    }

    public List<PhylogenyNode> getStrippedSpeciesTreeNodes() {
        return _stripped_species_tree_nodes;
    }

    private void stripGeneTree() {
        for( final PhylogenyNode g : _stripped_gene_tree_nodes ) {
            _gene_tree.deleteSubtree( g, true );
        }
        _gene_tree.clearHashIdToNodeMap();
        _gene_tree.externalNodesHaveChanged();
    }

    public Set<PhylogenyNode> getMappedExternalSpeciesTreeNodes() {
        return _mapped_species_tree_nodes;
    }

    public static TaxonomyComparisonBase determineTaxonomyComparisonBase( final Phylogeny gene_tree ) {
        int with_id_count = 0;
        int with_code_count = 0;
        int with_sn_count = 0;
        int max = 0;
        for( final PhylogenyNodeIterator iter = gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode g = iter.next();
            if ( g.getNodeData().isHasTaxonomy() ) {
                final Taxonomy tax = g.getNodeData().getTaxonomy();
                if ( ( tax.getIdentifier() != null ) && !ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) {
                    if ( ++with_id_count > max ) {
                        max = with_id_count;
                    }
                }
                if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                    if ( ++with_code_count > max ) {
                        max = with_code_count;
                    }
                }
                if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                    if ( ++with_sn_count > max ) {
                        max = with_sn_count;
                    }
                }
            }
        }
        if ( max == 0 ) {
            throw new IllegalArgumentException( "gene tree has no taxonomic data" );
        }
        else if ( max == 1 ) {
            throw new IllegalArgumentException( "gene tree has only one node with taxonomic data" );
        }
        else if ( max == with_id_count ) {
            return SDI.TaxonomyComparisonBase.ID;
        }
        else if ( max == with_sn_count ) {
            return SDI.TaxonomyComparisonBase.SCIENTIFIC_NAME;
        }
        else {
            return SDI.TaxonomyComparisonBase.CODE;
        }
    }

    public List<PhylogenyNode> getStrippedExternalGeneTreeNodes() {
        return _stripped_gene_tree_nodes;
    }

    @Override
    public final String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( "Most parsimonious duplication model: " + _most_parsimonious_duplication_model );
        sb.append( ForesterUtil.getLineSeparator() );
        sb.append( "Speciations sum                    : " + getSpeciationsSum() );
        sb.append( ForesterUtil.getLineSeparator() );
        sb.append( "Duplications sum                   : " + getDuplicationsSum() );
        sb.append( ForesterUtil.getLineSeparator() );
        if ( !_most_parsimonious_duplication_model ) {
            sb.append( "Speciation or duplications sum     : " + getSpeciationOrDuplicationEventsSum() );
            sb.append( ForesterUtil.getLineSeparator() );
        }
        sb.append( "mapping cost L                     : " + computeMappingCostL() );
        return sb.toString();
    }
}
