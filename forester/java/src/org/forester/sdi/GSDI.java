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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.sdi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;
import org.forester.util.ForesterUtil;

public final class GSDI implements GSDII {

    private final boolean                _most_parsimonious_duplication_model;
    private final int                    _speciation_or_duplication_events_sum;
    private final int                    _speciations_sum;
    private final int                    _duplications_sum;
    private final List<PhylogenyNode>    _stripped_gene_tree_nodes;
    private final List<PhylogenyNode>    _stripped_species_tree_nodes;
    private final Set<PhylogenyNode>     _mapped_species_tree_nodes;
    private final TaxonomyComparisonBase _tax_comp_base;
    private final SortedSet<String>      _scientific_names_mapped_to_reduced_specificity;

    public GSDI( final Phylogeny gene_tree,
                 final Phylogeny species_tree,
                 final boolean most_parsimonious_duplication_model,
                 final boolean strip_gene_tree,
                 final boolean strip_species_tree ) throws SDIException {
        _most_parsimonious_duplication_model = most_parsimonious_duplication_model;
        if ( gene_tree.getRoot().getNumberOfDescendants() == 3 ) {
            gene_tree.reRoot( gene_tree.getRoot().getChildNode( 2 ) );
        }
        final NodesLinkingResult nodes_linking_result = linkNodesOfG( gene_tree,
                                                                      species_tree,
                                                                      strip_gene_tree,
                                                                      strip_species_tree );
        _stripped_gene_tree_nodes = nodes_linking_result.getStrippedGeneTreeNodes();
        _stripped_species_tree_nodes = nodes_linking_result.getStrippedSpeciesTreeNodes();
        _mapped_species_tree_nodes = nodes_linking_result.getMappedSpeciesTreeNodes();
        _scientific_names_mapped_to_reduced_specificity = nodes_linking_result
                .getScientificNamesMappedToReducedSpecificity();
        _tax_comp_base = nodes_linking_result.getTaxCompBase();
        PhylogenyMethods.preOrderReId( species_tree );
        final GSDIsummaryResult gsdi_summary_result = geneTreePostOrderTraversal( gene_tree,
                                                                                  _most_parsimonious_duplication_model );
        _speciation_or_duplication_events_sum = gsdi_summary_result.getSpeciationOrDuplicationEventsSum();
        _speciations_sum = gsdi_summary_result.getSpeciationsSum();
        _duplications_sum = gsdi_summary_result.getDuplicationsSum();
    }

    public int getDuplicationsSum() {
        return _duplications_sum;
    }

    @Override
    public Set<PhylogenyNode> getMappedExternalSpeciesTreeNodes() {
        return _mapped_species_tree_nodes;
    }

    @Override
    public final SortedSet<String> getReMappedScientificNamesFromGeneTree() {
        return _scientific_names_mapped_to_reduced_specificity;
    }

    public final int getSpeciationOrDuplicationEventsSum() {
        return _speciation_or_duplication_events_sum;
    }

    @Override
    public final int getSpeciationsSum() {
        return _speciations_sum;
    }

    @Override
    public List<PhylogenyNode> getStrippedExternalGeneTreeNodes() {
        return _stripped_gene_tree_nodes;
    }

    @Override
    public List<PhylogenyNode> getStrippedSpeciesTreeNodes() {
        return _stripped_species_tree_nodes;
    }

    @Override
    public TaxonomyComparisonBase getTaxCompBase() {
        return _tax_comp_base;
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
        return sb.toString();
    }

    /**
     * Traverses the subtree of PhylogenyNode g in postorder, calculating the
     * mapping function M, and determines which nodes represent speciation
     * events and which ones duplication events.
     * <p>
     * Preconditions: Mapping M for external nodes must have been calculated and
     * the species tree must be labeled in preorder.
     * <p>
     * @return 
     * @throws SDIException 
     * 
     */
    final static GSDIsummaryResult geneTreePostOrderTraversal( final Phylogeny gene_tree,
                                                               final boolean most_parsimonious_duplication_model )
            throws SDIException {
        final GSDIsummaryResult res = new GSDIsummaryResult();
        for( final PhylogenyNodeIterator it = gene_tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode g = it.next();
            if ( g.isInternal() ) {
                if ( g.getNumberOfDescendants() != 2 ) {
                    throw new SDIException( "gene tree contains internal node with " + g.getNumberOfDescendants()
                            + " descendents" );
                }
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
                determineEvent( s1, g, most_parsimonious_duplication_model, res );
            }
        }
        return res;
    }

    final static GSDIsummaryResult geneTreePostOrderTraversal( final Phylogeny gene_tree,
                                                               final boolean most_parsimonious_duplication_model,
                                                               final int min_duplications ) throws SDIException {
        final GSDIsummaryResult res = new GSDIsummaryResult();
        for( final PhylogenyNodeIterator it = gene_tree.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode g = it.next();
            if ( g.isInternal() ) {
                if ( g.getNumberOfDescendants() != 2 ) {
                    throw new SDIException( "gene tree contains internal node with " + g.getNumberOfDescendants()
                            + " descendents" );
                }
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
                determineEvent( s1, g, most_parsimonious_duplication_model, res );
                if ( res.getDuplicationsSum() > min_duplications ) {
                    return null;
                }
            }
        }
        return res;
    }

    final static NodesLinkingResult linkNodesOfG( final Phylogeny gene_tree,
                                                  final Phylogeny species_tree,
                                                  final boolean strip_gene_tree,
                                                  final boolean strip_species_tree ) throws SDIException {
        final TaxonomyComparisonBase tax_comp_base = SDIutil.determineTaxonomyComparisonBase( gene_tree );
        if ( tax_comp_base == null ) {
            throw new RuntimeException( "failed to establish taxonomy linking base (taxonomy linking base is null)" );
        }
        return linkNodesOfG( gene_tree, species_tree, tax_comp_base, strip_gene_tree, strip_species_tree );
    }

    /**
     * This allows for linking of internal nodes of the species tree (as opposed
     * to just external nodes, as in the method it overrides.
     * If TaxonomyComparisonBase is null, it will try to determine it.
     * @throws SDIException 
     * 
     */
    final static NodesLinkingResult linkNodesOfG( final Phylogeny gene_tree,
                                                  final Phylogeny species_tree,
                                                  final TaxonomyComparisonBase tax_comp_base,
                                                  final boolean strip_gene_tree,
                                                  final boolean strip_species_tree ) throws SDIException {
        if ( tax_comp_base == null ) {
            throw new IllegalArgumentException( "taxonomy linking base is null" );
        }
        final Map<String, PhylogenyNode> species_to_node_map = new HashMap<String, PhylogenyNode>();
        final List<PhylogenyNode> species_tree_ext_nodes = new ArrayList<PhylogenyNode>();
        final NodesLinkingResult res = new NodesLinkingResult();
        res.setTaxCompBase( tax_comp_base );
        // Stringyfied taxonomy is the key, node is the value.
        for( final PhylogenyNodeIterator iter = species_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode s = iter.next();
            species_tree_ext_nodes.add( s );
            if ( s.getNodeData().isHasTaxonomy() ) {
                final String tax_str = SDIutil.taxonomyToString( s, res.getTaxCompBase() );
                if ( !ForesterUtil.isEmpty( tax_str ) ) {
                    if ( species_to_node_map.containsKey( tax_str ) ) {
                        throw new SDIException( "taxonomy \"" + tax_str + "\" is not unique in species tree (using "
                                + res.getTaxCompBase() + " for linking to gene tree)" );
                    }
                    species_to_node_map.put( tax_str, s );
                }
            }
        }
        // Retrieve the reference to the node with a matching stringyfied taxonomy.
        for( final PhylogenyNodeIterator iter = gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode g = iter.next();
            if ( !g.getNodeData().isHasTaxonomy() ) {
                if ( strip_gene_tree ) {
                    res.getStrippedGeneTreeNodes().add( g );
                }
                else {
                    throw new SDIException( "gene tree node \"" + g + "\" has no taxonomic data" );
                }
            }
            else {
                final String tax_str = SDIutil.taxonomyToString( g, res.getTaxCompBase() );
                if ( ForesterUtil.isEmpty( tax_str ) ) {
                    if ( strip_gene_tree ) {
                        res.getStrippedGeneTreeNodes().add( g );
                    }
                    else {
                        throw new SDIException( "gene tree node \"" + g + "\" has no appropriate taxonomic data" );
                    }
                }
                else {
                    PhylogenyNode s = species_to_node_map.get( tax_str );
                    if ( ( res.getTaxCompBase() == TaxonomyComparisonBase.SCIENTIFIC_NAME ) && ( s == null )
                            && ( ForesterUtil.countChars( tax_str, ' ' ) > 1 ) ) {
                        s = tryMapByRemovingOverlySpecificData( species_to_node_map,
                                                                tax_str,
                                                                res.getScientificNamesMappedToReducedSpecificity() );
                    }
                    if ( s == null ) {
                        if ( strip_gene_tree ) {
                            res.getStrippedGeneTreeNodes().add( g );
                        }
                        else {
                            throw new SDIException( "taxonomy \"" + g.getNodeData().getTaxonomy()
                                    + "\" not present in species tree" );
                        }
                    }
                    else {
                        g.setLink( s );
                        res.getMappedSpeciesTreeNodes().add( s );
                    }
                }
            }
        } // for loop
        if ( strip_gene_tree ) {
            stripTree( gene_tree, res.getStrippedGeneTreeNodes() );
            if ( gene_tree.isEmpty() || ( gene_tree.getNumberOfExternalNodes() < 2 ) ) {
                throw new SDIException( "species could not be mapped between gene tree and species tree (based on "
                        + res.getTaxCompBase() + ")" );
            }
        }
        if ( strip_species_tree ) {
            stripSpeciesTree( species_tree, species_tree_ext_nodes, res );
        }
        return res;
    }

    private final static void addScientificNamesMappedToReducedSpecificity( final String s1,
                                                                            final String s2,
                                                                            final SortedSet<String> scientific_names_mapped_to_reduced_specificity ) {
        scientific_names_mapped_to_reduced_specificity.add( s1 + " -> " + s2 );
    }

    private final static void determineEvent( final PhylogenyNode s,
                                              final PhylogenyNode g,
                                              final boolean most_parsimonious_duplication_model,
                                              final GSDIsummaryResult res ) {
        boolean oyako = false;
        if ( ( g.getChildNode1().getLink() == s ) || ( g.getChildNode2().getLink() == s ) ) {
            oyako = true;
        }
        if ( g.getLink().getNumberOfDescendants() == 2 ) {
            if ( oyako ) {
                g.getNodeData().setEvent( Event.createSingleDuplicationEvent() );
                res.increaseDuplicationsSum();
            }
            else {
                g.getNodeData().setEvent( Event.createSingleSpeciationEvent() );
                res.increaseSpeciationsSum();
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
                    g.getNodeData().setEvent( Event.createSingleDuplicationEvent() );
                    res.increaseDuplicationsSum();
                }
                else {
                    if ( most_parsimonious_duplication_model ) {
                        g.getNodeData().setEvent( Event.createSingleSpeciationEvent() );
                        res.increaseSpeciationsSum();
                    }
                    else {
                        g.getNodeData().setEvent( Event.createSingleSpeciationOrDuplicationEvent() );
                        res.increaseSpeciationOrDuplicationEventsSum();
                    }
                }
            }
            else {
                g.getNodeData().setEvent( Event.createSingleSpeciationEvent() );
                res.increaseSpeciationsSum();
            }
        }
    }

    private final static void stripSpeciesTree( final Phylogeny species_tree,
                                                final List<PhylogenyNode> species_tree_ext_nodes,
                                                final NodesLinkingResult res ) {
        for( final PhylogenyNode s : species_tree_ext_nodes ) {
            if ( !res.getMappedSpeciesTreeNodes().contains( s ) ) {
                species_tree.deleteSubtree( s, true );
                res.getStrippedSpeciesTreeNodes().add( s );
            }
        }
        species_tree.clearHashIdToNodeMap();
        species_tree.externalNodesHaveChanged();
    }

    private final static void stripTree( final Phylogeny phy, final List<PhylogenyNode> strip_nodes ) {
        for( final PhylogenyNode g : strip_nodes ) {
            phy.deleteSubtree( g, true );
        }
        phy.clearHashIdToNodeMap();
        phy.externalNodesHaveChanged();
    }

    private final static PhylogenyNode tryMapByRemovingOverlySpecificData( final Map<String, PhylogenyNode> species_to_node_map,
                                                                           final String tax_str,
                                                                           final SortedSet<String> scientific_names_mapped_to_reduced_specificity ) {
        PhylogenyNode s = tryMapByRemovingOverlySpecificData( species_to_node_map,
                                                              tax_str,
                                                              " (",
                                                              scientific_names_mapped_to_reduced_specificity );
        if ( s == null ) {
            if ( ForesterUtil.countChars( tax_str, ' ' ) == 2 ) {
                final String new_tax_str = tax_str.substring( 0, tax_str.lastIndexOf( ' ' ) ).trim();
                s = species_to_node_map.get( new_tax_str );
                if ( s != null ) {
                    addScientificNamesMappedToReducedSpecificity( tax_str,
                                                                  new_tax_str,
                                                                  scientific_names_mapped_to_reduced_specificity );
                }
            }
        }
        if ( s == null ) {
            for( final String t : new String[] { " subspecies ", " strain ", " variety ", " varietas ", " subvariety ",
                    " form ", " subform ", " cultivar ", " section ", " subsection " } ) {
                s = tryMapByRemovingOverlySpecificData( species_to_node_map,
                                                        tax_str,
                                                        t,
                                                        scientific_names_mapped_to_reduced_specificity );
                if ( s != null ) {
                    break;
                }
            }
        }
        return s;
    }

    private final static PhylogenyNode tryMapByRemovingOverlySpecificData( final Map<String, PhylogenyNode> species_to_node_map,
                                                                           final String tax_str,
                                                                           final String term,
                                                                           final SortedSet<String> scientific_names_mapped_to_reduced_specificity ) {
        final int i = tax_str.indexOf( term );
        if ( i > 4 ) {
            final String new_tax_str = tax_str.substring( 0, i ).trim();
            final PhylogenyNode s = species_to_node_map.get( new_tax_str );
            if ( s != null ) {
                addScientificNamesMappedToReducedSpecificity( tax_str,
                                                              new_tax_str,
                                                              scientific_names_mapped_to_reduced_specificity );
            }
            return s;
        }
        return null;
    }
}
