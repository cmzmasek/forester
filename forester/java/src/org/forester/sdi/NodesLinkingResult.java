
package org.forester.sdi;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;

final class NodesLinkingResult {

    private final List<PhylogenyNode> _stripped_gene_tree_nodes;
    private final List<PhylogenyNode> _stripped_species_tree_nodes;
    private final Set<PhylogenyNode>  _mapped_species_tree_nodes;
    private TaxonomyComparisonBase    _tax_comp_base;
    private final SortedSet<String>   _scientific_names_mapped_to_reduced_specificity;

    NodesLinkingResult() {
        _stripped_gene_tree_nodes = new ArrayList<PhylogenyNode>();
        _stripped_species_tree_nodes = new ArrayList<PhylogenyNode>();
        _mapped_species_tree_nodes = new HashSet<PhylogenyNode>();
        _scientific_names_mapped_to_reduced_specificity = new TreeSet<String>();
        _tax_comp_base = null;
    }

    final List<PhylogenyNode> getStrippedGeneTreeNodes() {
        return _stripped_gene_tree_nodes;
    }

    final List<PhylogenyNode> getStrippedSpeciesTreeNodes() {
        return _stripped_species_tree_nodes;
    }

    final Set<PhylogenyNode> getMappedSpeciesTreeNodes() {
        return _mapped_species_tree_nodes;
    }

    final TaxonomyComparisonBase getTaxCompBase() {
        return _tax_comp_base;
    }

    final void setTaxCompBase( final TaxonomyComparisonBase tax_comp_base ) {
        _tax_comp_base = tax_comp_base;
    }

    final SortedSet<String> getScientificNamesMappedToReducedSpecificity() {
        return _scientific_names_mapped_to_reduced_specificity;
    }
}
