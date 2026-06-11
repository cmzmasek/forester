// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

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

    final Set<PhylogenyNode> getMappedSpeciesTreeNodes() {
        return _mapped_species_tree_nodes;
    }

    final SortedSet<String> getScientificNamesMappedToReducedSpecificity() {
        return _scientific_names_mapped_to_reduced_specificity;
    }

    final List<PhylogenyNode> getStrippedGeneTreeNodes() {
        return _stripped_gene_tree_nodes;
    }

    final List<PhylogenyNode> getStrippedSpeciesTreeNodes() {
        return _stripped_species_tree_nodes;
    }

    final TaxonomyComparisonBase getTaxCompBase() {
        return _tax_comp_base;
    }

    final void setTaxCompBase( final TaxonomyComparisonBase tax_comp_base ) {
        _tax_comp_base = tax_comp_base;
    }
}
