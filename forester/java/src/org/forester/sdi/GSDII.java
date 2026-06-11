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

import java.util.List;
import java.util.Set;
import java.util.SortedSet;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.sdi.SDIutil.TaxonomyComparisonBase;

public interface GSDII {

    public abstract Set<PhylogenyNode> getMappedExternalSpeciesTreeNodes();

    public abstract SortedSet<String> getReMappedScientificNamesFromGeneTree();

    public abstract int getSpeciationsSum();

    public abstract List<PhylogenyNode> getStrippedExternalGeneTreeNodes();

    public abstract List<PhylogenyNode> getStrippedSpeciesTreeNodes();

    public abstract TaxonomyComparisonBase getTaxCompBase();
}