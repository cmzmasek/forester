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

package org.forester.go;

import java.util.List;

import org.forester.phylogeny.data.PhylogenyData;

public interface GoTerm extends PhylogenyData, Comparable<GoTerm> {

    public List<GoId> getAltIds();

    public String getComment();

    public String getDefinition();

    public GoId getGoId();

    public GoNameSpace getGoNameSpace();

    public List<GoRelationship> getGoRelationships();

    public List<GoSubset> getGoSubsets();

    public List<GoXRef> getGoXRefs();

    public String getName();

    public List<GoId> getSuperGoIds();

    public boolean isObsolete();
}
