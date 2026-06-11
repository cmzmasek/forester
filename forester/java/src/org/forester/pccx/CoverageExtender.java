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

package org.forester.pccx;

import java.io.PrintStream;
import java.util.List;

import org.forester.phylogeny.Phylogeny;

/*
 * @author Christian M. Zmasek
 */
public interface CoverageExtender {

    public abstract List<String> find( final List<Phylogeny> phylogenies,
                                       final List<String> already_covered,
                                       int number_names_to_find,
                                       final CoverageCalculationOptions options,
                                       final PrintStream out );
}