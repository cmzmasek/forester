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

import java.util.List;

import org.forester.phylogeny.Phylogeny;

/*
 * @author Christian M. Zmasek
 */
public class CoverageCalculator {

    private final CoverageCalculationMethod  _method;
    private final CoverageCalculationOptions _options;

    private CoverageCalculator( final CoverageCalculationMethod method, final CoverageCalculationOptions options ) {
        _method = method;
        _options = options;
    }

    public Coverage calculateCoverage( final List<Phylogeny> phylogenies,
                                       final List<String> names,
                                       final boolean annotate_phylogenies ) {
        return getMethod().calculateCoverage( phylogenies, names, getOptions(), annotate_phylogenies );
    }

    private CoverageCalculationMethod getMethod() {
        return _method;
    }

    private CoverageCalculationOptions getOptions() {
        return _options;
    }

    public static CoverageCalculator getInstance( final CoverageCalculationMethod method,
                                                  final CoverageCalculationOptions options ) {
        return new CoverageCalculator( method, options );
    }
}
