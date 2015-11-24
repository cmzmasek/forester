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
