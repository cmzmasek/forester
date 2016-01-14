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

package org.forester.protein;

public interface BinaryDomainCombination extends Comparable<BinaryDomainCombination> {

    public static final String SEPARATOR = "=";

    public String getId0();

    public String getId1();

    short getId0Code();

    short getId1Code();

    public StringBuffer toGraphDescribingLanguage( final OutputFormat format,
                                                   final String node_attribute,
                                                   String edge_attribute );

    public static enum DomainCombinationType {
        BASIC, DIRECTED, DIRECTED_ADJACTANT;
    }

    public static enum OutputFormat {
        DOT
    }
}