// $Id:
// forester -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2010 Christian M. Zmasek
// Copyright (C) 2008-2010 Burnham Institute for Medical Research
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

package org.forester.archaeopteryx.webservices;

// The interactive TreeBASE/TreeFam/Tree-of-Life tree readers were removed; these
// endpoint constants remain only for the URL-based tree-reading tests in org.forester.test.Test.
public final class WebserviceUtil {

    public static final String TOL_URL_BASE                    = "http://tolweb.org/onlinecontributors/app?service=external&page=xml/TreeStructureService&node_id=";
    public static final String TREE_FAM_URL_BASE               = "http://www.treefam.org/family/TF";
    public static final String TREEBASE_PHYLOWS_STUDY_URL_BASE = "https://treebase.org/treebase-web/phylows/study/TB2:S";
    public static final String TREEBASE_PHYLOWS_TREE_URL_BASE  = "https://treebase.org/treebase-web/phylows/tree/TB2:Tr";
}
