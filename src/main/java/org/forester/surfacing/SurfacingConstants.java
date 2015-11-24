// $Id:
//
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

package org.forester.surfacing;

import org.forester.util.ForesterUtil;

public class SurfacingConstants {

    public static final String AMIGO_LINK                     = "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&query=";
    public static final String EOL_LINK                       = "http://www.eol.org/search?q=";
    public static final String GO_LINK                        = "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&query=";
    public static final String GOOGLE_SCHOLAR_SEARCH          = "http://scholar.google.com/scholar?q=";
    public static final String GOOGLE_WEB_SEARCH_LINK         = "http://www.google.com/search?q=";
    public static final String NL                             = ForesterUtil.LINE_SEPARATOR;
    public static final String NONE                           = "[none]";
    public static final String PFAM_FAMILY_ID_LINK            = "http://pfam.xfam.org/family/";
    public static final String UNIPROT_TAXONOMY_ID_LINK       = "http://www.uniprot.org/taxonomy/";
    static final boolean       PRINT_MORE_DOM_SIMILARITY_INFO = false;
    static final boolean       SECONDARY_FEATURES_ARE_SCOP    = true;
    static final String        SECONDARY_FEATURES_SCOP_LINK   = "http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?key=";
}
