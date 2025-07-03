// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2000-2009 Christian M. Zmasek
// Copyright (C) 2007-2009 Burnham Institute for Medical Research
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

package org.forester.util;

public final class ForesterConstants {

    public final static String NH_COMMENT                          = "nh:comment";
    public final static String  FORESTER_VERSION                 = "1.051";
    public final static String  FORESTER_DATE                    = "190425";
    public final static String  PHYLO_XML_VERSION                = "1.20";
    public final static String  PHYLO_XML_LOCATION               = "http://www.phyloxml.org";
    public final static String  PHYLO_XML_XSD                    = "phyloxml.xsd";
    public final static String  XML_SCHEMA_INSTANCE              = "http://www.w3.org/2001/XMLSchema-instance";
    public final static String  LOCAL_PHYLOXML_XSD_RESOURCE      = "resources/phyloxml.xsd";
    public final static String  PHYLO_XML_SUFFIX                 = ".xml";
    public final static String  ID_NORMALIZED_FASTA_FILE_SUFFIX  = "_ni.fasta";
    public final static String  ID_NORMALIZED_NEXUS_FILE_SUFFIX  = "_ni.nexus";
    public final static String  ID_NORMALIZED_PHYLIP_FILE_SUFFIX = "_ni.phylip";
    public final static String  FASTA_FILE_SUFFIX                = ".fasta";
    public final static String  NEXUS_FILE_SUFFIX                = ".nexus";
    public final static String  PHYLIP_FILE_SUFFIX               = ".phylip";
    public final static String  ID_MAP_FILE_SUFFIX               = ".nim";
    public final static String  UTF_8                            = "UTF-8";
    public final static String  ISO_8859_1                       = "ISO-8859-1";
    public final static String GO_LINK                        = "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&query=";
    public final static String PFAM_FAMILY_ID_LINK            = "http://pfam.xfam.org/family/";

    
    public final static String  PHYLO_XML_REFERENCE              = "Han MV and Zmasek CM (2009): \"phyloXML: XML for evolutionary biology and comparative genomics\", BMC Bioinformatics 10:356";
    public final static boolean RELEASE                          = false;

    public enum PhylogeneticTreeFormats {
                                         NH,
                                         NHX,
                                         NEXUS,
                                         PHYLOXML
    }

    public static final String UNIPROT_TAXONOMY_ID_LINK       = "http://www.uniprot.org/taxonomy/";
    public static final String EOL_LINK                       = "http://www.eol.org/search?q=";
    public static final String GOOGLE_SCHOLAR_SEARCH          = "http://scholar.google.com/scholar?q=";
    public static final String GOOGLE_WEB_SEARCH_LINK         = "http://www.google.com/search?q=";
}
