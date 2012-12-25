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
// WWW: www.phylosoft.org/forester

package org.forester.io.parsers.phyloxml;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import org.forester.io.parsers.util.ParserUtils;

public final class PhyloXmlUtil {

    public static final String       OTHER                                      = "other";
    public static final String       UNKNOWN                                    = "unknown";
    public final static Pattern      SEQUENCE_SYMBOL_PATTERN                    = Pattern.compile( "\\S{1,20}" );
    public final static Pattern      TAXOMONY_CODE_PATTERN_STRICT               = ParserUtils.TAXOMONY_CODE_PATTERN_1;
    public final static Pattern      TAXOMONY_CODE_PATTERN_LAX                  = Pattern.compile( "[A-Z0-9]{3,6}" );
    public final static Pattern      LIT_REF_DOI_PATTERN                        = Pattern
                                                                                        .compile( "[a-zA-Z0-9_\\.]+\\S+" );
    public final static Set<String>  SEQUENCE_TYPES                             = new HashSet<String>();
    public final static List<String> TAXONOMY_RANKS_LIST                        = new ArrayList<String>();
    public final static Set<String>  TAXONOMY_RANKS_SET                         = new HashSet<String>();
    public static final int          ROUNDING_DIGITS_FOR_PHYLOXML_DOUBLE_OUTPUT = 9;
    public static final String       VECTOR_PROPERTY_REF                        = "vector:index=";
    public static final String       VECTOR_PROPERTY_TYPE                       = "xsd:decimal";
    public static final String       UNIPROT_TAX_PROVIDER                       = "uniprot";
    public static final String       SEQ_TYPE_RNA                               = "rna";
    public static final String       SEQ_TYPE_DNA                               = "dna";
    public static final String       SEQ_TYPE_PROTEIN                           = "protein";
    static {
        SEQUENCE_TYPES.add( SEQ_TYPE_RNA );
        SEQUENCE_TYPES.add( SEQ_TYPE_PROTEIN );
        SEQUENCE_TYPES.add( SEQ_TYPE_DNA );
        TAXONOMY_RANKS_LIST.add( "domain" );
        TAXONOMY_RANKS_LIST.add( "superkingdom" );
        TAXONOMY_RANKS_LIST.add( "kingdom" );
        TAXONOMY_RANKS_LIST.add( "subkingdom" );
        TAXONOMY_RANKS_LIST.add( "branch" );
        TAXONOMY_RANKS_LIST.add( "infrakingdom" );
        TAXONOMY_RANKS_LIST.add( "superphylum" );
        TAXONOMY_RANKS_LIST.add( "phylum" );
        TAXONOMY_RANKS_LIST.add( "subphylum" );
        TAXONOMY_RANKS_LIST.add( "infraphylum" );
        TAXONOMY_RANKS_LIST.add( "microphylum" );
        TAXONOMY_RANKS_LIST.add( "superdivision" );
        TAXONOMY_RANKS_LIST.add( "division" );
        TAXONOMY_RANKS_LIST.add( "subdivision" );
        TAXONOMY_RANKS_LIST.add( "infradivision" );
        TAXONOMY_RANKS_LIST.add( "superclass" );
        TAXONOMY_RANKS_LIST.add( "class" );
        TAXONOMY_RANKS_LIST.add( "subclass" );
        TAXONOMY_RANKS_LIST.add( "infraclass" );
        TAXONOMY_RANKS_LIST.add( "superlegion" );
        TAXONOMY_RANKS_LIST.add( "legion" );
        TAXONOMY_RANKS_LIST.add( "sublegion" );
        TAXONOMY_RANKS_LIST.add( "infralegion" );
        TAXONOMY_RANKS_LIST.add( "supercohort" );
        TAXONOMY_RANKS_LIST.add( "cohort" );
        TAXONOMY_RANKS_LIST.add( "subcohort" );
        TAXONOMY_RANKS_LIST.add( "infracohort" );
        TAXONOMY_RANKS_LIST.add( "superorder" );
        TAXONOMY_RANKS_LIST.add( "order" );
        TAXONOMY_RANKS_LIST.add( "suborder" );
        TAXONOMY_RANKS_LIST.add( "infraorder" );
        TAXONOMY_RANKS_LIST.add( "superfamily" );
        TAXONOMY_RANKS_LIST.add( "family" );
        TAXONOMY_RANKS_LIST.add( "subfamily" );
        TAXONOMY_RANKS_LIST.add( "supertribe" );
        TAXONOMY_RANKS_LIST.add( "tribe" );
        TAXONOMY_RANKS_LIST.add( "subtribe" );
        TAXONOMY_RANKS_LIST.add( "infratribe" );
        TAXONOMY_RANKS_LIST.add( "genus" );
        TAXONOMY_RANKS_LIST.add( "subgenus" );
        TAXONOMY_RANKS_LIST.add( "superspecies" );
        TAXONOMY_RANKS_LIST.add( "species" );
        TAXONOMY_RANKS_LIST.add( "subspecies" );
        TAXONOMY_RANKS_LIST.add( "variety" );
        TAXONOMY_RANKS_LIST.add( "varietas" );
        TAXONOMY_RANKS_LIST.add( "subvariety" );
        TAXONOMY_RANKS_LIST.add( "form" );
        TAXONOMY_RANKS_LIST.add( "subform" );
        TAXONOMY_RANKS_LIST.add( "cultivar" );
        TAXONOMY_RANKS_LIST.add( "strain" );
        TAXONOMY_RANKS_LIST.add( "section" );
        TAXONOMY_RANKS_LIST.add( "subsection" );
        TAXONOMY_RANKS_LIST.add( UNKNOWN );
        TAXONOMY_RANKS_LIST.add( OTHER );
        // same thing as set:
        TAXONOMY_RANKS_SET.add( "domain" );
        TAXONOMY_RANKS_SET.add( "superkingdom" );
        TAXONOMY_RANKS_SET.add( "kingdom" );
        TAXONOMY_RANKS_SET.add( "subkingdom" );
        TAXONOMY_RANKS_SET.add( "branch" );
        TAXONOMY_RANKS_SET.add( "infrakingdom" );
        TAXONOMY_RANKS_SET.add( "superphylum" );
        TAXONOMY_RANKS_SET.add( "phylum" );
        TAXONOMY_RANKS_SET.add( "subphylum" );
        TAXONOMY_RANKS_SET.add( "infraphylum" );
        TAXONOMY_RANKS_SET.add( "microphylum" );
        TAXONOMY_RANKS_SET.add( "superdivision" );
        TAXONOMY_RANKS_SET.add( "division" );
        TAXONOMY_RANKS_SET.add( "subdivision" );
        TAXONOMY_RANKS_SET.add( "infradivision" );
        TAXONOMY_RANKS_SET.add( "superclass" );
        TAXONOMY_RANKS_SET.add( "class" );
        TAXONOMY_RANKS_SET.add( "subclass" );
        TAXONOMY_RANKS_SET.add( "infraclass" );
        TAXONOMY_RANKS_SET.add( "superlegion" );
        TAXONOMY_RANKS_SET.add( "legion" );
        TAXONOMY_RANKS_SET.add( "sublegion" );
        TAXONOMY_RANKS_SET.add( "infralegion" );
        TAXONOMY_RANKS_SET.add( "supercohort" );
        TAXONOMY_RANKS_SET.add( "cohort" );
        TAXONOMY_RANKS_SET.add( "subcohort" );
        TAXONOMY_RANKS_SET.add( "infracohort" );
        TAXONOMY_RANKS_SET.add( "superorder" );
        TAXONOMY_RANKS_SET.add( "order" );
        TAXONOMY_RANKS_SET.add( "suborder" );
        TAXONOMY_RANKS_SET.add( "infraorder" );
        TAXONOMY_RANKS_SET.add( "superfamily" );
        TAXONOMY_RANKS_SET.add( "family" );
        TAXONOMY_RANKS_SET.add( "subfamily" );
        TAXONOMY_RANKS_SET.add( "supertribe" );
        TAXONOMY_RANKS_SET.add( "tribe" );
        TAXONOMY_RANKS_SET.add( "subtribe" );
        TAXONOMY_RANKS_SET.add( "infratribe" );
        TAXONOMY_RANKS_SET.add( "genus" );
        TAXONOMY_RANKS_SET.add( "subgenus" );
        TAXONOMY_RANKS_SET.add( "superspecies" );
        TAXONOMY_RANKS_SET.add( "species" );
        TAXONOMY_RANKS_SET.add( "subspecies" );
        TAXONOMY_RANKS_SET.add( "variety" );
        TAXONOMY_RANKS_SET.add( "varietas" );
        TAXONOMY_RANKS_SET.add( "subvariety" );
        TAXONOMY_RANKS_SET.add( "form" );
        TAXONOMY_RANKS_SET.add( "subform" );
        TAXONOMY_RANKS_SET.add( "cultivar" );
        TAXONOMY_RANKS_SET.add( "strain" );
        TAXONOMY_RANKS_SET.add( "section" );
        TAXONOMY_RANKS_SET.add( "subsection" );
        TAXONOMY_RANKS_SET.add( UNKNOWN );
        TAXONOMY_RANKS_SET.add( OTHER );
    };
}
