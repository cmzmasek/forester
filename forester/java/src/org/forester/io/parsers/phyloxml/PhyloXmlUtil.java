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

package org.forester.io.parsers.phyloxml;

import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;

import org.forester.io.parsers.util.ParserUtils;

public final class PhyloXmlUtil {

    public final static Pattern      SEQUENCE_SYMBOL_PATTERN                    = Pattern.compile( "\\S{1,20}" );
    public final static Pattern      TAXOMONY_CODE_PATTERN                      = Pattern
            .compile( ParserUtils.TAX_CODE );
    public final static Pattern      LIT_REF_DOI_PATTERN                        = Pattern
            .compile( "[a-zA-Z0-9_\\.]+\\S+" );
    public final static Set<String>  SEQUENCE_TYPES                             = new HashSet<String>();
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
       
    };
}
