// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// Copyright (C) 2003-2007 Ethalinda K.S. Cannon
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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.phylogeny.data.Identifier;

public final class SequenceIdParser {

    // gb_ADF31344_1_segmented_worms_
    // gb_AAA96518_1
    // gb_EHB07727_1_rodents_
    // dbj_BAF37827_1_turtles_
    // emb_CAA73223_1_primates_
    // lcl_91970_unknown_
    // mites|ref_XP_002434188_1
    // ref_XP_002434188_1_mites___ticks_
    // ref_NP_001121530_1_frogs___toads_
    //The format for GenBank Accession numbers are:
    //Nucleotide: 1 letter + 5 numerals OR 2 letters + 6 numerals
    //Protein:    3 letters + 5 numerals
    //http://www.ncbi.nlm.nih.gov/Sequin/acc.html
    private final static Pattern GENBANK_NUCLEOTIDE_AC_PATTERN_1 = Pattern
                                                                         .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]\\d{5}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)" );
    private final static Pattern GENBANK_NUCLEOTIDE_AC_PATTERN_2 = Pattern
                                                                         .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]{2}\\d{6}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)" );
    private final static Pattern GENBANK_PROTEIN_AC_PATTERN      = Pattern
                                                                         .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]{3}\\d{5}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)" );
    // RefSeq accession numbers can be distinguished from GenBank accessions 
    // by their distinct prefix format of 2 characters followed by an
    // underscore character ('_'). For example, a RefSeq protein accession is NP_015325. 
    private final static Pattern REFSEQ_PATTERN                  = Pattern
                                                                         .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]{2}_\\d{6,})(?:[^a-zA-Z0-9]|\\Z)" );
    // See: http://web.expasy.org/docs/userman.html#ID_line
    private final static Pattern TREMBL_PATTERN                  = Pattern
                                                                         .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z][0-9][A-Z0-9]{3}[0-9])(?:[^a-zA-Z0-9]|\\Z)" );
    private final static Pattern GI_PATTERN                      = Pattern
                                                                         .compile( "(?:\\b|_)(?:GI|gi)[|_=:](\\d+)(?:\\b|_)" );

    /**
     * Returns null if no match.
     * 
     */
    public final static Identifier parse( final String s ) {
        String v = parseGenbankAccessor( s );
        if ( !ForesterUtil.isEmpty( v ) ) {
            return new Identifier( v, Identifier.NCBI );
        }
        v = parseRefSeqAccessor( s );
        if ( !ForesterUtil.isEmpty( v ) ) {
            return new Identifier( v, Identifier.REFSEQ );
        }
        v = parseTrEMBLAccessor( s );
        if ( !ForesterUtil.isEmpty( v ) ) {
            return new Identifier( v, Identifier.SP );
        }
        return null;
    }

    public final static boolean isProtein( final String query ) {
        final String r1 = parseRefSeqAccessor( query );
        if ( !ForesterUtil.isEmpty( r1 ) && ( r1.charAt( 1 ) == 'P' ) ) {
            return true;
        }
        final String r2 = parseTrEMBLAccessor( query );
        if ( !ForesterUtil.isEmpty( r2 ) ) {
            return true;
        }
        return GENBANK_PROTEIN_AC_PATTERN.matcher( query ).lookingAt();
    }

    /**
     * Returns null if no match.
     * 
     */
    public static String parseGenbankAccessor( final String query ) {
        Matcher m = GENBANK_NUCLEOTIDE_AC_PATTERN_1.matcher( query );
        if ( m.lookingAt() ) {
            return m.group( 1 );
        }
        else {
            m = GENBANK_NUCLEOTIDE_AC_PATTERN_2.matcher( query );
            if ( m.lookingAt() ) {
                return m.group( 1 );
            }
            else {
                m = GENBANK_PROTEIN_AC_PATTERN.matcher( query );
                if ( m.lookingAt() ) {
                    return m.group( 1 );
                }
                else {
                    return null;
                }
            }
        }
    }

    /**
     * Returns null if no match.
     * 
     */
    public final static String parseRefSeqAccessor( final String query ) {
        final Matcher m = REFSEQ_PATTERN.matcher( query );
        if ( m.lookingAt() ) {
            return m.group( 1 );
        }
        return null;
    }

    /**
     * Returns null if no match.
     * 
     */
    private final static String parseTrEMBLAccessor( final String query ) {
        final Matcher m = TREMBL_PATTERN.matcher( query );
        if ( m.lookingAt() ) {
            return m.group( 1 );
        }
        return null;
    }

    private SequenceIdParser() {
        // Hiding the constructor.
    }

    public static String parseGInumber( final String query ) {
        final Matcher m = GI_PATTERN.matcher( query );
        if ( m.find() ) {
            return m.group( 1 );
        }
        return null;
    }
}
