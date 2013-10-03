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

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Sequence;

public final class SequenceAccessionTools {

    public final static Pattern  UNIPROT_KB_PATTERN_0            = Pattern
                                                                         .compile( "\\b([A-Z][0-9][A-Z0-9]{3}[0-9])\\b" );
    public final static Pattern  UNIPROT_KB_PATTERN_1            = Pattern
                                                                         .compile( "(?:\\b|_)(?:sp|tr)[\\.|\\-_=/\\\\]([A-Z][0-9][A-Z0-9]{3}[0-9])(?:\\b|_)" );
    public final static Pattern  UNIPROT_KB_PATTERN_2            = Pattern
                                                                         .compile( "(?:\\b|_)(?:[A-Z0-9]{2,5}|(?:[A-Z][0-9][A-Z0-9]{3}[0-9]))_(([A-Z9][A-Z]{2}[A-Z0-9]{2})|RAT|PIG|PEA)(?:\\b|_)" );
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
    private final static Pattern GI_PATTERN                      = Pattern
                                                                         .compile( "(?:\\b|_)(?:GI|gi)[|_=:](\\d+)(?:\\b|_)" );
    // RefSeq accession numbers can be distinguished from GenBank accessions 
    // by their distinct prefix format of 2 characters followed by an
    // underscore character ('_'). For example, a RefSeq protein accession is NP_015325. 
    private final static Pattern REFSEQ_PATTERN                  = Pattern
                                                                         .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]{2}_\\d{6,})(?:[^a-zA-Z0-9]|\\Z)" );

    private SequenceAccessionTools() {
        // Hiding the constructor.
    }

    public static String extractGenbankAccessor( final PhylogenyNode node ) {
        String v = null;
        if ( node.getNodeData().isHasSequence() ) {
            final Sequence seq = node.getNodeData().getSequence();
            if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
                v = parseGenbankAccessor( seq.getSymbol() );
            }
            if ( !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
                v = parseGenbankAccessor( seq.getGeneName() );
            }
            if ( ForesterUtil.isEmpty( v ) && !ForesterUtil.isEmpty( seq.getName() ) ) {
                v = parseGenbankAccessor( seq.getName() );
            }
            if ( ForesterUtil.isEmpty( v ) && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
                v = parseGenbankAccessor( seq.getAccession().getValue() );
            }
        }
        if ( ForesterUtil.isEmpty( v ) && !ForesterUtil.isEmpty( node.getName() ) ) {
            v = parseGenbankAccessor( node.getName() );
        }
        return v;
    }

    public static String extractGInumber( final PhylogenyNode node ) {
        String v = null;
        if ( node.getNodeData().isHasSequence() ) {
            final Sequence seq = node.getNodeData().getSequence();
            if ( ForesterUtil.isEmpty( v ) && !ForesterUtil.isEmpty( seq.getName() ) ) {
                v = parseGInumber( seq.getName() );
            }
            if ( ForesterUtil.isEmpty( v ) && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
                v = parseGInumber( seq.getAccession().getValue() );
            }
        }
        if ( ForesterUtil.isEmpty( v ) && !ForesterUtil.isEmpty( node.getName() ) ) {
            v = parseGInumber( node.getName() );
        }
        return v;
    }

    public static String extractRefSeqAccessor( final PhylogenyNode node ) {
        String v = null;
        if ( node.getNodeData().isHasSequence() ) {
            final Sequence seq = node.getNodeData().getSequence();
            if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
                v = parseRefSeqAccessor( seq.getSymbol() );
            }
            if ( !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
                v = parseRefSeqAccessor( seq.getGeneName() );
            }
            if ( ForesterUtil.isEmpty( v ) && !ForesterUtil.isEmpty( seq.getName() ) ) {
                v = parseRefSeqAccessor( seq.getName() );
            }
            if ( ForesterUtil.isEmpty( v ) && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
                v = parseRefSeqAccessor( seq.getAccession().getValue() );
            }
        }
        if ( ForesterUtil.isEmpty( v ) && !ForesterUtil.isEmpty( node.getName() ) ) {
            v = parseRefSeqAccessor( node.getName() );
        }
        return v;
    }

    public static String extractUniProtKbProteinSeqIdentifier( final PhylogenyNode node ) {
        String a = null;
        if ( node.getNodeData().isHasSequence() ) {
            final Sequence seq = node.getNodeData().getSequence();
            if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
                a = SequenceAccessionTools.extractUniProtKbProteinSeqIdentifier( seq.getSymbol() );
            }
            if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( seq.getName() ) ) {
                a = SequenceAccessionTools.extractUniProtKbProteinSeqIdentifier( seq.getName() );
            }
            if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
                a = SequenceAccessionTools.extractUniProtKbProteinSeqIdentifier( seq.getGeneName() );
            }
            if ( ForesterUtil.isEmpty( a ) && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
                a = SequenceAccessionTools.extractUniProtKbProteinSeqIdentifier( seq.getAccession().getValue() );
            }
        }
        if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( node.getName() ) ) {
            a = SequenceAccessionTools.extractUniProtKbProteinSeqIdentifier( node.getName() );
        }
        return a;
    }

    public static String extractUniProtKbProteinSeqIdentifier( final String str ) {
        Matcher m = UNIPROT_KB_PATTERN_0.matcher( str );
        if ( m.find() ) {
            return m.group( 1 );
        }
        m = UNIPROT_KB_PATTERN_1.matcher( str );
        if ( m.find() ) {
            return m.group( 1 );
        }
        m = UNIPROT_KB_PATTERN_2.matcher( str );
        if ( m.find() ) {
            return m.group();
        }
        return null;
    }

    public final static boolean isProtein( final String query ) {
        final String r1 = parseRefSeqAccessor( query );
        if ( !ForesterUtil.isEmpty( r1 ) && ( r1.charAt( 1 ) == 'P' ) ) {
            return true;
        }
        final String r2 = extractUniProtKbProteinSeqIdentifier( query );
        if ( !ForesterUtil.isEmpty( r2 ) ) {
            return true;
        }
        return GENBANK_PROTEIN_AC_PATTERN.matcher( query ).lookingAt();
    }

    public final static Accession parse( final PhylogenyNode n ) {
        String v = extractUniProtKbProteinSeqIdentifier( n );
        if ( !ForesterUtil.isEmpty( v ) ) {
            return new Accession( v, Accession.UNIPROT );
        }
        v = extractGenbankAccessor( n );
        if ( !ForesterUtil.isEmpty( v ) ) {
            return new Accession( v, Accession.NCBI );
        }
        v = extractRefSeqAccessor( n );
        if ( !ForesterUtil.isEmpty( v ) ) {
            return new Accession( v, Accession.REFSEQ );
        }
        v = extractGInumber( n );
        if ( !ForesterUtil.isEmpty( v ) ) {
            return new Accession( v, Accession.GI );
        }
        return null;
    }

    public final static Accession obtainFromSeqAccession( final PhylogenyNode node ) {
        if ( node.getNodeData().isHasSequence() && ( node.getNodeData().getSequence().getAccession() != null )
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() )
                && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() ) ) {
            final String source = node.getNodeData().getSequence().getAccession().getSource().toLowerCase();
            final String value = node.getNodeData().getSequence().getAccession().getValue();
            if ( ( source.startsWith( "uniprot" ) || source.equals( "swissprot" ) || source.equals( "trembl" ) || source
                    .equals( "sp" ) ) ) {
                return new Accession( value, Accession.UNIPROT );
            }
            else if ( source.equals( "embl" ) || source.equals( "ebi" ) ) {
                return new Accession( value, Accession.EMBL );
            }
            else if ( source.equals( "ncbi" ) || source.equals( "genbank" ) ) {
                return new Accession( value, Accession.NCBI );
            }
            else if ( source.equals( "refseq" ) ) {
                return new Accession( value, Accession.REFSEQ );
            }
            else if ( source.equals( "gi" ) ) {
                return new Accession( value, Accession.GI );
            }
        }
        return null;
    }

    /**
     * Returns null if no match.
     * 
     */
    public final static Accession parse( final String s ) {
        if ( !ForesterUtil.isEmpty( s ) ) {
            String v = extractUniProtKbProteinSeqIdentifier( s );
            if ( !ForesterUtil.isEmpty( v ) ) {
                return new Accession( v, Accession.UNIPROT );
            }
            v = parseGenbankAccessor( s );
            if ( !ForesterUtil.isEmpty( v ) ) {
                return new Accession( v, Accession.NCBI );
            }
            v = parseRefSeqAccessor( s );
            if ( !ForesterUtil.isEmpty( v ) ) {
                return new Accession( v, Accession.REFSEQ );
            }
            v = parseGInumber( s );
            if ( !ForesterUtil.isEmpty( v ) ) {
                return new Accession( v, Accession.GI );
            }
        }
        return null;
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

    public static String parseGenbankProteinAccessor( final String query ) {
        final Matcher m = GENBANK_PROTEIN_AC_PATTERN.matcher( query );
        if ( m.lookingAt() ) {
            return m.group( 1 );
        }
        else {
            return null;
        }
    }

    public static String parseGInumber( final String query ) {
        final Matcher m = GI_PATTERN.matcher( query );
        if ( m.find() ) {
            return m.group( 1 );
        }
        return null;
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
}
