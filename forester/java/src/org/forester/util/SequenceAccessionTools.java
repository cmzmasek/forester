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
    public final static Pattern  GENBANK_NUC_PATTERN_1 = Pattern
                                                               .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]\\d{5}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)" );
    public final static Pattern  GENBANK_NUC_PATTERN_2 = Pattern
                                                               .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]{2}\\d{6}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)" );
    public final static Pattern  GENBANK_PROT_PATTERN  = Pattern
                                                               .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]{3}\\d{5}(?:\\.\\d+)?)(?:[^a-zA-Z0-9]|\\Z)" );
    public final static Pattern  GI_PATTERN            = Pattern.compile( "(?:\\b|_)(?:GI|gi)[|_=:](\\d+)(?:\\b|_)" );
    public final static Pattern  UNIPROT_KB_PATTERN_0  = Pattern.compile( "\\b([A-Z][0-9][A-Z0-9]{3}[0-9])\\b" );
    public final static Pattern  UNIPROT_KB_PATTERN_1  = Pattern
                                                               .compile( "(?:\\b|_)(?:sp|tr)[\\.|\\-_=/\\\\]([A-Z][0-9][A-Z0-9]{3}[0-9])(?:\\b|_)" );
    public final static Pattern  UNIPROT_KB_PATTERN_2  = Pattern
                                                               .compile( "(?:\\b|_)(?:[A-Z0-9]{2,5}|(?:[A-Z][0-9][A-Z0-9]{3}[0-9]))_(([A-Z9][A-Z]{2}[A-Z0-9]{2})|RAT|PIG|PEA)(?:\\b|_)" );
    // RefSeq accession numbers can be distinguished from GenBank accessions 
    // by their distinct prefix format of 2 characters followed by an
    // underscore character ('_'). For example, a RefSeq protein accession is NP_015325. 
    private final static Pattern REFSEQ_PATTERN        = Pattern
                                                               .compile( "(?:\\A|.*[^a-zA-Z0-9])([A-Z]{2}_\\d{6,})(?:[^a-zA-Z0-9]|\\Z)" );

    private SequenceAccessionTools() {
        // Hiding the constructor.
    }

    public final static boolean isProteinDbQuery( final String query ) {
        final String r1 = parseRefSeqAccessorFromString( query );
        if ( !ForesterUtil.isEmpty( r1 ) && ( r1.charAt( 1 ) == 'P' ) ) {
            return true;
        }
        final String r2 = parseUniProtAccessorFromString( query );
        if ( !ForesterUtil.isEmpty( r2 ) ) {
            return true;
        }
        return GENBANK_PROT_PATTERN.matcher( query ).lookingAt();
    }

    public final static Accession obtainAccessorFromDataFields( final PhylogenyNode n ) {
        String a = obtainUniProtAccessorFromDataFields( n );
        if ( !ForesterUtil.isEmpty( a ) ) {
            return new Accession( a, Accession.UNIPROT );
        }
        a = obtainGenbankAccessorFromDataFields( n );
        if ( !ForesterUtil.isEmpty( a ) ) {
            return new Accession( a, Accession.NCBI );
        }
        a = obtainRefSeqAccessorFromDataFields( n );
        if ( !ForesterUtil.isEmpty( a ) ) {
            return new Accession( a, Accession.REFSEQ );
        }
        a = obtainGiNumberFromDataFields( n );
        if ( !ForesterUtil.isEmpty( a ) ) {
            return new Accession( a, Accession.GI );
        }
        return null;
    }

    public final static Accession obtainFromSeqAccession( final PhylogenyNode n ) {
        if ( n.getNodeData().isHasSequence() && ( n.getNodeData().getSequence().getAccession() != null )
                && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getAccession().getSource() )
                && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getAccession().getValue() ) ) {
            final String source = n.getNodeData().getSequence().getAccession().getSource().toLowerCase();
            final String value = n.getNodeData().getSequence().getAccession().getValue();
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

    public final static String obtainGenbankAccessorFromDataFields( final PhylogenyNode n ) {
        String a = null;
        if ( n.getNodeData().isHasSequence() ) {
            final Sequence seq = n.getNodeData().getSequence();
            if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
                a = parseGenbankAccessorFromString( seq.getSymbol() );
            }
            if ( !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
                a = parseGenbankAccessorFromString( seq.getGeneName() );
            }
            if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( seq.getName() ) ) {
                a = parseGenbankAccessorFromString( seq.getName() );
            }
            if ( ForesterUtil.isEmpty( a ) && ( n.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
                a = parseGenbankAccessorFromString( seq.getAccession().getValue() );
            }
        }
        if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( n.getName() ) ) {
            a = parseGenbankAccessorFromString( n.getName() );
        }
        return a;
    }

    public final static String obtainGiNumberFromDataFields( final PhylogenyNode n ) {
        String a = null;
        if ( n.getNodeData().isHasSequence() ) {
            final Sequence seq = n.getNodeData().getSequence();
            if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( seq.getName() ) ) {
                a = parseGInumberFromString( seq.getName() );
            }
            if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
                a = parseGInumberFromString( seq.getGeneName() );
            }
            if ( ForesterUtil.isEmpty( a ) && ( n.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
                a = parseGInumberFromString( seq.getAccession().getValue() );
            }
        }
        if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( n.getName() ) ) {
            a = parseGInumberFromString( n.getName() );
        }
        return a;
    }

    public final static String obtainRefSeqAccessorFromDataFields( final PhylogenyNode n ) {
        String a = null;
        if ( n.getNodeData().isHasSequence() ) {
            final Sequence seq = n.getNodeData().getSequence();
            if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
                a = parseRefSeqAccessorFromString( seq.getSymbol() );
            }
            if ( !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
                a = parseRefSeqAccessorFromString( seq.getGeneName() );
            }
            if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( seq.getName() ) ) {
                a = parseRefSeqAccessorFromString( seq.getName() );
            }
            if ( ForesterUtil.isEmpty( a ) && ( n.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
                a = parseRefSeqAccessorFromString( seq.getAccession().getValue() );
            }
        }
        if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( n.getName() ) ) {
            a = parseRefSeqAccessorFromString( n.getName() );
        }
        return a;
    }

    public final static String obtainUniProtAccessorFromDataFields( final PhylogenyNode n ) {
        String a = null;
        if ( n.getNodeData().isHasSequence() ) {
            final Sequence seq = n.getNodeData().getSequence();
            if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
                a = SequenceAccessionTools.parseUniProtAccessorFromString( seq.getSymbol() );
            }
            if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( seq.getName() ) ) {
                a = SequenceAccessionTools.parseUniProtAccessorFromString( seq.getName() );
            }
            if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
                a = SequenceAccessionTools.parseUniProtAccessorFromString( seq.getGeneName() );
            }
            if ( ForesterUtil.isEmpty( a ) && ( n.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
                a = SequenceAccessionTools.parseUniProtAccessorFromString( seq.getAccession().getValue() );
            }
        }
        if ( ForesterUtil.isEmpty( a ) && !ForesterUtil.isEmpty( n.getName() ) ) {
            a = SequenceAccessionTools.parseUniProtAccessorFromString( n.getName() );
        }
        return a;
    }

    public final static Accession parseAccessorFromString( final String s ) {
        if ( !ForesterUtil.isEmpty( s ) ) {
            String v = parseUniProtAccessorFromString( s );
            if ( !ForesterUtil.isEmpty( v ) ) {
                return new Accession( v, Accession.UNIPROT );
            }
            v = parseGenbankAccessorFromString( s );
            if ( !ForesterUtil.isEmpty( v ) ) {
                return new Accession( v, Accession.NCBI );
            }
            v = parseRefSeqAccessorFromString( s );
            if ( !ForesterUtil.isEmpty( v ) ) {
                return new Accession( v, Accession.REFSEQ );
            }
            v = parseGInumberFromString( s );
            if ( !ForesterUtil.isEmpty( v ) ) {
                return new Accession( v, Accession.GI );
            }
        }
        return null;
    }

    public final static String parseGenbankAccessorFromString( final String s ) {
        Matcher m = GENBANK_NUC_PATTERN_1.matcher( s );
        if ( m.lookingAt() ) {
            return m.group( 1 );
        }
        else {
            m = GENBANK_NUC_PATTERN_2.matcher( s );
            if ( m.lookingAt() ) {
                return m.group( 1 );
            }
            else {
                m = GENBANK_PROT_PATTERN.matcher( s );
                if ( m.lookingAt() ) {
                    return m.group( 1 );
                }
                else {
                    return null;
                }
            }
        }
    }

    public final static String parseGenbankProteinAccessorFromString( final String s ) {
        final Matcher m = GENBANK_PROT_PATTERN.matcher( s );
        if ( m.lookingAt() ) {
            return m.group( 1 );
        }
        else {
            return null;
        }
    }

    public final static String parseGInumberFromString( final String s ) {
        final Matcher m = GI_PATTERN.matcher( s );
        if ( m.find() ) {
            return m.group( 1 );
        }
        return null;
    }

    public final static String parseRefSeqAccessorFromString( final String s ) {
        final Matcher m = REFSEQ_PATTERN.matcher( s );
        if ( m.lookingAt() ) {
            return m.group( 1 );
        }
        return null;
    }

    public final static String parseUniProtAccessorFromString( final String s ) {
        Matcher m = UNIPROT_KB_PATTERN_0.matcher( s );
        if ( m.find() ) {
            return m.group( 1 );
        }
        m = UNIPROT_KB_PATTERN_1.matcher( s );
        if ( m.find() ) {
            return m.group( 1 );
        }
        m = UNIPROT_KB_PATTERN_2.matcher( s );
        if ( m.find() ) {
            return m.group();
        }
        return null;
    }
}
