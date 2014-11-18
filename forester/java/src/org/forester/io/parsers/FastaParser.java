// $Id:
//
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

package org.forester.io.parsers;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.msa.BasicMsa;
import org.forester.msa.Msa;
import org.forester.msa.MsaFormatException;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;

public class FastaParser {

    private static final Pattern NAME_REGEX      = Pattern.compile( "^\\s*>\\s*(.+)" );
    private static final Pattern SEQ_REGEX       = Pattern.compile( "^\\s*(.+)" );
    private static final Pattern ANYTHING_REGEX  = Pattern.compile( "[\\d\\s]+" );
    //>gi|71834668|ref|NP_001025424.1| Bcl2 [Danio rerio]
    public static final Pattern  FASTA_DESC_LINE = Pattern
            .compile( ">?\\s*([^|]+)\\|([^|]+)\\S*\\s+(.+)\\s+\\[(.+)\\]" );

    public static void main( final String[] args ) {
        final String a = ">gi|71834668|ref|NP_001025424.1| Bcl2 [Danio rerio]";
        final Matcher name_m = FASTA_DESC_LINE.matcher( a );
        if ( name_m.lookingAt() ) {
            System.out.println();
            System.out.println( name_m.group( 1 ) );
            System.out.println( name_m.group( 2 ) );
            System.out.println( name_m.group( 3 ) );
            System.out.println( name_m.group( 4 ) );
        }
        else {
            System.out.println( "Does not match." );
        }
    }

    static public boolean isLikelyFasta( final File f ) throws IOException {
        return isLikelyFasta( new FileInputStream( f ) );
    }

    static public boolean isLikelyFasta( final InputStream is ) throws IOException {
        final BufferedReader reader = new BufferedReader( new InputStreamReader( is, "UTF-8" ) );
        String line = null;
        while ( ( line = reader.readLine() ) != null ) {
            final boolean is_name_line = NAME_REGEX.matcher( line ).lookingAt();
            if ( canIgnore( line, true, false ) ) {
                continue;
            }
            else if ( is_name_line ) {
                reader.close();
                return true;
            }
            else if ( SEQ_REGEX.matcher( line ).lookingAt() ) {
                reader.close();
                return false;
            }
        }
        reader.close();
        return false;
    }

    static public Msa parseMsa( final File f ) throws IOException {
        return parseMsa( new FileInputStream( f ) );
    }

    static public Msa parseMsa( final InputStream is ) throws IOException {
        return BasicMsa.createInstance( parse( is ) );
    }

    static public Msa parseMsa( final String s ) throws IOException {
        return parseMsa( s.getBytes() );
    }

    static public Msa parseMsa( final byte[] bytes ) throws IOException {
        return parseMsa( new ByteArrayInputStream( bytes ) );
    }

    static public List<MolecularSequence> parse( final File f ) throws IOException {
        return parse( new FileInputStream( f ) );
    }

    static public List<MolecularSequence> parse( final InputStream is ) throws IOException {
        final BufferedReader reader = new BufferedReader( new InputStreamReader( is, "UTF-8" ) );
        String line = null;
        int line_counter = 0;
        boolean saw_first_seq = false;
        StringBuilder current_seq = null;
        StringBuilder name = null;
        final List<StringBuilder[]> temp_msa = new ArrayList<StringBuilder[]>();
        while ( ( line = reader.readLine() ) != null ) {
            ++line_counter;
            final Matcher name_m = NAME_REGEX.matcher( line );
            final boolean is_name_line = name_m.lookingAt();
            if ( canIgnore( line, saw_first_seq, is_name_line ) ) {
                continue;
            }
            final Matcher seq_m = SEQ_REGEX.matcher( line );
            if ( is_name_line ) {
                saw_first_seq = true;
                addSeq( name, current_seq, temp_msa );
                name = new StringBuilder( name_m.group( 1 ).trim() );
                current_seq = new StringBuilder();
            }
            else if ( seq_m.lookingAt() ) {
                if ( name.length() < 1 ) {
                    reader.close();
                    throw new MsaFormatException( "illegally formatted fasta msa (line: " + line_counter + "):\n\""
                            + trim( line ) + "\"" );
                }
                current_seq.append( seq_m.group( 1 ).replaceAll( "\\s+", "" ) );
            }
            else {
                reader.close();
                throw new MsaFormatException( "illegally formatted fasta msa (line: " + line_counter + "):\n\""
                        + trim( line ) + "\"" );
            }
        }
        addSeq( name, current_seq, temp_msa );
        reader.close();
        final List<MolecularSequence> seqs = new ArrayList<MolecularSequence>();
        for( int i = 0; i < temp_msa.size(); ++i ) {
            seqs.add( BasicSequence.createAaSequence( temp_msa.get( i )[ 0 ].toString(),
                                                      temp_msa.get( i )[ 1 ].toString() ) );
        }
        return seqs;
    }

    static private boolean canIgnore( final String line, final boolean saw_first_seq, final boolean is_name_line ) {
        if ( ( line.length() < 1 ) || ANYTHING_REGEX.matcher( line ).matches() ) {
            return true;
        }
        if ( !saw_first_seq && !is_name_line ) {
            return true;
        }
        return false;
    }

    private static void addSeq( final StringBuilder name, final StringBuilder seq, final List<StringBuilder[]> temp_msa ) {
        if ( ( name != null ) && ( seq != null ) && ( name.length() > 0 ) && ( seq.length() > 0 ) ) {
            final StringBuilder[] ary = new StringBuilder[ 2 ];
            ary[ 0 ] = name;
            ary[ 1 ] = seq;
            temp_msa.add( ary );
        }
    }

    private static String trim( final String line ) {
        if ( line.length() > 100 ) {
            return line.substring( 0, 100 ) + " ...";
        }
        return line;
    }
}
