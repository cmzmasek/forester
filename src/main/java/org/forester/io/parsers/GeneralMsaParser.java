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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.msa.BasicMsa;
import org.forester.msa.Msa;
import org.forester.msa.MsaFormatException;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;

public final class GeneralMsaParser {

    private static final Pattern NAME_SEQ_PATTERN          = Pattern.compile( "(\\S+)\\s+(\\S+)\\s*" );
    private static final Pattern INDENTED_SEQ_PATTERN      = Pattern.compile( "\\s+(\\S+)\\s*" );
    private static final Pattern NON_INDENTED_SEQ_PATTERN  = Pattern.compile( "(\\S+).*" );
    private static final Pattern PROBCONS_REGEX            = Pattern.compile( "^CLUSTAL" );
    private static final Pattern MUSCLE_REGEX              = Pattern.compile( "^MUSCLE\\s\\(" );
    private static final Pattern CLUSTAL_REGEX             = Pattern.compile( "^PROBCONS\\s" );
    private static final Pattern ANYTHING_REGEX            = Pattern.compile( "[\\d\\s]+" );
    private static final Pattern SELEX_SPECIAL_LINES_REGEX = Pattern.compile( "\\s+[*\\.:\\s]+" );
    private static final Pattern SPECIAL_LINES_REGEX       = Pattern.compile( "^\\s*(#|%|//|!!)" );
    private static final Pattern ERROR_REGEX               = Pattern.compile( "\\S+\\s+\\S+\\s+\\S+" );

    static private boolean canIgnore( final String line ) {
        if ( ( line.length() < 1 ) || ANYTHING_REGEX.matcher( line ).matches() ) {
            return true;
        }
        return ( SELEX_SPECIAL_LINES_REGEX.matcher( line ).matches() || SPECIAL_LINES_REGEX.matcher( line ).lookingAt() );
    }

    static private boolean isProgramNameLine( final String line ) {
        return ( PROBCONS_REGEX.matcher( line ).lookingAt() || CLUSTAL_REGEX.matcher( line ).lookingAt() || MUSCLE_REGEX
                .matcher( line ).lookingAt() );
    }

    static public Msa parse( final InputStream is ) throws IOException {
        int block = -1;
        int current_seq_index_per_block = -1;
        String current_name = null;
        boolean saw_ignorable = true;
        boolean is_first = true;
        final Map<String, StringBuilder> temp_msa = new HashMap<String, StringBuilder>();
        final List<String> names_in_order = new ArrayList<String>();
        final BufferedReader reader = new BufferedReader( new InputStreamReader( is, "UTF-8" ) );
        String line = null;
        int line_counter = 0;
        while ( ( line = reader.readLine() ) != null ) {
            ++line_counter;
            if ( canIgnore( line ) ) {
                saw_ignorable = true;
            }
            else if ( !( is_first && isProgramNameLine( line ) ) ) {
                if ( ERROR_REGEX.matcher( line ).lookingAt() ) {
                    throw new MsaFormatException( "unrecognized msa format (line: " + line_counter + "):\n\""
                            + trim( line ) + "\"" );
                }
                if ( canIgnore( line ) ) {
                    saw_ignorable = true;
                }
                final Matcher name_seq_m = NAME_SEQ_PATTERN.matcher( line );
                Matcher ind_seq_m = null;
                Matcher non_ind_seq_m = null;
                boolean ind_seq_m_matches = false;
                boolean non_ind_seq_m_matches = false;
                final boolean name_seq_m_matches = name_seq_m.matches();
                if ( !name_seq_m_matches ) {
                    ind_seq_m = INDENTED_SEQ_PATTERN.matcher( line );
                    ind_seq_m_matches = ind_seq_m.matches();
                    if ( !ind_seq_m_matches ) {
                        non_ind_seq_m = NON_INDENTED_SEQ_PATTERN.matcher( line );
                        non_ind_seq_m_matches = non_ind_seq_m.lookingAt();
                    }
                }
                if ( name_seq_m_matches || ind_seq_m_matches || non_ind_seq_m_matches ) {
                    if ( saw_ignorable ) {
                        ++block;
                        current_seq_index_per_block = -1;
                        saw_ignorable = false;
                    }
                    ++current_seq_index_per_block;
                    if ( name_seq_m_matches ) {
                        final String name = name_seq_m.group( 1 );
                        final String seq = name_seq_m.group( 2 );
                        if ( temp_msa.containsKey( name ) ) {
                            temp_msa.get( name ).append( seq );
                        }
                        else {
                            temp_msa.put( name, new StringBuilder( seq ) );
                            names_in_order.add( name );
                        }
                        current_name = name;
                    }
                    else if ( ind_seq_m_matches ) {
                        if ( temp_msa.containsKey( current_name ) ) {
                            temp_msa.get( current_name ).append( ind_seq_m.group( 1 ) );
                        }
                        else {
                            throw new MsaFormatException( "illegal msa format (line: " + line_counter + "):\n\""
                                    + trim( line ) + "\"" );
                        }
                    }
                    else if ( non_ind_seq_m_matches ) {
                        if ( block == 0 ) {
                            throw new MsaFormatException( "illegal msa format: first block cannot contain un-named sequence (line: "
                                    + line_counter + "):\n\"" + trim( line ) + "\"" );
                        }
                        else {
                            String name = "";
                            try {
                                name = names_in_order.get( current_seq_index_per_block );
                            }
                            catch ( final IndexOutOfBoundsException e ) {
                                throw new MsaFormatException( "illegalmsa format (line: " + line_counter + "):\n\""
                                        + trim( line ) + "\"" );
                            }
                            if ( temp_msa.containsKey( name ) ) {
                                temp_msa.get( name ).append( non_ind_seq_m.group( 1 ) );
                            }
                            else {
                                throw new MsaFormatException( "illegal msa format (line: " + line_counter + "):\n\""
                                        + trim( line ) + "\"" );
                            }
                        }
                        current_name = null;
                    }
                }
                else {
                    throw new MsaFormatException( "illegal msa format (line: " + line_counter + "):\n\"" + trim( line )
                                                  + "\"" );
                }
                if ( is_first ) {
                    is_first = false;
                }
            }
        } // while ( ( line = reader.readLine() ) != null )
        final List<MolecularSequence> seqs = new ArrayList<MolecularSequence>();
        for( int i = 0; i < names_in_order.size(); ++i ) {
            seqs.add( BasicSequence.createAaSequence( names_in_order.get( i ), temp_msa.get( names_in_order.get( i ) )
                                                      .toString() ) );
        }
        final Msa msa = BasicMsa.createInstance( seqs );
        return msa;
    }

    private static String trim( final String line ) {
        if ( line.length() > 100 ) {
            return line.substring( 0, 100 ) + " ...";
        }
        return line;
    }
}
