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
//
// "java -Xmx1024m -cp path\to\forester.jar org.forester.application.fasta_split
//
//

package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class fasta_split {

    final static private String PRG_NAME    = "fasta_split";
    final static private String PRG_VERSION = "1.1.0";
    final static private String PRG_DATE    = "2024-05-06";

    public static void main( final String[] args) {
        ForesterUtil.printProgramInformation( fasta_split.PRG_NAME, fasta_split.PRG_VERSION, fasta_split.PRG_DATE );
        System.out.println();
        if ( ( args.length != 3 && args.length != 4 ) ) {
            fasta_split.argumentsError();
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        final String pattern_str = cla.getName( 0 );
        final File infile = cla.getFile( 1 );
        final File outdir = cla.getFile( 2 );

        final String out_pattern;
        if ( args.length == 4 ) {
            out_pattern = cla.getName( 3 );
        }
        else {
            out_pattern = null;
        }

        Pattern pa = null;
        try {
            pa = Pattern.compile( pattern_str );
        }
        catch ( final Exception ex ) {
            ForesterUtil.fatalError( PRG_NAME, ex.getMessage() );
        }
        final String error = ForesterUtil.isReadableFile( infile );
        if ( !ForesterUtil.isEmpty( error ) ) {
            ForesterUtil.fatalError( PRG_NAME, error );
        }
        if ( !outdir.exists() ) {
            new File( outdir.toString() ).mkdir();
        }
        if ( !outdir.isDirectory() ) {
            ForesterUtil.fatalError( PRG_NAME, outdir + " is not a directory" );
        }
        List<MolecularSequence> seqs = null;
        try {
            seqs = FastaParser.parse( new FileInputStream( infile ) );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( ( seqs == null ) || seqs.isEmpty() ) {
            ForesterUtil.fatalError( PRG_NAME, infile + " appears empty" );
        }
        System.out.println( "Read " + seqs.size() + " sequences" );
        final Map<String, List<MolecularSequence>> output = new HashMap<String, List<MolecularSequence>>();
        for( final MolecularSequence seq : seqs ) {
            final String id = seq.getIdentifier().trim().replaceAll( "\\s*\\|\\s*", "|").replaceAll( "\\s+", "_");
            ( ( BasicSequence ) seq ).setIdentifier(id);
            final Matcher m = pa.matcher(seq.getIdentifier() );
            if ( m.find() ) {
                final String key = m.group( 1 ).toUpperCase();
                System.out.println( "Extracted: " + key );
                if ( !output.containsKey( key ) ) {
                    output.put( key, new ArrayList<MolecularSequence>() );
                }
                output.get( key ).add( seq );
            }
            else {
                System.out.println( "Warning: " + pattern_str + " not found in sequence \"" + seq.getIdentifier()
                        + "\"" );
                final String key = "UNKNOWN";
                if ( !output.containsKey( key ) ) {
                    output.put( key, new ArrayList<MolecularSequence>() );
                }
                output.get( key ).add( seq );
            }
        }
        int c = 0;
        int seqs_written = 0;
        for( final Map.Entry<String, List<MolecularSequence>> entry : output.entrySet() ) {
            if ( entry.getKey() == null) {
                ForesterUtil.fatalError( PRG_NAME, "regular expression appears faulty");
            }
            String s = entry.getKey().trim();

            s = s.replaceAll( "[\\./\\*\\s]+", "_" );
            s = s.replaceAll( "\\(", "~" );
            s = s.replaceAll( "\\)", "~" );
            final File of;
            if (out_pattern != null) {
                of = new File( outdir.getAbsolutePath() + "/" + out_pattern + s + ".fasta" );
            }
            else {
                of = new File( outdir.getAbsolutePath() + "/" + s + ".fasta" );
            }

            if ( of.exists() ) {
                ForesterUtil.fatalError( PRG_NAME, of + " already exists" );
            }
            System.out.println( ++c + ": writing " + of + " [" + entry.getValue().size() + " seqs]" );
            
            try {
                SequenceWriter.writeSeqs( entry.getValue(), of, SEQ_FORMAT.FASTA, 60 );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            }
            seqs_written += entry.getValue().size();
        }
        System.out.println( "Wrote " + seqs_written + " sequences" );
    }

    private static void argumentsError() {
        System.out.println( PRG_NAME + " <pattern> <infile> <outdir> [out base name]" );
        System.out.println();
        System.out.println( "Examples: " );
        System.out.println( "  " + PRG_NAME + " \"v-germ=(\\S+)\" tt.fasta outdir" );
        System.out.println( "  " + PRG_NAME + " \"(\\S+?)\\|\" seqs.fasta outdir" );
        System.out.println( "  " + PRG_NAME + " \"OS=(.+?)[A-Z]{2}=\" seqs.fasta outdir" );
        System.out.println( "  " + PRG_NAME + " \"^.*?\\\\|.*?\\\\|.*?\\\\|(.+?)\\\\|\" seqs.fasta outdir" );
        System.out.println( "  " + PRG_NAME + " \"segment.(\\d+)\" sequence.fasta outdir A_thick_billed_murre_Greenland_9045-2K_2014_Segment_" );


        System.out.println();
        System.exit( -1 );
    }
}
