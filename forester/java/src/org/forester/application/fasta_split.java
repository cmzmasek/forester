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
import org.forester.sequence.MolecularSequence;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class fasta_split {

    final static private String PRG_NAME    = "fasta_split";
    final static private String PRG_VERSION = "1.00";
    final static private String PRG_DATE    = "150325";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( fasta_split.PRG_NAME, fasta_split.PRG_VERSION, fasta_split.PRG_DATE );
        System.out.println();
        if ( ( args.length != 3 ) ) {
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
        final Map<String, List<MolecularSequence>> output = new HashMap<String, List<MolecularSequence>>();
        int cc = 0;
        for( final MolecularSequence seq : seqs ) {
            System.out.println( ++cc );
            final Matcher m = pa.matcher( seq.getIdentifier() );
            if ( m.find() ) {
                final String key = m.group( 1 );
                if ( !output.containsKey( key ) ) {
                    output.put( key, new ArrayList<MolecularSequence>() );
                }
                output.get( key ).add( seq );
            }
            else {
                //ForesterUtil.fatalError( PRG_NAME, pattern_str + " not found in sequence \"" + seq.getIdentifier() + "\"" );
                System.out.println( "warning: " + pattern_str + " not found in sequence \"" + seq.getIdentifier() + "\"" );
                final String key = "unknown";
                if ( !output.containsKey( key ) ) {
                    output.put( key, new ArrayList<MolecularSequence>() );
                }
                output.get( key ).add( seq );
            }
        }
        int c = 0;
        for( final Map.Entry<String, List<MolecularSequence>> entry : output.entrySet() ) {
            final File of = new File( outdir.getAbsolutePath().toString() + "/" + entry.getKey() + ".fasta" );
            if ( of.exists() ) {
                ForesterUtil.fatalError( PRG_NAME, of + " already exists" );
            }
            System.out.println( ++c + ": writing " + of );
            try {
                SequenceWriter.writeSeqs( entry.getValue(), of, SEQ_FORMAT.FASTA, 60 );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            }
        }
    }

    private static void argumentsError() {
        System.out.println( PRG_NAME + " <pattern> <infile> <outdir>" );
        System.out.println();
        System.exit( -1 );
    }
}
