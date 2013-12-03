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

package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.Sequence;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class check_fasta {

    final static private String PRG_NAME    = "check_fasta";
    final static private String PRG_VERSION = "1.00";
    final static private String PRG_DATE    = "131202";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( check_fasta.PRG_NAME, check_fasta.PRG_VERSION, check_fasta.PRG_DATE );
        System.out.println();
        if ( ( args.length != 2 ) ) {
            check_fasta.argumentsError();
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        final File indir = cla.getFile( 0 );
        final File outdir = cla.getFile( 1 );
        if ( !indir.isDirectory() ) {
            ForesterUtil.fatalError( PRG_NAME, indir + " is not a directory" );
        }
        if ( !outdir.isDirectory() ) {
            ForesterUtil.fatalError( PRG_NAME, outdir + " is not a directory" );
        }
        final File[] list_of_files = indir.listFiles();
        final List<File> infiles = new ArrayList<File>();
        for( final File file : list_of_files ) {
            if ( file.isFile()
                    && file.canRead()
                    && ( file.toString().toLowerCase().endsWith( ".fasta" ) || file.toString().toLowerCase()
                            .endsWith( ".fas" ) ) ) {
                infiles.add( file );
            }
        }
        Collections.sort( infiles );
        int c = 0;
        for( final File infile : infiles ) {
            System.out.println( ++c + "/" + infiles.size() + ": " + infile );
            execute( outdir, infile );
        }
    }

    private static void execute( final File outdir, final File infile ) {
        final File outfile = new File( outdir.getAbsolutePath().toString() + "/" + infile.getName() );
        if ( outfile.exists() ) {
            System.out.println( outfile + " already exists" );
        }
        else {
            try {
                final List<Sequence> seqs = FastaParser.parse( new FileInputStream( infile ) );
                final Map<String, Short> names = new HashMap<String, Short>();
                for( final Sequence seq : seqs ) {
                    procSeq( infile.toString(), names, seq );
                }
                SequenceWriter.writeSeqs( seqs, outfile, SEQ_FORMAT.FASTA, 60 );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            }
        }
    }

    private static void procSeq( final String infile, final Map<String, Short> names, final Sequence seq ) {
        final String name = seq.getIdentifier();
        if ( !names.containsKey( name ) ) {
            names.put( name, ( short ) 1 );
        }
        else {
            final short i = names.get( name );
            ( ( BasicSequence ) seq ).setIdentifier( name + "_" + i );
            names.put( name, ( short ) ( i + 1 ) );
            System.out.println( "  " + infile + i + ": " + seq.getIdentifier() );
        }
    }

    private static void argumentsError() {
        System.out.println( PRG_NAME + " <indir> <outdir>" );
        System.out.println();
        System.exit( -1 );
    }
}
