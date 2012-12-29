// $Id:
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

package org.forester.msa;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.Sequence;
import org.forester.util.SystemCommandExecutor;

public final class Mafft extends MsaInferrer {

    private final static String DEFAULT_PARAMETERS = "--maxiterate 1000 --localpair";
    private String              _error;
    private int                 _exit_code;
    private final String        _path_to_prg;

    public static MsaInferrer createInstance( final String path_to_prg ) throws IOException {
        return new Mafft( path_to_prg );
    }

    private Mafft( final String path_to_prg ) throws IOException {
        if ( !isInstalled( path_to_prg ) ) {
            throw new IOException( "cannot execute MAFFT with \"" + path_to_prg + "\"" );
        }
        _path_to_prg = new String( path_to_prg );
        init();
    }

    public static String getDefaultParameters() {
        return DEFAULT_PARAMETERS;
    }

    @Override
    public String getErrorDescription() {
        return _error;
    }

    @Override
    public int getExitCode() {
        return _exit_code;
    }

    @Override
    public Msa infer( final List<Sequence> seqs, final List<String> opts ) throws IOException, InterruptedException {
        final File file = File.createTempFile( "__mafft_input_", ".fasta" );
        file.deleteOnExit();
        final BufferedWriter writer = new BufferedWriter( new FileWriter( file ) );
        SequenceWriter.writeSeqs( seqs, writer, SEQ_FORMAT.FASTA, 100 );
        writer.close();
        final Msa msa = infer( file, opts );
        file.delete();
        return msa;
    }

    @Override
    public Msa infer( final File path_to_input_seqs, final List<String> opts ) throws IOException, InterruptedException {
        init();
        final List<String> my_opts = new ArrayList<String>();
        my_opts.add( _path_to_prg );
        for( int i = 0; i < opts.size(); i++ ) {
            my_opts.add( opts.get( i ) );
        }
        my_opts.add( path_to_input_seqs.getAbsolutePath() );
        final SystemCommandExecutor command_executor = new SystemCommandExecutor( my_opts );
        final int _exit_code = command_executor.executeCommand();
        final StringBuilder stderr = command_executor.getStandardErrorFromCommand();
        _error = stderr.toString();
        if ( _exit_code != 0 ) {
            throw new IOException( "MAFFT program failed, exit code: " + _exit_code + "\nCommand:\n" + my_opts
                    + "\nError:\n" + stderr );
        }
        final StringBuilder stdout = command_executor.getStandardOutputFromCommand();
        if ( ( stdout == null ) || ( stdout.length() < 2 ) ) {
            throw new IOException( "MAFFT program did not produce any output\nCommand:\n" + my_opts + "\nError:\n"
                    + stderr );
        }
        final Msa msa = FastaParser.parseMsa( stdout.toString() );
        return msa;
    }

    private void init() {
        _error = null;
        _exit_code = -100;
    }
}
