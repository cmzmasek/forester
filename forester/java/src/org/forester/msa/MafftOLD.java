
package org.forester.msa;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;

import org.forester.io.parsers.FastaParser;
import org.forester.util.ExternalProgram;
import org.forester.util.ForesterUtil;

public final class MafftOLD implements MsaInferrer {

    private String       _error;
    private int          _exit_code;
    private final String _path_to_prg;

    public static MsaInferrer createInstance( final String path_to_prg ) {
        return new MafftOLD( path_to_prg );
    }

    private MafftOLD( final String path_to_prg ) {
        _path_to_prg = new String( path_to_prg );
        init();
    }

    @Override
    public Object clone() {
        throw new NoSuchMethodError();
    }

    public String getErrorDescription() {
        return _error;
    }

    public int getExitCode() {
        return _exit_code;
    }

    public Msa infer( final File path_to_input_seqs, final List<String> opts ) throws IOException, InterruptedException {
        init();
        final String[] my_opts = new String[ opts.size() + 1 ];
        for( int i = 0; i < opts.size(); i++ ) {
            my_opts[ i ] = opts.get( i );
        }
        my_opts[ opts.size() ] = path_to_input_seqs.getAbsolutePath();
        final ExternalProgram mafft_prg = new ExternalProgram( _path_to_prg );
        mafft_prg.launch( my_opts );
        // _exit_code = mafft_prg.waitFor();
        // if ( _exit_code != 0 ) {
        //    throw new IOException( "MAFFT failed, exit code: " + _exit_code );
        // }
        final BufferedReader r = new BufferedReader( new InputStreamReader( mafft_prg.getErrorStream() ) );
        final StringBuffer error_sb = new StringBuffer();
        String line = null;
        while ( ( line = r.readLine() ) != null ) {
            error_sb.append( line );
            error_sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        r.close();
        if ( error_sb.length() > 0 ) {
            _error = error_sb.toString();
            throw new IOException( "MAFFT failed" );
        }
        final InputStream is = mafft_prg.getInputStream();
        final Msa msa = FastaParser.parseMsa( is );
        is.close();
        return msa;
    }

    private void init() {
        _error = null;
        _exit_code = -100;
    }
}
