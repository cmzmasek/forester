
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.msa.Msa;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.tools.MutationInference;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class m_label {

    private final static String PROPERTY_REF = "aptx:branch_event";
    private final static String PRG_NAME     = "m_label";
    private static final String MIN_SEQS_OPTION = "m";

    public static void main( final String args[] ) {
        
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( MIN_SEQS_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        if ( ( args.length < 2 ) || ( args.length > 3 ) ) {
            System.out.println( "\nWrong number of arguments.\n" );
            System.exit( -1 );
        }
        int min_seqs = -1;
        if ( cla.isOptionSet( MIN_SEQS_OPTION ) ) {
            try {
                min_seqs = cla.getOptionValueAsInt( MIN_SEQS_OPTION );
            }
            catch ( IOException e ) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        
        final File infile = new File( args[ 0 ] );
        final File alignment_file;
        final File outfile;
        if ( !infile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + infile + "] does not exist" );
        }
        if ( args.length > 2 ) {
            alignment_file = new File( args[ 1 ] );
            outfile = new File( args[ 2 ] );
            if ( !infile.exists() ) {
                ForesterUtil.fatalError( PRG_NAME, "[" + alignment_file + "] does not exist" );
            }
        }
        else {
            alignment_file = null;
            outfile = new File( args[ 1 ] );
        }
        
        
        
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        Msa msa = null;
        if ( alignment_file != null ) {
            try {
                msa = FastaParser.parseMsa( alignment_file );
            }
            catch ( final IOException e ) {
                e.printStackTrace();
            }
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( infile, true );
            p = factory.create( infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + infile + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        if ( msa != null ) {
            MutationInference.addMolSeqs( msa, p );
        }
        MutationInference.inferBranchEvents( p, PROPERTY_REF );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outfile );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outfile + "]: " + e.getMessage() );
        }
        System.out.println( "[" + PRG_NAME + "] wrote: [" + outfile + "]" );
        System.out.println( "[" + PRG_NAME + "] OK" );
        System.out.println();
    }
}
