
package org.forester.application;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.GeneralMsaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.msa.BasicMsa;
import org.forester.msa.Msa;
import org.forester.msa.Msa.MSA_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public class tap {

    final static private String PRG_NAME                = "tap";
    final static private String PRG_DATE                = "170327";
    final static private String PRG_DESC                = "Replacement of labels in multiple sequence files";
    final static private String PRG_VERSION             = "1.00";
    final static private String WWW                     = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String E_MAIL                  = "phyloxml@gmail.com";
    final static private String EXTRACT_TAXONOMY_OPTION = "t";
    final static private String ANNOTATION_OPTION       = "a";
    final static private String HELP_OPTION_1           = "help";
    final static private String HELP_OPTION_2           = "h";

    public static void main( final String args[] ) {
        try {
            ForesterUtil.printProgramInformation( PRG_NAME,
                                                  PRG_DESC,
                                                  PRG_VERSION,
                                                  PRG_DATE,
                                                  E_MAIL,
                                                  WWW,
                                                  ForesterUtil.getForesterLibraryInformation() );
            CommandLineArguments cla = null;
            try {
                cla = new CommandLineArguments( args );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            }
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) ) {
                System.out.println();
                print_help();
                System.exit( 0 );
            }
            String input = null;
            String output = null;
            String list_file = null;
            String i = null;
            if ( args.length == 3 ) {
                input = cla.getName( 0 );
                output = cla.getName( 1 );
                list_file = cla.getName( 2 );
            }
            else if ( args.length == 1 ) {
                input = cla.getName( 0 );
                i = null;
                if ( input.toLowerCase().endsWith( ".fasta" ) ) {
                    i = input.substring( 0, input.length() - 7 );
                }
                else if ( input.toLowerCase().endsWith( ".fsa" ) ) {
                    i = input.substring( 0, input.length() - 5 );
                }
                else {
                    i = input;
                }
                output = i + ForesterConstants.ID_NORMALIZED_FASTA_FILE_SUFFIX;
                list_file = i + ForesterConstants.ID_MAP_FILE_SUFFIX;
            }
            else {
                print_help();
                System.exit( -1 );
            }
            final List<String> allowed_options = new ArrayList<>();
            allowed_options.add( ANNOTATION_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            final File outfile_file = new File( output );
            final File listfile = new File( list_file );
            final File input_file = new File( input );
            final String error1 = ForesterUtil.isWritableFile( outfile_file );
            if ( !ForesterUtil.isEmpty( error1 ) ) {
                ForesterUtil.fatalError( PRG_NAME, error1 );
            }
            final String error2 = ForesterUtil.isWritableFile( listfile );
            if ( !ForesterUtil.isEmpty( error2 ) ) {
                ForesterUtil.fatalError( PRG_NAME, error2 );
            }
            final String error3 = ForesterUtil.isReadableFile( input_file );
            if ( !ForesterUtil.isEmpty( error3 ) ) {
                ForesterUtil.fatalError( PRG_NAME, error3 );
            }
            final boolean fasta_like = ForesterUtil.isLooksLikeFasta( input_file );
            final Msa.MSA_FORMAT output_format = MSA_FORMAT.FASTA;
            System.out.println();
            System.out.println( "Input alignment       : " + input );
            System.out.println( "Output alignment      : " + output );
            System.out.println( "Name list             : " + list_file );
            if ( fasta_like ) {
                System.out.println( "Input format          : Fasta" );
            }
            else {
                System.out.println( "Input format          : Phylip like" );
            }
            if ( output_format == MSA_FORMAT.FASTA ) {
                System.out.println( "Output format         : Fasta" );
            }
            else if ( output_format == MSA_FORMAT.NEXUS ) {
                System.out.println( "Output format         : Nexus" );
            }
            else if ( output_format == MSA_FORMAT.PHYLIP ) {
                System.out.println( "Output format         : Phylip" );
            }
            System.out.println();
            
            final List<MolecularSequence> seqs;
            final FileInputStream is = new FileInputStream( input_file );
            if ( FastaParser.isLikelyFasta( input_file ) ) {
                seqs = FastaParser.parse( is );
            }
            else {
                seqs = GeneralMsaParser.parseSeqs( is );
            }
            if ( seqs == null ) {
                ForesterUtil.fatalError( PRG_NAME, "failed to read MSA" );
            }
            if ( seqs.size() < 1 ) {
                ForesterUtil.fatalError( PRG_NAME, "MSA seems to be devoid of sequences" );
            }
           // TODO print number of seqs
           // TODO print number min length
           // TODO print max length
           // TODO OR
          //  TODO print length is aligned
          //  TODO if no aligned no phylip or nexus outpt
            //
           
            final List<MolecularSequence> seqs2 = new ArrayList<>();
            int counter = 0;
            final BufferedWriter writer = ForesterUtil.createBufferedWriter( list_file );
            for( final MolecularSequence seq : seqs ) {
                final String new_name = modify_name( seq.getIdentifier(), counter++, writer );
                final MolecularSequence ns = BasicSequence.createSequence( new_name,
                                                                           seq.getMolecularSequenceAsString() );
                seqs2.add( ns );
            }
            writer.flush();
            writer.close();
            final BufferedWriter seq_writer = ForesterUtil.createBufferedWriter( outfile_file );
            if ( ( output_format == MSA_FORMAT.NEXUS ) || ( output_format == MSA_FORMAT.PHYLIP ) ) {
                final Msa m = BasicMsa.createInstance( seqs2 );
                m.write( seq_writer, output_format );
            }
            else if ( output_format == MSA_FORMAT.FASTA ) {
                SequenceWriter.writeSeqs( seqs2, seq_writer, SEQ_FORMAT.FASTA, 60 );
            }
            seq_writer.flush();
            seq_writer.close();
            //                    Util.print_message( PRG_NAME, "wrote: " + list_file )
            //                    Util.print_message( PRG_NAME, "wrote: " + output )
        }
        catch ( final IllegalArgumentException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( PRG_NAME, "Unexpected errror!" );
        }
    }

    final static String modify_name( final String desc, final int counter, final Writer writer ) throws IOException {
        desc.replaceAll( "\\s+", " " );
        final String new_desc = Integer.toHexString( counter );
        if ( new_desc.length() > 9 ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "shortened identifier [" + new_desc + "] is too long (" + new_desc.length()
                                             + " characters)" );
        }
        writer.write( new_desc + "\t" + desc + "\n" );
        return new_desc;
    }

    private final static void print_help() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " [options] <gene tree file> <query>" );
        System.out.println();
        System.out.println( " options:" );
        //System.out.println( "  -" + SEP_OPTION + "=<separator>: the separator to be used" );
        System.out.println();
        System.out.println( "Example:" );
        System.out.println();
        System.out.println( " " + PRG_NAME + " -s=. my_tree.xml A.1.1.1" );
        System.out.println();
    }
}
