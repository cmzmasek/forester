
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
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

public class rid {

    final static private String PRG_NAME               = "rid";
    final static private String PRG_DATE               = "170902";
    final static private String PRG_DESC               = "sequence file reformatting and identifier normalization";
    final static private String PRG_VERSION            = "1.00";
    final static private String WWW                    = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String E_MAIL                 = "phyloxml@gmail.com";
    final static private String OUTPUT_FORMAT_OPTION   = "o";
    final static private String ID_NORM_OPTION         = "s";
    final static private String HELP_OPTION_1          = "help";
    final static private String HELP_OPTION_2          = "h";
    private static final String OUTPUT_FORMAT_FASTA    = "f";
    private static final String OUTPUT_FORMAT_PHYLIP   = "p";
    private static final String OUTPUT_FORMAT_NEXUS    = "n";
    private static final String OUTPUT_FORMAT_FASTA_L  = "fasta";
    private static final String OUTPUT_FORMAT_PHYLIP_L = "phylip";
    private static final String OUTPUT_FORMAT_NEXUS_L  = "nexus";

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
            if ( ( cla.getNumberOfNames() == 0 ) || cla.isOptionSet( HELP_OPTION_1 )
                    || cla.isOptionSet( HELP_OPTION_2 ) ) {
                System.out.println();
                print_help();
                System.exit( 0 );
            }
            String input_seqs_file_str = null;
            String output_seqs_file_str = null;
            String output_map_file_str = null;
            String input_seqs_name_wo_suffix = null;
            if ( ( cla.getNumberOfNames() == 2 ) || ( cla.getNumberOfNames() == 3 ) ) {
                input_seqs_file_str = cla.getName( 0 );
                output_seqs_file_str = cla.getName( 1 );
                if ( cla.getNumberOfNames() == 3 ) {
                    output_map_file_str = cla.getName( 2 );
                }
            }
            else if ( cla.getNumberOfNames() == 1 ) {
                input_seqs_file_str = cla.getName( 0 );
                input_seqs_name_wo_suffix = null;
                if ( input_seqs_file_str.toLowerCase().endsWith( ".fasta" ) ) {
                    input_seqs_name_wo_suffix = input_seqs_file_str.substring( 0, input_seqs_file_str.length() - 6 );
                }
                else if ( input_seqs_file_str.toLowerCase().endsWith( ".fsa" ) ) {
                    input_seqs_name_wo_suffix = input_seqs_file_str.substring( 0, input_seqs_file_str.length() - 4 );
                }
                else if ( input_seqs_file_str.toLowerCase().endsWith( ".phy" ) ) {
                    input_seqs_name_wo_suffix = input_seqs_file_str.substring( 0, input_seqs_file_str.length() - 4 );
                }
                else if ( input_seqs_file_str.toLowerCase().endsWith( ".aln" ) ) {
                    input_seqs_name_wo_suffix = input_seqs_file_str.substring( 0, input_seqs_file_str.length() - 4 );
                }
                else if ( input_seqs_file_str.toLowerCase().endsWith( ".phylip" ) ) {
                    input_seqs_name_wo_suffix = input_seqs_file_str.substring( 0, input_seqs_file_str.length() - 7 );
                }
                else if ( input_seqs_file_str.toLowerCase().endsWith( ".nex" ) ) {
                    input_seqs_name_wo_suffix = input_seqs_file_str.substring( 0, input_seqs_file_str.length() - 4 );
                }
                else if ( input_seqs_file_str.toLowerCase().endsWith( ".nexus" ) ) {
                    input_seqs_name_wo_suffix = input_seqs_file_str.substring( 0, input_seqs_file_str.length() - 5 );
                }
                else {
                    input_seqs_name_wo_suffix = input_seqs_file_str;
                }
                output_map_file_str = input_seqs_name_wo_suffix + ForesterConstants.ID_MAP_FILE_SUFFIX;
            }
            else {
                print_help();
                System.exit( -1 );
            }
            final List<String> allowed_options = new ArrayList<>();
            allowed_options.add( OUTPUT_FORMAT_OPTION );
            allowed_options.add( ID_NORM_OPTION );
            final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
            if ( dissallowed_options.length() > 0 ) {
                ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
            }
            final File input_seqs_file = new File( input_seqs_file_str );
            final String error0 = ForesterUtil.isReadableFile( input_seqs_file );
            if ( !ForesterUtil.isEmpty( error0 ) ) {
                ForesterUtil.fatalError( PRG_NAME, error0 );
            }
            final boolean input_seqs_fasta_like = ForesterUtil.isLooksLikeFasta( input_seqs_file );
            Msa.MSA_FORMAT output_format = MSA_FORMAT.FASTA;
            if ( cla.isOptionSet( OUTPUT_FORMAT_OPTION ) ) {
                if ( cla.isOptionValueSet( OUTPUT_FORMAT_OPTION ) ) {
                    final String output_format_str = cla.getOptionValue( OUTPUT_FORMAT_OPTION );
                    if ( output_format_str.equals( OUTPUT_FORMAT_FASTA )
                            || output_format_str.equalsIgnoreCase( OUTPUT_FORMAT_FASTA_L ) ) {
                        output_format = MSA_FORMAT.FASTA;
                    }
                    else if ( output_format_str.equals( OUTPUT_FORMAT_PHYLIP )
                            || output_format_str.equalsIgnoreCase( OUTPUT_FORMAT_PHYLIP_L ) ) {
                        output_format = MSA_FORMAT.PHYLIP;
                    }
                    else if ( output_format_str.equals( OUTPUT_FORMAT_NEXUS )
                            || output_format_str.equalsIgnoreCase( OUTPUT_FORMAT_NEXUS_L ) ) {
                        output_format = MSA_FORMAT.NEXUS;
                    }
                    else {
                        ForesterUtil.fatalError( PRG_NAME, "unknown output format option: " + output_format_str );
                    }
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "no value for output format option"  );
                }
            }
            final boolean normalize_identifiers;
            if ( cla.isOptionSet( ID_NORM_OPTION ) || ( cla.getNumberOfNames() == 3 ) ) {
                normalize_identifiers = true;
            }
            else {
                normalize_identifiers = false;
            }
            if ( normalize_identifiers && ForesterUtil.isEmpty( output_map_file_str ) ) {
                ForesterUtil.fatalError( PRG_NAME, "need to indicate name for output map file" );
            }
            final File output_map_file;
            if ( normalize_identifiers ) {
                output_map_file = new File( output_map_file_str );
                final String error = ForesterUtil.isWritableFile( output_map_file );
                if ( !ForesterUtil.isEmpty( error ) ) {
                    ForesterUtil.fatalError( PRG_NAME, error );
                }
            }
            else {
                output_map_file = null;
            }
            if ( cla.getNumberOfNames() == 1 ) {
                if ( normalize_identifiers ) {
                    if ( output_format == MSA_FORMAT.FASTA ) {
                        output_seqs_file_str = input_seqs_name_wo_suffix
                                + ForesterConstants.ID_NORMALIZED_FASTA_FILE_SUFFIX;
                    }
                    else if ( output_format == MSA_FORMAT.NEXUS ) {
                        output_seqs_file_str = input_seqs_name_wo_suffix
                                + ForesterConstants.ID_NORMALIZED_NEXUS_FILE_SUFFIX;
                    }
                    else if ( output_format == MSA_FORMAT.PHYLIP ) {
                        output_seqs_file_str = input_seqs_name_wo_suffix
                                + ForesterConstants.ID_NORMALIZED_PHYLIP_FILE_SUFFIX;
                    }
                }
                else {
                    if ( output_format == MSA_FORMAT.FASTA ) {
                        output_seqs_file_str = input_seqs_name_wo_suffix + ForesterConstants.FASTA_FILE_SUFFIX;
                        if ( ForesterUtil.isWritableFile( output_seqs_file_str ) != null ) {
                            output_seqs_file_str = input_seqs_name_wo_suffix + "_"
                                    + ForesterConstants.FASTA_FILE_SUFFIX;
                        }
                    }
                    else if ( output_format == MSA_FORMAT.NEXUS ) {
                        output_seqs_file_str = input_seqs_name_wo_suffix + ForesterConstants.NEXUS_FILE_SUFFIX;
                        if ( ForesterUtil.isWritableFile( output_seqs_file_str ) != null ) {
                            output_seqs_file_str = input_seqs_name_wo_suffix + "_"
                                    + ForesterConstants.NEXUS_FILE_SUFFIX;
                        }
                    }
                    else if ( output_format == MSA_FORMAT.PHYLIP ) {
                        output_seqs_file_str = input_seqs_name_wo_suffix + ForesterConstants.PHYLIP_FILE_SUFFIX;
                        if ( ForesterUtil.isWritableFile( output_seqs_file_str ) != null ) {
                            output_seqs_file_str = input_seqs_name_wo_suffix + "_"
                                    + ForesterConstants.PHYLIP_FILE_SUFFIX;
                        }
                    }
                }
            }
            final File outfile_seqs_file = new File( output_seqs_file_str );
            final String error1 = ForesterUtil.isWritableFile( outfile_seqs_file );
            if ( !ForesterUtil.isEmpty( error1 ) ) {
                ForesterUtil.fatalError( PRG_NAME, error1 );
            }
            System.out.println();
            if ( input_seqs_fasta_like ) {
                System.out.println( "Input format          : Fasta" );
            }
            else {
                System.out.println( "Input format          : Phylip like" );
            }
            System.out.println( "Input file            : " + input_seqs_file_str );
            if ( output_format == MSA_FORMAT.FASTA ) {
                System.out.println( "Output format         : Fasta" );
            }
            else if ( output_format == MSA_FORMAT.NEXUS ) {
                System.out.println( "Output format         : Nexus" );
            }
            else if ( output_format == MSA_FORMAT.PHYLIP ) {
                System.out.println( "Output format         : Phylip" );
            }
            System.out.println( "Output file           : " + output_seqs_file_str );
            System.out.println( "Shorten names         : " + normalize_identifiers );
            if ( normalize_identifiers ) {
                System.out.println( "Identifier map        : " + output_map_file_str );
            }
            final List<MolecularSequence> input_seqs;
            final FileInputStream is = new FileInputStream( input_seqs_file );
            if ( FastaParser.isLikelyFasta( input_seqs_file ) ) {
                input_seqs = FastaParser.parse( is );
            }
            else {
                input_seqs = GeneralMsaParser.parseSeqs( is );
            }
            if ( input_seqs == null ) {
                ForesterUtil.fatalError( PRG_NAME, "failed to read input sequences" );
            }
            if ( input_seqs.size() < 1 ) {
                ForesterUtil.fatalError( PRG_NAME, "input seems to be devoid of sequences" );
            }
            final BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
            for( final MolecularSequence seq : input_seqs ) {
                stats.addValue( seq.getLength() );
            }
            System.out.println( "Number of sequences   : " + input_seqs.size() );
            if ( !ForesterUtil.isEqual( stats.getMin(), stats.getMax() ) ) {
                System.out.println( "Sequence lenght min   : " + ( int ) stats.getMin() );
                System.out.println( "Sequence lenght max   : " + ( int ) stats.getMax() );
                if ( input_seqs.size() > 2 ) {
                    System.out.println( "Sequence length median: " + ( int ) stats.median() );
                }
                if ( ( output_format == MSA_FORMAT.NEXUS ) || ( output_format == MSA_FORMAT.PHYLIP ) ) {
                    ForesterUtil.fatalError( PRG_NAME,
                                             "Input is not an alignment, cannot write in Nexus or Phylip format" );
                }
            }
            else {
                System.out.println( "Alignment length      : " + ( int ) stats.getMax() );
            }
            final List<MolecularSequence> output_seqs = new ArrayList<>();
            int counter = 0;
            final BufferedWriter output_map_writer;
            if ( normalize_identifiers ) {
                output_map_writer = ForesterUtil.createBufferedWriter( output_map_file_str );
            }
            else {
                output_map_writer = null;
            }
            for( final MolecularSequence seq : input_seqs ) {
                final String new_name;
                if ( normalize_identifiers ) {
                    new_name = modify_name( seq.getIdentifier(), counter++, output_map_writer );
                }
                else {
                    new_name = seq.getIdentifier();
                }
                final MolecularSequence ns = BasicSequence.createGeneralSequence( new_name,
                                                                                  seq.getMolecularSequenceAsString() );
                output_seqs.add( ns );
            }
            System.out.println();
            if ( normalize_identifiers ) {
                output_map_writer.flush();
                output_map_writer.close();
                System.out.println( "Wrote                 : " + output_map_file );
            }
            final BufferedWriter seq_writer = ForesterUtil.createBufferedWriter( outfile_seqs_file );
            if ( ( output_format == MSA_FORMAT.NEXUS ) || ( output_format == MSA_FORMAT.PHYLIP ) ) {
                final Msa m = BasicMsa.createInstance( output_seqs );
                m.write( seq_writer, output_format );
            }
            else if ( output_format == MSA_FORMAT.FASTA ) {
                SequenceWriter.writeSeqs( output_seqs, seq_writer, SEQ_FORMAT.FASTA, 60 );
            }
            seq_writer.flush();
            seq_writer.close();
            System.out.println( "Wrote                 : " + outfile_seqs_file );
            System.out.println();
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
        System.out.println( PRG_NAME + " [options] <input sequences file> [output sequences file] [output map file]" );
        System.out.println();
        System.out.println( " options:" );
        System.out.println( "  -" + OUTPUT_FORMAT_OPTION + "=<format>: output format: " + OUTPUT_FORMAT_FASTA_L + " or "
                + OUTPUT_FORMAT_FASTA + " for Fasta (default), " + OUTPUT_FORMAT_PHYLIP_L + " or "
                + OUTPUT_FORMAT_PHYLIP + " for Phylip, " + OUTPUT_FORMAT_NEXUS_L + " or " + OUTPUT_FORMAT_NEXUS
                + " for Nexus" );
        System.out.println( "  -" + ID_NORM_OPTION + "         : to replace sequence names with short(er) identifiers" );
        System.out.println();
        System.out.println( "Example:" );
        System.out.println();
        System.out.println( " " + PRG_NAME + " -s -o=p my_seqs.fasta" );
        System.out.println();
    }
}
