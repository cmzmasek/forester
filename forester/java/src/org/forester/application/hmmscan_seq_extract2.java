
package org.forester.application;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser;
import org.forester.io.parsers.HmmscanPerDomainTableParser.INDIVIDUAL_SCORE_CUTOFF;
import org.forester.io.writers.SequenceWriter;
import org.forester.protein.Protein;
import org.forester.sequence.MolecularSequence;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class hmmscan_seq_extract2 {

    final static private String  PRG_NAME                              = "hmmscan_seq_extract";
    final static private String  PRG_VERSION                           = "1.0.0";
    final static private String  PRG_DATE                              = "20221101";
    final static private String  PRG_DESC                              = "hmmscan seq extract";
    final static private String  E_MAIL                                = "phyloxml@gmail.com";
    final static private String  WWW                                   = "https://github.com/cmzmasek/forester";
    final static private String  HELP_OPTION_1                         = "help";
    final static private String  HELP_OPTION_2                         = "h";
    private static final String  MAX_I_E_VALUE_OPTION                  = "ie";
    private static final String  MIN_REL_ENV_LENGTH_RATIO_OPTION       = "mrel";
    private static final int     HMMSCAN_PARSE_MAX_OVERLAP             = 5;
    private static final boolean HMMSCAN_PARSE_IGNORE_ENGULFED_DOMAINS = true;
    public static void main( final String args[] ) {
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
        if ( ( cla.getNumberOfNames() != 4 ) ) {
            print_help();
            System.exit( -1 );
        }
        final List<String> allowed_options = new ArrayList<>();
        allowed_options.add( MAX_I_E_VALUE_OPTION );
        allowed_options.add( MIN_REL_ENV_LENGTH_RATIO_OPTION );
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        double rel_env_length_ratio_cutoff = -1;
        double ie_value_max = -1;
        if ( cla.isOptionSet( MIN_REL_ENV_LENGTH_RATIO_OPTION ) ) {
            try {
                rel_env_length_ratio_cutoff = cla.getOptionValueAsDouble( MIN_REL_ENV_LENGTH_RATIO_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, "no acceptable value for min rel env length ratio" );
            }
        }
        if ( cla.isOptionSet( MAX_I_E_VALUE_OPTION ) ) {
            try {
                ie_value_max = cla.getOptionValueAsDouble( MAX_I_E_VALUE_OPTION );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, "no acceptable value for E-value maximum" );
            }
        }
        final File hmmscan_file = cla.getFile( 0 );
        final HmmscanPerDomainTableParser hmmscan_parser = new HmmscanPerDomainTableParser( hmmscan_file,
                                                                                            "unknown",
                                                                                            INDIVIDUAL_SCORE_CUTOFF.NONE );
        final File fasta_file = cla.getFile( 1 );
        final File output_dir = cla.getFile( 3 );
        final File inmap = cla.getFile( 2 );
        System.out.println( "Hmmscan output file                   : " + hmmscan_file );
        System.out.println( "Fasta file (used as input for hmmscan): " + fasta_file );
        System.out.println( "Maping file (DA SOG-name taxonomy)    : " + inmap );
        System.out.println( "Output directory                      : " + output_dir );
        System.out.println( "Ie value max cutoff                   : " + ie_value_max );
        System.out.println( "Rel env length ratio cutoff           : " + rel_env_length_ratio_cutoff );
        System.out.println( "Max domain overlap cutoff             : " + HMMSCAN_PARSE_MAX_OVERLAP );
        System.out.println( "Ignore engulfed domains               : " + HMMSCAN_PARSE_IGNORE_ENGULFED_DOMAINS );
        System.out.println();
        System.out.println();
        BasicTable<String> table = null;
        try {
            table = BasicTableParser.parse( inmap, '\t' );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read \"" + inmap + "\" [" + e.getMessage() + "]" );
        }
        if ( table.getNumberOfColumns() < 3 ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "table has " + table.getNumberOfColumns() + " columns, expected at least 3" );
        }
        if ( table.getNumberOfRows() < 1 ) {
            ForesterUtil.fatalError( PRG_NAME, "table has no rows (i.e. is empty)" );
        }
        ForesterUtil.programMessage( PRG_NAME,
                                     "Successfully read in table with " + table.getNumberOfColumns() + " columns and "
                                             + table.getNumberOfRows() + " rows" );
        final SortedMap<String, String> da_to_name_map = table.getColumnsAsMap( 0, 1 );
        final SortedMap<String, String> da_to_taxo_map = table.getColumnsAsMap( 0, 2 );
        if ( output_dir.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + output_dir + "] already exists" );
        }
        if ( !fasta_file.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + fasta_file + "] does not exists" );
        }
        if ( !hmmscan_file.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + hmmscan_file + "] does not exists" );
        }
        output_dir.mkdirs();
        final Map<String, String> id_to_seq_map = readInFastaFile( fasta_file );
        System.out.println( "FASTA file: " + id_to_seq_map.size() + " entries" );
        if ( ie_value_max > 0.0 ) {
            hmmscan_parser.setIEValueMaximum( ie_value_max );
        }
        if ( rel_env_length_ratio_cutoff > 0.0 ) {
            hmmscan_parser.setRelEnvLengthRatioCutoff( rel_env_length_ratio_cutoff );
        }
        hmmscan_parser.setIgnoreEngulfedDomains( HMMSCAN_PARSE_IGNORE_ENGULFED_DOMAINS );
        hmmscan_parser.setMaxAllowedOverlap( HMMSCAN_PARSE_MAX_OVERLAP );
        List<Protein> hmmscan_scan_result = null;
        try {
            hmmscan_scan_result = hmmscan_parser.parse();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        catch ( final Exception e ) {
            ForesterUtil.unexpectedFatalError( PRG_NAME, e.getMessage(), e );
        }
        System.out.println();
        System.out.println();
        final SortedMap<String, SortedSet<String>> da_to_ids_map = new TreeMap<>();
        for( final Protein protein : hmmscan_scan_result ) {
            final String da = protein.toDomainArchitectureString( "--", 1 );
            if ( !da_to_ids_map.containsKey( da ) ) {
                da_to_ids_map.put( da, new TreeSet<String>() );
            }
            da_to_ids_map.get( da ).add( protein.getProteinId().getId() );
            System.out.print( protein.getProteinId().getId() );
            System.out.print( '\t' );
            System.out.print( da );
            System.out.println();
        }
        /////
        final SortedMap<String, Integer> not_found = new TreeMap<>();
        final SortedMap<String, Integer> found = new TreeMap<>();
        System.out.println();
        System.out.println();
        for( final Protein protein : hmmscan_scan_result ) {
            final String da = protein.toDomainArchitectureString( "--", 1 );
            if ( da_to_name_map.containsKey( da ) && da_to_taxo_map.containsKey( da ) ) {
                final String da_with_name = da + "\t" + da_to_name_map.get( da ) + '_' + da_to_taxo_map.get( da );
            }
            else {
                if ( not_found.containsKey( da ) ) {
                    not_found.put( da, not_found.get( da ) + 1 );
                }
                else {
                    not_found.put( da, 1 );
                }
            }
        }
        /////
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println( "FOUND DAs [" + da_to_ids_map.size() + "]" );
        System.out.println();
        for( final Entry<String, SortedSet<String>> entry : da_to_ids_map.entrySet() ) {
            final String da = entry.getKey();
            System.out.print( da + '\t' + entry.getValue() );
            System.out.println();
        }
        //  System.out.println( "Found "  + protein_id);
        int counter = 0;
        for( final Entry<String, SortedSet<String>> entry : da_to_ids_map.entrySet() ) {
            final String da = entry.getKey();
            System.out.println( "DA " + da );
            try {
                String da_with_name = "";
                if ( da_to_name_map.containsKey( da ) && da_to_taxo_map.containsKey( da ) ) {
                    da_with_name = da_to_name_map.get( da ) + '_' + da_to_taxo_map.get( da );
                    da_with_name=  da_with_name.replace('/', '~').replace(' ', '_').replace('(', '_').replace(')', '_');
                }
                else {
                    //da_with_name = "__" + da;
                    da_with_name = da;
                }
                String out_name = output_dir + "/" + da_with_name + ".fasta";
                if ( da.length() > 80 ) {
                    ++counter;
                    out_name = output_dir + "/" + da.substring( 0, 80 ) + "__" + counter + ".fasta";
                }
                final FileWriter fw = new FileWriter( out_name );
                final BufferedWriter bw = new BufferedWriter( fw );
                final SortedSet<String> ids = entry.getValue();
                for( final String protein_id : ids ) {
                    if ( id_to_seq_map.containsKey( protein_id ) ) {
                        System.out.println( "Found " + protein_id );
                        bw.write( SequenceWriter.toFasta( protein_id, id_to_seq_map.get( protein_id ), 80 )
                                .toString() );
                        bw.newLine();
                    }
                    else {
                        ForesterUtil.fatalError( PRG_NAME, "NOT found: " + protein_id );
                    }
                }
                bw.close();
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, "failed to write fasta file:" + e.getLocalizedMessage() );
            }
        }
    }

    private static Map<String, String> readInFastaFile( final File fasta_file ) {
        List<MolecularSequence> seqs_from_fasta_file = null;
        try {
            seqs_from_fasta_file = FastaParser.parse( new FileInputStream( fasta_file ) );
        }
        catch ( final IOException e1 ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read fasta file:" + e1.getLocalizedMessage() );
        }
        final Map<String, String> id_to_seq_map = new HashMap<>();
        for( final MolecularSequence seq : seqs_from_fasta_file ) {
            final String id = seq.getIdentifier();
            final String proc_id = id.split( " " )[ 0 ];
            id_to_seq_map.put( proc_id, seq.getMolecularSequenceAsString() );
        }
        return id_to_seq_map;
    }

    private static void print_help() {
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME
                + " [options] <hmmscan output> <fasta file (used as input for hmmscan)> <maping file (DA SOG-name taxonomy), tab separated>" );
        System.out.println();
        System.out.println( " options:" );
        System.out.println( MAX_I_E_VALUE_OPTION + ": max (inclusive) iE-value" );
        System.out.println( MIN_REL_ENV_LENGTH_RATIO_OPTION + ": min (inclusive) relative envelope length ratio" );
    }
}