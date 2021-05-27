
package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.MolecularSequence;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

public class pango_select {

    private final static String PRG_NAME             = "pango_select";
    private static final String PRG_DATE             = "2021-05-24";
    private static final String PRG_VERSION          = "1.0.0";
    private final static int    MAX_SEQS_PER_LINEAGE = 3;
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( args.length != 3 ) {
            System.out.println( "\nWrong number of arguments, expected: lineage_file fasta_seqs_file outfile\n" );
            System.exit( -1 );
        }
        final File lineage_file = new File( args[ 0 ] );
        final File infile = new File( args[ 1 ] );
        final File outfile = new File( args[ 2 ] );
        if ( !lineage_file.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + lineage_file + "] does not exist" );
        }
        if ( !infile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + infile + "] does not exist" );
        }
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        BasicTable<String> lineage_table = null;
        try {
            lineage_table = BasicTableParser.parse( lineage_file, ',', false, false );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read [" + lineage_file + "] [" + e.getMessage() + "]" );
        }
        final SortedSet<String> all_lineages_sorted = new TreeSet<>();
        for( int i = 0; i < lineage_table.getNumberOfRows(); ++i ) {
            final String lineage = lineage_table.getValue( 1, i );
            if ( !lineage.equalsIgnoreCase( "lineage" ) ) {
                all_lineages_sorted.add( lineage );
            }
        }
        System.out.println( "All lineages: " + all_lineages_sorted.size() );
        System.out.println( "All lineages: " + all_lineages_sorted );
        List<MolecularSequence> seqs = null;
        try {
            seqs = FastaParser.parse( new FileInputStream( infile ) );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
        final SortedMap<String, MolecularSequence> id_to_seq = new TreeMap<>();
        for( final MolecularSequence seq : seqs ) {
            id_to_seq.put( seq.getIdentifier(), seq );
        }
        final List<MolecularSequence> out_seqs = new ArrayList<>();
        for( final String lineage : all_lineages_sorted ) {
            int counter = 0;
            for( int i = 0; i < lineage_table.getNumberOfRows(); ++i ) {
                final String lineage_from_table = lineage_table.getValue( 1, i );
                if ( ( counter < MAX_SEQS_PER_LINEAGE ) && lineage_from_table.equals( lineage ) ) {
                    System.out.println( lineage );
                    ++counter;
                    out_seqs.add( id_to_seq.get( lineage_table.getValue( 0, i ) ) );
                }
            }
        }
        try {
            SequenceWriter.writeSeqs( out_seqs, outfile, SEQ_FORMAT.FASTA, 80 );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
        System.out.println( "[" + PRG_NAME + "] wrote: [" + outfile + "] with " + out_seqs.size() + " sequences" );
        System.out.println( "[" + PRG_NAME + "] OK" );
        System.out.println();
    }
}
