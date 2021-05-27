
package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
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

public class pango_select_loc {

    private final static String PRG_NAME          = "pango_select_loc";
    private static final String PRG_DATE          = "2021-05-26";
    private static final String PRG_VERSION       = "1.0.0";
    private static final String TARGET_LINEAGES[] = new String[] { "A.23.1", "A.27", "B.1.1.318", "B.1.1.519",
            "B.1.1.7", "B.1.351", "B.1.427", "B.1.429", "B.1.525", "B.1.526", "B.1.526.1", "B.1.526.2", "B.1.526.3",
            "B.1.617.1", "B.1.617.2", "C.37", "P.1", "P.2", "P.3", "R.1" };
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( args.length != 3 ) {
            System.out.println( "\nWrong number of arguments, expected: lineage_file fasta_seqs_file outfile\n" );
            System.exit( -1 );
        }
        final File lineage_file = new File( args[ 0 ] );
        final File infile = new File( args[ 1 ] );
        final String outfile_base = args[ 2 ];
        if ( !lineage_file.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + lineage_file + "] does not exist" );
        }
        if ( !infile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + infile + "] does not exist" );
        }
        BasicTable<String> lineage_table = null;
        try {
            lineage_table = BasicTableParser.parse( lineage_file, ',', false, false );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read [" + lineage_file + "] [" + e.getMessage() + "]" );
        }
        final SortedSet<String> target_lineages = new TreeSet<>();
        for( final String l : TARGET_LINEAGES ) {
            target_lineages.add( l );
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
        final SortedMap<String, List<MolecularSequence>> lin_seqs_map = new TreeMap<>();
        for( final String lineage : all_lineages_sorted ) {
            if ( target_lineages.contains( lineage ) ) {
                for( int i = 0; i < lineage_table.getNumberOfRows(); ++i ) {
                    final String lineage_from_table = lineage_table.getValue( 1, i );
                    if ( lineage_from_table.equals( lineage ) ) {
                        if ( id_to_seq.containsKey( lineage_table.getValue( 0, i ) ) ) {
                            if ( !lin_seqs_map.containsKey( lineage ) ) {
                                lin_seqs_map.put( lineage, new ArrayList<MolecularSequence>() );
                            }
                            lin_seqs_map.get( lineage ).add( id_to_seq.get( lineage_table.getValue( 0, i ) ) );
                        }
                    }
                }
            }
        }
        writeResult( outfile_base, lin_seqs_map );
        System.out.println();
        System.out.println( "[" + PRG_NAME + "] OK" );
        System.out.println();
    }

    private static void writeResult( final String outfile_base,
                                     final SortedMap<String, List<MolecularSequence>> lin_seqs_map ) {
        for( final Entry<String, List<MolecularSequence>> entry : lin_seqs_map.entrySet() ) {
            final String lin = entry.getKey();
            final List<MolecularSequence> sequences = entry.getValue();
            final File outfile = new File( outfile_base + "_" + lin.replace( '.', '_' ) );
            if ( sequences.size() > 2 ) {
                if ( outfile.exists() ) {
                    ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
                }
                try {
                    SequenceWriter.writeSeqs( sequences, outfile, SEQ_FORMAT.FASTA, 80 );
                }
                catch ( final IOException e ) {
                    e.printStackTrace();
                    System.exit( -1 );
                }
                System.out.println( "wrote: [" + outfile + "] with " + sequences.size() + " sequences" );
            }
            else {
                System.out.println( "  ignored: [" + outfile + "] with " + sequences.size() + " sequences" );
            }
        }
    }
}
