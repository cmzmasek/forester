
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

import org.forester.evoinference.distance.PairwiseDistanceCalculator;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.io.parsers.FastaParser;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.msa.Mafft;
import org.forester.msa.Msa;
import org.forester.msa.MsaInferrer;
import org.forester.sequence.MolecularSequence;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

public class pango_select_2 {

    private final static String PRG_NAME             = "pango_select";
    private static final String PRG_DATE             = "2021-06-17";
    private static final String PRG_VERSION          = "0.0.1";
    private final static int    MAX_SEQS_PER_LINEAGE = 50;
    private final static String PATH_TO_MAFFT        = "/usr/local/bin/mafft";
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
            ForesterUtil.fatalError( PRG_NAME, "failed to parset [" + outfile + "]: " + e );
        }
        final SortedMap<String, MolecularSequence> id_to_seq = new TreeMap<>();
        for( final MolecularSequence seq : seqs ) {
            id_to_seq.put( seq.getIdentifier(), seq );
        }
        final List<MolecularSequence> out_seqs = new ArrayList<>();
        for( final String lineage : all_lineages_sorted ) {
            final List<MolecularSequence> seqs_per_lineage = new ArrayList<>();
            int counter = 0;
            for( int i = 0; i < lineage_table.getNumberOfRows(); ++i ) {
                final String lineage_from_table = lineage_table.getValue( 1, i );
                if ( ( counter < MAX_SEQS_PER_LINEAGE ) && lineage_from_table.equals( lineage ) ) {
                    System.out.println( lineage );
                    ++counter;
                    seqs_per_lineage.add( id_to_seq.get( lineage_table.getValue( 0, i ) ) );
                }
            }
            if ( seqs_per_lineage.size() > 3 ) {
                calcDistances( PATH_TO_MAFFT, out_seqs, seqs_per_lineage );
            }
            else {
                for( final MolecularSequence seq : seqs_per_lineage ) {
                    out_seqs.add( seq );
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

    private static void calcDistances( final String path_to_mafft,
                                       final List<MolecularSequence> out_seqs,
                                       final List<MolecularSequence> seqs_per_lineage ) {
        Msa msa = null;
        try {
            final MsaInferrer mafft = Mafft.createInstance( path_to_mafft );
            final List<String> opts = new ArrayList<>();
            try {
                msa = mafft.infer( seqs_per_lineage, opts );
            }
            catch ( final InterruptedException e ) {
                ForesterUtil.fatalError( PRG_NAME, "failed to execute MAFFT: " + e );
            }
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to execute MAFFT: " + e );
        }
        final BasicSymmetricalDistanceMatrix m = PairwiseDistanceCalculator.calcFractionalDissimilarities( msa );
        //System.out.println( m.toString() );
        //
        double min_sum = Double.MAX_VALUE;
        double max_1_sum = 0;
        double max_2_sum = 0;
        String min_id = "";
        String max_1_id = "";
        String max_2_id = "";
        for( int row = 0; row < m.getSize(); ++row ) {
            double sum = 0;
            for( int col = 0; col < m.getSize(); ++col ) {
                sum += m.getValue( col, row );
            }
            if ( sum < min_sum ) {
                min_sum = sum;
                min_id = m.getIdentifier( row );
            }
            if ( sum > max_1_sum ) {
                max_2_sum = max_1_sum;
                max_2_id = max_1_id;
                max_1_sum = sum;
                max_1_id = m.getIdentifier( row );
            }
            else if ( sum > max_2_sum ) {
                max_2_sum = sum;
                max_2_id = m.getIdentifier( row );
            }
        }
        System.out.println( "MIN : " + min_id + " = " + min_sum );
        System.out.println( "MAX1: " + max_1_id + " = " + max_1_sum );
        System.out.println( "MAX2: " + max_2_id + " = " + max_2_sum );
        for( final MolecularSequence seq : seqs_per_lineage ) {
            if ( seq.getIdentifier().equals( min_id ) || seq.getIdentifier().equals( max_1_id )
                    || seq.getIdentifier().equals( max_2_id ) ) {
                out_seqs.add( seq );
            }
        }
    }
}
