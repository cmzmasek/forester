
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

    private final static String PRG_NAME                     = "pango_select_2";
    private static final String PRG_DATE                     = "2021-06-26";
    private static final String PRG_VERSION                  = "2.0.0";
    private final static int    MAX_SEQS_PER_LINEAGE         = 1000;
    private final static int    MAX_SEQS_PER_LINEAGE_SPECIAL = 1500;
    private final static String PATH_TO_MAFFT                = "/usr/local/bin/mafft";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( args.length != 3 ) {
            System.out.println( "\nWrong number of arguments, expected: lineage_file fasta_seqs_file outfile\n" );
            System.exit( -1 );
        }
        System.out.println( "MAX_SEQS_PER_LINEAGE        : " + MAX_SEQS_PER_LINEAGE );
        System.out.println( "MAX_SEQS_PER_LINEAGE_SPECIAL: " + MAX_SEQS_PER_LINEAGE_SPECIAL );
        final File lineage_file = new File( args[ 0 ] );
        final File infile = new File( args[ 1 ] );
        final File outfile = new File( args[ 2 ] );
        final SortedSet<String> too_many = new TreeSet<>();
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
        final SortedSet<String> one_seqs_set = new TreeSet<>();
        final SortedSet<String> no_seqs_set = new TreeSet<>();
        for( final String lineage : all_lineages_sorted ) {
            final List<MolecularSequence> seqs_per_lineage = new ArrayList<>();
            int counter = 0;
            System.out.println();
            System.out.print( lineage + ": " );
            for( int i = 0; i < lineage_table.getNumberOfRows(); ++i ) {
                final String lineage_from_table = lineage_table.getValue( 1, i );
                if ( lineage_from_table.equals( lineage ) ) {
                    if ( ( lineage_from_table.equals( "B.1.1.7" ) && ( counter < MAX_SEQS_PER_LINEAGE_SPECIAL ) )
                            || ( counter < MAX_SEQS_PER_LINEAGE ) ) {
                        if ( id_to_seq.containsKey( lineage_table.getValue( 0, i ) ) ) {
                            ++counter;
                            seqs_per_lineage.add( id_to_seq.get( lineage_table.getValue( 0, i ) ) );
                        }
                    }
                    else {
                        too_many.add( lineage );
                    }
                }
            }
            System.out.println( counter );
            System.out.println( "seqs per lineage: " + seqs_per_lineage.size() );
            if ( seqs_per_lineage.size() > 4 ) {
                calcDistances_1_2_d( PATH_TO_MAFFT, out_seqs, seqs_per_lineage );
            }
            else {
                for( final MolecularSequence seq : seqs_per_lineage ) {
                    if ( seq != null ) {
                        out_seqs.add( seq );
                    }
                }
            }
            if ( seqs_per_lineage.size() == 1 ) {
                one_seqs_set.add( lineage );
            }
            else if ( seqs_per_lineage.size() == 0 ) {
                no_seqs_set.add( lineage );
            }
        }
        try {
            SequenceWriter.writeSeqs( out_seqs, outfile, SEQ_FORMAT.FASTA, 80 );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
        System.out.println();
        System.out.println( "Lineages with only one sequence: " );
        System.out.println( one_seqs_set );
        System.out.println();
        System.out.println( "Lineages with no sequence: " );
        System.out.println( no_seqs_set );
        System.out.println();
        System.out.println( "Lineages with more than " + MAX_SEQS_PER_LINEAGE + " sequences (" + too_many.size()
                + "): " );
        for( final String lin : too_many ) {
            System.out.println( lin );
        }
        System.out.println();
        System.out.println( "[" + PRG_NAME + "] wrote: [" + outfile + "] with " + out_seqs.size() + " sequences" );
        System.out.println( "[" + PRG_NAME + "] OK" );
        System.out.println();
    }

    private static void calcDistances_1_2_d( final String path_to_mafft,
                                             final List<MolecularSequence> out_seqs,
                                             final List<MolecularSequence> seqs_per_lineage ) {
        Msa msa = null;
        try {
            final MsaInferrer mafft = Mafft.createInstance( path_to_mafft );
            final List<String> opts = new ArrayList<>();
            //opts.add( "--retree 1" );
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
        final SortedMap<Double, String> sum_to_id_map = new TreeMap<>();
        for( int row = 0; row < m.getSize(); ++row ) {
            double sum = 0;
            for( int col = 0; col < m.getSize(); ++col ) {
                sum += m.getValue( col, row );
            }
            sum_to_id_map.put( sum, m.getIdentifier( row ) );
        }
        System.out.println( sum_to_id_map );
        final List<String> ids = new ArrayList<>( sum_to_id_map.values() );
        final String min_1_id = ids.get( 0 );
        final String min_2_id = ids.get( 1 );
        final String max_1_id = ids.get( ids.size() - 1 );
        String max_d_id = null;
        double max_d_d = 0.0;
        for( int row = 0; row < m.getSize(); ++row ) {
            if ( m.getIdentifier( row ).equals( max_1_id ) ) {
                for( int col = 0; col < m.getSize(); ++col ) {
                    if ( ( m.getValue( col, row ) > max_d_d ) && !m.getIdentifier( col ).equals( min_1_id )
                            && !m.getIdentifier( col ).equals( min_2_id ) ) {
                        max_d_d = m.getValue( col, row );
                        max_d_id = m.getIdentifier( col );
                    }
                }
            }
        }
        System.out.println( " min 1 id: " + min_1_id );
        System.out.println( " min 2 id: " + min_2_id );
        System.out.println( " max 1 id: " + max_1_id );
        System.out.println( " max d id: " + max_d_id );
        for( final MolecularSequence seq : seqs_per_lineage ) {
            if ( seq.getIdentifier().equals( min_1_id ) || seq.getIdentifier().equals( min_2_id )
                    || seq.getIdentifier().equals( max_1_id ) || seq.getIdentifier().equals( max_d_id ) ) {
                System.out.println( "     adding seq: " + seq.getIdentifier() );
                out_seqs.add( seq );
            }
        }
    }

    private static void calcDistances( final String path_to_mafft,
                                       final List<MolecularSequence> out_seqs,
                                       final List<MolecularSequence> seqs_per_lineage ) {
        Msa msa = null;
        try {
            final MsaInferrer mafft = Mafft.createInstance( path_to_mafft );
            final List<String> opts = new ArrayList<>();
            //opts.add( "--retree 1" );
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
