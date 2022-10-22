
package org.forester.application;

import java.io.File;
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
import org.forester.msa.Msa;
import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;

public class pango_select_4 {

    private final static String PRG_NAME    = "pango_select_4";
    private static final String PRG_DATE    = "2022-05-10";
    private static final String PRG_VERSION = "1.0.0";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( args.length != 2 ) {
            System.out.println( "\nWrong number of arguments, expected: indir outfile\n" );
            System.exit( -1 );
        }
        final File indir = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        final SortedSet<String> too_many = new TreeSet<>();
        if ( !indir.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + indir + "] does not exist" );
        }
        if ( !indir.isDirectory() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + indir + "] is not a directory" );
        }
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        final List<MolecularSequence> out_seqs = new ArrayList<>();
        for( final File f : indir.listFiles() ) {
            if ( !f.isDirectory() && f.getName().endsWith( ".fasta" ) ) {
                final File seqs_per_lineage_file = f.getAbsoluteFile();
                final String name = seqs_per_lineage_file.getName();
                final String afa = name.substring( 0, name.length() - 6 ) + ".afa";
                System.out.println( "fasta: " + seqs_per_lineage_file );
                System.out.println( "afa  : " + afa );
                System.out.println();
                calcDistances_1_2_d( new File( seqs_per_lineage_file.getParent() + "/" + afa ),
                                     seqs_per_lineage_file,
                                     out_seqs );
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
        System.out.println();
        System.out.println( "[" + PRG_NAME + "] wrote: [" + outfile + "] with " + out_seqs.size() + " sequences" );
        System.out.println( "[" + PRG_NAME + "] OK" );
        System.out.println();
    }

    private static void calcDistances_1_2_d( final File msa_file,
                                             final File seqs_per_lineage_file,
                                             final List<MolecularSequence> out_seqs ) {
        List<MolecularSequence> seqs_per_lineage = null;
        try {
            seqs_per_lineage = FastaParser.parse( seqs_per_lineage_file );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read fasta file: " + e );
        }
        if ( seqs_per_lineage.size() <= 4 ) {
            for( final MolecularSequence ms : seqs_per_lineage ) {
                out_seqs.add( ms );
            }
            return;
        }
        Msa msa = null;
        try {
            msa = FastaParser.parseMsa( msa_file );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read msa: " + e );
        }
        final BasicSymmetricalDistanceMatrix m = PairwiseDistanceCalculator.calcFractionalDissimilarities( msa );
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
}
