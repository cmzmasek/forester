//

package org.forester.application;

import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.forester.io.parsers.FastaParser;
import org.forester.msa.Msa;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.Sequence;
import org.forester.util.ForesterUtil;

public class aa {

    public static void main( final String args[] ) {
        try {
            System.out.println( "STARTING..." );
            final List<Sequence> orig = FastaParser
                    .parse( new FileInputStream( "C:\\Users\\zma\\Desktop\\RRMa_domains_ext_20.fasta" ) );
            final Msa msa = FastaParser.parseMsa( new FileInputStream( "C:\\Users\\zma\\Desktop\\test3_sorted.fasta" ) );
            final Set<Sequence> all_found_seqs = new HashSet<Sequence>();
            for( int i = 0; i < msa.getNumberOfSequences(); ++i ) {
                final String id = msa.getIdentifier( i );
                final String id_ = id.substring( 0, id.indexOf( "_" ) );
                final String range = id.substring( id.indexOf( "[" ) + 1, id.indexOf( "]" ) );
                //System.out.println( i + ": " + id + "=>" + id_ + " " + range );
                if ( ForesterUtil.isEmpty( id_ ) ) {
                    System.out.println( "ERROR: id is empty for: " + id );
                    System.exit( -1 );
                }
                if ( ForesterUtil.isEmpty( range ) ) {
                    System.out.println( "ERROR: range is empty for: " + id );
                    System.exit( -1 );
                }
                int found = 0;
                final List<Sequence> found_seqs = new ArrayList<Sequence>();
                for( final Sequence orig_seq : orig ) {
                    final String orig_seq_id = orig_seq.getIdentifier();
                    if ( ( orig_seq_id.indexOf( id_ ) >= 0 ) && ( orig_seq_id.indexOf( "[" + range + "]" ) >= 0 ) ) {
                        found++;
                        found_seqs.add( orig_seq );
                    }
                }
                if ( found > 0 ) {
                    for( final Sequence found_seq : found_seqs ) {
                        if ( found_seq.getLength() >= 85 ) {
                            all_found_seqs.add( BasicSequence.createAaSequence( id, found_seq
                                    .getMolecularSequenceAsString() ) );
                        }
                    }
                    if ( found > 1 ) {
                        System.out.println( i + ": " + id + "=>" + id_ + " " + range );
                        System.out.println( "  found: " + found );
                        for( final Sequence found_seq : found_seqs ) {
                            System.out.println( found_seq.toString() );
                        }
                    }
                }
                else {
                    System.out.println( "ERROR: not found: " + id );
                    System.exit( -1 );
                }
            }
            final String fasta_ary[] = new String[ all_found_seqs.size() ];
            int i = 0;
            for( final Sequence sequence : all_found_seqs ) {
                fasta_ary[ i ] = ">" + sequence.getIdentifier() + "\n" + sequence.getMolecularSequenceAsString();
                System.out.println( sequence );
                i++;
            }
            Arrays.sort( fasta_ary );
            for( final String element : fasta_ary ) {
                System.out.println( element );
            }
            System.out.println( "DONE." );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
        }
    }
}
