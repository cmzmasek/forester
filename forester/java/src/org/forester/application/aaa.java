
package org.forester.application;

import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.FastaParser;
import org.forester.sequence.Sequence;
import org.forester.util.ForesterUtil;

public class aaa {

    public final static Pattern GN_PATTERN    = Pattern.compile( "GN=(\\S+)\\s" );     //use w+ instead of S+ for more stringent setting.
    public final static Pattern RANGE_PATTERN = Pattern.compile( "\\[(\\d+-\\d+)\\]" ); //use w+ instead of S+ for more stringent setting.

    public static void main( final String args[] ) {
        try {
            System.out.println( "STARTING..." );
            final List<Sequence> orig = FastaParser
                    .parse( new FileInputStream( "C:\\Users\\zma\\Desktop\\RRMa_domains_ext_20_2.fasta" ) );
            final List<String> new_seqs = new ArrayList<String>();
            for( final Sequence seq : orig ) {
                final Matcher matcher = GN_PATTERN.matcher( seq.getIdentifier() );
                String gn = "";
                if ( matcher.find() ) {
                    gn = matcher.group( 1 );
                }
                else {
                    System.out.println( "ERROR: no gene for: " + seq.getIdentifier() );
                    System.exit( -1 );
                }
                new_seqs.add( ">" + gn + "|" + seq.getIdentifier() + "\n" + seq.getMolecularSequenceAsString() );
            }
            final Set<String> gn_ra_set = new HashSet<String>();
            final Set<String> mol_seq_set = new HashSet<String>();
            Collections.sort( new_seqs );
            int unique_counter = 0;
            int duplicate_counter_gn_ra = 0;
            int duplicate_counter_mol_seq = 0;
            String prev_gn = "____";
            for( final String seq : new_seqs ) {
                final Matcher matcher_ra = RANGE_PATTERN.matcher( seq );
                final Matcher matcher_gn = GN_PATTERN.matcher( seq );
                String range = "";
                if ( matcher_ra.find() ) {
                    range = matcher_ra.group( 1 );
                }
                else {
                    System.out.println( "ERROR: no range for: " + seq );
                    System.exit( -1 );
                }
                matcher_gn.find();
                final String gn = matcher_gn.group( 1 );
                final String gn_ra = gn + "_" + range;
                if ( !gn_ra_set.contains( gn_ra ) ) {
                    gn_ra_set.add( gn_ra );
                    final String mol_seq = seq.split( "\n" )[ 1 ];
                    if ( !mol_seq_set.contains( mol_seq ) ) {
                        mol_seq_set.add( mol_seq );
                        if ( prev_gn.equals( gn ) ) {
                            int count = same_protein_seqs.size();
                            if ( count == 1 ) {
                                System.out.println( seq );
                            }
                            else {
                                int c = 1;
                                for( final String s : same_protein_seqs ) {
                                    System.out.println( new StringBuffer( s ).insert( s.indexOf( "|" ),
                                                                                      "__" + c + "_OF_" + count )
                                            .toString() );
                                    c++;
                                }
                            }
                        }
                        prev_gn = gn;
                        System.out.println( seq );
                        unique_counter++;
                    }
                    else {
                        duplicate_counter_mol_seq++;
                    }
                }
                else {
                    duplicate_counter_gn_ra++;
                }
            }
            System.out.println( "unique   : " + unique_counter );
            System.out.println( "duplicate because gn and range same: " + duplicate_counter_gn_ra );
            System.out.println( "duplicate because mol seq same     : " + duplicate_counter_mol_seq );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
        }
    }
}
