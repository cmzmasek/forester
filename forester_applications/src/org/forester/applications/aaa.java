
package org.forester.applications;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.FastaParser;
import org.forester.sequence.Sequence;
import org.forester.util.EasyWriter;
import org.forester.util.ForesterUtil;

public class aaa {

    public final static Pattern GN_PATTERN    = Pattern.compile( "GN=(\\S+)\\s" );     //use w+ instead of S+ for more stringent setting.
    public final static Pattern RANGE_PATTERN = Pattern.compile( "\\[(\\d+-\\d+)\\]" ); //use w+ instead of S+ for more stringent setting.
    public final static int     MIN_LENGTH    = 85;

    public static void main( final String args[] ) {
        try {
            final EasyWriter out = ( EasyWriter ) ForesterUtil.createEasyWriter( "aaa_out" );
            System.out.println( "STARTING..." );
            final List<Sequence> too_short = new ArrayList<Sequence>();
            final List<Sequence> orig = FastaParser
                    .parse( new FileInputStream( "C:\\Users\\zma\\Desktop\\RRMa_domains_ext_20_2.fasta" ) );
            final int initial_number = orig.size();
            final List<String> new_seqs = new ArrayList<String>();
            for( final Sequence seq : orig ) {
                if ( seq.getLength() < MIN_LENGTH ) {
                    too_short.add( seq );
                    continue;
                }
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
            final List<String> duplicate_gn_ra = new ArrayList<String>();
            final List<String> duplicate_mol_seq = new ArrayList<String>();
            final List<String> new_seqs_unique = new ArrayList<String>();
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
                        new_seqs_unique.add( seq );
                        unique_counter++;
                    }
                    else {
                        duplicate_mol_seq.add( seq );
                    }
                }
                else {
                    duplicate_gn_ra.add( seq );
                }
            }
            String prev_gn = "___";
            boolean is_first = true;
            List<String> seqs_from_same_protein = new ArrayList<String>();
            for( final String seq : new_seqs_unique ) {
                final Matcher matcher_gn = GN_PATTERN.matcher( seq );
                matcher_gn.find();
                final String gn = matcher_gn.group( 1 );
                if ( !prev_gn.equals( gn ) && !is_first ) {
                    doit( seqs_from_same_protein, out );
                    seqs_from_same_protein = new ArrayList<String>();
                }
                prev_gn = gn;
                is_first = false;
                seqs_from_same_protein.add( seq );
            }
            doit( seqs_from_same_protein, out );
            out.println( "" );
            out.println( "" );
            out.println( "Removed because same GN and region:" );
            for( final String s : duplicate_gn_ra ) {
                out.println( s );
            }
            out.println( "" );
            out.println( "" );
            out.println( "Removed because identical mol sequence:" );
            for( final String s : duplicate_mol_seq ) {
                out.println( s );
            }
            out.println( "" );
            out.println( "" );
            out.println( "Removed because too short:" );
            for( final Sequence s : too_short ) {
                out.println( s.toString() );
            }
            out.println( "" );
            out.println( "" );
            out.println( "initial:" + initial_number );
            out.println( "ignored because shorter than " + MIN_LENGTH + "aa: " + too_short.size() );
            out.println( "unique   : " + unique_counter );
            out.println( "unique   : " + new_seqs_unique.size() );
            out.println( "duplicate because gn and range same: " + duplicate_gn_ra.size() );
            out.println( "duplicate because mol seq same     : " + duplicate_mol_seq.size() );
            out.flush();
            out.close();
            System.out.println( "DONE " );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
        }
    }

    private static void doit( final List<String> same_protein_seqs, final EasyWriter out ) throws IOException {
        final int count = same_protein_seqs.size();
        if ( count == 1 ) {
            out.println( same_protein_seqs.get( 0 ) );
        }
        else {
            int c = 1;
            for( final String s : same_protein_seqs ) {
                out.println( new StringBuffer( s ).insert( s.indexOf( "|" ), "__" + c + "_OF_" + count ).toString() );
                c++;
            }
        }
    }
}
