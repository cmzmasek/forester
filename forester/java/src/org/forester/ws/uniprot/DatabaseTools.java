package org.forester.ws.uniprot;

import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class DatabaseTools {
    //The format for GenBank Accession numbers are:
    //Nucleotide: 1 letter + 5 numerals OR 2 letters + 6 numerals
    //Protein:    3 letters + 5 numerals
    //http://www.ncbi.nlm.nih.gov/Sequin/acc.html
    
    private final static Pattern GENBANK_NUCLEOTIDE_AC_PATTERN_1 = Pattern
    .compile( "^.*[^a-zA-Z0-9]?([A-Z]\\d{5})[^a-zA-Z0-9]?" );
    
    private final static Pattern GENBANK_NUCLEOTIDE_AC_PATTERN_2 = Pattern
    .compile( "^.*[^a-zA-Z0-9]?([A-Z]{2}\\d{6})[^a-zA-Z0-9]?" );

    private final static Pattern GENBANK_PROTEIN_AC_PATTERN = Pattern
    .compile( "^.*[^a-zA-Z0-9]?([A-Z]{3}\\d{5})[^a-zA-Z0-9]?" );

    
    
    private final static boolean DEBUG              = false;

    /**
     * Returns null if no match.
     * 
     * @param query
     * @param db 
     * @return
     */
    static public String parseGenbankAccessor( final String query ) {
        Matcher m = GENBANK_NUCLEOTIDE_AC_PATTERN_1.matcher( query );
        if ( m.lookingAt() ) {
            return m.group( 1 );
        }
        else {
             m = GENBANK_NUCLEOTIDE_AC_PATTERN_2.matcher( query );
            if ( m.lookingAt() ) {
                return m.group( 1 );
            } 
            else {
                m = GENBANK_PROTEIN_AC_PATTERN.matcher( query );
                if ( m.lookingAt() ) {
                    return m.group( 1 );
                }
                else {
                    return null;
                }
            }
        }
    }

    static String extract( final String target, final String a, final String b ) {
        final int i_a = target.indexOf( a );
        final int i_b = target.indexOf( b );
        if ( ( i_a < 0 ) || ( i_b < i_a ) ) {
            throw new IllegalArgumentException( "attempt to extract from [" + target + "] between [" + a + "] and ["
                    + b + "]" );
        }
        return target.substring( i_a + a.length(), i_b ).trim();
    }



    static String extract( final String target, final String a ) {
        final int i_a = target.indexOf( a );
        return target.substring( i_a + a.length() ).trim();
    }

}
