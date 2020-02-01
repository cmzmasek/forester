
package org.forester.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class reform_for_bats {

    private final static Pattern VIPR_PATTERN = Pattern
            .compile( "(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)" );

    public static void main( final String args[] ) throws FileNotFoundException, IOException {
        if ( args.length != 2 ) {
            System.out.println( "\nWrong number of arguments, expected: infile map\n" );
            System.exit( -1 );
        }
        final File mapfile = new File( args[ 1 ] );
        final BufferedReader br = new BufferedReader( new FileReader( mapfile ) );
        String line;
        final SortedMap<String, String> id_to_name_map = new TreeMap<>();
        while ( ( line = br.readLine() ) != null ) {
            final String[] s = line.split( "\\s+" );
            final String key = s[ 0 ];
            final String value = s[ 1 ];
            //
            String new_value = "N";
            final Matcher m = VIPR_PATTERN.matcher( value );
            if ( m.matches() ) {
                final String gb_accession = m.group( 3 );
                System.out.println( "gb accession: " + gb_accession );
                if ( addSpecialData1( gb_accession ) || addSpecialData12( gb_accession )
                        || addSpecialData2( gb_accession ) ) {
                    new_value = "P";
                }
            }
            //
            id_to_name_map.put( key, new_value );
        }
        br.close();
        final File infile = new File( args[ 0 ] );
        final BufferedReader br2 = new BufferedReader( new FileReader( infile ) );
        String line2;
        boolean saw_translate = false;
        boolean saw_tree = false;
        while ( ( line2 = br2.readLine() ) != null ) {
            line2 = line2.trim();
            if ( line2.equals( "translate" ) ) {
                saw_translate = true;
            }
            else if ( line2.startsWith( "tree" ) ) {
                saw_tree = true;
            }
            else if ( saw_translate && !saw_tree ) {
                final String[] s = line2.split( "\\s+" );
                final String key = s[ 0 ];
                final String value = s[ 1 ].substring( 0, s[ 1 ].length() - 1 );
                final String new_value = id_to_name_map.get( value );
                // System.out.println( key + "->" + value + "->" + new_value );
                System.out.println( key + " " + new_value );
            }
        }
        br2.close();
    }

    private static boolean addSpecialData1( final String gb_accession ) {
        if ( gb_accession.equalsIgnoreCase( "KP126912" ) || gb_accession.equalsIgnoreCase( "KP100796" )
                || gb_accession.equalsIgnoreCase( "KP744827" ) || gb_accession.equalsIgnoreCase( "KP126911" )
                || gb_accession.equalsIgnoreCase( "KP100794" ) || gb_accession.equalsIgnoreCase( "KP100792" )
                || gb_accession.equalsIgnoreCase( "KP126910" ) || gb_accession.equalsIgnoreCase( "KP100793" )
                || gb_accession.equalsIgnoreCase( "KP322752" ) || gb_accession.equalsIgnoreCase( "KY358059" )
                || gb_accession.equalsIgnoreCase( "KX685078" ) || gb_accession.equalsIgnoreCase( "KX675263" )
                || gb_accession.equalsIgnoreCase( "KX675261" ) || gb_accession.equalsIgnoreCase( "KX675262" ) ) {
            return true;
        }
        else {
            return false;
        }
    }

    private static boolean addSpecialData12( final String gb_accession ) {
        if ( gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13964_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13965_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13966_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13967_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13968_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82282_1_2347_3273.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82283_1_2352_3278.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82284_1_2347_3273.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13999_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14000_1_1_825.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14002_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14003_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14004_1_1_828.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14005_1_2354_3280.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14006_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG14007_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82280_1_2350_3276.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEJ82281_1_2348_3274.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13980_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13981_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13982_1_2354_2537.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13983_1_2354_3280.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13994_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13995_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13996_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13997_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13969_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13970_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13973_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13975_1_2342_3268.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13976_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13977_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13978_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13979_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13957_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13959_1_2346_3272.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_QEG13962_1_2346_3272.1" ) ) {
            return true;
        }
        else {
            return false;
        }
    }

    private static boolean addSpecialData2( final String gb_accession ) {
        if ( gb_accession.equalsIgnoreCase( "VIPR_ALG4_694265759_2332_3258.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_694265750_2353_3279.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_694265773_1905_2831.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_694265778_2367_3296.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_217316401_2309_3187.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_217316471_2313_3182.1" )
                || gb_accession.equalsIgnoreCase( "NP_740518.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_348549271_2240_3148.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_333980908_2253_3137.1" )
                || gb_accession.equalsIgnoreCase( "VIPR_ALG4_348549277_2043_2906.1" ) ) {
            return true;
        }
        else {
            return false;
        }
    }
}
