
package org.forester.application;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.util.ForesterUtil;

public class reform_mb_for_bats {

    private final static String  PRG_NAME     = "reform_mb_for_bats";
    private final static Pattern VIPR_PATTERN = Pattern
            .compile( "(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)\\|(.*?)" );

    public static void main( final String args[] ) throws FileNotFoundException, IOException {
        if ( args.length != 4 ) {
            System.out
                    .println( "\nWrong number of arguments, expected: infile (=mr bayes run.t outfile) map outfile first_tree\n" );
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        final File mapfile = new File( args[ 1 ] );
        final File outfile = new File( args[ 2 ] );
        
        final int first_tree = Integer.parseInt( args[ 3 ] );
        
        if ( !infile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + infile + "] does not exist" );
        }
        if ( !mapfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + mapfile + "] does not exist" );
        }
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        final BufferedWriter writer = new BufferedWriter( new OutputStreamWriter( new FileOutputStream( outfile ),
                                                                                  "utf-8" ) );
        final SortedMap<String, String> id_to_name_map = x( mapfile );
        final BufferedReader br = new BufferedReader( new FileReader( infile ) );
        String line;
        boolean saw_translate = false;
        boolean saw_tree = false;
        int tree_counter = 0;
        while ( ( line = br.readLine() ) != null ) {
            line = line.trim();
            if ( line.startsWith( "#NEXUS" ) ) {
                writer.write( line );
                writer.write( "\n" );
            }
            else if ( line.startsWith( "[ID:" ) ) {
                writer.write( "\n" );
            }
            else if ( line.startsWith( "[Param:" ) ) {
                // ignore
            }
            else if ( line.startsWith( "begin trees" ) ) {
                // ignore
            }
            else if ( line.equals( "translate" ) ) {
                saw_translate = true;
                writer.write( "begin states;" );
                writer.write( "\n" );
            }
            else if ( line.startsWith( "tree" ) ) {
                if ( !saw_tree ) {
                    saw_tree = true;
                    writer.write( "End;" );
                    writer.write( "\n" );
                    writer.write( "\n" );
                    writer.write( "begin trees;" );
                    writer.write( "\n" );
                }
                ++tree_counter;
                if ( tree_counter >= first_tree ) {
                    line = line.replaceFirst( "\\[&U\\]", "[&R]" );
                    writer.write( line );
                    writer.write( "\n" );
                }
            }
            else if ( saw_translate && !saw_tree ) {
                final String[] s = line.split( "\\s+" );
                final String key = s[ 0 ];
                final String value = s[ 1 ].substring( 0, s[ 1 ].length() - 1 );
                final String new_value = id_to_name_map.get( value );
                writer.write( key + " " + new_value );
                writer.write( "\n" );
            }
            else if ( line.startsWith( "end" ) && saw_tree ) {
                writer.write( "End;" );
                writer.write( "\n" );
            }
        }
        br.close();
        writer.flush();
        writer.close();
    }

    private static SortedMap<String, String> x( final File mapfile ) throws FileNotFoundException, IOException {
        final BufferedReader br = new BufferedReader( new FileReader( mapfile ) );
        String line;
        final SortedMap<String, String> id_to_name_map = new TreeMap<>();
        while ( ( line = br.readLine() ) != null ) {
            final String[] s = line.split( "\\s+" );
            final String key = s[ 0 ];
            final String value = s[ 1 ];
            String new_value = "N";
            final Matcher m = VIPR_PATTERN.matcher( value );
            if ( m.matches() ) {
                final String gb_accession = m.group( 3 );
                if ( addSpecialData1( gb_accession ) || addSpecialData12( gb_accession )
                        || addSpecialData2( gb_accession ) ) {
                    new_value = "P";
                }
            }
            id_to_name_map.put( key, new_value );
        }
        br.close();
        return id_to_name_map;
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
