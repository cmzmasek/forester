// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.util;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.io.Writer;
import java.math.BigDecimal;
import java.net.URL;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.archaeopteryx.Constants;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.protein.BasicProtein;
import org.forester.protein.Domain;
import org.forester.protein.Protein;
import org.forester.surfacing.SurfacingUtil;

public final class ForesterUtil {

    public final static String       FILE_SEPARATOR                   = System.getProperty( "file.separator" );
    public final static String       LINE_SEPARATOR                   = System.getProperty( "line.separator" );
    public final static String       JAVA_VENDOR                      = System.getProperty( "java.vendor" );
    public final static String       JAVA_VERSION                     = System.getProperty( "java.version" );
    public final static String       OS_ARCH                          = System.getProperty( "os.arch" );
    public final static String       OS_NAME                          = System.getProperty( "os.name" );
    public final static String       OS_VERSION                       = System.getProperty( "os.version" );
    public final static Pattern      PARANTHESESABLE_NH_CHARS_PATTERN = Pattern.compile( "[(),;\\s:\\[\\]'\"]" );
    public final static double       ZERO_DIFF                        = 1.0E-9;
    public static final BigDecimal   NULL_BD                          = new BigDecimal( 0 );
    public static final NumberFormat FORMATTER_9;
    public static final NumberFormat FORMATTER_6;
    public static final NumberFormat FORMATTER_06;
    public static final NumberFormat FORMATTER_3;
    public static final String       NCBI_PROTEIN                     = "http://www.ncbi.nlm.nih.gov/protein/";
    public static final String       NCBI_NUCCORE                     = "http://www.ncbi.nlm.nih.gov/nuccore/";
    public final static String       UNIPROT_KB                       = "http://www.uniprot.org/uniprot/";
    public static final String       NCBI_GI                          = "http://www.ncbi.nlm.nih.gov/protein/gi:";
    public static final String       PDB                              = "http://www.pdb.org/pdb/explore/explore.do?pdbId=";
    static {
        final DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator( '.' );
        // dfs.setGroupingSeparator( ( char ) 0 );
        FORMATTER_9 = new DecimalFormat( "#.#########", dfs );
        FORMATTER_6 = new DecimalFormat( "#.######", dfs );
        FORMATTER_06 = new DecimalFormat( "0.######", dfs );
        FORMATTER_3 = new DecimalFormat( "#.###", dfs );
    }

    private ForesterUtil() {
    }

    public static int calculateOverlap( final Domain domain, final List<Boolean> covered_positions ) {
        int overlap_count = 0;
        for( int i = domain.getFrom(); i <= domain.getTo(); ++i ) {
            if ( ( i < covered_positions.size() ) && ( covered_positions.get( i ) == true ) ) {
                ++overlap_count;
            }
        }
        return overlap_count;
    }

    final public static void appendSeparatorIfNotEmpty( final StringBuffer sb, final char separator ) {
        if ( sb.length() > 0 ) {
            sb.append( separator );
        }
    }

    /**
     * 
     * Example regarding engulfment: ------------0.1 ----------0.2 --0.3 =>
     * domain with 0.3 is ignored
     * 
     * -----------0.1 ----------0.2 --0.3 => domain with 0.3 is ignored
     * 
     * 
     * ------------0.1 ----------0.3 --0.2 => domains with 0.3 and 0.2 are _not_
     * ignored
     * 
     * @param max_allowed_overlap
     *            maximal allowed overlap (inclusive) to be still considered not
     *            overlapping (zero or negative value to allow any overlap)
     * @param remove_engulfed_domains
     *            to remove domains which are completely engulfed by coverage of
     *            domains with better support
     * @param protein
     * @return
     */
    public static Protein removeOverlappingDomains( final int max_allowed_overlap,
                                                    final boolean remove_engulfed_domains,
                                                    final Protein protein ) {
        final Protein pruned_protein = new BasicProtein( protein.getProteinId().getId(), protein.getSpecies()
                .getSpeciesId(), protein.getLength() );
        final List<Domain> sorted = SurfacingUtil.sortDomainsWithAscendingConfidenceValues( protein );
        final List<Boolean> covered_positions = new ArrayList<Boolean>();
        for( final Domain domain : sorted ) {
            if ( ( ( max_allowed_overlap < 0 ) || ( ForesterUtil.calculateOverlap( domain, covered_positions ) <= max_allowed_overlap ) )
                    && ( !remove_engulfed_domains || !isEngulfed( domain, covered_positions ) ) ) {
                final int covered_positions_size = covered_positions.size();
                for( int i = covered_positions_size; i < domain.getFrom(); ++i ) {
                    covered_positions.add( false );
                }
                final int new_covered_positions_size = covered_positions.size();
                for( int i = domain.getFrom(); i <= domain.getTo(); ++i ) {
                    if ( i < new_covered_positions_size ) {
                        covered_positions.set( i, true );
                    }
                    else {
                        covered_positions.add( true );
                    }
                }
                pruned_protein.addProteinDomain( domain );
            }
        }
        return pruned_protein;
    }

    /**
     * Returns true is Domain domain falls in an uninterrupted stretch of
     * covered positions.
     * 
     * @param domain
     * @param covered_positions
     * @return
     */
    public static boolean isEngulfed( final Domain domain, final List<Boolean> covered_positions ) {
        for( int i = domain.getFrom(); i <= domain.getTo(); ++i ) {
            if ( ( i >= covered_positions.size() ) || ( covered_positions.get( i ) != true ) ) {
                return false;
            }
        }
        return true;
    }

    /**
     * This calculates a color. If value is equal to min the returned color is
     * minColor, if value is equal to max the returned color is maxColor,
     * otherwise a color 'proportional' to value is returned.
     * 
     * @param value
     *            the value 
     * @param min
     *            the smallest value 
     * @param max
     *            the largest value 
     * @param minColor
     *            the color for min
     * @param maxColor
     *            the color for max
     * @return a Color
     */
    final public static Color calcColor( double value,
                                         final double min,
                                         final double max,
                                         final Color minColor,
                                         final Color maxColor ) {
        if ( value < min ) {
            value = min;
        }
        if ( value > max ) {
            value = max;
        }
        final double x = ForesterUtil.calculateColorFactor( value, max, min );
        final int red = ForesterUtil.calculateColorComponent( minColor.getRed(), maxColor.getRed(), x );
        final int green = ForesterUtil.calculateColorComponent( minColor.getGreen(), maxColor.getGreen(), x );
        final int blue = ForesterUtil.calculateColorComponent( minColor.getBlue(), maxColor.getBlue(), x );
        return new Color( red, green, blue );
    }

    /**
     * This calculates a color. If value is equal to min the returned color is
     * minColor, if value is equal to max the returned color is maxColor, if
     * value is equal to mean the returned color is meanColor, otherwise a color
     * 'proportional' to value is returned -- either between min-mean or
     * mean-max
     * 
     * @param value
     *            the value
     * @param min
     *            the smallest value
     * @param max
     *            the largest value 
     * @param mean
     *            the mean/median value 
     * @param minColor
     *            the color for min
     * @param maxColor
     *            the color for max
     * @param meanColor
     *            the color for mean
     * @return a Color
     */
    final public static Color calcColor( double value,
                                         final double min,
                                         final double max,
                                         final double mean,
                                         final Color minColor,
                                         final Color maxColor,
                                         final Color meanColor ) {
        if ( value < min ) {
            value = min;
        }
        if ( value > max ) {
            value = max;
        }
        if ( value < mean ) {
            final double x = ForesterUtil.calculateColorFactor( value, mean, min );
            final int red = ForesterUtil.calculateColorComponent( minColor.getRed(), meanColor.getRed(), x );
            final int green = ForesterUtil.calculateColorComponent( minColor.getGreen(), meanColor.getGreen(), x );
            final int blue = ForesterUtil.calculateColorComponent( minColor.getBlue(), meanColor.getBlue(), x );
            return new Color( red, green, blue );
        }
        else if ( value > mean ) {
            final double x = ForesterUtil.calculateColorFactor( value, max, mean );
            final int red = ForesterUtil.calculateColorComponent( meanColor.getRed(), maxColor.getRed(), x );
            final int green = ForesterUtil.calculateColorComponent( meanColor.getGreen(), maxColor.getGreen(), x );
            final int blue = ForesterUtil.calculateColorComponent( meanColor.getBlue(), maxColor.getBlue(), x );
            return new Color( red, green, blue );
        }
        else {
            return meanColor;
        }
    }

    final public static String collapseWhiteSpace( final String s ) {
        return s.replaceAll( "[\\s]+", " " );
    }

    final public static void collection2file( final File file, final Collection<?> data, final String separator )
            throws IOException {
        final Writer writer = new BufferedWriter( new FileWriter( file ) );
        collection2writer( writer, data, separator );
        writer.close();
    }

    final public static void collection2writer( final Writer writer, final Collection<?> data, final String separator )
            throws IOException {
        boolean first = true;
        for( final Object object : data ) {
            if ( !first ) {
                writer.write( separator );
            }
            else {
                first = false;
            }
            writer.write( object.toString() );
        }
    }

    final public static String colorToHex( final Color color ) {
        final String rgb = Integer.toHexString( color.getRGB() );
        return rgb.substring( 2, rgb.length() );
    }

    synchronized public static void copyFile( final File in, final File out ) throws IOException {
        final FileInputStream in_s = new FileInputStream( in );
        final FileOutputStream out_s = new FileOutputStream( out );
        try {
            final byte[] buf = new byte[ 1024 ];
            int i = 0;
            while ( ( i = in_s.read( buf ) ) != -1 ) {
                out_s.write( buf, 0, i );
            }
        }
        catch ( final IOException e ) {
            throw e;
        }
        finally {
            if ( in_s != null ) {
                in_s.close();
            }
            if ( out_s != null ) {
                out_s.close();
            }
        }
    }

    final public static int countChars( final String str, final char c ) {
        int count = 0;
        for( int i = 0; i < str.length(); ++i ) {
            if ( str.charAt( i ) == c ) {
                ++count;
            }
        }
        return count;
    }

    final public static BufferedWriter createBufferedWriter( final File file ) throws IOException {
        if ( file.exists() ) {
            throw new IOException( "[" + file + "] already exists" );
        }
        return new BufferedWriter( new FileWriter( file ) );
    }

    final public static BufferedWriter createBufferedWriter( final String name ) throws IOException {
        return new BufferedWriter( new FileWriter( createFileForWriting( name ) ) );
    }

    final public static EasyWriter createEasyWriter( final File file ) throws IOException {
        return new EasyWriter( createBufferedWriter( file ) );
    }

    final public static BufferedWriter createEasyWriter( final String name ) throws IOException {
        return createEasyWriter( createFileForWriting( name ) );
    }

    final public static File createFileForWriting( final String name ) throws IOException {
        final File file = new File( name );
        if ( file.exists() ) {
            throw new IOException( "[" + name + "] already exists" );
        }
        return file;
    }

    final public static void ensurePresenceOfDate( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasDate() ) {
            node.getNodeData().setDate( new org.forester.phylogeny.data.Date() );
        }
    }

    final public static void ensurePresenceOfDistribution( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasDistribution() ) {
            node.getNodeData().setDistribution( new Distribution( "" ) );
        }
    }

    public static void ensurePresenceOfSequence( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasSequence() ) {
            node.getNodeData().setSequence( new Sequence() );
        }
    }

    public static void ensurePresenceOfTaxonomy( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasTaxonomy() ) {
            node.getNodeData().setTaxonomy( new Taxonomy() );
        }
    }

    public static void fatalError( final String message ) {
        System.err.println();
        System.err.println( "error: " + message );
        System.err.println();
        System.exit( -1 );
    }

    public static void fatalError( final String prg_name, final String message ) {
        System.err.println();
        System.err.println( "[" + prg_name + "] > " + message );
        System.err.println();
        System.exit( -1 );
    }

    public static void fatalErrorIfFileNotReadable( final File file ) {
        final String error = isReadableFile( file );
        if ( !isEmpty( error ) ) {
            System.err.println();
            System.err.println( "error: " + error );
            System.err.println();
            System.exit( -1 );
        }
    }

    public static void fatalErrorIfFileNotReadable( final String prg_name, final File file ) {
        final String error = isReadableFile( file );
        if ( !isEmpty( error ) ) {
            System.err.println();
            System.err.println( "[" + prg_name + "] > " + error );
            System.err.println();
            System.exit( -1 );
        }
    }

    public static String[] file2array( final File file ) throws IOException {
        final List<String> list = file2list( file );
        final String[] ary = new String[ list.size() ];
        int i = 0;
        for( final String s : list ) {
            ary[ i++ ] = s;
        }
        return ary;
    }

    public static String[][] file22dArray( final File file ) throws IOException {
        final List<String> list = new ArrayList<String>();
        final BufferedReader in = new BufferedReader( new FileReader( file ) );
        String str;
        while ( ( str = in.readLine() ) != null ) {
            str = str.trim();
            if ( ( str.length() > 0 ) && !str.startsWith( "#" ) ) {
                list.add( str );
            }
        }
        in.close();
        final String[][] ary = new String[ list.size() ][ 2 ];
        final Pattern pa = Pattern.compile( "(\\S+)\\s+(\\S+)" );
        int i = 0;
        for( final String s : list ) {
            final Matcher m = pa.matcher( s );
            if ( m.matches() ) {
                ary[ i ][ 0 ] = m.group( 1 );
                ary[ i ][ 1 ] = m.group( 2 );
                ++i;
            }
            else {
                throw new IOException( "unexpcted format: " + s );
            }
        }
        return ary;
    }

    final public static List<String> file2list( final File file ) throws IOException {
        final List<String> list = new ArrayList<String>();
        final BufferedReader in = new BufferedReader( new FileReader( file ) );
        String str;
        while ( ( str = in.readLine() ) != null ) {
            str = str.trim();
            if ( ( str.length() > 0 ) && !str.startsWith( "#" ) ) {
                for( final String s : splitString( str ) ) {
                    list.add( s );
                }
            }
        }
        in.close();
        return list;
    }

    final public static SortedSet<String> file2set( final File file ) throws IOException {
        final SortedSet<String> set = new TreeSet<String>();
        final BufferedReader in = new BufferedReader( new FileReader( file ) );
        String str;
        while ( ( str = in.readLine() ) != null ) {
            str = str.trim();
            if ( ( str.length() > 0 ) && !str.startsWith( "#" ) ) {
                for( final String s : splitString( str ) ) {
                    set.add( s );
                }
            }
        }
        in.close();
        return set;
    }

    final public static String getCurrentDateTime() {
        final DateFormat format = new SimpleDateFormat( "yyyy/MM/dd HH:mm:ss" );
        return format.format( new Date() );
    }

    final public static String getFileSeparator() {
        return ForesterUtil.FILE_SEPARATOR;
    }

    final public static String getFirstLine( final Object source ) throws FileNotFoundException, IOException {
        BufferedReader reader = null;
        if ( source instanceof File ) {
            final File f = ( File ) source;
            if ( !f.exists() ) {
                throw new IOException( "[" + f.getAbsolutePath() + "] does not exist" );
            }
            else if ( !f.isFile() ) {
                throw new IOException( "[" + f.getAbsolutePath() + "] is not a file" );
            }
            else if ( !f.canRead() ) {
                throw new IOException( "[" + f.getAbsolutePath() + "] is not a readable" );
            }
            reader = new BufferedReader( new FileReader( f ) );
        }
        else if ( source instanceof InputStream ) {
            reader = new BufferedReader( new InputStreamReader( ( InputStream ) source ) );
        }
        else if ( source instanceof String ) {
            reader = new BufferedReader( new StringReader( ( String ) source ) );
        }
        else if ( source instanceof StringBuffer ) {
            reader = new BufferedReader( new StringReader( source.toString() ) );
        }
        else if ( source instanceof URL ) {
            reader = new BufferedReader( new InputStreamReader( ( ( URL ) source ).openStream() ) );
        }
        else {
            throw new IllegalArgumentException( "dont know how to read [" + source.getClass() + "]" );
        }
        String line;
        while ( ( line = reader.readLine() ) != null ) {
            line = line.trim();
            if ( !ForesterUtil.isEmpty( line ) ) {
                if ( reader != null ) {
                    reader.close();
                }
                return line;
            }
        }
        if ( reader != null ) {
            reader.close();
        }
        return line;
    }

    final public static String getForesterLibraryInformation() {
        return "forester " + ForesterConstants.FORESTER_VERSION + " (" + ForesterConstants.FORESTER_DATE + ")";
    }

    final public static String getLineSeparator() {
        return ForesterUtil.LINE_SEPARATOR;
    }

    final public static void increaseCountingMap( final Map<String, Integer> counting_map, final String item_name ) {
        if ( !counting_map.containsKey( item_name ) ) {
            counting_map.put( item_name, 1 );
        }
        else {
            counting_map.put( item_name, counting_map.get( item_name ) + 1 );
        }
    }

    final public static boolean isContainsParanthesesableNhCharacter( final String nh ) {
        return PARANTHESESABLE_NH_CHARS_PATTERN.matcher( nh ).find();
    }

    final public static boolean isEmpty( final List<?> l ) {
        if ( ( l == null ) || l.isEmpty() ) {
            return true;
        }
        for( final Object o : l ) {
            if ( o != null ) {
                return false;
            }
        }
        return true;
    }

    final public static boolean isEmpty( final Set<?> s ) {
        if ( ( s == null ) || s.isEmpty() ) {
            return true;
        }
        for( final Object o : s ) {
            if ( o != null ) {
                return false;
            }
        }
        return true;
    }

    final public static boolean isEmpty( final String s ) {
        return ( ( s == null ) || ( s.length() < 1 ) );
    }

    final public static boolean isEqual( final double a, final double b ) {
        return ( ( Math.abs( a - b ) ) < ZERO_DIFF );
    }

    final public static boolean isEven( final int n ) {
        return ( n % 2 ) == 0;
    }

    /**
     * This determines whether String[] a and String[] b have at least one
     * String in common (intersect). Returns false if at least one String[] is
     * null or empty.
     * 
     * @param a
     *            a String[] b a String[]
     * @return true if both a and b or not empty or null and contain at least
     *         one element in common false otherwise
     */
    final public static boolean isIntersecting( final String[] a, final String[] b ) {
        if ( ( a == null ) || ( b == null ) ) {
            return false;
        }
        if ( ( a.length < 1 ) || ( b.length < 1 ) ) {
            return false;
        }
        for( final String ai : a ) {
            for( final String element : b ) {
                if ( ( ai != null ) && ( element != null ) && ai.equals( element ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    final public static double isLargerOrEqualToZero( final double d ) {
        if ( d > 0.0 ) {
            return d;
        }
        else {
            return 0.0;
        }
    }

    final public static boolean isNull( final BigDecimal s ) {
        return ( ( s == null ) || ( s.compareTo( NULL_BD ) == 0 ) );
    }

    final public static String isReadableFile( final File f ) {
        if ( !f.exists() ) {
            return "file [" + f + "] does not exist";
        }
        if ( f.isDirectory() ) {
            return "[" + f + "] is a directory";
        }
        if ( !f.isFile() ) {
            return "[" + f + "] is not a file";
        }
        if ( !f.canRead() ) {
            return "file [" + f + "] is not readable";
        }
        if ( f.length() < 1 ) {
            return "file [" + f + "] is empty";
        }
        return null;
    }

    final public static String isReadableFile( final String s ) {
        return isReadableFile( new File( s ) );
    }

    public final static boolean isWindows() {
        try {
            return OS_NAME.toLowerCase().indexOf( "win" ) > -1;
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "minor error: " + e );
            return false;
        }
    }

    public final static boolean isMac() {
        try {
            return OS_NAME.toLowerCase().startsWith( "mac" );
        }
        catch ( final Exception e ) {
            ForesterUtil.printWarningMessage( Constants.PRG_NAME, "minor error: " + e );
            return false;
        }
    }

    final public static String isWritableFile( final File f ) {
        if ( f.isDirectory() ) {
            return "[" + f + "] is a directory";
        }
        if ( f.exists() ) {
            return "[" + f + "] already exists";
        }
        return null;
    }

    /**
     * Helper for method "stringToColor".
     * <p>
     * (Last modified: 12/20/03)
     */
    final public static int limitRangeForColor( int i ) {
        if ( i > 255 ) {
            i = 255;
        }
        else if ( i < 0 ) {
            i = 0;
        }
        return i;
    }

    final public static SortedMap<Object, Integer> listToSortedCountsMap( final List<?> list ) {
        final SortedMap<Object, Integer> map = new TreeMap<Object, Integer>();
        for( final Object key : list ) {
            if ( !map.containsKey( key ) ) {
                map.put( key, 1 );
            }
            else {
                map.put( key, map.get( key ) + 1 );
            }
        }
        return map;
    }

    final public static void map2file( final File file,
                                       final Map<?, ?> data,
                                       final String entry_separator,
                                       final String data_separator ) throws IOException {
        final Writer writer = new BufferedWriter( new FileWriter( file ) );
        map2writer( writer, data, entry_separator, data_separator );
        writer.close();
    }

    final public static void map2writer( final Writer writer,
                                         final Map<?, ?> data,
                                         final String entry_separator,
                                         final String data_separator ) throws IOException {
        boolean first = true;
        for( final Entry<?, ?> entry : data.entrySet() ) {
            if ( !first ) {
                writer.write( data_separator );
            }
            else {
                first = false;
            }
            writer.write( entry.getKey().toString() );
            writer.write( entry_separator );
            writer.write( entry.getValue().toString() );
        }
    }

    final public static StringBuffer mapToStringBuffer( final Map<Object, Object> map, final String key_value_separator ) {
        final StringBuffer sb = new StringBuffer();
        for( final Object key : map.keySet() ) {
            sb.append( key.toString() );
            sb.append( key_value_separator );
            sb.append( map.get( key ).toString() );
            sb.append( ForesterUtil.getLineSeparator() );
        }
        return sb;
    }

    final public static String normalizeString( final String s,
                                                final int length,
                                                final boolean left_pad,
                                                final char pad_char ) {
        if ( s.length() > length ) {
            return s.substring( 0, length );
        }
        else {
            final StringBuffer pad = new StringBuffer( length - s.length() );
            for( int i = 0; i < ( length - s.length() ); ++i ) {
                pad.append( pad_char );
            }
            if ( left_pad ) {
                return pad + s;
            }
            else {
                return s + pad;
            }
        }
    }

    final public static BufferedReader obtainReader( final Object source ) throws IOException, FileNotFoundException {
        BufferedReader reader = null;
        if ( source instanceof File ) {
            final File f = ( File ) source;
            if ( !f.exists() ) {
                throw new IOException( "\"" + f.getAbsolutePath() + "\" does not exist" );
            }
            else if ( !f.isFile() ) {
                throw new IOException( "\"" + f.getAbsolutePath() + "\" is not a file" );
            }
            else if ( !f.canRead() ) {
                throw new IOException( "\"" + f.getAbsolutePath() + "\" is not a readable" );
            }
            reader = new BufferedReader( new FileReader( f ) );
        }
        else if ( source instanceof InputStream ) {
            reader = new BufferedReader( new InputStreamReader( ( InputStream ) source ) );
        }
        else if ( source instanceof String ) {
            reader = new BufferedReader( new StringReader( ( String ) source ) );
        }
        else if ( source instanceof StringBuffer ) {
            reader = new BufferedReader( new StringReader( source.toString() ) );
        }
        else {
            throw new IllegalArgumentException( "attempt to parse object of type [" + source.getClass()
                    + "] (can only parse objects of type File, InputStream, String, or StringBuffer)" );
        }
        return reader;
    }

    final public static StringBuffer pad( final double number, final int size, final char pad, final boolean left_pad ) {
        return pad( new StringBuffer( number + "" ), size, pad, left_pad );
    }

    final public static StringBuffer pad( final String string, final int size, final char pad, final boolean left_pad ) {
        return pad( new StringBuffer( string ), size, pad, left_pad );
    }

    final public static StringBuffer pad( final StringBuffer string,
                                          final int size,
                                          final char pad,
                                          final boolean left_pad ) {
        final StringBuffer padding = new StringBuffer();
        final int s = size - string.length();
        if ( s < 1 ) {
            return new StringBuffer( string.substring( 0, size ) );
        }
        for( int i = 0; i < s; ++i ) {
            padding.append( pad );
        }
        if ( left_pad ) {
            return padding.append( string );
        }
        else {
            return string.append( padding );
        }
    }

    final public static double parseDouble( final String str ) throws ParseException {
        if ( ForesterUtil.isEmpty( str ) ) {
            return 0.0;
        }
        return Double.parseDouble( str );
    }

    final public static int parseInt( final String str ) throws ParseException {
        if ( ForesterUtil.isEmpty( str ) ) {
            return 0;
        }
        return Integer.parseInt( str );
    }

    final public static void printArray( final Object[] a ) {
        for( int i = 0; i < a.length; ++i ) {
            System.out.println( "[" + i + "]=" + a[ i ] );
        }
    }

    final public static void printCountingMap( final Map<String, Integer> counting_map ) {
        for( final String key : counting_map.keySet() ) {
            System.out.println( key + ": " + counting_map.get( key ) );
        }
    }

    final public static void printErrorMessage( final String prg_name, final String message ) {
        System.err.println( "[" + prg_name + "] > error: " + message );
    }

    final public static void printProgramInformation( final String prg_name, final String prg_version, final String date ) {
        final int l = prg_name.length() + prg_version.length() + date.length() + 4;
        System.out.println();
        System.out.println( prg_name + " " + prg_version + " (" + date + ")" );
        for( int i = 0; i < l; ++i ) {
            System.out.print( "_" );
        }
        System.out.println();
    }

    final public static void printProgramInformation( final String prg_name,
                                                      final String prg_version,
                                                      final String date,
                                                      final String email,
                                                      final String www ) {
        printProgramInformation( prg_name, null, prg_version, date, email, www, null );
    }

    final public static void printProgramInformation( final String prg_name,
                                                      final String desc,
                                                      final String prg_version,
                                                      final String date,
                                                      final String email,
                                                      final String www,
                                                      final String based_on ) {
        String my_prg_name = new String( prg_name );
        if ( !ForesterUtil.isEmpty( desc ) ) {
            my_prg_name += ( " - " + desc );
        }
        final int l = my_prg_name.length() + prg_version.length() + date.length() + 4;
        System.out.println();
        System.out.println( my_prg_name + " " + prg_version + " (" + date + ")" );
        for( int i = 0; i < l; ++i ) {
            System.out.print( "_" );
        }
        System.out.println();
        System.out.println();
        System.out.println( "WWW     : " + www );
        System.out.println( "Contact : " + email );
        if ( !ForesterUtil.isEmpty( based_on ) ) {
            System.out.println( "Based on: " + based_on );
        }
        if ( !ForesterUtil.isEmpty( ForesterUtil.JAVA_VERSION ) && !ForesterUtil.isEmpty( ForesterUtil.JAVA_VENDOR ) ) {
            System.out.println();
            System.out.println( "[running on Java " + ForesterUtil.JAVA_VERSION + " " + ForesterUtil.JAVA_VENDOR + "]" );
        }
        System.out.println();
    }

    final public static void printWarningMessage( final String prg_name, final String message ) {
        System.out.println( "[" + prg_name + "] > warning: " + message );
    }

    final public static void programMessage( final String prg_name, final String message ) {
        System.out.println( "[" + prg_name + "] > " + message );
    }

    final public static String removeSuffix( final String file_name ) {
        final int i = file_name.lastIndexOf( '.' );
        if ( i > 1 ) {
            return file_name.substring( 0, i );
        }
        return file_name;
    }

    /**
     * Removes all white space from String s.
     * 
     * @return String s with white space removed
     */
    final public static String removeWhiteSpace( String s ) {
        int i;
        for( i = 0; i <= ( s.length() - 1 ); i++ ) {
            if ( ( s.charAt( i ) == ' ' ) || ( s.charAt( i ) == '\t' ) || ( s.charAt( i ) == '\n' )
                    || ( s.charAt( i ) == '\r' ) ) {
                s = s.substring( 0, i ) + s.substring( i + 1 );
                i--;
            }
        }
        return s;
    }

    final public static String replaceIllegalNhxCharacters( final String nhx ) {
        if ( nhx == null ) {
            return "";
        }
        return nhx.trim().replaceAll( "[\\[\\]']+", "_" );
    }

    final public static double round( final double value, final int decimal_place ) {
        BigDecimal bd = new BigDecimal( value );
        bd = bd.setScale( decimal_place, BigDecimal.ROUND_HALF_UP );
        return bd.doubleValue();
    }

    /**
     * Rounds d to an int.
     */
    final public static int roundToInt( final double d ) {
        return ( int ) ( d + 0.5 );
    }

    final public static int roundToInt( final float f ) {
        return ( int ) ( f + 0.5f );
    }

    final public static short roundToShort( final double d ) {
        return ( short ) ( d + 0.5 );
    }

    final public static String sanitizeString( final String s ) {
        if ( s == null ) {
            return "";
        }
        else {
            return s.trim();
        }
    }

    public static boolean seqIsLikelyToBeAa( final String s ) {
        final String seq = s.toLowerCase();
        if ( ( seq.indexOf( 'r' ) > -1 ) || ( seq.indexOf( 'd' ) > -1 ) || ( seq.indexOf( 'e' ) > -1 )
                || ( seq.indexOf( 'q' ) > -1 ) || ( seq.indexOf( 'h' ) > -1 ) || ( seq.indexOf( 'k' ) > -1 )
                || ( seq.indexOf( 'w' ) > -1 ) || ( seq.indexOf( 's' ) > -1 ) || ( seq.indexOf( 'm' ) > -1 )
                || ( seq.indexOf( 'p' ) > -1 ) || ( seq.indexOf( 'v' ) > -1 ) ) {
            return true;
        }
        return false;
    }

    final public static String stringArrayToString( final String[] a ) {
        return stringArrayToString( a, ", " );
    }

    final public static String stringArrayToString( final String[] a, final String separator ) {
        final StringBuilder sb = new StringBuilder();
        if ( ( a != null ) && ( a.length > 0 ) ) {
            for( int i = 0; i < ( a.length - 1 ); ++i ) {
                sb.append( a[ i ] + separator );
            }
            sb.append( a[ a.length - 1 ] );
        }
        return sb.toString();
    }

    final public static String[] stringListToArray( final List<String> list ) {
        if ( list != null ) {
            final String[] str = new String[ list.size() ];
            int i = 0;
            for( final String l : list ) {
                str[ i++ ] = l;
            }
            return str;
        }
        return null;
    }

    final public static String stringListToString( final List<String> l, final String separator ) {
        final StringBuilder sb = new StringBuilder();
        if ( ( l != null ) && ( l.size() > 0 ) ) {
            for( int i = 0; i < ( l.size() - 1 ); ++i ) {
                sb.append( l.get( i ) + separator );
            }
            sb.append( l.get( l.size() - 1 ) );
        }
        return sb.toString();
    }

    final public static String[] stringSetToArray( final Set<String> strings ) {
        final String[] str_array = new String[ strings.size() ];
        int i = 0;
        for( final String e : strings ) {
            str_array[ i++ ] = e;
        }
        return str_array;
    }

    final public static void unexpectedFatalError( final Exception e ) {
        System.err.println();
        System.err.println( "unexpected exception: should not have occured! Please contact program author(s)." );
        e.printStackTrace( System.err );
        System.err.println();
        System.exit( -1 );
    }

    final public static void unexpectedFatalError( final Error e ) {
        System.err.println();
        System.err.println( "unexpected error: should not have occured! Please contact program author(s)." );
        e.printStackTrace( System.err );
        System.err.println();
        System.exit( -1 );
    }

    final public static void unexpectedFatalError( final String message ) {
        System.err.println();
        System.err.println( "unexpected error: should not have occured! Please contact program author(s)." );
        System.err.println( message );
        System.err.println();
        System.exit( -1 );
    }

    final public static void unexpectedFatalError( final String prg_name, final Exception e ) {
        System.err.println();
        System.err.println( "[" + prg_name
                + "] > unexpected error; should not have occured! Please contact program author(s)." );
        e.printStackTrace( System.err );
        System.err.println();
        System.exit( -1 );
    }

    final public static void unexpectedFatalError( final String prg_name, final String message ) {
        System.err.println();
        System.err.println( "[" + prg_name
                + "] > unexpected error: should not have occured! Please contact program author(s)." );
        System.err.println( message );
        System.err.println();
        System.exit( -1 );
    }

    final public static void unexpectedFatalError( final String prg_name, final String message, final Exception e ) {
        System.err.println();
        System.err.println( "[" + prg_name
                + "] > unexpected error: should not have occured! Please contact program author(s)." );
        System.err.println( message );
        e.printStackTrace( System.err );
        System.err.println();
        System.exit( -1 );
    }

    public final static void updateProgress( final double progress_percentage ) {
        final int width = 50;
        System.out.print( "\r[" );
        int i = 0;
        for( ; i <= ForesterUtil.roundToInt( progress_percentage * width ); i++ ) {
            System.out.print( "." );
        }
        for( ; i < width; i++ ) {
            System.out.print( " " );
        }
        System.out.print( "]" );
    }

    public final static void updateProgress( final int i, final DecimalFormat f ) {
        System.out.print( "\r[" + f.format( i ) + "]" );
    }

    public final static String wordWrap( final String str, final int width ) {
        final StringBuilder sb = new StringBuilder( str );
        int start = 0;
        int ls = -1;
        int i = 0;
        while ( i < sb.length() ) {
            if ( sb.charAt( i ) == ' ' ) {
                ls = i;
            }
            if ( sb.charAt( i ) == '\n' ) {
                ls = -1;
                start = i + 1;
            }
            if ( i > ( ( start + width ) - 1 ) ) {
                if ( ls != -1 ) {
                    sb.setCharAt( ls, '\n' );
                    start = ls + 1;
                    ls = -1;
                }
                else {
                    sb.insert( i, '\n' );
                    start = i + 1;
                }
            }
            i++;
        }
        return sb.toString();
    }

    /**
     * Helper method for calcColor methods.
     * 
     * @param smallercolor_component_x
     *            color component the smaller color
     * @param largercolor_component_x
     *            color component the larger color
     * @param x
     *            factor
     * @return an int representing a color component
     */
    final private static int calculateColorComponent( final double smallercolor_component_x,
                                                      final double largercolor_component_x,
                                                      final double x ) {
        return ( int ) ( smallercolor_component_x + ( ( x * ( largercolor_component_x - smallercolor_component_x ) ) / 255.0 ) );
    }

    /**
     * Helper method for calcColor methods.
     * 
     * 
     * @param value
     *            the value
     * @param larger
     *            the largest value
     * @param smaller
     *            the smallest value
     * @return a normalized value between larger and smaller
     */
    final private static double calculateColorFactor( final double value, final double larger, final double smaller ) {
        return ( 255.0 * ( value - smaller ) ) / ( larger - smaller );
    }

    final private static String[] splitString( final String str ) {
        final String regex = "[\\s;,]+";
        return str.split( regex );
    }

    public final static void outOfMemoryError( final OutOfMemoryError e ) {
        System.err.println();
        System.err.println( "Java memory allocation might be too small, try \"-Xmx2048m\" java command line option" );
        System.err.println();
        e.printStackTrace( System.err );
        System.err.println();
        System.exit( -1 );
    }

    public final static Color obtainColorDependingOnTaxonomyGroup( final String tax_group ) {
        if ( !ForesterUtil.isEmpty( tax_group ) ) {
            if ( tax_group.equals( TaxonomyGroups.DEUTEROSTOMIA ) ) {
                return TaxonomyColors.DEUTEROSTOMIA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.PROTOSTOMIA ) ) {
                return TaxonomyColors.PROTOSTOMIA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.CNIDARIA ) ) {
                return TaxonomyColors.CNIDARIA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.PLACOZOA ) ) {
                return TaxonomyColors.PLACOZOA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.CTENOPHORA ) ) {
                return TaxonomyColors.CTENOPHORA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.PORIFERA ) ) {
                return TaxonomyColors.PORIFERA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.CHOANOFLAGELLIDA ) ) {
                return TaxonomyColors.CHOANOFLAGELLIDA;
            }
            else if ( tax_group.equals( TaxonomyGroups.ICHTHYOPHONIDA_FILASTEREA ) ) {
                return TaxonomyColors.ICHTHYOSPOREA_AND_FILASTEREA;
            }
            else if ( tax_group.equals( TaxonomyGroups.DIKARYA ) ) {
                return TaxonomyColors.DIKARYA_COLOR;
            }
            else if ( tax_group.equalsIgnoreCase( TaxonomyGroups.FUNGI )
                    || tax_group.equalsIgnoreCase( TaxonomyGroups.OTHER_FUNGI ) ) {
                return TaxonomyColors.OTHER_FUNGI_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.NUCLEARIIDAE_AND_FONTICULA_GROUP ) ) {
                return TaxonomyColors.NUCLEARIIDAE_AND_FONTICULA_GROUP_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.AMOEBOZOA ) ) {
                return TaxonomyColors.AMOEBOZOA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.EMBRYOPHYTA ) ) {
                return TaxonomyColors.EMBRYOPHYTA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.CHLOROPHYTA ) ) {
                return TaxonomyColors.CHLOROPHYTA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.RHODOPHYTA ) ) {
                return TaxonomyColors.RHODOPHYTA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.HACROBIA ) ) {
                return TaxonomyColors.HACROBIA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.GLAUCOCYSTOPHYCEAE ) ) {
                return TaxonomyColors.GLAUCOPHYTA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.STRAMENOPILES ) ) {
                return TaxonomyColors.STRAMENOPILES_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.ALVEOLATA ) ) {
                return TaxonomyColors.ALVEOLATA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.RHIZARIA ) ) {
                return TaxonomyColors.RHIZARIA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.EXCAVATA ) ) {
                return TaxonomyColors.EXCAVATA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.APUSOZOA ) ) {
                return TaxonomyColors.APUSOZOA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.ARCHAEA ) ) {
                return TaxonomyColors.ARCHAEA_COLOR;
            }
            else if ( tax_group.equals( TaxonomyGroups.BACTERIA ) ) {
                return TaxonomyColors.BACTERIA_COLOR;
            }
        }
        return null;
    }

    public final static String obtainNormalizedTaxonomyGroup( final String tax ) {
        if ( tax.equalsIgnoreCase( TaxonomyGroups.DEUTEROSTOMIA ) ) {
            return TaxonomyGroups.DEUTEROSTOMIA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.PROTOSTOMIA ) ) {
            return TaxonomyGroups.PROTOSTOMIA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.CNIDARIA ) ) {
            return TaxonomyGroups.CNIDARIA;
        }
        else if ( tax.toLowerCase().startsWith( "trichoplax" ) || tax.equalsIgnoreCase( TaxonomyGroups.PLACOZOA ) ) {
            return TaxonomyGroups.PLACOZOA;
        }
        else if ( tax.toLowerCase().startsWith( "mnemiopsis" ) || tax.equalsIgnoreCase( TaxonomyGroups.CTENOPHORA ) ) {
            return TaxonomyGroups.CTENOPHORA;
        }
        else if ( tax.toLowerCase().startsWith( "amphimedon" ) || tax.equalsIgnoreCase( TaxonomyGroups.PORIFERA ) ) {
            return TaxonomyGroups.PORIFERA;
        }
        else if ( tax.equalsIgnoreCase( "codonosigidae" ) || tax.equalsIgnoreCase( TaxonomyGroups.CHOANOFLAGELLIDA ) ) {
            return TaxonomyGroups.CHOANOFLAGELLIDA;
        }
        else if ( tax.toLowerCase().startsWith( TaxonomyGroups.ICHTHYOPHONIDA_FILASTEREA )
                || tax.toLowerCase().startsWith( "ichthyophonida and filasterea" )
                || tax.toLowerCase().startsWith( "ichthyosporea & filasterea" )
                || tax.toLowerCase().startsWith( "ichthyosporea and filasterea" ) ) {
            return TaxonomyGroups.ICHTHYOPHONIDA_FILASTEREA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.DIKARYA ) ) {
            return TaxonomyGroups.DIKARYA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.FUNGI ) || tax.equalsIgnoreCase( TaxonomyGroups.OTHER_FUNGI ) ) {
            return TaxonomyGroups.OTHER_FUNGI;
        }
        else if ( tax.toLowerCase().startsWith( "nucleariidae and fonticula" ) ) {
            return TaxonomyGroups.NUCLEARIIDAE_AND_FONTICULA_GROUP;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.AMOEBOZOA ) ) {
            return TaxonomyGroups.AMOEBOZOA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.EMBRYOPHYTA ) ) {
            return TaxonomyGroups.EMBRYOPHYTA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.CHLOROPHYTA ) ) {
            return TaxonomyGroups.CHLOROPHYTA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.RHODOPHYTA ) ) {
            return TaxonomyGroups.RHODOPHYTA;
        }
        else if ( tax.toLowerCase().startsWith( TaxonomyGroups.HACROBIA ) ) {
            return TaxonomyGroups.HACROBIA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.GLAUCOCYSTOPHYCEAE ) || tax.equalsIgnoreCase( "glaucophyta" ) ) {
            return TaxonomyGroups.GLAUCOCYSTOPHYCEAE;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.STRAMENOPILES ) ) {
            return TaxonomyGroups.STRAMENOPILES;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.ALVEOLATA ) ) {
            return TaxonomyGroups.ALVEOLATA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.RHIZARIA ) ) {
            return TaxonomyGroups.RHIZARIA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.EXCAVATA ) ) {
            return TaxonomyGroups.EXCAVATA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.APUSOZOA ) ) {
            return TaxonomyGroups.APUSOZOA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.ARCHAEA ) ) {
            return TaxonomyGroups.ARCHAEA;
        }
        else if ( tax.equalsIgnoreCase( TaxonomyGroups.BACTERIA ) ) {
            return TaxonomyGroups.BACTERIA;
        }
        return null;
    }
}
