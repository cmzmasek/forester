// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

package org.forester.ws.seqdb;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.go.GoTerm;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Accession.Source;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

public final class SequenceDbWsTools {

    public final static String   BASE_UNIPROT_URL        = "http://www.uniprot.org/";
    public final static int      DEFAULT_LINES_TO_RETURN = 4000;
    //public final static String   EMBL_DBS_EMBL           = "embl";
    public final static String   EMBL_DBS_REFSEQ_N       = "refseqn";
    public final static String   EMBL_DBS_REFSEQ_P       = "refseqp";
    public final static String   EMBL_GENBANK            = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=GENBANK&style=raw&id=";
    public final static String   EMBL_REFSEQ             = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=REFSEQ&style=raw&id=";
    public final static String   EMBL_EMBL               = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=EMBL&style=raw&id=";
    private final static boolean DEBUG                   = true;
    private final static String  URL_ENC                 = "UTF-8";

    public static List<UniProtTaxonomy> getTaxonomiesFromCommonNameStrict( final String cn,
                                                                           final int max_taxonomies_return )
            throws IOException {
        final List<UniProtTaxonomy> taxonomies = getTaxonomiesFromCommonName( cn, max_taxonomies_return );
        if ( ( taxonomies != null ) && ( taxonomies.size() > 0 ) ) {
            final List<UniProtTaxonomy> filtered_taxonomies = new ArrayList<UniProtTaxonomy>();
            for( final UniProtTaxonomy taxonomy : taxonomies ) {
                if ( taxonomy.getCommonName().equalsIgnoreCase( cn ) ) {
                    filtered_taxonomies.add( taxonomy );
                }
            }
            return filtered_taxonomies;
        }
        return null;
    }

    public static List<UniProtTaxonomy> getTaxonomiesFromId( final String id, final int max_taxonomies_return )
            throws IOException {
        final List<String> result = getTaxonomyStringFromId( id, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    /**
     * Does not return "sub-types".
     * For example, for "Mus musculus" only returns "Mus musculus"
     * and not "Mus musculus", "Mus musculus bactrianus", ...
     * 
     */
    public static List<UniProtTaxonomy> getTaxonomiesFromScientificNameStrict( final String sn,
                                                                               final int max_taxonomies_return )
            throws IOException {
        final List<UniProtTaxonomy> taxonomies = getTaxonomiesFromScientificName( sn, max_taxonomies_return );
        if ( ( taxonomies != null ) && ( taxonomies.size() > 0 ) ) {
            final List<UniProtTaxonomy> filtered_taxonomies = new ArrayList<UniProtTaxonomy>();
            for( final UniProtTaxonomy taxonomy : taxonomies ) {
                if ( taxonomy.getScientificName().equalsIgnoreCase( sn ) ) {
                    filtered_taxonomies.add( taxonomy );
                }
            }
            return filtered_taxonomies;
        }
        return null;
    }

    public static List<UniProtTaxonomy> getTaxonomiesFromTaxonomyCode( final String code,
                                                                       final int max_taxonomies_return )
            throws IOException {
        final String my_code = new String( code );
        final List<String> result = getTaxonomyStringFromTaxonomyCode( my_code, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    public static SequenceDatabaseEntry obtainEmblEntry( final Accession acc ) throws IOException {
        return obtainEmblEntry( acc, DEFAULT_LINES_TO_RETURN );
    }

    public static SequenceDatabaseEntry obtainEmblEntry( final Accession acc, final int max_lines_to_return )
            throws IOException {
        final List<String> lines = queryEmblDb( acc, max_lines_to_return );
        return EbiDbEntry.createInstanceFromPlainTextForRefSeq( lines );
    }

    public static SequenceDatabaseEntry obtainEntry( final String acc_str ) throws IOException {
        if ( ForesterUtil.isEmpty( acc_str ) ) {
            throw new IllegalArgumentException( "cannot not extract sequence db accessor from null or empty string" );
        }
        final Accession acc = SequenceAccessionTools.parseAccessorFromString( acc_str );
        if ( acc == null ) {
            throw new IllegalArgumentException( "could not extract acceptable sequence db accessor from \"" + acc_str
                    + "\"" );
        }
        if ( acc.getSource().equals( Source.REFSEQ.toString() ) || acc.getSource().equals( Source.EMBL.toString() )
                || acc.getSource().equals( Source.NCBI.toString() ) ) {
            return obtainEmblEntry( acc, DEFAULT_LINES_TO_RETURN );
        }
        else if ( acc.getSource().equals( Source.UNIPROT.toString() ) ) {
            return obtainUniProtEntry( acc.getValue(), DEFAULT_LINES_TO_RETURN );
        }
        else {
            throw new IllegalArgumentException( "don't know how to handle request for source \"" + acc.getSource()
                    + "\"" );
        }
    }

    public static SequenceDatabaseEntry obtainRefSeqEntryFromEmbl( final Accession acc ) throws IOException {
        return obtainRefSeqEntryFromEmbl( acc, DEFAULT_LINES_TO_RETURN );
    }

    public static SequenceDatabaseEntry obtainRefSeqEntryFromEmbl( final Accession acc, final int max_lines_to_return )
            throws IOException {
        final List<String> lines = queryEmblDbForRefSeqEntry( acc, max_lines_to_return );
        return EbiDbEntry.createInstanceFromPlainTextForRefSeq( lines );
    }

    public final static Accession obtainSeqAccession( final PhylogenyNode node ) {
        Accession acc = SequenceAccessionTools.obtainFromSeqAccession( node );
        if ( !isAccessionAcceptable( acc ) ) {
            acc = SequenceAccessionTools.obtainAccessorFromDataFields( node );
        }
        return acc;
    }

    public final static void obtainSeqInformation( final boolean allow_to_set_taxonomic_data,
                                                   final int lines_to_return,
                                                   final SortedSet<String> not_found,
                                                   final PhylogenyNode node ) throws IOException {
        final Accession acc = obtainSeqAccession( node );
        if ( !isAccessionAcceptable( acc ) ) {
            if ( node.isExternal() || !node.isEmpty() ) {
                not_found.add( node.toString() );
            }
        }
        else {
            addDataFromDbToNode( allow_to_set_taxonomic_data, lines_to_return, not_found, node, acc );
        }
    }

    public final static void obtainSeqInformation( final boolean allow_to_set_taxonomic_data,
                                                   final SortedSet<String> not_found,
                                                   final PhylogenyNode node ) throws IOException {
        obtainSeqInformation( allow_to_set_taxonomic_data, DEFAULT_LINES_TO_RETURN, not_found, node );
    }

    public final static SortedSet<String> obtainSeqInformation( final Phylogeny phy,
                                                                final boolean ext_nodes_only,
                                                                final boolean allow_to_set_taxonomic_data,
                                                                final int lines_to_return ) throws IOException {
        final SortedSet<String> not_found = new TreeSet<String>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() || !ext_nodes_only ) {
                obtainSeqInformation( allow_to_set_taxonomic_data, lines_to_return, not_found, node );
            }
        }
        return not_found;
    }

    public final static void obtainSeqInformation( final PhylogenyNode node ) throws IOException {
        obtainSeqInformation( true, DEFAULT_LINES_TO_RETURN, new TreeSet<String>(), node );
    }

    public static SequenceDatabaseEntry obtainUniProtEntry( final String query ) throws IOException {
        return obtainUniProtEntry( query, DEFAULT_LINES_TO_RETURN );
    }

    public static SequenceDatabaseEntry obtainUniProtEntry( final String query, final int max_lines_to_return )
            throws IOException {
        final List<String> lines = queryUniprot( "uniprot/" + query + ".txt", max_lines_to_return );
        return UniProtEntry.createInstanceFromPlainText( lines );
    }

    public static List<String> queryDb( final String query, int max_lines_to_return, final String base_url )
            throws IOException {
        if ( ForesterUtil.isEmpty( query ) ) {
            throw new IllegalArgumentException( "illegal attempt to use empty query " );
        }
        if ( max_lines_to_return < 1 ) {
            max_lines_to_return = 1;
        }
        final URL url = new URL( base_url + query );
        if ( DEBUG ) {
            System.out.println( "url: " + url.toString() );
        }
        final URLConnection urlc = url.openConnection();
        final BufferedReader in = new BufferedReader( new InputStreamReader( urlc.getInputStream() ) );
        String line;
        final List<String> result = new ArrayList<String>();
        while ( ( line = in.readLine() ) != null ) {
            if ( DEBUG ) {
                System.out.println( line );
            }
            result.add( line );
            if ( result.size() > max_lines_to_return ) {
                break;
            }
        }
        in.close();
        try {
            // To prevent accessing online dbs in too quick succession. 
            Thread.sleep( 20 );
        }
        catch ( final InterruptedException e ) {
            e.printStackTrace();
        }
        return result;
    }

    public static List<String> queryEmblDb( final Accession acc, final int max_lines_to_return ) throws IOException {
        final StringBuilder url_sb = new StringBuilder();
        //  url_sb.append( BASE_EMBL_DB_URL );
        System.out.println( "source: " + acc.getSource() );
        if ( acc.getSource().equals( Source.NCBI.toString() ) ) {
            url_sb.append( EMBL_GENBANK );
            //url_sb.append( '/' );
        }
        else if ( acc.getSource().equals( Source.REFSEQ.toString() ) ) {
            url_sb.append( EMBL_REFSEQ );
        }
        else if ( acc.getSource().equals( Source.EMBL.toString() ) ) {
            url_sb.append( EMBL_EMBL );
        }
        else {
            throw new IllegalArgumentException( "unable to handle source: " + acc.getSource() );
        }
        return queryDb( acc.getValue(), max_lines_to_return, url_sb.toString() );
    }

    public static List<String> queryEmblDbForRefSeqEntry( final Accession id, final int max_lines_to_return )
            throws IOException {
        final StringBuilder url_sb = new StringBuilder();
        url_sb.append( EMBL_REFSEQ );
        return queryDb( id.getValue(), max_lines_to_return, url_sb.toString() );
    }

    public static List<String> queryUniprot( final String query, final int max_lines_to_return ) throws IOException {
        return queryDb( query, max_lines_to_return, BASE_UNIPROT_URL );
    }

    final static String extractFrom( final String target, final String a ) {
        final int i_a = target.indexOf( a );
        return target.substring( i_a + a.length() ).trim();
    }

    final static String extractFromTo( final String target, final String a, final String b ) {
        final int i_a = target.indexOf( a );
        final int i_b = target.indexOf( b );
        if ( ( i_a < 0 ) || ( i_b < i_a ) ) {
            throw new IllegalArgumentException( "attempt to extract from \"" + target + "\" between \"" + a
                    + "\" and \"" + b + "\"" );
        }
        return target.substring( i_a + a.length(), i_b ).trim();
    }

    final static String extractTo( final String target, final String b ) {
        final int i_b = target.indexOf( b );
        return target.substring( 0, i_b ).trim();
    }

    private static void addDataFromDbToNode( final boolean allow_to_set_taxonomic_data,
                                             final int lines_to_return,
                                             final SortedSet<String> not_found,
                                             final PhylogenyNode node,
                                             final Accession acc ) throws IOException {
        SequenceDatabaseEntry db_entry = null;
        final String query = acc.getValue();
        if ( acc.getSource().equals( Source.UNIPROT.toString() ) ) {
            if ( DEBUG ) {
                System.out.println( "uniprot: " + query );
            }
            try {
                db_entry = obtainUniProtEntry( query, lines_to_return );
            }
            catch ( final FileNotFoundException e ) {
                // Eat this, and move to next.
            }
        }
        else if ( acc.getSource().equals( Source.REFSEQ.toString() ) ) {
            if ( DEBUG ) {
                System.out.println( "refseq: " + query );
            }
            try {
                db_entry = obtainRefSeqEntryFromEmbl( new Accession( query ), lines_to_return );
            }
            catch ( final FileNotFoundException e ) {
                // Eat this, and move to next.
            }
        }
        else if ( acc.getSource().equals( Source.EMBL.toString() ) || acc.getSource().equals( Source.NCBI.toString() )
                || acc.getSource().equals( Source.EMBL.toString() ) ) {
            if ( DEBUG ) {
                System.out.println( acc.toString() );
            }
            try {
                db_entry = obtainEmblEntry( acc, lines_to_return );
            }
            catch ( final FileNotFoundException e ) {
                // Eat this, and move to next.
            }
        }
        else if ( acc.getSource().equals( Source.GI.toString() ) ) {
            if ( DEBUG ) {
                System.out.println( "gi: " + query );
            }
            try {
                db_entry = obtainRefSeqEntryFromEmbl( new Accession( query ), lines_to_return );
            }
            catch ( final FileNotFoundException e ) {
                // Eat this, and move to next.
            }
        }
        if ( ( db_entry != null ) && !db_entry.isEmpty() ) {
            final Sequence seq = node.getNodeData().isHasSequence() ? node.getNodeData().getSequence() : new Sequence();
            if ( !ForesterUtil.isEmpty( db_entry.getAccession() ) ) {
                seq.setAccession( new Accession( db_entry.getAccession(), acc.getSource() ) );
            }
            if ( !ForesterUtil.isEmpty( db_entry.getSequenceName() ) ) {
                seq.setName( db_entry.getSequenceName() );
            }
            if ( !ForesterUtil.isEmpty( db_entry.getGeneName() ) ) {
                seq.setGeneName( db_entry.getGeneName() );
            }
            if ( !ForesterUtil.isEmpty( db_entry.getSequenceSymbol() ) ) {
                try {
                    seq.setSymbol( db_entry.getSequenceSymbol() );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    // Eat this exception.
                }
            }
            if ( ( db_entry.getGoTerms() != null ) && !db_entry.getGoTerms().isEmpty() ) {
                for( final GoTerm go : db_entry.getGoTerms() ) {
                    final Annotation ann = new Annotation( go.getGoId().getId() );
                    ann.setDesc( go.getName() );
                    seq.addAnnotation( ann );
                }
            }
            if ( ( db_entry.getCrossReferences() != null ) && !db_entry.getCrossReferences().isEmpty() ) {
                for( final Accession x : db_entry.getCrossReferences() ) {
                    seq.addCrossReference( x );
                }
            }
            if ( !ForesterUtil.isEmpty( db_entry.getChromosome() ) && !ForesterUtil.isEmpty( db_entry.getMap() ) ) {
                seq.setLocation( "chr " + db_entry.getChromosome() + ", " + db_entry.getMap() );
            }
            else if ( !ForesterUtil.isEmpty( db_entry.getChromosome() ) ) {
                seq.setLocation( "chr " + db_entry.getChromosome() );
            }
            else if ( !ForesterUtil.isEmpty( db_entry.getMap() ) ) {
                seq.setLocation( db_entry.getMap() );
            }
            final Taxonomy tax = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy() : new Taxonomy();
            if ( !ForesterUtil.isEmpty( db_entry.getTaxonomyScientificName() ) ) {
                tax.setScientificName( db_entry.getTaxonomyScientificName() );
            }
            if ( allow_to_set_taxonomic_data && !ForesterUtil.isEmpty( db_entry.getTaxonomyIdentifier() ) ) {
                tax.setIdentifier( new Identifier( db_entry.getTaxonomyIdentifier(), "uniprot" ) );
            }
            node.getNodeData().setTaxonomy( tax );
            node.getNodeData().setSequence( seq );
        }
        else {
            if ( node.isExternal() || !node.isEmpty() ) {
                not_found.add( node.toString() );
            }
        }
        try {
            Thread.sleep( 10 );// Sleep for 10 ms
        }
        catch ( final InterruptedException ie ) {
        }
    }

    private static String encode( final String str ) throws UnsupportedEncodingException {
        return URLEncoder.encode( str.trim(), URL_ENC );
    }

    private static List<UniProtTaxonomy> getTaxonomiesFromCommonName( final String cn, final int max_taxonomies_return )
            throws IOException {
        final List<String> result = getTaxonomyStringFromCommonName( cn, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    private static List<UniProtTaxonomy> getTaxonomiesFromScientificName( final String sn,
                                                                          final int max_taxonomies_return )
            throws IOException {
        final List<String> result = getTaxonomyStringFromScientificName( sn, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    private static List<String> getTaxonomyStringFromCommonName( final String cn, final int max_lines_to_return )
            throws IOException {
        return queryUniprot( "taxonomy/?query=common%3a%22" + encode( cn ) + "%22&format=tab", max_lines_to_return );
    }

    private static List<String> getTaxonomyStringFromId( final String id, final int max_lines_to_return )
            throws IOException {
        return queryUniprot( "taxonomy/?query=id%3a%22" + encode( id ) + "%22&format=tab", max_lines_to_return );
    }

    private static List<String> getTaxonomyStringFromScientificName( final String sn, final int max_lines_to_return )
            throws IOException {
        return queryUniprot( "taxonomy/?query=scientific%3a%22" + encode( sn ) + "%22&format=tab", max_lines_to_return );
    }

    private static List<String> getTaxonomyStringFromTaxonomyCode( final String code, final int max_lines_to_return )
            throws IOException {
        return queryUniprot( "taxonomy/?query=mnemonic%3a%22" + encode( code ) + "%22&format=tab", max_lines_to_return );
    }

    private final static boolean isAccessionAcceptable( final Accession acc ) {
        return ( !( ( acc == null ) || ForesterUtil.isEmpty( acc.getSource() ) || ForesterUtil.isEmpty( acc.getValue() ) || ( ( acc
                .getSource().equals( Source.UNIPROT.toString() ) )
                && ( acc.getSource().toString().equals( Source.EMBL.toString() ) ) && ( acc.getSource().toString()
                .equals( Source.REFSEQ.toString() ) ) ) ) );
    }

    private static List<UniProtTaxonomy> parseUniProtTaxonomy( final List<String> result ) throws IOException {
        final List<UniProtTaxonomy> taxonomies = new ArrayList<UniProtTaxonomy>();
        for( final String line : result ) {
            if ( ForesterUtil.isEmpty( line ) ) {
                // Ignore empty lines.
            }
            else if ( line.startsWith( "Taxon" ) ) {
                final String[] items = line.split( "\t" );
                if ( !( items[ 1 ].equalsIgnoreCase( "Mnemonic" ) && items[ 2 ].equalsIgnoreCase( "Scientific name" )
                        && items[ 3 ].equalsIgnoreCase( "Common name" ) && items[ 4 ].equalsIgnoreCase( "Synonym" )
                        && items[ 5 ].equalsIgnoreCase( "Other Names" ) && items[ 6 ].equalsIgnoreCase( "Reviewed" )
                        && items[ 7 ].equalsIgnoreCase( "Rank" ) && items[ 8 ].equalsIgnoreCase( "Lineage" ) ) ) {
                    throw new IOException( "Unreconized UniProt Taxonomy format: " + line );
                }
            }
            else {
                if ( line.split( "\t" ).length > 4 ) {
                    taxonomies.add( new UniProtTaxonomy( line ) );
                }
            }
        }
        return taxonomies;
    }
}
