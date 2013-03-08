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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceIdParser;

public final class SequenceDbWsTools {

    private static final boolean ALLOW_TAXONOMY_CODE_HACKS = true;                                         //TODO turn off for final realease!
    public final static String   BASE_UNIPROT_URL          = "http://www.uniprot.org/";
    public final static String   BASE_EMBL_DB_URL          = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/";
    public final static String   EMBL_DBS_EMBL             = "embl";
    public final static String   EMBL_DBS_REFSEQ_P         = "refseqp";
    public final static String   EMBL_DBS_REFSEQ_N         = "refseqn";
    private final static String  URL_ENC                   = "UTF-8";
    private final static boolean DEBUG                     = false;

    public static List<UniProtTaxonomy> getTaxonomiesFromCommonName( final String cn, final int max_taxonomies_return )
            throws IOException {
        final List<String> result = getTaxonomyStringFromCommonName( cn, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

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

    public static List<UniProtTaxonomy> getTaxonomiesFromScientificName( final String sn,
                                                                         final int max_taxonomies_return )
            throws IOException {
        // Hack!  Craniata? .. 
        if ( sn.equals( "Drosophila" ) ) {
            return uniProtTaxonomyToList( UniProtTaxonomy.DROSOPHILA_GENUS );
        }
        else if ( sn.equals( "Xenopus" ) ) {
            return uniProtTaxonomyToList( UniProtTaxonomy.XENOPUS_GENUS );
        }
        // else if ( sn.equals( "Nucleariidae and Fonticula group" ) ) {
        //     return hack( UniProtTaxonomy.NUCLEARIIDAE_AND_FONTICULA );
        // }
        final List<String> result = getTaxonomyStringFromScientificName( sn, max_taxonomies_return );
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
        if ( ALLOW_TAXONOMY_CODE_HACKS ) {
            final List<UniProtTaxonomy> l = resolveFakeTaxonomyCodes( max_taxonomies_return, my_code );
            if ( l != null ) {
                return l;
            }
        }
        final List<String> result = getTaxonomyStringFromTaxonomyCode( my_code, max_taxonomies_return );
        if ( result.size() > 0 ) {
            return parseUniProtTaxonomy( result );
        }
        return null;
    }

    public static SequenceDatabaseEntry obtainEmblEntry( final Identifier id, final int max_lines_to_return )
            throws IOException {
        final List<String> lines = queryEmblDb( id, max_lines_to_return );
        return EbiDbEntry.createInstanceFromPlainText( lines );
    }

    public static SequenceDatabaseEntry obtainRefSeqEntryFromEmbl( final Identifier id, final int max_lines_to_return )
            throws IOException {
        final List<String> lines = queryEmblDb( id, max_lines_to_return );
        return EbiDbEntry.createInstanceFromPlainTextForRefSeq( lines );
    }

    public static SortedSet<String> obtainSeqInformation( final Phylogeny phy,
                                                          final boolean ext_nodes_only,
                                                          final boolean allow_to_set_taxonomic_data,
                                                          final int lines_to_return ) throws IOException {
        final SortedSet<String> not_found = new TreeSet<String>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( ext_nodes_only && node.isInternal() ) {
                continue;
            }
            String query = null;
            Identifier id = null;
            Db db = Db.NONE;
            if ( node.getNodeData().isHasSequence() && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() )
                    && node.getNodeData().getSequence().getAccession().getValue().toLowerCase().startsWith( "uniprot" ) ) {
                query = node.getNodeData().getSequence().getAccession().getValue();
                db = Db.UNIPROT;
            }
            else if ( node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() )
                    && ( node.getNodeData().getSequence().getAccession().getValue().toLowerCase().startsWith( "embl" ) || node
                            .getNodeData().getSequence().getAccession().getValue().toLowerCase().startsWith( "ebi" ) ) ) {
                query = node.getNodeData().getSequence().getAccession().getValue();
                db = Db.EMBL;
            }
            else if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                if ( ( query = ForesterUtil.extractUniProtKbProteinSeqIdentifier( node ) ) != null ) {
                    db = Db.UNIPROT;
                }
                else if ( ( id = SequenceIdParser.parse( node.getName() ) ) != null ) {
                    if ( id.getProvider().equalsIgnoreCase( Identifier.NCBI ) ) {
                        db = Db.NCBI;
                    }
                    else if ( id.getProvider().equalsIgnoreCase( Identifier.REFSEQ ) ) {
                        db = Db.REFSEQ;
                    }
                }
            }
            if ( db == Db.NONE ) {
                not_found.add( node.getName() );
            }
            SequenceDatabaseEntry db_entry = null;
            if ( !ForesterUtil.isEmpty( query ) ) {
                if ( db == Db.UNIPROT ) {
                    if ( DEBUG ) {
                        System.out.println( "uniprot: " + query );
                    }
                    db_entry = obtainUniProtEntry( query, lines_to_return );
                }
                if ( ( db == Db.EMBL ) || ( ( db == Db.UNIPROT ) && ( db_entry == null ) ) ) {
                    if ( DEBUG ) {
                        System.out.println( "embl: " + query );
                    }
                    db_entry = obtainEmblEntry( new Identifier( query ), lines_to_return );
                    if ( ( db == Db.UNIPROT ) && ( db_entry != null ) ) {
                        db = Db.EMBL;
                    }
                }
            }
            else if ( ( db == Db.REFSEQ ) && ( id != null ) ) {
                db_entry = obtainRefSeqEntryFromEmbl( id, lines_to_return );
            }
            else if ( ( db == Db.NCBI ) && ( id != null ) ) {
                db_entry = obtainEmblEntry( id, lines_to_return ); //TODO ?
            }
            if ( ( db_entry != null ) && !db_entry.isEmpty() ) {
                final Sequence seq = node.getNodeData().isHasSequence() ? node.getNodeData().getSequence()
                        : new Sequence();
                if ( !ForesterUtil.isEmpty( db_entry.getAccession() ) ) {
                    String type = null;
                    if ( db == Db.EMBL ) {
                        type = "embl";
                    }
                    else if ( db == Db.UNIPROT ) {
                        type = "uniprot";
                    }
                    else if ( db == Db.NCBI ) {
                        type = "ncbi";
                    }
                    else if ( db == Db.REFSEQ ) {
                        type = "refseq";
                    }
                    seq.setAccession( new Accession( db_entry.getAccession(), type ) );
                }
                if ( !ForesterUtil.isEmpty( db_entry.getSequenceName() ) ) {
                    seq.setName( db_entry.getSequenceName() );
                }
                if ( !ForesterUtil.isEmpty( db_entry.getSequenceSymbol() ) ) {
                    seq.setSymbol( db_entry.getSequenceSymbol() );
                }
                final Taxonomy tax = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy()
                        : new Taxonomy();
                if ( !ForesterUtil.isEmpty( db_entry.getTaxonomyScientificName() ) ) {
                    tax.setScientificName( db_entry.getTaxonomyScientificName() );
                }
                if ( allow_to_set_taxonomic_data && !ForesterUtil.isEmpty( db_entry.getTaxonomyIdentifier() ) ) {
                    tax.setIdentifier( new Identifier( db_entry.getTaxonomyIdentifier(), "uniprot" ) );
                }
                node.getNodeData().setTaxonomy( tax );
                node.getNodeData().setSequence( seq );
            }
            else if ( db != Db.NONE ) {
                not_found.add( node.getName() );
            }
            try {
                Thread.sleep( 10 );// Sleep for 10 ms
            }
            catch ( final InterruptedException ie ) {
            }
        }
        return not_found;
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

    public static List<String> queryEmblDb( final Identifier id, final int max_lines_to_return ) throws IOException {
        final StringBuilder url_sb = new StringBuilder();
        url_sb.append( BASE_EMBL_DB_URL );
        if ( ForesterUtil.isEmpty( id.getProvider() ) || id.getProvider().equalsIgnoreCase( Identifier.NCBI ) ) {
            url_sb.append( SequenceDbWsTools.EMBL_DBS_EMBL );
            url_sb.append( '/' );
        }
        else if ( id.getProvider().equalsIgnoreCase( Identifier.REFSEQ ) ) {
            if ( id.getValue().toUpperCase().indexOf( 'P' ) == 1 ) {
                url_sb.append( SequenceDbWsTools.EMBL_DBS_REFSEQ_P );
                url_sb.append( '/' );
            }
            else {
                url_sb.append( SequenceDbWsTools.EMBL_DBS_REFSEQ_N );
                url_sb.append( '/' );
            }
        }
        return queryDb( id.getValue(), max_lines_to_return, url_sb.toString() );
    }

    public static List<String> queryUniprot( final String query, final int max_lines_to_return ) throws IOException {
        return queryDb( query, max_lines_to_return, BASE_UNIPROT_URL );
    }

    private static String encode( final String str ) throws UnsupportedEncodingException {
        return URLEncoder.encode( str.trim(), URL_ENC );
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

    private static List<UniProtTaxonomy> resolveFakeTaxonomyCodes( final int max_taxonomies_return, final String code )
            throws IOException {
        if ( code.equals( "CAP" ) ) {
            return getTaxonomiesFromId( "283909", max_taxonomies_return );
        }
        else if ( code.equals( "FUGRU" ) ) {
            return getTaxonomiesFromId( "31033", max_taxonomies_return );
        }
        else if ( code.equals( "GIALA" ) ) {
            return getTaxonomiesFromId( "5741", max_taxonomies_return );
        }
        else if ( code.equals( "TRIVE" ) ) {
            return getTaxonomiesFromId( "413071", max_taxonomies_return );
        }
        else if ( code.equals( "CAPOWC" ) ) {
            return getTaxonomiesFromId( "192875", max_taxonomies_return );
        }
        else if ( code.equals( "SPHARC" ) ) {
            return getTaxonomiesFromId( "667725", max_taxonomies_return );
        }
        else if ( code.equals( "THETRA" ) ) {
            return getTaxonomiesFromId( "529818", max_taxonomies_return );
        }
        else if ( code.equals( "CHLVUL" ) ) {
            return getTaxonomiesFromId( "574566", max_taxonomies_return );
        }
        else if ( code.equals( "CITCLE" ) ) {
            return getTaxonomiesFromId( "85681", max_taxonomies_return );
        }
        else if ( code.equals( "MYCPOP" ) ) {
            return getTaxonomiesFromId( "85929", max_taxonomies_return );
        }
        else if ( code.equals( "AGABB" ) ) {
            return getTaxonomiesFromId( "597362", max_taxonomies_return );
        }
        else if ( code.equals( "BAUCOM" ) ) {
            return getTaxonomiesFromId( "430998", max_taxonomies_return );
        }
        else if ( code.equals( "DICSQU" ) ) {
            return getTaxonomiesFromId( "114155", max_taxonomies_return );
        }
        else if ( code.equals( "FOMPIN" ) ) {
            return getTaxonomiesFromId( "40483", max_taxonomies_return );
        }
        else if ( code.equals( "HYDMA" ) ) {
            return getTaxonomiesFromId( "6085", max_taxonomies_return );
        }
        else if ( code.equals( "MYCFI" ) ) {
            return getTaxonomiesFromId( "83344", max_taxonomies_return );
        }
        else if ( code.equals( "OIDMAI" ) ) {
            return getTaxonomiesFromId( "78148", max_taxonomies_return );
        }
        else if ( code.equals( "OSTRC" ) ) {
            return getTaxonomiesFromId( "385169", max_taxonomies_return );
        }
        else if ( code.equals( "POSPL" ) ) {
            return getTaxonomiesFromId( "104341", max_taxonomies_return );
        }
        else if ( code.equals( "SAICOM" ) ) {
            return getTaxonomiesFromId( "5606", max_taxonomies_return );
        }
        else if ( code.equals( "SERLA" ) ) {
            return getTaxonomiesFromId( "85982", max_taxonomies_return );
        }
        else if ( code.equals( "SPORO" ) ) {
            return getTaxonomiesFromId( "40563", max_taxonomies_return );
        }
        else if ( code.equals( "ACRALC" ) ) {
            return getTaxonomiesFromId( "398408", max_taxonomies_return );
        }
        else if ( code.equals( "THITER" ) ) {
            return getTaxonomiesFromId( "35720", max_taxonomies_return );
        }
        else if ( code.equals( "MYCTHE" ) ) {
            return getTaxonomiesFromId( "78579", max_taxonomies_return );
        }
        else if ( code.equals( "CONPUT" ) ) {
            return getTaxonomiesFromId( "80637", max_taxonomies_return );
        }
        else if ( code.equals( "WOLCOC" ) ) {
            return getTaxonomiesFromId( "81056", max_taxonomies_return );
        }
        else if ( code.equals( "CLAGRA" ) ) {
            return getTaxonomiesFromId( "27339", max_taxonomies_return );
        }
        else if ( code.equals( "XANPAR" ) ) {
            return getTaxonomiesFromId( "107463", max_taxonomies_return );
        }
        else if ( code.equals( "HYDPIN" ) ) {
            return getTaxonomiesFromId( "388859", max_taxonomies_return );
        }
        else if ( code.equals( "SERLAC" ) ) {
            return getTaxonomiesFromId( "85982", max_taxonomies_return );
        }
        else {
            return null;
        }
    }

    private static List<UniProtTaxonomy> uniProtTaxonomyToList( final UniProtTaxonomy tax ) {
        final List<UniProtTaxonomy> l = new ArrayList<UniProtTaxonomy>();
        l.add( tax );
        return l;
    }

    public enum Db {
        UNIPROT, EMBL, NCBI, NONE, REFSEQ;
    }
}
