// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.archaeopteryx;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.ws.seqdb.TaxonLineage;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Tests for {@link CladeAnnotationReport}: one row per tip with its source and note, the ignored internal
 * annotations, the per-rank summary, the counts and the TSV. Headless; an in-memory service.
 */
public final class CladeAnnotationReportTest {

    public static void main( final String[] args ) {
        System.out.println( "CladeAnnotationReport: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return testRowsAndSummary() && testTreeSourceAndEdges() && testTsv();
    }

    /** The fixture: "order Rodentia" wrongly over a cat (id 9685 -> Carnivora in the cache) and a bare tip "u";
     *  a plain node over an ambiguous "Bacillus" and an un-looked-up "notyet"; a dog annotated at order itself. */
    private static Phylogeny fixture( final FakeService svc ) {
        svc.know( "9685", "order", "Carnivora" );
        svc.ambiguous( "Bacillus", new TaxonLineage( "55087", "genus", "Bacillus", null, null ),
                       new TaxonLineage( "1386", "genus", "Bacillus", null, null ) );
        final PhylogenyNode rodents = internal( "order", "Rodentia" );
        rodents.addAsChild( leaf( "cat", "Felis catus", "species", "9685" ) );
        rodents.addAsChild( leaf( "u", null, null, null ) );
        final PhylogenyNode plain = new PhylogenyNode();
        plain.addAsChild( leaf( "Bacillus", null, null, null ) );
        plain.addAsChild( leaf( "notyet", null, null, null ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( rodents );
        root.addAsChild( plain );
        root.addAsChild( leaf( "dog", "Carnivora", "order", null ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean testRowsAndSummary() {
        final FakeService svc = new FakeService();
        final Phylogeny phy = fixture( svc );
        final CladeAnnotationReport r = CladeAnnotationReport.of( phy, Collections.singletonList( "order" ), svc,
                                                                  TaxonSource.DATABASE );
        // one row per tip, in tree order, then the ignored internal annotation
        final List<String> nodes = new ArrayList<String>();
        for( final CladeAnnotationReport.Row row : r.getRows() ) {
            nodes.add( row.getNode() );
        }
        if ( !nodes.equals( Arrays.asList( "cat", "u", "Bacillus", "notyet", "dog", "[internal] Rodentia" ) ) ) {
            return fail( "rows must be the tips in tree order, then the ignored internal annotations; got " + nodes );
        }
        final Map<String, String> taxa = r.taxaByNodeForTest( "order" );
        if ( !"Carnivora".equals( taxa.get( "cat" ) ) || !"Rodentia".equals( taxa.get( "u" ) ) || !"".equals( taxa.get( "Bacillus" ) )
                || !"".equals( taxa.get( "notyet" ) ) || !"Carnivora".equals( taxa.get( "dog" ) ) ) {
            return fail( "wrong taxa: " + taxa );
        }
        final CladeAnnotationReport.Row cat = r.getRows().get( 0 );
        if ( !cat.getSource().equals( TaxonPlacement.Source.DATABASE_BY_ID.getLabel() )
                || !cat.getNote().contains( "the tree's ancestor says \"Rodentia\"" ) || !cat.getNote().contains( "9685" ) ) {
            return fail( "the cat row must name its source and the overridden ancestor; got " + cat.getSource() + " / " + cat.getNote() );
        }
        final CladeAnnotationReport.Row bacillus = r.getRows().get( 2 );
        if ( !bacillus.getSource().equals( TaxonPlacement.Source.UNPLACED_AMBIGUOUS.getLabel() )
                || !bacillus.getNote().contains( "2 taxa" ) || !bacillus.getNote().contains( "55087" ) ) {
            return fail( "an ambiguous tip lists its candidates; got " + bacillus.getNote() );
        }
        final CladeAnnotationReport.Row ignored = r.getRows().get( 5 );
        if ( !ignored.isInternal() || !"Rodentia".equals( ignored.getTaxon() )
                || !CladeAnnotationReport.IGNORED_SOURCE.equals( ignored.getSource() )
                || !ignored.getNote().contains( "Carnivora (1)" ) ) {
            return fail( "the contradicted internal annotation must be reported with who contradicts it; got "
                    + ignored.getSource() + " / " + ignored.getNote() );
        }
        if ( ( r.getConflictCount() != 1 ) || ( r.getIgnoredInternalCount() != 1 ) || ( r.getUnplacedCount() != 2 ) ) {
            return fail( "counts: conflicts " + r.getConflictCount() + " ignored " + r.getIgnoredInternalCount()
                    + " unplaced " + r.getUnplacedCount() );
        }
        final String line = r.getSummary();
        for( final String piece : new String[] { "order: 3 of 5 tips placed", "database, by its NCBI id 1",
                "ancestor's annotation 1", "own annotation 1", "2 unplaced", "ambiguous name 1", "not looked up yet 1",
                "the database overrode the tree on 1", "1 internal annotation ignored (Rodentia)" } ) {
            if ( !line.contains( piece ) ) {
                return fail( "summary lacks \"" + piece + "\": " + line );
            }
        }
        if ( r.getSummaryLines().size() != 1 ) {
            return fail( "one summary line per rank" );
        }
        // the same wrong annotation on several nodes is named once, with its count
        final String grouped = CladeAnnotationReport.summaryLine( "genus", 3, 3, new java.util.LinkedHashMap<TaxonPlacement.Source, Integer>(),
                                                                  0, Arrays.asList( "Bacteria", "Bacteria", "Bacteria", "Mus" ) );
        if ( !grouped.endsWith( "4 internal annotations ignored (Bacteria \u00d73, Mus)" ) ) {
            return fail( "repeated ignored names are grouped with a count; got " + grouped );
        }
        return true;
    }

    private static boolean testTreeSourceAndEdges() {
        final FakeService svc = new FakeService();
        final Phylogeny phy = fixture( svc );
        final CladeAnnotationReport t = CladeAnnotationReport.of( phy, Collections.singletonList( "order" ), svc,
                                                                  TaxonSource.TREE );
        if ( ( t.getConflictCount() != 0 ) || ( t.getIgnoredInternalCount() != 0 ) || ( t.getUnplacedCount() != 2 )
                || ( t.getRows().size() != 5 ) ) {
            return fail( "from the tree the cat is a rodent, nothing conflicts and nothing is ignored; got " + t.getSummary() );
        }
        if ( !"Rodentia".equals( t.taxaByNodeForTest( "order" ).get( "cat" ) ) ) {
            return fail( "the tree source places the cat where the tree says" );
        }
        // two ranks: two summary lines, rows of each rank tagged with it
        final CladeAnnotationReport two = CladeAnnotationReport.of( phy, Arrays.asList( "order", "genus" ), svc,
                                                                    TaxonSource.TREE );
        if ( ( two.getSummaryLines().size() != 2 ) || !two.getSummaryLines().get( 1 ).startsWith( "genus: 0 of 5" )
                || !"genus".equals( two.getRows().get( two.getRows().size() - 1 ).getRank() ) ) {
            return fail( "a report over two ranks has a line and rows per rank; got " + two.getSummary() );
        }
        // edges: an empty tree, a null rank list, a blank rank
        if ( !CladeAnnotationReport.of( new Phylogeny(), Collections.singletonList( "order" ), svc, TaxonSource.TREE ).getRows().isEmpty()
                || !CladeAnnotationReport.of( phy, null, svc, TaxonSource.TREE ).getRows().isEmpty()
                || !CladeAnnotationReport.of( phy, Collections.singletonList( "" ), svc, TaxonSource.TREE ).getRows().isEmpty() ) {
            return fail( "an empty tree, a null rank list or a blank rank yield an empty report" );
        }
        // a nameless tip is labelled by its taxonomy, a nameless taxonomy-less one "(unnamed)"
        final PhylogenyNode anon = leaf( null, "Felis catus", "species", null );
        if ( !"Felis catus".equals( CladeAnnotationReport.nodeLabel( anon ) )
                || !"(unnamed)".equals( CladeAnnotationReport.nodeLabel( new PhylogenyNode() ) ) ) {
            return fail( "node labels fall back from name to taxonomy to (unnamed)" );
        }
        return true;
    }

    private static boolean testTsv() {
        final FakeService svc = new FakeService();
        final Phylogeny phy = fixture( svc );
        phy.getNode( "notyet" ).setName( "not\tyet\nnow" ); // a tab and a line break in a tip name
        final CladeAnnotationReport r = CladeAnnotationReport.of( phy, Collections.singletonList( "order" ), svc,
                                                                  TaxonSource.DATABASE );
        final String[] lines = r.toTsv().split( "\n", -1 );
        if ( !"Rank\tNode\tTaxon\tSource\tNote".equals( lines[ 0 ] ) ) {
            return fail( "TSV header: " + lines[ 0 ] );
        }
        if ( ( lines.length != r.getRows().size() + 2 ) || !lines[ lines.length - 1 ].isEmpty() ) {
            return fail( "one line per row, newline-terminated; got " + lines.length + " pieces" );
        }
        if ( !lines[ 1 ].startsWith( "order\tcat\tCarnivora\t" + TaxonPlacement.Source.DATABASE_BY_ID.getLabel() + "\t" ) ) {
            return fail( "first data line: " + lines[ 1 ] );
        }
        for( final String l : lines ) {
            if ( l.split( "\t", -1 ).length != ( l.isEmpty() ? 1 : 5 ) ) {
                return fail( "every line must have exactly five cells (a cell's tab became a space): " + l );
            }
        }
        if ( !lines[ 4 ].contains( "not yet now" ) ) {
            return fail( "a tab or line break inside a cell becomes a space; got " + lines[ 4 ] );
        }
        return true;
    }

    private static PhylogenyNode internal( final String rank, final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        try {
            t.setRank( rank );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode leaf( final String name, final String sci, final String rank, final String id ) {
        final PhylogenyNode n = new PhylogenyNode();
        if ( name != null ) {
            n.setName( name );
        }
        if ( sci != null ) {
            final Taxonomy t = new Taxonomy();
            t.setScientificName( sci );
            if ( rank != null ) {
                try {
                    t.setRank( rank );
                }
                catch ( final Exception e ) {
                    throw new RuntimeException( e );
                }
            }
            if ( id != null ) {
                t.setIdentifier( new Identifier( id, "ncbi" ) );
            }
            n.getNodeData().setTaxonomy( t );
        }
        return n;
    }

    /** Cache-only in-memory service: {@code know} puts a lineage straight into the cache; {@code ambiguous} caches
     *  an EMPTY and the candidates. */
    private static final class FakeService implements TaxonomicLineageService {

        private final Map<String, TaxonLineage>       _cache     = new HashMap<String, TaxonLineage>();
        private final Map<String, List<TaxonLineage>> _ambiguous = new HashMap<String, List<TaxonLineage>>();

        void know( final String key, final String rank, final String name ) {
            _cache.put( key.toLowerCase( Locale.ROOT ), new TaxonLineage( null, null, null, null,
                    Collections.singletonList( new TaxonLineage.Ancestor( name, rank, null ) ) ) );
        }

        void ambiguous( final String name, final TaxonLineage... candidates ) {
            _cache.put( name.toLowerCase( Locale.ROOT ), TaxonLineage.EMPTY );
            _ambiguous.put( name.toLowerCase( Locale.ROOT ), Arrays.asList( candidates ) );
        }

        @Override
        public TaxonLineage lineageOf( final String taxon ) {
            return ( taxon == null ) ? null : _cache.get( taxon.trim().toLowerCase( Locale.ROOT ) );
        }

        @Override
        public TaxonLineage fetch( final String taxon ) {
            throw new UnsupportedOperationException( "the report must never fetch" );
        }

        @Override
        public List<TaxonLineage> ambiguousCandidates( final String taxon ) {
            final List<TaxonLineage> c = ( taxon == null ) ? null : _ambiguous.get( taxon.trim().toLowerCase( Locale.ROOT ) );
            return ( c == null ) ? Collections.<TaxonLineage> emptyList() : c;
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  CladeAnnotationReport test failed: " + msg );
        return false;
    }

    private CladeAnnotationReportTest() {
    }
}
