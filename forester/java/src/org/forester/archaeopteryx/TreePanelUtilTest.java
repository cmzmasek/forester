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

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.SortedSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.ws.seqdb.AccessionAwareLineageService;
import org.forester.ws.seqdb.Organism;
import org.forester.ws.seqdb.OrganismSource;
import org.forester.ws.seqdb.RankedLineage;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Unit tests for {@link TreePanelUtil}. Lives in the {@code org.forester.archaeopteryx} package
 * because the methods under test are package-private. Run standalone via {@link #main(String[])},
 * or as part of the suite via {@link #test()} (called from {@code org.forester.test.Test}).
 */
public final class TreePanelUtilTest {

    public static void main( final String[] args ) {
        System.out.println( "TreePanelUtil: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return testYDistanceToAvoidLabelOverlap() && testSupportSymbolMath() && testDetectConfidenceScaleMax()
                && testCapEntries() && testTaxonomyLabel() && testRankColorization() && testTipQueryName()
                && testCladeBands() && testRankColorizationViaSequenceIds() && testInternalLabelAboveBranchLayout()
                && testAbbreviateScientificName() && testSupportColor() && testScaleGridLines();
    }

    /**
     * Distance grid lines step by {@code spacing} from one spacing right of the root up to and including
     * the deepest tip; non-positive spacing or a zero-depth tree yields none.
     */
    private static boolean testScaleGridLines() {
        // origin 10, spacing 5, max 30 -> 15,20,25,30 (max included exactly); assert EVERY element so a
        // regression to accumulation or an off-by-one start (origin + i*spacing) is caught, not just count.
        final float[] a = TreePanelUtil.scaleGridLineXs( 10f, 5f, 30f );
        if ( ( a.length != 4 ) || ( a[ 0 ] != 15f ) || ( a[ 1 ] != 20f ) || ( a[ 2 ] != 25f ) || ( a[ 3 ] != 30f ) ) {
            return fail( "expected lines 15,20,25,30; got " + java.util.Arrays.toString( a ) );
        }
        // non-integer origin (exactly-representable .25 steps to avoid float-equality flak): lines are
        // origin + (i+1)*spacing -- the first is 10.75, NOT 10.5 (which an off-by-one start would give).
        final float[] b = TreePanelUtil.scaleGridLineXs( 10.5f, 0.25f, 11.5f );
        if ( ( b.length != 4 ) || ( b[ 0 ] != 10.75f ) || ( b[ 1 ] != 11.0f ) || ( b[ 2 ] != 11.25f )
                || ( b[ 3 ] != 11.5f ) ) {
            return fail( "non-integer origin: expected 10.75,11.0,11.25,11.5; got " + java.util.Arrays.toString( b ) );
        }
        // max not on a boundary -> last line is the largest multiple <= max (no overshoot past 30)
        if ( TreePanelUtil.scaleGridLineXs( 10f, 5f, 32f ).length != 4 ) {
            return fail( "must not place a line beyond max_x" );
        }
        // no room for even one line, and degenerate spacings -> empty
        if ( ( TreePanelUtil.scaleGridLineXs( 10f, 5f, 12f ).length != 0 )
                || ( TreePanelUtil.scaleGridLineXs( 10f, 0f, 30f ).length != 0 )
                || ( TreePanelUtil.scaleGridLineXs( 10f, -1f, 30f ).length != 0 ) ) {
            return fail( "too-small / non-positive spacing must yield no grid lines" );
        }
        return true;
    }

    /**
     * COLOR_BRANCHES gradient: full strong color at fraction 1, fading toward the background as support
     * drops, monotonic in between, fraction clamped. Theme-agnostic (works from any strong/background pair).
     */
    private static boolean testSupportColor() {
        final java.awt.Color black = java.awt.Color.BLACK;
        final java.awt.Color white = java.awt.Color.WHITE;
        // strong support -> full branch color (no fade)
        if ( !black.equals( TreePanelUtil.supportColor( 1.0, black, white ) ) ) {
            return fail( "fraction 1 must be the full strong color; got " + TreePanelUtil.supportColor( 1.0, black, white ) );
        }
        // weakest support fades 80% toward the background: 0 + 0.8*255 = 204
        final java.awt.Color weak = TreePanelUtil.supportColor( 0.0, black, white );
        if ( ( weak.getRed() != 204 ) || ( weak.getGreen() != 204 ) || ( weak.getBlue() != 204 ) ) {
            return fail( "fraction 0 must fade 80% toward the background; got " + weak );
        }
        // monotonic: stronger support is darker (closer to the strong color) on a light background
        if ( !( TreePanelUtil.supportColor( 0.25, black, white ).getRed()
                > TreePanelUtil.supportColor( 0.75, black, white ).getRed() ) ) {
            return fail( "stronger support must be closer to the strong color" );
        }
        // theme-aware: on a dark background the weak color fades toward dark, not light
        final java.awt.Color dark_bg = new java.awt.Color( 30, 30, 30 );
        final java.awt.Color light_branch = new java.awt.Color( 230, 230, 230 );
        if ( TreePanelUtil.supportColor( 0.0, light_branch, dark_bg ).getRed() >= light_branch.getRed() ) {
            return fail( "on a dark theme, weak support must fade toward the (dark) background" );
        }
        // out-of-range fractions clamp
        if ( !black.equals( TreePanelUtil.supportColor( 1.5, black, white ) )
                || !weak.equals( TreePanelUtil.supportColor( -0.5, black, white ) ) ) {
            return fail( "fraction must clamp to 0..1" );
        }
        return true;
    }

    /**
     * The binomial abbreviation ("Homo sapiens" -&gt; "H. sapiens", genus initial + full species, extra
     * epithets kept) and -- the point of the guard -- that malformed names (single token, trailing/leading
     * whitespace) are returned verbatim instead of throwing an {@link ArrayIndexOutOfBoundsException}.
     */
    private static boolean testAbbreviateScientificName() {
        if ( !"H. sapiens".equals( TreePanelUtil.abbreviateScientificName( "Homo sapiens" ) ) ) {
            return fail( "binomial: 'Homo sapiens' -> 'H. sapiens'; got "
                    + TreePanelUtil.abbreviateScientificName( "Homo sapiens" ) );
        }
        if ( !"H. sapiens neanderthalensis"
                .equals( TreePanelUtil.abbreviateScientificName( "Homo sapiens neanderthalensis" ) ) ) {
            return fail( "trinomial: genus abbreviated, further epithets kept verbatim" );
        }
        if ( !"E. coli".equals( TreePanelUtil.abbreviateScientificName( "Escherichia  coli" ) ) ) {
            return fail( "collapsed internal whitespace: 'Escherichia  coli' -> 'E. coli'" );
        }
        // malformed inputs must not throw and must come back unchanged
        if ( !"Homo ".equals( TreePanelUtil.abbreviateScientificName( "Homo " ) ) ) {
            return fail( "single-token-with-trailing-space must be returned verbatim, not throw" );
        }
        if ( !"Homo".equals( TreePanelUtil.abbreviateScientificName( "Homo" ) ) ) {
            return fail( "single token must be returned verbatim" );
        }
        if ( !" sapiens".equals( TreePanelUtil.abbreviateScientificName( " sapiens" ) ) ) {
            return fail( "leading whitespace (empty first token) must be returned verbatim" );
        }
        return true;
    }

    /**
     * The publication-style internal label sits to the LEFT of the node, right-aligned so the
     * (rightmost) node-data segment ends just left of the node, with the taxonomy segment to its left
     * and a baseline just above the branch. Verifies the right-alignment invariant for both-present,
     * data-only and taxonomy-only cases, and the left-edge clamp for a label that would overflow.
     */
    private static boolean testInternalLabelAboveBranchLayout() {
        final float node_x = 100f;
        final float node_y = 50f;
        final int hbs = 3;       // half box size
        final int gap = 5;
        final int descent = 4;
        final float min_x = 2f;
        final float right = node_x - hbs - 2;                  // 95: the node's left edge
        // both segments present: taxo (40) | gap (5) | data (30), data ending at the node
        final float[] both = TreePanelUtil.internalLabelAboveBranchLayout( node_x, node_y, hbs, 40, 30, gap, descent, min_x );
        if ( ( both[ 0 ] != 20f ) || ( both[ 1 ] != 65f ) || ( both[ 2 ] != 45f ) ) {
            return fail( "both segments: expected {20,65,45}; got " + both[ 0 ] + "," + both[ 1 ] + "," + both[ 2 ] );
        }
        if ( both[ 1 ] + 30 != right ) {
            return fail( "node-data segment must end at the node's left edge" );
        }
        if ( both[ 0 ] + 40 != ( both[ 1 ] - gap ) ) {
            return fail( "taxonomy segment must end one gap left of the node-data segment" );
        }
        // data only: no gap applied, data still right-aligned to the node
        final float[] data_only = TreePanelUtil.internalLabelAboveBranchLayout( node_x, node_y, hbs, 0, 30, gap, descent, min_x );
        if ( ( data_only[ 1 ] != 65f ) || ( data_only[ 1 ] + 30 != right ) ) {
            return fail( "data-only: node-data must end at the node's left edge with no gap; got " + data_only[ 1 ] );
        }
        // taxonomy only: no gap applied, taxonomy right-aligned to the node
        final float[] taxo_only = TreePanelUtil.internalLabelAboveBranchLayout( node_x, node_y, hbs, 40, 0, gap, descent, min_x );
        if ( ( taxo_only[ 0 ] != 55f ) || ( taxo_only[ 0 ] + 40 != right ) ) {
            return fail( "taxo-only: taxonomy must end at the node's left edge with no gap; got " + taxo_only[ 0 ] );
        }
        // baseline sits above the branch line (smaller y is higher on screen)
        if ( both[ 2 ] >= node_y ) {
            return fail( "label baseline must be above the branch line" );
        }
        // left-edge clamp: a wide label on a node near the root would right-align to a negative x
        // (taxo_x = 30-3-2-30-5-40 = -50); it must shift right so the leftmost glyph is at min_x, while
        // the internal taxo->data spacing (gap + taxo_width = 45) is preserved.
        final float[] clamped = TreePanelUtil.internalLabelAboveBranchLayout( 30f, node_y, hbs, 40, 30, gap, descent, min_x );
        if ( clamped[ 0 ] != min_x ) {
            return fail( "clamp: leftmost (taxonomy) segment must start at min_x; got " + clamped[ 0 ] );
        }
        if ( ( clamped[ 1 ] - clamped[ 0 ] ) != ( gap + 40 ) ) {
            return fail( "clamp: taxo->data spacing must be preserved after the shift; got "
                    + ( clamped[ 1 ] - clamped[ 0 ] ) );
        }
        return true;
    }

    /**
     * End-to-end (no network) of the very-common case the user hit: a tree whose tips are UniProt/SwissProt
     * <i>sequence</i> identifiers (no taxonomy on the nodes). The default service is an
     * {@link AccessionAwareLineageService}, so the fetch pass bridges each accession to its organism's
     * lineage and the cache-only colorize then places every tip -- exactly the production flow, with the
     * network replaced by in-memory fakes.
     */
    private static boolean testRankColorizationViaSequenceIds() {
        // delegate keyed by NCBI tax-id (what the bridge resolves an accession's organism to)
        final FakeLineageService delegate = new FakeLineageService();
        delegate.know( "9606", lineage( "order", "Primates" ) );  // human
        delegate.know( "10090", lineage( "order", "Rodentia" ) ); // mouse
        // organism source: accession -> organism NCBI tax-id
        final FakeOrganismSource seqs = new FakeOrganismSource();
        seqs.know( "P12345", "9606" );   // UniProt accession, human
        seqs.know( "RL7_HUMAN", "9606" );// SwissProt entry name, human
        seqs.know( "P63017", "10090" );  // UniProt accession, mouse
        final AccessionAwareLineageService svc = new AccessionAwareLineageService( delegate, seqs );

        // ((P12345, RL7_HUMAN), Q9MOUS): two human tips form one clade, the mouse tip is a second
        final PhylogenyNode humans = new PhylogenyNode();
        humans.addAsChild( bareLeaf( "P12345" ) );
        humans.addAsChild( bareLeaf( "RL7_HUMAN" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( humans );
        root.addAsChild( bareLeaf( "P63017" ) );
        final Phylogeny tree = new Phylogeny();
        tree.setRoot( root );
        tree.externalNodesHaveChanged();

        // every tip is unplaceable from the tree alone (sequence ids, not taxa)
        final SortedSet<String> unresolved = TreePanelUtil.unresolvedTipTaxa( tree, "order", svc );
        if ( unresolved.size() != 3 ) {
            return fail( "all three sequence-id tips must need online resolution; got " + unresolved );
        }
        try {
            for( final String name : unresolved ) {
                svc.fetch( name ); // the off-EDT fetch pass: bridges acc -> organism -> ranked lineage
            }
        }
        catch ( final java.io.IOException e ) {
            return fail( "fake fetch must not throw: " + e );
        }
        // cache-only colorize now places all three tips by order
        final Map<PhylogenyNode, String> assignment = TreePanelUtil.assignTipsToRankTaxon( tree, "order", svc );
        if ( assignment.size() != 3 ) {
            return fail( "all three tips must be placed at rank order after the fetch; got " + assignment );
        }
        final Map<String, Color> legend = new LinkedHashMap<String, Color>();
        final int colorizations = TreePanelUtil.colorPhylogenyAccordingToRanks( tree, "order", svc, legend );
        if ( colorizations != 2 ) {
            return fail( "expected 2 colorized clades (Primates clade + Rodentia tip); got " + colorizations );
        }
        if ( ( legend.size() != 2 ) || !legend.containsKey( "Primates" ) || !legend.containsKey( "Rodentia" ) ) {
            return fail( "legend must hold the two resolved orders; got " + legend.keySet() );
        }
        return true;
    }

    /** A leaf with only a node name (no taxonomy) -- so tipQueryName returns the name (a sequence id). */
    private static PhylogenyNode bareLeaf( final String node_name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( node_name );
        return n;
    }

    /** In-memory organism source: maps an accession value to an organism NCBI tax-id. */
    private static final class FakeOrganismSource implements OrganismSource {

        private final Map<String, String> _organism_id = new HashMap<String, String>();

        void know( final String accession, final String organism_id ) {
            _organism_id.put( accession.toUpperCase( Locale.ROOT ), organism_id );
        }

        @Override
        public Organism organismOf( final Accession acc ) {
            final String id = _organism_id.get( acc.getValue().toUpperCase( Locale.ROOT ) );
            return ( id == null ) ? Organism.EMPTY : new Organism( id, null );
        }
    }

    /**
     * tipQueryName prefers the TIP's own identity (scientific name, else code, else common name, else
     * node name) and only falls back to an ancestor's scientific name when the tip carries no identity
     * at all -- so a code-only tip under a sci-name ancestor is queried by its own (more specific) code,
     * not the ancestor's coarser name.
     */
    private static boolean testTipQueryName() {
        // tip's own scientific name wins
        final PhylogenyNode sci = new PhylogenyNode();
        final Taxonomy st = new Taxonomy();
        st.setScientificName( "Felis" );
        sci.getNodeData().setTaxonomy( st );
        if ( !"Felis".equals( TreePanelUtil.tipQueryName( sci ) ) ) {
            return fail( "tip's own scientific name must be used" );
        }
        // a code-only tip under a sci-name-bearing ancestor returns the TIP's code, not the ancestor's name
        final PhylogenyNode ancestor = new PhylogenyNode();
        final Taxonomy at = new Taxonomy();
        at.setScientificName( "Mammalia" );
        ancestor.getNodeData().setTaxonomy( at );
        final PhylogenyNode code_leaf = new PhylogenyNode();
        final Taxonomy ct = new Taxonomy();
        try {
            ct.setTaxonomyCode( "FELCA" ); // validated against the taxonomy-code format
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        code_leaf.getNodeData().setTaxonomy( ct );
        ancestor.addAsChild( code_leaf );
        if ( !"FELCA".equals( TreePanelUtil.tipQueryName( code_leaf ) ) ) {
            return fail( "a code-only tip must be queried by its own code, not the ancestor's scientific name" );
        }
        // a tip with no taxonomy but a node name uses the node name
        final PhylogenyNode named = new PhylogenyNode();
        named.setName( "Homo_sapiens" );
        if ( !"Homo_sapiens".equals( TreePanelUtil.tipQueryName( named ) ) ) {
            return fail( "a tip with only a node name must use the node name" );
        }
        // a tip with no identity at all falls back to the nearest ancestor's scientific name
        final PhylogenyNode anc2 = new PhylogenyNode();
        final Taxonomy a2 = new Taxonomy();
        a2.setScientificName( "Carnivora" );
        anc2.getNodeData().setTaxonomy( a2 );
        final PhylogenyNode bare = new PhylogenyNode();
        anc2.addAsChild( bare );
        if ( !"Carnivora".equals( TreePanelUtil.tipQueryName( bare ) ) ) {
            return fail( "an identity-less tip must fall back to the ancestor's scientific name" );
        }
        return true;
    }

    /**
     * The assignment-based rank colorizer (the major-flaw fix): every tip is placed under its taxon
     * at the chosen rank -- from an in-tree rank annotation when present, else the lineage service's
     * cache -- and each maximal monophyletic run of one taxon is colored. The test tree mixes a
     * Rodentia clade annotated at rank=order with two cats/dogs that are only genus-annotated and
     * sit in different parts of the tree (paraphyletic Carnivora), plus one tip the DB can't resolve.
     */
    private static boolean testRankColorization() {
        // a fake taxonomy DB: Felis and Canis resolve to order Carnivora; "Nonexistus" is unknown
        final FakeLineageService svc = new FakeLineageService();
        svc.know( "Felis", lineage( "class", "Mammalia", "order", "Carnivora", "genus", "Felis" ) );
        svc.know( "Canis", lineage( "class", "Mammalia", "order", "Carnivora", "genus", "Canis" ) );

        final Phylogeny tree = mammalTree();
        final PhylogenyNode felis = findLeaf( tree, "Felis" );
        final PhylogenyNode canis = findLeaf( tree, "Canis" );

        // --- phase 1: cache empty, so only the in-tree-annotated Rodentia tips are placeable ---
        Map<PhylogenyNode, String> assignment = TreePanelUtil.assignTipsToRankTaxon( tree, "order", svc );
        if ( assignment.size() != 2 ) {
            return fail( "with an empty cache only the 2 Rodentia tips should be placeable, got " + assignment.size() );
        }
        final SortedSet<String> unresolved = TreePanelUtil.unresolvedTipTaxa( tree, "order", svc );
        // Felis, Canis, and the unknown tip's query name must all be flagged for fetching
        if ( !unresolved.contains( "Felis" ) || !unresolved.contains( "Canis" ) || !unresolved.contains( "Nonexistus" ) ) {
            return fail( "unresolved tip taxa must include Felis, Canis and Nonexistus; got " + unresolved );
        }

        // --- background fetch (simulated): resolve every unresolved name once ---
        try {
            for( final String name : unresolved ) {
                svc.fetch( name );
            }
        }
        catch ( final Exception e ) {
            return fail( "fake fetch must not throw: " + e );
        }
        // after fetching, nothing is left to fetch (no repeated prompts)
        if ( !TreePanelUtil.unresolvedTipTaxa( tree, "order", svc ).isEmpty() ) {
            return fail( "after fetching all names, unresolvedTipTaxa must be empty" );
        }

        // --- phase 2: now Felis and Canis place at Carnivora; the unknown tip stays unplaced ---
        assignment = TreePanelUtil.assignTipsToRankTaxon( tree, "order", svc );
        if ( assignment.size() != 4 ) {
            return fail( "after fetch, 4 tips (2 Rodentia + Felis + Canis) should be placed, got " + assignment.size() );
        }
        if ( !"Carnivora".equals( assignment.get( felis ) ) || !"Carnivora".equals( assignment.get( canis ) ) ) {
            return fail( "Felis and Canis must be placed under order Carnivora" );
        }

        // maximal roots: the Rodentia clade + Felis + Canis (paraphyletic Carnivora -> two roots, one taxon)
        final Map<PhylogenyNode, String> roots = TreePanelUtil.maximalMonochromaticRoots( tree, assignment );
        if ( roots.size() != 3 ) {
            return fail( "expected 3 maximal monochromatic roots (Rodentia clade, Felis, Canis), got " + roots.size() );
        }
        int carnivora_roots = 0;
        for( final String t : roots.values() ) {
            if ( "Carnivora".equals( t ) ) {
                ++carnivora_roots;
            }
        }
        if ( carnivora_roots != 2 ) {
            return fail( "paraphyletic Carnivora must yield two same-taxon roots, got " + carnivora_roots );
        }

        // full colorize: 3 clades colored, legend has the 2 distinct taxa, Felis and Canis share a color
        final Map<String, Color> legend = new LinkedHashMap<String, Color>();
        final int colorized = TreePanelUtil.colorPhylogenyAccordingToRanks( tree, "order", svc, legend );
        if ( colorized != 3 ) {
            return fail( "colorize should report 3 colored clades, got " + colorized );
        }
        if ( ( legend.size() != 2 ) || !legend.containsKey( "Carnivora" ) || !legend.containsKey( "Rodentia" ) ) {
            return fail( "legend must list exactly Carnivora and Rodentia; got " + legend.keySet() );
        }
        final Color felis_c = felis.getBranchData().getBranchColor().getValue();
        final Color canis_c = canis.getBranchData().getBranchColor().getValue();
        if ( ( felis_c == null ) || !felis_c.equals( canis_c ) ) {
            return fail( "Felis and Canis (same order) must get the same color" );
        }
        if ( felis_c.equals( legend.get( "Rodentia" ) ) ) {
            return fail( "Carnivora and Rodentia must get distinct colors" );
        }
        // the unresolvable tip is never colored
        if ( findLeaf( tree, "x_unknown" ).getBranchData().getBranchColor() != null ) {
            return fail( "an unplaceable tip must be left uncolored, not swept into a neighbor's color" );
        }
        return true;
    }

    /** Clade bands reuse the same assignment as the colorizer: one band per maximal-monophyletic clade. */
    private static boolean testCladeBands() {
        final FakeLineageService svc = new FakeLineageService();
        svc.know( "Felis", lineage( "class", "Mammalia", "order", "Carnivora", "genus", "Felis" ) );
        svc.know( "Canis", lineage( "class", "Mammalia", "order", "Carnivora", "genus", "Canis" ) );
        final Phylogeny tree = mammalTree();
        try {
            for( final String name : TreePanelUtil.unresolvedTipTaxa( tree, "order", svc ) ) {
                svc.fetch( name );
            }
        }
        catch ( final Exception e ) {
            return fail( "fake fetch must not throw: " + e );
        }
        final java.util.List<CladeBand> bands = TreePanelUtil.cladeBands( tree, "order", svc );
        // the same three maximal clades the colorizer finds: the Rodentia clade + Felis + Canis
        if ( bands.size() != 3 ) {
            return fail( "expected 3 clade bands (Rodentia, Felis, Canis), got " + bands.size() );
        }
        final java.util.List<Color> carnivora_colors = new java.util.ArrayList<Color>();
        Color rodentia_color = null;
        for( final CladeBand b : bands ) {
            if ( ( b.getColor() == null ) || ( b.getRoot() == null ) ) {
                return fail( "every band must carry a color and a clade root" );
            }
            if ( "Carnivora".equals( b.getTaxon() ) ) {
                carnivora_colors.add( b.getColor() );
            }
            else if ( "Rodentia".equals( b.getTaxon() ) ) {
                rodentia_color = b.getColor();
            }
            else {
                return fail( "unexpected band taxon: " + b.getTaxon() );
            }
        }
        if ( ( carnivora_colors.size() != 2 ) || ( rodentia_color == null ) ) {
            return fail( "paraphyletic Carnivora must yield 2 bands + 1 Rodentia band" );
        }
        if ( !carnivora_colors.get( 0 ).equals( carnivora_colors.get( 1 ) ) ) {
            return fail( "the two Carnivora bands (same taxon) must share a color" );
        }
        if ( carnivora_colors.get( 0 ).equals( rodentia_color ) ) {
            return fail( "Carnivora and Rodentia bands must have distinct colors" );
        }
        // The branch colorizer ("Colorize Subtrees via Taxonomic Rank") handles the SAME polyphyly /
        // gene-duplication shape: all 3 monophyletic clades (2 separate Carnivora + 1 Rodentia) are
        // colorized, but they collapse to just 2 colors (one per taxon NAME), so the two Carnivora
        // clades match. This is the property that matters for large gene trees.
        final java.util.Map<String, Color> branch_legend = new java.util.HashMap<String, Color>();
        final int colorizations = TreePanelUtil.colorPhylogenyAccordingToRanks( mammalTree(), "order", svc,
                                                                                branch_legend );
        if ( colorizations != 3 ) {
            return fail( "branch colorizer: expected 3 colorized clades (2 Carnivora + 1 Rodentia), got "
                    + colorizations );
        }
        if ( branch_legend.size() != 2 ) {
            return fail( "branch colorizer: 3 polyphyletic clades must collapse to 2 colors (by taxon name), got "
                    + branch_legend.size() );
        }
        // a user color override (taxon -> color) replaces the auto-assigned color for that taxon only
        final Color override = new Color( 12, 34, 56 );
        final java.util.Map<String, Color> overrides = new java.util.HashMap<String, Color>();
        overrides.put( "Carnivora", override );
        int overridden = 0;
        for( final CladeBand b : TreePanelUtil.cladeBands( tree, "order", svc, overrides ) ) {
            if ( "Carnivora".equals( b.getTaxon() ) ) {
                if ( !override.equals( b.getColor() ) ) {
                    return fail( "Carnivora band did not pick up the color override" );
                }
                ++overridden;
            }
            else if ( override.equals( b.getColor() ) ) {
                return fail( "override leaked onto a non-overridden band" );
            }
        }
        if ( overridden != 2 ) {
            return fail( "both Carnivora bands must use the override color" );
        }
        // the branch colorizer honors the same override (it surfaces in the legend it fills)
        final java.util.Map<String, Color> legend = new java.util.HashMap<String, Color>();
        TreePanelUtil.colorPhylogenyAccordingToRanks( mammalTree(), "order", svc, legend, overrides );
        if ( !override.equals( legend.get( "Carnivora" ) ) ) {
            return fail( "colorPhylogenyAccordingToRanks did not apply the override to its legend" );
        }
        // degenerate inputs yield no bands (never throw)
        if ( !TreePanelUtil.cladeBands( tree, "", svc ).isEmpty()
                || !TreePanelUtil.cladeBands( null, "order", svc ).isEmpty() ) {
            return fail( "empty rank / null tree must yield no bands" );
        }
        return true;
    }

    /** Builds a LinkedHashMap of rank-&gt;name pairs from a flat rank1,name1,rank2,name2,... list. */
    private static Map<String, String> lineage( final String... rank_name_pairs ) {
        final Map<String, String> m = new LinkedHashMap<String, String>();
        for( int i = 0; ( i + 1 ) < rank_name_pairs.length; i += 2 ) {
            m.put( rank_name_pairs[ i ], rank_name_pairs[ i + 1 ] );
        }
        return m;
    }

    /**
     * root
     *  +-- A
     *  |    +-- Rodentia (internal, rank=order)
     *  |    |    +-- Mus     (genus)
     *  |    |    +-- Rattus  (genus)
     *  |    +-- Felis        (genus; resolves to order Carnivora only via the DB)
     *  +-- B
     *       +-- Canis        (genus; order Carnivora via the DB) -- not a clade with Felis (paraphyly)
     *       +-- x_unknown    (scientific name the DB cannot resolve)
     */
    private static Phylogeny mammalTree() {
        final PhylogenyNode rodentia = internalOrder( "Rodentia" );
        rodentia.addAsChild( genusLeaf( "Mus" ) );
        rodentia.addAsChild( genusLeaf( "Rattus" ) );
        final PhylogenyNode a = new PhylogenyNode();
        a.addAsChild( rodentia );
        a.addAsChild( genusLeaf( "Felis" ) );
        final PhylogenyNode b = new PhylogenyNode();
        b.addAsChild( genusLeaf( "Canis" ) );
        b.addAsChild( namedLeaf( "x_unknown", "Nonexistus" ) ); // queryable, but unknown to the DB
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( a );
        root.addAsChild( b );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode internalOrder( final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        try {
            t.setRank( "order" );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode genusLeaf( final String sci ) {
        final PhylogenyNode n = namedLeaf( sci, sci );
        try {
            n.getNodeData().getTaxonomy().setRank( "genus" );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        return n;
    }

    private static PhylogenyNode namedLeaf( final String node_name, final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( node_name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode findLeaf( final Phylogeny tree, final String node_name ) {
        for( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = tree
                .iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( node_name.equals( n.getName() ) ) {
                return n;
            }
        }
        return null;
    }

    /** In-memory {@link TaxonomicLineageService}: {@code lineageOf} is cache-only; {@code fetch} copies from the "DB". */
    private static final class FakeLineageService implements TaxonomicLineageService {

        private final Map<String, RankedLineage> _db    = new HashMap<String, RankedLineage>();
        private final Map<String, RankedLineage> _cache = new HashMap<String, RankedLineage>();

        void know( final String name, final Map<String, String> rank_to_name ) {
            _db.put( name.toLowerCase( Locale.ROOT ), new RankedLineage( rank_to_name ) );
        }

        @Override
        public RankedLineage lineageOf( final String taxon ) {
            return ( taxon == null ) ? null : _cache.get( taxon.toLowerCase( Locale.ROOT ) );
        }

        @Override
        public RankedLineage fetch( final String taxon ) {
            final String k = taxon.toLowerCase( Locale.ROOT );
            final RankedLineage rl = _db.containsKey( k ) ? _db.get( k ) : RankedLineage.EMPTY;
            _cache.put( k, rl );
            return rl;
        }
    }

    /**
     * taxonomyLabel picks scientific name, else common name, else taxonomy code -- so a rank
     * colorization legend gets a row even for taxa identified only by a common name or code (not just
     * a scientific name).
     */
    private static boolean testTaxonomyLabel() {
        final Taxonomy sci_and_common = new Taxonomy();
        sci_and_common.setScientificName( "Diptera" );
        sci_and_common.setCommonName( "flies" );
        if ( !"Diptera".equals( TreePanelUtil.taxonomyLabel( sci_and_common ) ) ) {
            return fail( "scientific name must win over common name" );
        }
        final Taxonomy common_only = new Taxonomy();
        common_only.setCommonName( "flies" );
        if ( !"flies".equals( TreePanelUtil.taxonomyLabel( common_only ) ) ) {
            return fail( "common name must be used when no scientific name (the legend-empty bug fix)" );
        }
        if ( !"".equals( TreePanelUtil.taxonomyLabel( new Taxonomy() ) )
                || !"".equals( TreePanelUtil.taxonomyLabel( null ) ) ) {
            return fail( "an empty or null taxonomy must yield an empty label" );
        }
        return true;
    }

    /** capEntries bounds a legend to its first N entries (iteration order), preserving keys+colors. */
    private static boolean testCapEntries() {
        if ( !TreePanelUtil.capEntries( null, 5 ).isEmpty() ) {
            return fail( "null input must yield an empty cap" );
        }
        final Map<String, Color> in = new LinkedHashMap<String, Color>();
        in.put( "a", Color.RED );
        in.put( "b", Color.GREEN );
        in.put( "c", Color.BLUE );
        // cap >= size keeps everything, in order
        final Map<String, Color> all = TreePanelUtil.capEntries( in, 5 );
        if ( !new ArrayList<String>( all.keySet() ).equals( new ArrayList<String>( in.keySet() ) ) ) {
            return fail( "cap >= size must keep all entries in iteration order; got " + all.keySet() );
        }
        if ( all.get( "a" ) != Color.RED ) {
            return fail( "cap must preserve the colors" );
        }
        // cap < size keeps exactly the first N entries
        final Map<String, Color> capped = TreePanelUtil.capEntries( in, 2 );
        if ( ( capped.size() != 2 ) || !capped.containsKey( "a" ) || !capped.containsKey( "b" )
                || capped.containsKey( "c" ) ) {
            return fail( "cap < size must keep exactly the first N entries; got " + capped.keySet() );
        }
        // cap of 0 -> empty (the legend's "+N more" footer then accounts for all of them)
        if ( !TreePanelUtil.capEntries( in, 0 ).isEmpty() ) {
            return fail( "cap of 0 must yield an empty map" );
        }
        return true;
    }

    /**
     * The node-symbol support math (see TreePanel.paintNodeSupportSymbol): scale detection picks the
     * absolute ceiling (1 or 100) so symbols mean the same across trees; the support fraction clamps
     * to 0..1; size interpolates min..max by that fraction; the threshold test compares the fraction
     * to the cutoff.
     */
    private static boolean testSupportSymbolMath() {
        // scale ceiling: anything above 1 implies the 0..100 family, otherwise 0..1
        if ( TreePanelUtil.confidenceScaleMaxFor( 0.0 ) != 1.0 ) {
            return fail( "empty/zero support implies the 0..1 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 0.95 ) != 1.0 ) {
            return fail( "max 0.95 (posterior probability) implies the 0..1 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 1.0 ) != 1.0 ) {
            return fail( "max exactly 1.0 stays on the 0..1 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 1.01 ) != 100.0 ) {
            return fail( "any value above 1 implies the 0..100 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 70.0 ) != 100.0 ) {
            return fail( "max 70 (bootstrap) implies the 0..100 scale" );
        }
        // support fraction: clamped to 0..1, scale-relative
        if ( TreePanelUtil.supportFraction( 50.0, 100.0 ) != 0.5 ) {
            return fail( "50 on a 0..100 scale is fraction 0.5" );
        }
        if ( TreePanelUtil.supportFraction( 0.8, 1.0 ) != 0.8 ) {
            return fail( "0.8 on a 0..1 scale is fraction 0.8" );
        }
        if ( TreePanelUtil.supportFraction( 150.0, 100.0 ) != 1.0 ) {
            return fail( "fraction must clamp to 1.0 above the scale ceiling" );
        }
        if ( TreePanelUtil.supportFraction( -1.0, 100.0 ) != 0.0 ) {
            return fail( "negative support clamps to fraction 0.0" );
        }
        if ( TreePanelUtil.supportFraction( 50.0, 0.0 ) != 0.0 ) {
            return fail( "a non-positive scale must not divide; fraction 0.0" );
        }
        // size interpolation: min at 0, max at full support, monotonic
        final float min = 2.0f;
        final float max = 8.0f;
        if ( TreePanelUtil.supportSymbolSize( 0.0, 100.0, min, max ) != min ) {
            return fail( "zero support gives the minimum symbol size" );
        }
        if ( TreePanelUtil.supportSymbolSize( 100.0, 100.0, min, max ) != max ) {
            return fail( "full support gives the maximum symbol size" );
        }
        if ( TreePanelUtil.supportSymbolSize( 50.0, 100.0, min, max ) != 5.0f ) {
            return fail( "half support gives the midpoint symbol size" );
        }
        if ( TreePanelUtil.supportSymbolSize( 30.0, 100.0, min, max ) >= TreePanelUtil
                .supportSymbolSize( 60.0, 100.0, min, max ) ) {
            return fail( "symbol size must grow with support" );
        }
        // threshold test (cutoff is a fraction of the scale)
        if ( !TreePanelUtil.isSupportAtOrAboveThreshold( 95.0, 100.0, 0.95 ) ) {
            return fail( "95/100 must meet the 0.95 cutoff" );
        }
        if ( TreePanelUtil.isSupportAtOrAboveThreshold( 94.0, 100.0, 0.95 ) ) {
            return fail( "94/100 must fall below the 0.95 cutoff" );
        }
        if ( !TreePanelUtil.isSupportAtOrAboveThreshold( 0.96, 1.0, 0.95 ) ) {
            return fail( "0.96 on a 0..1 scale must meet the 0.95 cutoff (scale-independent)" );
        }
        // symbol center: the middle of the branch (support is a branch property), not the node
        // rectangular: x is the branch midpoint; y is the node's y (horizontal branch)
        final float[] rect = TreePanelUtil.supportSymbolCenter( 10f, 30f, 5f, 20f, false );
        if ( ( rect[ 0 ] != 20f ) || ( rect[ 1 ] != 20f ) ) {
            return fail( "rectangular support symbol must sit at branch-midpoint x and node y; got " + rect[ 0 ] + ","
                    + rect[ 1 ] );
        }
        // radial (unrooted/circular): the branch is slanted, so y is the segment midpoint too
        final float[] radial = TreePanelUtil.supportSymbolCenter( 10f, 30f, 5f, 25f, true );
        if ( ( radial[ 0 ] != 20f ) || ( radial[ 1 ] != 15f ) ) {
            return fail( "radial support symbol must sit at the 2-D branch midpoint; got " + radial[ 0 ] + ","
                    + radial[ 1 ] );
        }
        return true;
    }

    /** detectConfidenceScaleMax scans only internal-node confidences and infers the absolute scale. */
    private static boolean testDetectConfidenceScaleMax() {
        // a bootstrap tree (support 90 at the internal node) -> 0..100 scale
        if ( TreePanelUtil.detectConfidenceScaleMax( treeWithInternalConfidence( 90.0, "bootstrap" ) ) != 100.0 ) {
            return fail( "a tree with bootstrap support 90 must be detected as the 0..100 scale" );
        }
        // a Bayesian tree (posterior probability 0.99) -> 0..1 scale
        if ( TreePanelUtil.detectConfidenceScaleMax( treeWithInternalConfidence( 0.99, "posterior" ) ) != 1.0 ) {
            return fail( "a tree with posterior probability 0.99 must be detected as the 0..1 scale" );
        }
        // no confidences anywhere -> defaults to the 0..1 scale (harmless; nothing is drawn)
        final Phylogeny bare = new Phylogeny();
        final PhylogenyNode r = new PhylogenyNode();
        r.addAsChild( new PhylogenyNode() );
        r.addAsChild( new PhylogenyNode() );
        bare.setRoot( r );
        bare.externalNodesHaveChanged();
        if ( TreePanelUtil.detectConfidenceScaleMax( bare ) != 1.0 ) {
            return fail( "a tree without confidences defaults to the 0..1 scale" );
        }
        return true;
    }

    /** (internal:conf, leaf, leaf) under a root -- the confidence sits on an internal branch. */
    private static Phylogeny treeWithInternalConfidence( final double value, final String type ) {
        final Phylogeny p = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode internal = new PhylogenyNode();
        internal.getBranchData().addConfidence( new Confidence( value, type ) );
        internal.addAsChild( new PhylogenyNode() );
        internal.addAsChild( new PhylogenyNode() );
        root.addAsChild( internal );
        root.addAsChild( new PhylogenyNode() );
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    /**
     * The y-distance returned for a given label height must (a) space adjacent leaf rows
     * (2 * y-distance apart) at least one label-height apart so labels do not overlap, and
     * (b) drive the dynamic-hiding factor -- the same formula TreePanel.calcDynamicHidingFactor
     * uses -- down to <= 1 so the "Dyna Hide" indicator clears.
     */
    private static boolean testYDistanceToAvoidLabelOverlap() {
        final int[] heights = { 2, 8, 10, 11, 12, 14, 16, 20, 27, 40 };
        float previous = -1.0f;
        for( final int h : heights ) {
            final float y_dist = TreePanelUtil.yDistanceToAvoidLabelOverlap( h );
            if ( y_dist <= 0.0f ) {
                return fail( "y-distance must be positive for height " + h + " (got " + y_dist + ")" );
            }
            // (a) leaf rows are 2 * y-distance apart; that must be >= the label height
            if ( ( 2.0f * y_dist ) < h ) {
                return fail( "labels would overlap at height " + h + ": spacing " + ( 2.0f * y_dist ) + " < " + h );
            }
            // (b) same formula as TreePanel.calcDynamicHidingFactor: round( h / (1.5 * y-distance) )
            final int hiding_factor = (int) ( 0.5 + ( h / ( 1.5 * y_dist ) ) );
            if ( hiding_factor > 1 ) {
                return fail( "dynamic-hiding factor should be <= 1 at height " + h + " (got " + hiding_factor + ")" );
            }
            // monotonic in the label height (taller labels never need less spacing)
            if ( y_dist < previous ) {
                return fail( "y-distance should grow with label height; broke at height " + h );
            }
            previous = y_dist;
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [TreePanelUtilTest] " + message );
        return false;
    }

    private TreePanelUtilTest() {
        // not instantiable
    }
}
