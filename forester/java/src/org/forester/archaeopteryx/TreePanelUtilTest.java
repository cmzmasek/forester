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
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.List;
import java.util.SortedSet;
import java.util.function.Function;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.archaeopteryx.tools.AncestralTaxonomyInferrer;
import org.forester.phylogeny.data.Identifier;
import org.forester.ws.seqdb.AccessionAwareLineageService;
import org.forester.ws.seqdb.Organism;
import org.forester.ws.seqdb.OrganismSource;
import org.forester.ws.seqdb.TaxonLineage;
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
                && testRankTaxonCounts() && testTaxonomyLabel() && testRankColorization() && testTipQueryName()
                && testTipIdentityWins() && testRankColorizationTaxIdKeying() && testInternalTaxonGapFill()
                && testRankResolutionIdThenNameCache() && testNamelessTaxonomyTipResolvesByNodeName()
                && testInferenceFeedsRankAssignment() && testWriteCladeTaxonomies() && testInternalTaxaByRank()
                && testCladeBands() && testRankColorizationViaSequenceIds() && testInternalLabelAboveBranchLayout()
                && testAbbreviateScientificName() && testSupportColor() && testBreakLongBranches() && testScaleGridLines()
                && testScaleAxisTickValues() && testCalendarTickYears() && testMaAxisTicks() && testFormatCompactNumber() && testHpdBarXRange() && testSpindleHalfHeight()
                && testOrientationTransform() && testInternalLabelAlignWidth() && testAutoTipLabelDirection()
                && testUserVisiblePropertiesText() && testTipLineagesAndUnresolved() && testInferenceStrings()
                && testIsDuplicateOfAncestorTaxon() && testScaleAxisFloating() && testDomainBoxHeight()
                && testAncestralPieData();
    }

    /**
     * The ancestral-state pie data layer: parseBraceList (top-level commas, quotes, braces), ancestralStateTraits
     * (detect beast:<trait>_set + _set_prob), stateDistribution (zip + normalize; tip single state; length mismatch
     * / malformed -> empty), and collectAncestralStates.
     */
    private static boolean testAncestralPieData() {
        final List<String> l = TreePanelUtil.parseBraceList( "{Africa,Europe,Asia}" );
        if ( ( l.size() != 3 ) || !l.get( 0 ).equals( "Africa" ) || !l.get( 2 ).equals( "Asia" ) ) {
            return fail( "parseBraceList must split {a,b,c} into 3 elements, got " + l );
        }
        final List<String> q = TreePanelUtil.parseBraceList( "{ \"x,y\" , z }" );
        if ( ( q.size() != 2 ) || !q.get( 0 ).equals( "x,y" ) || !q.get( 1 ).equals( "z" ) ) {
            return fail( "parseBraceList must keep a quoted comma and trim, got " + q );
        }
        if ( !TreePanelUtil.parseBraceList( "{}" ).isEmpty() || !TreePanelUtil.parseBraceList( "" ).isEmpty() ) {
            return fail( "parseBraceList must be empty for {} / empty input" );
        }
        // root(location distribution) -> { tipA(location=Africa), tipB(location=Europe) }
        final PhylogenyNode root = new PhylogenyNode();
        withProp( root, "beast:location_set", "{Africa,Europe,Asia}" );
        withProp( root, "beast:location_set_prob", "{0.5,0.3,0.2}" );
        final PhylogenyNode tipA = new PhylogenyNode();
        tipA.setName( "A" );
        withProp( tipA, "beast:location", "Africa" );
        final PhylogenyNode tipB = new PhylogenyNode();
        tipB.setName( "B" );
        withProp( tipB, "beast:location", "Europe" );
        root.addAsChild( tipA );
        root.addAsChild( tipB );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        if ( ( TreePanelUtil.ancestralStateTraits( phy ).size() != 1 )
                || !TreePanelUtil.ancestralStateTraits( phy ).contains( "location" ) ) {
            return fail( "ancestralStateTraits must detect the 'location' distribution trait" );
        }
        final List<TreePanelUtil.StateProbability> d = TreePanelUtil.stateDistribution( root, "location" );
        double sum = 0.0;
        for( final TreePanelUtil.StateProbability sp : d ) {
            sum += sp.getProbability();
        }
        if ( ( d.size() != 3 ) || ( Math.abs( sum - 1.0 ) > 1e-9 ) || !d.get( 0 ).getState().equals( "Africa" )
                || ( Math.abs( d.get( 0 ).getProbability() - 0.5 ) > 1e-9 ) ) {
            return fail( "internal-node distribution must be 3 states summing to 1 with Africa=0.5" );
        }
        final List<TreePanelUtil.StateProbability> t = TreePanelUtil.stateDistribution( tipA, "location" );
        if ( ( t.size() != 1 ) || !t.get( 0 ).getState().equals( "Africa" ) || ( t.get( 0 ).getProbability() != 1.0 ) ) {
            return fail( "a tip's single observed state must be one full-probability wedge" );
        }
        // Auspice / Nextstrain traits (nextstrain: prefix) flow through the SAME pie layer as BEAST's
        final PhylogenyNode ns_root = new PhylogenyNode();
        withProp( ns_root, "nextstrain:region_set", "{Asia,Europe}" );
        withProp( ns_root, "nextstrain:region_set_prob", "{0.7,0.3}" );
        final PhylogenyNode ns_tip = new PhylogenyNode();
        ns_tip.setName( "ns" );
        withProp( ns_tip, "nextstrain:region", "Asia" );
        ns_root.addAsChild( ns_tip );
        final Phylogeny ns_phy = new Phylogeny();
        ns_phy.setRoot( ns_root );
        ns_phy.externalNodesHaveChanged();
        if ( !TreePanelUtil.ancestralStateTraits( ns_phy ).contains( "region" ) ) {
            return fail( "ancestralStateTraits must detect a nextstrain:<trait> distribution too" );
        }
        final List<TreePanelUtil.StateProbability> nd = TreePanelUtil.stateDistribution( ns_root, "region" );
        if ( ( nd.size() != 2 ) || !nd.get( 0 ).getState().equals( "Asia" )
                || ( Math.abs( nd.get( 0 ).getProbability() - 0.7 ) > 1e-9 ) ) {
            return fail( "a nextstrain: distribution must be 2 states with Asia=0.7, got " + nd );
        }
        if ( TreePanelUtil.stateDistribution( ns_tip, "region" ).size() != 1 ) {
            return fail( "a nextstrain: tip's single observed state must be one wedge" );
        }
        // unnormalized probabilities normalize
        final PhylogenyNode n2 = new PhylogenyNode();
        withProp( n2, "beast:host_set", "{Bat,Human}" );
        withProp( n2, "beast:host_set_prob", "{3,1}" );
        final List<TreePanelUtil.StateProbability> d2 = TreePanelUtil.stateDistribution( n2, "host" );
        if ( ( d2.size() != 2 ) || ( Math.abs( d2.get( 0 ).getProbability() - 0.75 ) > 1e-9 ) ) {
            return fail( "unnormalized {3,1} must normalize to {0.75,0.25}" );
        }
        // set/prob length mismatch -> no distribution
        final PhylogenyNode n3 = new PhylogenyNode();
        withProp( n3, "beast:x_set", "{A,B}" );
        withProp( n3, "beast:x_set_prob", "{0.5}" );
        if ( !TreePanelUtil.stateDistribution( n3, "x" ).isEmpty() ) {
            return fail( "a set/prob length mismatch must yield no distribution" );
        }
        // malformed probability -> no distribution
        final PhylogenyNode n4 = new PhylogenyNode();
        withProp( n4, "beast:y_set", "{A,B}" );
        withProp( n4, "beast:y_set_prob", "{0.5,bad}" );
        if ( !TreePanelUtil.stateDistribution( n4, "y" ).isEmpty() ) {
            return fail( "a malformed probability must yield no distribution" );
        }
        final SortedSet<String> all = TreePanelUtil.collectAncestralStates( phy, "location" );
        if ( ( all.size() != 3 ) || !all.contains( "Africa" ) || !all.contains( "Asia" ) || !all.contains( "Europe" ) ) {
            return fail( "collectAncestralStates must gather all distinct states, got " + all );
        }
        // parseBraceList also strips [..] brackets and splits only TOP-LEVEL commas (a nested list stays one element)
        final List<String> br = TreePanelUtil.parseBraceList( "[a,b]" );
        if ( ( br.size() != 2 ) || !br.get( 1 ).equals( "b" ) ) {
            return fail( "parseBraceList must strip [..] brackets, got " + br );
        }
        final List<String> nested = TreePanelUtil.parseBraceList( "{a,{b,c},d}" );
        if ( ( nested.size() != 3 ) || !nested.get( 1 ).equals( "{b,c}" ) ) {
            return fail( "parseBraceList must split only top-level commas (nested stays one element), got " + nested );
        }
        // an apostrophe in a state name must NOT swallow the rest (real place names: Xi'an, Cote d'Ivoire)
        final List<String> apos = TreePanelUtil.parseBraceList( "{Xi'an,Beijing}" );
        if ( ( apos.size() != 2 ) || !apos.get( 0 ).equals( "Xi'an" ) || !apos.get( 1 ).equals( "Beijing" ) ) {
            return fail( "parseBraceList must treat an apostrophe as literal, got " + apos );
        }
        // grayShade (B&W wedge/legend ramp): distinct grays for distinct indices; a single element is mid-gray
        if ( TreePanelUtil.grayShade( 0, 3 ).equals( TreePanelUtil.grayShade( 2, 3 ) ) ) {
            return fail( "grayShade must give distinct grays for distinct wedge indices" );
        }
        if ( TreePanelUtil.grayShade( 0, 1 ).getRed() != 128 ) {
            return fail( "grayShade of a single element must be a neutral mid-gray" );
        }
        // a NEGATIVE probability is rejected like NaN (no plausible-but-wrong pie from corrupt data)
        final PhylogenyNode n5 = new PhylogenyNode();
        withProp( n5, "beast:z_set", "{A,B,C}" );
        withProp( n5, "beast:z_set_prob", "{0.9,-0.4,0.5}" );
        if ( !TreePanelUtil.stateDistribution( n5, "z" ).isEmpty() ) {
            return fail( "a negative probability must yield no distribution" );
        }
        // a MALFORMED distribution must NOT fall through to the bare single-state disc (which would overstate an
        // uncertain internal node as a confident 100%): a length-mismatched set + a bare beast:trait -> NO pie
        final PhylogenyNode n6 = new PhylogenyNode();
        withProp( n6, "beast:w_set", "{A,B}" );
        withProp( n6, "beast:w_set_prob", "{0.5}" );
        withProp( n6, "beast:w", "A" );
        if ( !TreePanelUtil.stateDistribution( n6, "w" ).isEmpty() ) {
            return fail( "a malformed distribution must not fall through to a misleading single-state disc" );
        }
        return true;
    }

    private static void withProp( final PhylogenyNode node, final String ref, final String value ) {
        PropertiesList pl = node.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            node.getNodeData().setProperties( pl );
        }
        pl.addProperty( new Property( ref, value, "", "xsd:string", AppliesTo.NODE ) );
    }

    /**
     * The scale axis FLOATS at the viewport bottom / breadth edge on screen (so it never scrolls out of view when
     * zoomed), but stays anchored to the tree/export extent for a file export (WYSIWYG) and to the whole canvas for
     * the direct File&gt;Print path.
     */
    private static boolean testScaleAxisFloating() {
        // horizontal: screen -> viewport bottom (ignores the taller canvas); export -> extent; print(height 0) -> canvas
        if ( TreePanelUtil.scaleAxisFloatingBottom( false, false, 0, 0, 1000, 300 ) != 300 ) {
            return fail( "screen axis must float at the viewport bottom (300), not the canvas bottom" );
        }
        if ( TreePanelUtil.scaleAxisFloatingBottom( false, true, 10, 500, 1000, 300 ) != 510 ) {
            return fail( "a file export must anchor the axis to the export extent (graphics_file_y + height = 510)" );
        }
        if ( TreePanelUtil.scaleAxisFloatingBottom( true, false, 0, 400, 1000, 300 ) != 400 ) {
            return fail( "a PDF export with a real height must anchor to that extent (400)" );
        }
        if ( TreePanelUtil.scaleAxisFloatingBottom( true, false, 0, 0, 1000, 300 ) != 1000 ) {
            return fail( "File>Print (export flag, height 0) must anchor to the whole canvas (1000)" );
        }
        // vertical ruler: screen -> viewport edge on the ruler's side (opposite the tree); export -> tree-anchored x
        if ( TreePanelUtil.scaleAxisRulerX( false, false, 700, 1, 50, 400 ) != 51 ) {
            return fail( "screen ruler with the tree to the right must float to the LEFT viewport edge (viewport_x+1)" );
        }
        if ( TreePanelUtil.scaleAxisRulerX( false, false, 700, -1, 50, 400 ) != 449 ) {
            return fail( "screen ruler with the tree to the left must float to the RIGHT viewport edge (x+width-1)" );
        }
        if ( ( TreePanelUtil.scaleAxisRulerX( false, true, 700, 1, 50, 400 ) != 700 )
                || ( TreePanelUtil.scaleAxisRulerX( true, false, 700, -1, 50, 400 ) != 700 ) ) {
            return fail( "an export must keep the tree-anchored breadth position (700)" );
        }
        return true;
    }

    /**
     * isDuplicateOfAncestorTaxon flags a non-root internal node whose taxon matches its nearest annotated ancestor
     * (so the redundant label can be suppressed at draw time), skipping un-annotated intermediates; tips, the root,
     * and un-annotated nodes never qualify.
     */
    private static boolean testIsDuplicateOfAncestorTaxon() {
        // the label the paint path would show; here = scientific name (what the default checkboxes render)
        final Function<PhylogenyNode, String> sci = n -> TreePanelUtil.taxonomyLabel( n.getNodeData().getTaxonomy() );
        // root(Carnivora) -> mid(Carnivora) -> {Felis, Panthera}; root -> Canis
        final PhylogenyNode mid = taxInternal( "Carnivora" );
        mid.addAsChild( genusLeaf( "Felis" ) );
        mid.addAsChild( genusLeaf( "Panthera" ) );
        final PhylogenyNode root = taxInternal( "Carnivora" );
        root.addAsChild( mid );
        root.addAsChild( genusLeaf( "Canis" ) );
        if ( TreePanelUtil.isDuplicateOfAncestorTaxon( root, sci ) ) {
            return fail( "the root can never be a duplicate" );
        }
        if ( !TreePanelUtil.isDuplicateOfAncestorTaxon( mid, sci ) ) {
            return fail( "a Carnivora node inside a Carnivora node must be a duplicate" );
        }
        if ( TreePanelUtil.isDuplicateOfAncestorTaxon( mid.getChildNode( 0 ), sci ) ) {
            return fail( "an external tip is never a duplicate" );
        }
        // a DIFFERENT taxon than the ancestor is not a duplicate
        final PhylogenyNode fam = taxInternal( "Felidae" );
        fam.addAsChild( genusLeaf( "Lynx" ) );
        final PhylogenyNode root2 = taxInternal( "Carnivora" );
        root2.addAsChild( fam );
        if ( TreePanelUtil.isDuplicateOfAncestorTaxon( fam, sci ) ) {
            return fail( "Felidae inside Carnivora is NOT a duplicate" );
        }
        // skip an un-annotated intermediate: A(Carnivora) -> B(no taxonomy) -> C(Carnivora) duplicates A
        final PhylogenyNode c = taxInternal( "Carnivora" );
        c.addAsChild( genusLeaf( "Felis" ) );
        final PhylogenyNode b = new PhylogenyNode();
        b.addAsChild( c );
        final PhylogenyNode a = taxInternal( "Carnivora" );
        a.addAsChild( b );
        if ( TreePanelUtil.isDuplicateOfAncestorTaxon( b, sci ) ) {
            return fail( "an internal node without a taxonomy is never a duplicate" );
        }
        if ( !TreePanelUtil.isDuplicateOfAncestorTaxon( c, sci ) ) {
            return fail( "C(Carnivora) must duplicate A(Carnivora) across an un-annotated B" );
        }
        // ... and NOT when the nearest annotated ancestor differs
        final PhylogenyNode c2 = taxInternal( "Rodentia" );
        c2.addAsChild( genusLeaf( "Mus" ) );
        final PhylogenyNode b2 = new PhylogenyNode();
        b2.addAsChild( c2 );
        final PhylogenyNode a2 = taxInternal( "Mammalia" );
        a2.addAsChild( b2 );
        if ( TreePanelUtil.isDuplicateOfAncestorTaxon( c2, sci ) ) {
            return fail( "Rodentia inside Mammalia is NOT a duplicate" );
        }
        // the label BASIS is the caller's: two nodes sharing a scientific name but showing a different CODE are
        // duplicates under the scientific-name labeler yet NOT under a code-showing labeler (the checkbox-aware fix)
        final Function<PhylogenyNode, String> code = n -> {
            final String cd = n.getNodeData().getTaxonomy().getTaxonomyCode();
            return ( cd == null ) ? "" : cd;
        };
        final PhylogenyNode cc = taxWithCode( "Bacteria", "BBBBB" );
        cc.addAsChild( genusLeaf( "x" ) );
        final PhylogenyNode aa = taxWithCode( "Bacteria", "AAAAA" );
        aa.addAsChild( cc );
        if ( !TreePanelUtil.isDuplicateOfAncestorTaxon( cc, sci ) ) {
            return fail( "same scientific name IS a duplicate under the scientific-name labeler" );
        }
        if ( TreePanelUtil.isDuplicateOfAncestorTaxon( cc, code ) ) {
            return fail( "differing codes are NOT a duplicate under a code-showing labeler" );
        }
        // an ancestor whose label renders EMPTY (only a tax-id) is skipped: A(Carnivora) -> B(empty) -> C(Carnivora)
        final PhylogenyNode c3 = taxInternal( "Carnivora" );
        c3.addAsChild( genusLeaf( "Felis" ) );
        final PhylogenyNode b3 = taxIdOnly();
        b3.addAsChild( c3 );
        final PhylogenyNode a3 = taxInternal( "Carnivora" );
        a3.addAsChild( b3 );
        if ( !TreePanelUtil.isDuplicateOfAncestorTaxon( c3, sci ) ) {
            return fail( "C(Carnivora) must duplicate A(Carnivora) across an empty-label B" );
        }
        return true;
    }

    /** An internal node carrying a taxonomy with just a scientific name (no rank). */
    private static PhylogenyNode taxInternal( final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode taxWithCode( final String sci, final String code ) {
        final PhylogenyNode n = taxInternal( sci );
        try {
            n.getNodeData().getTaxonomy().setTaxonomyCode( code );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        return n;
    }

    /** An internal node with a Taxonomy that renders NO visible label (only a tax-id). */
    private static PhylogenyNode taxIdOnly() {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setIdentifier( new Identifier( "9999", "ncbi" ) );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    /**
     * tipLineages resolves each tip from the cache (preferred) or a reconstruction of its stored taxonomy, and
     * tipsWithoutLineage flags only the tips that carry no usable lineage in-tree and are not yet cached (the set
     * to fetch online for ancestral-taxonomy inference).
     */
    private static boolean testTipLineagesAndUnresolved() {
        final FakeLineageService svc = new FakeLineageService();
        final PhylogenyNode felis = genusLeaf( "Felis" ); // resolved via the cache after a fetch
        svc.know( "Felis", lineage( "order", "Carnivora", "genus", "Felis" ) );
        svc.fetch( "Felis" );
        final PhylogenyNode mus = genusLeaf( "Mus" ); // carries an in-memory getLineage() + an NCBI id (prior Fetch)
        mus.getNodeData().getTaxonomy().setLineage( namesOf( "Eukaryota", "Rodentia", "Mus" ) );
        mus.getNodeData().getTaxonomy().setIdentifier( new Identifier( "10090", "ncbi" ) );
        final PhylogenyNode bos = genusLeaf( "Bos" ); // stored lineage but a NON-ncbi identifier (must be ignored)
        bos.getNodeData().getTaxonomy().setLineage( namesOf( "Eukaryota", "Laurasiatheria", "Bos" ) );
        bos.getNodeData().getTaxonomy().setIdentifier( new Identifier( "P0", "uniprot" ) );
        final PhylogenyNode solo = genusLeaf( "Solo" ); // a stored lineage that is ONLY the tip's own taxon (no ancestor)
        solo.getNodeData().getTaxonomy().setLineage( namesOf( "Solo" ) );
        final PhylogenyNode plain = new PhylogenyNode(); // no taxonomy at all
        plain.setName( "iso_1" );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( felis );
        root.addAsChild( mus );
        root.addAsChild( bos );
        root.addAsChild( solo );
        root.addAsChild( plain );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        // needs an online fetch: the plain tip AND 'Solo' (whose only stored lineage entry is its own taxon, so it
        // has no ancestor to intersect) -- but NOT Felis (cached) or Mus/Bos (real stored ancestors)
        final SortedSet<String> unresolved = TreePanelUtil.tipsWithoutLineage( phy, svc );
        if ( ( unresolved.size() != 2 ) || !unresolved.contains( "iso_1" ) || !unresolved.contains( "Solo" ) ) {
            return fail( "tipsWithoutLineage should be {Solo, iso_1}, got " + unresolved );
        }
        final Map<PhylogenyNode, TaxonLineage> map = TreePanelUtil.tipLineages( phy, svc );
        if ( !"Carnivora".equals( map.get( felis ).at( "order" ) ) ) {
            return fail( "Felis lineage should come from the cache (order Carnivora)" );
        }
        final TaxonLineage mus_l = map.get( mus );
        if ( mus_l.isEmpty() || !mus_l.lineageNames().contains( "Rodentia" )
                || !"Mus".equals( mus_l.getScientificName() ) ) {
            return fail( "Mus lineage should be reconstructed from its stored taxonomy (ancestors incl. Rodentia)" );
        }
        if ( !"10090".equals( mus_l.getTaxId() ) ) {
            return fail( "Mus reconstruction should carry its NCBI tax-id 10090, got " + mus_l.getTaxId() );
        }
        if ( map.get( bos ).getTaxId() != null ) {
            return fail( "a non-'ncbi' identifier must be ignored (Bos tax-id should be null)" );
        }
        if ( !map.get( plain ).isEmpty() ) {
            return fail( "a tip with no taxonomy should map to an EMPTY lineage" );
        }
        return true;
    }

    private static ArrayList<String> namesOf( final String... names ) {
        final ArrayList<String> l = new ArrayList<String>();
        for( final String n : names ) {
            l.add( n );
        }
        return l;
    }

    /** The pure inference provenance + completion-summary strings are grammatically singular/plural-correct. */
    private static boolean testInferenceStrings() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( genusLeaf( "A" ) );
        root.addAsChild( genusLeaf( "B" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setName( "mytree" );
        phy.externalNodesHaveChanged();
        final String p1 = AncestralTaxonomyInferrer.inferenceProvenance( phy, 1, false );
        if ( !p1.contains( "1 internal node in" ) || p1.contains( "overwriting" ) || !p1.contains( "\"mytree\"" )
                || !p1.contains( "2 tips" ) ) {
            return fail( "singular no-overwrite provenance is wrong: " + p1 );
        }
        final String p3 = AncestralTaxonomyInferrer.inferenceProvenance( phy, 3, true );
        if ( !p3.contains( "3 internal nodes" ) || !p3.contains( "overwriting existing internal taxa" ) ) {
            return fail( "plural overwrite provenance is wrong: " + p3 );
        }
        if ( !"Assigned a taxonomy to 1 internal node.".equals( AncestralTaxonomyInferrer.inferenceSummary( 1, 0, 0 ) ) ) {
            return fail( "singular summary is wrong" );
        }
        final String s = AncestralTaxonomyInferrer.inferenceSummary( 3, 2, 1 );
        if ( !s.contains( "Assigned taxonomies to 3 internal nodes." )
                || !s.contains( "Kept 2 existing internal-node taxonomies" )
                || !s.contains( "1 tip had no ancestor lineage" ) ) {
            return fail( "plural summary is wrong: " + s );
        }
        final String s0 = AncestralTaxonomyInferrer.inferenceSummary( 0, 1, 2 );
        if ( !s0.contains( "Kept 1 existing internal-node taxonomy " ) || !s0.contains( "2 tips had no ancestor" ) ) {
            return fail( "mixed-count summary is wrong: " + s0 );
        }
        return true;
    }

    /** The node-data display text hides internal {@code aptx:*} metadata (the persisted Re-import profile) but keeps
     *  real user properties -- so hovering / inspecting the root of an annotation-imported tree shows no machinery. */
    private static boolean testUserVisiblePropertiesText() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( new Property( "data:host", "cat", "", "xsd:string", AppliesTo.NODE ) );
        props.addProperty( new Property( "aptx:import_profile", "v1;/data/x.csv;0;TIP_NAME;name", "", "xsd:string",
                AppliesTo.NODE ) );
        props.addProperty( new Property( "data:reads", "42", "", "xsd:decimal", AppliesTo.NODE ) );
        final String text = TreePanelUtil.userVisiblePropertiesText( props ).toString();
        if ( text.contains( "aptx:import_profile" ) || text.contains( "v1;" ) ) {
            return fail( "the internal aptx:import_profile property must not surface in the node-data display: " + text );
        }
        if ( !text.contains( "data:host" ) || !text.contains( "cat" ) || !text.contains( "data:reads" )
                || !text.contains( "42" ) ) {
            return fail( "real user properties must still be shown: " + text );
        }
        // the two kept rows are newline-joined, with NO blank line where the internal row was skipped (order is the
        // list's own ref-sorted order, so assert the clean-join shape rather than a fixed order)
        if ( text.startsWith( "\n" ) || text.endsWith( "\n" ) || text.contains( "\n\n" ) ) {
            return fail( "filtered text must not leave a blank line where the internal row was skipped: [" + text + "]" );
        }
        if ( text.split( "\n", -1 ).length != 2 ) {
            return fail( "filtered text should be exactly the two visible rows: [" + text + "]" );
        }
        if ( TreePanelUtil.userVisiblePropertiesText( null ).length() != 0 ) {
            return fail( "a null property list should yield empty text" );
        }
        if ( !TreePanelUtil.isInternalPropertyRef( "aptx:import_profile" )
                || TreePanelUtil.isInternalPropertyRef( "data:host" )
                || TreePanelUtil.isInternalPropertyRef( null ) ) {
            return fail( "isInternalPropertyRef must match only the aptx: namespace" );
        }
        return true;
    }

    /** Auto tip-label angle: upright when labels fit the tip spacing, 45° while their diagonal fits, else vertical. */
    private static boolean testAutoTipLabelDirection() {
        // longest label 40px wide
        if ( TreePanelUtil.autoTipLabelDirection( 50.0, 40.0 ) != Options.TIP_LABEL_DIRECTION.HORIZONTAL ) {
            return fail( "a label narrower than the tip spacing should stay upright (0°)" );
        }
        if ( TreePanelUtil.autoTipLabelDirection( 40.0, 40.0 ) != Options.TIP_LABEL_DIRECTION.HORIZONTAL ) {
            return fail( "a label exactly at the tip spacing should still be upright" );
        }
        // 40px label, spacing 35: too wide for upright, but 40*0.707=28.3 <= 35 -> diagonal
        if ( TreePanelUtil.autoTipLabelDirection( 35.0, 40.0 ) != Options.TIP_LABEL_DIRECTION.ANGLED ) {
            return fail( "a label whose diagonal projection fits should be 45°" );
        }
        // 40px label, spacing 20: even the diagonal (28.3) overflows -> vertical
        if ( TreePanelUtil.autoTipLabelDirection( 20.0, 40.0 ) != Options.TIP_LABEL_DIRECTION.VERTICAL ) {
            return fail( "a dense layout should fall back to vertical (90°)" );
        }
        // degenerate (no layout yet) -> vertical, no divide-by-zero
        if ( ( TreePanelUtil.autoTipLabelDirection( 0.0, 40.0 ) != Options.TIP_LABEL_DIRECTION.VERTICAL )
                || ( TreePanelUtil.autoTipLabelDirection( 50.0, 0.0 ) != Options.TIP_LABEL_DIRECTION.VERTICAL ) ) {
            return fail( "a non-positive spacing/width should fall back to vertical" );
        }
        return true;
    }

    /** The vertical internal-label right-alignment width: the taxonomy's trailing space is dropped only when the
     *  taxonomy is the rightmost element (no node data), so a taxonomy-only label ends flush at the branch. */
    private static boolean testInternalLabelAlignWidth() {
        final int space = 4;
        // taxonomy only -> drop the trailing space so it right-aligns flush at the branch
        if ( TreePanelUtil.internalLabelAlignWidth( 30, 0, space, true, true ) != ( 30 - space ) ) {
            System.out.println( "  internalLabelAlignWidth: taxonomy-only must drop the trailing space" );
            return false;
        }
        // taxonomy + node data -> node data (no trailing space) is rightmost, so keep the full width
        if ( TreePanelUtil.internalLabelAlignWidth( 30, 20, space, true, false ) != 50 ) {
            System.out.println( "  internalLabelAlignWidth: taxonomy+data must keep the full width" );
            return false;
        }
        // node data only (no taxonomy) -> nothing to drop
        if ( TreePanelUtil.internalLabelAlignWidth( 0, 20, space, false, false ) != 20 ) {
            System.out.println( "  internalLabelAlignWidth: data-only must keep the full width" );
            return false;
        }
        // never negative (a taxonomy narrower than one space)
        if ( TreePanelUtil.internalLabelAlignWidth( 2, 0, space, true, true ) != 0 ) {
            System.out.println( "  internalLabelAlignWidth: must clamp to >= 0" );
            return false;
        }
        return true;
    }

    /** The vertical-orientation rotation R: exact corner mappings, a pure rotation (determinant +1, no mirror), all
     *  four logical corners landing in the positive quadrant, and an R^-1 round-trip. Pure math -> headless. */
    private static boolean testOrientationTransform() {
        final double w = 800; // logical depth extent (root x=0 .. tip x=w)
        final double h = 300; // logical breadth/tip-spread extent (y=0 .. y=h)
        if ( !TreePanelUtil.orientationTransformFor( Options.TREE_ORIENTATION.ROOT_LEFT, w, h ).isIdentity() ) {
            return fail( "ROOT_LEFT transform should be the identity" );
        }
        // ROOT_TOP: (x,y) -> (h - y, x)
        final AffineTransform top = TreePanelUtil.orientationTransformFor( Options.TREE_ORIENTATION.ROOT_TOP, w, h );
        if ( !mapsTo( top, 0, 0, h, 0 ) || !mapsTo( top, w, 0, h, w ) || !mapsTo( top, 0, h, 0, 0 )
                || !mapsTo( top, w, h, 0, w ) ) {
            return fail( "ROOT_TOP corner mapping wrong" );
        }
        // ROOT_BOTTOM: (x,y) -> (y, w - x)
        final AffineTransform bot = TreePanelUtil.orientationTransformFor( Options.TREE_ORIENTATION.ROOT_BOTTOM, w, h );
        if ( !mapsTo( bot, 0, 0, 0, w ) || !mapsTo( bot, w, 0, 0, 0 ) || !mapsTo( bot, 0, h, h, w )
                || !mapsTo( bot, w, h, h, 0 ) ) {
            return fail( "ROOT_BOTTOM corner mapping wrong" );
        }
        for ( final AffineTransform r : new AffineTransform[] { top, bot } ) {
            if ( Math.abs( r.getDeterminant() - 1.0 ) > 1.0e-9 ) {
                return fail( "orientation transform must be a pure rotation (determinant +1)" );
            }
            for ( final double[] c : new double[][] { { 0, 0 }, { w, 0 }, { 0, h }, { w, h } } ) {
                final Point2D.Double p = new Point2D.Double( c[ 0 ], c[ 1 ] );
                r.transform( p, p );
                if ( ( p.x < -1.0e-6 ) || ( p.y < -1.0e-6 ) ) {
                    return fail( "orientation transform leaves a corner in a negative quadrant: " + p );
                }
            }
            try {
                final Point2D.Double p = new Point2D.Double( 123.5, 47.25 );
                final Point2D.Double q = new Point2D.Double();
                r.transform( p, q );
                r.createInverse().transform( q, q );
                if ( ( Math.abs( q.x - p.x ) > 1.0e-6 ) || ( Math.abs( q.y - p.y ) > 1.0e-6 ) ) {
                    return fail( "R^-1(R(p)) != p" );
                }
            }
            catch ( final java.awt.geom.NoninvertibleTransformException e ) {
                return fail( "orientation transform should be invertible: " + e );
            }
        }
        return true;
    }

    private static boolean mapsTo( final AffineTransform t, final double x, final double y, final double ex,
                                   final double ey ) {
        final Point2D.Double p = new Point2D.Double( x, y );
        t.transform( p, p );
        return ( Math.abs( p.x - ex ) < 1.0e-6 ) && ( Math.abs( p.y - ey ) < 1.0e-6 );
    }

    /** A node-age HPD bar is anchored to the node's own x and offset by the age deltas (older->left, younger->right);
     *  the bar always straddles the node's x -- independent of tree height / ultrametricity. */
    private static boolean testHpdBarXRange() {
        // node at x=100, age 90, interval [82,98], corr 2: left = 100-(98-90)*2 = 84, right = 100+(90-82)*2 = 116
        final float[] r = TreePanelUtil.hpdBarXRange( 100f, 90.0, 82.0, 98.0, 2.0 );
        if ( ( Math.abs( r[ 0 ] - 84f ) > 0.001f ) || ( Math.abs( r[ 1 ] - 116f ) > 0.001f ) ) {
            return fail( "hpdBarXRange: expected [84,116], got " + java.util.Arrays.toString( r ) );
        }
        // it must STRADDLE the node's x (100), and doing so must NOT depend on any tree height: an asymmetric node
        // (x=50, age 5 in [2,20]) still brackets x=50: left=50-(20-5)*2=20, right=50+(5-2)*2=56
        final float[] a = TreePanelUtil.hpdBarXRange( 50f, 5.0, 2.0, 20.0, 2.0 );
        if ( !( ( r[ 0 ] < 100f ) && ( 100f < r[ 1 ] ) && ( a[ 0 ] < 50f ) && ( 50f < a[ 1 ] ) ) ) {
            return fail( "the bar must straddle the node's x: " + java.util.Arrays.toString( r ) + " / "
                    + java.util.Arrays.toString( a ) );
        }
        return true;
    }

    /** The node-age SPINDLE profile: 0 at both interval ends, peaks at the point estimate, asymmetric when the estimate
     *  is off-centre, and degenerate-safe. */
    /**
     * The protein-domain box-height clamp: track the tip-row spacing (getYdistance), but clamp into [min,max] so the
     * boxes stay readable when zoomed out (floor) and stay bars-not-blocks when zoomed in (ceiling). The user tuned
     * this across rounds 3-5 to MIN=6/MAX=16; shared by both domain-height sites in TreePanel.
     */
    private static boolean testDomainBoxHeight() {
        if ( TreePanelUtil.domainBoxHeight( 1.79f, 6, 16 ) != 6 ) { // tightly packed (yd below the floor) -> floor
            return fail( "a tightly-packed row (spacing below min) must floor at min" );
        }
        if ( TreePanelUtil.domainBoxHeight( 40f, 6, 16 ) != 16 ) { // very expanded (yd above the ceiling) -> ceiling
            return fail( "a very expanded row (spacing above max) must cap at max" );
        }
        if ( TreePanelUtil.domainBoxHeight( 9f, 6, 16 ) != 9 ) { // in range -> the (rounded) spacing
            return fail( "an in-range row must use the rounded spacing" );
        }
        if ( TreePanelUtil.domainBoxHeight( 9.6f, 6, 16 ) != 10 ) { // ROUNDED, not truncated
            return fail( "the spacing must be rounded, not truncated" );
        }
        // RESPONDS to vertical zoom once above the floor: a wider row spacing gives a taller box
        if ( TreePanelUtil.domainBoxHeight( 12f, 6, 16 ) <= TreePanelUtil.domainBoxHeight( 8f, 6, 16 ) ) {
            return fail( "domain box height must grow with row spacing above the floor" );
        }
        return true;
    }

    private static boolean testSpindleHalfHeight() {
        // interval [0,10], estimate (peak) at 3, max half-height 5
        if ( TreePanelUtil.spindleHalfHeightAt( 0, 0, 10, 3, 5 ) != 0 ) {
            return fail( "spindle must be 0 at the low end" );
        }
        if ( TreePanelUtil.spindleHalfHeightAt( 10, 0, 10, 3, 5 ) != 0 ) {
            return fail( "spindle must be 0 at the high end" );
        }
        if ( Math.abs( TreePanelUtil.spindleHalfHeightAt( 3, 0, 10, 3, 5 ) - 5 ) > 1e-9 ) {
            return fail( "spindle must reach h_max at the peak" );
        }
        // ASYMMETRIC: at the SAME absolute distance (1.5) from the off-centre peak (3), the short lobe [0,3] is lower
        // than the long lobe [3,10] -- because the short lobe rises/falls more steeply
        final double left = TreePanelUtil.spindleHalfHeightAt( 1.5, 0, 10, 3, 5 );  // 1.5 below the peak (short lobe)
        final double right = TreePanelUtil.spindleHalfHeightAt( 4.5, 0, 10, 3, 5 ); // 1.5 above the peak (long lobe)
        if ( !( ( left > 0 ) && ( left < 5 ) && ( right > 0 ) && ( right < 5 ) ) ) {
            return fail( "spindle lobe heights must be between 0 and h_max" );
        }
        if ( left >= right ) {
            return fail( "the spindle must be asymmetric: the short lobe is lower at the same distance from the peak "
                    + "(left=" + left + " right=" + right + ")" );
        }
        // degenerate: zero-width interval or non-positive h_max -> 0
        if ( ( TreePanelUtil.spindleHalfHeightAt( 5, 5, 5, 5, 5 ) != 0 )
                || ( TreePanelUtil.spindleHalfHeightAt( 3, 0, 10, 3, 0 ) != 0 ) ) {
            return fail( "spindle must be 0 for a degenerate interval or non-positive h_max" );
        }
        return true;
    }

    /**
     * Labeled-axis tick VALUES: 0 (the root) then stepping by spacing up to and including max_distance; degenerate
     * inputs yield none. (Distinct from the grid lines, which skip 0.)
     */
    private static boolean testCalendarTickYears() {
        // a SARS-CoV-2-scale span (2019.9 .. 2022.5) -> whole-year ticks 2020, 2021, 2022
        final double[] a = TreePanelUtil.calendarTickYears( 2019.9, 2022.5 );
        if ( ( a.length != 3 ) || ( a[ 0 ] != 2020.0 ) || ( a[ 1 ] != 2021.0 ) || ( a[ 2 ] != 2022.0 ) ) {
            return fail( "expected calendar ticks 2020,2021,2022; got " + java.util.Arrays.toString( a ) );
        }
        // every tick is a whole year (no decimal-year labels)
        for ( final double y : a ) {
            if ( y != Math.rint( y ) ) {
                return fail( "calendar ticks must be whole years, got " + y );
            }
        }
        // a ~15-year span must use a WHOLE-year step (regression: a 2.5-year step gave half-year ticks like 2007.5)
        final double[] mid = TreePanelUtil.calendarTickYears( 2005.0, 2020.0 );
        for ( final double y : mid ) {
            if ( y != Math.rint( y ) ) {
                return fail( "a 15-year span must use whole-year ticks (no 2.5-year step), got " + y );
            }
        }
        // a wider span uses a coarser (still whole-year) step and stays within [from, to]
        final double[] b = TreePanelUtil.calendarTickYears( 1850.0, 2020.0 ); // span 170 -> ~25-year step
        if ( b.length < 3 ) {
            return fail( "a 170-year span should yield several ticks; got " + java.util.Arrays.toString( b ) );
        }
        for ( final double y : b ) {
            if ( ( y < 1850.0 ) || ( y > 2020.0 ) || ( y != Math.rint( y ) ) ) {
                return fail( "wide-span ticks must be whole years within [1850,2020], got " + y );
            }
        }
        // niceYearStep is at least 1 (whole years) and grows with the target
        if ( ( TreePanelUtil.niceYearStep( 0.3 ) != 1.0 ) || ( TreePanelUtil.niceYearStep( 3.0 ) < 1.0 )
                || ( TreePanelUtil.niceYearStep( 40.0 ) < 10.0 ) ) {
            return fail( "niceYearStep must be >= 1 and scale up with the target span" );
        }
        // degenerate: non-positive span -> no ticks
        if ( ( TreePanelUtil.calendarTickYears( 2021.0, 2021.0 ).length != 0 )
                || ( TreePanelUtil.calendarTickYears( 2022.0, 2020.0 ).length != 0 ) ) {
            return fail( "a non-positive calendar span must yield no ticks" );
        }
        return true;
    }

    /** The numeric "Ma before present" geologic-axis ticks: niceAxisStep allows SUB-1 steps (unlike whole-year
     *  niceYearStep), so a shallow tree gets fine ticks and a deep one gets 50/500-Ma ticks. */
    private static boolean testMaAxisTicks() {
        if ( Math.abs( TreePanelUtil.niceAxisStep( 0.3 ) - 0.5 ) > 1e-12 ) {
            return fail( "niceAxisStep(0.3) should be 0.5 (sub-1 allowed), got " + TreePanelUtil.niceAxisStep( 0.3 ) );
        }
        if ( ( Math.abs( TreePanelUtil.niceAxisStep( 35.0 ) - 50.0 ) > 1e-9 )
                || ( Math.abs( TreePanelUtil.niceAxisStep( 543.0 ) - 1000.0 ) > 1e-9 ) ) {
            return fail( "niceAxisStep should round UP to a nice 1/2/5 x 10^k step" );
        }
        if ( ( TreePanelUtil.niceAxisStep( 0.0 ) != 0.0 ) || ( TreePanelUtil.niceAxisStep( -5.0 ) != 0.0 ) ) {
            return fail( "niceAxisStep of a non-positive raw must be 0" );
        }
        final double[] dino = TreePanelUtil.maAxisTickValues( 250.0 ); // 0..250 by 50
        if ( ( dino.length != 6 ) || ( dino[ 0 ] != 0.0 ) || ( Math.abs( dino[ 1 ] - 50.0 ) > 1e-9 )
                || ( Math.abs( dino[ 5 ] - 250.0 ) > 1e-9 ) ) {
            return fail( "maAxisTickValues(250) expected 0..250 by 50, got " + java.util.Arrays.toString( dino ) );
        }
        final double[] shallow = TreePanelUtil.maAxisTickValues( 5.0 ); // fine whole-Ma ticks 0..5
        if ( ( shallow.length != 6 ) || ( Math.abs( shallow[ 1 ] - 1.0 ) > 1e-9 ) ) {
            return fail( "maAxisTickValues(5) expected 1-Ma ticks 0..5, got " + java.util.Arrays.toString( shallow ) );
        }
        if ( TreePanelUtil.maAxisTickValues( 0.0 ).length != 0 ) {
            return fail( "maAxisTickValues of a non-positive root age must be empty" );
        }
        return true;
    }

    private static boolean testScaleAxisTickValues() {
        final double[] a = TreePanelUtil.scaleAxisTickValues( 0.7, 0.1 );
        if ( ( a.length != 8 ) || ( a[ 0 ] != 0.0 ) || ( Math.abs( a[ 1 ] - 0.1 ) > 1e-9 )
                || ( Math.abs( a[ 7 ] - 0.7 ) > 1e-9 ) ) {
            return fail( "expected ticks 0..0.7 (8 values incl. 0); got " + java.util.Arrays.toString( a ) );
        }
        // integer spacing, max exactly on a boundary
        final double[] b = TreePanelUtil.scaleAxisTickValues( 30.0, 10.0 );
        if ( ( b.length != 4 ) || ( b[ 0 ] != 0.0 ) || ( b[ 3 ] != 30.0 ) ) {
            return fail( "expected ticks 0,10,20,30; got " + java.util.Arrays.toString( b ) );
        }
        // float tolerance: 3*0.1 = 0.30000000000000004, so max=0.3 must still include the 0.3 tick (not drop it)
        if ( TreePanelUtil.scaleAxisTickValues( 0.3, 0.1 ).length != 4 ) {
            return fail( "a tick whose value rounds to max_distance must be included despite float error" );
        }
        // degenerate: no depth / non-positive spacing -> no ticks
        if ( ( TreePanelUtil.scaleAxisTickValues( 0.0, 0.1 ).length != 0 )
                || ( TreePanelUtil.scaleAxisTickValues( 0.5, 0.0 ).length != 0 )
                || ( TreePanelUtil.scaleAxisTickValues( 0.5, -1.0 ).length != 0 ) ) {
            return fail( "zero-depth or non-positive spacing must yield no axis ticks" );
        }
        // pathological depth/spacing ratio (corrupt/raw-count branch lengths) -> NO axis, not a giant array / int
        // overflow / hang. 1e12 / 100 = 1e10 ticks; well past the ceiling.
        if ( TreePanelUtil.scaleAxisTickValues( 1.0e12, 100.0 ).length != 0 ) {
            return fail( "an absurd tick count must yield no axis (guard against OOM / int overflow)" );
        }
        if ( TreePanelUtil.scaleAxisTickValues( 1.0e15, 1.0e-3 ).length != 0 ) { // would overflow int(n)
            return fail( "a count that overflows int must yield no axis, not a NegativeArraySizeException" );
        }
        return true;
    }

    /** The shared compact figure-label formatter: whole numbers as integers, small magnitudes keep ~3 sig digits. */
    private static boolean testFormatCompactNumber() {
        if ( !"120".equals( TreePanelUtil.formatCompactNumber( 120.0 ) )
                || !"0".equals( TreePanelUtil.formatCompactNumber( 0.0 ) )
                || !"2016.5".equals( TreePanelUtil.formatCompactNumber( 2016.5 ) )
                || !"0.3".equals( TreePanelUtil.formatCompactNumber( 0.30000000000000004 ) ) // float noise -> clean
                || !"0.005".equals( TreePanelUtil.formatCompactNumber( 0.005 ) ) ) { // small magnitude kept legible
            return fail( "formatCompactNumber: got " + TreePanelUtil.formatCompactNumber( 0.30000000000000004 ) + "/"
                    + TreePanelUtil.formatCompactNumber( 0.005 ) );
        }
        return true;
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

        // --- phase 1: cache empty, so only the Rodentia tips are placeable -- via the ancestor FALLBACK, since
        //     Mus/Rattus are genus-annotated (not order) and unknown to the DB ---
        Map<PhylogenyNode, String> assignment = TreePanelUtil.assignTipsToRankTaxon( tree, "order", svc );
        if ( assignment.size() != 2 ) {
            return fail( "with an empty cache only the 2 Rodentia tips should be placeable, got " + assignment.size() );
        }
        final SortedSet<String> unresolved = TreePanelUtil.unresolvedTipTaxa( tree, "order", svc );
        // Felis, Canis, and the unknown tip's query name must all be flagged for fetching
        if ( !unresolved.contains( "Felis" ) || !unresolved.contains( "Canis" ) || !unresolved.contains( "Nonexistus" ) ) {
            return fail( "unresolved tip taxa must include Felis, Canis and Nonexistus; got " + unresolved );
        }
        // "tip identity wins": Mus/Rattus no longer self-resolve (they're genus-annotated, under the Rodentia
        // ANCESTOR), so they too are now flagged for fetch -- an ancestor annotation no longer suppresses it
        if ( ( unresolved.size() != 5 ) || !unresolved.contains( "Mus" ) || !unresolved.contains( "Rattus" ) ) {
            return fail( "a tip under an annotated ancestor must be flagged for fetch (tip identity wins); got "
                    + unresolved );
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
        final Map<PhylogenyNode, TreePanelUtil.RankTaxon> roots = TreePanelUtil
                .maximalMonochromaticRoots( tree, TreePanelUtil.assignNodesToRankTaxon( tree, "order", svc ) );
        if ( roots.size() != 3 ) {
            return fail( "expected 3 maximal monochromatic roots (Rodentia clade, Felis, Canis), got " + roots.size() );
        }
        int carnivora_roots = 0;
        for( final TreePanelUtil.RankTaxon t : roots.values() ) {
            if ( "Carnivora".equals( t.getName() ) ) {
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

    /**
     * "Tip identity wins, internal annotations only fill gaps" (the fix for imperfect manual internal taxonomy):
     * a tip that can resolve its OWN rank taxon overrides a wrong ancestor annotation; an unresolvable tip falls
     * back to the ancestor annotation; and a tip under an annotated ancestor is still flagged for fetching so its
     * own identity CAN be resolved (and thereby override the ancestor).
     */
    private static boolean testTipIdentityWins() {
        // fixture: an internal node WRONGLY annotated order=Rodentia over a single cat tip (Felis)
        final FakeLineageService svc = new FakeLineageService();
        svc.know( "Felis", lineage( "order", "Carnivora", "genus", "Felis" ) );
        Phylogeny tree = wrongAncestorTree();
        PhylogenyNode felis = findLeaf( tree, "Felis" );
        // (c) even under the annotated ancestor, the tip must be flagged for fetch (empty cache)
        if ( !TreePanelUtil.unresolvedTipTaxa( tree, "order", svc ).contains( "Felis" ) ) {
            return fail( "a tip under an annotated ancestor must still be flagged for fetch" );
        }
        try {
            svc.fetch( "Felis" ); // warm the cache with the tip's own resolution
        }
        catch ( final Exception e ) {
            return fail( "fake fetch must not throw: " + e );
        }
        // (a) tip identity wins: Felis resolves its OWN order (Carnivora) and is NOT overridden by the wrong ancestor
        if ( !"Carnivora".equals( TreePanelUtil.assignTipsToRankTaxon( tree, "order", svc ).get( felis ) ) ) {
            return fail( "a tip that resolves its own order must win over a wrong ancestor annotation" );
        }
        // (b) fallback: a tip that cannot resolve itself (unknown name, empty cache) uses the ancestor annotation
        final FakeLineageService empty = new FakeLineageService();
        tree = wrongAncestorTree();
        felis = findLeaf( tree, "Felis" );
        if ( !"Rodentia".equals( TreePanelUtil.assignTipsToRankTaxon( tree, "order", empty ).get( felis ) ) ) {
            return fail( "an unresolvable tip must fall back to the nearest ancestor annotation" );
        }
        return true;
    }

    /** A cat tip (genus Felis) whose immediate ancestor is deliberately mis-annotated as order Rodentia. */
    private static Phylogeny wrongAncestorTree() {
        final PhylogenyNode wrong = internalOrder( "Rodentia" );
        wrong.addAsChild( genusLeaf( "Felis" ) );
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( wrong );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /**
     * Spine B tax-id keying: HOMONYMS (two distinct taxa sharing a scientific name at a rank) are NOT merged into
     * one clade, while a SYNONYM (two names sharing a tax-id) groups as one -- the grouping keys on the NCBI tax-id
     * ({@link TaxonLineage#taxIdAt}), read here straight off each tip's own {@code <taxonomy>} identifier.
     */
    private static boolean testRankColorizationTaxIdKeying() {
        // homonym: two "Drosophila" tips (fly genus 7215 vs a homonymous genus 999999), one in each child clade.
        // Name-keying would make the parent uniform ("Drosophila") = ONE merged clade; tax-id keying keeps them apart.
        final PhylogenyNode da = genusLeafId( "da", "Drosophila", "7215" );
        final PhylogenyNode db = genusLeafId( "db", "Drosophila", "999999" );
        final PhylogenyNode root_h = node( da, db );
        final Phylogeny homonyms = phylogeny( root_h );
        final Map<PhylogenyNode, TreePanelUtil.RankTaxon> h_roots = TreePanelUtil
                .maximalMonochromaticRoots( homonyms, TreePanelUtil.assignNodesToRankTaxon( homonyms, "genus", null ) );
        if ( ( h_roots.size() != 2 ) || h_roots.containsKey( root_h ) ) {
            return fail( "homonyms (same name, different tax-id) must NOT merge: expected 2 roots and no merged parent, "
                    + "got " + h_roots.size() + " roots (parent merged: " + h_roots.containsKey( root_h ) + ")" );
        }

        // synonym: two names ("Drosophila" / "Sophophora") sharing tax-id 7215 group as ONE clade at their parent
        final PhylogenyNode root_s = node( genusLeafId( "sa", "Drosophila", "7215" ),
                                           genusLeafId( "sb", "Sophophora", "7215" ) );
        final Phylogeny synonyms = phylogeny( root_s );
        final Map<PhylogenyNode, TreePanelUtil.RankTaxon> s_roots = TreePanelUtil
                .maximalMonochromaticRoots( synonyms, TreePanelUtil.assignNodesToRankTaxon( synonyms, "genus", null ) );
        if ( ( s_roots.size() != 1 ) || !s_roots.containsKey( root_s ) ) {
            return fail( "a synonym (different names, same tax-id) must group as ONE clade at the parent; got "
                    + s_roots.size() + " roots" );
        }
        return true;
    }

    /**
     * Spine B internal-node assignment + gap-fill (what the tip-only, ancestor-fallback path could not do): an
     * INTERNAL node resolved to the coloring rank -- here a GENUS-annotated node whose ORDER comes from the service,
     * so the tips' ancestor-fallback (which keys on the exact coloring rank) misses it -- (a) colors a clade of
     * otherwise-unresolvable tips, and (c) sweeps an unplaced sibling into its taxon; (b) without such an authorizing
     * ancestor an unplaced tip still breaks the clade ("unplaced tips aren't silently swept").
     */
    private static boolean testInternalTaxonGapFill() {
        final FakeLineageService svc = new FakeLineageService();
        svc.know( "Felis", lineage( "order", "Carnivora", "genus", "Felis" ) );
        try {
            svc.fetch( "Felis" ); // warm the cache-only lineageOf
        }
        catch ( final Exception e ) {
            return fail( "fake fetch must not throw: " + e );
        }

        // (a) a GENUS-annotated internal node (order comes from the service) gap-fills a clade of unresolvable tips
        final PhylogenyNode felis_clade = genusInternal( "Felis" );
        felis_clade.addAsChild( bareLeaf( "unknown_a" ) );
        felis_clade.addAsChild( bareLeaf( "unknown_b" ) );
        final Phylogeny gap = phylogeny( node( felis_clade, orderLeaf( "rat", "Rodentia", null ) ) );
        final Map<PhylogenyNode, TreePanelUtil.RankTaxon> g_roots = TreePanelUtil
                .maximalMonochromaticRoots( gap, TreePanelUtil.assignNodesToRankTaxon( gap, "order", svc ) );
        if ( ( g_roots.get( felis_clade ) == null ) || !"Carnivora".equals( g_roots.get( felis_clade ).getName() ) ) {
            return fail( "an internal node resolved to the rank must gap-fill a clade of unresolvable tips" );
        }

        // (b) an unplaced tip with NO authorizing ancestor still breaks the clade -- only the placed tip is a root
        final PhylogenyNode felis1 = orderLeaf( "felis1", "Carnivora", null );
        final PhylogenyNode plain = node( felis1, bareLeaf( "u1" ) ); // a plain, UN-annotated internal node
        final Phylogeny broken = phylogeny( node( plain, orderLeaf( "rat1", "Rodentia", null ) ) );
        final Map<PhylogenyNode, TreePanelUtil.RankTaxon> b_roots = TreePanelUtil
                .maximalMonochromaticRoots( broken, TreePanelUtil.assignNodesToRankTaxon( broken, "order", svc ) );
        if ( b_roots.containsKey( plain ) || ( b_roots.get( felis1 ) == null ) ) {
            return fail( "an unplaced tip with no authorizing ancestor must break the clade (only the placed tip is a root)" );
        }

        // (c) the same shape, but the clade IS resolved to Carnivora (genus Felis -> order via service) -> the
        //     unplaced sibling is swept in: one Carnivora root at the clade, the placed tip is not its own root
        final PhylogenyNode felis2 = orderLeaf( "felis2", "Carnivora", null );
        final PhylogenyNode auth = genusInternal( "Felis" );
        auth.addAsChild( felis2 );
        auth.addAsChild( bareLeaf( "swept" ) );
        final Phylogeny fill = phylogeny( node( auth, orderLeaf( "rat2", "Rodentia", null ) ) );
        final Map<PhylogenyNode, TreePanelUtil.RankTaxon> c_roots = TreePanelUtil
                .maximalMonochromaticRoots( fill, TreePanelUtil.assignNodesToRankTaxon( fill, "order", svc ) );
        if ( ( c_roots.get( auth ) == null ) || !"Carnivora".equals( c_roots.get( auth ).getName() )
                || c_roots.containsKey( felis2 ) ) {
            return fail( "an authorizing ancestor (resolved to the rank) must sweep an unplaced tip into its taxon" );
        }
        return true;
    }

    /** A leaf with own {@code <taxonomy>} = {sci, rank=genus, ncbi:tax_id}. */
    private static PhylogenyNode genusLeafId( final String node_name, final String sci, final String tax_id ) {
        return rankLeaf( node_name, sci, "genus", tax_id );
    }

    /** A leaf with own {@code <taxonomy>} = {sci, rank=order, optional ncbi:tax_id}. */
    private static PhylogenyNode orderLeaf( final String node_name, final String sci, final String tax_id ) {
        return rankLeaf( node_name, sci, "order", tax_id );
    }

    private static PhylogenyNode rankLeaf( final String node_name, final String sci, final String rank,
                                           final String tax_id ) {
        final PhylogenyNode n = namedLeaf( node_name, sci );
        try {
            n.getNodeData().getTaxonomy().setRank( rank );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        if ( tax_id != null ) {
            n.getNodeData().getTaxonomy().setIdentifier( new Identifier( tax_id, "ncbi" ) );
        }
        return n;
    }

    /** An INTERNAL node with own {@code <taxonomy>} = {sci, rank=genus} (its order comes from the service). */
    private static PhylogenyNode genusInternal( final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        try {
            t.setRank( "genus" );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode node( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode c : children ) {
            n.addAsChild( c );
        }
        return n;
    }

    private static Phylogeny phylogeny( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /**
     * Regression: a tip annotated below the coloring rank WITH an NCBI id must still colorize when the lineage cache
     * was primed by NAME (as the colorize fetch flow does) -- the id lookup misses the name-keyed cache, so resolution
     * must fall back to the name (id-then-name), else id-bearing tips silently go uncolored.
     */
    private static boolean testRankResolutionIdThenNameCache() {
        final FakeLineageService svc = new FakeLineageService();
        svc.know( "Felis catus", lineage( "order", "Carnivora", "species", "Felis catus" ) );
        try {
            svc.fetch( "Felis catus" ); // primes the cache under the NAME, like unresolvedTipTaxa/fetch
        }
        catch ( final Exception e ) {
            return fail( "fake fetch must not throw: " + e );
        }
        final PhylogenyNode cat = rankLeaf( "cat", "Felis catus", "species", "9685" ); // species + NCBI id
        final Phylogeny tree = phylogeny( node( cat, orderLeaf( "rat", "Rodentia", null ) ) );
        final TreePanelUtil.RankTaxon rt = TreePanelUtil.assignNodesToRankTaxon( tree, "order", svc ).get( cat );
        if ( ( rt == null ) || !"Carnivora".equals( rt.getName() ) ) {
            return fail( "an id-annotated tip must resolve via the name-keyed cache (id lookup falls back to name)" );
        }
        return true;
    }

    /**
     * Regression: a tip carrying an effectively NAMELESS {@code <taxonomy>} (rank-only, no sci/code/common, no NCBI
     * id) but whose NODE NAME is the resolvable organism must still be placed by its node name -- the has-taxonomy
     * path finds no own name, so the node-name fallback (tipQueryName) must not be gated away.
     */
    private static boolean testNamelessTaxonomyTipResolvesByNodeName() {
        final FakeLineageService svc = new FakeLineageService();
        svc.know( "Mus", lineage( "order", "Rodentia", "genus", "Mus" ) );
        try {
            svc.fetch( "Mus" );
        }
        catch ( final Exception e ) {
            return fail( "fake fetch must not throw: " + e );
        }
        final PhylogenyNode tip = namelessTaxonomyLeaf( "Mus" ); // node name "Mus"; taxonomy present but nameless
        final Phylogeny tree = phylogeny( node( tip, orderLeaf( "felis", "Carnivora", null ) ) );
        final TreePanelUtil.RankTaxon rt = TreePanelUtil.assignNodesToRankTaxon( tree, "order", svc ).get( tip );
        if ( ( rt == null ) || !"Rodentia".equals( rt.getName() ) ) {
            return fail( "a tip with a nameless taxonomy must still resolve by its node name" );
        }
        return true;
    }

    /**
     * The full A-&gt;B-&gt;C loop: ancestral-taxonomy INFERENCE (Spine C) writes a rank-annotated internal
     * {@code <taxonomy>}, and the node-&gt;taxon ASSIGNMENT (Spine B) reads it -- so an inferred internal taxon FEEDS
     * the rank colorize / clade bands. This is the payoff of the internal-node-taxonomy overhaul, guarded end-to-end.
     */
    private static boolean testInferenceFeedsRankAssignment() {
        final FakeLineageService svc = new FakeLineageService();
        svc.know( "Homo", lineage( "class", "Mammalia", "order", "Primates", "genus", "Homo" ) );
        svc.know( "Pan", lineage( "class", "Mammalia", "order", "Primates", "genus", "Pan" ) );
        try {
            svc.fetch( "Homo" );
            svc.fetch( "Pan" );
        }
        catch ( final Exception e ) {
            return fail( "fake fetch must not throw: " + e );
        }
        // two bare-named primate tips under one internal node; NO internal taxonomy yet
        final PhylogenyNode mrca = node( bareLeaf( "Homo" ), bareLeaf( "Pan" ) );
        final Phylogeny tree = phylogeny( mrca );
        // (C) infer: the tips resolve to the same order, so the MRCA is assigned order=Primates
        final Map<PhylogenyNode, TaxonLineage> tl = TreePanelUtil.tipLineages( tree, svc );
        final org.forester.analysis.AncestralTaxonomyInference.InferenceResult r = org.forester.analysis.AncestralTaxonomyInference
                .inferInternalTaxonomies( tree, tl, false );
        if ( r.getAssigned() < 1 ) {
            return fail( "inference must assign the shared order to the internal node; assigned=" + r.getAssigned() );
        }
        if ( !mrca.getNodeData().isHasTaxonomy() || !"order".equalsIgnoreCase( mrca.getNodeData().getTaxonomy().getRank() )
                || !"Primates".equals( mrca.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "the inferred internal taxonomy must carry rank=order, name=Primates" );
        }
        // (B) the node->taxon assignment READS that inferred internal taxon -- inference now feeds colorize
        final TreePanelUtil.RankTaxon rt = TreePanelUtil.assignNodesToRankTaxon( tree, "order", svc ).get( mrca );
        if ( ( rt == null ) || !"Primates".equals( rt.getName() ) ) {
            return fail( "the node->taxon assignment must read the inferred internal taxon at its rank" );
        }
        return true;
    }

    /**
     * "Annotate Clades by Rank -&gt; write into the tree": each maximal-monochromatic clade's taxon at the rank is
     * WRITTEN onto that clade's INTERNAL root {@code <taxonomy>} (name + rank + NCBI id), a tip keeps its own (more
     * specific) taxonomy, and an already-annotated root is preserved unless overwrite.
     */
    private static boolean testWriteCladeTaxonomies() {
        // a plain internal clade of two order-Carnivora tips (tax-id 33554); a SPECIES-level outgroup tip whose order
        // (Rodentia) comes from the service -- so if a tip were wrongly written, its species would be DOWNGRADED
        final FakeLineageService svc = new FakeLineageService();
        svc.know( "Rattus", lineage( "order", "Rodentia", "genus", "Rattus" ) );
        try {
            svc.fetch( "Rattus" );
        }
        catch ( final Exception e ) {
            return fail( "fake fetch must not throw: " + e );
        }
        final PhylogenyNode carn = node( orderLeaf( "felis", "Carnivora", "33554" ),
                                         orderLeaf( "panthera", "Carnivora", "33554" ) );
        final PhylogenyNode rat = rankLeaf( "rat", "Rattus", "species", null ); // species-level; order via the service
        final Phylogeny tree = phylogeny( node( carn, rat ) );
        final int wrote = TreePanelUtil.writeCladeTaxonomies( tree, "order", svc, false );
        if ( wrote != 1 ) {
            return fail( "write should annotate exactly the 1 internal clade root, got " + wrote );
        }
        if ( !carn.getNodeData().isHasTaxonomy() || !"order".equalsIgnoreCase( carn.getNodeData().getTaxonomy().getRank() )
                || !"Carnivora".equals( carn.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "the internal clade root must be written with rank=order, name=Carnivora" );
        }
        final Identifier id = carn.getNodeData().getTaxonomy().getIdentifier();
        if ( ( id == null ) || !"33554".equals( id.getValue() ) || !"ncbi".equalsIgnoreCase( id.getProvider() ) ) {
            return fail( "the written taxonomy must carry the NCBI tax-id 33554, got " + id );
        }
        if ( !"Rattus".equals( rat.getNodeData().getTaxonomy().getScientificName() )
                || !"species".equalsIgnoreCase( rat.getNodeData().getTaxonomy().getRank() ) ) {
            return fail( "a tip's own (species) taxonomy must NOT be downgraded to its order by the clade write" );
        }
        // preserve: an already-annotated clade root is left alone unless overwrite
        if ( TreePanelUtil.writeCladeTaxonomies( tree, "order", svc, false ) != 0 ) {
            return fail( "an already-annotated clade root must be preserved when overwrite is off" );
        }
        if ( TreePanelUtil.writeCladeTaxonomies( tree, "order", svc, true ) != 1 ) {
            return fail( "overwrite must re-write the clade root" );
        }
        // provenance sentence (pure; the caller APPENDS it)
        tree.setName( "mammals" );
        final String prov = TreePanelUtil.cladeTaxonomyProvenance( tree, "order", 1, false );
        if ( !prov.contains( "rank order" ) || !prov.contains( "1 internal clade node" ) || !prov.contains( "mammals" ) ) {
            return fail( "provenance sentence malformed: " + prov );
        }
        return true;
    }

    /**
     * {@code internalTaxaByRank} groups the DISTINCT internal-node taxa by rank (canonical Linnaean order), counts the
     * internal nodes per taxon, orders taxa by count-desc, and EXCLUDES the tips -- the data for the internal-taxonomy
     * backbone key.
     */
    private static boolean testInternalTaxaByRank() {
        final PhylogenyNode carn1 = internalTax( "order", "Carnivora" );
        carn1.addAsChild( bareLeaf( "felis" ) );
        final PhylogenyNode carn2 = internalTax( "order", "Carnivora" ); // Carnivora on a SECOND internal node -> count 2
        carn2.addAsChild( bareLeaf( "canis" ) );
        final PhylogenyNode fel = internalTax( "family", "Felidae" );
        fel.addAsChild( bareLeaf( "panthera" ) );
        final PhylogenyNode rod = internalTax( "order", "Rodentia" );
        rod.addAsChild( orderLeaf( "mus", "Rodentia", null ) ); // a TIP annotated at order -> excluded from the key
        final Phylogeny tree = phylogeny( node( node( carn1, carn2, fel ), rod ) );
        final java.util.LinkedHashMap<String, java.util.LinkedHashMap<String, Integer>> m = TreePanelUtil
                .internalTaxaByRank( tree );
        final java.util.List<String> ranks = new java.util.ArrayList<>( m.keySet() );
        if ( ( ranks.size() != 2 ) || !"order".equals( ranks.get( 0 ) ) || !"family".equals( ranks.get( 1 ) ) ) {
            return fail( "internalTaxaByRank must group by rank in canonical order (order, family); got " + ranks );
        }
        final java.util.LinkedHashMap<String, Integer> orders = m.get( "order" );
        if ( !Integer.valueOf( 2 ).equals( orders.get( "Carnivora" ) )
                || !Integer.valueOf( 1 ).equals( orders.get( "Rodentia" ) ) ) {
            return fail( "Carnivora must count 2 internal nodes, Rodentia 1; got " + orders );
        }
        if ( !"Carnivora".equals( orders.keySet().iterator().next() ) ) {
            return fail( "taxa within a rank must be ordered by count desc (Carnivora first)" );
        }
        if ( ( m.get( "family" ).size() != 1 ) || !m.get( "family" ).containsKey( "Felidae" ) ) {
            return fail( "the family rank must list exactly Felidae; got " + m.get( "family" ) );
        }
        return true;
    }

    /** An INTERNAL node with own {@code <taxonomy>} = {sci, rank}. */
    private static PhylogenyNode internalTax( final String rank, final String sci ) {
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

    /** A leaf whose node name is {@code node_name} and whose {@code <taxonomy>} is present but NAMELESS (rank only). */
    private static PhylogenyNode namelessTaxonomyLeaf( final String node_name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( node_name );
        final Taxonomy t = new Taxonomy();
        try {
            t.setRank( "genus" );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
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

        private final Map<String, TaxonLineage> _db    = new HashMap<String, TaxonLineage>();
        private final Map<String, TaxonLineage> _cache = new HashMap<String, TaxonLineage>();

        void know( final String name, final Map<String, String> rank_to_name ) {
            // model a resolved lineage as rank-carrying ancestors (own taxon left blank); TaxonLineage.at(rank)
            // then returns the name, exactly as the old RankedLineage did
            final ArrayList<TaxonLineage.Ancestor> anc = new ArrayList<TaxonLineage.Ancestor>();
            for( final Map.Entry<String, String> e : rank_to_name.entrySet() ) {
                anc.add( new TaxonLineage.Ancestor( e.getValue(), e.getKey(), null ) );
            }
            _db.put( name.toLowerCase( Locale.ROOT ), new TaxonLineage( null, null, null, null, anc ) );
        }

        @Override
        public TaxonLineage lineageOf( final String taxon ) {
            return ( taxon == null ) ? null : _cache.get( taxon.toLowerCase( Locale.ROOT ) );
        }

        @Override
        public TaxonLineage fetch( final String taxon ) {
            final String k = taxon.toLowerCase( Locale.ROOT );
            final TaxonLineage rl = _db.containsKey( k ) ? _db.get( k ) : TaxonLineage.EMPTY;
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

    /** The counts_out overloads report each taxon's tip count (from the same assignment used to color). */
    private static boolean testRankTaxonCounts() {
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
        // 2 Rodentia + Felis + Canis, with Felis & Canis both order Carnivora -> Carnivora:2, Rodentia:2
        final Map<String, Integer> counts = new LinkedHashMap<String, Integer>();
        TreePanelUtil.colorPhylogenyAccordingToRanks( tree, "order", svc, new LinkedHashMap<String, Color>(), null,
                                                      counts );
        if ( ( counts.size() != 2 ) || ( intval( counts, "Carnivora" ) != 2 ) || ( intval( counts, "Rodentia" ) != 2 ) ) {
            return fail( "rank counts should be Carnivora=2, Rodentia=2; got " + counts );
        }
        // the clade-band overload reports the same per-taxon counts
        final Map<String, Integer> band_counts = new LinkedHashMap<String, Integer>();
        TreePanelUtil.cladeBands( tree, "order", svc, null, band_counts );
        if ( !counts.equals( band_counts ) ) {
            return fail( "cladeBands counts must match the colorize counts; got " + band_counts );
        }
        // a null counts_out is safe (the pre-counts callers pass null)
        TreePanelUtil.colorPhylogenyAccordingToRanks( tree, "order", svc, new LinkedHashMap<String, Color>(), null,
                                                      null );
        return true;
    }

    private static int intval( final Map<String, Integer> m, final String k ) {
        final Integer v = m.get( k );
        return ( v == null ) ? -1 : v.intValue();
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
        if ( TreePanelUtil.confidenceScaleMaxFor( 100.0 ) != 100.0 ) {
            return fail( "max exactly 100 stays on the 0..100 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 100.5 ) != 1000.0 ) {
            return fail( "any value above 100 implies the 0..1000 scale" );
        }
        if ( TreePanelUtil.confidenceScaleMaxFor( 800.0 ) != 1000.0 ) {
            return fail( "max 800 implies the 0..1000 scale" );
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

    /** medianPositiveBranchLength (positive-only, so zeros don't poison it), longBranchBreakCap, and cappedTreeHeight
     *  (a long branch contributes only the cap; a cap larger than every branch == the ordinary tree height). */
    private static boolean testBreakLongBranches() {
        // ((A:1,B:1):1, OUT:40) -- one outlier branch, the rest ~1
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode clade = blChild( root, null, 1 );
        blChild( clade, "A", 1 );
        blChild( clade, "B", 1 );
        blChild( root, "OUT", 40 );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        // positive dtps = {1,1,1,40} -> median 1 (the root's own default negative distance is excluded)
        if ( Math.abs( TreePanelUtil.medianPositiveBranchLength( phy ) - 1.0 ) > 1e-9 ) {
            return fail( "median positive branch length must be 1, got " + TreePanelUtil.medianPositiveBranchLength( phy ) );
        }
        final double cap = TreePanelUtil.longBranchBreakCap( phy, 8.0 );
        if ( Math.abs( cap - 8.0 ) > 1e-9 ) {
            return fail( "cap must be 8 * median = 8, got " + cap );
        }
        // cappedTreeHeight: root->clade->A = 2, root->OUT = min(40,8) = 8 -> deepest capped path 8 (< true height 40)
        if ( Math.abs( TreePanelUtil.cappedTreeHeight( phy, cap ) - 8.0 ) > 1e-9 ) {
            return fail( "capped height must be 8, got " + TreePanelUtil.cappedTreeHeight( phy, cap ) );
        }
        // a cap larger than every branch reproduces the ordinary tree height (no rescale for a well-behaved tree)
        final double huge = TreePanelUtil.cappedTreeHeight( phy, 1.0e9 );
        if ( ( Math.abs( huge - 40.0 ) > 1e-9 ) || ( Math.abs( huge - phy.calculateHeight( true ) ) > 1e-9 ) ) {
            return fail( "capped height with a huge cap must equal the ordinary height 40, got " + huge );
        }
        // a long INTERNAL branch (a whole clade behind it) must be capped too, pulling its subtree in: ((A:1,B:1):40, C:1)
        // -> the internal branch caps to 8, so the deepest path root->clade(8)->A(1) = 9 (not 41), NOT off at the true 41
        final PhylogenyNode ir = new PhylogenyNode();
        final PhylogenyNode iclade = blChild( ir, null, 40 ); // the long internal branch
        final PhylogenyNode ia = blChild( iclade, "A", 1 );
        blChild( iclade, "B", 1 );
        final PhylogenyNode ic = blChild( ir, "C", 1 );
        final Phylogeny iphy = new Phylogeny();
        iphy.setRoot( ir );
        iphy.externalNodesHaveChanged();
        final double icap = TreePanelUtil.longBranchBreakCap( iphy, 8.0 ); // median {1,1,1,40} = 1 -> cap 8
        if ( ( Math.abs( icap - 8.0 ) > 1e-9 )
                || ( Math.abs( TreePanelUtil.cappedTreeHeight( iphy, icap ) - 9.0 ) > 1e-9 ) ) {
            return fail( "a long internal branch must cap (cap 8, capped height 9), got cap=" + icap + " height="
                    + TreePanelUtil.cappedTreeHeight( iphy, icap ) );
        }
        // cappedDistanceToRoot places a node's radius/spoke radially: it sums CAPPED branch lengths to the root, so a tip
        // behind the long internal branch pulls in (A: min(1,8)+min(40,8) = 9, NOT the true 41) while a normal tip is exact
        // (C: 1); the root is 0. This is the radial (circular/unrooted) analogue of the rectangular capped depth.
        if ( ( Math.abs( TreePanelUtil.cappedDistanceToRoot( ia, icap ) - 9.0 ) > 1e-9 )
                || ( Math.abs( TreePanelUtil.cappedDistanceToRoot( ic, icap ) - 1.0 ) > 1e-9 )
                || ( TreePanelUtil.cappedDistanceToRoot( ir, icap ) != 0.0 ) ) {
            return fail( "cappedDistanceToRoot wrong (A must be 9, C must be 1, root 0), got A="
                    + TreePanelUtil.cappedDistanceToRoot( ia, icap ) + " C=" + TreePanelUtil.cappedDistanceToRoot( ic, icap )
                    + " root=" + TreePanelUtil.cappedDistanceToRoot( ir, icap ) );
        }
        // zeros (polytomy branches) are excluded from the median: {0,0,2,4,6} -> median 4, not 2
        final PhylogenyNode r2 = new PhylogenyNode();
        blChild( r2, "a", 0 );
        blChild( r2, "b", 0 );
        blChild( r2, "c", 2 );
        blChild( r2, "d", 4 );
        blChild( r2, "e", 6 );
        final Phylogeny phy2 = new Phylogeny();
        phy2.setRoot( r2 );
        phy2.externalNodesHaveChanged();
        if ( Math.abs( TreePanelUtil.medianPositiveBranchLength( phy2 ) - 4.0 ) > 1e-9 ) {
            return fail( "median must ignore zero-length branches (expect 4), got "
                    + TreePanelUtil.medianPositiveBranchLength( phy2 ) );
        }
        // a POSITIVE root branch is included in the capped height (parity with calculateHeight, which adds it), so a
        // well-behaved tree (nothing over the cap) is drawn unchanged instead of overshooting
        final PhylogenyNode r4 = new PhylogenyNode();
        r4.setDistanceToParent( 3.0 ); // the root's own incoming branch
        final PhylogenyNode c4 = blChild( r4, null, 1 );
        blChild( c4, "a", 1 );
        blChild( c4, "b", 1 );
        blChild( r4, "c", 2 );
        final Phylogeny phy4 = new Phylogeny();
        phy4.setRoot( r4 );
        phy4.externalNodesHaveChanged();
        if ( Math.abs( TreePanelUtil.cappedTreeHeight( phy4, 1.0e9 ) - phy4.calculateHeight( true ) ) > 1e-9 ) {
            return fail( "capped height (huge cap) must equal calculateHeight incl. the root branch, got "
                    + TreePanelUtil.cappedTreeHeight( phy4, 1.0e9 ) + " vs " + phy4.calculateHeight( true ) );
        }
        // the RADIAL normalizer EXCLUDES the root branch (matching cappedDistanceToRoot / getMaxDistanceToRoot), so the
        // deepest capped tip fills the ring exactly even with a positive root branch. phy4's deepest tip is 2 from the
        // root (a: 1+1, or c: 2); cappedMaxDistanceToRoot must be 2, i.e. cappedTreeHeight (5, root-included) minus the
        // 3.0 root branch -- else a circular/unrooted tree with a root branch would under-fill the ring/diameter.
        if ( ( Math.abs( TreePanelUtil.cappedMaxDistanceToRoot( phy4, 1.0e9 ) - 2.0 ) > 1e-9 )
                || ( Math.abs( TreePanelUtil.cappedMaxDistanceToRoot( phy4, 1.0e9 )
                        - ( TreePanelUtil.cappedTreeHeight( phy4, 1.0e9 ) - 3.0 ) ) > 1e-9 ) ) {
            return fail( "cappedMaxDistanceToRoot must EXCLUDE the root branch (expect 2 = cappedTreeHeight 5 - root 3), "
                    + "got " + TreePanelUtil.cappedMaxDistanceToRoot( phy4, 1.0e9 ) );
        }
        // collapse-aware: cappedTreeHeight mirrors calculateHeight under display-collapse (a collapsed clade counts only
        // to its root). Deep clade (tips at depth 6) + a shallow tip (depth 2); collapsing the deep clade drops the
        // deepest DRAWN tip to 2 -- so the depth scale/extent must NOT keep measuring through the collapse.
        final PhylogenyNode cr = new PhylogenyNode();
        final PhylogenyNode deep = blChild( cr, null, 1 );
        blChild( deep, "t1", 5 );
        blChild( deep, "t2", 5 );
        blChild( cr, "shallow", 2 );
        final Phylogeny cphy = new Phylogeny();
        cphy.setRoot( cr );
        cphy.externalNodesHaveChanged();
        deep.setCollapse( true );
        if ( Math.abs( cphy.calculateHeight( true ) - cphy.calculateHeight( false ) ) < 1e-9 ) {
            return fail( "test setup: collapsing the deep clade must change calculateHeight" );
        }
        for( final boolean tc : new boolean[] { true, false } ) {
            // with a huge cap (nothing capped), capped height must EXACTLY equal calculateHeight for BOTH collapse modes
            if ( Math.abs( TreePanelUtil.cappedTreeHeight( cphy, 1.0e9, tc ) - cphy.calculateHeight( tc ) ) > 1e-9 ) {
                return fail( "capped height (huge cap, collapse=" + tc + ") must equal calculateHeight, got "
                        + TreePanelUtil.cappedTreeHeight( cphy, 1.0e9, tc ) + " vs " + cphy.calculateHeight( tc ) );
            }
        }
        // scale-bar interval buckets: sized from the DRAWN extent, so a capped tree's bar reflects the ingroup scale
        // (e.g. a deep-outlier tree with maxDist 40 -> bucket 1, but capped extent 1.6 -> bucket 0.1)
        if ( ( TreePanelUtil.niceScaleBarDistance( 0.4 ) != 0.01 ) || ( TreePanelUtil.niceScaleBarDistance( 1.6 ) != 0.1 )
                || ( TreePanelUtil.niceScaleBarDistance( 4.0 ) != 0.1 ) || ( TreePanelUtil.niceScaleBarDistance( 40 ) != 1 )
                || ( TreePanelUtil.niceScaleBarDistance( 400 ) != 10 )
                || ( TreePanelUtil.niceScaleBarDistance( 0 ) != 0.0 ) ) {
            return fail( "niceScaleBarDistance buckets wrong (0.4->0.01, 1.6/4->0.1, 40->1, 400->10, 0->0)" );
        }
        // a branch-length-less tree: no positive branch -> cap 0 -> capping inactive
        final PhylogenyNode r3 = new PhylogenyNode();
        blChild( r3, "x", 0 );
        blChild( r3, "y", 0 );
        final Phylogeny phy3 = new Phylogeny();
        phy3.setRoot( r3 );
        phy3.externalNodesHaveChanged();
        if ( ( TreePanelUtil.medianPositiveBranchLength( phy3 ) != 0 )
                || ( TreePanelUtil.longBranchBreakCap( phy3, 8.0 ) != 0 )
                || ( TreePanelUtil.cappedTreeHeight( phy3, 0 ) != 0 ) ) {
            return fail( "a branch-length-less tree must yield median/cap/capped-height 0 (capping inactive)" );
        }
        return true;
    }

    private static PhylogenyNode blChild( final PhylogenyNode parent, final String name, final double dtp ) {
        final PhylogenyNode n = new PhylogenyNode();
        if ( name != null ) {
            n.setName( name );
        }
        n.setDistanceToParent( dtp );
        parent.addAsChild( n );
        return n;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [TreePanelUtilTest] " + message );
        return false;
    }

    private TreePanelUtilTest() {
        // not instantiable
    }
}
