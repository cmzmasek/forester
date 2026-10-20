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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.forester.archaeopteryx.TreePanelUtil.RankTaxon;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The joint clades-by-rank rules, on the SAME fixtures Archaeopteryx.js pins in test/rank_test.js, with its
 * hand-computed numbers (Christian's decisions, 2026-10-19, ported 0.11.186): a tip's taxon at a rank from its own
 * taxonomy, else a property named for the rank in ANY namespace, else the nearest ancestor naming one; a tip the
 * tree places at no taxon does NOT vote (a clade is a maximal subtree whose placed tips all name one taxon, the
 * unplaced swept in); an internal annotation its tips contradict is reported, never corrected; strain, substrain
 * and isolate are no ranks; a clade / lineage scheme of the file's own is. The 12-tip tree, drawn:
 * <pre>
 *   root
 *   +-- A   taxonomy genus Felis
 *   |   +-- a1  taxonomy species "Felis catus"      genus: Felis, by the ANCESTOR A
 *   |   +-- a2  property vipr:Genus = "felis"       genus: Felis, by the PROPERTY (spelled lower)
 *   |   +-- a3  nothing                             genus: Felis, by the ancestor A
 *   +-- B   nothing
 *   |   +-- b1  property repseq:genus = "Panthera"  genus: Panthera
 *   |   +-- b2  nothing, no annotated ancestor      genus: UNPLACED -- does not vote
 *   |   +-- b3  property x:Genus = "  Panthera  "   genus: Panthera (trimmed)
 *   +-- C   nothing
 *   |   +-- c1  property vipr:Genus = "Felis"       genus: Felis -- a second Felis clade
 *   |   +-- c2  property vipr:Genus = "Lynx"
 *   +-- D   taxonomy genus Canis
 *   |   +-- d1  nothing                             genus: Canis, by the ancestor
 *   |   +-- d2  property foo:Host = "dog"           genus: Canis (Host is no rank)
 *   +-- E   taxonomy genus Vulpes
 *       +-- e1  taxonomy genus Nyctereutes AND property vipr:Genus = "Vulpes": Nyctereutes (own taxonomy first)
 *       +-- e2  property vipr:Genus = "Vulpes"; taxonomy rank other "x": Vulpes (property over ancestor)
 * </pre>
 * Clades at genus: A (Felis), B (Panthera; b2 swept in), c1, c2, D (Canis), e1, e2 -- seven; the root, C and E are
 * mixed. Taxa by placed tips: Felis 4 in 2 clades, Canis 2, Panthera 2, Lynx 1, Nyctereutes 1, Vulpes 1; 12 tips,
 * 11 placed. E's own Vulpes is contradicted by e1. At species only a1 is placed, and the whole tree is one clade of
 * it (the voting rule's reach, pinned on purpose). Then the cross-run tree both programs ship,
 * forester/test_data/rank/Caliciviridae_100.xml, at Archaeopteryx.js's numbers.
 */
public final class RankCladesTest {

    public static void main( final String[] args ) {
        System.out.println( test() ? "RankCladesTest: OK." : "RankCladesTest: FAILED." );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        try {
            sources( ok );
            cladesAtGenus( ok );
            conflicts( ok );
            reach( ok );
            offered( ok );
            caliciviridae( ok );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ok[ 0 ] = fail( "unexpected: " + e );
        }
        return ok[ 0 ];
    }

    // ---- the three sources of a tip's taxon, and their precedence ----
    private static void sources( final boolean[] ok ) {
        final Phylogeny t = fixture();
        final Map<PhylogenyNode, TaxonPlacement> p = TreePanelUtil.placeNodes( t, "genus", null, TaxonSource.TREE );
        check( ok, placed( p, t, "e1", "Nyctereutes" ) && p.get( node( t, "e1" ) ).getSource() == TaxonPlacement.Source.OWN_ANNOTATION,
               "a tip's own taxonomy at the rank, before its property" );
        check( ok, placed( p, t, "b3", "Panthera" ) && placed( p, t, "b1", "Panthera" )
                && p.get( node( t, "b3" ) ).getSource() == TaxonPlacement.Source.OWN_METADATA,
               "a property named for the rank, any namespace, trimmed" );
        check( ok, placed( p, t, "a1", "Felis" ) && p.get( node( t, "a1" ) ).getSource() == TaxonPlacement.Source.ANCESTOR_ANNOTATION
                && placed( p, t, "d2", "Canis" ) && placed( p, t, "e2", "Vulpes" )
                && p.get( node( t, "e2" ) ).getSource() == TaxonPlacement.Source.OWN_METADATA,
               "the nearest ancestor's own taxon after the tip's own sources; Host is no rank; a property over an ancestor" );
        check( ok, !p.get( node( t, "b2" ) ).isPlaced(), "a tip with nothing and no annotated ancestor is unplaced" );
        check( ok, p.get( node( t, "A" ) ) != null && "Felis".equals( p.get( node( t, "A" ) ).getTaxon().getName() )
                && ( p.get( node( t, "B" ) ) == null ), "an internal node is read for its own taxon only" );
        final Map<PhylogenyNode, TaxonPlacement> other = TreePanelUtil.placeNodes( t, "other", null, TaxonSource.TREE );
        check( ok, other.isEmpty() || !other.containsKey( node( t, "e2" ) ) || !other.get( node( t, "e2" ) ).isPlaced(),
               "a rank of 'other' is no rank" );
        final Map<PhylogenyNode, TaxonPlacement> host = TreePanelUtil.placeNodes( t, "host", null, TaxonSource.TREE );
        check( ok, !host.get( node( t, "d2" ) ).isPlaced(), "a property that is no rank places nothing" );
        check( ok, placed( TreePanelUtil.placeNodes( t, "Genus", null, TaxonSource.TREE ), t, "a2", "felis" ),
               "the rank is matched by its key (Genus); a placement keeps the tip's own spelling (felis)" );
        check( ok, "Felis".equals( TreePanelUtil.assignNodesToRankTaxon( t, "genus", null, TaxonSource.TREE ).get( node( t, "a2" ) ).getName() ),
               "the assignment shows one spelling per taxon, the first met in the tree's order (Felis)" );
    }

    // ---- the clades at genus ----
    private static void cladesAtGenus( final boolean[] ok ) {
        final Phylogeny t = fixture();
        final Map<PhylogenyNode, RankTaxon> a = TreePanelUtil.assignNodesToRankTaxon( t, "genus", null, TaxonSource.TREE );
        final Map<PhylogenyNode, RankTaxon> roots = TreePanelUtil.maximalMonochromaticRoots( t, a );
        final List<String> names = new ArrayList<String>();
        final Map<String, String> taxon_of = new HashMap<String, String>();
        for( final Entry<PhylogenyNode, RankTaxon> e : roots.entrySet() ) {
            names.add( e.getKey().getName() );
            taxon_of.put( e.getKey().getName(), e.getValue().getName() );
        }
        names.sort( null );
        check( ok, names.equals( Arrays.asList( "A", "B", "D", "c1", "c2", "e1", "e2" ) ), "seven maximal clades: A B D c1 c2 e1 e2, got " + names );
        check( ok, "Felis".equals( taxon_of.get( "A" ) ) && "Panthera".equals( taxon_of.get( "B" ) ) && "Felis".equals( taxon_of.get( "c1" ) )
                && "Lynx".equals( taxon_of.get( "c2" ) ) && "Canis".equals( taxon_of.get( "D" ) ) && "Nyctereutes".equals( taxon_of.get( "e1" ) )
                && "Vulpes".equals( taxon_of.get( "e2" ) ), "each clade's taxon: " + taxon_of );
        check( ok, !roots.containsKey( node( t, "C" ) ) && !roots.containsKey( node( t, "E" ) ) && !roots.containsKey( t.getRoot() ),
               "C, E and the root are mixed" );
        check( ok, !a.containsKey( node( t, "b2" ) ) && node( t, "b2" ).getParent() == node( t, "B" ),
               "b2 is swept into B's clade and names no taxon" );
        // the legend: taxa by placed tips, one spelling, the unplaced count
        final Map<String, Integer> counts = new TreeMap<String, Integer>();
        final Map<String, java.awt.Color> legend = new HashMap<String, java.awt.Color>();
        TreePanelUtil.colorPhylogenyAccordingToRanks( t, "genus", null, legend, null, counts, TaxonSource.TREE );
        check( ok, counts.equals( map( "Felis", 4, "Canis", 2, "Panthera", 2, "Lynx", 1, "Nyctereutes", 1, "Vulpes", 1 ) ),
               "taxa by placed tips, grouped case-insensitively, the first spelling shown: " + counts );
        check( ok, legend.size() == 6 && legend.containsKey( "Felis" ) && !legend.containsKey( "felis" ), "the legend lists Felis once" );
        check( ok, TreePanelUtil.unplacedTipCount( t, a ) == 1, "one tip unplaced (the legend's no-value row)" );
        check( ok, node( t, "b2" ).getBranchData().getBranchColor() == null
                && node( t, "b1" ).getBranchData().getBranchColor() != null,
               "colorized: the swept-in tip's own branch stays uncoloured, its placed sibling's is coloured" );
    }

    // ---- an internal annotation its tips contradict is reported, never corrected ----
    private static void conflicts( final boolean[] ok ) {
        final Phylogeny t = fixture();
        final CladeAnnotationReport r = CladeAnnotationReport.of( t, Arrays.asList( "genus" ), null, TaxonSource.TREE );
        int internal = 0;
        String e_note = null;
        for( final CladeAnnotationReport.Row row : r.getRows() ) {
            if ( row.isInternal() ) {
                ++internal;
                if ( row.getNode().contains( "E" ) ) {
                    e_note = row.getNote();
                }
            }
        }
        check( ok, internal == 1 && ( e_note != null ) && e_note.contains( "Nyctereutes" ), "one contradicted internal node, E (Vulpes) against e1 (Nyctereutes); A's felis is the same taxon; D's tips name nothing: " + internal + " " + e_note );
        check( ok, "Vulpes".equals( node( t, "E" ).getNodeData().getTaxonomy().getScientificName() ), "never corrected" );
        check( ok, CladeAnnotationReport.of( t, Arrays.asList( "species", "family" ), null, TaxonSource.TREE ).getRows().stream()
                .noneMatch( CladeAnnotationReport.Row::isInternal ), "no internal annotation at the rank, no conflict" );
        // a nearer internal annotation contradicts a farther one through the tips it places: F says Felis, its
        // child G says Lynx, the tips name nothing of their own, so they take G's Lynx -- F is contradicted, G not
        final PhylogenyNode g = clade( "G", tax( "genus", "Lynx" ), tip( "g1" ), tip( "g2" ) );
        final PhylogenyNode f = clade( "F", tax( "genus", "Felis" ), g );
        final Phylogeny n = tree( clade( "root", null, f, tip( "h1", prop( "x:genus", "Felis" ) ) ) );
        int f_rows = 0, g_rows = 0;
        for( final CladeAnnotationReport.Row row : CladeAnnotationReport.of( n, Arrays.asList( "genus" ), null, TaxonSource.TREE ).getRows() ) {
            if ( row.isInternal() ) {
                if ( row.getNode().contains( "F" ) && row.getNote().contains( "Lynx" ) ) {
                    ++f_rows;
                }
                if ( row.getNode().contains( "G" ) ) {
                    ++g_rows;
                }
            }
        }
        check( ok, f_rows == 1 && g_rows == 0, "F contradicted by g1 (Lynx, placed by G), G not: " + f_rows + "/" + g_rows );
        // an INTERNAL node's own rank FIELD counts as its annotation does: it places the bare tips below it, and a tip
        // that names another genus contradicts it, reported as a node row
        final PhylogenyNode m = clade( "M", prop( "x:genus", "Felis" ), tip( "m1" ), tip( "m2" ), tip( "m3", prop( "vipr:Genus", "Lynx" ) ) );
        final Phylogeny fields = tree( clade( "root", null, m, tip( "o1", prop( "x:genus", "Canis" ) ) ) );
        final Map<PhylogenyNode, TaxonPlacement> fp = TreePanelUtil.placeNodes( fields, "genus", null, TaxonSource.TREE );
        check( ok, placed( fp, fields, "m1", "Felis" ) && fp.get( node( fields, "m1" ) ).getSource() == TaxonPlacement.Source.ANCESTOR_ANNOTATION
                && placed( fp, fields, "m3", "Lynx" ) && ( fp.get( node( fields, "M" ) ) != null )
                && fp.get( node( fields, "M" ) ).getSource() == TaxonPlacement.Source.OWN_METADATA,
               "an internal node's own field places the bare tips below it, and is the node's own placement" );
        int m_rows = 0;
        for( final CladeAnnotationReport.Row row : CladeAnnotationReport.of( fields, Arrays.asList( "genus" ), null, TaxonSource.TREE ).getRows() ) {
            if ( row.isInternal() && row.getNode().contains( "M" ) && "Felis".equals( row.getTaxon() ) && row.getNote().contains( "Lynx" ) ) {
                ++m_rows;
            }
        }
        check( ok, m_rows == 1, "a node whose own FIELD a tip contradicts is reported like an annotated one: " + m_rows );
    }

    // ---- the voting rule's reach, a pure tree, nameless taxonomies ----
    private static void reach( final boolean[] ok ) {
        final Phylogeny t = fixture();
        final Map<PhylogenyNode, RankTaxon> sp = TreePanelUtil.assignNodesToRankTaxon( t, "species", null, TaxonSource.TREE );
        final Map<PhylogenyNode, RankTaxon> roots = TreePanelUtil.maximalMonochromaticRoots( t, sp );
        check( ok, roots.size() == 1 && roots.containsKey( t.getRoot() ) && "Felis catus".equals( roots.get( t.getRoot() ).getName() )
                && TreePanelUtil.unplacedTipCount( t, sp ) == 11, "species: the one placed tip makes the whole tree its clade" );
        check( ok, TreePanelUtil.maximalMonochromaticRoots( t, TreePanelUtil.assignNodesToRankTaxon( t, "family", null, TaxonSource.TREE ) ).isEmpty(),
               "a rank placing nothing gives no clade" );
        final Phylogeny pure = tree( clade( "r", null, tip( "x", prop( "a:family", "F" ) ),
                                            clade( "y", null, tip( "y1", prop( "b:Family", "f" ) ), tip( "y2" ) ) ) );
        final Map<PhylogenyNode, RankTaxon> fa = TreePanelUtil.assignNodesToRankTaxon( pure, "family", null, TaxonSource.TREE );
        final Map<PhylogenyNode, RankTaxon> pr = TreePanelUtil.maximalMonochromaticRoots( pure, fa );
        check( ok, pr.size() == 1 && pr.containsKey( pure.getRoot() ) && fa.size() == 2 && "F".equals( pr.get( pure.getRoot() ).getName() ),
               "a tree pure for one taxon is one clade at the root (3 tips, 2 placed, F and f one taxon)" );
        final Phylogeny nameless = tree( clade( "r", null, tip( "x", tax( "genus", null ) ),
                                                tip( "y", code( "genus", "CANFA" ) ), tip( "z", common( "genus", "dog" ) ) ) );
        final Map<PhylogenyNode, TaxonPlacement> np = TreePanelUtil.placeNodes( nameless, "genus", null, TaxonSource.TREE );
        check( ok, !np.get( node( nameless, "x" ) ).isPlaced() && placed( np, nameless, "y", "CANFA" ) && placed( np, nameless, "z", "dog" ),
               "a taxonomy without a usable name places nothing; a code or a common name does" );
    }

    // ---- what is offered: the repeat rule, the refused names, the file's own schemes, the default ----
    private static void offered( final boolean[] ok ) {
        final Phylogeny t = fixture();
        final List<RankSourcePanel.Coverage> cov = RankSourcePanel.coverage( t, null, TaxonSource.TREE );
        check( ok, cov.size() == 1 && "genus".equals( cov.get( 0 ).key ) && cov.get( 0 ).placed == 11 && cov.get( 0 ).tips == 12
                && cov.get( 0 ).taxa == 6, "genus offered (11 of 12 placed, 6 taxa); species places one tip of one taxon, not offered; other is no rank: "
                + describe( cov ) );
        // the strains repeat (two tips each), so the repeat rule would let them through: refused by name. (The JS
        // fixture also puts "Sub_Strain" and "isolate" into <taxonomy> ranks; phyloXML's rank vocabulary, which
        // Taxonomy.setRank enforces here, has neither, so those two are properties only on this side.)
        final Phylogeny strains = tree( clade( "r", null,
                tip( "s1", prop( "vipr:Strain", "A" ), prop( "x:sub-strain", "A1" ), prop( "y:Isolate", "I1" ), tax( "strain", "A" ) ),
                tip( "s2", prop( "vipr:Strain", "A" ), prop( "x:sub_strain", "A1" ), prop( "y:Isolate", "I1" ), tax( "strain", "A" ) ),
                tip( "s3", prop( "vipr:Strain", "B" ), prop( "x:substrain", "B1" ), prop( "y:Isolate", "I2" ), tax( "strain", "B" ) ),
                tip( "s4", prop( "vipr:Strain", "B" ), prop( "x:substrain", "B1" ), prop( "y:Isolate", "I2" ), tax( "other", "o" ) ) ) );
        check( ok, RankSourcePanel.coverage( strains, null, TaxonSource.TREE ).isEmpty(), "strain, substrain, isolate, other: nothing offered" );
        // a property is a rank by the list, or by naming a clade or lineage; a field is not
        final Phylogeny schemes = tree( clade( "r", null, schemeTip( "t1" ), schemeTip( "t2" ) ) );
        final List<RankSourcePanel.Coverage> sc = RankSourcePanel.coverage( schemes, null, TaxonSource.TREE );
        final List<String> keys = new ArrayList<String>(), labels = new ArrayList<String>();
        for( final RankSourcePanel.Coverage c : sc ) {
            keys.add( c.key );
            labels.add( c.label );
        }
        check( ok, String.join( " ", keys ).equals( "serotype subtype clademembership h5clade pangolineagel1" )
                && String.join( " ", labels ).equals( "serotype subtype clade_membership H5Clade PANGO_Lineage_L1" ),
               "listed ranks first in their order, then the file's own schemes alphabetically under the file's spelling; host_group and genome_type are fields: "
                + keys + " " + labels );
        check( ok, placed( TreePanelUtil.placeNodes( schemes, "H5-Clade", null, TaxonSource.TREE ), schemes, "t1", "2.3.4.4" ),
               "the file's own scheme places by its folded name" );
        // a scheme only INTERNAL nodes carry (a branch label, as the Auspice reader writes them) is no rank to offer;
        // a listed rank on an internal node is (it reaches the tips below)
        final Phylogeny internal_only = tree( clade( "r", null,
                clade( "i", prop( "ird:H5Clade", "2.3.4.4" ), tip( "u1" ), tip( "u2" ) ),
                clade( "j", prop( "x:family", "Fam" ), tip( "u3" ), tip( "u4" ) ) ) );
        final List<String> io = new ArrayList<String>();
        for( final RankSourcePanel.Coverage c : RankSourcePanel.coverage( internal_only, null, TaxonSource.TREE ) ) {
            io.add( c.key );
        }
        check( ok, io.equals( Arrays.asList( "family" ) ), "an internal-only scheme is not offered; an internal listed rank is: " + io );
        // the default: the finest listed rank naming 2..20 taxa, else the most taxa
        check( ok, "subtype".equals( RankSourcePanel.defaultRank( Arrays.asList( cv( "species", 3 ), cv( "subtype", 2 ), cv( "subclade", 30 ) ), 20 ) )
                && "genus".equals( RankSourcePanel.defaultRank( Arrays.asList( cv( "family", 1 ), cv( "genus", 4 ), cv( "species", 234 ) ), 20 ) )
                && "species".equals( RankSourcePanel.defaultRank( Arrays.asList( cv( "class", 1 ), cv( "family", 1 ), cv( "species", 63 ) ), 20 ) )
                && "family".equals( RankSourcePanel.defaultRank( Arrays.asList( cv( "family", 2 ), cv( "genus", 10 ), cv( "species", 66 ) ), 5 ) )
                && ( RankSourcePanel.defaultRank( new ArrayList<RankSourcePanel.Coverage>(), 20 ) == null ),
               "the default rank" );
    }

    // ---- the cross-run tree both programs ship ----
    private static void caliciviridae( final boolean[] ok ) throws Exception {
        final File f = new File( System.getProperty( "user.dir" ), "forester/test_data/rank/Caliciviridae_100.xml" );
        if ( !f.exists() ) {
            ok[ 0 ] = fail( "fixture missing: " + f );
            return;
        }
        final Phylogeny t = ParserBasedPhylogenyFactory.getInstance().create( f, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        final Map<PhylogenyNode, RankTaxon> g = TreePanelUtil.assignNodesToRankTaxon( t, "genus", null, TaxonSource.TREE );
        final Map<PhylogenyNode, RankTaxon> roots = TreePanelUtil.maximalMonochromaticRoots( t, g );
        int placed = 0;
        final Map<String, Integer> clades_of = new TreeMap<String, Integer>();
        final java.util.Set<String> taxa = new java.util.TreeSet<String>();
        for( final Entry<PhylogenyNode, RankTaxon> e : g.entrySet() ) {
            if ( e.getKey().isExternal() ) {
                ++placed;
                taxa.add( e.getValue().getName() );
            }
        }
        int single = 0;
        final List<Integer> nebo = new ArrayList<Integer>();
        for( final Entry<PhylogenyNode, RankTaxon> e : roots.entrySet() ) {
            clades_of.merge( e.getValue().getName(), 1, Integer::sum );
            if ( e.getKey().isExternal() && "Norovirus".equals( e.getValue().getName() ) ) {
                ++single;
            }
            if ( "Nebovirus".equals( e.getValue().getName() ) ) {
                nebo.add( e.getKey().isExternal() ? 1 : e.getKey().getAllExternalDescendants().size() );
            }
        }
        nebo.sort( null );
        check( ok, t.getNumberOfExternalNodes() == 97 && placed == 87 && taxa.size() == 10 && roots.size() == 17,
               "Caliciviridae at genus: 97 tips, 87 placed, 10 taxa, 17 clades under the no-vote rule; got " + t.getNumberOfExternalNodes()
                + "/" + placed + "/" + taxa.size() + "/" + roots.size() );
        check( ok, Integer.valueOf( 2 ).equals( clades_of.get( "Norovirus" ) ) && single == 1,
               "Norovirus is two clades, one of them a single tip: " + clades_of );
        // Nebovirus is FOUR clades, and the tree forces it: a Vesivirus tip (Bovine_calicivirus|Mukti) stands beside
        // the Kirklareli tip, so their parent is mixed; BoNeV16 is that mixed node's sibling, and the 3-tip "Bovine
        // nebovirus" species node is a third branch -- 9 (one swept-in tip), 1, 1 and 3 tips
        check( ok, nebo.equals( Arrays.asList( 1, 1, 3, 9 ) ), "Nebovirus: four clades of 9, 1, 1 and 3 tips, got " + nebo );
        check( ok, clades_of.equals( map( "Bavovirus", 1, "Lagovirus", 1, "Nacovirus", 1, "Nebovirus", 4, "Norovirus", 2, "Recovirus", 1,
                                         "Salovirus", 1, "Sapovirus", 3, "Valovirus", 1, "Vesivirus", 2 ) ),
               "the clades per genus: " + clades_of );
        final CladeAnnotationReport fam = CladeAnnotationReport.of( t, Arrays.asList( "family" ), null, TaxonSource.TREE );
        final java.util.Set<String> fam_taxa = new java.util.TreeSet<String>();
        for( final RankTaxon rt : TreePanelUtil.assignNodesToRankTaxon( t, "family", null, TaxonSource.TREE ).values() ) {
            fam_taxa.add( rt.getName() );
        }
        long internal = fam.getRows().stream().filter( CladeAnnotationReport.Row::isInternal ).count();
        check( ok, fam_taxa.size() == 2 && internal == 1 && fam.getRows().stream().anyMatch( r -> r.isInternal() && "Caliciviridae".equals( r.getTaxon() )
                && r.getNote().contains( "unclassified Caliciviridae" ) ),
               "at family: 2 taxa and one conflict, the root's Caliciviridae against a tip's \"unclassified Caliciviridae\" (a different string is a different taxon): "
                + fam_taxa + ", " + internal + " internal" );
    }

    // ---- fixture ----

    static Phylogeny fixture() {
        final PhylogenyNode a = clade( "A", tax( "genus", "Felis" ), tip( "a1", tax( "species", "Felis catus" ) ), tip( "a2", prop( "vipr:Genus", "felis" ) ), tip( "a3" ) );
        final PhylogenyNode b = clade( "B", null, tip( "b1", prop( "repseq:genus", "Panthera" ) ), tip( "b2" ), tip( "b3", prop( "x:Genus", "  Panthera  " ) ) );
        final PhylogenyNode c = clade( "C", null, tip( "c1", prop( "vipr:Genus", "Felis" ) ), tip( "c2", prop( "vipr:Genus", "Lynx" ) ) );
        final PhylogenyNode d = clade( "D", tax( "genus", "Canis" ), tip( "d1" ), tip( "d2", prop( "foo:Host", "dog" ) ) );
        final PhylogenyNode e = clade( "E", tax( "genus", "Vulpes" ), tip( "e1", tax( "genus", "Nyctereutes" ), prop( "vipr:Genus", "Vulpes" ) ),
                                       tip( "e2", prop( "vipr:Genus", "Vulpes" ), tax( "other", "x" ) ) );
        return tree( clade( "root", null, a, b, c, d, e ) );
    }

    static PhylogenyNode schemeTip( final String name ) {
        return tip( name, prop( "ird:H5Clade", "2.3.4.4" ), prop( "vipr:PANGO_Lineage_L1", "B.1" ), prop( "ns:clade_membership", "20A" ),
                    prop( "BVBRC:host_group", "Human" ), prop( "x:genome_type", "ssRNA" ), prop( "ird:Subtype", "H5N1" ), prop( "q:Serotype", "O1" ) );
    }

    /** A datum to attach: a taxonomy or a property. */
    interface Datum {
        void attach( PhylogenyNode n );
    }

    static Datum tax( final String rank, final String scientific_name ) {
        return n -> taxonomy( n, rank, scientific_name, null, null );
    }

    static Datum code( final String rank, final String code ) {
        return n -> taxonomy( n, rank, null, code, null );
    }

    static Datum common( final String rank, final String common_name ) {
        return n -> taxonomy( n, rank, null, null, common_name );
    }

    static void taxonomy( final PhylogenyNode n, final String rank, final String sci, final String code, final String common ) {
        try {
            final Taxonomy t = new Taxonomy();
            if ( sci != null ) {
                t.setScientificName( sci );
            }
            if ( code != null ) {
                t.setTaxonomyCode( code );
            }
            if ( common != null ) {
                t.setCommonName( common );
            }
            t.setRank( rank );
            n.getNodeData().addTaxonomy( t );
        }
        catch ( final org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException e ) {
            throw new IllegalArgumentException( e );
        }
    }

    static Datum prop( final String ref, final String value ) {
        return n -> {
            if ( n.getNodeData().getProperties() == null ) {
                n.getNodeData().setProperties( new PropertiesList() );
            }
            n.getNodeData().getProperties().addProperty( new Property( ref, value, "", "xsd:string", AppliesTo.NODE ) );
        };
    }

    static PhylogenyNode tip( final String name, final Datum... data ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( 1 );
        for( final Datum d : data ) {
            d.attach( n );
        }
        return n;
    }

    static PhylogenyNode clade( final String name, final Datum own, final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( 1 );
        if ( own != null ) {
            own.attach( n );
        }
        for( final PhylogenyNode c : children ) {
            n.addAsChild( c );
        }
        return n;
    }

    static Phylogeny tree( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    static PhylogenyNode node( final Phylogeny t, final String name ) {
        return t.getNode( name );
    }

    static boolean placed( final Map<PhylogenyNode, TaxonPlacement> p, final Phylogeny t, final String name, final String taxon ) {
        final TaxonPlacement tp = p.get( node( t, name ) );
        return ( tp != null ) && tp.isPlaced() && taxon.equals( tp.getTaxon().getName() );
    }

    static RankSourcePanel.Coverage cv( final String rank, final int taxa ) {
        return new RankSourcePanel.Coverage( RankVocabulary.key( rank ), rank, 0, taxa, taxa, taxa, 0 );
    }

    static String describe( final List<RankSourcePanel.Coverage> cov ) {
        final StringBuilder sb = new StringBuilder();
        for( final RankSourcePanel.Coverage c : cov ) {
            sb.append( c.label ).append( '(' ).append( c.placed ).append( '/' ).append( c.tips ).append( ';' ).append( c.taxa ).append( ") " );
        }
        return sb.toString();
    }

    static Map<String, Integer> map( final Object... kv ) {
        final Map<String, Integer> m = new TreeMap<String, Integer>();
        for( int i = 0; i < kv.length; i += 2 ) {
            m.put( (String) kv[ i ], (Integer) kv[ i + 1 ] );
        }
        return m;
    }

    private static void check( final boolean[] ok, final boolean c, final String what ) {
        if ( !c ) {
            System.out.println( "  [RankCladesTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private static boolean fail( final String what ) {
        System.out.println( "  [RankCladesTest] " + what );
        return false;
    }

    private RankCladesTest() {
    }
}
