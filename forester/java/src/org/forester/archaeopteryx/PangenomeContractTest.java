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
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.forester.archaeopteryx.tools.NodeDataImporter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The JOINT pangenome fixture: a 100-strain x 40-gene presence/absence matrix, shared byte-for-byte with
 * Archaeopteryx.js so both programs can be measured on identical input.
 * <p>
 * Each cell is an ordinal certainty 0..4 -- evidence strength about a hidden binary truth, not a count -- and a
 * blank cell means NOT ASSESSED, which is not the same as 0. Gene presence was simulated by gain/loss on the tree,
 * so the clade-specific genes form real blocks and the 40 genes deliberately differ in distinctness (core /
 * clade-specific / patchy / rare / ambiguous). That heterogeneity is the point: it gives the 40 fields 40 distinct
 * visualization scores, so a ranking cannot be confused with an alphabetical tiebreak.
 * <p>
 * Three things are pinned here:
 * <ol>
 * <li><b>The table join</b> -- the tree plus the TSV through {@link NodeDataImporter}, the path a user walks, with
 * the missing cells left unfilled rather than zero-filled.</li>
 * <li><b>The visualization candidates</b> -- rank, ref, label, tier and score against
 * {@code desktop_candidates.tsv}. Archaeopteryx.js reproduced this table independently, in a different language,
 * agreeing on all 40 fields to under 1e-6; that agreement is the strongest cross-program evidence either side has,
 * and it is only worth anything while both sides read the SAME bytes. Do not regenerate these files.</li>
 * <li><b>The column layout order</b> -- {@code View > Clustergram} lays the matrix out in the TSV's own column
 * order, not in candidate order (see {@link TreePanelUtil#propertyRefsInSourceOrder}).</li>
 * </ol>
 * The fixture files are NOT reproducible from the generator that made them (it drew a random number per element of
 * a set, whose iteration order is hash-randomized per process), which is exactly why the artifact is pinned rather
 * than a recipe.
 */
public final class PangenomeContractTest {

    private static final String DIR       = "forester/test_data/pangenome_contract/";
    private static final String TREE      = "pangenome_tree.nwk";
    private static final String MATRIX    = "pangenome_matrix.tsv";
    private static final String EXPECTED  = "desktop_candidates.tsv";
    private static final int    TIPS      = 100;
    private static final int    GENES     = 40;
    /** Non-empty cells: 4000 minus the 169 deliberately unassessed ones. */
    private static final int    FILLED    = 3831;
    private static final double SCORE_EPS = 1.0e-6;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "PangenomeContract: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final Phylogeny phy = loadJoined();
        if ( phy == null ) {
            return false;
        }
        return joinOk( phy ) && candidatesMatchExpected( phy ) && columnsFollowTableOrder( phy );
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [PangenomeContractTest] " + msg );
        return false;
    }

    private static File file( final String name ) {
        return new File( System.getProperty( "user.dir" ), DIR + name );
    }

    /** The tree as the application opens it, with the TSV joined onto it exactly as the import tool would. */
    private static Phylogeny loadJoined() {
        try {
            for( final String n : new String[] { TREE, MATRIX, EXPECTED } ) {
                if ( !file( n ).exists() ) {
                    fail( "the joint fixture is missing: " + file( n ).getAbsolutePath() );
                    return null;
                }
            }
            final Phylogeny[] phys = FigureRenderer.readTrees( file( TREE ) );
            if ( ( phys.length != 1 ) || ( phys[ 0 ] == null ) || phys[ 0 ].isEmpty() ) {
                fail( TREE + " did not yield exactly one non-empty tree" );
                return null;
            }
            final Phylogeny phy = phys[ 0 ];
            final String tsv = Files.readString( file( MATRIX ).toPath() );
            final NodeDataImporter.Table table = NodeDataImporter.parseTable( tsv );
            // key column 0 is "strain"; matching by tip name is what a user picks in the dialog
            final NodeDataImporter.ImportResult res = NodeDataImporter
                    .apply( phy, table, 0, NodeDataImporter.MatchBy.TIP_NAME );
            if ( res.getTipsAnnotated() != TIPS ) {
                fail( "the join should annotate all " + TIPS + " tips, got " + res.getTipsAnnotated() );
                return null;
            }
            if ( res.getRowsMatched() != TIPS ) {
                fail( "every table row should match a tip, got " + res.getRowsMatched() );
                return null;
            }
            if ( res.getPropertyColumns().size() != GENES ) {
                fail( "the join should bring " + GENES + " gene columns, got " + res.getPropertyColumns().size() );
                return null;
            }
            return phy;
        }
        catch ( final Exception e ) {
            fail( "loading the joint fixture threw " + e );
            return null;
        }
    }

    /**
     * The join itself: every tip annotated, and -- the part that matters -- a blank cell left UNFILLED rather than
     * turned into a 0. 169 of the 4000 cells are deliberately unassessed, including one draft genome missing most
     * of its genes; if they became zeros the matrix would claim 169 confident absences it has no evidence for.
     */
    private static boolean joinOk( final Phylogeny phy ) {
        if ( phy.getNumberOfExternalNodes() != TIPS ) {
            return fail( "the fixture tree should have " + TIPS + " tips, got " + phy.getNumberOfExternalNodes() );
        }
        int filled = 0;
        int sparsest = GENES;
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( ( n.getNodeData() == null ) || ( n.getNodeData().getProperties() == null ) ) {
                return fail( "tip " + n.getName() + " carries no properties after the join" );
            }
            int here = 0;
            for( final org.forester.phylogeny.data.Property p : n.getNodeData().getProperties().getProperties() ) {
                if ( p.getRef().startsWith( "meta:" ) ) {
                    ++here;
                }
            }
            filled += here;
            sparsest = Math.min( sparsest, here );
        }
        if ( filled != FILLED ) {
            return fail( "a blank cell must stay UNFILLED, not become a 0: expected " + FILLED
                    + " gene properties across the tree, got " + filled );
        }
        if ( sparsest >= GENES ) {
            return fail( "the draft genome should carry FEWER than " + GENES
                    + " genes -- without it the fixture cannot show that missing is not zero" );
        }
        return true;
    }

    /**
     * Rank, ref, label, tier and score against the pinned table. Archaeopteryx.js matches this same table
     * independently; a change here is a change to a cross-program contract, not just to a local expectation.
     */
    private static boolean candidatesMatchExpected( final Phylogeny phy ) {
        final List<String[]> expected = new ArrayList<String[]>();
        try {
            final List<String> lines = Files.readAllLines( file( EXPECTED ).toPath() );
            for( int i = 1; i < lines.size(); ++i ) { // row 0 is the header
                final String line = lines.get( i );
                if ( line.trim().length() > 0 ) {
                    expected.add( line.split( "\t" ) );
                }
            }
        }
        catch ( final Exception e ) {
            return fail( "could not read " + EXPECTED + ": " + e );
        }
        if ( expected.size() != GENES ) {
            return fail( EXPECTED + " should hold " + GENES + " rows, got " + expected.size() );
        }
        // GUARD: a table that happened to be alphabetical by label could not tell a ranking from a sort, so this
        // fixture would pin nothing. Assert it is NOT alphabetical before trusting it.
        boolean alphabetical = true;
        for( int i = 1; i < expected.size(); ++i ) {
            if ( expected.get( i - 1 )[ 2 ].compareToIgnoreCase( expected.get( i )[ 2 ] ) > 0 ) {
                alphabetical = false;
                break;
            }
        }
        if ( alphabetical ) {
            return fail( "the expected table is in alphabetical label order, so it cannot distinguish the"
                    + " candidate RANKING from a sort -- the fixture would pin nothing" );
        }
        final List<PropertyColorScheme.VisCandidate> got = PropertyColorScheme.visualizationCandidates( phy );
        if ( got.size() != expected.size() ) {
            return fail( "expected " + expected.size() + " candidates, got " + got.size() );
        }
        for( int i = 0; i < expected.size(); ++i ) {
            final String[] want = expected.get( i );
            final PropertyColorScheme.VisCandidate c = got.get( i );
            if ( !want[ 1 ].equals( c._ref ) ) {
                return fail( "rank " + i + ": expected ref " + want[ 1 ] + ", got " + c._ref );
            }
            if ( !want[ 2 ].equals( c._label ) ) {
                return fail( "rank " + i + " (" + c._ref + "): expected label \"" + want[ 2 ] + "\", got \""
                        + c._label + "\"" );
            }
            if ( Integer.parseInt( want[ 3 ] ) != c.tier() ) {
                return fail( "rank " + i + " (" + c._ref + "): expected tier " + want[ 3 ] + ", got " + c.tier() );
            }
            final double want_score = Double.parseDouble( want[ 4 ] );
            if ( Math.abs( want_score - c._score ) > SCORE_EPS ) {
                return fail( "rank " + i + " (" + c._ref + "): expected score " + want_score + ", got " + c._score );
            }
        }
        return true;
    }

    /**
     * The matrix lays out in the TSV's own column order. The expected candidate table is ranked by score, so its
     * order is demonstrably NOT the table's order -- which is what makes this assertion worth making.
     */
    private static boolean columnsFollowTableOrder( final Phylogeny phy ) {
        final List<String> headers = new ArrayList<String>();
        try {
            final String first = Files.readAllLines( file( MATRIX ).toPath() ).get( 0 );
            final String[] cols = first.split( "\t" );
            for( int i = 1; i < cols.length; ++i ) { // column 0 is the "strain" key
                headers.add( "meta:" + cols[ i ] );
            }
        }
        catch ( final Exception e ) {
            return fail( "could not read the header of " + MATRIX + ": " + e );
        }
        if ( headers.size() != GENES ) {
            return fail( MATRIX + " should have " + GENES + " gene columns, got " + headers.size() );
        }
        final List<AnnotationColumns.ColumnSpec> specs = MainFrame.clustergramColumnSpecs( phy );
        final List<String> matrix = new ArrayList<String>();
        for( final AnnotationColumns.ColumnSpec s : specs ) {
            if ( s._type == AnnotationColumns.Type.MATRIX ) {
                matrix.add( s._ref );
            }
        }
        if ( matrix.size() != GENES ) {
            return fail( "all " + GENES + " numeric gene fields should become MATRIX columns, got " + matrix.size() );
        }
        if ( !String.join( " ", headers ).equals( String.join( " ", matrix ) ) ) {
            return fail( "the matrix must lay out in the TABLE's column order.\n    expected: "
                    + String.join( " ", headers ).substring( 0, 90 ) + "...\n    got     : "
                    + String.join( " ", matrix ).substring( 0, 90 ) + "..." );
        }
        // and the neighbouring wrong answer, named: candidate order is a real, different order for this fixture
        final List<String> candidate_order = PropertyColorScheme.colorableRefs( phy );
        if ( String.join( " ", candidate_order ).equals( String.join( " ", headers ) ) ) {
            return fail( "candidate order equals table order for this fixture, so the assertion above proves"
                    + " nothing -- the fixture has lost the property that made it useful" );
        }
        return true;
    }
}
