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

import org.forester.archaeopteryx.tools.AlignmentImporter;
import org.forester.io.parsers.FastaParser;
import org.forester.msa.Msa;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/** Headless tests for joining an aligned FASTA onto tips by name ({@link AlignmentImporter}). */
public final class AlignmentImporterTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "AlignmentImporter: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            // tips A, B, C; alignment for A, B and an unmatched row X -> C stays without a sequence
            final Phylogeny phy = threeTips();
            final Msa msa = FastaParser.parseMsa( ">A\nMK-L\n>B\nMKQL\n>X\nQQQQ\n" );
            final AlignmentImporter.Result r = AlignmentImporter.apply( phy, msa );

            if ( r.getTipsAligned() != 2 ) {
                return fail( "two tips (A, B) should be aligned, got " + r.getTipsAligned() );
            }
            if ( r.getUnmatchedRows() != 1 ) {
                return fail( "one alignment row (X) should be unmatched, got " + r.getUnmatchedRows() );
            }
            if ( r.getTipsWithoutSequence() != 1 ) {
                return fail( "one tip (C) should have no sequence, got " + r.getTipsWithoutSequence() );
            }
            if ( r.getAlignmentLength() != 4 ) {
                return fail( "alignment length should be 4, got " + r.getAlignmentLength() );
            }
            final PhylogenyNode a = byName( phy, "A" );
            if ( !a.getNodeData().isHasSequence() || !"MK-L".equals( a.getNodeData().getSequence().getMolecularSequence() )
                    || !a.getNodeData().getSequence().isMolecularSequenceAligned() ) {
                return fail( "tip A must carry its aligned row 'MK-L' with the aligned flag set" );
            }
            if ( byName( phy, "C" ).getNodeData().isHasSequence() ) {
                return fail( "tip C (no matching row) must have no sequence" );
            }
            if ( !AptxUtil.hasAlignedSequences( phy ) ) {
                return fail( "hasAlignedSequences must be true after the import" );
            }

            // an UNEQUAL-length FASTA is not an alignment -> the parser rejects it
            boolean rejected = false;
            try {
                FastaParser.parseMsa( ">A\nMKL\n>B\nMK\n" );
            }
            catch ( final Exception ex ) {
                rejected = true;
            }
            if ( !rejected ) {
                return fail( "an unequal-length FASTA must be rejected as not an alignment" );
            }
            return true;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static Phylogeny threeTips() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( named( "A" ) );
        root.addAsChild( named( "B" ) );
        root.addAsChild( named( "C" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode named( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode byName( final Phylogeny phy, final String name ) {
        for( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( name.equals( ext.getName() ) ) {
                return ext;
            }
        }
        return null;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [AlignmentImporterTest] " + message );
        return false;
    }

    private AlignmentImporterTest() {
        // not instantiable
    }
}
