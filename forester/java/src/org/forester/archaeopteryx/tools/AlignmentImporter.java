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

package org.forester.archaeopteryx.tools;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.forester.msa.Msa;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;

/**
 * Joins a multiple sequence alignment ({@link Msa}, e.g. parsed from an aligned FASTA) onto a tree's external tips by
 * NAME, writing each matched tip's gapped row as its molecular sequence (aligned). So the alignment becomes part of the
 * tree data -- it round-trips to phyloXML and the display reads it back per tip. Pure and GUI-free (headless-testable);
 * mutates {@code phy} in place, so the caller should snapshot for undo BEFORE calling {@link #apply}.
 */
public final class AlignmentImporter {

    /** Outcome counts for the result dialog. */
    public static final class Result {

        private final int _tips_aligned;
        private final int _unmatched_rows;
        private final int _tips_without_sequence;
        private final int _alignment_length;

        Result( final int tips_aligned, final int unmatched_rows, final int tips_without_sequence,
                final int alignment_length ) {
            _tips_aligned = tips_aligned;
            _unmatched_rows = unmatched_rows;
            _tips_without_sequence = tips_without_sequence;
            _alignment_length = alignment_length;
        }

        public int getTipsAligned() {
            return _tips_aligned;
        }

        public int getUnmatchedRows() {
            return _unmatched_rows;
        }

        public int getTipsWithoutSequence() {
            return _tips_without_sequence;
        }

        public int getAlignmentLength() {
            return _alignment_length;
        }

        public String summary() {
            final StringBuilder sb = new StringBuilder();
            sb.append( "Aligned " ).append( _tips_aligned ).append( _tips_aligned == 1 ? " tip" : " tips" );
            sb.append( " (" ).append( _alignment_length ).append( " columns)." );
            if ( _unmatched_rows > 0 ) {
                sb.append( "\n" ).append( _unmatched_rows )
                        .append( _unmatched_rows == 1 ? " alignment sequence had" : " alignment sequences had" )
                        .append( " no matching tip (by name)." );
            }
            if ( _tips_without_sequence > 0 ) {
                sb.append( "\n" ).append( _tips_without_sequence )
                        .append( _tips_without_sequence == 1 ? " tip has" : " tips have" )
                        .append( " no alignment sequence." );
            }
            return sb.toString();
        }
    }

    /**
     * Writes each MSA row onto the external tip whose (trimmed) name matches the row identifier, as an aligned
     * molecular sequence. A row matching several same-named tips annotates all of them; a blank name never matches.
     */
    public static Result apply( final Phylogeny phy, final Msa msa ) {
        final Map<String, List<PhylogenyNode>> by_name = indexTipsByName( phy );
        final Set<PhylogenyNode> aligned = Collections.newSetFromMap( new IdentityHashMap<PhylogenyNode, Boolean>() );
        int unmatched_rows = 0;
        for( int row = 0; row < msa.getNumberOfSequences(); ++row ) {
            final String id = trimmed( msa.getIdentifier( row ) );
            final List<PhylogenyNode> tips = id.isEmpty() ? null : by_name.get( id );
            if ( ( tips == null ) || tips.isEmpty() ) {
                ++unmatched_rows;
                continue;
            }
            final String seq = msa.getSequenceAsString( row ).toString();
            for( final PhylogenyNode tip : tips ) {
                writeAlignedSequence( tip, seq );
                aligned.add( tip );
            }
        }
        final int tips_without = phy.getNumberOfExternalNodes() - aligned.size();
        return new Result( aligned.size(), unmatched_rows, tips_without, msa.getLength() );
    }

    private static Map<String, List<PhylogenyNode>> indexTipsByName( final Phylogeny phy ) {
        final Map<String, List<PhylogenyNode>> by_name = new HashMap<String, List<PhylogenyNode>>();
        for( final PhylogenyNode ext : phy.getExternalNodes() ) {
            final String key = trimmed( ext.getName() );
            if ( !key.isEmpty() ) {
                List<PhylogenyNode> list = by_name.get( key );
                if ( list == null ) {
                    list = new ArrayList<PhylogenyNode>();
                    by_name.put( key, list );
                }
                list.add( ext );
            }
        }
        return by_name;
    }

    private static void writeAlignedSequence( final PhylogenyNode tip, final String seq ) {
        Sequence sequence;
        if ( tip.getNodeData().isHasSequence() ) {
            sequence = tip.getNodeData().getSequence();
        }
        else {
            sequence = new Sequence();
            tip.getNodeData().addSequence( sequence );
        }
        sequence.setMolecularSequence( seq );
        sequence.setMolecularSequenceAligned( true );
    }

    private static String trimmed( final String s ) {
        return ( s == null ) ? "" : s.trim();
    }

    private AlignmentImporter() {
        // not instantiable
    }
}
