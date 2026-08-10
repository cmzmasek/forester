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

import org.forester.archaeopteryx.MainFrame;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

/**
 * The undoable "install" half of the Analysis-menu "Infer Ancestor Taxonomies" operation. The taxa themselves
 * are computed by the pure {@link org.forester.analysis.AncestralTaxonomyInference} engine (driven from
 * {@code MainFrame.executeLineageInference}, which also handles the online resolve); this class just commits the
 * enriched tree to the panel so the commit path stays testable in isolation, without running the (possibly
 * online) resolve. EDT-only.
 */
public final class AncestralTaxonomyInferrer {

    private final MainFrame  _mf;
    private final TreePanel  _treepanel;
    private final Phylogeny  _phy;
    private final boolean    _overwrite;

    public AncestralTaxonomyInferrer( final MainFrame mf,
                                      final TreePanel treepanel,
                                      final Phylogeny phy,
                                      final boolean overwrite ) {
        _mf = mf;
        _treepanel = treepanel;
        _phy = phy;
        _overwrite = overwrite;
    }

    /**
     * EDT-only: checkpoint the live (pre-inference) tree so the inferred taxonomies can be undone, append a
     * provenance sentence to the description, then install the enriched tree. {@code assigned} is the number of
     * internal nodes that received a taxonomy (for the provenance sentence).
     */
    public void commit( final int assigned ) {
        _treepanel.pushUndoCheckpoint( "Infer Ancestor Taxonomies" );
        final String sentence = inferenceProvenance( _phy, assigned, _overwrite );
        final String existing = _phy.getDescription();
        _phy.setDescription( ForesterUtil.isEmpty( existing ) ? sentence : existing + " " + sentence );
        _phy.setRerootable( false );
        _treepanel.setTree( _phy );
        _mf.showWhole();
        _treepanel.setEdited( true );
    }

    /** The user-facing completion summary (pure/testable). Grammatically singular/plural-correct. */
    public static String inferenceSummary( final int assigned, final int skipped_existing, final int unresolved_tips ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( ( assigned == 1 ) ? "Assigned a taxonomy to 1 internal node."
                : "Assigned taxonomies to " + assigned + " internal nodes." );
        if ( skipped_existing > 0 ) {
            sb.append( "\nKept " ).append( skipped_existing ).append( " existing internal-node taxonom" )
                    .append( ( skipped_existing == 1 ) ? "y" : "ies" ).append( " (enable overwrite to replace)." );
        }
        if ( unresolved_tips > 0 ) {
            sb.append( "\n" ).append( unresolved_tips ).append( " tip" ).append( ( unresolved_tips == 1 ) ? "" : "s" )
                    .append( " had no ancestor lineage, so some internal nodes could not be inferred." );
        }
        return sb.toString();
    }

    /** The provenance sentence appended to the tree description (pure/testable). Append, never overwrite. */
    public static String inferenceProvenance( final Phylogeny phy, final int assigned, final boolean overwrite ) {
        final String name = ForesterUtil.isEmpty( phy.getName() ) ? "" : phy.getName();
        return "Used ancestral-taxonomy inference (deepest shared descendant lineage"
                + ( overwrite ? ", overwriting existing internal taxa" : "" ) + ") to assign taxonomies to "
                + assigned + " internal node" + ( assigned == 1 ? "" : "s" ) + " in tree named \"" + name
                + "\" with " + phy.getNumberOfExternalNodes() + " tips.";
    }
}
