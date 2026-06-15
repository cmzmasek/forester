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

import java.util.SortedSet;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;
import org.forester.ws.seqdb.NcbiTaxonomyLineageService;
import org.forester.ws.seqdb.SequenceTaxonomyResolver;
import org.forester.ws.seqdb.WebSequenceFetcher;

/**
 * The Tools-menu "Fetch Sequence &amp; Taxonomic Data" operation: enriches the current tree from
 * external web services in one background pass. It runs the clean, timeout/cancel-aware
 * {@link SequenceTaxonomyResolver} (UniProt REST + NCBI E-utilities) and commits + reports on the EDT.
 *
 * <p>As a {@link RunnableProcess} it is a {@code CancelFlag}: the user can cancel it from the process
 * menu, and the resolver stops promptly at its next per-node check.
 */
public final class SequenceAndTaxonomyDataObtainer extends RunnableProcess {

    private static final int           MAX_NODES_TO_LIST = 20;
    private final Phylogeny            _phy;
    private final MainFrameApplication _mf;
    private final TreePanel            _treepanel;

    public SequenceAndTaxonomyDataObtainer( final MainFrameApplication mf,
                                            final TreePanel treepanel,
                                            final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    @Override
    public void run() {
        start( _mf, "sequence & taxonomy data" );
        final SequenceTaxonomyResolver.Result result;
        try {
            final SequenceTaxonomyResolver resolver = new SequenceTaxonomyResolver( new WebSequenceFetcher(),
                                                                                    new NcbiTaxonomyLineageService() );
            result = resolver.resolve( _phy, this ); // `this` is the CancelFlag
        }
        finally {
            end( _mf );
        }
        // Swing is not thread-safe: commit the enriched tree and show the report on the EDT.
        SwingUtilities.invokeLater( () -> {
            if ( shouldCommit( result ) ) {
                _treepanel.setTree( _phy );
                _mf.showWhole();
                _treepanel.setEdited( true );
            }
            final int type = hasIssues( result ) ? JOptionPane.WARNING_MESSAGE : JOptionPane.INFORMATION_MESSAGE;
            try {
                JOptionPane.showMessageDialog( _mf, buildCompletionMessage( result ), "Sequence & Taxonomy Tool Completed",
                                               type );
            }
            catch ( final Exception e ) {
                // not important if the dialog fails
            }
        } );
    }

    /** Commit the enriched tree whenever something was written (even a partial/cancelled run). */
    static boolean shouldCommit( final SequenceTaxonomyResolver.Result r ) {
        return ( r.getSequencesWritten() + r.getTaxonomiesWritten() ) > 0;
    }

    /** True if the run was cancelled, errored, or left any node unresolved. */
    static boolean hasIssues( final SequenceTaxonomyResolver.Result r ) {
        return r.wasCancelled() || ( r.getError() != null ) || !r.getNotFound().isEmpty();
    }

    /** Builds the report shown after the run. Never returns null. */
    static String buildCompletionMessage( final SequenceTaxonomyResolver.Result r ) {
        if ( !hasIssues( r ) ) {
            return "Sequence and taxonomy data successfully obtained (" + r.getSequencesWritten() + " sequence(s), "
                    + r.getTaxonomiesWritten() + " taxonomy/taxonomies).";
        }
        final StringBuilder sb = new StringBuilder();
        if ( r.wasCancelled() ) {
            sb.append( "Cancelled before completion.\n" );
        }
        if ( r.getError() != null ) {
            sb.append( "Error: " ).append( r.getError() ).append( "\n" );
        }
        sb.append( "Obtained " ).append( r.getSequencesWritten() ).append( " sequence(s) and " )
                .append( r.getTaxonomiesWritten() ).append( " taxonomy/taxonomies.\n" );
        final SortedSet<String> not_found = r.getNotFound();
        if ( !not_found.isEmpty() ) {
            sb.append( "\nNo data was found for the following " );
            sb.append( ( not_found.size() == 1 ) ? "node:\n" : ( "nodes (total: " + not_found.size() + "):\n" ) );
            int i = 0;
            for( final String s : not_found ) {
                if ( i >= MAX_NODES_TO_LIST ) {
                    sb.append( "...\n" );
                    break;
                }
                sb.append( s ).append( "\n" );
                ++i;
            }
        }
        return sb.toString().trim();
    }
}
