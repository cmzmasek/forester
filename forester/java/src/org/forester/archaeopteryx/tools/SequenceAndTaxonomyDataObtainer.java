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

import java.io.IOException;
import java.net.UnknownHostException;
import java.util.SortedSet;

import javax.swing.JOptionPane;

import org.forester.analysis.AncestralTaxonomyInferenceException;
import org.forester.analysis.TaxonomyDataManager;
import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;
import org.forester.ws.seqdb.SequenceDbWsTools;

/**
 * Single Tools-menu operation that enriches a tree from external web services in one pass:
 * it first obtains sequence information (UniProt / EMBL-GenBank) and then resolves detailed
 * taxonomic information (UniProt Taxonomy). This merges what used to be the two separate
 * "Obtain Sequence Information" and "Obtain Detailed Taxonomic Information" tools into one.
 *
 * <p>Both phases run best-effort on the same phylogeny copy in a single background thread:
 * a failure or unresolved node in one phase does not abort the other, and whatever could be
 * obtained is committed to the displayed tree, with a single combined report dialog at the end.
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
        execute();
    }

    private void execute() {
        start( _mf, "sequence & taxonomy data" );
        SortedSet<String> seq_not_found = null;
        SortedSet<String> tax_not_found = null;
        String seq_error = null;
        String tax_error = null;
        try {
            try {
                seq_not_found = SequenceDbWsTools.obtainSeqInformation( _phy,
                                                                       false,
                                                                       true,
                                                                       SequenceDbWsTools.DEFAULT_LINES_TO_RETURN );
            }
            catch ( final UnknownHostException e ) {
                seq_error = "could not connect to the sequence database";
            }
            catch ( final IOException e ) {
                seq_error = e.getLocalizedMessage();
            }
            try {
                tax_not_found = TaxonomyDataManager.obtainDetailedTaxonomicInformation( _phy, true );
            }
            catch ( final UnknownHostException e ) {
                tax_error = "could not connect to the taxonomy database";
            }
            catch ( final IOException e ) {
                tax_error = e.getLocalizedMessage();
            }
            catch ( final AncestralTaxonomyInferenceException e ) {
                tax_error = e.getLocalizedMessage();
            }
        }
        finally {
            end( _mf );
        }
        // Only commit (replace the displayed tree and mark it edited) when at least one phase
        // completed without error; on a total failure nothing was obtained, so leave the tree
        // and its edited-state untouched -- just report the errors below.
        if ( shouldCommit( seq_error, tax_error ) ) {
            _treepanel.setTree( _phy );
            _mf.showWhole();
            _treepanel.setEdited( true );
        }
        final String message = buildCompletionMessage( seq_not_found, seq_error, tax_not_found, tax_error );
        final int type = hasIssues( seq_not_found, seq_error, tax_not_found, tax_error ) ? JOptionPane.WARNING_MESSAGE
                : JOptionPane.INFORMATION_MESSAGE;
        try {
            JOptionPane.showMessageDialog( _mf, message, "Sequence & Taxonomy Tool Completed", type );
        }
        catch ( final Exception e ) {
            // Not important if this fails, do nothing.
        }
    }

    /**
     * True if at least one phase completed without error, so whatever it obtained should be
     * committed to the displayed tree. When both phases fail, nothing was obtained and the
     * tree (and its edited-state) is left untouched.
     */
    static boolean shouldCommit( final String seq_error, final String tax_error ) {
        return ( seq_error == null ) || ( tax_error == null );
    }

    /** True if either phase produced an error or left at least one node unresolved. */
    static boolean hasIssues( final SortedSet<String> seq_not_found,
                              final String seq_error,
                              final SortedSet<String> tax_not_found,
                              final String tax_error ) {
        return ( seq_error != null ) || ( tax_error != null ) || ( ( seq_not_found != null ) && !seq_not_found.isEmpty() )
                || ( ( tax_not_found != null ) && !tax_not_found.isEmpty() );
    }

    /** Builds the combined report shown after both phases have run. Never returns null. */
    static String buildCompletionMessage( final SortedSet<String> seq_not_found,
                                          final String seq_error,
                                          final SortedSet<String> tax_not_found,
                                          final String tax_error ) {
        if ( !hasIssues( seq_not_found, seq_error, tax_not_found, tax_error ) ) {
            return "Sequence and taxonomy data successfully obtained.";
        }
        final StringBuilder sb = new StringBuilder();
        appendSection( sb, "Sequence data", "sequence information", seq_not_found, seq_error );
        appendSection( sb, "Taxonomy data", "taxonomy", tax_not_found, tax_error );
        return sb.toString().trim();
    }

    private static void appendSection( final StringBuilder sb,
                                       final String title,
                                       final String what,
                                       final SortedSet<String> not_found,
                                       final String error ) {
        if ( error != null ) {
            sb.append( title ).append( ": error - " ).append( error ).append( "\n\n" );
        }
        else if ( ( not_found != null ) && !not_found.isEmpty() ) {
            sb.append( title ).append( ": no " ).append( what );
            if ( not_found.size() == 1 ) {
                sb.append( " was found for the following node:\n" );
            }
            else {
                sb.append( " was found for the following nodes (total: " ).append( not_found.size() ).append( "):\n" );
            }
            int i = 0;
            for( final String s : not_found ) {
                if ( i >= MAX_NODES_TO_LIST ) {
                    sb.append( "...\n" );
                    break;
                }
                sb.append( s ).append( "\n" );
                ++i;
            }
            sb.append( "\n" );
        }
    }
}
