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

import java.io.IOException;
import java.net.UnknownHostException;
import java.util.SortedSet;

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.tools.RunnableProcess;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Resolves the tip taxa that could not be placed from the tree's own data, by fetching their lineages
 * from the taxonomy database (off the EDT) into the shared {@link TaxonomicLineageService} cache, then
 * runs a {@link Completion} on the EDT. Used by both "Colorize Subtrees via Taxonomic Rank" and
 * "Annotate Clades by Rank": the EDT completion only reads the now-populated cache, so it does the
 * colorize / clade-band build without ever blocking the event thread on the network.
 *
 * <p>Lives in {@code org.forester.archaeopteryx} (not {@code .tools}) so a {@link Completion} can drive
 * package-private {@link TreePanel} methods; it extends the public {@link RunnableProcess} for the
 * process pool, so a long fetch is cancellable from the Process menu.
 */
public final class OnlineTaxonResolver extends RunnableProcess {

    /** Runs on the EDT after the off-EDT fetch; {@code error} is non-null when the network failed. */
    interface Completion {
        void finished( String error );
    }

    private final MainFrame               _mf;
    private final String                  _label;
    private final SortedSet<String>       _names;
    private final TaxonomicLineageService _service;
    private final Completion              _completion;

    OnlineTaxonResolver( final MainFrame mf,
                         final String label,
                         final SortedSet<String> names,
                         final Completion completion ) {
        _mf = mf;
        _label = label;
        _names = names;
        _completion = completion;
        _service = TreePanelUtil.getDefaultLineageService();
    }

    @Override
    public void run() {
        start( _mf, _label );
        String error = null;
        try {
            for( final String name : _names ) {
                if ( isCancelled() ) {
                    break; // user cancelled from the process menu -- finish with what was resolved so far
                }
                try {
                    _service.fetch( name );
                }
                catch ( final UnknownHostException e ) {
                    error = "could not connect to the taxonomy database (no internet connection?)";
                    break; // a connection failure will hit every remaining taxon -- stop now
                }
                catch ( final IOException e ) {
                    error = e.getLocalizedMessage();
                    // a single failed taxon should not abort the rest -- keep going
                }
            }
        }
        finally {
            end( _mf );
        }
        final String err = error;
        SwingUtilities.invokeLater( () -> _completion.finished( err ) );
    }
}
