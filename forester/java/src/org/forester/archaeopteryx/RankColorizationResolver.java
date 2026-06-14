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

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.tools.RunnableProcess;
import org.forester.ws.seqdb.TaxonomicLineageService;

/**
 * Background half of "Colorize Subtrees via Taxonomic Rank": resolves the tip taxa that could not be
 * placed at the chosen rank from the tree's own data by fetching their lineages from the taxonomy
 * database (off the EDT), then recolors and reports the result on the EDT. The fetched lineages land
 * in the shared {@link TaxonomicLineageService} cache, so the EDT recolor (which only reads the
 * cache) now places those tips -- never blocking the event thread on the network.
 *
 * <p>Lives in the {@code org.forester.archaeopteryx} package (not {@code .tools}) so it can drive the
 * package-private {@link TreePanel#colorByRank}/{@link TreePanel#reportRankColorization} without
 * widening their visibility; it still extends the public {@link RunnableProcess} for the process pool.
 */
public final class RankColorizationResolver extends RunnableProcess {

    private final MainFrame                _mf;
    private final TreePanel                _tp;
    private final String                   _rank;
    private final SortedSet<String>        _names;
    private final TaxonomicLineageService  _service;

    RankColorizationResolver( final MainFrame mf,
                              final TreePanel tp,
                              final String rank,
                              final SortedSet<String> names ) {
        _mf = mf;
        _tp = tp;
        _rank = rank;
        _names = names;
        _service = TreePanelUtil.getDefaultLineageService();
    }

    @Override
    public void run() {
        start( _mf, "rank colorization (" + _rank + ")" );
        String error = null;
        try {
            for( final String name : _names ) {
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
        SwingUtilities.invokeLater( () -> {
            final int colorized = _tp.colorByRank( _rank );
            if ( err != null ) {
                // colorByRank mutated branch colors but does not set the edited flag (only
                // reportRankColorization does); on the error path we skip that dialog, so mark the
                // tree edited here -- otherwise the colorization is lost with no "save changes?" prompt.
                if ( colorized > 0 ) {
                    _tp.setEdited( true );
                }
                JOptionPane.showMessageDialog( _mf,
                                               "Colorized " + colorized + " clade(s), but some taxa could not be resolved:\n"
                                                       + err,
                                               "Taxonomy Rank-Colorization (" + _rank + ")",
                                               JOptionPane.WARNING_MESSAGE );
            }
            else {
                _tp.reportRankColorization( _rank, colorized );
            }
        } );
    }
}
