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

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import org.forester.analysis.AncestralTaxonomyInference;
import org.forester.analysis.AncestralTaxonomyInferenceException;
import org.forester.archaeopteryx.MainFrame;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;

/**
 * The Analysis-menu "Infer Ancestor Taxonomies" operation: assigns each internal node the deepest
 * common lineage of its descendants. This is a <b>local</b> computation (no network -- it intersects the
 * lineages already present on the descendants), run on a background thread; the result is committed to
 * the tree and reported on the EDT (Swing is not thread-safe).
 */
public class AncestralTaxonomyInferrer extends RunnableProcess {

    private final Phylogeny _phy;
    private final MainFrame _mf;
    private final TreePanel _treepanel;

    public AncestralTaxonomyInferrer( final MainFrame mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    private void inferTaxonomies() {
        start( _mf, "ancestral taxonomy" );
        String error = null;
        try {
            AncestralTaxonomyInference.inferTaxonomyFromDescendents( _phy );
        }
        catch ( final AncestralTaxonomyInferenceException e ) {
            error = e.getMessage();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            error = e.toString();
        }
        catch ( final Error e ) {
            error = e.toString();
        }
        finally {
            end( _mf );
        }
        final String err = error;
        // commit + report on the EDT
        SwingUtilities.invokeLater( () -> {
            if ( err != null ) {
                JOptionPane.showMessageDialog( _mf, err, "Error during ancestral taxonomy inference",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            _phy.setRerootable( false );
            _treepanel.setTree( _phy );
            _mf.showWhole();
            _treepanel.setEdited( true );
            try {
                JOptionPane.showMessageDialog( _mf, "Ancestral taxonomy inference successfully completed",
                                               "Ancestral Taxonomy Inference Completed", JOptionPane.INFORMATION_MESSAGE );
            }
            catch ( final Exception e ) {
                // not important if the dialog fails
            }
        } );
    }

    @Override
    public void run() {
        inferTaxonomies();
    }
}
