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

import java.net.UnknownHostException;

import javax.swing.JOptionPane;

import org.forester.analysis.AncestralTaxonomyInference;
import org.forester.analysis.AncestralTaxonomyInferenceException;
import org.forester.archaeopteryx.MainFrame;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;
import org.forester.ws.seqdb.SequenceDbWsTools;

public class AncestralTaxonomyInferrer extends RunnableProcess {

    private final Phylogeny            _phy;
    private final MainFrame _mf;
    private final TreePanel            _treepanel;

    public AncestralTaxonomyInferrer( final MainFrame mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    public static String getBaseUrl() {
        return SequenceDbWsTools.BASE_UNIPROT_REST_URL;
    }

    private void inferTaxonomies() {
        start( _mf, "ancestral taxonomy" );
        try {
            AncestralTaxonomyInference.inferTaxonomyFromDescendents( _phy );
        }
        catch ( final AncestralTaxonomyInferenceException e ) {
            end( _mf );
            JOptionPane.showMessageDialog( _mf,
                                           e.getMessage(),
                                           "Error during ancestral taxonomy inference",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final UnknownHostException e ) {
            end( _mf );
            JOptionPane.showMessageDialog( _mf,
                                           "Could not connect to \"" + getBaseUrl() + "\"",
                                           "Network error during ancestral taxonomy inference",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final Exception e ) {
            end( _mf );
            e.printStackTrace();
            JOptionPane.showMessageDialog( _mf,
                                           e.toString(),
                                           "Unexpected exception during ancestral taxonomy inference",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final Error e ) {
            end( _mf );
            JOptionPane.showMessageDialog( _mf,
                                           e.toString(),
                                           "Unexpected error during ancestral taxonomy inference",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        _phy.setRerootable( false );
        _treepanel.setTree( _phy );
        _mf.showWhole();
        _treepanel.setEdited( true );
        end( _mf );
        try {
            JOptionPane.showMessageDialog( _mf,
                                           "Ancestral taxonomy inference successfully completed",
                                           "Ancestral Taxonomy Inference Completed",
                                           JOptionPane.INFORMATION_MESSAGE );
        }
        catch ( final Exception e ) {
            // Not important if this fails, do nothing.
        }
    }

    @Override
    public void run() {
        inferTaxonomies();
    }
}
