// Exp $
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx.tools;

import java.net.UnknownHostException;

import javax.swing.JOptionPane;

import org.forester.analysis.AncestralTaxonomyInference;
import org.forester.analysis.AncestralTaxonomyInferenceException;
import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;
import org.forester.ws.uniprot.UniProtWsTools;

public class AncestralTaxonomyInferrer implements Runnable {

    private final Phylogeny            _phy;
    private final MainFrameApplication _mf;
    private final TreePanel            _treepanel;

    public AncestralTaxonomyInferrer( final MainFrameApplication mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    private String getBaseUrl() {
        return UniProtWsTools.BASE_URL;
    }

    private void inferTaxonomies() {
        _mf.getMainPanel().getCurrentTreePanel().setWaitCursor();
        try {
            AncestralTaxonomyInference.inferTaxonomyFromDescendents( _phy );
        }
        catch ( final AncestralTaxonomyInferenceException e ) {
            _mf.getMainPanel().getCurrentTreePanel().setArrowCursor();
            JOptionPane.showMessageDialog( _mf,
                                           e.getMessage(),
                                           "Error during ancestral taxonomy inference",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final UnknownHostException e ) {
            _mf.getMainPanel().getCurrentTreePanel().setArrowCursor();
            JOptionPane.showMessageDialog( _mf,
                                           "Could not connect to \"" + getBaseUrl() + "\"",
                                           "Network error during ancestral taxonomy inference",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final Exception e ) {
            _mf.getMainPanel().getCurrentTreePanel().setArrowCursor();
            e.printStackTrace();
            JOptionPane.showMessageDialog( _mf,
                                           e.toString(),
                                           "Unexpected error during ancestral taxonomy inference",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        _mf.getMainPanel().getCurrentTreePanel().setArrowCursor();
        _phy.setRerootable( false );
        _treepanel.setTree( _phy );
        _mf.showWhole();
        _treepanel.setEdited( true );
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
