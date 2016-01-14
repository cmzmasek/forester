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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx.tools;

import java.io.IOException;
import java.net.UnknownHostException;
import java.util.SortedSet;

import javax.swing.JOptionPane;

import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;
import org.forester.ws.seqdb.SequenceDbWsTools;

public final class SequenceDataRetriver extends RunnableProcess {

    private final Phylogeny            _phy;
    private final MainFrameApplication _mf;
    private final TreePanel            _treepanel;
    public final static boolean        DEBUG = false;

    public SequenceDataRetriver( final MainFrameApplication mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    @Override
    public void run() {
        execute();
    }

    private void execute() {
        start( _mf, "sequence data" );
        SortedSet<String> not_found = null;
        try {
            not_found = SequenceDbWsTools.obtainSeqInformation( _phy,
                                                                false,
                                                                true,
                                                                SequenceDbWsTools.DEFAULT_LINES_TO_RETURN );
        }
        catch ( final UnknownHostException e ) {
            JOptionPane.showMessageDialog( _mf,
                                           e.getLocalizedMessage(),
                                           "Network error during sequence data gathering",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            JOptionPane.showMessageDialog( _mf,
                                           e.toString(),
                                           "Failed to obtain sequence data",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        finally {
            end( _mf );
        }
        _treepanel.setTree( _phy );
        _mf.showWhole();
        _treepanel.setEdited( true );
        if ( ( not_found != null ) && ( not_found.size() > 0 ) ) {
            int max = not_found.size();
            boolean more = false;
            if ( max > 20 ) {
                more = true;
                max = 20;
            }
            final StringBuffer sb = new StringBuffer();
            if ( not_found.size() == 1 ) {
                sb.append( "For the following node no data was found:\n" );
            }
            else {
                sb.append( "For the following nodes no data was found (total: " + not_found.size() + "):\n" );
            }
            int i = 0;
            for( final String string : not_found ) {
                if ( i > 19 ) {
                    break;
                }
                sb.append( string );
                sb.append( "\n" );
                ++i;
            }
            if ( more ) {
                sb.append( "..." );
            }
            try {
                JOptionPane.showMessageDialog( _mf,
                                               sb.toString(),
                                               "Sequence Tool Completed",
                                               JOptionPane.WARNING_MESSAGE );
            }
            catch ( final Exception e ) {
                // Not important if this fails, do nothing.
            }
        }
        else {
            try {
                JOptionPane.showMessageDialog( _mf,
                                               "Sequence tool successfully completed",
                                               "Sequence Tool Completed",
                                               JOptionPane.INFORMATION_MESSAGE );
            }
            catch ( final Exception e ) {
                // Not important if this fails, do nothing.
            }
        }
    }
}
