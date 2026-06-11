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

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

final class NodeFrame extends javax.swing.JFrame {

    private static final long serialVersionUID = -6943510233968557246L;
    private final TreePanel   _reepanel;
    private int               _index           = -1;

    NodeFrame( final PhylogenyNode n, final Phylogeny tree, final TreePanel tp, final int x ) {
        super( "Node " + ( ForesterUtil.isEmpty( n.getName() ) ? n.getId() : n.getName() ) );
        _reepanel = tp;
        setSize( AptxConstants.NODE_FRAME_SIZE );
        _index = x;
        final Container contentPane = getContentPane();
        final NodePanel nodepanel = new NodePanel( n, this );
        contentPane.add( nodepanel );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                close();
            }
        } );
        setResizable( true );
        nodepanel.setVisible( true );
        setVisible( true );
    }

    
    void close() {
        remove(); // to release slot in array
        dispose();
    }
    
    NodeFrame( final PhylogenyNode n, final Phylogeny tree, final TreePanel tp, final int x, final String dummy ) {
        super( "Editable Node " + ( ForesterUtil.isEmpty( n.getName() ) ? n.getId() : n.getName() ) );
        _reepanel = tp;
        setSize( AptxConstants.NODE_FRAME_SIZE );
        _index = x;
        final Container contentPane = getContentPane();
        final NodeEditPanel nodepanel = new NodeEditPanel( n, tp, this );
        contentPane.add( nodepanel, BorderLayout.CENTER );
        addWindowListener( new WindowAdapter() {

            @Override
            public void windowClosing( final WindowEvent e ) {
                try {
                    nodepanel.writeAll();
                }
                catch ( final Exception ex ) {
                    // Do nothing.
                }
                remove(); // to release slot in array
                dispose();
            }
        } );
        setResizable( false );
        nodepanel.setVisible( true );
        setVisible( true );
    }

    TreePanel getTreePanel() {
        return _reepanel;
    }

    void remove() {
        if ( _index > -1 ) {
            _reepanel.removeEditNodeFrame( _index ); // to release slot in array
        }
    }
}
