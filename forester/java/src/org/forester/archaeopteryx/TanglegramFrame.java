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

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JToolBar;
import javax.swing.JViewport;

import org.forester.phylogeny.Phylogeny;

/**
 * A stand-alone window hosting a {@link TanglegramPanel}. Deliberately a separate top-level frame (not a tab): the
 * main workspace's tab machinery is index-coupled to {@code TreePanel}, so a non-{@code TreePanel} tab would desync
 * it. A separate window also fits a read-only comparison view -- it carries only a minimal zoom/fit toolbar plus a
 * summary of the link (field, connectors, crossings, unmatched), and none of the tree-editing control panel.
 */
final class TanglegramFrame extends JFrame {

    private static final long     serialVersionUID = 1L;
    private final TanglegramPanel _panel;

    TanglegramFrame( final Phylogeny left, final Phylogeny right, final TanglegramLinker.LinkField field,
                     final String left_name, final String right_name ) {
        super( "Tanglegram: " + left_name + "  ↔  " + right_name );
        _panel = new TanglegramPanel( left, right, field );
        final JScrollPane scroll = new JScrollPane( _panel );
        scroll.getViewport().setScrollMode( JViewport.SIMPLE_SCROLL_MODE );
        setLayout( new BorderLayout() );
        add( buildToolbar( _panel, left_name, right_name ), BorderLayout.NORTH );
        add( scroll, BorderLayout.CENTER );
        setDefaultCloseOperation( DISPOSE_ON_CLOSE );
        setSize( 1000, 700 );
        setLocationRelativeTo( null );
    }

    TanglegramPanel getTanglegramPanel() {
        return _panel;
    }

    private static JToolBar buildToolbar( final TanglegramPanel panel, final String left_name,
                                          final String right_name ) {
        final JToolBar bar = new JToolBar();
        bar.setFloatable( false );
        final JButton zoom_in = new JButton( "Zoom In" );
        final JButton zoom_out = new JButton( "Zoom Out" );
        final JButton fit = new JButton( "Fit" );
        zoom_in.addActionListener( e -> panel.zoomIn() );
        zoom_out.addActionListener( e -> panel.zoomOut() );
        fit.addActionListener( e -> panel.fit() );
        bar.add( zoom_in );
        bar.add( zoom_out );
        bar.add( fit );
        bar.addSeparator();
        bar.add( new JLabel( summary( panel, left_name, right_name ) ) );
        bar.add( Box.createHorizontalGlue() );
        return bar;
    }

    private static String summary( final TanglegramPanel panel, final String left_name, final String right_name ) {
        return "Link: " + panel.getField().label() + "   •   " + panel.getResult().getLinks().size()
                + " connectors, " + panel.getCrossingCount() + " crossings   •   " + left_name + " ("
                + panel.getLeftTipCount() + " tips) ↔ " + right_name + " (" + panel.getRightTipCount()
                + " tips)   •   " + panel.getUnmatchedCount() + " unmatched";
    }
}
