// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2014 Christian M. Zmasek
// Copyright (C) 2014 Sanford-Burnham Medical Research Institute
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.msa_compactor;

import java.awt.BorderLayout;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.UIManager;
import javax.swing.WindowConstants;

import org.forester.util.ForesterUtil;

import com.approximatrix.charting.coordsystem.BoxCoordSystem;
import com.approximatrix.charting.model.MultiScatterDataModel;
import com.approximatrix.charting.render.MultiScatterChartRenderer;
import com.approximatrix.charting.swing.ChartPanel;

public final class Chart extends JDialog implements ActionListener {

    private static final long   serialVersionUID = -5292420246132943515L;
    private ChartPanel          _chart_panel     = null;
    private final JMenuItem     _m_exit          = new JMenuItem();
    private List<MsaProperties> _msa_props;
    private final boolean       _show_msa_qual;
    private final int           _initial_number_of_seqs;
    private final String        _title;

    private Chart( final List<MsaProperties> msa_props,
                   final int initial_number_of_seqs,
                   final boolean show_msa_qual,
                   final String title ) {
        super();
        _msa_props = msa_props;
        _title = title;
        _initial_number_of_seqs = initial_number_of_seqs;
        _show_msa_qual = show_msa_qual;
        setTitle( "msa compactor" );
        setSize( 600, 500 );
        setResizable( true );
        final JPanel content_pane = new JPanel();
        content_pane.setLayout( new BorderLayout() );
        setContentPane( content_pane );
        final JMenuBar menu_bar = new JMenuBar();
        final JMenu file_menu = new JMenu();
        file_menu.setText( "File" );
        _m_exit.setText( "Exit" );
        file_menu.add( _m_exit );
        menu_bar.add( file_menu );
        setJMenuBar( menu_bar );
        setDefaultCloseOperation( WindowConstants.DISPOSE_ON_CLOSE );
        _m_exit.addActionListener( this );
        content_pane.add( obtainChartPanel(), BorderLayout.CENTER );
    }

    @Override
    public void actionPerformed( final java.awt.event.ActionEvent e ) {
        if ( e.getSource() == _m_exit ) {
            dispose();
        }
    }

    private ChartPanel obtainChartPanel() {
        if ( _chart_panel == null ) {
            final MultiScatterDataModel model = new MultiScatterDataModel();
            if ( ( _msa_props == null ) || _msa_props.isEmpty() ) {
                _msa_props = new ArrayList<MsaProperties>();
                final MsaProperties p0 = new MsaProperties( 10, 200, 0.5, 0.1 );
                final MsaProperties p1 = new MsaProperties( 9, 190, 0.49, 0.2 );
                final MsaProperties p2 = new MsaProperties( 8, 150, 0.2, 0.3 );
                final MsaProperties p3 = new MsaProperties( 7, 145, 0.2, 0.4 );
                _msa_props.add( p0 );
                _msa_props.add( p1 );
                _msa_props.add( p2 );
                _msa_props.add( p3 );
            }
            final double[][] seqs_length = new double[ _msa_props.size() ][ 2 ];
            for( int i = 0; i < _msa_props.size(); ++i ) {
                seqs_length[ i ][ 0 ] = _initial_number_of_seqs - _msa_props.get( i ).getNumberOfSequences();
                seqs_length[ i ][ 1 ] = _msa_props.get( i ).getLength();
            }
            model.addData( seqs_length, "Length" );
            model.setSeriesLine( "Series " + "Length", true );
            model.setSeriesMarker( "Series " + "Length", false );
            final double[][] seqs_gaps = new double[ _msa_props.size() ][ 2 ];
            for( int i = 0; i < _msa_props.size(); ++i ) {
                seqs_gaps[ i ][ 0 ] = _initial_number_of_seqs - _msa_props.get( i ).getNumberOfSequences();
                seqs_gaps[ i ][ 1 ] = ForesterUtil.roundToInt( _msa_props.get( i ).getGapRatio() * 200 );
            }
            model.addData( seqs_gaps, "Gap ratio" );
            model.setSeriesLine( "Series " + "Gap ratio", true );
            model.setSeriesMarker( "Series " + "Gap ratio", false );
            if ( _show_msa_qual ) {
                final double[][] seqs_identity = new double[ _msa_props.size() ][ 2 ];
                for( int i = 0; i < _msa_props.size(); ++i ) {
                    seqs_identity[ i ][ 0 ] = _initial_number_of_seqs - _msa_props.get( i ).getNumberOfSequences();
                    seqs_identity[ i ][ 1 ] = ForesterUtil
                            .roundToInt( _msa_props.get( i ).getAverageIdentityRatio() * 200 );
                }
                model.addData( seqs_identity, "mean MSA column identity" );
                model.setSeriesLine( "Series " + "mean MSA column identity", true );
                model.setSeriesMarker( "Series " + "mean MSA column identity", false );
            }
            final BoxCoordSystem coord = new BoxCoordSystem( model );
            coord.setUnitFont( coord.getUnitFont().deriveFont( 16.0f ) );
            coord.setXAxisUnit( "Number of Removed Sequences" );
            coord.setPaintGrid( true );
            coord.setYAxisUnit( "MSA Length" );
            _chart_panel = new ChartPanel( model, _title );
            _chart_panel.setCoordSystem( coord );
            final MultiScatterChartRenderer renderer = new MultiScatterChartRenderer( coord, model );
            renderer.setAllowBuffer( false );
            _chart_panel.addChartRenderer( renderer, 0 );
        }
        return _chart_panel;
    }

    public static void display( final List<MsaProperties> msa_props,
                                final int initial_number_of_seqs,
                                final boolean show_msa_qual,
                                final String title ) {
        try {
            UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
        }
        final Chart chart = new Chart( msa_props, initial_number_of_seqs, show_msa_qual, title );
        chart.setVisible( true );
    }

    public static void main( final String[] args ) {
        try {
            UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName() );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
        }
        final Chart temp = new Chart( null, 0, true, "title" );
        temp.setVisible( true );
    }
}
