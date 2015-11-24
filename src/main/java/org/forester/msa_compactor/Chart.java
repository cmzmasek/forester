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
import java.text.DecimalFormat;
import java.text.NumberFormat;
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

    final private static NumberFormat NF_1             = new DecimalFormat( "0.##" );
    private static final long         serialVersionUID = -5292420246132943515L;
    private ChartPanel                _chart_panel     = null;
    private final int                 _initial_number_of_seqs;
    private final JMenuItem           _m_exit          = new JMenuItem();
    private final List<MsaProperties> _msa_props;
    private final boolean             _show_msa_qual;
    private final String              _title;

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
            final double[][] seqs_length = new double[ _msa_props.size() ][ 2 ];
            int max_length = -1;
            int min_length = Integer.MAX_VALUE;
            double max_gap_ratio = -1;
            double min_gap_ratio = Double.MAX_VALUE;
            double max_avg_gap_count = -1;
            double min_avg_gap_count = Double.MAX_VALUE;
            for( int i = 0; i < _msa_props.size(); ++i ) {
                seqs_length[ i ][ 0 ] = _initial_number_of_seqs - _msa_props.get( i ).getNumberOfSequences();
                //
                final int length = _msa_props.get( i ).getLength();
                seqs_length[ i ][ 1 ] = length;
                if ( length > max_length ) {
                    max_length = length;
                }
                if ( length < min_length ) {
                    min_length = length;
                }
                //
                final double gap_ratio = _msa_props.get( i ).getGapRatio();
                if ( gap_ratio > max_gap_ratio ) {
                    max_gap_ratio = gap_ratio;
                }
                if ( gap_ratio < min_gap_ratio ) {
                    min_gap_ratio = gap_ratio;
                }
                //
                final double avg_gap_count = _msa_props.get( i ).getAvgNumberOfGaps();
                if ( avg_gap_count > max_avg_gap_count ) {
                    max_avg_gap_count = avg_gap_count;
                }
                if ( avg_gap_count < min_avg_gap_count ) {
                    min_avg_gap_count = avg_gap_count;
                }
            }
            model.addData( seqs_length, "Length" + " (" + minMaxToString( min_length, max_length ) + ")" );
            model.setSeriesLine( "Series " + "Length", true );
            model.setSeriesMarker( "Series " + "Length", false );
            final double[][] seqs_gaps = new double[ _msa_props.size() ][ 2 ];
            double max_ent7 = -1;
            double max_ent21 = -1;
            double min_ent7 = Double.MAX_VALUE;
            double min_ent21 = Double.MAX_VALUE;
            if ( _show_msa_qual ) {
                for( int i = 0; i < _msa_props.size(); ++i ) {
                    final double ent7 = _msa_props.get( i ).getEntropy7();
                    if ( ent7 > max_ent7 ) {
                        max_ent7 = ent7;
                    }
                    if ( ent7 < max_ent7 ) {
                        min_ent7 = ent7;
                    }
                    final double ent21 = _msa_props.get( i ).getEntropy21();
                    if ( ent21 > min_ent21 ) {
                        max_ent21 = ent21;
                    }
                    if ( ent21 < min_ent21 ) {
                        min_ent21 = ent21;
                    }
                }
            }
            final double gap_ratio_factor = ( max_length / 2.0 ) / max_gap_ratio;
            final double avg_gaps_counts_factor = ( max_length / 2.0 ) / max_avg_gap_count;
            final double ent7_factor = ( max_length / 2.0 ) / max_ent7;
            final double ent21_factor = ( max_length / 2.0 ) / max_ent21;
            for( int i = 0; i < _msa_props.size(); ++i ) {
                seqs_gaps[ i ][ 0 ] = _initial_number_of_seqs - _msa_props.get( i ).getNumberOfSequences();
                seqs_gaps[ i ][ 1 ] = ForesterUtil.roundToInt( _msa_props.get( i ).getGapRatio() * gap_ratio_factor );
            }
            model.addData( seqs_gaps, "Gap Ratio" + " (" + minMaxToString( min_gap_ratio, max_gap_ratio ) + ")" );
            model.setSeriesLine( "Series " + "Gap Ratio", true );
            model.setSeriesMarker( "Series " + "Gap Ratio", false );
            final double[][] gap_counts = new double[ _msa_props.size() ][ 2 ];
            for( int i = 0; i < _msa_props.size(); ++i ) {
                gap_counts[ i ][ 0 ] = _initial_number_of_seqs - _msa_props.get( i ).getNumberOfSequences();
                gap_counts[ i ][ 1 ] = ForesterUtil.roundToInt( _msa_props.get( i ).getAvgNumberOfGaps()
                                                                * avg_gaps_counts_factor );
            }
            model.addData( gap_counts, "Mean Gap Count" + " (" + minMaxToString( min_avg_gap_count, max_avg_gap_count )
                    + ")" );
            model.setSeriesLine( "Series " + "Mean Gap Count", true );
            model.setSeriesMarker( "Series " + "Mean Gap Count", false );
            if ( _show_msa_qual ) {
                final double[][] entropy7 = new double[ _msa_props.size() ][ 2 ];
                for( int i = 0; i < _msa_props.size(); ++i ) {
                    entropy7[ i ][ 0 ] = _initial_number_of_seqs - _msa_props.get( i ).getNumberOfSequences();
                    entropy7[ i ][ 1 ] = ForesterUtil.roundToInt( _msa_props.get( i ).getEntropy7() * ent7_factor );
                }
                model.addData( entropy7, "Entropy norm 7" + " (" + minMaxToString( min_ent7, max_ent7 ) + ")" );
                model.setSeriesLine( "Series " + "Entropy norm 7", true );
                model.setSeriesMarker( "Series " + "Entropy norm 7", false );
                //
                final double[][] entropy21 = new double[ _msa_props.size() ][ 2 ];
                for( int i = 0; i < _msa_props.size(); ++i ) {
                    entropy21[ i ][ 0 ] = _initial_number_of_seqs - _msa_props.get( i ).getNumberOfSequences();
                    entropy21[ i ][ 1 ] = ForesterUtil.roundToInt( _msa_props.get( i ).getEntropy21() * ent21_factor );
                }
                model.addData( entropy21, "Entropy norm 21" + " (" + minMaxToString( min_ent21, max_ent21 ) + ")" );
                model.setSeriesLine( "Series " + "Entropy norm 21", true );
                model.setSeriesMarker( "Series " + "Entropy norm 21", false );
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

    private final static String minMaxToString( final double min, final double max ) {
        return NF_1.format( min ) + "-" + NF_1.format( max );
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
