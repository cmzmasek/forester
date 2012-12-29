// $Id:
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

package org.forester.development;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;

import javax.swing.JComponent;

import org.forester.msa.Msa;

public class MsaRenderer extends JComponent {

    private static final long serialVersionUID = -68078011081748093L;

    public static boolean isMouseEventAltered( final MouseEvent event ) {
        return event.isShiftDown() || event.isAltDown() || event.isControlDown() || event.isAltGraphDown()
                || event.isMetaDown();
    }
    //private PlateDisplayPanel _heatMapPanel;
    private final int        _rows;
    private final int        _columns;
    private int              _wellSize;
    private AbstractRenderer _wells[][];
    private double           _min;
    private double           _max;
    private double           _mean;
    private Color            _minColor;
    private Color            _maxColor;
    private Color            _meanColor;
    private boolean          _useMean;
    // private Rubberband        _rubberband;
    private final Msa        _msa;
    private final JComponent _parent;

    public MsaRenderer( final Msa msa, final int unit_size, final JComponent parent ) {
        _parent = parent;
        _msa = msa;
        _rows = _msa.getNumberOfSequences();
        _columns = _msa.getLength();
        setWellSize( unit_size );
        addMouseListeners();
        initializeWells();
        //setRubberband( new RubberbandRectangle( this ) );
    }

    private void addMouseListeners() {
        addMouseMotionListener( new MouseMotionAdapter() {

            @Override
            public void mouseDragged( final MouseEvent event ) {
                //  if ( ( ( event.getModifiers() & 0x10 ) != 0 ) && getRubberband().isActive()
                //          && !PlateRenderer.isMouseEventAltered( event ) ) {
                //     getRubberband().stretch( event.getPoint() );
                // }
            }
        } );
        addMouseListener( new MouseAdapter() {

            @Override
            public void mouseClicked( final MouseEvent event ) {
                //  if ( ( ( event.getModifiers() & 4 ) != 0 ) || PlateRenderer.isMouseEventAltered( event ) ) {
                //    getInfo( new Rectangle( event.getX(), event.getY(), 1, 1 ) );
                //   }
                //   else {
                //   changeSelected( new Rectangle( event.getX(), event.getY(), 1, 1 ), true );
                //      event.consume();
                //  }
            }

            @Override
            public void mousePressed( final MouseEvent event ) {
                // if ( ( ( event.getModifiers() & 4 ) != 0 ) || PlateRenderer.isMouseEventAltered( event ) ) {
                //    getInfo( new Rectangle( event.getX(), event.getY(), 1, 1 ) );
                // }
                // if ( ( ( event.getModifiers() & 0x10 ) != 0 ) && getRubberband().isActive() ) {
                //     getRubberband().anchor( event.getPoint() );
                // }
            }

            @Override
            public void mouseReleased( final MouseEvent event ) {
                // if ( ( ( event.getModifiers() & 0x10 ) != 0 ) && getRubberband().isActive() ) {
                //     getRubberband().end( event.getPoint() );
                //     rubberbandEnded( getRubberband() );
                // }
            }
        } );
    }

    public AbstractRenderer getAbstractRenderer( final int row, final int col ) {
        return _wells[ row ][ col ];
    }

    public int getColumns() {
        return _columns;
    }

    private double getMax() {
        return _max;
    }

    private Color getMaxColor() {
        return _maxColor;
    }

    private double getMean() {
        return _mean;
    }

    private Color getMeanColor() {
        return _meanColor;
    }

    private double getMin() {
        return _min;
    }

    private Color getMinColor() {
        return _minColor;
    }

    @Override
    public Dimension getPreferredSize() {
        final int width = ( getWellSize() + 1 ) * ( getColumns() + 1 );
        final int hight = ( ( getWellSize() + 1 ) * ( getRows() + 1 ) ) + 30;
        return new Dimension( width, hight );
    }

    public int getRows() {
        return _rows;
    }

    //  private Rubberband getRubberband() {
    //      return _rubberband;
    //  }
    private int getWellSize() {
        return _wellSize;
    }

    private void initializeWells() {
        _wells = new AbstractRenderer[ getRows() + 1 ][ getColumns() + 1 ];
        for( int row = 0; row < getRows(); row++ ) {
            for( int col = 0; col < ( getColumns() + 1 ); col++ ) {
                AbstractRenderer r;
                if ( col == getColumns() ) {
                    //  r = new LabelRenderer( PlateData.ALPHABET[ row % PlateData.ALPHABET.length ] + "", this );
                }
                //else if ( getPlateData().getData( row, col ) == null ) {
                //    r = new WellRenderer( new WellData(), this );
                // }
                else {
                    r = new ResidueRenderer( getMsa().getResidueAt( row, col ), this );
                }
                //  r.setVisible( true );
                //  setAbstractRenderer( r, row, col );
            }
        }
        for( int col = 0; col < ( getColumns() + 1 ); col++ ) {
            //  AbstractRenderer r;
            if ( col == getColumns() ) {
                //      r = new LabelRenderer( "", this );
            }
            else {
                //       r = new LabelRenderer( ( col + 1 ) + "", this );
            }
            //  r.setVisible( true );
            // setAbstractRenderer( r, getRows(), col );
        }
    }

    private Msa getMsa() {
        return _msa;
    }

    public void inverseMarkedOfWell( final int well_row, final int well_col ) {
        final ResidueRenderer rend = ( ResidueRenderer ) getAbstractRenderer( well_row, well_col );
        if ( rend.isMarked() ) {
            rend.setIsMarked( false );
        }
        else {
            rend.setIsMarked( true );
        }
    }

    private boolean isUseMean() {
        return _useMean;
    }

    @Override
    public void paint( final Graphics g ) {
        g.setColor( Color.white );
        //   g.setFont( getPlateDisplayPanel().getPlateTitleFont() );
        // g
        //         .drawString( "Number:" + getPlateData().getName() + " Replicate:"
        //                 + ( getPlateData().getReplicateNumber() + 1 ), 10, 20 );
        for( int row = 0; row < ( getRows() + 1 ); row++ ) {
            for( int col = 0; col < ( getColumns() + 1 ); col++ ) {
                getAbstractRenderer( row, col ).paint( g );
            }
        }
    }

    public void resetWellColors() {
        for( int row = 0; row < getRows(); row++ ) {
            for( int col = 0; col < getColumns(); col++ ) {
                final ResidueRenderer r = ( ResidueRenderer ) getAbstractRenderer( row, col );
                if ( isUseMean() ) {
                    r.resetWellColor( getMin(), getMax(), getMean(), getMinColor(), getMaxColor(), getMeanColor() );
                }
                else {
                    r.resetWellColor( getMin(), getMax(), getMinColor(), getMaxColor() );
                }
            }
        }
    }

    public void resetWellSize( final int well_size ) {
        setWellSize( well_size );
        final int factor = well_size + 0;
        for( int row = 0; row < ( getRows() + 1 ); row++ ) {
            for( int col = 0; col < ( getColumns() + 1 ); col++ ) {
                final AbstractRenderer r = getAbstractRenderer( row, col );
                r.setX( 10 + ( factor * col ) );
                r.setY( ( factor * row ) + 30 );
                r.setWellSize( well_size );
            }
        }
    }

    //  private void rubberbandEnded( final Rubberband rb ) {
    // changeSelected( rb.getBounds(), false );
    //     repaint();
    // }
    private void setAbstractRenderer( final AbstractRenderer ar, final int row, final int col ) {
        _wells[ row ][ col ] = ar;
    }

    public void setFlaggedStatusOfOutlierWells( final boolean set_flagged_to_this ) {
        for( int row = 0; row < getRows(); row++ ) {
            for( int col = 0; col < getColumns(); col++ ) {
                final ResidueRenderer wr = ( ResidueRenderer ) getAbstractRenderer( row, col );
                // if ( wr.isFlaggedByStatisticalAnalysis() ) {
                //     wr.setIsUserFlagged( set_flagged_to_this );
                //  }
            }
        }
    }

    public void setFlaggedStatusOfSelectedWells( final boolean set_flagged_to_this ) {
        for( int row = 0; row < getRows(); row++ ) {
            for( int col = 0; col < getColumns(); col++ ) {
                final ResidueRenderer wr = ( ResidueRenderer ) getAbstractRenderer( row, col );
                if ( wr.isSelected() ) {
                    //      wr.setIsUserFlagged( set_flagged_to_this );
                    wr.setIsSelected( false );
                }
            }
        }
    }

    public void setIsFlaggingStatusChangedToFalse() {
        for( int row = 0; row < getRows(); row++ ) {
            for( int col = 0; col < getColumns(); col++ ) {
                final ResidueRenderer wr = ( ResidueRenderer ) getAbstractRenderer( row, col );
                //      wr.setIsFlaggingStatusChanged( false );
            }
        }
    }

    private void setIsSelectedOfAll( final boolean isSelected ) {
        for( int col = 0; col < ( getColumns() + 1 ); col++ ) {
            setIsSelectedOfColumn( col, isSelected );
        }
    }

    private void setIsSelectedOfColumn( final int column, final boolean isSelected ) {
        for( int row = 0; row < ( getRows() + 1 ); row++ ) {
            getAbstractRenderer( row, column ).setIsSelected( isSelected );
        }
    }

    private void setIsSelectedOfRow( final int row, final boolean isSelected ) {
        for( int col = 0; col < ( getColumns() + 1 ); col++ ) {
            getAbstractRenderer( row, col ).setIsSelected( isSelected );
        }
    }

    private void setIsSelectedOfRowAlternating( final int row, final boolean even ) {
        boolean selected = even;
        for( int col = 0; col < getColumns(); col++ ) {
            getAbstractRenderer( row, col ).setIsSelected( selected );
            selected = !selected;
        }
    }

    private void setIsSelectedToQuarter( final int quarter ) {
        boolean evenRow = false;
        boolean evenColumn = false;
        if ( quarter <= 1 ) {
            evenRow = true;
            evenColumn = true;
        }
        else if ( quarter == 2 ) {
            evenColumn = true;
        }
        else if ( quarter == 3 ) {
            evenRow = true;
        }
        for( int row = 0; row < getRows(); row++ ) {
            if ( evenColumn ) {
                setIsSelectedOfRowAlternating( row, evenRow );
            }
            else {
                setIsSelectedOfRow( row, false );
            }
            evenColumn = !evenColumn;
        }
    }

    public void setMarkedOfAllWellsToFalse() {
        for( int row = 0; row < getRows(); row++ ) {
            for( int col = 0; col < getColumns(); col++ ) {
                final ResidueRenderer wr = ( ResidueRenderer ) getAbstractRenderer( row, col );
                //  rend.setIsMarked( false );
            }
        }
    }

    public void setMax( final double max ) {
        _max = max;
    }

    void setMaxColor( final Color maxColor ) {
        _maxColor = maxColor;
    }

    void setMean( final double mean ) {
        _mean = mean;
    }

    public void setMeanColor( final Color meanColor ) {
        _meanColor = meanColor;
    }

    public void setMin( final double min ) {
        _min = min;
    }

    void setMinColor( final Color minColor ) {
        _minColor = minColor;
    }

    //    private void setRubberband( final Rubberband rb ) {
    //        if ( _rubberband != null ) {
    //            _rubberband.setActive( false );
    //        }
    //        _rubberband = rb;
    //        if ( _rubberband != null ) {
    //            _rubberband.setComponent( this );
    //            _rubberband.setActive( true );
    //        }
    //    }
    public void setUseMean( final boolean useMean ) {
        _useMean = useMean;
    }

    private void setWellSize( final int wellSize ) {
        _wellSize = wellSize;
    }

    public void unSelectUnMarkAll() {
        for( int row = 0; row < getRows(); row++ ) {
            for( int col = 0; col < getColumns(); col++ ) {
                final ResidueRenderer wr = ( ResidueRenderer ) getAbstractRenderer( row, col );
                wr.setIsSelected( false );
                wr.setIsMarked( false );
            }
        }
    }
}
