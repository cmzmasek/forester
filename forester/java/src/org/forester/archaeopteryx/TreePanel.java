// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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

package org.forester.archaeopteryx;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.InputEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.geom.AffineTransform;
import java.awt.geom.Arc2D;
import java.awt.geom.CubicCurve2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.QuadCurve2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URLEncoder;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;

import javax.swing.BorderFactory;
import javax.swing.JApplet;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JTextArea;
import javax.swing.Popup;
import javax.swing.PopupFactory;

import org.forester.archaeopteryx.Configuration.EXT_NODE_DATA_RETURN_ON;
import org.forester.archaeopteryx.ControlPanel.NodeClickAction;
import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.archaeopteryx.phylogeny.data.RenderableVector;
import org.forester.archaeopteryx.tools.Blast;
import org.forester.archaeopteryx.tools.ImageLoader;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.io.writers.SequenceWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.NodeData.NODE_DATA;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.PropertiesMap;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.SequenceRelation;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;
import org.forester.util.TaxonomyUtil;

public final class TreePanel extends JPanel implements ActionListener, MouseWheelListener, Printable {

    final static Cursor                  ARROW_CURSOR                                       = Cursor.getPredefinedCursor( Cursor.DEFAULT_CURSOR );
    final static Cursor                  CUT_CURSOR                                         = Cursor.getPredefinedCursor( Cursor.CROSSHAIR_CURSOR );
    final static Cursor                  HAND_CURSOR                                        = Cursor.getPredefinedCursor( Cursor.HAND_CURSOR );
    final static Cursor                  MOVE_CURSOR                                        = Cursor.getPredefinedCursor( Cursor.MOVE_CURSOR );
    final static Cursor                  WAIT_CURSOR                                        = Cursor.getPredefinedCursor( Cursor.WAIT_CURSOR );
    final private static double          _180_OVER_PI                                       = 180.0 / Math.PI;
    private static final float           ANGLE_ROTATION_UNIT                                = ( float ) ( Math.PI / 32 );
    private final static int             CONFIDENCE_LEFT_MARGIN                             = 4;
    private final static int             EURO_D                                             = 10;
    private final static NumberFormat    FORMATTER_BRANCH_LENGTH;
    private final static NumberFormat    FORMATTER_CONFIDENCE;
    private static final float           HALF_PI                                            = ( float ) ( Math.PI / 2.0 );
    private final static int             LIMIT_FOR_HQ_RENDERING                             = 2000;
    private final static int             MAX_NODE_FRAMES                                    = 10;
    private final static int             MAX_SUBTREES                                       = 100;
    private final static int             MIN_ROOT_LENGTH                                    = 3;
    private final static int             MOVE                                               = 20;
    private final static String          NODE_POPMENU_NODE_CLIENT_PROPERTY                  = "node";
    private static final float           ONEHALF_PI                                         = ( float ) ( 1.5 * Math.PI );
    private static final short           OV_BORDER                                          = 10;
    private final static double          OVERVIEW_FOUND_NODE_BOX_SIZE                       = 2;
    private final static double          OVERVIEW_FOUND_NODE_BOX_SIZE_HALF                  = 1;
    private static final float           PI                                                 = ( float ) ( Math.PI );
    final private static Font            POPUP_FONT                                         = new Font( Configuration.getDefaultFontFamilyName(),
                                                                                                        Font.PLAIN,
                                                                                                        12 );
    private static final float           ROUNDED_D                                          = 8;
    private final static long            serialVersionUID                                   = -978349745916505029L;
    private static final BasicStroke     STROKE_005                                         = new BasicStroke( 0.05f );
    private static final BasicStroke     STROKE_01                                          = new BasicStroke( 0.1f );
    private static final BasicStroke     STROKE_025                                         = new BasicStroke( 0.25f );
    private static final BasicStroke     STROKE_05                                          = new BasicStroke( 0.5f );
    private static final BasicStroke     STROKE_075                                         = new BasicStroke( 0.75f );
    private static final BasicStroke     STROKE_1                                           = new BasicStroke( 1f );
    private static final BasicStroke     STROKE_2                                           = new BasicStroke( 2f );
    private static final double          TWO_PI                                             = 2 * Math.PI;
    private final static int             WIGGLE                                             = 2;
    HashMap<Long, Short>                 _nodeid_dist_to_leaf                               = new HashMap<Long, Short>();
    final private Arc2D                  _arc                                               = new Arc2D.Double();
    private AffineTransform              _at;
    private int                          _circ_max_depth;
    final private Set<Long>              _collapsed_external_nodeid_set                     = new HashSet<Long>();
    private JColorChooser                _color_chooser                                     = null;
    private Configuration                _configuration                                     = null;
    private ControlPanel                 _control_panel                                     = null;
    private final CubicCurve2D           _cubic_curve                                       = new CubicCurve2D.Float();
    private Set<Long>                    _current_external_nodes                            = null;
    private StringBuilder                _current_external_nodes_data_buffer                = new StringBuilder();
    private int                          _current_external_nodes_data_buffer_change_counter = 0;
    private int                          _domain_structure_e_value_thr_exp                  = Constants.DOMAIN_STRUCTURE_E_VALUE_THR_DEFAULT_EXP;
    private double                       _domain_structure_width                            = Constants.DOMAIN_STRUCTURE_DEFAULT_WIDTH;
    private int                          _dynamic_hiding_factor                             = 0;
    private boolean                      _edited                                            = false;
    private final Ellipse2D              _ellipse                                           = new Ellipse2D.Float();
    private int                          _external_node_index                               = 0;
    private Set<Long>                    _found_nodes_0                                     = null;
    private Set<Long>                    _found_nodes_1                                     = null;
    private final FontRenderContext      _frc                                               = new FontRenderContext( null,
                                                                                                                     false,
                                                                                                                     false );
    private PHYLOGENY_GRAPHICS_TYPE      _graphics_type                                     = PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR;
    private PhylogenyNode                _highlight_node                                    = null;
    private boolean                      _in_ov                                             = false;
    private boolean                      _in_ov_rect                                        = false;
    private float                        _last_drag_point_x                                 = 0;
    private float                        _last_drag_point_y                                 = 0;
    private final Line2D                 _line                                              = new Line2D.Float();
    private int                          _longest_ext_node_info                             = 0;
    private MainPanel                    _main_panel                                        = null;
    private double                       _max_distance_to_root                              = -1;
    private Popup                        _node_desc_popup;
    private int                          _node_frame_index                                  = 0;
    private final NodeFrame[]            _node_frames                                       = new NodeFrame[ TreePanel.MAX_NODE_FRAMES ];
    private JPopupMenu                   _node_popup_menu                                   = null;
    private JMenuItem                    _node_popup_menu_items[]                           = null;
    private PhylogenyNode[]              _nodes_in_preorder                                 = null;
    private Options                      _options                                           = null;
    private float                        _ov_max_height                                     = 0;
    private float                        _ov_max_width                                      = 0;
    private boolean                      _ov_on                                             = false;
    private final Rectangle2D            _ov_rectangle                                      = new Rectangle2D.Float();
    private final Rectangle              _ov_virtual_rectangle                              = new Rectangle();
    private float                        _ov_x_correction_factor                            = 0.0f;
    private float                        _ov_x_distance                                     = 0;
    private int                          _ov_x_position                                     = 0;
    private float                        _ov_y_distance                                     = 0;
    private int                          _ov_y_position                                     = 0;
    private int                          _ov_y_start                                        = 0;
    private final boolean                _phy_has_branch_lengths;
    private Phylogeny                    _phylogeny                                         = null;
    private final Path2D.Float           _polygon                                           = new Path2D.Float();
    private final StringBuffer           _popup_buffer                                      = new StringBuffer();
    private final QuadCurve2D            _quad_curve                                        = new QuadCurve2D.Float();
    private Sequence                     _query_sequence                                    = null;
    private final Rectangle2D            _rectangle                                         = new Rectangle2D.Float();
    private final RenderingHints         _rendering_hints                                   = new RenderingHints( RenderingHints.KEY_RENDERING,
                                                                                                                  RenderingHints.VALUE_RENDER_DEFAULT );
    private JTextArea                    _rollover_popup;
    private PhylogenyNode                _root;
    private final StringBuilder          _sb                                                = new StringBuilder();
    private double                       _scale_distance                                    = 0.0;
    private String                       _scale_label                                       = null;
    // expression values menu:
    private DescriptiveStatistics        _statistics_for_vector_data;
    private final Phylogeny[]            _sub_phylogenies                                   = new Phylogeny[ TreePanel.MAX_SUBTREES ];
    private final PhylogenyNode[]        _sub_phylogenies_temp_roots                        = new PhylogenyNode[ TreePanel.MAX_SUBTREES ];
    private int                          _subtree_index                                     = 0;
    private File                         _treefile                                          = null;
    private float                        _urt_factor                                        = 1;
    private float                        _urt_factor_ov                                     = 1;
    final private HashMap<Long, Double>  _urt_nodeid_angle_map                              = new HashMap<Long, Double>();
    final private HashMap<Long, Integer> _urt_nodeid_index_map                              = new HashMap<Long, Integer>();
    private double                       _urt_starting_angle                                = ( float ) ( Math.PI / 2 );
    private float                        _x_correction_factor                               = 0.0f;
    private float                        _x_distance                                        = 0.0f;
    private float                        _y_distance                                        = 0.0f;
    //  private Image                           offscreenImage;
    //  private Graphics                        offscreenGraphics;
    //  private Dimension                       offscreenDimension;
    static {
        final DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator( '.' );
        FORMATTER_CONFIDENCE = new DecimalFormat( "#.###", dfs );
        FORMATTER_BRANCH_LENGTH = new DecimalFormat( "#.###", dfs );
    }

    TreePanel( final Phylogeny t, final Configuration configuration, final MainPanel tjp ) {
        requestFocusInWindow();
        addKeyListener( new KeyAdapter() {

            @Override
            public void keyPressed( final KeyEvent key_event ) {
                keyPressedCalls( key_event );
                requestFocusInWindow();
            }
        } );
        addFocusListener( new FocusAdapter() {

            @Override
            public void focusGained( final FocusEvent e ) {
                requestFocusInWindow();
            }
        } );
        if ( ( t == null ) || t.isEmpty() ) {
            throw new IllegalArgumentException( "attempt to draw phylogeny which is null or empty" );
        }
        _graphics_type = tjp.getOptions().getPhylogenyGraphicsType();
        _main_panel = tjp;
        _configuration = configuration;
        _phylogeny = t;
        _phy_has_branch_lengths = AptxUtil.isHasAtLeastOneBranchLengthLargerThanZero( _phylogeny );
        init();
        // if ( !_phylogeny.isEmpty() ) {
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        checkForVectorProperties( _phylogeny );
        // }
        setBackground( getTreeColorSet().getBackgroundColor() );
        final MouseListener mouse_listener = new MouseListener( this );
        addMouseListener( mouse_listener );
        addMouseMotionListener( mouse_listener );
        addMouseWheelListener( this );
        calculateScaleDistance();
        FORMATTER_CONFIDENCE.setMaximumFractionDigits( configuration.getNumberOfDigitsAfterCommaForConfidenceValues() );
        FORMATTER_BRANCH_LENGTH.setMaximumFractionDigits( configuration
                .getNumberOfDigitsAfterCommaForBranchLengthValues() );
    }

    @Override
    final public void actionPerformed( final ActionEvent e ) {
        boolean done = false;
        final JMenuItem node_popup_menu_item = ( JMenuItem ) e.getSource();
        for( int index = 0; ( index < _node_popup_menu_items.length ) && !done; index++ ) {
            // NOTE: index corresponds to the indices of click-to options
            // in the control panel.
            if ( node_popup_menu_item == _node_popup_menu_items[ index ] ) {
                // Set this as the new default click-to action
                _main_panel.getControlPanel().setClickToAction( index );
                final PhylogenyNode node = ( PhylogenyNode ) _node_popup_menu
                        .getClientProperty( NODE_POPMENU_NODE_CLIENT_PROPERTY );
                handleClickToAction( _control_panel.getActionWhenNodeClicked(), node );
                done = true;
            }
        }
        repaint();
        requestFocusInWindow();
    }

    public synchronized Hashtable<String, BufferedImage> getImageMap() {
        return getMainPanel().getImageMap();
    }

    final public MainPanel getMainPanel() {
        return _main_panel;
    }

    /**
     * Get a pointer to the phylogeny 
     * 
     * @return a pointer to the phylogeny
     */
    public final Phylogeny getPhylogeny() {
        return _phylogeny;
    }

    @Override
    final public void mouseWheelMoved( final MouseWheelEvent e ) {
        final int notches = e.getWheelRotation();
        if ( inOvVirtualRectangle( e ) ) {
            if ( !isInOvRect() ) {
                setInOvRect( true );
                repaint();
            }
        }
        else {
            if ( isInOvRect() ) {
                setInOvRect( false );
                repaint();
            }
        }
        if ( e.isControlDown() ) {
            if ( notches < 0 ) {
                getTreeFontSet().increaseFontSize();
                getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
            else {
                getTreeFontSet().decreaseFontSize( 1, false );
                getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
        }
        else if ( e.isShiftDown() ) {
            if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                    || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                if ( notches < 0 ) {
                    for( int i = 0; i < ( -notches ); ++i ) {
                        setStartingAngle( ( getStartingAngle() % TWO_PI ) + ANGLE_ROTATION_UNIT );
                        getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
                else {
                    for( int i = 0; i < notches; ++i ) {
                        setStartingAngle( ( getStartingAngle() % TWO_PI ) - ANGLE_ROTATION_UNIT );
                        if ( getStartingAngle() < 0 ) {
                            setStartingAngle( TWO_PI + getStartingAngle() );
                        }
                        getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
            }
            else {
                if ( notches < 0 ) {
                    for( int i = 0; i < ( -notches ); ++i ) {
                        getControlPanel().zoomInY( Constants.WHEEL_ZOOM_IN_FACTOR );
                        getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
                else {
                    for( int i = 0; i < notches; ++i ) {
                        getControlPanel().zoomOutY( Constants.WHEEL_ZOOM_OUT_FACTOR );
                        getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
            }
        }
        else {
            if ( notches < 0 ) {
                for( int i = 0; i < ( -notches ); ++i ) {
                    getControlPanel().zoomInX( Constants.WHEEL_ZOOM_IN_FACTOR,
                                               Constants.WHEEL_ZOOM_IN_X_CORRECTION_FACTOR );
                    getControlPanel().zoomInY( Constants.WHEEL_ZOOM_IN_FACTOR );
                    getControlPanel().displayedPhylogenyMightHaveChanged( false );
                }
            }
            else {
                for( int i = 0; i < notches; ++i ) {
                    getControlPanel().zoomOutY( Constants.WHEEL_ZOOM_OUT_FACTOR );
                    getControlPanel().zoomOutX( Constants.WHEEL_ZOOM_OUT_FACTOR,
                                                Constants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR );
                    getControlPanel().displayedPhylogenyMightHaveChanged( false );
                }
            }
        }
        requestFocus();
        requestFocusInWindow();
        requestFocus();
    }

    @Override
    final public void paintComponent( final Graphics g ) {
        // Dimension currentSize = getSize();
        //  if ( offscreenImage == null || !currentSize.equals( offscreenDimension ) ) {
        // call the 'java.awt.Component.createImage(...)' method to get an
        // image
        //   offscreenImage = createImage( currentSize.width, currentSize.height );
        //  offscreenGraphics = offscreenImage.getGraphics();
        //  offscreenDimension = currentSize;
        // }
        // super.paintComponent( g ); //why?
        //final Graphics2D g2d = ( Graphics2D ) offscreenGraphics;
        final Graphics2D g2d = ( Graphics2D ) g;
        g2d.setRenderingHints( _rendering_hints );
        paintPhylogeny( g2d, false, false, 0, 0, 0, 0 );
        //g.drawImage( offscreenImage, 0, 0, this );
    }

    @Override
    final public int print( final Graphics g, final PageFormat page_format, final int page_index )
            throws PrinterException {
        if ( page_index > 0 ) {
            return ( NO_SUCH_PAGE );
        }
        else {
            final Graphics2D g2d = ( Graphics2D ) g;
            g2d.translate( page_format.getImageableX(), page_format.getImageableY() );
            // Turn off double buffering !?
            paintPhylogeny( g2d, true, false, 0, 0, 0, 0 );
            // Turn double buffering back on !?
            return ( PAGE_EXISTS );
        }
    }

    public final void setEdited( final boolean edited ) {
        _edited = edited;
    }

    public synchronized void setImageMap( final Hashtable<String, BufferedImage> image_map ) {
        getMainPanel().setImageMap( image_map );
    }

    /**
     * Set a phylogeny tree.
     * 
     * @param t
     *            an instance of a Phylogeny
     */
    public final void setTree( final Phylogeny t ) {
        setNodeInPreorderToNull();
        _phylogeny = t;
    }

    public final void setWaitCursor() {
        setCursor( WAIT_CURSOR );
        repaint();
    }

    @Override
    public void update( final Graphics g ) {
        paint( g );
    }

    final void calcMaxDepth() {
        if ( _phylogeny != null ) {
            _circ_max_depth = PhylogenyMethods.calculateMaxDepth( _phylogeny );
        }
    }

    /**
     * Set parameters for printing the displayed tree
     * 
     */
    final void calcParametersForPainting( final int x, final int y, final boolean recalc_longest_ext_node_info ) {
        // updateStyle(); not needed?
        if ( ( _phylogeny != null ) && !_phylogeny.isEmpty() ) {
            initNodeData();
            if ( recalc_longest_ext_node_info ) {
                calculateLongestExtNodeInfo();
                if ( getOptions().isAllowFontSizeChange() ) {
                    if ( ( getLongestExtNodeInfo() > ( x * 0.6 ) )
                            && ( getTreeFontSet().getLargeFont().getSize() > 2 + TreeFontSet.FONT_SIZE_CHANGE_STEP ) ) {
                        while ( ( getLongestExtNodeInfo() > ( x * 0.7 ) )
                                && ( getTreeFontSet().getLargeFont().getSize() > 2 ) ) {
                            getMainPanel().getTreeFontSet().decreaseFontSize( getConfiguration().getMinBaseFontSize(),
                                                                              true );
                            calculateLongestExtNodeInfo();
                        }
                    }
                    else {
                        while ( ( getLongestExtNodeInfo() < ( x * 0.6 ) )
                                && ( getTreeFontSet().getLargeFont().getSize() <= getTreeFontSet().getLargeFontMemory()
                                        .getSize() - TreeFontSet.FONT_SIZE_CHANGE_STEP ) ) {
                            getMainPanel().getTreeFontSet().increaseFontSize();
                            calculateLongestExtNodeInfo();
                        }
                    }
                }
            }
            int ext_nodes = _phylogeny.getRoot().getNumberOfExternalNodes();
            final int max_depth = PhylogenyMethods.calculateMaxDepth( _phylogeny );
            if ( ext_nodes == 1 ) {
                ext_nodes = max_depth;
                if ( ext_nodes < 1 ) {
                    ext_nodes = 1;
                }
            }
            updateOvSizes();
            float xdist = 0;
            float ov_xdist = 0;
            if ( !isNonLinedUpCladogram() && !isUniformBranchLengthsForCladogram() ) {
                xdist = ( float ) ( ( x - getLongestExtNodeInfo() - TreePanel.MOVE ) / ( ext_nodes + 3.0 ) );
                ov_xdist = ( float ) ( getOvMaxWidth() / ( ext_nodes + 3.0 ) );
            }
            else {
                xdist = ( ( x - getLongestExtNodeInfo() - TreePanel.MOVE ) / ( max_depth + 1 ) );
                ov_xdist = ( getOvMaxWidth() / ( max_depth + 1 ) );
            }
            float ydist = ( float ) ( ( y - TreePanel.MOVE ) / ( ext_nodes * 2.0 ) );
            if ( xdist < 0.0 ) {
                xdist = 0.0f;
            }
            if ( ov_xdist < 0.0 ) {
                ov_xdist = 0.0f;
            }
            if ( ydist < 0.0 ) {
                ydist = 0.0f;
            }
            setXdistance( xdist );
            setYdistance( ydist );
            setOvXDistance( ov_xdist );
            final double height = _phylogeny.getHeight();
            if ( height > 0 ) {
                final float corr = ( float ) ( ( x - TreePanel.MOVE - getLongestExtNodeInfo() - getXdistance() ) / height );
                setXcorrectionFactor( corr > 0 ? corr : 0 );
                final float ov_corr = ( float ) ( ( getOvMaxWidth() - getOvXDistance() ) / height );
                setOvXcorrectionFactor( ov_corr > 0 ? ov_corr : 0 );
            }
            else {
                setXcorrectionFactor( 0 );
                setOvXcorrectionFactor( 0 );
            }
            _circ_max_depth = max_depth;
            setUpUrtFactor();
            //
            if ( getOptions().isAllowFontSizeChange() ) {
                if ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                        && ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                    //                int dynamic_hiding_factor = calcDynamicHidingFactor();
                    //                if ( dynamic_hiding_factor > 1 ) {
                    //                    while ( dynamic_hiding_factor > 1
                    //                            && getTreeFontSet()._fm_large.getHeight() > TreeFontSet.SMALL_FONTS_BASE ) {
                    //                        getTreeFontSet().decreaseFontSize( 1, true );
                    //                        dynamic_hiding_factor = calcDynamicHidingFactor();
                    //                    }
                    //                }
                    //                else if ( getTreeFontSet().isDecreasedSizeBySystem() ) {
                    //                    while ( dynamic_hiding_factor < 1 && getTreeFontSet()._fm_large.getHeight() < 12 ) {
                    //                        getTreeFontSet().increaseFontSize();
                    //                        dynamic_hiding_factor = calcDynamicHidingFactor();
                    //                    }
                    //                }
                }
            }
            //
        }
    }

    final void calculateLongestExtNodeInfo() {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return;
        }
        int max_length = ForesterUtil.roundToInt( ( getSize().getWidth() - MOVE )
                * Constants.EXT_NODE_INFO_LENGTH_MAX_RATIO );
        if ( max_length < 40 ) {
            max_length = 40;
        }
        int longest = 30;
        for( final PhylogenyNode node : _phylogeny.getExternalNodes() ) {
            int sum = 0;
            if ( node.isCollapse() ) {
                continue;
            }
            final StringBuilder sb = new StringBuilder();
            nodeDataAsSB( node, sb );
            if ( node.getNodeData().isHasTaxonomy() ) {
                nodeTaxonomyDataAsSB( node.getNodeData().getTaxonomy(), sb );
            }
            boolean use_vis = false;
            final Graphics2D g = ( Graphics2D ) getGraphics();
            if ( getControlPanel().isUseVisualStyles() ) {
                use_vis = setFont( g, node, false );
            }
            if ( !use_vis ) {
                sum = getFontMetricsForLargeDefaultFont().stringWidth( sb.toString() );
            }
            else {
                sum = getFontMetrics( g.getFont() ).stringWidth( sb.toString() );
            }
            if ( getControlPanel().isShowBinaryCharacters() && node.getNodeData().isHasBinaryCharacters() ) {
                sum += getFontMetricsForLargeDefaultFont().stringWidth( node.getNodeData().getBinaryCharacters()
                        .getGainedCharactersAsStringBuffer().toString() );
            }
            if ( sum >= max_length ) {
                setLongestExtNodeInfo( max_length );
                return;
            }
            if ( sum > longest ) {
                longest = sum;
            }
        }
        if ( longest >= max_length ) {
            setLongestExtNodeInfo( max_length );
        }
        else {
            setLongestExtNodeInfo( longest );
        }
    }

    final void calculateScaleDistance() {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return;
        }
        final double height = getMaxDistanceToRoot();
        if ( height > 0 ) {
            if ( ( height <= 0.5 ) ) {
                setScaleDistance( 0.01 );
            }
            else if ( height <= 5.0 ) {
                setScaleDistance( 0.1 );
            }
            else if ( height <= 50.0 ) {
                setScaleDistance( 1 );
            }
            else if ( height <= 500.0 ) {
                setScaleDistance( 10 );
            }
            else {
                setScaleDistance( 100 );
            }
        }
        else {
            setScaleDistance( 0.0 );
        }
        String scale_label = String.valueOf( getScaleDistance() );
        if ( !ForesterUtil.isEmpty( _phylogeny.getDistanceUnit() ) ) {
            scale_label += " [" + _phylogeny.getDistanceUnit() + "]";
        }
        setScaleLabel( scale_label );
    }

    final Color calculateTaxonomyBasedColor( final Taxonomy tax ) {
        if ( getOptions().isColorByTaxonomicGroup() ) {
            if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                boolean ex = false;
                String group = null;
                try {
                    group = TaxonomyUtil.getTaxGroupByTaxCode( tax.getTaxonomyCode() );
                }
                catch ( final Exception e ) {
                    ex = true;
                }
                if ( !ex && !ForesterUtil.isEmpty( group ) ) {
                    final Color c = ForesterUtil.obtainColorDependingOnTaxonomyGroup( group );
                    if ( c != null ) {
                        return c;
                    }
                }
            }
            return getTreeColorSet().getTaxonomyColor();
        }
        else {
            if ( ForesterUtil.isEmpty( tax.getTaxonomyCode() ) && ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                return getTreeColorSet().getTaxonomyColor();
            }
            Color c = null;
            if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                c = getControlPanel().getSpeciesColors().get( tax.getTaxonomyCode() );
            }
            if ( ( c == null ) && !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                c = getControlPanel().getSpeciesColors().get( tax.getScientificName() );
            }
            if ( c == null ) {
                if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                    c = TreePanelUtil.calculateColorFromString( tax.getTaxonomyCode(), true );
                    getControlPanel().getSpeciesColors().put( tax.getTaxonomyCode(), c );
                }
                else {
                    c = TreePanelUtil.calculateColorFromString( tax.getScientificName(), true );
                    getControlPanel().getSpeciesColors().put( tax.getScientificName(), c );
                }
            }
            return c;
        }
    }
    
    final Color calculateSequenceBasedColor( final Sequence seq ) {
            if ( ForesterUtil.isEmpty( seq.getName() ) ) {
                return getTreeColorSet().getSequenceColor();
            }
            Color c = null;
            final String seq_name = seq.getName();
            c = getControlPanel().getSequenceColors().get( seq_name  );
            if ( c == null ) {
                    c = TreePanelUtil.calculateColorFromString( seq_name, false );
                    getControlPanel().getSequenceColors().put( seq_name, c );
            }
            return c;
    }

    void checkForVectorProperties( final Phylogeny phy ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final PhylogenyNodeIterator iter = phy.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.getNodeData().getProperties() != null ) {
                final PropertiesMap pm = node.getNodeData().getProperties();
                final double[] vector = new double[ pm.getProperties().size() ];
                int counter = 0;
                for( final String ref : pm.getProperties().keySet() ) {
                    if ( ref.startsWith( PhyloXmlUtil.VECTOR_PROPERTY_REF ) ) {
                        final Property p = pm.getProperty( ref );
                        final String value_str = p.getValue();
                        final String index_str = ref
                                .substring( PhyloXmlUtil.VECTOR_PROPERTY_REF.length(), ref.length() );
                        double d = -100;
                        try {
                            d = Double.parseDouble( value_str );
                        }
                        catch ( final NumberFormatException e ) {
                            JOptionPane.showMessageDialog( this, "Could not parse \"" + value_str
                                    + "\" into a decimal value", "Problem with Vector Data", JOptionPane.ERROR_MESSAGE );
                            return;
                        }
                        int i = -1;
                        try {
                            i = Integer.parseInt( index_str );
                        }
                        catch ( final NumberFormatException e ) {
                            JOptionPane.showMessageDialog( this,
                                                           "Could not parse \"" + index_str
                                                                   + "\" into index for vector data",
                                                           "Problem with Vector Data",
                                                           JOptionPane.ERROR_MESSAGE );
                            return;
                        }
                        if ( i < 0 ) {
                            JOptionPane.showMessageDialog( this,
                                                           "Attempt to use negative index for vector data",
                                                           "Problem with Vector Data",
                                                           JOptionPane.ERROR_MESSAGE );
                            return;
                        }
                        vector[ i ] = d;
                        ++counter;
                        stats.addValue( d );
                    }
                }
                final List<Double> vector_l = new ArrayList<Double>( counter );
                for( int i = 0; i < counter; ++i ) {
                    vector_l.add( vector[ i ] );
                }
                node.getNodeData().setVector( vector_l );
            }
        }
        if ( stats.getN() > 0 ) {
            _statistics_for_vector_data = stats;
        }
    }

    void clearCurrentExternalNodesDataBuffer() {
        setCurrentExternalNodesDataBuffer( new StringBuilder() );
    }

    /**
     * Collapse the tree from the given node
     * 
     * @param node
     *            a PhylogenyNode
     */
    final void collapse( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot collapse in unrooted display type",
                                           "Attempt to collapse in unrooted display",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( !node.isExternal() && !node.isRoot() ) {
            final boolean collapse = !node.isCollapse();
            TreePanelUtil.collapseSubtree( node, collapse );
            updateSetOfCollapsedExternalNodes();
            _phylogeny.recalculateNumberOfExternalDescendants( true );
            resetNodeIdToDistToLeafMap();
            calculateLongestExtNodeInfo();
            setNodeInPreorderToNull();
            _control_panel.displayedPhylogenyMightHaveChanged( true );
            resetPreferredSize();
            updateOvSizes();
            _main_panel.adjustJScrollPane();
            repaint();
        }
    }

    final void collapseSpeciesSpecificSubtrees() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setWaitCursor();
        TreePanelUtil.collapseSpeciesSpecificSubtrees( _phylogeny );
        updateSetOfCollapsedExternalNodes();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        calculateLongestExtNodeInfo();
        setNodeInPreorderToNull();
        resetPreferredSize();
        _main_panel.adjustJScrollPane();
        setArrowCursor();
        repaint();
    }

    final void colorRank( final String rank ) {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setWaitCursor();
        AptxUtil.removeBranchColors( _phylogeny );
        final int colorizations = TreePanelUtil.colorPhylogenyAccordingToRanks( _phylogeny, rank, this );
        if ( colorizations > 0 ) {
            _control_panel.setColorBranches( true );
            if ( _control_panel.getUseVisualStylesCb() != null ) {
                _control_panel.getUseVisualStylesCb().setSelected( true );
            }
            if ( _control_panel.getColorAccSpeciesCb() != null ) {
                _control_panel.getColorAccSpeciesCb().setSelected( false );
            }
           
            _options.setColorLabelsSameAsParentBranch( true );
            _control_panel.repaint();
        }
        setArrowCursor();
        repaint();
        if ( colorizations > 0 ) {
            String msg = "Taxonomy colorization via " + rank + " completed:\n";
            if ( colorizations > 1 ) {
                msg += "colorized " + colorizations + " subtrees";
            }
            else {
                msg += "colorized one subtree";
            }
            setEdited( true );
            JOptionPane.showMessageDialog( this,
                                           msg,
                                           "Taxonomy Colorization Completed (" + rank + ")",
                                           JOptionPane.INFORMATION_MESSAGE );
        }
        else {
            String msg = "Could not taxonomy colorize any subtree via " + rank + ".\n";
            msg += "Possible solutions (given that suitable taxonomic information is present):\n";
            msg += "select a different rank (e.g. phylum, genus, ...)\n";
            msg += "  and/or\n";
            msg += "execute:\n";
            msg += "1. \"" + MainFrameApplication.OBTAIN_DETAILED_TAXONOMIC_INFORMATION + "\" (Tools)\n";
            msg += "2. \"" + MainFrameApplication.INFER_ANCESTOR_TAXONOMIES + "\" (Analysis)";
            JOptionPane.showMessageDialog( this, msg, "Taxonomy Colorization Failed", JOptionPane.WARNING_MESSAGE );
        }
    }

    final void confColor() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setWaitCursor();
        AptxUtil.removeBranchColors( _phylogeny );
        TreePanelUtil.colorPhylogenyAccordingToConfidenceValues( _phylogeny, this );
        _control_panel.setColorBranches( true );
        if ( _control_panel.getUseVisualStylesCb() != null ) {
            _control_panel.getUseVisualStylesCb().setSelected( true );
        }
        setArrowCursor();
        repaint();
    }

    final void decreaseDomainStructureEvalueThreshold() {
        if ( _domain_structure_e_value_thr_exp > -20 ) {
            _domain_structure_e_value_thr_exp -= 1;
        }
    }

    /**
     * Find the node, if any, at the given location
     * 
     * @param x
     * @param y
     * @return pointer to the node at x,y, null if not found
     */
    final PhylogenyNode findNode( final int x, final int y ) {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return null;
        }
        final int half_box_size_plus_wiggle = ( getOptions().getDefaultNodeShapeSize() / 2 ) + WIGGLE;
        for( final PhylogenyNodeIterator iter = _phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( ( _phylogeny.isRooted() || !node.isRoot() || ( node.getNumberOfDescendants() > 2 ) )
                    && ( ( node.getXcoord() - half_box_size_plus_wiggle ) <= x )
                    && ( ( node.getXcoord() + half_box_size_plus_wiggle ) >= x )
                    && ( ( node.getYcoord() - half_box_size_plus_wiggle ) <= y )
                    && ( ( node.getYcoord() + half_box_size_plus_wiggle ) >= y ) ) {
                return node;
            }
        }
        return null;
    }

    final Configuration getConfiguration() {
        return _configuration;
    }

    final ControlPanel getControlPanel() {
        return _control_panel;
    }

    String getCurrentExternalNodesDataBufferAsString() {
        return _current_external_nodes_data_buffer.toString();
    }

    int getCurrentExternalNodesDataBufferChangeCounter() {
        return _current_external_nodes_data_buffer_change_counter;
    }

    final int getDomainStructureEvalueThreshold() {
        return _domain_structure_e_value_thr_exp;
    }

    final Set<Long> getFoundNodes0() {
        return _found_nodes_0;
    }

    final Set<Long> getFoundNodes1() {
        return _found_nodes_1;
    }

    final Color getGraphicsForNodeBoxWithColorForParentBranch( final PhylogenyNode node ) {
        if ( getControlPanel().isUseVisualStyles() && ( PhylogenyMethods.getBranchColorValue( node ) != null ) ) {
            return ( PhylogenyMethods.getBranchColorValue( node ) );
        }
        else {
            return ( getTreeColorSet().getBranchColor() );
        }
    }

    final int getLongestExtNodeInfo() {
        return _longest_ext_node_info;
    }

    final Options getOptions() {
        if ( _options == null ) {
            _options = getControlPanel().getOptions();
        }
        return _options;
    }

    final Rectangle2D getOvRectangle() {
        return _ov_rectangle;
    }

    final Rectangle getOvVirtualRectangle() {
        return _ov_virtual_rectangle;
    }

    final PHYLOGENY_GRAPHICS_TYPE getPhylogenyGraphicsType() {
        return _graphics_type;
    }

    final double getStartingAngle() {
        return _urt_starting_angle;
    }

    DescriptiveStatistics getStatisticsForExpressionValues() {
        return _statistics_for_vector_data;
    }

  
    final Color getTaxonomyBasedColor( final PhylogenyNode node ) {
        if ( node.getNodeData().isHasTaxonomy() ) {
            return calculateTaxonomyBasedColor( node.getNodeData().getTaxonomy() );
        }
        // return non-colorized color
        return getTreeColorSet().getTaxonomyColor();
    }
    

    final Color getSequenceBasedColor( final PhylogenyNode node ) {
        if ( node.getNodeData().isHasSequence() ) {
            return calculateSequenceBasedColor( node.getNodeData().getSequence() );
        }
        // return non-colorized color
        return getTreeColorSet().getSequenceColor();
    }
    

    /**
     * @return pointer to colorset for tree drawing
     */
    final TreeColorSet getTreeColorSet() {
        return getMainPanel().getTreeColorSet();
    }

    final File getTreeFile() {
        return _treefile;
    }

    final float getXcorrectionFactor() {
        return _x_correction_factor;
    }

    final float getXdistance() {
        return _x_distance;
    }

    final float getYdistance() {
        return _y_distance;
    }

    final void increaseDomainStructureEvalueThreshold() {
        if ( _domain_structure_e_value_thr_exp < 3 ) {
            _domain_structure_e_value_thr_exp += 1;
        }
    }

    final void initNodeData() {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return;
        }
        double max_original_domain_structure_width = 0.0;
        for( final PhylogenyNode node : _phylogeny.getExternalNodes() ) {
            if ( node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                RenderableDomainArchitecture rds = null;
                if ( !( node.getNodeData().getSequence().getDomainArchitecture() instanceof RenderableDomainArchitecture ) ) {
                    rds = new RenderableDomainArchitecture( node.getNodeData().getSequence().getDomainArchitecture(),
                                                            getConfiguration() );
                    node.getNodeData().getSequence().setDomainArchitecture( rds );
                }
                else {
                    rds = ( RenderableDomainArchitecture ) node.getNodeData().getSequence().getDomainArchitecture();
                }
                if ( getControlPanel().isShowDomainArchitectures() ) {
                    final double dsw = rds.getOriginalSize().getWidth();
                    if ( dsw > max_original_domain_structure_width ) {
                        max_original_domain_structure_width = dsw;
                    }
                }
            }
        }
        if ( getControlPanel().isShowDomainArchitectures() ) {
            final double ds_factor_width = _domain_structure_width / max_original_domain_structure_width;
            for( final PhylogenyNode node : _phylogeny.getExternalNodes() ) {
                if ( node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                    final RenderableDomainArchitecture rds = ( RenderableDomainArchitecture ) node.getNodeData()
                            .getSequence().getDomainArchitecture();
                    rds.setRenderingFactorWidth( ds_factor_width );
                    rds.setParameter( _domain_structure_e_value_thr_exp );
                }
            }
        }
    }

    final boolean inOv( final MouseEvent e ) {
        return ( ( e.getX() > ( getVisibleRect().x + getOvXPosition() + 1 ) )
                && ( e.getX() < ( ( getVisibleRect().x + getOvXPosition() + getOvMaxWidth() ) - 1 ) )
                && ( e.getY() > ( getVisibleRect().y + getOvYPosition() + 1 ) ) && ( e.getY() < ( ( getVisibleRect().y
                + getOvYPosition() + getOvMaxHeight() ) - 1 ) ) );
    }

    final boolean inOvRectangle( final MouseEvent e ) {
        return ( ( e.getX() >= ( getOvRectangle().getX() - 1 ) )
                && ( e.getX() <= ( getOvRectangle().getX() + getOvRectangle().getWidth() + 1 ) )
                && ( e.getY() >= ( getOvRectangle().getY() - 1 ) ) && ( e.getY() <= ( getOvRectangle().getY()
                + getOvRectangle().getHeight() + 1 ) ) );
    }

    final boolean isApplet() {
        return getMainPanel() instanceof MainPanelApplets;
    }

    final boolean isCanCollapse() {
        return ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
    }

    final boolean isCanColorSubtree() {
        return ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
    }

    final boolean isCanCopy() {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && getOptions().isEditable() );
    }

    final boolean isCanCut( final PhylogenyNode node ) {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && getOptions().isEditable() && !node
                .isRoot() );
    }

    final boolean isCanDelete() {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && getOptions().isEditable() );
    }

    final boolean isCanPaste() {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && getOptions().isEditable()
                && ( getCutOrCopiedTree() != null ) && !getCutOrCopiedTree().isEmpty() );
    }

    final boolean isCanReroot() {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && ( _subtree_index < 1 ) );
    }

    final boolean isCanSubtree( final PhylogenyNode node ) {
        return ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) && !node.isExternal() && ( !node
                .isRoot() || ( _subtree_index > 0 ) ) );
    }

    final boolean isCurrentTreeIsSubtree() {
        return ( _subtree_index > 0 );
    }

    final boolean isEdited() {
        return _edited;
    }

    final boolean isInOvRect() {
        return _in_ov_rect;
    }

    final boolean isOvOn() {
        return _ov_on;
    }

    final boolean isPhyHasBranchLengths() {
        return _phy_has_branch_lengths;
    }

    final void midpointRoot() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        if ( !_phylogeny.isRerootable() ) {
            JOptionPane.showMessageDialog( this,
                                           "This is not rerootable",
                                           "Not rerootable",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        setNodeInPreorderToNull();
        setWaitCursor();
        PhylogenyMethods.midpointRoot( _phylogeny );
        resetNodeIdToDistToLeafMap();
        setArrowCursor();
        setEdited( true );
        repaint();
    }

    final void mouseClicked( final MouseEvent e ) {
        if ( getOptions().isShowOverview() && isOvOn() && isInOv() ) {
            final double w_ratio = getVisibleRect().width / getOvRectangle().getWidth();
            final double h_ratio = getVisibleRect().height / getOvRectangle().getHeight();
            double x = ( e.getX() - getVisibleRect().x - getOvXPosition() - ( getOvRectangle().getWidth() / 2.0 ) )
                    * w_ratio;
            double y = ( e.getY() - getVisibleRect().y - getOvYPosition() - ( getOvRectangle().getHeight() / 2.0 ) )
                    * h_ratio;
            if ( x < 0 ) {
                x = 0;
            }
            if ( y < 0 ) {
                y = 0;
            }
            final double max_x = getWidth() - getVisibleRect().width;
            final double max_y = getHeight() - getVisibleRect().height;
            if ( x > max_x ) {
                x = max_x;
            }
            if ( y > max_y ) {
                y = max_y;
            }
            getMainPanel().getCurrentScrollPane().getViewport()
                    .setViewPosition( new Point( ForesterUtil.roundToInt( x ), ForesterUtil.roundToInt( y ) ) );
            setInOvRect( true );
            repaint();
        }
        else {
            final PhylogenyNode node = findNode( e.getX(), e.getY() );
            if ( node != null ) {
                if ( !node.isRoot() && node.getParent().isCollapse() ) {
                    return;
                }
                _highlight_node = node;
                // Check if shift key is down
                if ( ( e.getModifiers() & InputEvent.SHIFT_MASK ) != 0 ) {
                    // Yes, so add to _found_nodes
                    if ( getFoundNodes0() == null ) {
                        setFoundNodes0( new HashSet<Long>() );
                    }
                    getFoundNodes0().add( node.getId() );
                    // Check if control key is down
                }
                else if ( ( e.getModifiers() & InputEvent.CTRL_MASK ) != 0 ) {
                    // Yes, so pop-up menu
                    displayNodePopupMenu( node, e.getX(), e.getY() );
                    // Handle unadorned click
                }
                else {
                    // Check for right mouse button
                    if ( e.getModifiers() == 4 ) {
                        displayNodePopupMenu( node, e.getX(), e.getY() );
                    }
                    else {
                        // if not in _found_nodes, clear _found_nodes
                        handleClickToAction( _control_panel.getActionWhenNodeClicked(), node );
                    }
                }
            }
            else {
                // no node was clicked
                _highlight_node = null;
            }
        }
        repaint();
    }

    final void mouseDragInBrowserPanel( final MouseEvent e ) {
        setCursor( MOVE_CURSOR );
        final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
        scroll_position.x -= ( e.getX() - getLastDragPointX() );
        scroll_position.y -= ( e.getY() - getLastDragPointY() );
        if ( scroll_position.x < 0 ) {
            scroll_position.x = 0;
        }
        else {
            final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
            if ( scroll_position.x > max_x ) {
                scroll_position.x = max_x;
            }
        }
        if ( scroll_position.y < 0 ) {
            scroll_position.y = 0;
        }
        else {
            final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
            if ( scroll_position.y > max_y ) {
                scroll_position.y = max_y;
            }
        }
        if ( isOvOn() || getOptions().isShowScale() ) {
            repaint();
        }
        getMainPanel().getCurrentScrollPane().getViewport().setViewPosition( scroll_position );
    }

    final void mouseDragInOvRectangle( final MouseEvent e ) {
        setCursor( HAND_CURSOR );
        final double w_ratio = getVisibleRect().width / getOvRectangle().getWidth();
        final double h_ratio = getVisibleRect().height / getOvRectangle().getHeight();
        final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
        double dx = ( ( w_ratio * e.getX() ) - ( w_ratio * getLastDragPointX() ) );
        double dy = ( ( h_ratio * e.getY() ) - ( h_ratio * getLastDragPointY() ) );
        scroll_position.x = ForesterUtil.roundToInt( scroll_position.x + dx );
        scroll_position.y = ForesterUtil.roundToInt( scroll_position.y + dy );
        if ( scroll_position.x <= 0 ) {
            scroll_position.x = 0;
            dx = 0;
        }
        else {
            final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
            if ( scroll_position.x >= max_x ) {
                dx = 0;
                scroll_position.x = max_x;
            }
        }
        if ( scroll_position.y <= 0 ) {
            dy = 0;
            scroll_position.y = 0;
        }
        else {
            final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
            if ( scroll_position.y >= max_y ) {
                dy = 0;
                scroll_position.y = max_y;
            }
        }
        repaint();
        getMainPanel().getCurrentScrollPane().getViewport().setViewPosition( scroll_position );
        setLastMouseDragPointX( ( float ) ( e.getX() + dx ) );
        setLastMouseDragPointY( ( float ) ( e.getY() + dy ) );
    }

    final void mouseMoved( final MouseEvent e ) {
        requestFocusInWindow();
        if ( _current_external_nodes != null ) {
            _current_external_nodes = null;
            repaint();
        }
        if ( getControlPanel().isNodeDescPopup() ) {
            if ( _node_desc_popup != null ) {
                _node_desc_popup.hide();
                _node_desc_popup = null;
            }
        }
        if ( getOptions().isShowOverview() && isOvOn() ) {
            if ( inOvVirtualRectangle( e ) ) {
                if ( !isInOvRect() ) {
                    setInOvRect( true );
                    repaint();
                }
            }
            else {
                if ( isInOvRect() ) {
                    setInOvRect( false );
                    repaint();
                }
            }
        }
        if ( inOv( e ) && getOptions().isShowOverview() && isOvOn() ) {
            if ( !isInOv() ) {
                setInOv( true );
            }
        }
        else {
            if ( isInOv() ) {
                setInOv( false );
            }
            final PhylogenyNode node = findNode( e.getX(), e.getY() );
            if ( ( node != null ) && ( node.isRoot() || !node.getParent().isCollapse() ) ) {
                if ( ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.GET_EXT_DESC_DATA ) ) {
                    for( final PhylogenyNode n : node.getAllExternalDescendants() ) {
                        addToCurrentExternalNodes( n.getId() );
                    }
                    setCursor( HAND_CURSOR );
                    repaint();
                }
                else if ( ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.CUT_SUBTREE )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.COPY_SUBTREE )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.PASTE_SUBTREE )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.DELETE_NODE_OR_SUBTREE )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.REROOT )
                        || ( getControlPanel().getActionWhenNodeClicked() == NodeClickAction.ADD_NEW_NODE ) ) {
                    setCursor( CUT_CURSOR );
                }
                else {
                    setCursor( HAND_CURSOR );
                    if ( getControlPanel().isNodeDescPopup() ) {
                        showNodeDataPopup( e, node );
                    }
                }
            }
            else {
                setCursor( ARROW_CURSOR );
            }
        }
    }

    final void mouseReleasedInBrowserPanel( final MouseEvent e ) {
        setCursor( ARROW_CURSOR );
    }

    final void multiplyUrtFactor( final float f ) {
        _urt_factor *= f;
    }

    final JApplet obtainApplet() {
        return ( ( MainPanelApplets ) getMainPanel() ).getApplet();
    }

    final void paintBranchCircular( final PhylogenyNode p,
                                    final PhylogenyNode c,
                                    final Graphics2D g,
                                    final boolean radial_labels,
                                    final boolean to_pdf,
                                    final boolean to_graphics_file ) {
        final double angle = _urt_nodeid_angle_map.get( c.getId() );
        final double root_x = _root.getXcoord();
        final double root_y = _root.getYcoord();
        final double dx = root_x - p.getXcoord();
        final double dy = root_y - p.getYcoord();
        final double parent_radius = Math.sqrt( ( dx * dx ) + ( dy * dy ) );
        final double arc = ( _urt_nodeid_angle_map.get( p.getId() ) ) - angle;
        assignGraphicsForBranchWithColorForParentBranch( c, false, g, to_pdf, to_graphics_file );
        if ( ( c.isFirstChildNode() || c.isLastChildNode() )
                && ( ( Math.abs( parent_radius * arc ) > 1.5 ) || to_pdf || to_graphics_file ) ) {
            final double r2 = 2.0 * parent_radius;
            drawArc( root_x - parent_radius, root_y - parent_radius, r2, r2, ( -angle - arc ), arc, g );
        }
        drawLine( c.getXcoord(),
                  c.getYcoord(),
                  root_x + ( Math.cos( angle ) * parent_radius ),
                  root_y + ( Math.sin( angle ) * parent_radius ),
                  g );
        paintNodeBox( c.getXcoord(), c.getYcoord(), c, g, to_pdf, to_graphics_file );
        if ( c.isExternal() ) {
            final boolean is_in_found_nodes = isInFoundNodes0( c ) || isInFoundNodes1( c )
                    || isInCurrentExternalNodes( c );
            if ( ( _dynamic_hiding_factor > 1 ) && !is_in_found_nodes
                    && ( ( _urt_nodeid_index_map.get( c.getId() ) % _dynamic_hiding_factor ) != 1 ) ) {
                return;
            }
            paintNodeDataUnrootedCirc( g, c, to_pdf, to_graphics_file, radial_labels, 0, is_in_found_nodes );
        }
    }

    final void paintBranchCircularLite( final PhylogenyNode p, final PhylogenyNode c, final Graphics2D g ) {
        final double angle = _urt_nodeid_angle_map.get( c.getId() );
        final double root_x = _root.getXSecondary();
        final double root_y = _root.getYSecondary();
        final double dx = root_x - p.getXSecondary();
        final double dy = root_y - p.getYSecondary();
        final double arc = ( _urt_nodeid_angle_map.get( p.getId() ) ) - angle;
        final double parent_radius = Math.sqrt( ( dx * dx ) + ( dy * dy ) );
        g.setColor( getTreeColorSet().getOvColor() );
        if ( ( c.isFirstChildNode() || c.isLastChildNode() ) && ( Math.abs( arc ) > 0.02 ) ) {
            final double r2 = 2.0 * parent_radius;
            drawArc( root_x - parent_radius, root_y - parent_radius, r2, r2, ( -angle - arc ), arc, g );
        }
        drawLine( c.getXSecondary(),
                  c.getYSecondary(),
                  root_x + ( Math.cos( angle ) * parent_radius ),
                  root_y + ( Math.sin( angle ) * parent_radius ),
                  g );
        if ( isInFoundNodes( c ) || isInCurrentExternalNodes( c ) ) {
            g.setColor( getColorForFoundNode( c ) );
            drawRectFilled( c.getXSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF, c.getYSecondary()
                    - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF, OVERVIEW_FOUND_NODE_BOX_SIZE, OVERVIEW_FOUND_NODE_BOX_SIZE, g );
        }
    }

    final void paintCircular( final Phylogeny phy,
                              final double starting_angle,
                              final int center_x,
                              final int center_y,
                              final int radius,
                              final Graphics2D g,
                              final boolean to_pdf,
                              final boolean to_graphics_file ) {
        final int circ_num_ext_nodes = phy.getNumberOfExternalNodes() - _collapsed_external_nodeid_set.size();
        System.out.println( "# collapsed external = " + _collapsed_external_nodeid_set.size() );
        _root = phy.getRoot();
        _root.setXcoord( center_x );
        _root.setYcoord( center_y );
        final boolean radial_labels = getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL;
        double current_angle = starting_angle;
        int i = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isCollapse() ) {
                n.setXcoord( ( float ) ( center_x + ( radius * Math.cos( current_angle ) ) ) );
                n.setYcoord( ( float ) ( center_y + ( radius * Math.sin( current_angle ) ) ) );
                _urt_nodeid_angle_map.put( n.getId(), current_angle );
                _urt_nodeid_index_map.put( n.getId(), i++ );
                current_angle += ( TWO_PI / circ_num_ext_nodes );
            }
            else {
                //TODO remove me
                System.out.println( "is collapse" + n.getName() );
            }
        }
        paintCirculars( phy.getRoot(), phy, center_x, center_y, radius, radial_labels, g, to_pdf, to_graphics_file );
        paintNodeBox( _root.getXcoord(), _root.getYcoord(), _root, g, to_pdf, to_graphics_file );
    }

    final void paintCircularLite( final Phylogeny phy,
                                  final double starting_angle,
                                  final int center_x,
                                  final int center_y,
                                  final int radius,
                                  final Graphics2D g ) {
        final int circ_num_ext_nodes = phy.getNumberOfExternalNodes();
        _root = phy.getRoot();
        _root.setXSecondary( center_x );
        _root.setYSecondary( center_y );
        double current_angle = starting_angle;
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.setXSecondary( ( float ) ( center_x + ( radius * Math.cos( current_angle ) ) ) );
            n.setYSecondary( ( float ) ( center_y + ( radius * Math.sin( current_angle ) ) ) );
            _urt_nodeid_angle_map.put( n.getId(), current_angle );
            current_angle += ( TWO_PI / circ_num_ext_nodes );
        }
        paintCircularsLite( phy.getRoot(), phy, center_x, center_y, radius, g );
    }

    final void paintPhylogeny( final Graphics2D g,
                               final boolean to_pdf,
                               final boolean to_graphics_file,
                               final int graphics_file_width,
                               final int graphics_file_height,
                               final int graphics_file_x,
                               final int graphics_file_y ) {
        if ( ( _phylogeny == null ) || _phylogeny.isEmpty() ) {
            return;
        }
        if ( _control_panel.isShowSequenceRelations() ) {
            _query_sequence = _control_panel.getSelectedQuerySequence();
        }
        // Color the background
        if ( !to_pdf ) {
            final Rectangle r = getVisibleRect();
            if ( !getOptions().isBackgroundColorGradient() || getOptions().isPrintBlackAndWhite() ) {
                g.setColor( getTreeColorSet().getBackgroundColor() );
                if ( !to_graphics_file ) {
                    g.fill( r );
                }
                else {
                    if ( getOptions().isPrintBlackAndWhite() ) {
                        g.setColor( Color.WHITE );
                    }
                    g.fillRect( graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height );
                }
            }
            else {
                if ( !to_graphics_file ) {
                    g.setPaint( new GradientPaint( r.x, r.y, getTreeColorSet().getBackgroundColor(), r.x, r.y
                            + r.height, getTreeColorSet().getBackgroundColorGradientBottom() ) );
                    g.fill( r );
                }
                else {
                    g.setPaint( new GradientPaint( graphics_file_x,
                                                   graphics_file_y,
                                                   getTreeColorSet().getBackgroundColor(),
                                                   graphics_file_x,
                                                   graphics_file_y + graphics_file_height,
                                                   getTreeColorSet().getBackgroundColorGradientBottom() ) );
                    g.fillRect( graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height );
                }
            }
            setupStroke( g );
        }
        else {
            g.setStroke( new BasicStroke( getOptions().getPrintLineWidth() ) );
        }
        if ( ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                && ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
            _external_node_index = 0;
            // Position starting X of tree
            if ( !_phylogeny.isRooted() /*|| ( _subtree_index > 0 )*/) {
                _phylogeny.getRoot().setXcoord( TreePanel.MOVE );
            }
            else if ( ( _phylogeny.getRoot().getDistanceToParent() > 0.0 ) && getControlPanel().isDrawPhylogram() ) {
                _phylogeny.getRoot().setXcoord( ( float ) ( TreePanel.MOVE + ( _phylogeny.getRoot()
                        .getDistanceToParent() * getXcorrectionFactor() ) ) );
            }
            else {
                _phylogeny.getRoot().setXcoord( TreePanel.MOVE + getXdistance() );
            }
            // Position starting Y of tree
            _phylogeny.getRoot().setYcoord( ( getYdistance() * _phylogeny.getRoot().getNumberOfExternalNodes() )
                    + ( TreePanel.MOVE / 2.0f ) );
            final int dynamic_hiding_factor = calcDynamicHidingFactor();
            if ( getControlPanel().isDynamicallyHideData() ) {
                if ( dynamic_hiding_factor > 1 ) {
                    getControlPanel().setDynamicHidingIsOn( true );
                }
                else {
                    getControlPanel().setDynamicHidingIsOn( false );
                }
            }
            if ( _nodes_in_preorder == null ) {
                _nodes_in_preorder = new PhylogenyNode[ _phylogeny.getNodeCount() ];
                int i = 0;
                for( final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
                    _nodes_in_preorder[ i++ ] = it.next();
                }
            }
            //final PhylogenyNodeIterator it;
            //for( it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
            //    paintNodeRectangular( g, it.next(), to_pdf, getControlPanel().isDynamicallyHideData()
            //            && ( dynamic_hiding_factor > 1 ), dynamic_hiding_factor, to_graphics_file );
            //}
            for( final PhylogenyNode element : _nodes_in_preorder ) {
                paintNodeRectangular( g, element, to_pdf, getControlPanel().isDynamicallyHideData()
                        && ( dynamic_hiding_factor > 1 ), dynamic_hiding_factor, to_graphics_file );
            }
            if ( getOptions().isShowScale() && getControlPanel().isDrawPhylogram() && ( getScaleDistance() > 0.0 ) ) {
                if ( !( to_graphics_file || to_pdf ) ) {
                    paintScale( g,
                                getVisibleRect().x,
                                getVisibleRect().y + getVisibleRect().height,
                                to_pdf,
                                to_graphics_file );
                }
                else {
                    paintScale( g, graphics_file_x, graphics_file_y + graphics_file_height, to_pdf, to_graphics_file );
                }
            }
            if ( getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf ) {
                paintPhylogenyLite( g );
            }
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            if ( getControlPanel().getDynamicallyHideData() != null ) {
                getControlPanel().setDynamicHidingIsOn( false );
            }
            final double angle = getStartingAngle();
            final boolean radial_labels = getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL;
            _dynamic_hiding_factor = 0;
            if ( getControlPanel().isDynamicallyHideData() ) {
                _dynamic_hiding_factor = ( int ) ( ( getFontMetricsForLargeDefaultFont().getHeight() * 1.5 * getPhylogeny()
                        .getNumberOfExternalNodes() ) / ( TWO_PI * 10 ) );
            }
            if ( getControlPanel().getDynamicallyHideData() != null ) {
                if ( _dynamic_hiding_factor > 1 ) {
                    getControlPanel().setDynamicHidingIsOn( true );
                }
                else {
                    getControlPanel().setDynamicHidingIsOn( false );
                }
            }
            paintUnrooted( _phylogeny.getRoot(),
                           angle,
                           ( float ) ( angle + ( 2 * Math.PI ) ),
                           radial_labels,
                           g,
                           to_pdf,
                           to_graphics_file );
            if ( getOptions().isShowScale() ) {
                if ( !( to_graphics_file || to_pdf ) ) {
                    paintScale( g,
                                getVisibleRect().x,
                                getVisibleRect().y + getVisibleRect().height,
                                to_pdf,
                                to_graphics_file );
                }
                else {
                    paintScale( g, graphics_file_x, graphics_file_y + graphics_file_height, to_pdf, to_graphics_file );
                }
            }
            if ( getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf ) {
                g.setColor( getTreeColorSet().getOvColor() );
                paintUnrootedLite( _phylogeny.getRoot(),
                                   angle,
                                   angle + ( 2 * Math.PI ),
                                   g,
                                   ( getUrtFactorOv() / ( getVisibleRect().width / getOvMaxWidth() ) ) );
                paintOvRectangle( g );
            }
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            final int radius = ( int ) ( ( Math.min( getPreferredSize().getWidth(), getPreferredSize().getHeight() ) / 2 ) - ( MOVE + getLongestExtNodeInfo() ) );
            final int d = radius + MOVE + getLongestExtNodeInfo();
            _dynamic_hiding_factor = 0;
            if ( getControlPanel().isDynamicallyHideData() && ( radius > 0 ) ) {
                _dynamic_hiding_factor = ( int ) ( ( getFontMetricsForLargeDefaultFont().getHeight() * 1.5 * getPhylogeny()
                        .getNumberOfExternalNodes() ) / ( TWO_PI * radius ) );
            }
            if ( getControlPanel().getDynamicallyHideData() != null ) {
                if ( _dynamic_hiding_factor > 1 ) {
                    getControlPanel().setDynamicHidingIsOn( true );
                }
                else {
                    getControlPanel().setDynamicHidingIsOn( false );
                }
            }
            paintCircular( _phylogeny, getStartingAngle(), d, d, radius > 0 ? radius : 0, g, to_pdf, to_graphics_file );
            if ( getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf ) {
                final int radius_ov = ( int ) ( getOvMaxHeight() < getOvMaxWidth() ? getOvMaxHeight() / 2
                        : getOvMaxWidth() / 2 );
                double x_scale = 1.0;
                double y_scale = 1.0;
                int x_pos = getVisibleRect().x + getOvXPosition();
                int y_pos = getVisibleRect().y + getOvYPosition();
                if ( getWidth() > getHeight() ) {
                    x_scale = ( double ) getHeight() / getWidth();
                    x_pos = ForesterUtil.roundToInt( x_pos / x_scale );
                }
                else {
                    y_scale = ( double ) getWidth() / getHeight();
                    y_pos = ForesterUtil.roundToInt( y_pos / y_scale );
                }
                _at = g.getTransform();
                g.scale( x_scale, y_scale );
                paintCircularLite( _phylogeny,
                                   getStartingAngle(),
                                   x_pos + radius_ov,
                                   y_pos + radius_ov,
                                   ( int ) ( radius_ov - ( getLongestExtNodeInfo() / ( getVisibleRect().width / getOvRectangle()
                                           .getWidth() ) ) ),
                                   g );
                g.setTransform( _at );
                paintOvRectangle( g );
            }
        }
    }

    final void recalculateMaxDistanceToRoot() {
        _max_distance_to_root = PhylogenyMethods.calculateMaxDistanceToRoot( getPhylogeny() );
    }

    /**
     * Remove all edit-node frames
     */
    final void removeAllEditNodeJFrames() {
        for( int i = 0; i <= ( TreePanel.MAX_NODE_FRAMES - 1 ); i++ ) {
            if ( _node_frames[ i ] != null ) {
                _node_frames[ i ].dispose();
                _node_frames[ i ] = null;
            }
        }
        _node_frame_index = 0;
    }

    /**
     * Remove a node-edit frame.
     */
    final void removeEditNodeFrame( final int i ) {
        _node_frame_index--;
        _node_frames[ i ] = null;
        if ( i < _node_frame_index ) {
            for( int j = 0; j < ( _node_frame_index - 1 ); j++ ) {
                _node_frames[ j ] = _node_frames[ j + 1 ];
            }
            _node_frames[ _node_frame_index ] = null;
        }
    }

    final void reRoot( final PhylogenyNode node ) {
        if ( !getPhylogeny().isRerootable() ) {
            JOptionPane.showMessageDialog( this,
                                           "This is not rerootable",
                                           "Not rerootable",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot reroot in unrooted display type",
                                           "Attempt to reroot tree in unrooted display",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        getPhylogeny().reRoot( node );
        getPhylogeny().recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setNodeInPreorderToNull();
        resetPreferredSize();
        getMainPanel().adjustJScrollPane();
        setEdited( true );
        repaint();
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            getControlPanel().showWhole();
        }
    }

    final void resetNodeIdToDistToLeafMap() {
        _nodeid_dist_to_leaf = new HashMap<Long, Short>();
    }

    final void resetPreferredSize() {
        if ( ( getPhylogeny() == null ) || getPhylogeny().isEmpty() ) {
            return;
        }
        int x = 0;
        int y = 0;
        y = TreePanel.MOVE
                + ForesterUtil.roundToInt( getYdistance() * getPhylogeny().getRoot().getNumberOfExternalNodes() * 2 );
        if ( getControlPanel().isDrawPhylogram() ) {
            x = TreePanel.MOVE
                    + getLongestExtNodeInfo()
                    + ForesterUtil
                            .roundToInt( ( getXcorrectionFactor() * getPhylogeny().getHeight() ) + getXdistance() );
        }
        else {
            if ( !isNonLinedUpCladogram() && !isUniformBranchLengthsForCladogram() ) {
                x = TreePanel.MOVE
                        + getLongestExtNodeInfo()
                        + ForesterUtil.roundToInt( getXdistance()
                                * ( getPhylogeny().getRoot().getNumberOfExternalNodes() + 2 ) );
            }
            else {
                x = TreePanel.MOVE
                        + getLongestExtNodeInfo()
                        + ForesterUtil.roundToInt( getXdistance()
                                * ( PhylogenyMethods.calculateMaxDepth( getPhylogeny() ) + 1 ) );
            }
        }
        setPreferredSize( new Dimension( x, y ) );
    }

    final void selectNode( final PhylogenyNode node ) {
        if ( ( getFoundNodes0() != null ) && getFoundNodes0().contains( node.getId() ) ) {
            getFoundNodes0().remove( node.getId() );
            getControlPanel().setSearchFoundCountsOnLabel0( getFoundNodes0().size() );
            if ( getFoundNodes0().size() < 1 ) {
                getControlPanel().searchReset0();
            }
        }
        else {
            getControlPanel().getSearchFoundCountsLabel0().setVisible( true );
            getControlPanel().getSearchResetButton0().setEnabled( true );
            getControlPanel().getSearchResetButton0().setVisible( true );
            if ( getFoundNodes0() == null ) {
                setFoundNodes0( new HashSet<Long>() );
            }
            getFoundNodes0().add( node.getId() );
            getControlPanel().setSearchFoundCountsOnLabel0( getFoundNodes0().size() );
        }
    }

    final void setArrowCursor() {
        setCursor( ARROW_CURSOR );
        repaint();
    }

    final void setControlPanel( final ControlPanel atv_control ) {
        _control_panel = atv_control;
    }

    void setCurrentExternalNodesDataBuffer( final StringBuilder sb ) {
        increaseCurrentExternalNodesDataBufferChangeCounter();
        _current_external_nodes_data_buffer = sb;
    }

    final void setFoundNodes0( final Set<Long> found_nodes ) {
        _found_nodes_0 = found_nodes;
    }

    final void setFoundNodes1( final Set<Long> found_nodes ) {
        _found_nodes_1 = found_nodes;
    }

    final void setInOvRect( final boolean in_ov_rect ) {
        _in_ov_rect = in_ov_rect;
    }

    final void setLargeFonts() {
        getTreeFontSet().largeFonts();
    }

    final void setLastMouseDragPointX( final float x ) {
        _last_drag_point_x = x;
    }

    final void setLastMouseDragPointY( final float y ) {
        _last_drag_point_y = y;
    }

    final void setLongestExtNodeInfo( final int i ) {
        _longest_ext_node_info = i;
    }

    final void setMediumFonts() {
        getTreeFontSet().mediumFonts();
    }

    final void setNodeInPreorderToNull() {
        _nodes_in_preorder = null;
    }

    final void setOvOn( final boolean ov_on ) {
        _ov_on = ov_on;
    }

    final void setPhylogenyGraphicsType( final PHYLOGENY_GRAPHICS_TYPE graphics_type ) {
        _graphics_type = graphics_type;
        setTextAntialias();
    }

    final void setSmallFonts() {
        getTreeFontSet().smallFonts();
    }

    final void setStartingAngle( final double starting_angle ) {
        _urt_starting_angle = starting_angle;
    }

    void setStatisticsForExpressionValues( final DescriptiveStatistics statistics_for_expression_values ) {
        _statistics_for_vector_data = statistics_for_expression_values;
    }

    final void setSuperTinyFonts() {
        getTreeFontSet().superTinyFonts();
    }

    final void setTextAntialias() {
        if ( ( _phylogeny != null ) && !_phylogeny.isEmpty() ) {
            if ( _phylogeny.getNumberOfExternalNodes() <= LIMIT_FOR_HQ_RENDERING ) {
                _rendering_hints.put( RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY );
            }
            else {
                _rendering_hints.put( RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED );
            }
        }
        if ( getMainPanel().getOptions().isAntialiasScreen() ) {
            _rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
            // try {
            _rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB );
            // }
            // catch ( final Throwable e ) {
            //    _rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
            //}
        }
        else {
            _rendering_hints.put( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF );
            _rendering_hints.put( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF );
        }
    }

    final void setTinyFonts() {
        getTreeFontSet().tinyFonts();
    }

    final void setTreeFile( final File treefile ) {
        _treefile = treefile;
    }

    final void setXcorrectionFactor( final float f ) {
        _x_correction_factor = f;
    }

    final void setXdistance( final float x ) {
        _x_distance = x;
    }

    final void setYdistance( final float y ) {
        _y_distance = y;
    }

    final void sortDescendants( final PhylogenyNode node ) {
        if ( !node.isExternal() ) {
            DESCENDANT_SORT_PRIORITY pri = DESCENDANT_SORT_PRIORITY.TAXONOMY;
            if ( ( !getControlPanel().isShowTaxonomyScientificNames() && !getControlPanel().isShowTaxonomyCode() && !getControlPanel()
                    .isShowTaxonomyCommonNames() ) ) {
                if ( ( getControlPanel().isShowSequenceAcc() || getControlPanel().isShowSeqNames() || getControlPanel()
                        .isShowSeqSymbols() ) ) {
                    pri = DESCENDANT_SORT_PRIORITY.SEQUENCE;
                }
                else if ( getControlPanel().isShowNodeNames() ) {
                    pri = DESCENDANT_SORT_PRIORITY.NODE_NAME;
                }
            }
            PhylogenyMethods.sortNodeDescendents( node, pri );
            setNodeInPreorderToNull();
            _phylogeny.externalNodesHaveChanged();
            _phylogeny.clearHashIdToNodeMap();
            _phylogeny.recalculateNumberOfExternalDescendants( true );
            resetNodeIdToDistToLeafMap();
            setEdited( true );
        }
        repaint();
    }

    final void subTree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot get a sub/super tree in unrooted display",
                                           "Attempt to get sub/super tree in unrooted display",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( node.isExternal() ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot get a subtree of a external node",
                                           "Attempt to get subtree of external node",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        if ( node.isRoot() && !isCurrentTreeIsSubtree() ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot get a subtree of the root node",
                                           "Attempt to get subtree of root node",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        setNodeInPreorderToNull();
        if ( !node.isExternal() && !node.isRoot() && ( _subtree_index <= ( TreePanel.MAX_SUBTREES - 1 ) ) ) {
            _sub_phylogenies[ _subtree_index ] = _phylogeny;
            _sub_phylogenies_temp_roots[ _subtree_index ] = node;
            ++_subtree_index;
            _phylogeny = TreePanelUtil.subTree( node, _phylogeny );
            updateSubSuperTreeButton();
        }
        else if ( node.isRoot() && isCurrentTreeIsSubtree() ) {
            superTree();
        }
        _main_panel.getControlPanel().showWhole();
        repaint();
    }

    final void superTree() {
        setNodeInPreorderToNull();
        final PhylogenyNode temp_root = _sub_phylogenies_temp_roots[ _subtree_index - 1 ];
        for( final PhylogenyNode n : temp_root.getDescendants() ) {
            n.setParent( temp_root );
        }
        _sub_phylogenies[ _subtree_index ] = null;
        _sub_phylogenies_temp_roots[ _subtree_index ] = null;
        _phylogeny = _sub_phylogenies[ --_subtree_index ];
        updateSubSuperTreeButton();
    }

    final void swap( final PhylogenyNode node ) {
        if ( node.isExternal() || ( node.getNumberOfDescendants() < 2 ) ) {
            return;
        }
        if ( node.getNumberOfDescendants() > 2 ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot swap descendants of nodes with more than 2 descendants",
                                           "Cannot swap descendants",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        if ( !node.isExternal() ) {
            node.swapChildren();
            setNodeInPreorderToNull();
            _phylogeny.externalNodesHaveChanged();
            _phylogeny.clearHashIdToNodeMap();
            _phylogeny.recalculateNumberOfExternalDescendants( true );
            resetNodeIdToDistToLeafMap();
            setEdited( true );
        }
        repaint();
    }

    final void taxColor() {
        if ( ( _phylogeny == null ) || ( _phylogeny.getNumberOfExternalNodes() < 2 ) ) {
            return;
        }
        setWaitCursor();
        TreePanelUtil.colorPhylogenyAccordingToExternalTaxonomy( _phylogeny, this );
        _control_panel.setColorBranches( true );
        if ( _control_panel.getUseVisualStylesCb() != null ) {
            _control_panel.getUseVisualStylesCb().setSelected( true );
        }
        setEdited( true );
        setArrowCursor();
        repaint();
    }

    final void updateOvSettings() {
        switch ( getOptions().getOvPlacement() ) {
            case LOWER_LEFT:
                setOvXPosition( OV_BORDER );
                setOvYPosition( ForesterUtil.roundToInt( getVisibleRect().height - OV_BORDER - getOvMaxHeight() ) );
                setOvYStart( ForesterUtil.roundToInt( getOvYPosition() + ( getOvMaxHeight() / 2 ) ) );
                break;
            case LOWER_RIGHT:
                setOvXPosition( ForesterUtil.roundToInt( getVisibleRect().width - OV_BORDER - getOvMaxWidth() ) );
                setOvYPosition( ForesterUtil.roundToInt( getVisibleRect().height - OV_BORDER - getOvMaxHeight() ) );
                setOvYStart( ForesterUtil.roundToInt( getOvYPosition() + ( getOvMaxHeight() / 2 ) ) );
                break;
            case UPPER_RIGHT:
                setOvXPosition( ForesterUtil.roundToInt( getVisibleRect().width - OV_BORDER - getOvMaxWidth() ) );
                setOvYPosition( OV_BORDER );
                setOvYStart( ForesterUtil.roundToInt( OV_BORDER + ( getOvMaxHeight() / 2 ) ) );
                break;
            default:
                setOvXPosition( OV_BORDER );
                setOvYPosition( OV_BORDER );
                setOvYStart( ForesterUtil.roundToInt( OV_BORDER + ( getOvMaxHeight() / 2 ) ) );
                break;
        }
    }

    final void updateOvSizes() {
        if ( ( getWidth() > ( 1.05 * getVisibleRect().width ) ) || ( getHeight() > ( 1.05 * getVisibleRect().height ) ) ) {
            setOvOn( true );
            float l = getLongestExtNodeInfo();
            final float w_ratio = getOvMaxWidth() / getWidth();
            l *= w_ratio;
            final int ext_nodes = _phylogeny.getRoot().getNumberOfExternalNodes();
            setOvYDistance( getOvMaxHeight() / ( 2 * ext_nodes ) );
            float ov_xdist = 0;
            if ( !isNonLinedUpCladogram() && !isUniformBranchLengthsForCladogram() ) {
                ov_xdist = ( ( getOvMaxWidth() - l ) / ( ext_nodes ) );
            }
            else {
                ov_xdist = ( ( getOvMaxWidth() - l ) / ( PhylogenyMethods.calculateMaxDepth( _phylogeny ) ) );
            }
            float ydist = ( float ) ( ( getOvMaxWidth() / ( ext_nodes * 2.0 ) ) );
            if ( ov_xdist < 0.0 ) {
                ov_xdist = 0.0f;
            }
            if ( ydist < 0.0 ) {
                ydist = 0.0f;
            }
            setOvXDistance( ov_xdist );
            final double height = _phylogeny.getHeight();
            if ( height > 0 ) {
                final float ov_corr = ( float ) ( ( ( getOvMaxWidth() - l ) - getOvXDistance() ) / height );
                setOvXcorrectionFactor( ov_corr > 0 ? ov_corr : 0 );
            }
            else {
                setOvXcorrectionFactor( 0 );
            }
        }
        else {
            setOvOn( false );
        }
    }

    void updateSetOfCollapsedExternalNodes() {
        final Phylogeny phy = getPhylogeny();
        _collapsed_external_nodeid_set.clear();
        if ( phy != null ) {
            E: for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
                final PhylogenyNode ext_node = it.next();
                PhylogenyNode n = ext_node;
                while ( !n.isRoot() ) {
                    if ( n.isCollapse() ) {
                        _collapsed_external_nodeid_set.add( ext_node.getId() );
                        ext_node.setCollapse( true );
                        continue E;
                    }
                    n = n.getParent();
                }
            }
        }
    }

    final void updateSubSuperTreeButton() {
        if ( _subtree_index < 1 ) {
            getControlPanel().deactivateButtonToReturnToSuperTree();
        }
        else {
            getControlPanel().activateButtonToReturnToSuperTree( _subtree_index );
        }
    }

    final void zoomInDomainStructure() {
        if ( _domain_structure_width < 2000 ) {
            _domain_structure_width *= 1.2;
        }
    }

    final void zoomOutDomainStructure() {
        if ( _domain_structure_width > 20 ) {
            _domain_structure_width *= 0.8;
        }
    }

    private void abbreviateScientificName( final String sn, final StringBuilder sb ) {
        final String[] a = sn.split( "\\s+" );
        sb.append( a[ 0 ].substring( 0, 1 ) );
        sb.append( a[ 1 ].substring( 0, 2 ) );
        if ( a.length > 2 ) {
            for( int i = 2; i < a.length; i++ ) {
                sb.append( " " );
                sb.append( a[ i ] );
            }
        }
    }

    final private void addEmptyNode( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        final String label = createASimpleTextRepresentationOfANode( node );
        String msg = "";
        if ( ForesterUtil.isEmpty( label ) ) {
            msg = "How to add the new, empty node?";
        }
        else {
            msg = "How to add the new, empty node to node" + label + "?";
        }
        final Object[] options = { "As sibling", "As descendant", "Cancel" };
        final int r = JOptionPane.showOptionDialog( this,
                                                    msg,
                                                    "Addition of Empty New Node",
                                                    JOptionPane.CLOSED_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE,
                                                    null,
                                                    options,
                                                    options[ 2 ] );
        boolean add_as_sibling = true;
        if ( r == 1 ) {
            add_as_sibling = false;
        }
        else if ( r != 0 ) {
            return;
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( new PhylogenyNode() );
        phy.setRooted( true );
        if ( add_as_sibling ) {
            if ( node.isRoot() ) {
                JOptionPane.showMessageDialog( this,
                                               "Cannot add sibling to root",
                                               "Attempt to add sibling to root",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            phy.addAsSibling( node );
        }
        else {
            phy.addAsChild( node );
        }
        setNodeInPreorderToNull();
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setEdited( true );
        repaint();
    }

    final private void addToCurrentExternalNodes( final long i ) {
        if ( _current_external_nodes == null ) {
            _current_external_nodes = new HashSet<Long>();
        }
        _current_external_nodes.add( i );
    }

    final private void assignGraphicsForBranchWithColorForParentBranch( final PhylogenyNode node,
                                                                        final boolean is_vertical,
                                                                        final Graphics g,
                                                                        final boolean to_pdf,
                                                                        final boolean to_graphics_file ) {
        final NodeClickAction action = _control_panel.getActionWhenNodeClicked();
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else if ( ( ( action == NodeClickAction.COPY_SUBTREE ) || ( action == NodeClickAction.CUT_SUBTREE )
                || ( action == NodeClickAction.DELETE_NODE_OR_SUBTREE ) || ( action == NodeClickAction.PASTE_SUBTREE ) || ( action == NodeClickAction.ADD_NEW_NODE ) )
                && ( getCutOrCopiedTree() != null )
                && ( getCopiedAndPastedNodes() != null )
                && !to_pdf
                && !to_graphics_file && getCopiedAndPastedNodes().contains( node.getId() ) ) {
            g.setColor( getTreeColorSet().getFoundColor0() );
        }
        else if ( getControlPanel().isUseVisualStyles() && ( PhylogenyMethods.getBranchColorValue( node ) != null ) ) {
            g.setColor( PhylogenyMethods.getBranchColorValue( node ) );
        }
        else if ( to_pdf ) {
            g.setColor( getTreeColorSet().getBranchColorForPdf() );
        }
        else {
            g.setColor( getTreeColorSet().getBranchColor() );
        }
    }

    final private void blast( final PhylogenyNode node ) {
        if ( !isCanBlast( node ) ) {
            JOptionPane.showMessageDialog( this,
                                           "Insufficient information present",
                                           "Cannot Blast",
                                           JOptionPane.INFORMATION_MESSAGE );
            return;
        }
        else {
            final String query = Blast.obtainQueryForBlast( node );
            System.out.println( "query for BLAST is: " + query );
            char type = '?';
            if ( !ForesterUtil.isEmpty( query ) ) {
                if ( node.getNodeData().isHasSequence() ) {
                    if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getType() ) ) {
                        if ( node.getNodeData().getSequence().getType().toLowerCase()
                                .equals( PhyloXmlUtil.SEQ_TYPE_PROTEIN ) ) {
                            type = 'p';
                        }
                        else {
                            type = 'n';
                        }
                    }
                    else if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) {
                        if ( ForesterUtil.seqIsLikelyToBeAa( node.getNodeData().getSequence().getMolecularSequence() ) ) {
                            type = 'p';
                        }
                        else {
                            type = 'n';
                        }
                    }
                }
                if ( type == '?' ) {
                    if ( SequenceAccessionTools.isProteinDbQuery( query ) ) {
                        type = 'p';
                    }
                    else {
                        type = 'n';
                    }
                }
                JApplet applet = null;
                if ( isApplet() ) {
                    applet = obtainApplet();
                }
                try {
                    Blast.openNcbiBlastWeb( query, type == 'n', applet, this );
                }
                catch ( final Exception e ) {
                    e.printStackTrace();
                }
                if ( Constants.ALLOW_DDBJ_BLAST ) {
                    try {
                        System.out.println( "trying: " + query );
                        final Blast s = new Blast();
                        s.ddbjBlast( query );
                    }
                    catch ( final Exception e ) {
                        e.printStackTrace();
                    }
                }
            }
        }
    }

    private final int calcDynamicHidingFactor() {
        return ( int ) ( 0.5 + ( getFontMetricsForLargeDefaultFont().getHeight() / ( 1.5 * getYdistance() ) ) );
    }

    /**
     * Calculate the length of the distance between the given node and its
     * parent.
     * 
     * @param node
     * @param ext_node_x
     * @factor
     * @return the distance value
     */
    final private float calculateBranchLengthToParent( final PhylogenyNode node, final float factor ) {
        if ( getControlPanel().isDrawPhylogram() ) {
            if ( node.getDistanceToParent() < 0.0 ) {
                return 0.0f;
            }
            return ( float ) ( getXcorrectionFactor() * node.getDistanceToParent() );
        }
        else {
            if ( ( factor == 0 ) || isNonLinedUpCladogram() ) {
                return getXdistance();
            }
            return getXdistance() * factor;
        }
    }

    final private Color calculateColorForAnnotation( final SortedSet<Annotation> ann ) {
        Color c = getTreeColorSet().getAnnotationColor();
        if ( getControlPanel().isColorAccordingToAnnotation() && ( getControlPanel().getAnnotationColors() != null ) ) {
            final StringBuilder sb = new StringBuilder();
            for( final Annotation a : ann ) {
                sb.append( !ForesterUtil.isEmpty( a.getRefValue() ) ? a.getRefValue() : a.getDesc() );
            }
            final String ann_str = sb.toString();
            if ( !ForesterUtil.isEmpty( ann_str ) ) {
                c = getControlPanel().getAnnotationColors().get( ann_str );
                if ( c == null ) {
                    c = TreePanelUtil.calculateColorFromString( ann_str, false );
                    getControlPanel().getAnnotationColors().put( ann_str, c );
                }
                if ( c == null ) {
                    c = getTreeColorSet().getAnnotationColor();
                }
            }
        }
        
        return c;
    }

    final private float calculateOvBranchLengthToParent( final PhylogenyNode node, final int factor ) {
        if ( getControlPanel().isDrawPhylogram() ) {
            if ( node.getDistanceToParent() < 0.0 ) {
                return 0.0f;
            }
            return ( float ) ( getOvXcorrectionFactor() * node.getDistanceToParent() );
        }
        else {
            if ( ( factor == 0 ) || isNonLinedUpCladogram() ) {
                return getOvXDistance();
            }
            return getOvXDistance() * factor;
        }
    }

    final private void cannotOpenBrowserWarningMessage( final String type_type ) {
        JOptionPane.showMessageDialog( this,
                                       "Cannot launch web browser for " + type_type + " data of this node",
                                       "Cannot launch web browser",
                                       JOptionPane.WARNING_MESSAGE );
    }

    private void changeNodeFont( final PhylogenyNode node ) {
        final FontChooser fc = new FontChooser();
        Font f = null;
        if ( ( node.getNodeData().getNodeVisualData() != null ) && !node.getNodeData().getNodeVisualData().isEmpty() ) {
            f = node.getNodeData().getNodeVisualData().getFont();
        }
        if ( f != null ) {
            fc.setFont( f );
        }
        else {
            fc.setFont( getMainPanel().getTreeFontSet().getLargeFont() );
        }
        fc.showDialog( this, "Select Font" );
        if ( ( fc.getFont() != null ) && !ForesterUtil.isEmpty( fc.getFont().getFamily().trim() ) ) {
            List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
            if ( ( getFoundNodes0() != null ) || ( getFoundNodes1() != null ) ) {
                nodes = getFoundNodesAsListOfPhylogenyNodes();
            }
            nodes.add( node );
            for( final PhylogenyNode n : nodes ) {
                if ( n.getNodeData().getNodeVisualData() == null ) {
                    n.getNodeData().setNodeVisualData( new NodeVisualData() );
                }
                final NodeVisualData vd = n.getNodeData().getNodeVisualData();
                final Font ff = fc.getFont();
                vd.setFontName( ff.getFamily().trim() );
                int s = ff.getSize();
                if ( s < 0 ) {
                    s = 0;
                }
                if ( s > Byte.MAX_VALUE ) {
                    s = Byte.MAX_VALUE;
                }
                vd.setFontSize( s );
                vd.setFontStyle( ff.getStyle() );
            }
            if ( _control_panel.getUseVisualStylesCb() != null ) {
                getControlPanel().getUseVisualStylesCb().setSelected( true );
            }
        }
        repaint();
    }

    final private void colorizeNodes( final Color c,
                                      final PhylogenyNode node,
                                      final List<PhylogenyNode> additional_nodes ) {
        _control_panel.setColorBranches( true );
        if ( _control_panel.getUseVisualStylesCb() != null ) {
            _control_panel.getUseVisualStylesCb().setSelected( true );
        }
        if ( node != null ) {
            colorizeNodesHelper( c, node );
        }
        if ( additional_nodes != null ) {
            for( final PhylogenyNode n : additional_nodes ) {
                colorizeNodesHelper( c, n );
            }
        }
        repaint();
    }

    final private void colorizeSubtree( final Color c,
                                        final PhylogenyNode node,
                                        final List<PhylogenyNode> additional_nodes ) {
        _control_panel.setColorBranches( true );
        if ( _control_panel.getUseVisualStylesCb() != null ) {
            _control_panel.getUseVisualStylesCb().setSelected( true );
        }
        if ( node != null ) {
            for( final PreorderTreeIterator it = new PreorderTreeIterator( node ); it.hasNext(); ) {
                it.next().getBranchData().setBranchColor( new BranchColor( c ) );
            }
        }
        if ( additional_nodes != null ) {
            for( final PhylogenyNode n : additional_nodes ) {
                n.getBranchData().setBranchColor( new BranchColor( c ) );
            }
        }
        repaint();
    }

    private void colorNodeFont( final PhylogenyNode node ) {
        _color_chooser.setPreviewPanel( new JPanel() );
        NodeColorizationActionListener al;
        if ( ( getFoundNodes0() != null ) || ( getFoundNodes1() != null ) ) {
            final List<PhylogenyNode> additional_nodes = getFoundNodesAsListOfPhylogenyNodes();
            al = new NodeColorizationActionListener( _color_chooser, node, additional_nodes );
        }
        else {
            al = new NodeColorizationActionListener( _color_chooser, node );
        }
        final JDialog dialog = JColorChooser.createDialog( this, "Node colorization", true, _color_chooser, al, null );
        dialog.setVisible( true );
    }

    final private void colorSubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot colorize subtree in unrooted display type",
                                           "Attempt to colorize subtree in unrooted display",
                                           JOptionPane.WARNING_MESSAGE );
            return;
        }
        _color_chooser.setPreviewPanel( new JPanel() );
        SubtreeColorizationActionListener al;
        if ( ( getFoundNodes0() != null ) || ( getFoundNodes1() != null ) ) {
            final List<PhylogenyNode> additional_nodes = getFoundNodesAsListOfPhylogenyNodes();
            al = new SubtreeColorizationActionListener( _color_chooser, node, additional_nodes );
        }
        else {
            al = new SubtreeColorizationActionListener( _color_chooser, node );
        }
        final JDialog dialog = JColorChooser
                .createDialog( this, "Subtree colorization", true, _color_chooser, al, null );
        dialog.setVisible( true );
    }

    final private void copySubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        setNodeInPreorderToNull();
        setCutOrCopiedTree( _phylogeny.copy( node ) );
        final List<PhylogenyNode> nodes = PhylogenyMethods.getAllDescendants( node );
        final Set<Long> node_ids = new HashSet<Long>( nodes.size() );
        for( final PhylogenyNode n : nodes ) {
            node_ids.add( n.getId() );
        }
        node_ids.add( node.getId() );
        setCopiedAndPastedNodes( node_ids );
        repaint();
    }

    final private String createASimpleTextRepresentationOfANode( final PhylogenyNode node ) {
        final String tax = PhylogenyMethods.getSpecies( node );
        String label = node.getName();
        if ( !ForesterUtil.isEmpty( label ) && !ForesterUtil.isEmpty( tax ) ) {
            label = label + " " + tax;
        }
        else if ( !ForesterUtil.isEmpty( tax ) ) {
            label = tax;
        }
        else {
            label = "";
        }
        if ( !ForesterUtil.isEmpty( label ) ) {
            label = " [" + label + "]";
        }
        return label;
    }

    final private void cutSubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if ( node.isRoot() ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot cut entire tree as subtree",
                                           "Attempt to cut entire tree",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final String label = createASimpleTextRepresentationOfANode( node );
        final int r = JOptionPane.showConfirmDialog( null,
                                                     "Cut subtree" + label + "?",
                                                     "Confirm Cutting of Subtree",
                                                     JOptionPane.YES_NO_OPTION );
        if ( r != JOptionPane.OK_OPTION ) {
            return;
        }
        setNodeInPreorderToNull();
        setCopiedAndPastedNodes( null );
        setCutOrCopiedTree( _phylogeny.copy( node ) );
        _phylogeny.deleteSubtree( node, true );
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setEdited( true );
        repaint();
    }

    final private void cycleColors() {
        getMainPanel().getTreeColorSet().cycleColorScheme();
        for( final TreePanel tree_panel : getMainPanel().getTreePanels() ) {
            tree_panel.setBackground( getMainPanel().getTreeColorSet().getBackgroundColor() );
        }
    }

    final private void decreaseOvSize() {
        if ( ( getOvMaxWidth() > 20 ) && ( getOvMaxHeight() > 20 ) ) {
            setOvMaxWidth( getOvMaxWidth() - 5 );
            setOvMaxHeight( getOvMaxHeight() - 5 );
            updateOvSettings();
            getControlPanel().displayedPhylogenyMightHaveChanged( false );
        }
    }

    final private void deleteNodeOrSubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if ( node.isRoot() && ( node.getNumberOfDescendants() != 1 ) ) {
            JOptionPane.showMessageDialog( this,
                                           "Cannot delete entire tree",
                                           "Attempt to delete entire tree",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final String label = createASimpleTextRepresentationOfANode( node );
        final Object[] options = { "Node only", "Entire subtree", "Cancel" };
        final int r = JOptionPane.showOptionDialog( this,
                                                    "Delete" + label + "?",
                                                    "Delete Node/Subtree",
                                                    JOptionPane.CLOSED_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE,
                                                    null,
                                                    options,
                                                    options[ 2 ] );
        setNodeInPreorderToNull();
        boolean node_only = true;
        if ( r == 1 ) {
            node_only = false;
        }
        else if ( r != 0 ) {
            return;
        }
        if ( node_only ) {
            PhylogenyMethods.removeNode( node, _phylogeny );
        }
        else {
            _phylogeny.deleteSubtree( node, true );
        }
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setEdited( true );
        repaint();
    }

    final private void displayNodePopupMenu( final PhylogenyNode node, final int x, final int y ) {
        makePopupMenus( node );
        _node_popup_menu.putClientProperty( NODE_POPMENU_NODE_CLIENT_PROPERTY, node );
        _node_popup_menu.show( this, x, y );
    }

    final private void drawArc( final double x,
                                final double y,
                                final double width,
                                final double heigth,
                                final double start_angle,
                                final double arc_angle,
                                final Graphics2D g ) {
        _arc.setArc( x, y, width, heigth, _180_OVER_PI * start_angle, _180_OVER_PI * arc_angle, Arc2D.OPEN );
        g.draw( _arc );
    }

    final private void drawLine( final double x1, final double y1, final double x2, final double y2, final Graphics2D g ) {
        if ( ( x1 == x2 ) && ( y1 == y2 ) ) {
            return;
        }
        _line.setLine( x1, y1, x2, y2 );
        g.draw( _line );
    }

    final private void drawOval( final double x,
                                 final double y,
                                 final double width,
                                 final double heigth,
                                 final Graphics2D g ) {
        _ellipse.setFrame( x, y, width, heigth );
        g.draw( _ellipse );
    }

    final private void drawOvalFilled( final double x,
                                       final double y,
                                       final double width,
                                       final double heigth,
                                       final Graphics2D g ) {
        _ellipse.setFrame( x, y, width, heigth );
        g.fill( _ellipse );
    }

    final private void drawOvalGradient( final float x,
                                         final float y,
                                         final float width,
                                         final float heigth,
                                         final Graphics2D g,
                                         final Color color_1,
                                         final Color color_2,
                                         final Color color_border ) {
        _ellipse.setFrame( x, y, width, heigth );
        g.setPaint( new GradientPaint( x, y, color_1, ( x + width ), ( y + heigth ), color_2, false ) );
        g.fill( _ellipse );
        if ( color_border != null ) {
            g.setPaint( color_border );
            g.draw( _ellipse );
        }
    }

    final private void drawRect( final float x, final float y, final float width, final float heigth, final Graphics2D g ) {
        _rectangle.setFrame( x, y, width, heigth );
        g.draw( _rectangle );
    }

    final private void drawRectFilled( final double x,
                                       final double y,
                                       final double width,
                                       final double heigth,
                                       final Graphics2D g ) {
        _rectangle.setFrame( x, y, width, heigth );
        g.fill( _rectangle );
    }

    final private void drawRectGradient( final float x,
                                         final float y,
                                         final float width,
                                         final float heigth,
                                         final Graphics2D g,
                                         final Color color_1,
                                         final Color color_2,
                                         final Color color_border ) {
        _rectangle.setFrame( x, y, width, heigth );
        g.setPaint( new GradientPaint( x, y, color_1, ( x + width ), ( y + heigth ), color_2, false ) );
        g.fill( _rectangle );
        if ( color_border != null ) {
            g.setPaint( color_border );
            g.draw( _rectangle );
        }
    }

    private double drawTaxonomyImage( final double x, final double y, final PhylogenyNode node, final Graphics2D g ) {
        final List<Uri> us = new ArrayList<Uri>();
        for( final Taxonomy t : node.getNodeData().getTaxonomies() ) {
            for( final Uri uri : t.getUris() ) {
                us.add( uri );
            }
        }
        double offset = 0;
        for( final Uri uri : us ) {
            if ( uri != null ) {
                final String uri_str = uri.getValue().toString().toLowerCase();
                if ( getImageMap().containsKey( uri_str ) ) {
                    final BufferedImage bi = getImageMap().get( uri_str );
                    if ( ( bi != null ) && ( bi.getHeight() > 5 ) && ( bi.getWidth() > 5 ) ) {
                        double scaling_factor = 1;
                        if ( getOptions().isAllowMagnificationOfTaxonomyImages()
                                || ( bi.getHeight() > ( 1.8 * getYdistance() ) ) ) {
                            scaling_factor = ( 1.8 * getYdistance() ) / bi.getHeight();
                        }
                        // y = y - ( 0.9 * getYdistance() );
                        final double hs = bi.getHeight() * scaling_factor;
                        double ws = ( bi.getWidth() * scaling_factor ) + offset;
                        final double my_y = y - ( 0.5 * hs );
                        final int x_w = ( int ) ( x + ws + 0.5 );
                        final int y_h = ( int ) ( my_y + hs + 0.5 );
                        if ( ( ( x_w - x ) > 7 ) && ( ( y_h - my_y ) > 7 ) ) {
                            g.drawImage( bi,
                                         ( int ) ( x + 0.5 + offset ),
                                         ( int ) ( my_y + 0.5 ),
                                         x_w,
                                         y_h,
                                         0,
                                         0,
                                         bi.getWidth(),
                                         bi.getHeight(),
                                         null );
                            ws += 8;
                        }
                        else {
                            ws = 0.0;
                        }
                        offset = ws;
                    }
                }
            }
        }
        return offset;
    }

    final private void errorMessageNoCutCopyPasteInUnrootedDisplay() {
        JOptionPane.showMessageDialog( this,
                                       "Cannot cut, copy, paste, add, or delete subtrees/nodes in unrooted display",
                                       "Attempt to cut/copy/paste/add/delete in unrooted display",
                                       JOptionPane.ERROR_MESSAGE );
    }

    private final Color getColorForFoundNode( final PhylogenyNode n ) {
        if ( isInCurrentExternalNodes( n ) ) {
            return getTreeColorSet().getFoundColor0();
        }
        else if ( isInFoundNodes0( n ) && !isInFoundNodes1( n ) ) {
            return getTreeColorSet().getFoundColor0();
        }
        else if ( !isInFoundNodes0( n ) && isInFoundNodes1( n ) ) {
            return getTreeColorSet().getFoundColor1();
        }
        else {
            return getTreeColorSet().getFoundColor0and1();
        }
    }

    final private Set<Long> getCopiedAndPastedNodes() {
        return getMainPanel().getCopiedAndPastedNodes();
    }

    final private Set<Long> getCurrentExternalNodes() {
        return _current_external_nodes;
    }

    final private Phylogeny getCutOrCopiedTree() {
        return getMainPanel().getCutOrCopiedTree();
    }

    private FontMetrics getFontMetricsForLargeDefaultFont() {
        return getTreeFontSet().getFontMetricsLarge();
    }

    private List<PhylogenyNode> getFoundNodesAsListOfPhylogenyNodes() {
        final List<PhylogenyNode> additional_nodes = new ArrayList<PhylogenyNode>();
        if ( getFoundNodes0() != null ) {
            for( final Long id : getFoundNodes0() ) {
                additional_nodes.add( _phylogeny.getNode( id ) );
            }
        }
        if ( getFoundNodes1() != null ) {
            for( final Long id : getFoundNodes1() ) {
                if ( ( getFoundNodes0() == null ) || !getFoundNodes0().contains( id ) ) {
                    additional_nodes.add( _phylogeny.getNode( id ) );
                }
            }
        }
        return additional_nodes;
    }

    final private float getLastDragPointX() {
        return _last_drag_point_x;
    }

    final private float getLastDragPointY() {
        return _last_drag_point_y;
    }

    final private short getMaxBranchesToLeaf( final PhylogenyNode node ) {
        if ( !_nodeid_dist_to_leaf.containsKey( node.getId() ) ) {
            final short m = PhylogenyMethods.calculateMaxBranchesToLeaf( node );
            _nodeid_dist_to_leaf.put( node.getId(), m );
            return m;
        }
        else {
            return _nodeid_dist_to_leaf.get( node.getId() );
        }
    }

    final private double getMaxDistanceToRoot() {
        if ( _max_distance_to_root < 0 ) {
            recalculateMaxDistanceToRoot();
        }
        return _max_distance_to_root;
    }

    final private float getOvMaxHeight() {
        return _ov_max_height;
    }

    final private float getOvMaxWidth() {
        return _ov_max_width;
    }

    final private float getOvXcorrectionFactor() {
        return _ov_x_correction_factor;
    }

    final private float getOvXDistance() {
        return _ov_x_distance;
    }

    final private int getOvXPosition() {
        return _ov_x_position;
    }

    final private float getOvYDistance() {
        return _ov_y_distance;
    }

    final private int getOvYPosition() {
        return _ov_y_position;
    }

    final private int getOvYStart() {
        return _ov_y_start;
    }

    final private List<Accession> getPdbAccs( final PhylogenyNode node ) {
        final List<Accession> pdb_ids = new ArrayList<Accession>();
        if ( node.getNodeData().isHasSequence() ) {
            final Sequence seq = node.getNodeData().getSequence();
            if ( !ForesterUtil.isEmpty( seq.getCrossReferences() ) ) {
                final SortedSet<Accession> cross_refs = seq.getCrossReferences();
                for( final Accession acc : cross_refs ) {
                    if ( acc.getSource().equalsIgnoreCase( "pdb" ) ) {
                        pdb_ids.add( acc );
                    }
                }
            }
        }
        return pdb_ids;
    }

    final private double getScaleDistance() {
        return _scale_distance;
    }

    final private String getScaleLabel() {
        return _scale_label;
    }

    final private TreeFontSet getTreeFontSet() {
        return getMainPanel().getTreeFontSet();
    }

    final private float getUrtFactor() {
        return _urt_factor;
    }

    final private float getUrtFactorOv() {
        return _urt_factor_ov;
    }

    final private void handleClickToAction( final NodeClickAction action, final PhylogenyNode node ) {
        switch ( action ) {
            case SHOW_DATA:
                showNodeFrame( node );
                break;
            case COLLAPSE:
                collapse( node );
                break;
            case REROOT:
                reRoot( node );
                break;
            case SUBTREE:
                subTree( node );
                break;
            case SWAP:
                swap( node );
                break;
            case COLOR_SUBTREE:
                colorSubtree( node );
                break;
            case COLOR_NODE_FONT:
                colorNodeFont( node );
                break;
            case CHANGE_NODE_FONT:
                changeNodeFont( node );
                break;
            case OPEN_SEQ_WEB:
                openSeqWeb( node );
                break;
            case BLAST:
                blast( node );
                break;
            case OPEN_TAX_WEB:
                openTaxWeb( node );
                break;
            case OPEN_PDB_WEB:
                openPdbWeb( node );
                break;
            case CUT_SUBTREE:
                cutSubtree( node );
                break;
            case COPY_SUBTREE:
                copySubtree( node );
                break;
            case PASTE_SUBTREE:
                pasteSubtree( node );
                break;
            case DELETE_NODE_OR_SUBTREE:
                deleteNodeOrSubtree( node );
                break;
            case ADD_NEW_NODE:
                addEmptyNode( node );
                break;
            case EDIT_NODE_DATA:
                showNodeEditFrame( node );
                break;
            case SELECT_NODES:
                selectNode( node );
                break;
            case SORT_DESCENDENTS:
                sortDescendants( node );
                break;
            case GET_EXT_DESC_DATA:
                showExtDescNodeData( node );
                break;
            default:
                throw new IllegalArgumentException( "unknown action: " + action );
        }
    }

    final private void increaseCurrentExternalNodesDataBufferChangeCounter() {
        _current_external_nodes_data_buffer_change_counter++;
    }

    final private void increaseOvSize() {
        if ( ( getOvMaxWidth() < ( getMainPanel().getCurrentScrollPane().getViewport().getVisibleRect().getWidth() / 2 ) )
                && ( getOvMaxHeight() < ( getMainPanel().getCurrentScrollPane().getViewport().getVisibleRect()
                        .getHeight() / 2 ) ) ) {
            setOvMaxWidth( getOvMaxWidth() + 5 );
            setOvMaxHeight( getOvMaxHeight() + 5 );
            updateOvSettings();
            getControlPanel().displayedPhylogenyMightHaveChanged( false );
        }
    }

    final private void init() {
        _color_chooser = new JColorChooser();
        _rollover_popup = new JTextArea();
        _rollover_popup.setFont( POPUP_FONT );
        resetNodeIdToDistToLeafMap();
        setTextAntialias();
        setTreeFile( null );
        setEdited( false );
        initializeOvSettings();
        setStartingAngle( ( TWO_PI * 3 ) / 4 );
        final ImageLoader il = new ImageLoader( this );
        new Thread( il ).start();
    }

    final private void initializeOvSettings() {
        setOvMaxHeight( getConfiguration().getOvMaxHeight() );
        setOvMaxWidth( getConfiguration().getOvMaxWidth() );
    }

    final private boolean inOvVirtualRectangle( final int x, final int y ) {
        return ( ( x >= ( getOvVirtualRectangle().x - 1 ) )
                && ( x <= ( getOvVirtualRectangle().x + getOvVirtualRectangle().width + 1 ) )
                && ( y >= ( getOvVirtualRectangle().y - 1 ) ) && ( y <= ( getOvVirtualRectangle().y
                + getOvVirtualRectangle().height + 1 ) ) );
    }

    final private boolean inOvVirtualRectangle( final MouseEvent e ) {
        return ( inOvVirtualRectangle( e.getX(), e.getY() ) );
    }

    final private boolean isCanBlast( final PhylogenyNode node ) {
        if ( !node.getNodeData().isHasSequence() && ForesterUtil.isEmpty( node.getName() ) ) {
            return false;
        }
        return Blast.isContainsQueryForBlast( node );
    }

    final private String isCanOpenSeqWeb( final PhylogenyNode node ) {
        final Accession a = SequenceAccessionTools.obtainAccessorFromDataFields( node );
        if ( a != null ) {
            return a.getValue();
        }
        return null;
    }

    final private boolean isCanOpenTaxWeb( final PhylogenyNode node ) {
        if ( node.getNodeData().isHasTaxonomy()
                && ( ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) )
                        || ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) )
                        || ( !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getCommonName() ) ) || ( ( node
                        .getNodeData().getTaxonomy().getIdentifier() != null ) && !ForesterUtil.isEmpty( node
                        .getNodeData().getTaxonomy().getIdentifier().getValue() ) ) ) ) {
            return true;
        }
        else {
            return false;
        }
    }

    final private boolean isInCurrentExternalNodes( final PhylogenyNode node ) {
        return ( ( getCurrentExternalNodes() != null ) && getCurrentExternalNodes().contains( node.getId() ) );
    }

    private boolean isInFoundNodes( final PhylogenyNode n ) {
        return isInFoundNodes0( n ) || isInFoundNodes1( n );
    }

    final private boolean isInFoundNodes0( final PhylogenyNode node ) {
        return ( ( getFoundNodes0() != null ) && getFoundNodes0().contains( node.getId() ) );
    }

    final private boolean isInFoundNodes1( final PhylogenyNode node ) {
        return ( ( getFoundNodes1() != null ) && getFoundNodes1().contains( node.getId() ) );
    }

    final private boolean isInOv() {
        return _in_ov;
    }

    final private boolean isNodeDataInvisible( final PhylogenyNode node ) {
        int y_dist = 40;
        if ( getControlPanel().isShowTaxonomyImages() ) {
            y_dist = 40 + ( int ) getYdistance();
        }
        return ( ( node.getYcoord() < ( getVisibleRect().getMinY() - y_dist ) )
                || ( node.getYcoord() > ( getVisibleRect().getMaxY() + y_dist ) ) || ( ( node.getParent() != null ) && ( node
                .getParent().getXcoord() > getVisibleRect().getMaxX() ) ) );
    }

    final private boolean isNodeDataInvisibleUnrootedCirc( final PhylogenyNode node ) {
        return ( ( node.getYcoord() < ( getVisibleRect().getMinY() - 20 ) )
                || ( node.getYcoord() > ( getVisibleRect().getMaxY() + 20 ) )
                || ( node.getXcoord() < ( getVisibleRect().getMinX() - 20 ) ) || ( node.getXcoord() > ( getVisibleRect()
                .getMaxX() + 20 ) ) );
    }

    final private boolean isNonLinedUpCladogram() {
        return getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP;
    }

    final private boolean isUniformBranchLengthsForCladogram() {
        return getOptions().getCladogramType() == CLADOGRAM_TYPE.TOTAL_NODE_SUM_DEP;
    }

    final private void keyPressedCalls( final KeyEvent e ) {
        if ( isOvOn() && ( getMousePosition() != null ) && ( getMousePosition().getLocation() != null ) ) {
            if ( inOvVirtualRectangle( getMousePosition().x, getMousePosition().y ) ) {
                if ( !isInOvRect() ) {
                    setInOvRect( true );
                }
            }
            else if ( isInOvRect() ) {
                setInOvRect( false );
            }
        }
        if ( e.getModifiersEx() == InputEvent.CTRL_DOWN_MASK ) {
            if ( ( e.getKeyCode() == KeyEvent.VK_DELETE ) || ( e.getKeyCode() == KeyEvent.VK_HOME )
                    || ( e.getKeyCode() == KeyEvent.VK_F ) ) {
                getMainPanel().getTreeFontSet().mediumFonts();
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
            else if ( ( e.getKeyCode() == KeyEvent.VK_SUBTRACT ) || ( e.getKeyCode() == KeyEvent.VK_MINUS ) ) {
                getMainPanel().getTreeFontSet().decreaseFontSize( 1, false );
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
            else if ( plusPressed( e.getKeyCode() ) ) {
                getMainPanel().getTreeFontSet().increaseFontSize();
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
            }
        }
        else {
            if ( ( e.getKeyCode() == KeyEvent.VK_DELETE ) || ( e.getKeyCode() == KeyEvent.VK_HOME )
                    || ( e.getKeyCode() == KeyEvent.VK_F ) ) {
                getControlPanel().showWhole();
            }
            else if ( ( e.getKeyCode() == KeyEvent.VK_UP ) || ( e.getKeyCode() == KeyEvent.VK_DOWN )
                    || ( e.getKeyCode() == KeyEvent.VK_LEFT ) || ( e.getKeyCode() == KeyEvent.VK_RIGHT ) ) {
                if ( e.getModifiersEx() == InputEvent.SHIFT_DOWN_MASK ) {
                    if ( e.getKeyCode() == KeyEvent.VK_UP ) {
                        getMainPanel().getControlPanel().zoomInY( Constants.WHEEL_ZOOM_IN_FACTOR );
                        getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_DOWN ) {
                        getMainPanel().getControlPanel().zoomOutY( Constants.WHEEL_ZOOM_OUT_FACTOR );
                        getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_LEFT ) {
                        getMainPanel().getControlPanel().zoomOutX( Constants.WHEEL_ZOOM_OUT_FACTOR,
                                                                   Constants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR );
                        getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_RIGHT ) {
                        getMainPanel().getControlPanel().zoomInX( Constants.WHEEL_ZOOM_IN_FACTOR,
                                                                  Constants.WHEEL_ZOOM_IN_FACTOR );
                        getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
                    }
                }
                else {
                    final int d = 80;
                    int dx = 0;
                    int dy = -d;
                    if ( e.getKeyCode() == KeyEvent.VK_DOWN ) {
                        dy = d;
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_LEFT ) {
                        dx = -d;
                        dy = 0;
                    }
                    else if ( e.getKeyCode() == KeyEvent.VK_RIGHT ) {
                        dx = d;
                        dy = 0;
                    }
                    final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
                    scroll_position.x = scroll_position.x + dx;
                    scroll_position.y = scroll_position.y + dy;
                    if ( scroll_position.x <= 0 ) {
                        scroll_position.x = 0;
                    }
                    else {
                        final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                                - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
                        if ( scroll_position.x >= max_x ) {
                            scroll_position.x = max_x;
                        }
                    }
                    if ( scroll_position.y <= 0 ) {
                        scroll_position.y = 0;
                    }
                    else {
                        final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                                - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
                        if ( scroll_position.y >= max_y ) {
                            scroll_position.y = max_y;
                        }
                    }
                    repaint();
                    getMainPanel().getCurrentScrollPane().getViewport().setViewPosition( scroll_position );
                }
            }
            else if ( ( e.getKeyCode() == KeyEvent.VK_SUBTRACT ) || ( e.getKeyCode() == KeyEvent.VK_MINUS ) ) {
                getMainPanel().getControlPanel().zoomOutY( Constants.WHEEL_ZOOM_OUT_FACTOR );
                getMainPanel().getControlPanel().zoomOutX( Constants.WHEEL_ZOOM_OUT_FACTOR,
                                                           Constants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR );
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
            }
            else if ( plusPressed( e.getKeyCode() ) ) {
                getMainPanel().getControlPanel().zoomInX( Constants.WHEEL_ZOOM_IN_FACTOR,
                                                          Constants.WHEEL_ZOOM_IN_FACTOR );
                getMainPanel().getControlPanel().zoomInY( Constants.WHEEL_ZOOM_IN_FACTOR );
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( false );
            }
            else if ( e.getKeyCode() == KeyEvent.VK_S ) {
                if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                        || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                    setStartingAngle( ( getStartingAngle() % TWO_PI ) + ANGLE_ROTATION_UNIT );
                    getControlPanel().displayedPhylogenyMightHaveChanged( false );
                }
            }
            else if ( e.getKeyCode() == KeyEvent.VK_A ) {
                if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                        || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
                    setStartingAngle( ( getStartingAngle() % TWO_PI ) - ANGLE_ROTATION_UNIT );
                    if ( getStartingAngle() < 0 ) {
                        setStartingAngle( TWO_PI + getStartingAngle() );
                    }
                    getControlPanel().displayedPhylogenyMightHaveChanged( false );
                }
            }
            else if ( e.getKeyCode() == KeyEvent.VK_D ) {
                boolean selected = false;
                if ( getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.HORIZONTAL ) {
                    getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
                    selected = true;
                }
                else {
                    getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.HORIZONTAL );
                }
                if ( getMainPanel().getMainFrame() == null ) {
                    // Must be "E" applet version.
                    final ArchaeopteryxE ae = ( ArchaeopteryxE ) ( ( MainPanelApplets ) getMainPanel() ).getApplet();
                    if ( ae.getlabelDirectionCbmi() != null ) {
                        ae.getlabelDirectionCbmi().setSelected( selected );
                    }
                }
                else {
                    getMainPanel().getMainFrame().getlabelDirectionCbmi().setSelected( selected );
                }
                repaint();
            }
            else if ( e.getKeyCode() == KeyEvent.VK_X ) {
                switchDisplaygetPhylogenyGraphicsType();
                repaint();
            }
            else if ( e.getKeyCode() == KeyEvent.VK_C ) {
                cycleColors();
                repaint();
            }
            else if ( getOptions().isShowOverview() && isOvOn() && ( e.getKeyCode() == KeyEvent.VK_O ) ) {
                MainFrame.cycleOverview( getOptions(), this );
                repaint();
            }
            else if ( getOptions().isShowOverview() && isOvOn() && ( e.getKeyCode() == KeyEvent.VK_I ) ) {
                increaseOvSize();
            }
            else if ( getOptions().isShowOverview() && isOvOn() && ( e.getKeyCode() == KeyEvent.VK_U ) ) {
                decreaseOvSize();
            }
            e.consume();
        }
    }

    final private void makePopupMenus( final PhylogenyNode node ) {
        _node_popup_menu = new JPopupMenu();
        final List<String> clickto_names = _main_panel.getControlPanel().getSingleClickToNames();
        _node_popup_menu_items = new JMenuItem[ clickto_names.size() ];
        for( int i = 0; i < clickto_names.size(); i++ ) {
            final String title = clickto_names.get( i );
            _node_popup_menu_items[ i ] = new JMenuItem( title );
            if ( title.equals( Configuration.clickto_options[ Configuration.open_seq_web ][ 0 ] ) ) {
                final String id = isCanOpenSeqWeb( node );
                if ( !ForesterUtil.isEmpty( id ) ) {
                    _node_popup_menu_items[ i ].setText( _node_popup_menu_items[ i ].getText() + " [" + id + "]" );
                    _node_popup_menu_items[ i ].setEnabled( true );
                }
                else {
                    _node_popup_menu_items[ i ].setEnabled( false );
                }
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.open_pdb_web ][ 0 ] ) ) {
                final List<Accession> accs = getPdbAccs( node );
                _node_popup_menu_items[ i ] = new JMenuItem( title );
                if ( !ForesterUtil.isEmpty( accs ) ) {
                    if ( accs.size() == 1 ) {
                        _node_popup_menu_items[ i ].setText( _node_popup_menu_items[ i ].getText() + " ["
                                + TreePanelUtil.pdbAccToString( accs, 0 ) + "]" );
                        _node_popup_menu_items[ i ].setEnabled( true );
                    }
                    else if ( accs.size() == 2 ) {
                        _node_popup_menu_items[ i ].setText( _node_popup_menu_items[ i ].getText() + " ["
                                + TreePanelUtil.pdbAccToString( accs, 0 ) + ", "
                                + TreePanelUtil.pdbAccToString( accs, 1 ) + "]" );
                        _node_popup_menu_items[ i ].setEnabled( true );
                    }
                    else if ( accs.size() == 3 ) {
                        _node_popup_menu_items[ i ].setText( _node_popup_menu_items[ i ].getText() + " ["
                                + TreePanelUtil.pdbAccToString( accs, 0 ) + ", "
                                + TreePanelUtil.pdbAccToString( accs, 1 ) + ", "
                                + TreePanelUtil.pdbAccToString( accs, 2 ) + "]" );
                        _node_popup_menu_items[ i ].setEnabled( true );
                    }
                    else {
                        _node_popup_menu_items[ i ].setText( _node_popup_menu_items[ i ].getText() + " ["
                                + TreePanelUtil.pdbAccToString( accs, 0 ) + ", "
                                + TreePanelUtil.pdbAccToString( accs, 1 ) + ", "
                                + TreePanelUtil.pdbAccToString( accs, 2 ) + ", + " + ( accs.size() - 3 ) + " more]" );
                        _node_popup_menu_items[ i ].setEnabled( true );
                    }
                }
                else {
                    _node_popup_menu_items[ i ].setEnabled( false );
                }
                //
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.open_tax_web ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanOpenTaxWeb( node ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.blast ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanBlast( node ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.delete_subtree_or_node ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanDelete() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.cut_subtree ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanCut( node ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.copy_subtree ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanCopy() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.paste_subtree ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
                _node_popup_menu_items[ i ].setEnabled( isCanPaste() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.edit_node_data ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.add_new_node ][ 0 ] ) ) {
                if ( !getOptions().isEditable() ) {
                    continue;
                }
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.reroot ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanReroot() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.collapse_uncollapse ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( ( isCanCollapse() && !node.isExternal() ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.color_subtree ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanColorSubtree() );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.subtree ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( isCanSubtree( node ) );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.swap ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( node.getNumberOfDescendants() == 2 );
            }
            else if ( title.equals( Configuration.clickto_options[ Configuration.sort_descendents ][ 0 ] ) ) {
                _node_popup_menu_items[ i ].setEnabled( node.getNumberOfDescendants() > 1 );
            }
            _node_popup_menu_items[ i ].addActionListener( this );
            _node_popup_menu.add( _node_popup_menu_items[ i ] );
        }
    }

    private final void nodeDataAsSB( final PhylogenyNode node, final StringBuilder sb ) {
        if ( getControlPanel().isShowNodeNames() && ( node.getName().length() > 0 ) ) {
            if ( sb.length() > 0 ) {
                sb.append( " " );
            }
            sb.append( node.getName() );
        }
        if ( node.getNodeData().isHasSequence() ) {
            if ( getControlPanel().isShowSeqSymbols() && ( node.getNodeData().getSequence().getSymbol().length() > 0 ) ) {
                if ( sb.length() > 0 ) {
                    sb.append( " " );
                }
                sb.append( node.getNodeData().getSequence().getSymbol() );
            }
            if ( getControlPanel().isShowGeneNames() && ( node.getNodeData().getSequence().getGeneName().length() > 0 ) ) {
                if ( sb.length() > 0 ) {
                    sb.append( " " );
                }
                sb.append( node.getNodeData().getSequence().getGeneName() );
            }
            if ( getControlPanel().isShowSeqNames() && ( node.getNodeData().getSequence().getName().length() > 0 ) ) {
                if ( sb.length() > 0 ) {
                    sb.append( " " );
                }
                sb.append( node.getNodeData().getSequence().getName() );
            }
            if ( getControlPanel().isShowSequenceAcc() && ( node.getNodeData().getSequence().getAccession() != null ) ) {
                if ( sb.length() > 0 ) {
                    sb.append( " " );
                }
                if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() ) ) {
                    sb.append( node.getNodeData().getSequence().getAccession().getSource() );
                    sb.append( ":" );
                }
                sb.append( node.getNodeData().getSequence().getAccession().getValue() );
            }
        }
        if ( getControlPanel().isShowProperties() && node.getNodeData().isHasProperties() ) {
            if ( sb.length() > 0 ) {
                sb.append( " " );
            }
            sb.append( propertiesToString( node ) );
        }
    }

    private final void nodeTaxonomyDataAsSB( final Taxonomy taxonomy, final StringBuilder sb ) {
        if ( _control_panel.isShowTaxonomyCode() && !ForesterUtil.isEmpty( taxonomy.getTaxonomyCode() ) ) {
            sb.append( taxonomy.getTaxonomyCode() );
            sb.append( " " );
        }
        if ( _control_panel.isShowTaxonomyScientificNames() && _control_panel.isShowTaxonomyCommonNames() ) {
            if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() )
                    && !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                if ( getOptions().isAbbreviateScientificTaxonNames()
                        && ( taxonomy.getScientificName().indexOf( ' ' ) > 0 ) ) {
                    abbreviateScientificName( taxonomy.getScientificName(), sb );
                }
                else {
                    sb.append( taxonomy.getScientificName() );
                }
                sb.append( " (" );
                sb.append( taxonomy.getCommonName() );
                sb.append( ") " );
            }
            else if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() ) ) {
                if ( getOptions().isAbbreviateScientificTaxonNames()
                        && ( taxonomy.getScientificName().indexOf( ' ' ) > 0 ) ) {
                    abbreviateScientificName( taxonomy.getScientificName(), sb );
                }
                else {
                    sb.append( taxonomy.getScientificName() );
                }
                sb.append( " " );
            }
            else if ( !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                sb.append( taxonomy.getCommonName() );
                sb.append( " " );
            }
        }
        else if ( _control_panel.isShowTaxonomyScientificNames() ) {
            if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() ) ) {
                if ( getOptions().isAbbreviateScientificTaxonNames()
                        && ( taxonomy.getScientificName().indexOf( ' ' ) > 0 ) ) {
                    abbreviateScientificName( taxonomy.getScientificName(), sb );
                }
                else {
                    sb.append( taxonomy.getScientificName() );
                }
                sb.append( " " );
            }
        }
        else if ( _control_panel.isShowTaxonomyCommonNames() ) {
            if ( !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                sb.append( taxonomy.getCommonName() );
                sb.append( " " );
            }
        }
    }

    private final String obtainTitleForExtDescNodeData() {
        switch ( getOptions().getExtDescNodeDataToReturn() ) {
            case NODE_NAME:
                return "Node Names";
            case GENE_NAME:
                return "Gene Names";
            case SEQUENCE_NAME:
                return "Sequence Names";
            case SEQUENCE_SYMBOL:
                return "Sequence Symbols";
            case SEQUENCE_MOL_SEQ:
                return "Molecular Sequences";
            case SEQUENCE_MOL_SEQ_FASTA:
                return "Molecular Sequences (Fasta)";
            case SEQUENCE_ACC:
                return "Sequence Accessors";
            case TAXONOMY_SCIENTIFIC_NAME:
                return "Scientific Names";
            case TAXONOMY_CODE:
                return "Taxonomy Codes";
            case TAXONOMY_COMM0N_NAME:
                return "Taxonomy Common Names";
            case UNKNOWN:
                return "User Selected Data";
            default:
                throw new IllegalArgumentException( "unknown data element: "
                        + getOptions().getExtDescNodeDataToReturn() );
        }
    }

    final private void openPdbWeb( final PhylogenyNode node ) {
        final List<Accession> pdb_ids = getPdbAccs( node );
        if ( ForesterUtil.isEmpty( pdb_ids ) ) {
            cannotOpenBrowserWarningMessage( "PDB" );
            return;
        }
        final List<String> uri_strs = TreePanelUtil.createUrisForPdbWeb( node, pdb_ids, getConfiguration(), this );
        if ( !ForesterUtil.isEmpty( uri_strs ) ) {
            for( final String uri_str : uri_strs ) {
                try {
                    AptxUtil.launchWebBrowser( new URI( uri_str ),
                                               isApplet(),
                                               isApplet() ? obtainApplet() : null,
                                               "_aptx_seq" );
                }
                catch ( final IOException e ) {
                    AptxUtil.showErrorMessage( this, e.toString() );
                    e.printStackTrace();
                }
                catch ( final URISyntaxException e ) {
                    AptxUtil.showErrorMessage( this, e.toString() );
                    e.printStackTrace();
                }
            }
        }
        else {
            cannotOpenBrowserWarningMessage( "PDB" );
        }
    }

    final private void openSeqWeb( final PhylogenyNode node ) {
        if ( ForesterUtil.isEmpty( isCanOpenSeqWeb( node ) ) ) {
            cannotOpenBrowserWarningMessage( "sequence" );
            return;
        }
        final String uri_str = TreePanelUtil.createUriForSeqWeb( node, getConfiguration(), this );
        if ( !ForesterUtil.isEmpty( uri_str ) ) {
            try {
                AptxUtil.launchWebBrowser( new URI( uri_str ),
                                           isApplet(),
                                           isApplet() ? obtainApplet() : null,
                                           "_aptx_seq" );
            }
            catch ( final IOException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
            catch ( final URISyntaxException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else {
            cannotOpenBrowserWarningMessage( "sequence" );
        }
    }

    final private void openTaxWeb( final PhylogenyNode node ) {
        if ( !isCanOpenTaxWeb( node ) ) {
            cannotOpenBrowserWarningMessage( "taxonomic" );
            return;
        }
        String uri_str = null;
        final Taxonomy tax = node.getNodeData().getTaxonomy();
        if ( ( tax.getIdentifier() != null ) && !ForesterUtil.isEmpty( tax.getIdentifier().getValue() )
                && tax.getIdentifier().getValue().startsWith( "http://" ) ) {
            try {
                uri_str = new URI( tax.getIdentifier().getValue() ).toString();
            }
            catch ( final URISyntaxException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                uri_str = null;
                e.printStackTrace();
            }
        }
        else if ( ( tax.getIdentifier() != null )
                && !ForesterUtil.isEmpty( tax.getIdentifier().getValue() )
                && !ForesterUtil.isEmpty( tax.getIdentifier().getProvider() )
                && ( tax.getIdentifier().getProvider().equalsIgnoreCase( "ncbi" ) || tax.getIdentifier().getProvider()
                        .equalsIgnoreCase( "uniprot" ) ) ) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/"
                        + URLEncoder.encode( tax.getIdentifier().getValue(), ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/?query="
                        + URLEncoder.encode( tax.getScientificName(), ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/?query="
                        + URLEncoder.encode( tax.getTaxonomyCode(), ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else if ( !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/?query="
                        + URLEncoder.encode( tax.getCommonName(), ForesterConstants.UTF8 );
            }
            catch ( final UnsupportedEncodingException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        if ( !ForesterUtil.isEmpty( uri_str ) ) {
            try {
                AptxUtil.launchWebBrowser( new URI( uri_str ),
                                           isApplet(),
                                           isApplet() ? obtainApplet() : null,
                                           "_aptx_tax" );
            }
            catch ( final IOException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
            catch ( final URISyntaxException e ) {
                AptxUtil.showErrorMessage( this, e.toString() );
                e.printStackTrace();
            }
        }
        else {
            cannotOpenBrowserWarningMessage( "taxonomic" );
        }
    }

    final private void paintBranchLength( final Graphics2D g,
                                          final PhylogenyNode node,
                                          final boolean to_pdf,
                                          final boolean to_graphics_file ) {
        g.setFont( getTreeFontSet().getSmallFont() );
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else {
            g.setColor( getTreeColorSet().getBranchLengthColor() );
        }
        if ( !node.isRoot() ) {
            if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                TreePanel.drawString( FORMATTER_BRANCH_LENGTH.format( node.getDistanceToParent() ), node.getParent()
                        .getXcoord() + EURO_D, node.getYcoord() - getTreeFontSet().getSmallMaxDescent(), g );
            }
            else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                TreePanel.drawString( FORMATTER_BRANCH_LENGTH.format( node.getDistanceToParent() ), node.getParent()
                        .getXcoord() + ROUNDED_D, node.getYcoord() - getTreeFontSet().getSmallMaxDescent(), g );
            }
            else {
                TreePanel.drawString( FORMATTER_BRANCH_LENGTH.format( node.getDistanceToParent() ), node.getParent()
                        .getXcoord() + 3, node.getYcoord() - getTreeFontSet().getSmallMaxDescent(), g );
            }
        }
        else {
            TreePanel.drawString( FORMATTER_BRANCH_LENGTH.format( node.getDistanceToParent() ), 3, node.getYcoord()
                    - getTreeFontSet().getSmallMaxDescent(), g );
        }
    }

    final private void paintBranchLite( final Graphics2D g,
                                        final float x1,
                                        final float x2,
                                        final float y1,
                                        final float y2,
                                        final PhylogenyNode node ) {
        g.setColor( getTreeColorSet().getOvColor() );
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR ) {
            drawLine( x1, y1, x2, y2, g );
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CONVEX ) {
            _quad_curve.setCurve( x1, y1, x1, y2, x2, y2 );
            ( g ).draw( _quad_curve );
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CURVED ) {
            final float dx = x2 - x1;
            final float dy = y2 - y1;
            _cubic_curve.setCurve( x1, y1, x1 + ( dx * 0.4f ), y1 + ( dy * 0.2f ), x1 + ( dx * 0.6f ), y1
                    + ( dy * 0.8f ), x2, y2 );
            ( g ).draw( _cubic_curve );
        }
        else {
            final float x2a = x2;
            final float x1a = x1;
            // draw the vertical line
            if ( node.isFirstChildNode() || node.isLastChildNode() ) {
                drawLine( x1, y1, x1, y2, g );
            }
            // draw the horizontal line
            drawLine( x1a, y2, x2a, y2, g );
        }
    }

    /**
     * Paint a branch which consists of a vertical and a horizontal bar
     * @param is_ind_found_nodes 
     */
    final private void paintBranchRectangular( final Graphics2D g,
                                               final float x1,
                                               final float x2,
                                               final float y1,
                                               final float y2,
                                               final PhylogenyNode node,
                                               final boolean to_pdf,
                                               final boolean to_graphics_file ) {
        assignGraphicsForBranchWithColorForParentBranch( node, false, g, to_pdf, to_graphics_file );
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR ) {
            drawLine( x1, y1, x2, y2, g );
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CONVEX ) {
            _quad_curve.setCurve( x1, y1, x1, y2, x2, y2 );
            g.draw( _quad_curve );
        }
        else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CURVED ) {
            final float dx = x2 - x1;
            final float dy = y2 - y1;
            _cubic_curve.setCurve( x1, y1, x1 + ( dx * 0.4f ), y1 + ( dy * 0.2f ), x1 + ( dx * 0.6f ), y1
                    + ( dy * 0.8f ), x2, y2 );
            g.draw( _cubic_curve );
        }
        else {
            final float x2a = x2;
            final float x1a = x1;
            float y2_r = 0;
            if ( node.isFirstChildNode() || node.isLastChildNode()
                    || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE )
                    || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) ) {
                if ( !to_graphics_file
                        && !to_pdf
                        && ( ( ( y2 < ( getVisibleRect().getMinY() - 20 ) ) && ( y1 < ( getVisibleRect().getMinY() - 20 ) ) ) || ( ( y2 > ( getVisibleRect()
                                .getMaxY() + 20 ) ) && ( y1 > ( getVisibleRect().getMaxY() + 20 ) ) ) ) ) {
                    // Do nothing.
                }
                else {
                    if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                        float x2c = x1 + EURO_D;
                        if ( x2c > x2a ) {
                            x2c = x2a;
                        }
                        drawLine( x1, y1, x2c, y2, g );
                    }
                    else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                        if ( y2 > y1 ) {
                            y2_r = y2 - ROUNDED_D;
                            if ( y2_r < y1 ) {
                                y2_r = y1;
                            }
                            drawLine( x1, y1, x1, y2_r, g );
                        }
                        else {
                            y2_r = y2 + ROUNDED_D;
                            if ( y2_r > y1 ) {
                                y2_r = y1;
                            }
                            drawLine( x1, y1, x1, y2_r, g );
                        }
                    }
                    else {
                        drawLine( x1, y1, x1, y2, g );
                    }
                }
            }
            // draw the horizontal line
            if ( !to_graphics_file && !to_pdf
                    && ( ( y2 < ( getVisibleRect().getMinY() - 20 ) ) || ( y2 > ( getVisibleRect().getMaxY() + 20 ) ) ) ) {
                return;
            }
            float x1_r = 0;
            if ( !getControlPanel().isWidthBranches() || ( PhylogenyMethods.getBranchWidthValue( node ) == 1 ) ) {
                if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                    x1_r = x1a + ROUNDED_D;
                    if ( x1_r < x2a ) {
                        drawLine( x1_r, y2, x2a, y2, g );
                    }
                }
                else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                    final float x1c = x1a + EURO_D;
                    if ( x1c < x2a ) {
                        drawLine( x1c, y2, x2a, y2, g );
                    }
                }
                else {
                    drawLine( x1a, y2, x2a, y2, g );
                }
            }
            else {
                final double w = PhylogenyMethods.getBranchWidthValue( node );
                if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                    x1_r = x1a + ROUNDED_D;
                    if ( x1_r < x2a ) {
                        drawRectFilled( x1_r, y2 - ( w / 2 ), x2a - x1_r, w, g );
                    }
                }
                else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                    final float x1c = x1a + EURO_D;
                    if ( x1c < x2a ) {
                        drawRectFilled( x1c, y2 - ( w / 2 ), x2a - x1c, w, g );
                    }
                }
                else {
                    drawRectFilled( x1a, y2 - ( w / 2 ), x2a - x1a, w, g );
                }
            }
            if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) ) {
                if ( x1_r > x2a ) {
                    x1_r = x2a;
                }
                if ( y2 > y2_r ) {
                    final double diff = y2 - y2_r;
                    _arc.setArc( x1, y2_r - diff, 2 * ( x1_r - x1 ), 2 * diff, 180, 90, Arc2D.OPEN );
                }
                else {
                    _arc.setArc( x1, y2, 2 * ( x1_r - x1 ), 2 * ( y2_r - y2 ), 90, 90, Arc2D.OPEN );
                }
                g.draw( _arc );
            }
        }
        if ( node.isExternal() ) {
            paintNodeBox( x2, y2, node, g, to_pdf, to_graphics_file );
        }
    }

    final private double paintCirculars( final PhylogenyNode n,
                                         final Phylogeny phy,
                                         final float center_x,
                                         final float center_y,
                                         final double radius,
                                         final boolean radial_labels,
                                         final Graphics2D g,
                                         final boolean to_pdf,
                                         final boolean to_graphics_file ) {
        if ( n.isExternal() || n.isCollapse() ) { //~~circ collapse
            if ( !_urt_nodeid_angle_map.containsKey( n.getId() ) ) {
                System.out.println( "no " + n + " =====>>>>>>> ERROR!" );//TODO
            }
            return _urt_nodeid_angle_map.get( n.getId() );
        }
        else {
            final List<PhylogenyNode> descs = n.getDescendants();
            double sum = 0;
            for( final PhylogenyNode desc : descs ) {
                sum += paintCirculars( desc,
                                       phy,
                                       center_x,
                                       center_y,
                                       radius,
                                       radial_labels,
                                       g,
                                       to_pdf,
                                       to_graphics_file );
            }
            double r = 0;
            if ( !n.isRoot() ) {
                r = 1 - ( ( ( double ) _circ_max_depth - n.calculateDepth() ) / _circ_max_depth );
            }
            final double theta = sum / descs.size();
            n.setXcoord( ( float ) ( center_x + ( r * radius * Math.cos( theta ) ) ) );
            n.setYcoord( ( float ) ( center_y + ( r * radius * Math.sin( theta ) ) ) );
            _urt_nodeid_angle_map.put( n.getId(), theta );
            for( final PhylogenyNode desc : descs ) {
                paintBranchCircular( n, desc, g, radial_labels, to_pdf, to_graphics_file );
            }
            return theta;
        }
    }

    final private void paintCircularsLite( final PhylogenyNode n,
                                           final Phylogeny phy,
                                           final int center_x,
                                           final int center_y,
                                           final int radius,
                                           final Graphics2D g ) {
        if ( n.isExternal() ) {
            return;
        }
        else {
            final List<PhylogenyNode> descs = n.getDescendants();
            for( final PhylogenyNode desc : descs ) {
                paintCircularsLite( desc, phy, center_x, center_y, radius, g );
            }
            float r = 0;
            if ( !n.isRoot() ) {
                r = 1 - ( ( ( float ) _circ_max_depth - n.calculateDepth() ) / _circ_max_depth );
            }
            final double theta = _urt_nodeid_angle_map.get( n.getId() );
            n.setXSecondary( ( float ) ( center_x + ( radius * r * Math.cos( theta ) ) ) );
            n.setYSecondary( ( float ) ( center_y + ( radius * r * Math.sin( theta ) ) ) );
            for( final PhylogenyNode desc : descs ) {
                paintBranchCircularLite( n, desc, g );
            }
        }
    }

    final private void paintCollapsedNode( final Graphics2D g,
                                           final PhylogenyNode node,
                                           final boolean to_graphics_file,
                                           final boolean to_pdf,
                                           final boolean is_in_found_nodes ) {
        Color c = null;
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            c = Color.BLACK;
        }
        else if ( is_in_found_nodes ) {
            c = getColorForFoundNode( node );
        }
        else if ( getControlPanel().isColorAccordingToSequence() ) {
            c = getSequenceBasedColor( node );
        }
        else if ( getControlPanel().isColorAccordingToTaxonomy() ) {
            c = getTaxonomyBasedColor( node );
        }
        else if ( getOptions().isColorLabelsSameAsParentBranch() && getControlPanel().isUseVisualStyles()
                && ( PhylogenyMethods.getBranchColorValue( node ) != null ) ) {
            c = PhylogenyMethods.getBranchColorValue( node );
        }
        else {
            c = getTreeColorSet().getCollapseFillColor();
        }
        double d = node.getAllExternalDescendants().size();
        if ( d > 1000 ) {
            d = ( 3 * _y_distance ) / 3;
        }
        else {
            d = ( Math.log10( d ) * _y_distance ) / 2.5;
        }
        final int box_size = getOptions().getDefaultNodeShapeSize() + 1;
        if ( d < box_size ) {
            d = box_size;
        }
        final float xx = node.getXcoord() - ( 2 * box_size );
        final float xxx = xx > node.getParent().getXcoord() + 1 ? xx : node.getParent().getXcoord() + 1;
        _polygon.reset();
        _polygon.moveTo( xxx, node.getYcoord() );
        _polygon.lineTo( node.getXcoord() + 1, node.getYcoord() - d );
        _polygon.lineTo( node.getXcoord() + 1, node.getYcoord() + d );
        _polygon.closePath();
        if ( getOptions().getDefaultNodeFill() == NodeVisualData.NodeFill.SOLID ) {
            g.setColor( c );
            g.fill( _polygon );
        }
        else if ( getOptions().getDefaultNodeFill() == NodeVisualData.NodeFill.NONE ) {
            g.setColor( getBackground() );
            g.fill( _polygon );
            g.setColor( c );
            g.draw( _polygon );
        }
        else if ( getOptions().getDefaultNodeFill() == NodeFill.GRADIENT ) {
            g.setPaint( new GradientPaint( xxx, node.getYcoord(), getBackground(), node.getXcoord(), ( float ) ( node
                    .getYcoord() - d ), c, false ) );
            g.fill( _polygon );
            g.setPaint( c );
            g.draw( _polygon );
        }
        paintNodeData( g, node, to_graphics_file, to_pdf, is_in_found_nodes );
    }

    final private void paintConfidenceValues( final Graphics2D g,
                                              final PhylogenyNode node,
                                              final boolean to_pdf,
                                              final boolean to_graphics_file ) {
        final List<Confidence> confidences = node.getBranchData().getConfidences();
        boolean not_first = false;
        Collections.sort( confidences );
        final StringBuilder sb = new StringBuilder();
        for( final Confidence confidence : confidences ) {
            final double value = confidence.getValue();
            if ( value != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                if ( value < getOptions().getMinConfidenceValue() ) {
                    return;
                }
                if ( not_first ) {
                    sb.append( "/" );
                }
                else {
                    not_first = true;
                }
                sb.append( FORMATTER_CONFIDENCE.format( ForesterUtil.round( value, getOptions()
                        .getNumberOfDigitsAfterCommaForConfidenceValues() ) ) );
                if ( getOptions().isShowConfidenceStddev() ) {
                    if ( confidence.getStandardDeviation() != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                        sb.append( "(" );
                        sb.append( FORMATTER_CONFIDENCE.format( ForesterUtil.round( confidence.getStandardDeviation(),
                                                                                    getOptions()
                                                                                            .getNumberOfDigitsAfterCommaForConfidenceValues() ) ) );
                        sb.append( ")" );
                    }
                }
            }
        }
        if ( sb.length() > 0 ) {
            final float parent_x = node.getParent().getXcoord();
            float x = node.getXcoord();
            g.setFont( getTreeFontSet().getSmallFont() );
            if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
                x += EURO_D;
            }
            else if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
                x += ROUNDED_D;
            }
            if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
                g.setColor( Color.BLACK );
            }
            else {
                g.setColor( getTreeColorSet().getConfidenceColor() );
            }
            final String conf_str = sb.toString();
            TreePanel.drawString( conf_str,
                                  parent_x
                                          + ( ( x - parent_x - getTreeFontSet().getFontMetricsSmall()
                                                  .stringWidth( conf_str ) ) / 2 ),
                                  ( node.getYcoord() + getTreeFontSet().getSmallMaxAscent() ) - 1,
                                  g );
        }
    }

    final private void paintGainedAndLostCharacters( final Graphics2D g,
                                                     final PhylogenyNode node,
                                                     final String gained,
                                                     final String lost ) {
        if ( node.getParent() != null ) {
            final float parent_x = node.getParent().getXcoord();
            final float x = node.getXcoord();
            g.setFont( getTreeFontSet().getLargeFont() );
            g.setColor( getTreeColorSet().getGainedCharactersColor() );
            if ( Constants.SPECIAL_CUSTOM ) {
                g.setColor( Color.BLUE );
            }
            TreePanel
                    .drawString( gained,
                                 parent_x
                                         + ( ( x - parent_x - getFontMetricsForLargeDefaultFont().stringWidth( gained ) ) / 2 ),
                                 ( node.getYcoord() - getFontMetricsForLargeDefaultFont().getMaxDescent() ),
                                 g );
            g.setColor( getTreeColorSet().getLostCharactersColor() );
            TreePanel
                    .drawString( lost,
                                 parent_x
                                         + ( ( x - parent_x - getFontMetricsForLargeDefaultFont().stringWidth( lost ) ) / 2 ),
                                 ( node.getYcoord() + getFontMetricsForLargeDefaultFont().getMaxAscent() ),
                                 g );
        }
    }

    /**
     * Draw a box at the indicated node.
     * 
     * @param x
     * @param y
     * @param node
     * @param g
     */
    final private void paintNodeBox( final float x,
                                     final float y,
                                     final PhylogenyNode node,
                                     final Graphics2D g,
                                     final boolean to_pdf,
                                     final boolean to_graphics_file ) {
        if ( node.isCollapse() ) {
            return;
        }
        // if this node should be highlighted, do so
        if ( ( _highlight_node == node ) && !to_pdf && !to_graphics_file ) {
            g.setColor( getTreeColorSet().getFoundColor0() );
            drawOval( x - 8, y - 8, 16, 16, g );
            drawOval( x - 9, y - 8, 17, 17, g );
            drawOval( x - 9, y - 9, 18, 18, g );
        }
        if ( ( isInFoundNodes( node ) || isInCurrentExternalNodes( node ) )
                || ( getOptions().isShowDefaultNodeShapesExternal() && node.isExternal() )
                || ( getOptions().isShowDefaultNodeShapesInternal() && node.isInternal() )
                || ( getControlPanel().isUseVisualStyles() && ( ( node.getNodeData().getNodeVisualData() != null ) && ( ( node
                        .getNodeData().getNodeVisualData().getNodeColor() != null )
                        || ( node.getNodeData().getNodeVisualData().getSize() != NodeVisualData.DEFAULT_SIZE )
                        || ( node.getNodeData().getNodeVisualData().getFillType() != NodeFill.DEFAULT ) || ( node
                        .getNodeData().getNodeVisualData().getShape() != NodeShape.DEFAULT ) ) ) )
                || ( getControlPanel().isEvents() && node.isHasAssignedEvent() && ( node.getNodeData().getEvent()
                        .isDuplication()
                        || node.getNodeData().getEvent().isSpeciation() || node.getNodeData().getEvent()
                        .isSpeciationOrDuplication() ) ) ) {
            NodeVisualData vis = null;
            if ( getControlPanel().isUseVisualStyles() && ( node.getNodeData().getNodeVisualData() != null )
                    && ( !node.getNodeData().getNodeVisualData().isEmpty() ) ) {
                vis = node.getNodeData().getNodeVisualData();
            }
            float box_size = getOptions().getDefaultNodeShapeSize();
            if ( ( vis != null ) && ( vis.getSize() != NodeVisualData.DEFAULT_SIZE ) ) {
                box_size = vis.getSize();
            }
            final float half_box_size = box_size / 2.0f;
            Color outline_color = null;
            if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
                outline_color = Color.BLACK;
            }
            else if ( isInFoundNodes( node ) || isInCurrentExternalNodes( node ) ) {
                outline_color = getColorForFoundNode( node );
            }
            else if ( vis != null ) {
                if ( vis.getNodeColor() != null ) {
                    outline_color = vis.getNodeColor();
                }
                else if ( vis.getFontColor() != null ) {
                    outline_color = vis.getFontColor();
                }
            }
            else if ( getControlPanel().isEvents() && TreePanelUtil.isHasAssignedEvent( node ) ) {
                final Event event = node.getNodeData().getEvent();
                if ( event.isDuplication() ) {
                    outline_color = getTreeColorSet().getDuplicationBoxColor();
                }
                else if ( event.isSpeciation() ) {
                    outline_color = getTreeColorSet().getSpecBoxColor();
                }
                else if ( event.isSpeciationOrDuplication() ) {
                    outline_color = getTreeColorSet().getDuplicationOrSpeciationColor();
                }
            }
            if ( outline_color == null ) {
                outline_color = getGraphicsForNodeBoxWithColorForParentBranch( node );
                if ( to_pdf && ( outline_color == getTreeColorSet().getBranchColor() ) ) {
                    outline_color = getTreeColorSet().getBranchColorForPdf();
                }
            }
            NodeShape shape = null;
            if ( vis != null ) {
                if ( vis.getShape() == NodeShape.CIRCLE ) {
                    shape = NodeShape.CIRCLE;
                }
                else if ( vis.getShape() == NodeShape.RECTANGLE ) {
                    shape = NodeShape.RECTANGLE;
                }
            }
            if ( shape == null ) {
                if ( getOptions().getDefaultNodeShape() == NodeShape.CIRCLE ) {
                    shape = NodeShape.CIRCLE;
                }
                else if ( getOptions().getDefaultNodeShape() == NodeShape.RECTANGLE ) {
                    shape = NodeShape.RECTANGLE;
                }
            }
            NodeFill fill = null;
            if ( vis != null ) {
                if ( vis.getFillType() == NodeFill.SOLID ) {
                    fill = NodeFill.SOLID;
                }
                else if ( vis.getFillType() == NodeFill.NONE ) {
                    fill = NodeFill.NONE;
                }
                else if ( vis.getFillType() == NodeFill.GRADIENT ) {
                    fill = NodeFill.GRADIENT;
                }
            }
            if ( fill == null ) {
                if ( getOptions().getDefaultNodeFill() == NodeFill.SOLID ) {
                    fill = NodeFill.SOLID;
                }
                else if ( getOptions().getDefaultNodeFill() == NodeFill.NONE ) {
                    fill = NodeFill.NONE;
                }
                else if ( getOptions().getDefaultNodeFill() == NodeFill.GRADIENT ) {
                    fill = NodeFill.GRADIENT;
                }
            }
            Color vis_fill_color = null;
            if ( ( vis != null ) && ( vis.getNodeColor() != null ) ) {
                vis_fill_color = vis.getNodeColor();
            }
            if ( shape == NodeShape.CIRCLE ) {
                if ( fill == NodeFill.GRADIENT ) {
                    drawOvalGradient( x - half_box_size, y - half_box_size, box_size, box_size, g, to_pdf ? Color.WHITE
                            : outline_color, to_pdf ? outline_color : getBackground(), outline_color );
                }
                else if ( fill == NodeFill.NONE ) {
                    Color background = getBackground();
                    if ( to_pdf ) {
                        background = Color.WHITE;
                    }
                    drawOvalGradient( x - half_box_size,
                                      y - half_box_size,
                                      box_size,
                                      box_size,
                                      g,
                                      background,
                                      background,
                                      outline_color );
                }
                else if ( fill == NodeVisualData.NodeFill.SOLID ) {
                    if ( vis_fill_color != null ) {
                        g.setColor( vis_fill_color );
                    }
                    else {
                        g.setColor( outline_color );
                    }
                    drawOvalFilled( x - half_box_size, y - half_box_size, box_size, box_size, g );
                }
            }
            else if ( shape == NodeVisualData.NodeShape.RECTANGLE ) {
                if ( fill == NodeVisualData.NodeFill.GRADIENT ) {
                    drawRectGradient( x - half_box_size, y - half_box_size, box_size, box_size, g, to_pdf ? Color.WHITE
                            : outline_color, to_pdf ? outline_color : getBackground(), outline_color );
                }
                else if ( fill == NodeVisualData.NodeFill.NONE ) {
                    Color background = getBackground();
                    if ( to_pdf ) {
                        background = Color.WHITE;
                    }
                    drawRectGradient( x - half_box_size,
                                      y - half_box_size,
                                      box_size,
                                      box_size,
                                      g,
                                      background,
                                      background,
                                      outline_color );
                }
                else if ( fill == NodeVisualData.NodeFill.SOLID ) {
                    if ( vis_fill_color != null ) {
                        g.setColor( vis_fill_color );
                    }
                    else {
                        g.setColor( outline_color );
                    }
                    drawRectFilled( x - half_box_size, y - half_box_size, box_size, box_size, g );
                }
            }
        }
    }

    final private int paintNodeData( final Graphics2D g,
                                     final PhylogenyNode node,
                                     final boolean to_graphics_file,
                                     final boolean to_pdf,
                                     final boolean is_in_found_nodes ) {
        if ( isNodeDataInvisible( node ) && !to_graphics_file && !to_pdf ) {
            return 0;
        }
        if ( getOptions().isShowBranchLengthValues()
                && ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR )
                        || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) )
                && ( !node.isRoot() ) && ( node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) ) {
            paintBranchLength( g, node, to_pdf, to_graphics_file );
        }
        if ( !getControlPanel().isShowInternalData() && !node.isExternal() && !node.isCollapse() ) {
            return 0;
        }
        _sb.setLength( 0 );
        int x = 0;
        final int half_box_size = getOptions().getDefaultNodeShapeSize() / 2;
        if ( getControlPanel().isShowTaxonomyImages()
                && ( getImageMap() != null )
                && !getImageMap().isEmpty()
                && node.getNodeData().isHasTaxonomy()
                && ( ( node.getNodeData().getTaxonomy().getUris() != null ) && !node.getNodeData().getTaxonomy()
                        .getUris().isEmpty() ) ) {
            x += drawTaxonomyImage( node.getXcoord() + 2 + half_box_size, node.getYcoord(), node, g );
        }
        if ( ( getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyScientificNames() || getControlPanel()
                .isShowTaxonomyCommonNames() ) && node.getNodeData().isHasTaxonomy() ) {
            x += paintTaxonomy( g, node, is_in_found_nodes, to_pdf, to_graphics_file, x );
        }
        setColor( g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getSequenceColor() );
        if ( node.isCollapse() && ( ( !node.isRoot() && !node.getParent().isCollapse() ) || node.isRoot() ) ) {
            if ( _sb.length() > 0 ) {
                _sb.setLength( 0 );
                _sb.append( " (" );
                _sb.append( node.getAllExternalDescendants().size() );
                _sb.append( ")" );
            }
        }
        else {
            _sb.setLength( 0 );
        }
        nodeDataAsSB( node, _sb );
        final boolean using_visual_font = setFont( g, node, is_in_found_nodes );
        float down_shift_factor = 3.0f;
        if ( !node.isExternal() && ( node.getNumberOfDescendants() == 1 ) ) {
            down_shift_factor = 1;
        }
        final float pos_x = node.getXcoord() + x + 2 + half_box_size;
        float pos_y;
        if ( !using_visual_font ) {
            pos_y = ( node.getYcoord() + ( getFontMetricsForLargeDefaultFont().getAscent() / down_shift_factor ) );
        }
        else {
            pos_y = ( node.getYcoord() + ( getFontMetrics( g.getFont() ).getAscent() / down_shift_factor ) );
        }
        final String sb_str = _sb.toString();
        // GUILHEM_BEG ______________
        if ( _control_panel.isShowSequenceRelations() && node.getNodeData().isHasSequence()
                && ( _query_sequence != null ) ) {
            int nodeTextBoundsWidth = 0;
            if ( sb_str.length() > 0 ) {
                final Rectangle2D node_text_bounds = new TextLayout( sb_str, g.getFont(), _frc ).getBounds(); //would like to remove this 'new', but how...
                nodeTextBoundsWidth = ( int ) node_text_bounds.getWidth();
            }
            if ( node.getNodeData().getSequence().equals( _query_sequence ) ) {
                if ( nodeTextBoundsWidth > 0 ) { // invert font color and background color to show that this is the query sequence
                    g.fillRect( ( int ) pos_x - 1, ( int ) pos_y - 8, nodeTextBoundsWidth + 5, 11 );
                    g.setColor( getTreeColorSet().getBackgroundColor() );
                }
            }
            else {
                final List<SequenceRelation> seqRelations = node.getNodeData().getSequence().getSequenceRelations();
                for( final SequenceRelation seqRelation : seqRelations ) {
                    final boolean fGotRelationWithQuery = ( seqRelation.getRef0().isEqual( _query_sequence ) || seqRelation
                            .getRef1().isEqual( _query_sequence ) )
                            && seqRelation.getType().equals( getControlPanel().getSequenceRelationTypeBox()
                                    .getSelectedItem() );
                    if ( fGotRelationWithQuery ) { // we will underline the text to show that this sequence is ortholog to the query
                        final double linePosX = node.getXcoord() + 2 + half_box_size;
                        final String sConfidence = ( !getControlPanel().isShowSequenceRelationConfidence() || ( seqRelation
                                .getConfidence() == null ) ) ? null : " (" + seqRelation.getConfidence().getValue()
                                + ")";
                        if ( sConfidence != null ) {
                            float confidenceX = pos_x;
                            if ( sb_str.length() > 0 ) {
                                confidenceX += new TextLayout( sb_str, g.getFont(), _frc ).getBounds().getWidth()
                                        + CONFIDENCE_LEFT_MARGIN;
                            }
                            if ( confidenceX > linePosX ) { // let's only display confidence value if we are already displaying at least one of Prot/Gene Name and Taxonomy Code 
                                final int confidenceWidth = ( int ) new TextLayout( sConfidence, g.getFont(), _frc )
                                        .getBounds().getWidth();
                                TreePanel.drawString( sConfidence, confidenceX, pos_y, g );
                                x += CONFIDENCE_LEFT_MARGIN + confidenceWidth;
                            }
                        }
                        if ( ( x + nodeTextBoundsWidth ) > 0 ) /* we only underline if there is something displayed */
                        {
                            if ( nodeTextBoundsWidth == 0 ) {
                                nodeTextBoundsWidth -= 3; /* the gap between taxonomy code and node name should not be underlined if nothing comes after it */
                            }
                            else {
                                nodeTextBoundsWidth += 2;
                            }
                            g.drawLine( ( int ) linePosX + 1, 3 + ( int ) pos_y, ( int ) linePosX + x
                                    + nodeTextBoundsWidth, 3 + ( int ) pos_y );
                            break;
                        }
                    }
                }
            }
        }
        if ( sb_str.length() > 0 ) {
            TreePanel.drawString( sb_str, pos_x, pos_y, g );
        }
        // GUILHEM_END _____________
        if ( _sb.length() > 0 ) {
            if ( !using_visual_font && !is_in_found_nodes ) {
                x += getFontMetricsForLargeDefaultFont().stringWidth( _sb.toString() ) + 5;
            }
            else {
                x += getFontMetrics( g.getFont() ).stringWidth( _sb.toString() ) + 5;
            }
        }
        if ( getControlPanel().isShowAnnotation() && node.getNodeData().isHasSequence()
                && ( node.getNodeData().getSequence().getAnnotations() != null )
                && ( !node.getNodeData().getSequence().getAnnotations().isEmpty() ) ) {
            final SortedSet<Annotation> ann = node.getNodeData().getSequence().getAnnotations();
            if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
                g.setColor( Color.BLACK );
            }
            else if ( getControlPanel().isColorAccordingToAnnotation() ) {
                g.setColor( calculateColorForAnnotation( ann ) );
            }
            final String ann_str = TreePanelUtil.createAnnotationString( ann, getOptions().isShowAnnotationRefSource() );
            TreePanel.drawString( ann_str, node.getXcoord() + x + 3 + half_box_size, node.getYcoord()
                    + ( getFontMetricsForLargeDefaultFont().getAscent() / down_shift_factor ), g );
            _sb.setLength( 0 );
            _sb.append( ann_str );
            if ( _sb.length() > 0 ) {
                if ( !using_visual_font && !is_in_found_nodes ) {
                    x += getFontMetricsForLargeDefaultFont().stringWidth( _sb.toString() ) + 5;
                }
                else {
                    x += getFontMetrics( g.getFont() ).stringWidth( _sb.toString() ) + 5;
                }
            }
        }
        if ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR )
                || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE )
                || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) ) {
            if ( ( getControlPanel().isShowBinaryCharacters() || getControlPanel().isShowBinaryCharacterCounts() )
                    && node.getNodeData().isHasBinaryCharacters() ) {
                if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
                    g.setColor( Color.BLACK );
                }
                else {
                    g.setColor( getTreeColorSet().getBinaryDomainCombinationsColor() );
                }
                if ( getControlPanel().isShowBinaryCharacters() ) {
                    TreePanel.drawString( node.getNodeData().getBinaryCharacters().getPresentCharactersAsStringBuffer()
                            .toString(), node.getXcoord() + x + 1 + half_box_size, node.getYcoord()
                            + ( getFontMetricsForLargeDefaultFont().getAscent() / down_shift_factor ), g );
                    paintGainedAndLostCharacters( g, node, node.getNodeData().getBinaryCharacters()
                            .getGainedCharactersAsStringBuffer().toString(), node.getNodeData().getBinaryCharacters()
                            .getLostCharactersAsStringBuffer().toString() );
                }
                else {
                    TreePanel
                            .drawString( " " + node.getNodeData().getBinaryCharacters().getPresentCount(),
                                         node.getXcoord() + x + 4 + half_box_size,
                                         node.getYcoord()
                                                 + ( getFontMetricsForLargeDefaultFont().getAscent() / down_shift_factor ),
                                         g );
                    paintGainedAndLostCharacters( g, node, "+"
                            + node.getNodeData().getBinaryCharacters().getGainedCount(), "-"
                            + node.getNodeData().getBinaryCharacters().getLostCount() );
                }
            }
        }
        return x;
    }

    final private void paintNodeDataUnrootedCirc( final Graphics2D g,
                                                  final PhylogenyNode node,
                                                  final boolean to_pdf,
                                                  final boolean to_graphics_file,
                                                  final boolean radial_labels,
                                                  final double ur_angle,
                                                  final boolean is_in_found_nodes ) {
        if ( isNodeDataInvisibleUnrootedCirc( node ) && !to_graphics_file && !to_pdf ) {
            return;
        }
        _sb.setLength( 0 );
        _sb.append( " " );
        if ( node.getNodeData().isHasTaxonomy()
                && ( getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyScientificNames() || getControlPanel()
                        .isShowTaxonomyCommonNames() ) ) {
            final Taxonomy taxonomy = node.getNodeData().getTaxonomy();
            if ( _control_panel.isShowTaxonomyCode() && !ForesterUtil.isEmpty( taxonomy.getTaxonomyCode() ) ) {
                _sb.append( taxonomy.getTaxonomyCode() );
                _sb.append( " " );
            }
            if ( _control_panel.isShowTaxonomyScientificNames() && _control_panel.isShowTaxonomyCommonNames() ) {
                if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() )
                        && !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                    _sb.append( taxonomy.getScientificName() );
                    _sb.append( " (" );
                    _sb.append( taxonomy.getCommonName() );
                    _sb.append( ") " );
                }
                else if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() ) ) {
                    _sb.append( taxonomy.getScientificName() );
                    _sb.append( " " );
                }
                else if ( !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                    _sb.append( taxonomy.getCommonName() );
                    _sb.append( " " );
                }
            }
            else if ( _control_panel.isShowTaxonomyScientificNames() ) {
                if ( !ForesterUtil.isEmpty( taxonomy.getScientificName() ) ) {
                    _sb.append( taxonomy.getScientificName() );
                    _sb.append( " " );
                }
            }
            else if ( _control_panel.isShowTaxonomyCommonNames() ) {
                if ( !ForesterUtil.isEmpty( taxonomy.getCommonName() ) ) {
                    _sb.append( taxonomy.getCommonName() );
                    _sb.append( " " );
                }
            }
        }
        if ( node.isCollapse() && ( ( !node.isRoot() && !node.getParent().isCollapse() ) || node.isRoot() ) ) {
            _sb.append( " [" );
            _sb.append( node.getAllExternalDescendants().size() );
            _sb.append( "]" );
        }
        if ( getControlPanel().isShowNodeNames() && ( node.getName().length() > 0 ) ) {
            if ( _sb.length() > 0 ) {
                _sb.append( " " );
            }
            _sb.append( node.getName() );
        }
        if ( node.getNodeData().isHasSequence() ) {
            if ( getControlPanel().isShowSequenceAcc() && ( node.getNodeData().getSequence().getAccession() != null ) ) {
                if ( _sb.length() > 0 ) {
                    _sb.append( " " );
                }
                if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() ) ) {
                    _sb.append( node.getNodeData().getSequence().getAccession().getSource() );
                    _sb.append( ":" );
                }
                _sb.append( node.getNodeData().getSequence().getAccession().getValue() );
            }
            if ( getControlPanel().isShowSeqNames() && ( node.getNodeData().getSequence().getName().length() > 0 ) ) {
                if ( _sb.length() > 0 ) {
                    _sb.append( " " );
                }
                _sb.append( node.getNodeData().getSequence().getName() );
            }
        }
        //g.setFont( getTreeFontSet().getLargeFont() );
        //if ( is_in_found_nodes ) {
        //    g.setFont( getTreeFontSet().getLargeFont().deriveFont( Font.BOLD ) );
        // }
        if ( _sb.length() > 1 ) {
            setColor( g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getSequenceColor() );
            final boolean using_visual_font = setFont( g, node, is_in_found_nodes );
            final String sb_str = _sb.toString();
            double m = 0;
            if ( _graphics_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                m = _urt_nodeid_angle_map.get( node.getId() ) % TWO_PI;
            }
            else {
                m = ( float ) ( ur_angle % TWO_PI );
            }
            _at = g.getTransform();
            boolean need_to_reset = false;
            final float x_coord = node.getXcoord();
            float y_coord;
            if ( !using_visual_font ) {
                y_coord = node.getYcoord() + ( getFontMetricsForLargeDefaultFont().getAscent() / 3.0f );
            }
            else {
                y_coord = node.getYcoord() + ( getFontMetrics( g.getFont() ).getAscent() / 3.0f );
            }
            if ( radial_labels ) {
                need_to_reset = true;
                boolean left = false;
                if ( ( m > HALF_PI ) && ( m < ONEHALF_PI ) ) {
                    m -= PI;
                    left = true;
                }
                g.rotate( m, x_coord, node.getYcoord() );
                if ( left ) {
                    if ( !using_visual_font ) {
                        g.translate( -( getFontMetricsForLargeDefaultFont().getStringBounds( sb_str, g ).getWidth() ),
                                     0 );
                    }
                    else {
                        g.translate( -( getFontMetrics( g.getFont() ).getStringBounds( sb_str, g ).getWidth() ), 0 );
                    }
                }
            }
            else {
                if ( ( m > HALF_PI ) && ( m < ONEHALF_PI ) ) {
                    need_to_reset = true;
                    if ( !using_visual_font ) {
                        g.translate( -getFontMetricsForLargeDefaultFont().getStringBounds( sb_str, g ).getWidth(), 0 );
                    }
                    else {
                        g.translate( -getFontMetrics( g.getFont() ).getStringBounds( sb_str, g ).getWidth(), 0 );
                    }
                }
            }
            TreePanel.drawString( sb_str, x_coord, y_coord, g );
            if ( need_to_reset ) {
                g.setTransform( _at );
            }
        }
    }

    final private void paintNodeLite( final Graphics2D g, final PhylogenyNode node ) {
        if ( node.isCollapse() ) {
            if ( !node.isRoot() && !node.getParent().isCollapse() ) {
                paintCollapsedNode( g, node, false, false, false );
            }
            return;
        }
        if ( isInFoundNodes( node ) || isInCurrentExternalNodes( node ) ) {
            g.setColor( getColorForFoundNode( node ) );
            drawRectFilled( node.getXSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF, node.getYSecondary()
                    - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF, OVERVIEW_FOUND_NODE_BOX_SIZE, OVERVIEW_FOUND_NODE_BOX_SIZE, g );
        }
        float new_x = 0;
        if ( !node.isExternal() && !node.isCollapse() ) {
            boolean first_child = true;
            float y2 = 0.0f;
            final int parent_max_branch_to_leaf = getMaxBranchesToLeaf( node );
            for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
                final PhylogenyNode child_node = node.getChildNode( i );
                int factor_x;
                if ( !isUniformBranchLengthsForCladogram() ) {
                    factor_x = node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes();
                }
                else {
                    factor_x = parent_max_branch_to_leaf - getMaxBranchesToLeaf( child_node );
                }
                if ( first_child ) {
                    first_child = false;
                    y2 = node.getYSecondary()
                            - ( getOvYDistance() * ( node.getNumberOfExternalNodes() - child_node
                                    .getNumberOfExternalNodes() ) );
                }
                else {
                    y2 += getOvYDistance() * child_node.getNumberOfExternalNodes();
                }
                final float x2 = calculateOvBranchLengthToParent( child_node, factor_x );
                new_x = x2 + node.getXSecondary();
                final float diff_y = node.getYSecondary() - y2;
                final float diff_x = node.getXSecondary() - new_x;
                if ( ( diff_y > 2 ) || ( diff_y < -2 ) || ( diff_x > 2 ) || ( diff_x < -2 ) ) {
                    paintBranchLite( g, node.getXSecondary(), new_x, node.getYSecondary(), y2, child_node );
                }
                child_node.setXSecondary( new_x );
                child_node.setYSecondary( y2 );
                y2 += getOvYDistance() * child_node.getNumberOfExternalNodes();
            }
        }
    }

    final private void paintNodeRectangular( final Graphics2D g,
                                             final PhylogenyNode node,
                                             final boolean to_pdf,
                                             final boolean dynamically_hide,
                                             final int dynamic_hiding_factor,
                                             final boolean to_graphics_file ) {
        final boolean is_in_found_nodes = isInFoundNodes( node ) || isInCurrentExternalNodes( node );
        if ( node.isCollapse() ) {
            if ( ( !node.isRoot() && !node.getParent().isCollapse() ) ) {
                paintCollapsedNode( g, node, to_graphics_file, to_pdf, is_in_found_nodes );
            }
            return;
        }
        if ( node.isExternal() ) {
            ++_external_node_index;
        }
        // Confidence values
        if ( getControlPanel().isShowConfidenceValues()
                && !node.isExternal()
                && !node.isRoot()
                && ( ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED )
                        || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) || ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) )
                && node.getBranchData().isHasConfidences() ) {
            paintConfidenceValues( g, node, to_pdf, to_graphics_file );
        }
        // Draw a line to root:
        if ( node.isRoot() && _phylogeny.isRooted() ) {
            paintRootBranch( g, node.getXcoord(), node.getYcoord(), node, to_pdf, to_graphics_file );
        }
        float new_x = 0;
        float new_x_min = Float.MAX_VALUE;
        final boolean disallow_shortcutting = dynamic_hiding_factor < 40;
        float min_dist = 1.5f;
        if ( !disallow_shortcutting ) {
            //   System.out.println( dynamic_hiding_factor );
            if ( dynamic_hiding_factor > 4000 ) {
                min_dist = 4;
            }
            else if ( dynamic_hiding_factor > 1000 ) {
                min_dist = 3;
            }
            else if ( dynamic_hiding_factor > 100 ) {
                min_dist = 2;
            }
        }
        if ( !node.isExternal() && !node.isCollapse() ) {
            boolean first_child = true;
            float y2 = 0.0f;
            final int parent_max_branch_to_leaf = getMaxBranchesToLeaf( node );
            for( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
                final PhylogenyNode child_node = node.getChildNode( i );
                int factor_x;
                if ( !isUniformBranchLengthsForCladogram() ) {
                    factor_x = node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes();
                }
                else {
                    factor_x = parent_max_branch_to_leaf - getMaxBranchesToLeaf( child_node );
                }
                if ( first_child ) {
                    first_child = false;
                    y2 = node.getYcoord()
                            - ( _y_distance * ( node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes() ) );
                }
                else {
                    y2 += _y_distance * child_node.getNumberOfExternalNodes();
                }
                final float x2 = calculateBranchLengthToParent( child_node, factor_x );
                new_x = x2 + node.getXcoord();
                if ( dynamically_hide && ( x2 < new_x_min ) ) {
                    new_x_min = x2;
                }
                final float diff_y = node.getYcoord() - y2;
                final float diff_x = node.getXcoord() - new_x;
                if ( disallow_shortcutting || ( diff_y > min_dist ) || ( diff_y < -min_dist ) || ( diff_x > min_dist )
                        || ( diff_x < -min_dist ) || to_graphics_file || to_pdf ) {
                    paintBranchRectangular( g,
                                            node.getXcoord(),
                                            new_x,
                                            node.getYcoord(),
                                            y2,
                                            child_node,
                                            to_pdf,
                                            to_graphics_file );
                }
                child_node.setXcoord( new_x );
                child_node.setYcoord( y2 );
                y2 += _y_distance * child_node.getNumberOfExternalNodes();
            }
            paintNodeBox( node.getXcoord(), node.getYcoord(), node, g, to_pdf, to_graphics_file );
        }
        if ( dynamically_hide
                && !is_in_found_nodes
                && ( ( node.isExternal() && ( ( _external_node_index % dynamic_hiding_factor ) != 1 ) ) || ( !node
                        .isExternal() && ( ( new_x_min < 20 ) || ( ( _y_distance * node.getNumberOfExternalNodes() ) < getFontMetricsForLargeDefaultFont()
                        .getHeight() ) ) ) ) ) {
            return;
        }
        final int x = paintNodeData( g, node, to_graphics_file, to_pdf, is_in_found_nodes );
        paintNodeWithRenderableData( x, g, node, to_graphics_file, to_pdf );
    }

    final private void paintNodeWithRenderableData( final int x,
                                                    final Graphics2D g,
                                                    final PhylogenyNode node,
                                                    final boolean to_graphics_file,
                                                    final boolean to_pdf ) {
        if ( isNodeDataInvisible( node ) && !to_graphics_file ) {
            return;
        }
        if ( ( !getControlPanel().isShowInternalData() && !node.isExternal() ) ) {
            return;
        }
        if ( getControlPanel().isShowDomainArchitectures() && node.getNodeData().isHasSequence()
                && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
            RenderableDomainArchitecture rds = null;
            if ( node.getNodeData().getSequence().getDomainArchitecture() instanceof RenderableDomainArchitecture ) {
                try {
                    rds = ( RenderableDomainArchitecture ) node.getNodeData().getSequence().getDomainArchitecture();
                }
                catch ( final ClassCastException cce ) {
                    cce.printStackTrace();
                }
                if ( rds != null ) {
                    rds.setRenderingHeight( 6 );
                    rds.render( node.getXcoord() + x, node.getYcoord() - 3, g, this, to_pdf );
                }
            }
        }
        //////////////
        if ( getControlPanel().isShowVectorData() && ( node.getNodeData().getVector() != null )
                && ( node.getNodeData().getVector().size() > 0 ) && ( getStatisticsForExpressionValues() != null ) ) {
            final RenderableVector rv = RenderableVector.createInstance( node.getNodeData().getVector(),
                                                                         getStatisticsForExpressionValues(),
                                                                         getConfiguration() );
            if ( rv != null ) {
                int xx = 0;
                PhylogenyNode my_node = node;
                if ( !getControlPanel().isDrawPhylogram() ) {
                    my_node = getPhylogeny().getFirstExternalNode();
                }
                if ( getControlPanel().isShowTaxonomyCode() && ( PhylogenyMethods.getSpecies( my_node ).length() > 0 ) ) {
                    xx += getFontMetricsForLargeDefaultFont()
                            .stringWidth( PhylogenyMethods.getSpecies( my_node ) + " " );
                }
                if ( getControlPanel().isShowNodeNames() && ( my_node.getName().length() > 0 ) ) {
                    xx += getFontMetricsForLargeDefaultFont().stringWidth( my_node.getName() + " " );
                }
                rv.render( my_node.getXcoord() + xx, node.getYcoord() - 5, g, this, to_pdf );
            }
        }
        //////////////
    }

    final private void paintOvRectangle( final Graphics2D g ) {
        final float w_ratio = ( ( float ) getWidth() ) / getVisibleRect().width;
        final float h_ratio = ( ( float ) getHeight() ) / getVisibleRect().height;
        final float x_ratio = ( ( float ) getWidth() ) / getVisibleRect().x;
        final float y_ratio = ( ( float ) getHeight() ) / getVisibleRect().y;
        final float width = getOvMaxWidth() / w_ratio;
        final float height = getOvMaxHeight() / h_ratio;
        final float x = getVisibleRect().x + getOvXPosition() + ( getOvMaxWidth() / x_ratio );
        final float y = getVisibleRect().y + getOvYPosition() + ( getOvMaxHeight() / y_ratio );
        g.setColor( getTreeColorSet().getFoundColor0() );
        getOvRectangle().setRect( x, y, width, height );
        final Stroke s = g.getStroke();
        g.setStroke( STROKE_1 );
        if ( ( width < 6 ) && ( height < 6 ) ) {
            drawRectFilled( x, y, 6, 6, g );
            getOvVirtualRectangle().setRect( x, y, 6, 6 );
        }
        else if ( width < 6 ) {
            drawRectFilled( x, y, 6, height, g );
            getOvVirtualRectangle().setRect( x, y, 6, height );
        }
        else if ( height < 6 ) {
            drawRectFilled( x, y, width, 6, g );
            getOvVirtualRectangle().setRect( x, y, width, 6 );
        }
        else {
            drawRect( x, y, width, height, g );
            if ( isInOvRect() ) {
                drawRect( x + 1, y + 1, width - 2, height - 2, g );
            }
            getOvVirtualRectangle().setRect( x, y, width, height );
        }
        g.setStroke( s );
    }

    final private void paintPhylogenyLite( final Graphics2D g ) {
        _phylogeny
                .getRoot()
                .setXSecondary( ( float ) ( getVisibleRect().x + getOvXPosition() + ( MOVE / ( getVisibleRect().width / getOvRectangle()
                        .getWidth() ) ) ) );
        _phylogeny.getRoot().setYSecondary( ( getVisibleRect().y + getOvYStart() ) );
        final Stroke s = g.getStroke();
        g.setStroke( STROKE_05 );
        for( final PhylogenyNode element : _nodes_in_preorder ) {
            paintNodeLite( g, element );
        }
        g.setStroke( s );
        paintOvRectangle( g );
    }

    /**
     * Paint the root branch. (Differs from others because it will always be a
     * single horizontal line).
     * @param to_graphics_file 
     * 
     * @return new x1 value
     */
    final private void paintRootBranch( final Graphics2D g,
                                        final float x1,
                                        final float y1,
                                        final PhylogenyNode root,
                                        final boolean to_pdf,
                                        final boolean to_graphics_file ) {
        assignGraphicsForBranchWithColorForParentBranch( root, false, g, to_pdf, to_graphics_file );
        float d = getXdistance();
        if ( getControlPanel().isDrawPhylogram() && ( root.getDistanceToParent() > 0.0 ) ) {
            d = ( float ) ( getXcorrectionFactor() * root.getDistanceToParent() );
        }
        if ( d < MIN_ROOT_LENGTH ) {
            d = MIN_ROOT_LENGTH;
        }
        if ( !getControlPanel().isWidthBranches() || ( PhylogenyMethods.getBranchWidthValue( root ) == 1 ) ) {
            drawLine( x1 - d, root.getYcoord(), x1, root.getYcoord(), g );
        }
        else {
            final double w = PhylogenyMethods.getBranchWidthValue( root );
            drawRectFilled( x1 - d, root.getYcoord() - ( w / 2 ), d, w, g );
        }
        paintNodeBox( x1, root.getYcoord(), root, g, to_pdf, to_graphics_file );
    }

    final private void paintScale( final Graphics2D g,
                                   int x1,
                                   int y1,
                                   final boolean to_pdf,
                                   final boolean to_graphics_file ) {
        x1 += MOVE;
        final double x2 = x1 + ( getScaleDistance() * getXcorrectionFactor() );
        y1 -= 12;
        final int y2 = y1 - 8;
        final int y3 = y1 - 4;
        g.setFont( getTreeFontSet().getSmallFont() );
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else {
            g.setColor( getTreeColorSet().getBranchLengthColor() );
        }
        final Stroke s = g.getStroke();
        g.setStroke( STROKE_1 );
        drawLine( x1, y1, x1, y2, g );
        drawLine( x2, y1, x2, y2, g );
        drawLine( x1, y3, x2, y3, g );
        if ( getScaleLabel() != null ) {
            g.drawString( getScaleLabel(), ( x1 + 2 ), y3 - 2 );
        }
        g.setStroke( s );
    }

    final private int paintTaxonomy( final Graphics2D g,
                                     final PhylogenyNode node,
                                     final boolean is_in_found_nodes,
                                     final boolean to_pdf,
                                     final boolean to_graphics_file,
                                     final float x_shift ) {
        final Taxonomy taxonomy = node.getNodeData().getTaxonomy();
        final boolean using_visual_font = setFont( g, node, is_in_found_nodes );
        setColor( g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getTaxonomyColor() );
        final float start_x = node.getXcoord() + 3 + ( getOptions().getDefaultNodeShapeSize() / 2 ) + x_shift;
        float start_y;
        if ( !using_visual_font ) {
            start_y = node.getYcoord()
                    + ( getFontMetricsForLargeDefaultFont().getAscent() / ( node.getNumberOfDescendants() == 1 ? 1
                            : 3.0f ) );
        }
        else {
            start_y = node.getYcoord()
                    + ( getFontMetrics( g.getFont() ).getAscent() / ( node.getNumberOfDescendants() == 1 ? 1 : 3.0f ) );
        }
        _sb.setLength( 0 );
        nodeTaxonomyDataAsSB( taxonomy, _sb );
        final String label = _sb.toString();
        /* GUILHEM_BEG */
        if ( _control_panel.isShowSequenceRelations() && ( label.length() > 0 )
                && ( node.getNodeData().isHasSequence() ) && node.getNodeData().getSequence().equals( _query_sequence ) ) {
            // invert font color and background color to show that this is the query sequence
            final Rectangle2D nodeTextBounds = new TextLayout( label, g.getFont(), new FontRenderContext( null,
                                                                                                          false,
                                                                                                          false ) )
                    .getBounds();
            g.fillRect( ( int ) start_x - 1, ( int ) start_y - 8, ( int ) nodeTextBounds.getWidth() + 4, 11 );
            g.setColor( getTreeColorSet().getBackgroundColor() );
        }
        /* GUILHEM_END */
        TreePanel.drawString( label, start_x, start_y, g );
        if ( !using_visual_font && !is_in_found_nodes ) {
            return getFontMetricsForLargeDefaultFont().stringWidth( label );
        }
        return getFontMetrics( g.getFont() ).stringWidth( label );
    }

    final private void paintUnrooted( final PhylogenyNode n,
                                      final double low_angle,
                                      final double high_angle,
                                      final boolean radial_labels,
                                      final Graphics2D g,
                                      final boolean to_pdf,
                                      final boolean to_graphics_file ) {
        if ( n.isRoot() ) {
            n.setXcoord( getWidth() / 2 );
            n.setYcoord( getHeight() / 2 );
        }
        if ( n.isExternal() ) {
            paintNodeDataUnrootedCirc( g,
                                       n,
                                       to_pdf,
                                       to_graphics_file,
                                       radial_labels,
                                       ( high_angle + low_angle ) / 2,
                                       isInFoundNodes( n ) || isInCurrentExternalNodes( n ) );
            return;
        }
        final float num_enclosed = n.getNumberOfExternalNodes();
        final float x = n.getXcoord();
        final float y = n.getYcoord();
        double current_angle = low_angle;
        // final boolean n_below = n.getYcoord() < getVisibleRect().getMinY() - 20;
        // final boolean n_above = n.getYcoord() > getVisibleRect().getMaxY() + 20;
        // final boolean n_left = n.getXcoord() < getVisibleRect().getMinX() - 20;
        // final boolean n_right = n.getXcoord() > getVisibleRect().getMaxX() + 20;
        for( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
            final PhylogenyNode desc = n.getChildNode( i );
            ///  if ( ( ( n_below ) & ( desc.getYcoord() < getVisibleRect().getMinY() - 20 ) )
            //          || ( ( n_above ) & ( desc.getYcoord() > getVisibleRect().getMaxY() + 20 ) )
            //         || ( ( n_left ) & ( desc.getXcoord() < getVisibleRect().getMinX() - 20 ) )
            //          || ( ( n_right ) & ( desc.getXcoord() > getVisibleRect().getMaxX() + 20 ) ) ) {
            //     continue;
            // }
            //if ( ( desc.getYcoord() > n.getYcoord() ) && ( n.getYcoord() > getVisibleRect().getMaxY() - 20 ) ) {
            //    continue;
            //}
            //if ( ( desc.getYcoord() < n.getYcoord() ) && ( n.getYcoord() < getVisibleRect().getMinY() + 20 ) ) {
            //    continue;
            // }
            final int desc_num_enclosed = desc.getNumberOfExternalNodes();
            final double arc_size = ( desc_num_enclosed / num_enclosed ) * ( high_angle - low_angle );
            float length;
            if ( isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram() ) {
                if ( desc.getDistanceToParent() < 0 ) {
                    length = 0;
                }
                else {
                    length = ( float ) ( desc.getDistanceToParent() * getUrtFactor() );
                }
            }
            else {
                length = getUrtFactor();
            }
            final double mid_angle = current_angle + ( arc_size / 2 );
            final float new_x = ( float ) ( x + ( Math.cos( mid_angle ) * length ) );
            final float new_y = ( float ) ( y + ( Math.sin( mid_angle ) * length ) );
            desc.setXcoord( new_x );
            desc.setYcoord( new_y );
            paintUnrooted( desc, current_angle, current_angle + arc_size, radial_labels, g, to_pdf, to_graphics_file );
            current_angle += arc_size;
            assignGraphicsForBranchWithColorForParentBranch( desc, false, g, to_pdf, to_graphics_file );
            drawLine( x, y, new_x, new_y, g );
            paintNodeBox( new_x, new_y, desc, g, to_pdf, to_graphics_file );
        }
        if ( n.isRoot() ) {
            paintNodeBox( n.getXcoord(), n.getYcoord(), n, g, to_pdf, to_graphics_file );
        }
    }

    final private void paintUnrootedLite( final PhylogenyNode n,
                                          final double low_angle,
                                          final double high_angle,
                                          final Graphics2D g,
                                          final float urt_ov_factor ) {
        if ( n.isRoot() ) {
            final int x_pos = ( int ) ( getVisibleRect().x + getOvXPosition() + ( getOvMaxWidth() / 2 ) );
            final int y_pos = ( int ) ( getVisibleRect().y + getOvYPosition() + ( getOvMaxHeight() / 2 ) );
            n.setXSecondary( x_pos );
            n.setYSecondary( y_pos );
        }
        if ( n.isExternal() ) {
            return;
        }
        final float num_enclosed = n.getNumberOfExternalNodes();
        final float x = n.getXSecondary();
        final float y = n.getYSecondary();
        double current_angle = low_angle;
        for( int i = 0; i < n.getNumberOfDescendants(); ++i ) {
            final PhylogenyNode desc = n.getChildNode( i );
            final int desc_num_enclosed = desc.getNumberOfExternalNodes();
            final double arc_size = ( desc_num_enclosed / num_enclosed ) * ( high_angle - low_angle );
            float length;
            if ( isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram() ) {
                if ( desc.getDistanceToParent() < 0 ) {
                    length = 0;
                }
                else {
                    length = ( float ) ( desc.getDistanceToParent() * urt_ov_factor );
                }
            }
            else {
                length = urt_ov_factor;
            }
            final double mid_angle = current_angle + ( arc_size / 2 );
            final float new_x = ( float ) ( x + ( Math.cos( mid_angle ) * length ) );
            final float new_y = ( float ) ( y + ( Math.sin( mid_angle ) * length ) );
            desc.setXSecondary( new_x );
            desc.setYSecondary( new_y );
            if ( isInFoundNodes( desc ) || isInCurrentExternalNodes( desc ) ) {
                g.setColor( getColorForFoundNode( desc ) );
                drawRectFilled( desc.getXSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                                desc.getYSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                                OVERVIEW_FOUND_NODE_BOX_SIZE,
                                OVERVIEW_FOUND_NODE_BOX_SIZE,
                                g );
                g.setColor( getTreeColorSet().getOvColor() );
            }
            paintUnrootedLite( desc, current_angle, current_angle + arc_size, g, urt_ov_factor );
            current_angle += arc_size;
            drawLine( x, y, new_x, new_y, g );
        }
    }

    final private void pasteSubtree( final PhylogenyNode node ) {
        if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if ( ( getCutOrCopiedTree() == null ) || getCutOrCopiedTree().isEmpty() ) {
            JOptionPane.showMessageDialog( this,
                                           "No tree in buffer (need to copy or cut a subtree first)",
                                           "Attempt to paste with empty buffer",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        final String label = createASimpleTextRepresentationOfANode( getCutOrCopiedTree().getRoot() );
        final Object[] options = { "As sibling", "As descendant", "Cancel" };
        final int r = JOptionPane.showOptionDialog( this,
                                                    "How to paste subtree" + label + "?",
                                                    "Paste Subtree",
                                                    JOptionPane.CLOSED_OPTION,
                                                    JOptionPane.QUESTION_MESSAGE,
                                                    null,
                                                    options,
                                                    options[ 2 ] );
        boolean paste_as_sibling = true;
        if ( r == 1 ) {
            paste_as_sibling = false;
        }
        else if ( r != 0 ) {
            return;
        }
        final Phylogeny buffer_phy = getCutOrCopiedTree().copy();
        buffer_phy.setAllNodesToNotCollapse();
        PhylogenyMethods.preOrderReId( buffer_phy );
        buffer_phy.setRooted( true );
        boolean need_to_show_whole = false;
        if ( paste_as_sibling ) {
            if ( node.isRoot() ) {
                JOptionPane.showMessageDialog( this,
                                               "Cannot paste sibling to root",
                                               "Attempt to paste sibling to root",
                                               JOptionPane.ERROR_MESSAGE );
                return;
            }
            buffer_phy.addAsSibling( node );
        }
        else {
            if ( ( node.getNumberOfExternalNodes() == 1 ) && node.isRoot() ) {
                need_to_show_whole = true;
                _phylogeny = buffer_phy;
            }
            else {
                buffer_phy.addAsChild( node );
            }
        }
        if ( getCopiedAndPastedNodes() == null ) {
            setCopiedAndPastedNodes( new HashSet<Long>() );
        }
        final List<PhylogenyNode> nodes = PhylogenyMethods.obtainAllNodesAsList( buffer_phy );
        final Set<Long> node_ids = new HashSet<Long>( nodes.size() );
        for( final PhylogenyNode n : nodes ) {
            node_ids.add( n.getId() );
        }
        node_ids.add( node.getId() );
        getCopiedAndPastedNodes().addAll( node_ids );
        setNodeInPreorderToNull();
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants( true );
        resetNodeIdToDistToLeafMap();
        setEdited( true );
        if ( need_to_show_whole ) {
            getControlPanel().showWhole();
        }
        repaint();
    }

    private final StringBuffer propertiesToString( final PhylogenyNode node ) {
        final PropertiesMap properties = node.getNodeData().getProperties();
        final StringBuffer sb = new StringBuffer();
        boolean first = true;
        for( final String ref : properties.getPropertyRefs() ) {
            if ( first ) {
                first = false;
            }
            else {
                sb.append( " " );
            }
            final Property p = properties.getProperty( ref );
            sb.append( TreePanelUtil.getPartAfterColon( p.getRef() ) );
            sb.append( "=" );
            sb.append( p.getValue() );
            if ( !ForesterUtil.isEmpty( p.getUnit() ) ) {
                sb.append( TreePanelUtil.getPartAfterColon( p.getUnit() ) );
            }
        }
        return sb;
    }

    private void setColor( final Graphics2D g,
                           final PhylogenyNode node,
                           final boolean to_graphics_file,
                           final boolean to_pdf,
                           final boolean is_in_found_nodes,
                           final Color default_color ) {
        if ( ( to_pdf || to_graphics_file ) && getOptions().isPrintBlackAndWhite() ) {
            g.setColor( Color.BLACK );
        }
        else if ( is_in_found_nodes ) {
            g.setColor( getColorForFoundNode( node ) );
        }
        else if ( getControlPanel().isUseVisualStyles() && ( node.getNodeData().getNodeVisualData() != null )
                && ( node.getNodeData().getNodeVisualData().getFontColor() != null ) ) {
            g.setColor( node.getNodeData().getNodeVisualData().getFontColor() );
        }
        else if ( getControlPanel().isColorAccordingToSequence() ) {
            g.setColor( getSequenceBasedColor( node ) );
        }
        else if ( getControlPanel().isColorAccordingToTaxonomy() ) {
            g.setColor( getTaxonomyBasedColor( node ) );
        }
        else if ( getControlPanel().isColorAccordingToAnnotation()
                && ( node.getNodeData().isHasSequence() && ( node.getNodeData().getSequence().getAnnotations() != null ) && ( !node
                        .getNodeData().getSequence().getAnnotations().isEmpty() ) ) ) {
            g.setColor( calculateColorForAnnotation( node.getNodeData().getSequence().getAnnotations() ) );
        }
        else if ( getOptions().isColorLabelsSameAsParentBranch() && getControlPanel().isUseVisualStyles()
                && ( PhylogenyMethods.getBranchColorValue( node ) != null ) ) {
            g.setColor( PhylogenyMethods.getBranchColorValue( node ) );
        }
        else if ( to_pdf ) {
            g.setColor( Color.BLACK );
        }
        else {
            g.setColor( default_color );
        }
    }

    final private void setCopiedAndPastedNodes( final Set<Long> nodeIds ) {
        getMainPanel().setCopiedAndPastedNodes( nodeIds );
    }

    final private void setCutOrCopiedTree( final Phylogeny cut_or_copied_tree ) {
        getMainPanel().setCutOrCopiedTree( cut_or_copied_tree );
    }

    private boolean setFont( final Graphics2D g, final PhylogenyNode node, final boolean is_in_found_nodes ) {
        Font visual_font = null;
        if ( getControlPanel().isUseVisualStyles() && ( node.getNodeData().getNodeVisualData() != null ) ) {
            visual_font = node.getNodeData().getNodeVisualData().getFont();
            g.setFont( visual_font != null ? visual_font : getTreeFontSet().getLargeFont() );
        }
        else {
            g.setFont( getTreeFontSet().getLargeFont() );
        }
        if ( is_in_found_nodes ) {
            g.setFont( g.getFont().deriveFont( Font.BOLD ) );
        }
        return visual_font != null;
    }

    final private void setInOv( final boolean in_ov ) {
        _in_ov = in_ov;
    }

    final private void setOvMaxHeight( final float ov_max_height ) {
        _ov_max_height = ov_max_height;
    }

    final private void setOvMaxWidth( final float ov_max_width ) {
        _ov_max_width = ov_max_width;
    }

    final private void setOvXcorrectionFactor( final float f ) {
        _ov_x_correction_factor = f;
    }

    final private void setOvXDistance( final float ov_x_distance ) {
        _ov_x_distance = ov_x_distance;
    }

    final private void setOvXPosition( final int ov_x_position ) {
        _ov_x_position = ov_x_position;
    }

    final private void setOvYDistance( final float ov_y_distance ) {
        _ov_y_distance = ov_y_distance;
    }

    final private void setOvYPosition( final int ov_y_position ) {
        _ov_y_position = ov_y_position;
    }

    final private void setOvYStart( final int ov_y_start ) {
        _ov_y_start = ov_y_start;
    }

    final private void setScaleDistance( final double scale_distance ) {
        _scale_distance = scale_distance;
    }

    final private void setScaleLabel( final String scale_label ) {
        _scale_label = scale_label;
    }

    private final void setupStroke( final Graphics2D g ) {
        if ( getYdistance() < 0.001 ) {
            g.setStroke( STROKE_005 );
        }
        else if ( getYdistance() < 0.01 ) {
            g.setStroke( STROKE_01 );
        }
        else if ( getYdistance() < 0.5 ) {
            g.setStroke( STROKE_025 );
        }
        else if ( getYdistance() < 1 ) {
            g.setStroke( STROKE_05 );
        }
        else if ( getYdistance() < 2 ) {
            g.setStroke( STROKE_075 );
        }
        else if ( getYdistance() < 20 ) {
            g.setStroke( STROKE_1 );
        }
        else {
            g.setStroke( STROKE_2 );
        }
    }

    final private void setUpUrtFactor() {
        final int d = getVisibleRect().width < getVisibleRect().height ? getVisibleRect().width
                : getVisibleRect().height;
        if ( isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram() ) {
            setUrtFactor( ( float ) ( d / ( 2 * getMaxDistanceToRoot() ) ) );
        }
        else {
            final int max_depth = _circ_max_depth;
            if ( max_depth > 0 ) {
                setUrtFactor( d / ( 2 * max_depth ) );
            }
            else {
                setUrtFactor( d / 2 );
            }
        }
        setUrtFactorOv( getUrtFactor() );
    }

    final private void setUrtFactor( final float urt_factor ) {
        _urt_factor = urt_factor;
    }

    final private void setUrtFactorOv( final float urt_factor_ov ) {
        _urt_factor_ov = urt_factor_ov;
    }

    private void showExtDescNodeData( final PhylogenyNode node ) {
        final List<String> data = new ArrayList<String>();
        final List<PhylogenyNode> nodes = node.getAllExternalDescendants();
        if ( ( getFoundNodes0() != null ) || ( getFoundNodes1() != null ) ) {
            for( final PhylogenyNode n : getFoundNodesAsListOfPhylogenyNodes() ) {
                if ( !nodes.contains( n ) ) {
                    nodes.add( n );
                }
            }
        }
        for( final PhylogenyNode n : nodes ) {
            switch ( getOptions().getExtDescNodeDataToReturn() ) {
                case NODE_NAME:
                    if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                        data.add( n.getName() );
                    }
                    break;
                case SEQUENCE_NAME:
                    if ( n.getNodeData().isHasSequence()
                            && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getName() ) ) {
                        data.add( n.getNodeData().getSequence().getName() );
                    }
                    break;
                case GENE_NAME:
                    if ( n.getNodeData().isHasSequence()
                            && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getGeneName() ) ) {
                        data.add( n.getNodeData().getSequence().getGeneName() );
                    }
                    break;
                case SEQUENCE_SYMBOL:
                    if ( n.getNodeData().isHasSequence()
                            && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getSymbol() ) ) {
                        data.add( n.getNodeData().getSequence().getSymbol() );
                    }
                    break;
                case SEQUENCE_MOL_SEQ:
                    if ( n.getNodeData().isHasSequence()
                            && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getMolecularSequence() ) ) {
                        data.add( n.getNodeData().getSequence().getMolecularSequence() );
                    }
                    break;
                case SEQUENCE_MOL_SEQ_FASTA:
                    final StringBuilder sb = new StringBuilder();
                    if ( n.getNodeData().isHasSequence()
                            && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getMolecularSequence() ) ) {
                        final StringBuilder ann = new StringBuilder();
                        if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                            ann.append( n.getName() );
                            ann.append( "|" );
                        }
                        if ( !ForesterUtil.isEmpty( n.getNodeData().getSequence().getSymbol() ) ) {
                            ann.append( "SYM=" );
                            ann.append( n.getNodeData().getSequence().getSymbol() );
                            ann.append( "|" );
                        }
                        if ( !ForesterUtil.isEmpty( n.getNodeData().getSequence().getName() ) ) {
                            ann.append( "NAME=" );
                            ann.append( n.getNodeData().getSequence().getName() );
                            ann.append( "|" );
                        }
                        if ( !ForesterUtil.isEmpty( n.getNodeData().getSequence().getGeneName() ) ) {
                            ann.append( "GN=" );
                            ann.append( n.getNodeData().getSequence().getGeneName() );
                            ann.append( "|" );
                        }
                        if ( n.getNodeData().getSequence().getAccession() != null ) {
                            ann.append( "ACC=" );
                            ann.append( n.getNodeData().getSequence().getAccession().asText() );
                            ann.append( "|" );
                        }
                        if ( n.getNodeData().isHasTaxonomy() ) {
                            if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
                                ann.append( "TAXID=" );
                                ann.append( n.getNodeData().getTaxonomy().getTaxonomyCode() );
                                ann.append( "|" );
                            }
                            if ( !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                                ann.append( "SN=" );
                                ann.append( n.getNodeData().getTaxonomy().getScientificName() );
                                ann.append( "|" );
                            }
                        }
                        String ann_str;
                        if ( ann.charAt( ann.length() - 1 ) == '|' ) {
                            ann_str = ann.substring( 0, ann.length() - 1 );
                        }
                        else {
                            ann_str = ann.toString();
                        }
                        sb.append( SequenceWriter.toFasta( ann_str, n.getNodeData().getSequence()
                                .getMolecularSequence(), 60 ) );
                        data.add( sb.toString() );
                    }
                    break;
                case SEQUENCE_ACC:
                    if ( n.getNodeData().isHasSequence() && ( n.getNodeData().getSequence().getAccession() != null )
                            && !ForesterUtil.isEmpty( n.getNodeData().getSequence().getAccession().toString() ) ) {
                        data.add( n.getNodeData().getSequence().getAccession().toString() );
                    }
                    break;
                case TAXONOMY_SCIENTIFIC_NAME:
                    if ( n.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                        data.add( n.getNodeData().getTaxonomy().getScientificName() );
                    }
                    break;
                case TAXONOMY_COMM0N_NAME:
                    if ( n.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getCommonName() ) ) {
                        data.add( n.getNodeData().getTaxonomy().getCommonName() );
                    }
                    break;
                case TAXONOMY_CODE:
                    if ( n.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
                        data.add( n.getNodeData().getTaxonomy().getTaxonomyCode() );
                    }
                    break;
                case UNKNOWN:
                    TreePanelUtil.showExtDescNodeDataUserSelectedHelper( getControlPanel(), n, data );
                    break;
                default:
                    throw new IllegalArgumentException( "unknown data element: "
                            + getOptions().getExtDescNodeDataToReturn() );
            }
        } // for loop
        final StringBuilder sb = new StringBuilder();
        final int size = TreePanelUtil.makeSB( data, getOptions(), sb );
        if ( ( getConfiguration().getExtNodeDataReturnOn() == EXT_NODE_DATA_RETURN_ON.CONSOLE )
                || ( getConfiguration().getExtNodeDataReturnOn() == EXT_NODE_DATA_RETURN_ON.BUFFER_ONLY ) ) {
            if ( getConfiguration().getExtNodeDataReturnOn() == EXT_NODE_DATA_RETURN_ON.CONSOLE ) {
                System.out.println( sb );
            }
            if ( sb.length() < 1 ) {
                clearCurrentExternalNodesDataBuffer();
            }
            else {
                setCurrentExternalNodesDataBuffer( sb );
            }
        }
        else if ( getConfiguration().getExtNodeDataReturnOn() == EXT_NODE_DATA_RETURN_ON.WINODW ) {
            if ( sb.length() < 1 ) {
                TreePanelUtil.showInformationMessage( this, "No Appropriate Data (" + obtainTitleForExtDescNodeData()
                        + ")", "Descendants of selected node do not contain selected data" );
                clearCurrentExternalNodesDataBuffer();
            }
            else {
                setCurrentExternalNodesDataBuffer( sb );
                String title;
                if ( ( getFoundNodes0() != null ) && !getFoundNodes0().isEmpty() ) {
                    title = ( getOptions().getExtDescNodeDataToReturn() == NODE_DATA.UNKNOWN ? "Data"
                            : obtainTitleForExtDescNodeData() )
                            + " for "
                            + data.size()
                            + " nodes, unique entries: "
                            + size;
                }
                else {
                    title = ( getOptions().getExtDescNodeDataToReturn() == NODE_DATA.UNKNOWN ? "Data"
                            : obtainTitleForExtDescNodeData() )
                            + " for "
                            + data.size()
                            + "/"
                            + node.getNumberOfExternalNodes()
                            + " external descendats of node "
                            + node
                            + ", unique entries: " + size;
                }
                final String s = sb.toString().trim();
                if ( getMainPanel().getMainFrame() == null ) {
                    // Must be "E" applet version.
                    final ArchaeopteryxE ae = ( ArchaeopteryxE ) ( ( MainPanelApplets ) getMainPanel() ).getApplet();
                    ae.showTextFrame( s, title );
                }
                else {
                    getMainPanel().getMainFrame().showTextFrame( s, title );
                }
            }
        }
    }

    final private void showNodeDataPopup( final MouseEvent e, final PhylogenyNode node ) {
        try {
            if ( ( node.getName().length() > 0 )
                    || ( node.getNodeData().isHasTaxonomy() && !TreePanelUtil.isTaxonomyEmpty( node.getNodeData()
                            .getTaxonomy() ) )
                    || ( node.getNodeData().isHasSequence() && !TreePanelUtil.isSequenceEmpty( node.getNodeData()
                            .getSequence() ) ) || ( node.getNodeData().isHasDate() )
                    || ( node.getNodeData().isHasDistribution() ) || node.getBranchData().isHasConfidences() ) {
                _popup_buffer.setLength( 0 );
                short lines = 0;
                if ( node.getName().length() > 0 ) {
                    lines++;
                    _popup_buffer.append( node.getName() );
                }
                if ( node.getNodeData().isHasTaxonomy()
                        && !TreePanelUtil.isTaxonomyEmpty( node.getNodeData().getTaxonomy() ) ) {
                    lines++;
                    boolean enc_data = false;
                    final Taxonomy tax = node.getNodeData().getTaxonomy();
                    if ( _popup_buffer.length() > 0 ) {
                        _popup_buffer.append( "\n" );
                    }
                    if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                        _popup_buffer.append( "[" );
                        _popup_buffer.append( tax.getTaxonomyCode() );
                        _popup_buffer.append( "]" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " " );
                        }
                        _popup_buffer.append( tax.getScientificName() );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( tax.getCommonName() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " (" );
                        }
                        else {
                            _popup_buffer.append( "(" );
                        }
                        _popup_buffer.append( tax.getCommonName() );
                        _popup_buffer.append( ")" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( tax.getAuthority() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " (" );
                        }
                        else {
                            _popup_buffer.append( "(" );
                        }
                        _popup_buffer.append( tax.getAuthority() );
                        _popup_buffer.append( ")" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( tax.getRank() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " [" );
                        }
                        else {
                            _popup_buffer.append( "[" );
                        }
                        _popup_buffer.append( tax.getRank() );
                        _popup_buffer.append( "]" );
                        enc_data = true;
                    }
                    if ( tax.getSynonyms().size() > 0 ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " " );
                        }
                        _popup_buffer.append( "[" );
                        int counter = 1;
                        for( final String syn : tax.getSynonyms() ) {
                            if ( !ForesterUtil.isEmpty( syn ) ) {
                                enc_data = true;
                                _popup_buffer.append( syn );
                                if ( counter < tax.getSynonyms().size() ) {
                                    _popup_buffer.append( ", " );
                                }
                            }
                            counter++;
                        }
                        _popup_buffer.append( "]" );
                    }
                    if ( !enc_data ) {
                        if ( ( tax.getIdentifier() != null ) && !ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) {
                            if ( !ForesterUtil.isEmpty( tax.getIdentifier().getProvider() ) ) {
                                _popup_buffer.append( "[" );
                                _popup_buffer.append( tax.getIdentifier().getProvider() );
                                _popup_buffer.append( "] " );
                            }
                            _popup_buffer.append( tax.getIdentifier().getValue() );
                        }
                    }
                }
                if ( node.getNodeData().isHasSequence()
                        && !TreePanelUtil.isSequenceEmpty( node.getNodeData().getSequence() ) ) {
                    lines++;
                    boolean enc_data = false;
                    if ( _popup_buffer.length() > 0 ) {
                        _popup_buffer.append( "\n" );
                    }
                    final Sequence seq = node.getNodeData().getSequence();
                    if ( seq.getAccession() != null ) {
                        _popup_buffer.append( "[" );
                        if ( !ForesterUtil.isEmpty( seq.getAccession().getSource() ) ) {
                            _popup_buffer.append( seq.getAccession().getSource() );
                            _popup_buffer.append( ":" );
                        }
                        _popup_buffer.append( seq.getAccession().getValue() );
                        _popup_buffer.append( "]" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " [" );
                        }
                        else {
                            _popup_buffer.append( "[" );
                        }
                        _popup_buffer.append( seq.getSymbol() );
                        _popup_buffer.append( "]" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " [" );
                        }
                        else {
                            _popup_buffer.append( "[" );
                        }
                        _popup_buffer.append( seq.getGeneName() );
                        _popup_buffer.append( "]" );
                        enc_data = true;
                    }
                    if ( !ForesterUtil.isEmpty( seq.getName() ) ) {
                        if ( enc_data ) {
                            _popup_buffer.append( " " );
                        }
                        _popup_buffer.append( seq.getName() );
                    }
                }
                if ( node.getNodeData().isHasDate() ) {
                    lines++;
                    if ( _popup_buffer.length() > 0 ) {
                        _popup_buffer.append( "\n" );
                    }
                    _popup_buffer.append( node.getNodeData().getDate().asSimpleText() );
                }
                if ( node.getNodeData().isHasDistribution() ) {
                    lines++;
                    if ( _popup_buffer.length() > 0 ) {
                        _popup_buffer.append( "\n" );
                    }
                    _popup_buffer.append( node.getNodeData().getDistribution().asSimpleText() );
                }
                if ( node.getBranchData().isHasConfidences() ) {
                    final List<Confidence> confs = node.getBranchData().getConfidences();
                    for( final Confidence confidence : confs ) {
                        lines++;
                        if ( _popup_buffer.length() > 0 ) {
                            _popup_buffer.append( "\n" );
                        }
                        if ( !ForesterUtil.isEmpty( confidence.getType() ) ) {
                            _popup_buffer.append( "[" );
                            _popup_buffer.append( confidence.getType() );
                            _popup_buffer.append( "] " );
                        }
                        _popup_buffer
                                .append( FORMATTER_CONFIDENCE.format( ForesterUtil.round( confidence.getValue(),
                                                                                          getOptions()
                                                                                                  .getNumberOfDigitsAfterCommaForConfidenceValues() ) ) );
                        if ( confidence.getStandardDeviation() != Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                            _popup_buffer.append( " (sd=" );
                            _popup_buffer.append( FORMATTER_CONFIDENCE.format( ForesterUtil.round( confidence
                                    .getStandardDeviation(), getOptions()
                                    .getNumberOfDigitsAfterCommaForConfidenceValues() ) ) );
                            _popup_buffer.append( ")" );
                        }
                    }
                }
                if ( node.getNodeData().isHasProperties() ) {
                    final PropertiesMap properties = node.getNodeData().getProperties();
                    for( final String ref : properties.getPropertyRefs() ) {
                        _popup_buffer.append( "\n" );
                        final Property p = properties.getProperty( ref );
                        _popup_buffer.append( TreePanelUtil.getPartAfterColon( p.getRef() ) );
                        _popup_buffer.append( "=" );
                        _popup_buffer.append( p.getValue() );
                        if ( !ForesterUtil.isEmpty( p.getUnit() ) ) {
                            _popup_buffer.append( TreePanelUtil.getPartAfterColon( p.getUnit() ) );
                        }
                    }
                }
                if ( _popup_buffer.length() > 0 ) {
                    if ( !getConfiguration().isUseNativeUI() ) {
                        _rollover_popup
                                .setBorder( BorderFactory.createLineBorder( getTreeColorSet().getBranchColor() ) );
                        _rollover_popup.setBackground( getTreeColorSet().getBackgroundColor() );
                        if ( isInFoundNodes0( node ) && !isInFoundNodes1( node ) ) {
                            _rollover_popup.setForeground( getTreeColorSet().getFoundColor0() );
                        }
                        else if ( !isInFoundNodes0( node ) && isInFoundNodes1( node ) ) {
                            _rollover_popup.setForeground( getTreeColorSet().getFoundColor1() );
                        }
                        else if ( isInFoundNodes0( node ) && isInFoundNodes1( node ) ) {
                            _rollover_popup.setForeground( getTreeColorSet().getFoundColor0and1() );
                        }
                        else {
                            _rollover_popup.setForeground( getTreeColorSet().getSequenceColor() );
                        }
                    }
                    else {
                        _rollover_popup.setBorder( BorderFactory.createLineBorder( Color.BLACK ) );
                    }
                    _rollover_popup.setText( _popup_buffer.toString() );
                    _node_desc_popup = PopupFactory.getSharedInstance().getPopup( null,
                                                                                  _rollover_popup,
                                                                                  e.getLocationOnScreen().x + 10,
                                                                                  e.getLocationOnScreen().y
                                                                                          - ( lines * 20 ) );
                    _node_desc_popup.show();
                }
            }
        }
        catch ( final Exception ex ) {
            // Do nothing.
        }
    }

    final private void showNodeEditFrame( final PhylogenyNode n ) {
        if ( _node_frame_index < TreePanel.MAX_NODE_FRAMES ) {
            // pop up edit box for single node
            _node_frames[ _node_frame_index ] = new NodeFrame( n, _phylogeny, this, _node_frame_index, "" );
            _node_frame_index++;
        }
        else {
            JOptionPane.showMessageDialog( this, "too many node windows are open" );
        }
    }

    final private void showNodeFrame( final PhylogenyNode n ) {
        if ( _node_frame_index < TreePanel.MAX_NODE_FRAMES ) {
            // pop up edit box for single node
            _node_frames[ _node_frame_index ] = new NodeFrame( n, _phylogeny, this, _node_frame_index );
            _node_frame_index++;
        }
        else {
            JOptionPane.showMessageDialog( this, "too many node windows are open" );
        }
    }

    final private void switchDisplaygetPhylogenyGraphicsType() {
        switch ( getPhylogenyGraphicsType() ) {
            case RECTANGULAR:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
                break;
            case EURO_STYLE:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
                break;
            case ROUNDED:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CURVED );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CURVED );
                break;
            case CURVED:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR );
                break;
            case TRIANGULAR:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CONVEX );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CONVEX );
                break;
            case CONVEX:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                break;
            case UNROOTED:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                break;
            case CIRCULAR:
                setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                getOptions().setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                break;
            default:
                throw new RuntimeException( "unkwnown display type: " + getPhylogenyGraphicsType() );
        }
        if ( getControlPanel().getDynamicallyHideData() != null ) {
            if ( getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
                getControlPanel().getDynamicallyHideData().setEnabled( false );
            }
            else {
                getControlPanel().getDynamicallyHideData().setEnabled( true );
            }
        }
        if ( isPhyHasBranchLengths() && ( getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) ) {
            getControlPanel().setDrawPhylogramEnabled( true );
        }
        else {
            getControlPanel().setDrawPhylogramEnabled( false );
        }
        if ( getMainPanel().getMainFrame() == null ) {
            // Must be "E" applet version.
            ( ( ArchaeopteryxE ) ( ( MainPanelApplets ) getMainPanel() ).getApplet() )
                    .setSelectedTypeInTypeMenu( getPhylogenyGraphicsType() );
        }
        else {
            getMainPanel().getMainFrame().setSelectedTypeInTypeMenu( getPhylogenyGraphicsType() );
        }
    }

    private final static void colorizeNodesHelper( final Color c, final PhylogenyNode node ) {
        if ( node.getNodeData().getNodeVisualData() == null ) {
            node.getNodeData().setNodeVisualData( new NodeVisualData() );
        }
        node.getNodeData().getNodeVisualData().setFontColor( new Color( c.getRed(), c.getGreen(), c.getBlue() ) );
    }

    final private static void drawString( final String str, final float x, final float y, final Graphics2D g ) {
        g.drawString( str, x, y );
    }

    final private static boolean plusPressed( final int key_code ) {
        return ( ( key_code == KeyEvent.VK_ADD ) || ( key_code == KeyEvent.VK_PLUS )
                || ( key_code == KeyEvent.VK_EQUALS ) || ( key_code == KeyEvent.VK_SEMICOLON ) || ( key_code == KeyEvent.VK_1 ) );
    }

    final private class NodeColorizationActionListener implements ActionListener {

        List<PhylogenyNode> _additional_nodes = null;
        JColorChooser       _chooser          = null;
        PhylogenyNode       _node             = null;

        NodeColorizationActionListener( final JColorChooser chooser, final PhylogenyNode node ) {
            _chooser = chooser;
            _node = node;
        }

        NodeColorizationActionListener( final JColorChooser chooser,
                                        final PhylogenyNode node,
                                        final List<PhylogenyNode> additional_nodes ) {
            _chooser = chooser;
            _node = node;
            _additional_nodes = additional_nodes;
        }

        @Override
        public void actionPerformed( final ActionEvent e ) {
            final Color c = _chooser.getColor();
            if ( c != null ) {
                colorizeNodes( c, _node, _additional_nodes );
            }
        }
    }

    final private class SubtreeColorizationActionListener implements ActionListener {

        List<PhylogenyNode> _additional_nodes = null;
        JColorChooser       _chooser          = null;
        PhylogenyNode       _node             = null;

        SubtreeColorizationActionListener( final JColorChooser chooser, final PhylogenyNode node ) {
            _chooser = chooser;
            _node = node;
        }

        SubtreeColorizationActionListener( final JColorChooser chooser,
                                           final PhylogenyNode node,
                                           final List<PhylogenyNode> additional_nodes ) {
            _chooser = chooser;
            _node = node;
            _additional_nodes = additional_nodes;
        }

        @Override
        public void actionPerformed( final ActionEvent e ) {
            final Color c = _chooser.getColor();
            if ( c != null ) {
                colorizeSubtree( c, _node, _additional_nodes );
            }
        }
    }
}
