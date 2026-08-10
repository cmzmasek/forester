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
import java.awt.Toolkit;
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
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
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
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;

import javax.swing.BorderFactory;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.Popup;
import javax.swing.PopupFactory;
import javax.swing.Timer;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.ControlPanel.NodeClickAction;
import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.archaeopteryx.phylogeny.data.RenderableMsaSequence;
import org.forester.archaeopteryx.phylogeny.data.RenderableVector;
import org.forester.archaeopteryx.tools.Blast;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.SequenceRelation;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;
import org.forester.util.TaxonomyUtil;

public final class TreePanel extends JPanel implements ActionListener, MouseWheelListener, Printable {

    final private class NodeColorizationActionListener implements ActionListener {

        List<PhylogenyNode> _additional_nodes = null;
        JColorChooser _chooser = null;
        PhylogenyNode _node = null;

        NodeColorizationActionListener(final JColorChooser chooser, final PhylogenyNode node) {
            _chooser = chooser;
            _node = node;
        }

        NodeColorizationActionListener(final JColorChooser chooser,
                                       final PhylogenyNode node,
                                       final List<PhylogenyNode> additional_nodes) {
            _chooser = chooser;
            _node = node;
            _additional_nodes = additional_nodes;
        }

        @Override
        public void actionPerformed(final ActionEvent e) {
            final Color c = _chooser.getColor();
            if (c != null) {
                colorizeNodes(c, _node, _additional_nodes);
            }
        }
    }

    final private class SubtreeColorizationActionListener implements ActionListener {

        List<PhylogenyNode> _additional_nodes = null;
        JColorChooser _chooser = null;
        PhylogenyNode _node = null;

        SubtreeColorizationActionListener(final JColorChooser chooser, final PhylogenyNode node) {
            _chooser = chooser;
            _node = node;
        }

        SubtreeColorizationActionListener(final JColorChooser chooser,
                                          final PhylogenyNode node,
                                          final List<PhylogenyNode> additional_nodes) {
            _chooser = chooser;
            _node = node;
            _additional_nodes = additional_nodes;
        }

        @Override
        public void actionPerformed(final ActionEvent e) {
            final Color c = _chooser.getColor();
            if (c != null) {
                colorizeSubtree(c, _node, _additional_nodes);
            }
        }
    }

    public final static boolean SPECIAL_DOMAIN_COLORING = true;
    final static Cursor ARROW_CURSOR = Cursor
            .getPredefinedCursor(Cursor.DEFAULT_CURSOR);
    final static Cursor CUT_CURSOR = Cursor
            .getPredefinedCursor(Cursor.CROSSHAIR_CURSOR);
    final static Cursor HAND_CURSOR = Cursor
            .getPredefinedCursor(Cursor.HAND_CURSOR);
    final static Cursor MOVE_CURSOR = Cursor
            .getPredefinedCursor(Cursor.MOVE_CURSOR);
    final static Cursor WAIT_CURSOR = Cursor
            .getPredefinedCursor(Cursor.WAIT_CURSOR);
    final private static double _180_OVER_PI = 180.0 / Math.PI;
    private static final float ANGLE_ROTATION_UNIT = (float) (Math.PI
            / 32);
    private final static int CONFIDENCE_LEFT_MARGIN = 4;
    // Gap (px) between the node shape and the first label glyph. The taxonomy label and the node-data label that
    // follows it MUST use the same value (via labelSegmentStartX) so the node data begins exactly where the
    // taxonomy label ends. When the two differed (taxonomy 3, node data 2), the node data started a pixel inside
    // the taxonomy box and an italic scientific name's right overhang overlapped the following node name.
    private final static int LABEL_GAP_AFTER_NODE_SHAPE = 2;
    // Extra horizontal gap between the taxonomy and node-data segments of an above-the-branch internal label, ON
    // TOP of the trailing space the taxonomy segment already carries (taxonomyLabel emits a trailing " "). Kept at
    // 0 so an internal label's taxonomy->name spacing matches the external path (which separates the two with only
    // that trailing space); a non-zero value here double-counts the gap and reads as too wide.
    private final static int INTERNAL_LABEL_SEGMENT_GAP = 0;
    // Smallest x an above-the-branch internal label may start at: a long label on an internal node near
    // the root would otherwise grow leftward off-canvas and lose its leading characters, so it is shifted
    // right to keep its leftmost glyph visible at this margin instead of being clipped.
    private final static float INTERNAL_LABEL_MIN_LEFT_MARGIN = 2;
    // Fractional-metrics context for measuring label advance on the vector-export path; see
    // fractionalAdvanceWidth().
    private final static FontRenderContext FRC_FRACTIONAL = new FontRenderContext(null, true, true);
    private final static int EURO_D = 10;
    private final static NumberFormat FORMATTER_BRANCH_LENGTH;
    private final static NumberFormat FORMATTER_CONFIDENCE;
    // MAD support values lie in [0,1], so they get their own fixed (up to 2-decimal) format, kept
    // independent of the "digits after comma" setting that the regular confidences use.
    private final static NumberFormat FORMATTER_MAD;
    private static final float HALF_PI = (float) (Math.PI
            / 2.0);
    private final static int LIMIT_FOR_HQ_RENDERING = 2000;
    private final static int MAX_NODE_FRAMES = 10;
    private final static int MAX_SUBTREES = 100;
    private final static int MIN_ROOT_LENGTH = 3;
    private final static int MOVE = 20;
    private final static String NODE_POPMENU_NODE_CLIENT_PROPERTY = "node";
    private static final float ONEHALF_PI = (float) (1.5
            * Math.PI);
    private static final short OV_BORDER = 10;
    private final static double OVERVIEW_FOUND_NODE_BOX_SIZE = 2;
    private final static double OVERVIEW_FOUND_NODE_BOX_SIZE_HALF = 1;
    private static final float PI = (float) (Math.PI);
    final private static Font POPUP_FONT = new Font(Configuration
            .getDefaultFontFamilyName(), Font.PLAIN, 12);
    private static final float ROUNDED_D = 8;
    private final static long serialVersionUID = -978349745916505029L;
    private static final BasicStroke STROKE_0025 = new BasicStroke(0.025f);
    private static final BasicStroke STROKE_005 = new BasicStroke(0.05f);
    private static final BasicStroke STROKE_01 = new BasicStroke(0.1f);
    private static final BasicStroke STROKE_025 = new BasicStroke(0.25f);
    private static final BasicStroke STROKE_05 = new BasicStroke(0.5f);
    private static final BasicStroke STROKE_075 = new BasicStroke(0.75f);
    private static final BasicStroke STROKE_1 = new BasicStroke(1f);
    private static final BasicStroke STROKE_2 = new BasicStroke(2f);


    // Fine dotted "leader" that links each tip to its lined-up label/data in aligned-phylogram mode (the
    // iTOL-style guide line). Round-capped dots at a visible width -- the old per-branch strokes were
    // sub-pixel (0.01-0.1f) and rendered effectively invisible.
    private static final BasicStroke LEADER_STROKE = new BasicStroke(0.7f,
            BasicStroke.CAP_ROUND,
            BasicStroke.JOIN_ROUND,
            1.0f,
            new float[]{
                    1.0f, 2.0f},
            0f);
    // Neutral guide color for the lined-up-data connector (see connectorColor()/drawConnection).
    private static final Color CONNECTOR_GUIDE_COLOR = new Color(200, 200, 200);
    // How far the (faint) scale-grid color is blended from the background toward the branch-length color.
    private static final double SCALE_GRID_BLEND = 0.18;
    private static final int    SCALE_AXIS_TICK_LEN = 4;  // length (px) of the labeled scale-axis tick marks
    private static final int    SCALE_AXIS_LABEL_GAP = 4; // min px between adjacent tick labels (else the label is decimated)
    private static final int    SCALE_AXIS_UNIT_GAP = 5;  // gap before the trailing [unit] label
    private static final int    VERTICAL_DOMAIN_GAP = 8;  // gap (px) between a tilted tip label and its domain track (vertical)
    private static final int    VERTICAL_BREADTH_PAD = 12; // extra breadth margin (px) per side after fitting a vertical tree
    // Node age (HPD) bars: the bar thickness and its translucent fill (blue on screen/color export; a neutral gray in
    // a black-and-white export so it is not the only colored element).
    private static final int    HPD_BAR_HEIGHT = 7;
    private static final Color  HPD_BAR_COLOR = new Color(70, 130, 220, 90);  // translucent blue, FigTree-like
    private static final Color  HPD_BAR_COLOR_BW = new Color(90, 90, 90, 70); // translucent gray for B&W export
    // Zebra row stripes: a faint translucent band, darkening a light background / lightening a dark one.
    private static final Color  ZEBRA_STRIPE_ON_LIGHT = new Color(0, 0, 0, 16);
    private static final Color  ZEBRA_STRIPE_ON_DARK = new Color(255, 255, 255, 20);
    // "Dim Non-Matches": how far a non-hit label is blended toward the background while a search is active.
    private static final double DIM_NON_MATCH_FRACTION = 0.72;
    // "Pulse Found Nodes": a translucent found-color halo disc around each hit; on screen its radius+alpha breathe
    // over PULSE_PERIOD_MS, in exports it renders once at peak (static glow). Repainted at ~PULSE_INTERVAL_MS.
    private static final int    PULSE_INTERVAL_MS = 55;
    private static final int    PULSE_PERIOD_MS = 1300;
    private static final long   PULSE_PERIOD_NS = PULSE_PERIOD_MS * 1_000_000L; // monotonic clock base for the phase
    private static final float  HALO_BASE_RADIUS = 4f;   // radius at the trough of the pulse
    private static final float  HALO_AMP_RADIUS = 6f;    // added at the peak (so radius breathes 4 -> 10 px)
    private static final int    HALO_MIN_ALPHA = 18;
    private static final int    HALO_MAX_ALPHA = 70;
    private static final double TWO_PI = 2 * Math.PI;
    private final static int WIGGLE = 2;
    HashMap<Long, Short> _nodeid_dist_to_leaf = new HashMap<>();
    final private Arc2D _arc = new Arc2D.Double();
    private AffineTransform _at;
    // Root-on-top / root-on-bottom orientation: the tree is laid out LOGICALLY (depth->Xcoord, breadth->Ycoord)
    // exactly as always, then the whole rectangular canvas is rotated by _orientation_R (rebuilt each paint from
    // the logical extents). Geometry rides R for free; text is re-anchored upright/45deg by withNodeTextFrame, and
    // the viewport-fixed chrome is drawn after restoring _orientation_base_transform (the pre-rotation CTM).
    // _orientation_R_inverse maps device points back to logical space for hit-testing. All null in ROOT_LEFT.
    private AffineTransform _orientation_R;
    private AffineTransform _orientation_R_inverse;
    private AffineTransform _orientation_base_transform;
    // R depends only on the orientation + the logical extents (a function of the layout params + tree structure); those
    // change only via calcParametersForPainting / resetPreferredSize (never during a plain repaint), so R is cached and
    // rebuilt lazily -- a hover/pulse/scroll repaint reuses it instead of re-walking the tree (calculateHeight) for R.
    private boolean _orientation_transform_dirty = true;
    private Options.TREE_ORIENTATION _orientation_R_built_for;
    private int _circ_max_depth;
    final private Set<Long> _collapsed_external_nodeid_set = new HashSet<>();
    private JColorChooser _color_chooser = null;
    private Configuration _configuration = null;
    private ControlPanel _control_panel = null;
    private int _domain_structure_e_value_thr_exp = AptxConstants.DOMAIN_STRUCTURE_E_VALUE_THR_DEFAULT_EXP;
    private double _domain_structure_width = AptxConstants.DOMAIN_STRUCTURE_DEFAULT_WIDTH;
    private int _dynamic_hiding_factor = 0;
    private boolean _edited = false;
    private final Ellipse2D _ellipse = new Ellipse2D.Float();
    private int _external_node_index = 0;
    private Set<Long> _found_nodes_0 = null;
    private Set<Long> _found_nodes_1 = null;
    // The target currently previewed on hover in Select-Node(s) mode (or null): a single node (_hover_subtree
    // false -> one marker on the node) or a branch (_hover_subtree true -> markers on the subtree's tips), drawn
    // translucent as "will be added / will be removed" so a click's effect is visible before committing.
    private PhylogenyNode _hover_node = null;
    private boolean       _hover_subtree = false;
    // The target just clicked: its hover preview stays suppressed until the pointer moves off it, so it does not
    // instantly flip to the "will be removed" grey over what the click just selected.
    private PhylogenyNode _click_suppressed = null;
    private final FontRenderContext _frc = new FontRenderContext(null,
            false,
            false);
    private PHYLOGENY_GRAPHICS_TYPE _graphics_type = PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR;
    private boolean _in_ov = false;
    private boolean _in_ov_rect = false;
    private float _last_drag_point_x = 0;
    private float _last_drag_point_y = 0;
    private final Line2D _line = new Line2D.Float();
    private int _longest_ext_node_info = 0;
    private PhylogenyNode _ext_node_with_longest_txt_info = null;
    private MainPanel _main_panel = null;
    private double _max_distance_to_root = -1;
    private Popup _node_desc_popup;
    private int _node_frame_index = 0;
    private final NodeFrame[] _node_frames = new NodeFrame[TreePanel.MAX_NODE_FRAMES];
    private JPopupMenu _node_popup_menu = null;
    private JMenuItem _node_popup_menu_items[] = null;
    private PhylogenyNode[] _nodes_in_preorder = null;
    private Options _options = null;
    private float _ov_max_height = 0;
    private float _ov_max_width = 0;
    private boolean _ov_on = false;
    private final Rectangle2D _ov_rectangle = new Rectangle2D.Float();
    private final Rectangle _ov_virtual_rectangle = new Rectangle();
    private float _ov_x_correction_factor = 0.0f;
    private float _ov_x_distance = 0;
    private int _ov_x_position = 0;
    private float _ov_y_distance = 0;
    private int _ov_y_position = 0;
    private int _ov_y_start = 0;
    private final boolean _phy_has_branch_lengths;
    private Phylogeny _phylogeny = null;
    private final Path2D.Float _polygon = new Path2D.Float();
    private final StringBuffer _popup_buffer = new StringBuffer();
    private Sequence _query_sequence = null;
    private final Rectangle2D _rectangle = new Rectangle2D.Float();
    private final RenderingHints _rendering_hints = new RenderingHints(RenderingHints.KEY_RENDERING,
            RenderingHints.VALUE_RENDER_DEFAULT);
    private JTextArea _rollover_popup;
    private PhylogenyNode _root;
    private final StringBuilder _sb = new StringBuilder();
    // Set transiently by AptxUtil around a PNG export so paintPhylogeny leaves the background unfilled
    // (transparent). Off for screen and every other export format.
    private boolean _export_transparent_background = false;
    // Cache for the italic-derived scientific-name font: deriveFont() allocates, and taxonomyLabel runs
    // per node per repaint, so re-derive only when the base font actually changes (see italicOf()).
    private Font _italic_base_font;
    private Font _italic_derived_font;
    // Same cache for the bold-derived found-label font (see boldOf() / setFont). Also composes with italicOf (a
    // found scientific name derives bold-italic); the italic cache's single slot still re-derives at each
    // found<->non-found boundary, but that is bounded per-boundary, not per-label, and bold is off by default.
    private Font _bold_base_font;
    private Font _bold_derived_font;
    // "Dim Non-Matches": computed once per paint -- true only when the dim option is on AND at least one found
    // node is actually DRAWN (not all hidden under a collapse / absent from the displayed subtree), so the tree
    // never washes out with nothing emphasised. Read per-label in setColor.
    private boolean _has_visible_found_node = false;
    // "Pulse Found Nodes": the EDT timer that drives the on-screen animation, the drawn halo bounds it repaints
    // (only those small regions, not the whole canvas), and whether the last screen paint drew any halo.
    private Timer   _pulse_timer;
    private final List<Rectangle> _found_halo_bounds = new ArrayList<>();
    private boolean _has_visible_found_halo = false;
    // Next/previous "step through search hits": the current position in the ordered hit list (-1 = not positioned)
    // and the last node centered on (a collapsed hit's drawn triangle, else the hit itself) -- for tests.
    private int          _search_hit_index = -1;
    private PhylogenyNode _last_step_target;
    private double _scale_distance = 0.0;
    private String _scale_label = null;
    private DescriptiveStatistics _statistics_for_vector_data;
    private final Phylogeny[] _sub_phylogenies = new Phylogeny[TreePanel.MAX_SUBTREES];
    private final PhylogenyNode[] _sub_phylogenies_temp_roots = new PhylogenyNode[TreePanel.MAX_SUBTREES];
    private int _subtree_index = 0;
    private File _treefile = null;
    private float _urt_factor = 1;
    private float _urt_factor_ov = 1;
    final private HashMap<Long, Double> _urt_nodeid_angle_map = new HashMap<>();
    final private HashMap<Long, Integer> _urt_nodeid_index_map = new HashMap<>();
    private double _urt_starting_angle = (float) (Math.PI
            / 2);
    private float _x_correction_factor = 0.0f;
    private float _x_distance = 0.0f;
    private float _y_distance = 0.0f;
    // support-scale ceiling (1 or 100) for node-symbol support visualization; recomputed per paint
    private double _confidence_scale_max = 1.0;
    private int _length_of_longest_text;
    private int _length_of_longest_text_only; // longest tip TEXT label (name+taxonomy) only, excluding domain/vector tracks
    // AUTO tip-label direction resolved ONCE per calcParametersForPainting pass (null = option is not AUTO). Caching it
    // keeps the breadth RESERVES and the PAINT in agreement within a pass -- both read this instead of re-resolving AUTO
    // against a y-distance that the reserves see stale (pre-setYdistance) but the paint sees fresh, which clipped the
    // outermost upright labels when the fresh spacing flipped the resolved angle. See effectiveTipLabelDirection().
    private Options.TIP_LABEL_DIRECTION _resolved_auto_tip_dir = null;
    private int _longest_domain;
    private double _longest_rendered_domain; // widest domain track in px (rendered width), for the vertical depth reserve

    static {
        final DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator('.');
        FORMATTER_CONFIDENCE = new DecimalFormat("#.###", dfs);
        FORMATTER_BRANCH_LENGTH = new DecimalFormat("#.###", dfs);
        FORMATTER_MAD = new DecimalFormat("0.##", dfs); // e.g. 0 -> "0", 0.05 -> "0.05", 0.123 -> "0.12"
    }

    TreePanel(final Phylogeny t, final Configuration configuration, final MainPanel tjp) {
        requestFocusInWindow();
        addKeyListener(new KeyAdapter() {

            @Override
            public void keyPressed(final KeyEvent key_event) {
                keyPressedCalls(key_event);
                requestFocusInWindow();
            }
        });
        addFocusListener(new FocusAdapter() {

            @Override
            public void focusGained(final FocusEvent e) {
                requestFocusInWindow();
            }
        });
        if ((t == null) || t.isEmpty()) {
            throw new IllegalArgumentException("attempt to draw phylogeny which is null or empty");
        }
        _graphics_type = tjp.getOptions().getPhylogenyGraphicsType();
        // seed the "Color by" palette from the shared default so a new tab inherits the last-chosen palette
        if (!ForesterUtil.isEmpty(tjp.getOptions().getColorPaletteName())) {
            _color_palette_name = tjp.getOptions().getColorPaletteName();
        }
        _main_panel = tjp;
        _configuration = configuration;
        _phylogeny = t;
        _phy_has_branch_lengths = AptxUtil.isHasAtLeastOneBranchLengthLargerThanZero(_phylogeny);
        init();

        _phylogeny.recalculateNumberOfExternalDescendants(true);

        setBackground(getTreeColorSet().getBackgroundColor());
        final MouseListener mouse_listener = new MouseListener(this);
        addMouseListener(mouse_listener);
        addMouseMotionListener(mouse_listener);
        addMouseWheelListener(this);
        calculateScaleDistance();
        FORMATTER_CONFIDENCE.setMaximumFractionDigits(configuration.getNumberOfDigitsAfterCommaForConfidenceValues());
        FORMATTER_BRANCH_LENGTH
                .setMaximumFractionDigits(configuration.getNumberOfDigitsAfterCommaForBranchLengthValues());
    }

    @Override
    final public void actionPerformed(final ActionEvent e) {
        boolean done = false;
        final JMenuItem node_popup_menu_item = (JMenuItem) e.getSource();
        for (int index = 0; (index < _node_popup_menu_items.length) && !done; index++) {
            // NOTE: index corresponds to the indices of click-to options
            // in the control panel.
            if (node_popup_menu_item == _node_popup_menu_items[index]) {
                // Set this as the new default click-to action
                _main_panel.getControlPanel().setClickToAction(index);
                final PhylogenyNode node = (PhylogenyNode) _node_popup_menu
                        .getClientProperty(NODE_POPMENU_NODE_CLIENT_PROPERTY);
                handleClickToAction(_control_panel.getActionWhenNodeClicked(), node);
                done = true;
            }
        }
        repaint();
        requestFocusInWindow();
    }

    private static BasicStroke makeStroke(final float width) {
        return new BasicStroke(width);
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

    public final TreeColorSet getTreeColorSet() {
        return getMainPanel().getTreeColorSet();
    }

    @Override
    final public void mouseWheelMoved(final MouseWheelEvent e) {
        final int notches = e.getWheelRotation();
        if (inOvVirtualRectangle(e)) {
            if (!isInOvRect()) {
                setInOvRect(true);
                repaint();
            }
        } else {
            if (isInOvRect()) {
                setInOvRect(false);
                repaint();
            }
        }
        if (e.isControlDown() && e.isShiftDown()) {
            if (notches < 0) {
                getTreeFontSet().increaseUserFontSize();
            } else {
                getTreeFontSet().decreaseUserFontSize();
            }
            getControlPanel().displayedPhylogenyMightHaveChanged(true);
            resetPreferredSize();
            updateOvSizes();
            repaint();
        } else if (e.isShiftDown() && e.isAltDown()) {
            // horizontal-on-screen zoom (flips to the tip-spread in a vertical orientation)
            if (notches < 0) {
                for (int i = 0; i < (-notches); ++i) {
                    getControlPanel().zoomInScreenX(AptxConstants.WHEEL_ZOOM_IN_FACTOR,
                            AptxConstants.WHEEL_ZOOM_IN_FACTOR);
                    getControlPanel().displayedPhylogenyMightHaveChanged(false);
                }
            } else {
                for (int i = 0; i < notches; ++i) {
                    getControlPanel().zoomOutScreenX(AptxConstants.WHEEL_ZOOM_OUT_FACTOR,
                            AptxConstants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR);
                    getControlPanel().displayedPhylogenyMightHaveChanged(false);
                }
            }
        } else if (e.isShiftDown()) {
            if ((getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                    || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
                if (notches < 0) {
                    for (int i = 0; i < (-notches); ++i) {
                        setStartingAngle((getStartingAngle() % TWO_PI) + ANGLE_ROTATION_UNIT);
                        getControlPanel().displayedPhylogenyMightHaveChanged(false);
                    }
                } else {
                    for (int i = 0; i < notches; ++i) {
                        setStartingAngle((getStartingAngle() % TWO_PI) - ANGLE_ROTATION_UNIT);
                        if (getStartingAngle() < 0) {
                            setStartingAngle(TWO_PI + getStartingAngle());
                        }
                        getControlPanel().displayedPhylogenyMightHaveChanged(false);
                    }
                }
            } else {
                // vertical-on-screen zoom (flips to the depth axis in a vertical orientation)
                if (notches < 0) {
                    for (int i = 0; i < (-notches); ++i) {
                        getControlPanel().zoomInScreenY(AptxConstants.WHEEL_ZOOM_IN_FACTOR,
                                AptxConstants.WHEEL_ZOOM_IN_X_CORRECTION_FACTOR);
                        getControlPanel().displayedPhylogenyMightHaveChanged(false);
                    }
                } else {
                    for (int i = 0; i < notches; ++i) {
                        getControlPanel().zoomOutScreenY(AptxConstants.WHEEL_ZOOM_OUT_FACTOR,
                                AptxConstants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR);
                        getControlPanel().displayedPhylogenyMightHaveChanged(false);
                    }
                }
            }
        } else {
            if (notches < 0) {
                for (int i = 0; i < (-notches); ++i) {
                    getControlPanel().zoomInX(AptxConstants.WHEEL_ZOOM_IN_FACTOR,
                            AptxConstants.WHEEL_ZOOM_IN_X_CORRECTION_FACTOR);
                    getControlPanel().zoomInY(AptxConstants.WHEEL_ZOOM_IN_FACTOR);
                    getControlPanel().displayedPhylogenyMightHaveChanged(false);
                }
            } else {
                for (int i = 0; i < notches; ++i) {
                    getControlPanel().zoomOutY(AptxConstants.WHEEL_ZOOM_OUT_FACTOR);
                    getControlPanel().zoomOutX(AptxConstants.WHEEL_ZOOM_OUT_FACTOR,
                            AptxConstants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR);
                    getControlPanel().displayedPhylogenyMightHaveChanged(false);
                }
            }
        }
        requestFocus();
        requestFocusInWindow();
        requestFocus();
    }

    @Override
    final public void paintComponent(final Graphics g) {
        final Graphics2D g2d = (Graphics2D) g;
        g2d.setRenderingHints(_rendering_hints);
        paintPhylogeny(g2d, false, false, 0, 0, 0, 0);
    }

    @Override
    final public int print(final Graphics g, final PageFormat page_format, final int page_index)
            throws PrinterException {
        if (page_index > 0) {
            return (NO_SUCH_PAGE);
        } else {
            final Graphics2D g2d = (Graphics2D) g;
            g2d.translate(page_format.getImageableX(), page_format.getImageableY());
            // Turn off double buffering !?
            paintPhylogeny(g2d, true, false, 0, 0, 0, 0);
            // Turn double buffering back on !?
            return (PAGE_EXISTS);
        }
    }

    public final void setEdited(final boolean edited) {
        _edited = edited;
        // Undo/redo safety net: any edit made outside of an undo/redo restore invalidates the redo history
        // (it no longer describes a reachable future). Checkpointed ops already clear redo via checkpoint();
        // this also covers mutations that were NOT checkpointed, so a later Redo can never install a tree
        // unrelated to the current one. (The undo stack is left intact -- worst case an undo reverts one such
        // uncheckpointed edit too, which is recoverable via redo, rather than corrupting the tree.)
        if (edited && !_restoring_snapshot && (_history != null)) {
            _history.clearRedo();
            notifyEditMenu();
        }
    }


    /**
     * Set a phylogeny tree.
     *
     * @param t an instance of a Phylogeny
     */
    public final void setTree(final Phylogeny t) {
        setNodeInPreorderToNull(); // also clears any stale branch-hover preview
        _phylogeny = t;
    }

    public final void setWaitCursor() {
        setCursor(WAIT_CURSOR);
        repaint();
    }

    @Override
    public void update(final Graphics g) {
        paint(g);
    }

    final private void addEmptyNode(final PhylogenyNode node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        final String label = createASimpleTextRepresentationOfANode(node);
        String msg = "";
        if (ForesterUtil.isEmpty(label)) {
            msg = "How to add the new, empty node?";
        } else {
            msg = "How to add the new, empty node to node" + label + "?";
        }
        final Object[] options = {"As sibling", "As descendant", "Cancel"};
        final int r = JOptionPane.showOptionDialog(this,
                msg,
                "Addition of Empty New Node",
                JOptionPane.CLOSED_OPTION,
                JOptionPane.QUESTION_MESSAGE,
                null,
                options,
                options[2]);
        boolean add_as_sibling = true;
        if (r == 1) {
            add_as_sibling = false;
        } else if (r != 0) {
            return;
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot(new PhylogenyNode());
        phy.setRooted(true);
        if (add_as_sibling) {
            if (node.isRoot()) {
                JOptionPane.showMessageDialog(this,
                        "Cannot add sibling to root",
                        "Attempt to add sibling to root",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }
            pushUndoCheckpoint("Add Node");
            phy.addAsSibling(node);
        } else {
            pushUndoCheckpoint("Add Node");
            phy.addAsChild(node);
        }
        setNodeInPreorderToNull();
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants(true);
        resetNodeIdToDistToLeafMap();
        setEdited(true);
        repaint();
    }

    final private void assignGraphicsForBranchWithColorForParentBranch(final PhylogenyNode node,
                                                                       final boolean is_vertical,
                                                                       final Graphics g,
                                                                       final boolean to_pdf,
                                                                       final boolean to_graphics_file) {
        final NodeClickAction action = _control_panel.getActionWhenNodeClicked();
        if ((to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite()) {
            g.setColor(Color.BLACK);
        } else if (((action == NodeClickAction.COPY_SUBTREE) || (action == NodeClickAction.CUT_SUBTREE)
                || (action == NodeClickAction.DELETE_NODE_OR_SUBTREE) || (action == NodeClickAction.PASTE_SUBTREE)
                || (action == NodeClickAction.ADD_NEW_NODE)) && (getCutOrCopiedTree() != null)
                && (getCopiedAndPastedNodes() != null) && !to_pdf && !to_graphics_file
                && getCopiedAndPastedNodes().contains(node.getId())) {
            g.setColor(getTreeColorSet().getFoundColor0());
        } else if ((getOptions().getSupportVisualization() == Options.SUPPORT_VISUALIZATION.COLOR_BRANCHES)
                && !node.isExternal() && node.getBranchData().isHasConfidences()
                && (PhylogenyMethods.getConfidenceValue(node) >= 0.0)) {
            g.setColor(supportBranchColor(node, to_pdf));
        } else if (getControlPanel().isUseVisualStyles() && (PhylogenyMethods.getBranchColorValue(node) != null)) {
            g.setColor(PhylogenyMethods.getBranchColorValue(node));
        } else if (to_pdf) {
            g.setColor(getTreeColorSet().getBranchColorForPdf());
        } else {
            g.setColor(getTreeColorSet().getBranchColor());
        }
    }

    /**
     * The branch (parent&rarr;node) color when "Show support as: branch color" is on: the full branch
     * color for strong support, fading toward the background for weak support (theme-aware -- works on
     * both light and dark backgrounds). The branch's confidence and the tree's support scale drive the
     * fraction, exactly as the support symbols do.
     */
    private Color supportBranchColor(final PhylogenyNode node, final boolean to_pdf) {
        final double fraction = TreePanelUtil.supportFraction(PhylogenyMethods.getConfidenceValue(node),
                _confidence_scale_max);
        final Color strong = to_pdf ? getTreeColorSet().getBranchColorForPdf() : getTreeColorSet().getBranchColor();
        return TreePanelUtil.supportColor(fraction, strong, getTreeColorSet().getBackgroundColor());
    }

    final private void blast(final PhylogenyNode node) {
        if (!isCanBlast(node)) {
            JOptionPane.showMessageDialog(this,
                    "Insufficient information present",
                    "Cannot Blast",
                    JOptionPane.INFORMATION_MESSAGE);
            return;
        } else {
            final String query = Blast.obtainQueryForBlast(node);
            System.out.println("query for BLAST is: " + query);
            char type = '?';
            if (!ForesterUtil.isEmpty(query)) {
                if (node.getNodeData().isHasSequence()) {
                    if (!ForesterUtil.isEmpty(node.getNodeData().getSequence().getType())) {
                        if (node.getNodeData().getSequence().getType().toLowerCase()
                                .equals(PhyloXmlUtil.SEQ_TYPE_PROTEIN)) {
                            type = 'p';
                        } else {
                            type = 'n';
                        }
                    } else if (!ForesterUtil.isEmpty(node.getNodeData().getSequence().getMolecularSequence())) {
                        if (ForesterUtil
                                .seqIsLikelyToBeAa(node.getNodeData().getSequence().getMolecularSequence())) {
                            type = 'p';
                        } else {
                            type = 'n';
                        }
                    }
                }
                if (type == '?') {
                    if (SequenceAccessionTools.isProteinDbQuery(query)) {
                        type = 'p';
                    } else {
                        type = 'n';
                    }
                }
                try {
                    Blast.openNcbiBlastWeb(query, type == 'n', this);
                } catch (final Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private final int calcDynamicHidingFactor() {
        return (int) (0.5 + (getFontMetricsForLargeDefaultFont().getHeight() / (1.5 * getYdistance())));
    }

    final private int calcLengthOfLongestText() {
        if (_ext_node_with_longest_txt_info == null) {
            return 0;
        }
        final StringBuilder sb = new StringBuilder();
        nodeDataAsSB(_ext_node_with_longest_txt_info, sb);
        int sum = getFontMetricsForLargeDefaultFont().stringWidth(sb.toString());
        // Measure the taxonomy part italic-aware (the scientific name may be drawn in italics), so the
        // reserved width matches what is painted instead of a roman over/under-estimate.
        if (_ext_node_with_longest_txt_info.getNodeData().isHasTaxonomy()) {
            sum += taxonomyLabelWidth(_ext_node_with_longest_txt_info.getNodeData().getTaxonomy(),
                    getTreeFontSet().getLargeFont());
        }
        return sum;
    }

    /**
     * Calculate the length of the distance between the given node and its
     * parent.
     *
     * @param node
     * @return the distance value
     * @factor
     */
    final private float calculateBranchLengthToParent(final PhylogenyNode node, final float factor) {
        if (getControlPanel().isDrawPhylogram()) {
            if (node.getDistanceToParent() < 0.0) {
                return 0.0f;
            }
            return (float) (getXcorrectionFactor() * node.getDistanceToParent());
        } else {
            if ((factor == 0) || isNonLinedUpCladogram()) {
                return getXdistance();
            }
            return getXdistance() * factor;
        }
    }

    final private float calculateOvBranchLengthToParent(final PhylogenyNode node, final int factor) {
        if (getControlPanel().isDrawPhylogram()) {
            if (node.getDistanceToParent() < 0.0) {
                return 0.0f;
            }
            return (float) (getOvXcorrectionFactor() * node.getDistanceToParent());
        } else {
            if ((factor == 0) || isNonLinedUpCladogram()) {
                return getOvXDistance();
            }
            return getOvXDistance() * factor;
        }
    }

    final private void cannotOpenBrowserWarningMessage(final String type_type) {
        JOptionPane.showMessageDialog(this,
                "Cannot launch web browser for " + type_type + " data of this node",
                "Cannot launch web browser",
                JOptionPane.WARNING_MESSAGE);
    }

    final private void colorizeNodes(final Color c,
                                     final PhylogenyNode node,
                                     final List<PhylogenyNode> additional_nodes) {
        _control_panel.setColorBranches(true);
        if (_control_panel.getUseVisualStylesCb() != null) {
            _control_panel.getUseVisualStylesCb().setSelected(true);
        }
        if (node != null) {
            colorizeNodesHelper(c, node);
        }
        if (additional_nodes != null) {
            for (final PhylogenyNode n : additional_nodes) {
                colorizeNodesHelper(c, n);
            }
        }
        repaint();
    }

    final private void colorizeSubtree(final Color c,
                                       final PhylogenyNode node,
                                       final List<PhylogenyNode> additional_nodes) {
        _control_panel.setColorBranches(true);
        if (_control_panel.getUseVisualStylesCb() != null) {
            _control_panel.getUseVisualStylesCb().setSelected(true);
        }
        if (node != null) {
            for (final PreorderTreeIterator it = new PreorderTreeIterator(node); it.hasNext(); ) {
                it.next().getBranchData().setBranchColor(new BranchColor(c));
            }
        }
        if (additional_nodes != null) {
            for (final PhylogenyNode an : additional_nodes) {
                for (final PreorderTreeIterator it = new PreorderTreeIterator(an); it.hasNext(); ) {
                    it.next().getBranchData().setBranchColor(new BranchColor(c));
                }
            }
        }
        repaint();
    }

    private void colorNodeFont(final PhylogenyNode node) {
        _color_chooser.setPreviewPanel(new JPanel());
        NodeColorizationActionListener al;
        int count = 1;
        if ((getFoundNodes0() != null) || (getFoundNodes1() != null)) {
            final List<PhylogenyNode> additional_nodes = getFoundNodesAsListOfPhylogenyNodes();
            al = new NodeColorizationActionListener(_color_chooser, node, additional_nodes);
            count = additional_nodes.size();
            if (!additional_nodes.contains(node)) {
                count++;
            }
        } else {
            al = new NodeColorizationActionListener(_color_chooser, node);
        }
        String title = "Change the (node and font) color for ";
        if (count == 1) {
            title += "one node";
        } else {
            title += (count + " nodes");
        }
        final JDialog dialog = JColorChooser.createDialog(this, title, true, _color_chooser, al, null);
        setEdited(true);
        dialog.setVisible(true);
    }

    final private void colorSubtree(final PhylogenyNode node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            JOptionPane.showMessageDialog(this,
                    "Cannot colorize subtree in unrooted display type",
                    "Attempt to colorize subtree in unrooted display",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }
        _color_chooser.setPreviewPanel(new JPanel());
        final SubtreeColorizationActionListener al;
        final boolean color_found = getOptions().isColorAllFoundNodesWhenColoringSubtree();
        if (color_found && ((getFoundNodes0() != null) || (getFoundNodes1() != null))) {
            final List<PhylogenyNode> additional_nodes = getFoundNodesAsListOfPhylogenyNodes();
            al = new SubtreeColorizationActionListener(_color_chooser, node, additional_nodes);
        } else {
            al = new SubtreeColorizationActionListener(_color_chooser, node);
        }
        final JDialog dialog = JColorChooser
                .createDialog(this, "Subtree colorization", true, _color_chooser, al, null);
        setEdited(true);
        dialog.setVisible(true);
    }

    final private void copySubtree(final PhylogenyNode node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        setNodeInPreorderToNull();
        setCutOrCopiedTree(_phylogeny.copy(node));
        final List<PhylogenyNode> nodes = PhylogenyMethods.getAllDescendants(node);
        final Set<Long> node_ids = new HashSet<>(nodes.size());
        for (final PhylogenyNode n : nodes) {
            node_ids.add(n.getId());
        }
        node_ids.add(node.getId());
        setCopiedAndPastedNodes(node_ids);
        repaint();
    }

    final private String createASimpleTextRepresentationOfANode(final PhylogenyNode node) {
        final String tax = PhylogenyMethods.getSpecies(node);
        String label = node.getName();
        if (!ForesterUtil.isEmpty(label) && !ForesterUtil.isEmpty(tax)) {
            label = label + " " + tax;
        } else if (!ForesterUtil.isEmpty(tax)) {
            label = tax;
        } else {
            label = "";
        }
        if (!ForesterUtil.isEmpty(label)) {
            label = " [" + label + "]";
        }
        return label;
    }

    final private void cutSubtree(final PhylogenyNode node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if (node.isRoot()) {
            JOptionPane.showMessageDialog(this,
                    "Cannot cut entire tree as subtree",
                    "Attempt to cut entire tree",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        final String label = createASimpleTextRepresentationOfANode(node);
        final int r = JOptionPane.showConfirmDialog(null,
                "Cut subtree" + label + "?",
                "Confirm Cutting of Subtree",
                JOptionPane.YES_NO_OPTION);
        if (r != JOptionPane.OK_OPTION) {
            return;
        }
        pushUndoCheckpoint("Cut Subtree");
        setNodeInPreorderToNull();
        setCopiedAndPastedNodes(null);
        setCutOrCopiedTree(_phylogeny.copy(node));
        _phylogeny.deleteSubtree(node, true);
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants(true);
        resetNodeIdToDistToLeafMap();
        setEdited(true);
        repaint();
    }

    final private void decreaseOvSize() {
        if ((getOvMaxWidth() > 20) && (getOvMaxHeight() > 20)) {
            setOvMaxWidth(getOvMaxWidth() - 5);
            setOvMaxHeight(getOvMaxHeight() - 5);
            updateOvSettings();
            getControlPanel().displayedPhylogenyMightHaveChanged(false);
        }
    }

    final private void deleteNodeOrSubtree(final PhylogenyNode node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if (node.isRoot() && (node.getNumberOfDescendants() != 1)) {
            JOptionPane.showMessageDialog(this,
                    "Cannot delete entire tree",
                    "Attempt to delete entire tree",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        final String label = createASimpleTextRepresentationOfANode(node);
        final Object[] options = {"Node only", "Entire subtree", "Cancel"};
        final int r = JOptionPane.showOptionDialog(this,
                "Delete" + label + "?",
                "Delete Node/Subtree",
                JOptionPane.CLOSED_OPTION,
                JOptionPane.QUESTION_MESSAGE,
                null,
                options,
                options[2]);
        setNodeInPreorderToNull();
        boolean node_only = true;
        if (r == 1) {
            node_only = false;
        } else if (r != 0) {
            return;
        }
        pushUndoCheckpoint(node_only ? "Delete Node" : "Delete Subtree");
        if (node_only) {
            PhylogenyMethods.removeNode(node, _phylogeny);
        } else {
            _phylogeny.deleteSubtree(node, true);
        }
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants(true);
        resetNodeIdToDistToLeafMap();
        setEdited(true);
        repaint();
    }

    final private void displayNodePopupMenu(final PhylogenyNode node, final int x, final int y) {
        makePopupMenus(node);
        _node_popup_menu.putClientProperty(NODE_POPMENU_NODE_CLIENT_PROPERTY, node);
        _node_popup_menu.show(this, x, y);
    }

    final private void drawArc(final double x,
                               final double y,
                               final double width,
                               final double heigth,
                               final double start_angle,
                               final double arc_angle,
                               final Graphics2D g) {
        _arc.setArc(x, y, width, heigth, _180_OVER_PI * start_angle, _180_OVER_PI * arc_angle, Arc2D.OPEN);
        g.draw(_arc);
    }

    final private void drawLine(final double x1,
                                final double y1,
                                final double x2,
                                final double y2,
                                final Graphics2D g) {
        if ((x1 == x2) && (y1 == y2)) {
            return;
        }
        _line.setLine(x1, y1, x2, y2);
        g.draw(_line);
    }

    final private void drawOval(final double x,
                                final double y,
                                final double width,
                                final double heigth,
                                final Graphics2D g) {
        _ellipse.setFrame(x, y, width, heigth);
        g.draw(_ellipse);
    }

    final private void drawOvalFilled(final double x,
                                      final double y,
                                      final double width,
                                      final double heigth,
                                      final Graphics2D g) {
        _ellipse.setFrame(x, y, width, heigth);
        g.fill(_ellipse);
    }

    final private void drawOvalGradient(final float x,
                                        final float y,
                                        final float width,
                                        final float heigth,
                                        final Graphics2D g,
                                        final Color color_1,
                                        final Color color_2,
                                        final Color color_border) {
        _ellipse.setFrame(x, y, width, heigth);
        g.setPaint(new GradientPaint(x, y, color_1, (x + width), (y + heigth), color_2, false));
        g.fill(_ellipse);
        if (color_border != null) {
            g.setPaint(color_border);
            g.draw(_ellipse);
        }
    }

    final private void drawRect(final float x,
                                final float y,
                                final float width,
                                final float heigth,
                                final Graphics2D g) {
        _rectangle.setFrame(x, y, width, heigth);
        g.draw(_rectangle);
    }

    final private void drawRectFilled(final double x,
                                      final double y,
                                      final double width,
                                      final double heigth,
                                      final Graphics2D g) {
        _rectangle.setFrame(x, y, width, heigth);
        g.fill(_rectangle);
    }

    final private void drawRectGradient(final float x,
                                        final float y,
                                        final float width,
                                        final float heigth,
                                        final Graphics2D g,
                                        final Color color_1,
                                        final Color color_2,
                                        final Color color_border) {
        _rectangle.setFrame(x, y, width, heigth);
        g.setPaint(new GradientPaint(x, y, color_1, (x + width), (y + heigth), color_2, false));
        g.fill(_rectangle);
        if (color_border != null) {
            g.setPaint(color_border);
            g.draw(_rectangle);
        }
    }



    final private void errorMessageNoCutCopyPasteInUnrootedDisplay() {
        JOptionPane.showMessageDialog(this,
                "Cannot cut, copy, paste, add, or delete subtrees/nodes in unrooted display",
                "Attempt to cut/copy/paste/add/delete in unrooted display",
                JOptionPane.ERROR_MESSAGE);
    }

    private final Color getColorForFoundNode(final PhylogenyNode n) {
        if (isInFoundNodes0(n) && !isInFoundNodes1(n)) {
            return getTreeColorSet().getFoundColor0();
        } else if (!isInFoundNodes0(n) && isInFoundNodes1(n)) {
            return getTreeColorSet().getFoundColor1();
        } else {
            return getTreeColorSet().getFoundColor0and1();
        }
    }

    final private Set<Long> getCopiedAndPastedNodes() {
        return getMainPanel().getCopiedAndPastedNodes();
    }

    final private Phylogeny getCutOrCopiedTree() {
        return getMainPanel().getCutOrCopiedTree();
    }

    private FontMetrics getFontMetricsForLargeDefaultFont() {
        return getTreeFontSet().getFontMetricsLarge();
    }

    /** Pixel height of the (large) font used for leaf labels. */
    final int getLargeFontHeight() {
        return getFontMetricsForLargeDefaultFont().getHeight();
    }

    final private float getLastDragPointX() {
        return _last_drag_point_x;
    }

    final private float getLastDragPointY() {
        return _last_drag_point_y;
    }

    final private double getMaxDistanceToRoot() {
        if (_max_distance_to_root < 0) {
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

    final private List<Accession> getPdbAccs(final PhylogenyNode node) {
        final List<Accession> pdb_ids = new ArrayList<>();
        if (node.getNodeData().isHasSequence()) {
            final Sequence seq = node.getNodeData().getSequence();
            if (!ForesterUtil.isEmpty(seq.getCrossReferences())) {
                final SortedSet<Accession> cross_refs = seq.getCrossReferences();
                for (final Accession acc : cross_refs) {
                    if (acc.getSource().equalsIgnoreCase("pdb")) {
                        pdb_ids.add(acc);
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

    final private void handleClickToAction(final NodeClickAction action, final PhylogenyNode node) {
        switch (action) {
            case SHOW_DATA:
                showNodeFrame(node);
                break;
            case COLLAPSE:
                collapse(node);
                break;
            case REROOT:
                reRoot(node);
                break;
            case SUBTREE:
                subTree(node);
                break;
            case SWAP:
                swap(node);
                break;
            case COLOR_SUBTREE:
                colorSubtree(node);
                break;
            case COLOR_NODE_FONT:
                colorNodeFont(node);
                break;
            case OPEN_SEQ_WEB:
                openSeqWeb(node);
                break;
            case BLAST:
                blast(node);
                break;
            case OPEN_TAX_WEB:
                openTaxWeb(node);
                break;
            case OPEN_PDB_WEB:
                openPdbWeb(node);
                break;
            case CUT_SUBTREE:
                cutSubtree(node);
                break;
            case COPY_SUBTREE:
                copySubtree(node);
                break;
            case PASTE_SUBTREE:
                pasteSubtree(node);
                break;
            case DELETE_NODE_OR_SUBTREE:
                deleteNodeOrSubtree(node);
                break;
            case ADD_NEW_NODE:
                addEmptyNode(node);
                break;
            case EDIT_NODE_DATA:
                showNodeEditFrame(node);
                break;
            case SELECT_NODES:
                selectNode(node);
                break;
            case SORT_DESCENDENTS:
                sortDescendants(node);
                break;
            case UNCOLLAPSE_ALL:
                uncollapseAll(node);
                break;
            case ORDER_SUBTREE:
                orderSubtree(node);
                break;
            default:
                throw new IllegalArgumentException("unknown action: " + action);
        }
    }

    final private void increaseOvSize() {
        if ((getOvMaxWidth() < (getMainPanel().getCurrentScrollPane().getViewport().getVisibleRect().getWidth()
                / 2))
                && (getOvMaxHeight() < (getMainPanel().getCurrentScrollPane().getViewport().getVisibleRect()
                .getHeight() / 2))) {
            setOvMaxWidth(getOvMaxWidth() + 5);
            setOvMaxHeight(getOvMaxHeight() + 5);
            updateOvSettings();
            getControlPanel().displayedPhylogenyMightHaveChanged(false);
        }
    }

    final private void init() {
        _color_chooser = new JColorChooser();
        _rollover_popup = new JTextArea();
        _rollover_popup.setFont(POPUP_FONT);
        resetNodeIdToDistToLeafMap();
        setTextAntialias();
        setTreeFile(null);
        setEdited(false);
        initializeOvSettings();
        setStartingAngle((TWO_PI * 3) / 4);
    }

    final private void initializeOvSettings() {
        setOvMaxHeight(getConfiguration().getOvMaxHeight());
        setOvMaxWidth(getConfiguration().getOvMaxWidth());
    }

    final private boolean inOvVirtualRectangle(final int x, final int y) {
        return ((x >= (getOvVirtualRectangle().x - 1))
                && (x <= (getOvVirtualRectangle().x + getOvVirtualRectangle().width + 1))
                && (y >= (getOvVirtualRectangle().y - 1))
                && (y <= (getOvVirtualRectangle().y + getOvVirtualRectangle().height + 1)));
    }

    final private boolean inOvVirtualRectangle(final MouseEvent e) {
        return (inOvVirtualRectangle(e.getX(), e.getY()));
    }

    final private boolean isCanBlast(final PhylogenyNode node) {
        if (!node.getNodeData().isHasSequence() && ForesterUtil.isEmpty(node.getName())) {
            return false;
        }
        return Blast.isContainsQueryForBlast(node);
    }

    final private String isCanOpenSeqWeb(final PhylogenyNode node) {
        final Accession a = SequenceAccessionTools.obtainAccessorFromDataFields(node);
        if (a != null) {
            return a.getValue();
        }
        return null;
    }

    final private boolean isCanOpenTaxWeb(final PhylogenyNode node) {
        if (node.getNodeData().isHasTaxonomy() && ((!ForesterUtil
                .isEmpty(node.getNodeData().getTaxonomy().getScientificName()))
                || (!ForesterUtil.isEmpty(node.getNodeData().getTaxonomy().getTaxonomyCode()))
                || (!ForesterUtil.isEmpty(node.getNodeData().getTaxonomy().getCommonName()))
                || ((node.getNodeData().getTaxonomy().getIdentifier() != null)
                && !ForesterUtil.isEmpty(node.getNodeData().getTaxonomy().getIdentifier().getValue())))) {
            return true;
        } else {
            return false;
        }
    }

    private boolean isInFoundNodes(final PhylogenyNode n) {
        return isInFoundNodes0(n) || isInFoundNodes1(n);
    }

    /** Whether a search/selection is currently active (either found set has at least one node). Deliberately
     *  STRICTER than the bare {@code getFoundNodes0() != null} idiom used elsewhere: a 0-hit search can leave a
     *  non-null empty set, and dimming must NOT engage then. */
    private boolean hasFoundNodes() {
        return ((getFoundNodes0() != null) && !getFoundNodes0().isEmpty())
                || ((getFoundNodes1() != null) && !getFoundNodes1().isEmpty());
    }

    /** True iff at least one found node is currently DRAWN (in the displayed (sub)tree and not hidden under a
     *  collapsed ancestor) -- the gate for the "Dim Non-Matches" wash, so the whole tree never fades with no hit
     *  on screen (all hits collapsed, or navigated into a subtree that contains none). Early-terminates. */
    private boolean anyVisibleFoundNode() {
        if (!hasFoundNodes()) {
            return false;
        }
        for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            // cheap found-set lookup first, so the O(depth) hidden-under-collapse walk runs only for actual hits
            if (isInFoundNodes(n) && !isHiddenUnderCollapse(n)) {
                return true;
            }
        }
        return false;
    }

    /** All found (search box A + B, and manual selection) nodes in tree (preorder) order -- the step-through list. */
    private List<PhylogenyNode> orderedFoundNodes() {
        final List<PhylogenyNode> hits = new ArrayList<>();
        if ((_phylogeny == null) || !hasFoundNodes()) {
            return hits;
        }
        for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if (isInFoundNodes(n)) {
                hits.add(n);
            }
        }
        return hits;
    }

    /** Centers the NEXT ({@code dir}=+1) or PREVIOUS ({@code dir}=-1) search/selection hit in the viewport, wrapping
     *  around. A hit hidden under a collapse scrolls to its (drawn) collapsed-clade triangle instead. Updates the
     *  left panel's "k / N" navigator. Called by the panel's arrow buttons and the View-menu Find Next/Previous. */
    final void stepToFoundNode(final int dir) {
        final List<PhylogenyNode> hits = orderedFoundNodes();
        if (hits.isEmpty()) {
            _search_hit_index = -1;
            _last_step_target = null;
            getControlPanel().updateSearchHitNavigation();
            return;
        }
        int idx = _search_hit_index;
        if (idx < 0) {
            idx = (dir >= 0) ? 0 : (hits.size() - 1); // first step lands on the first (or last) hit
        } else {
            idx = ((((idx + dir) % hits.size()) + hits.size()) % hits.size()); // wrap in both directions
        }
        _search_hit_index = idx;
        final PhylogenyNode hit = hits.get(idx);
        _last_step_target = isHiddenUnderCollapse(hit) ? outermostCollapsedAncestor(hit) : hit;
        centerOnNode(_last_step_target);
        getControlPanel().updateSearchHitNavigation();
    }

    /** The outermost collapsed ancestor of {@code n} (the drawn triangle a hidden hit belongs to); {@code n} itself
     *  if it is not hidden under any collapse. */
    private static PhylogenyNode outermostCollapsedAncestor(final PhylogenyNode n) {
        PhylogenyNode result = n;
        for (PhylogenyNode p = n.getParent(); p != null; p = p.getParent()) {
            if (p.isCollapse()) {
                result = p;
            }
        }
        return result;
    }

    /** Scrolls the tree's viewport so {@code node} is centered (clamped to the scrollable extent). */
    private void centerOnNode(final PhylogenyNode node) {
        final JScrollPane sp = getMainPanel().getCurrentScrollPane();
        if (sp == null) {
            return;
        }
        final Rectangle view = sp.getViewport().getViewRect();
        // use the node's ON-SCREEN position (rotated by R in a vertical orientation), not its logical coords
        final Point2D.Double screen = screenPointFor(node);
        int x = (int) Math.round(screen.x) - (view.width / 2);
        int y = (int) Math.round(screen.y) - (view.height / 2);
        x = Math.max(0, Math.min(x, Math.max(0, getWidth() - view.width)));
        y = Math.max(0, Math.min(y, Math.max(0, getHeight() - view.height)));
        sp.getViewport().setViewPosition(new Point(x, y));
        repaint();
    }

    /** Number of search/selection hits (for the "k / N" navigator). */
    final int getSearchHitCount() {
        return orderedFoundNodes().size();
    }

    /** Current 0-based step-through position, or -1 if not positioned yet. */
    final int getSearchHitIndex() {
        return _search_hit_index;
    }

    /** Test hook: the node last centered on by {@link #stepToFoundNode(int)}. */
    PhylogenyNode getLastStepTargetForTest() {
        return _last_step_target;
    }

    /** Test hook: center the viewport on a node (drives the orientation-aware {@link #centerOnNode}). */
    void centerOnNodeForTest(final PhylogenyNode node) {
        centerOnNode(node);
    }

    final private boolean isInFoundNodes0(final PhylogenyNode node) {
        return ((getFoundNodes0() != null) && getFoundNodes0().contains(node.getId()));
    }

    final private boolean isInFoundNodes1(final PhylogenyNode node) {
        return ((getFoundNodes1() != null) && getFoundNodes1().contains(node.getId()));
    }

    final private boolean isInOv() {
        return _in_ov;
    }

    final private boolean isNodeDataInvisible(final PhylogenyNode node) {
        // work entirely in LOGICAL space (node coords are logical): map the viewport back through R once, then cull the
        // node's incoming-branch elbow bounding box, grown by the label reach so a node whose (possibly tilted) label
        // pokes into the viewport is never culled. Correct in EVERY orientation (vertical included -- see the disabled
        // per-orientation special-case this replaced).
        final Rectangle2D vis = logicalVisibleRect();
        if (vis == null) {
            return false; // no meaningful viewport -> draw everything (correct, just not thrifty)
        }
        double minx = node.getXcoord(), maxx = node.getXcoord(), miny = node.getYcoord(), maxy = node.getYcoord();
        if (node.getParent() != null) {
            minx = Math.min(minx, node.getParent().getXcoord());
            maxx = Math.max(maxx, node.getParent().getXcoord());
            miny = Math.min(miny, node.getParent().getYcoord());
            maxy = Math.max(maxy, node.getParent().getYcoord());
        }
        final double margin = 40 + Math.ceil(getLongestExtNodeInfo() / Math.sqrt(2.0));
        return ((maxx < (vis.getMinX() - margin)) || (minx > (vis.getMaxX() + margin))
                || (maxy < (vis.getMinY() - margin)) || (miny > (vis.getMaxY() + margin)));
    }

    final private boolean isNodeDataInvisibleUnrootedCirc(final PhylogenyNode node) {
        return ((node.getYcoord() < (getVisibleRect().getMinY() - 20))
                || (node.getYcoord() > (getVisibleRect().getMaxY() + 20))
                || (node.getXcoord() < (getVisibleRect().getMinX() - 20))
                || (node.getXcoord() > (getVisibleRect().getMaxX() + 20)));
    }

    final private boolean isNonLinedUpCladogram() {
        return getOptions().getCladogramType() == CLADOGRAM_TYPE.NON_LINED_UP;
    }

    final private void keyPressedCalls(final KeyEvent e) {
        // Keys carrying the platform menu-shortcut modifier (Cmd on macOS, Ctrl elsewhere) belong to the menu
        // accelerators (Save, Open, Close Tab, Copy Image, Undo/Redo, ...). This focused-canvas KeyListener must
        // NOT handle or consume them, or it would shadow the accelerator entirely (Swing runs a focused
        // component's KeyListeners before its WHEN_IN_FOCUSED_WINDOW menu accelerators, and skips the latter once
        // the event is consumed) -- and, since the single-letter branches below match by key code alone, it would
        // read e.g. Cmd+S as a bare 'S' (rotate) or Cmd+O as 'O' (cycle overview).
        if ((e.getModifiersEx() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx()) != 0) {
            return;
        }
        if (isOvOn() && (getMousePosition() != null) && (getMousePosition().getLocation() != null)) {
            if (inOvVirtualRectangle(getMousePosition().x, getMousePosition().y)) {
                if (!isInOvRect()) {
                    setInOvRect(true);
                }
            } else if (isInOvRect()) {
                setInOvRect(false);
            }
        }
        // Only keys this handler actually acts on are consumed; anything else propagates (so e.g. focus
        // traversal and other bindings aren't silently swallowed by the focus-grabbing canvas).
        boolean handled = false;
        if (e.isAltDown()) {
            if ((e.getKeyCode() == KeyEvent.VK_DELETE) || (e.getKeyCode() == KeyEvent.VK_HOME)
                    || (e.getKeyCode() == KeyEvent.VK_C) || (e.getKeyCode() == KeyEvent.VK_BACK_SPACE)) {
                getControlPanel().showWhole();
                handled = true;
            } else if (e.isShiftDown()
                    && ((e.getKeyCode() == KeyEvent.VK_SUBTRACT) || (e.getKeyCode() == KeyEvent.VK_MINUS))) {
                getMainPanel().getTreeFontSet().decreaseUserFontSize();
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(true);
                handled = true;
            } else if (e.isShiftDown() && plusPressed(e.getKeyCode())) {
                getMainPanel().getTreeFontSet().increaseUserFontSize();
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(true);
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_O) {
                getControlPanel().orderPressed(this);
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_R) {
                if (e.isShiftDown()) {
                    getControlPanel().returnedToWholeTreePressed();
                } else {
                    getControlPanel().returnedToSuperTreePressed();
                }
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_U) {
                getControlPanel().uncollapseAll(this);
                getControlPanel().displayedPhylogenyMightHaveChanged(false);
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_E) {
                getControlPanel().expandYToFitLabels();
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_W) {
                if (isVerticalOrientation()) {
                    getControlPanel().fitHeight();
                } else {
                    getControlPanel().fitWidth();
                }
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_UP) {
                getMainPanel().getControlPanel().zoomInScreenY(AptxConstants.WHEEL_ZOOM_IN_FACTOR,
                        AptxConstants.WHEEL_ZOOM_IN_X_CORRECTION_FACTOR);
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(false);
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_DOWN) {
                getMainPanel().getControlPanel().zoomOutScreenY(AptxConstants.WHEEL_ZOOM_OUT_FACTOR,
                        AptxConstants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR);
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(false);
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_LEFT) {
                getMainPanel().getControlPanel().zoomOutScreenX(AptxConstants.WHEEL_ZOOM_OUT_FACTOR,
                        AptxConstants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR);
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(false);
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
                getMainPanel().getControlPanel().zoomInScreenX(AptxConstants.WHEEL_ZOOM_IN_FACTOR,
                        AptxConstants.WHEEL_ZOOM_IN_FACTOR);
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(false);
                handled = true;
            } else if ((e.getKeyCode() == KeyEvent.VK_SUBTRACT) || (e.getKeyCode() == KeyEvent.VK_MINUS)) {
                getMainPanel().getControlPanel().zoomOutY(AptxConstants.WHEEL_ZOOM_OUT_FACTOR);
                getMainPanel().getControlPanel().zoomOutX(AptxConstants.WHEEL_ZOOM_OUT_FACTOR,
                        AptxConstants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR);
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(false);
                handled = true;
            } else if (plusPressed(e.getKeyCode())) {
                getMainPanel().getControlPanel().zoomInX(AptxConstants.WHEEL_ZOOM_IN_FACTOR,
                        AptxConstants.WHEEL_ZOOM_IN_FACTOR);
                getMainPanel().getControlPanel().zoomInY(AptxConstants.WHEEL_ZOOM_IN_FACTOR);
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(false);
                handled = true;
            }
        } else {
            if ((e.getKeyCode() == KeyEvent.VK_UP) || (e.getKeyCode() == KeyEvent.VK_DOWN)
                    || (e.getKeyCode() == KeyEvent.VK_LEFT) || (e.getKeyCode() == KeyEvent.VK_RIGHT)) {
                final int d = 80;
                int dx = 0;
                int dy = -d;
                if (e.getKeyCode() == KeyEvent.VK_DOWN) {
                    dy = d;
                } else if (e.getKeyCode() == KeyEvent.VK_LEFT) {
                    dx = -d;
                    dy = 0;
                } else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
                    dx = d;
                    dy = 0;
                }
                final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
                scroll_position.x = scroll_position.x + dx;
                scroll_position.y = scroll_position.y + dy;
                if (scroll_position.x <= 0) {
                    scroll_position.x = 0;
                } else {
                    final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                            - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
                    if (scroll_position.x >= max_x) {
                        scroll_position.x = max_x;
                    }
                }
                if (scroll_position.y <= 0) {
                    scroll_position.y = 0;
                } else {
                    final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                            - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
                    if (scroll_position.y >= max_y) {
                        scroll_position.y = max_y;
                    }
                }
                repaint();
                getMainPanel().getCurrentScrollPane().getViewport().setViewPosition(scroll_position);
                handled = true;
            } else if (e.getKeyCode() == KeyEvent.VK_S) {
                if ((getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                        || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
                    setStartingAngle((getStartingAngle() % TWO_PI) + ANGLE_ROTATION_UNIT);
                    getControlPanel().displayedPhylogenyMightHaveChanged(false);
                    handled = true;
                }
            } else if (e.getKeyCode() == KeyEvent.VK_A) {
                if ((getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                        || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
                    setStartingAngle((getStartingAngle() % TWO_PI) - ANGLE_ROTATION_UNIT);
                    if (getStartingAngle() < 0) {
                        setStartingAngle(TWO_PI + getStartingAngle());
                    }
                    getControlPanel().displayedPhylogenyMightHaveChanged(false);
                    handled = true;
                }
            } else if (getOptions().isShowOverview() && isOvOn() && (e.getKeyCode() == KeyEvent.VK_O)) {
                MainFrame.cycleOverview(getOptions(), this);
                repaint();
                handled = true;
            } else if (getOptions().isShowOverview() && isOvOn() && (e.getKeyCode() == KeyEvent.VK_I)) {
                increaseOvSize();
                handled = true;
            } else if (getOptions().isShowOverview() && isOvOn() && (e.getKeyCode() == KeyEvent.VK_U)) {
                decreaseOvSize();
                handled = true;
            }
        }
        if ((e.getKeyCode() == KeyEvent.VK_HOME) || (e.getKeyCode() == KeyEvent.VK_ESCAPE)) {
            getControlPanel().showWhole();
            handled = true;
        } else if (e.getKeyCode() == KeyEvent.VK_PAGE_UP) {
            getMainPanel().getTreeFontSet().increaseUserFontSize();
            getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(true);
            handled = true;
        } else if (e.getKeyCode() == KeyEvent.VK_PAGE_DOWN) {
            getMainPanel().getTreeFontSet().decreaseUserFontSize();
            getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(true);
            handled = true;
        }
        if (handled) {
            e.consume();
        }
    }

    final private void makePopupMenus(final PhylogenyNode node) {
        _node_popup_menu = new JPopupMenu();
        final List<String> clickto_names = _main_panel.getControlPanel().getSingleClickToNames();
        _node_popup_menu_items = new JMenuItem[clickto_names.size()];
        for (int i = 0; i < clickto_names.size(); i++) {
            final String title = clickto_names.get(i);
            _node_popup_menu_items[i] = new JMenuItem(title);
            if (title.equals(Configuration.clickto_options[Configuration.open_seq_web][0])) {
                final List<Sequence> seqs = node.getNodeData().getSequences();
                boolean addedPerSeqItems = false;
                if (seqs != null) {
                    for (final Sequence seq : seqs) {
                        final Accession acc = SequenceAccessionTools.obtainAccessorFromSequence(seq);
                        if (acc != null) {
                            final String seqName = seq.getName();
                            final String trimmedName = (!ForesterUtil.isEmpty(seqName) && seqName.length() > 20)
                                    ? seqName.substring(0, 20) + "..."
                                    : seqName;
                            final String label = ForesterUtil.isEmpty(trimmedName)
                                    ? acc.getValue()
                                    : trimmedName + " [" + acc.getValue() + "]";
                            final JMenuItem seqItem = new JMenuItem(title + " [" + label + "]");
                            final Sequence capturedSeq = seq;
                            seqItem.addActionListener(new ActionListener() {
                                @Override
                                public void actionPerformed(final ActionEvent e) {
                                    openSeqWebForSequence(node, capturedSeq);
                                }
                            });
                            _node_popup_menu.add(seqItem);
                            addedPerSeqItems = true;
                        }
                    }
                }
                if (addedPerSeqItems) {
                    continue;
                }
                final String id = isCanOpenSeqWeb(node);
                if (!ForesterUtil.isEmpty(id)) {
                    _node_popup_menu_items[i].setText(_node_popup_menu_items[i].getText() + " [" + id + "]");
                    _node_popup_menu_items[i].setEnabled(true);
                } else {
                    _node_popup_menu_items[i].setEnabled(false);
                }
            } else if (title.equals(Configuration.clickto_options[Configuration.open_pdb_web][0])) {
                final List<Accession> accs = getPdbAccs(node);
                _node_popup_menu_items[i] = new JMenuItem(title);
                if (!ForesterUtil.isEmpty(accs)) {
                    if (accs.size() == 1) {
                        _node_popup_menu_items[i].setText(_node_popup_menu_items[i].getText() + " ["
                                + TreePanelUtil.pdbAccToString(accs, 0) + "]");
                        _node_popup_menu_items[i].setEnabled(true);
                    } else if (accs.size() == 2) {
                        _node_popup_menu_items[i].setText(_node_popup_menu_items[i].getText() + " ["
                                + TreePanelUtil.pdbAccToString(accs, 0) + ", "
                                + TreePanelUtil.pdbAccToString(accs, 1) + "]");
                        _node_popup_menu_items[i].setEnabled(true);
                    } else if (accs.size() == 3) {
                        _node_popup_menu_items[i].setText(_node_popup_menu_items[i].getText() + " ["
                                + TreePanelUtil.pdbAccToString(accs, 0) + ", "
                                + TreePanelUtil.pdbAccToString(accs, 1) + ", "
                                + TreePanelUtil.pdbAccToString(accs, 2) + "]");
                        _node_popup_menu_items[i].setEnabled(true);
                    } else {
                        _node_popup_menu_items[i].setText(_node_popup_menu_items[i].getText() + " ["
                                + TreePanelUtil.pdbAccToString(accs, 0) + ", "
                                + TreePanelUtil.pdbAccToString(accs, 1) + ", "
                                + TreePanelUtil.pdbAccToString(accs, 2) + ", + " + (accs.size() - 3) + " more]");
                        _node_popup_menu_items[i].setEnabled(true);
                    }
                } else {
                    _node_popup_menu_items[i].setEnabled(false);
                }
            } else if (title.equals(Configuration.clickto_options[Configuration.open_tax_web][0])) {
                _node_popup_menu_items[i].setEnabled(isCanOpenTaxWeb(node));
            } else if (title.equals(Configuration.clickto_options[Configuration.blast][0])) {
                _node_popup_menu_items[i].setEnabled(isCanBlast(node));
            } else if (title.equals(Configuration.clickto_options[Configuration.delete_subtree_or_node][0])) {
                if (!getOptions().isEditable()) {
                    continue;
                }
                _node_popup_menu_items[i].setEnabled(isCanDelete());
            } else if (title.equals(Configuration.clickto_options[Configuration.cut_subtree][0])) {
                if (!getOptions().isEditable()) {
                    continue;
                }
                _node_popup_menu_items[i].setEnabled(isCanCut(node));
            } else if (title.equals(Configuration.clickto_options[Configuration.copy_subtree][0])) {
                if (!getOptions().isEditable()) {
                    continue;
                }
                _node_popup_menu_items[i].setEnabled(isCanCopy());
            } else if (title.equals(Configuration.clickto_options[Configuration.paste_subtree][0])) {
                if (!getOptions().isEditable()) {
                    continue;
                }
                _node_popup_menu_items[i].setEnabled(isCanPaste());
            } else if (title.equals(Configuration.clickto_options[Configuration.edit_node_data][0])) {
                if (!getOptions().isEditable()) {
                    continue;
                }
            } else if (title.equals(Configuration.clickto_options[Configuration.add_new_node][0])) {
                if (!getOptions().isEditable()) {
                    continue;
                }
            } else if (title.equals(Configuration.clickto_options[Configuration.reroot][0])) {
                _node_popup_menu_items[i].setEnabled(isCanReroot());
            } else if (title.equals(Configuration.clickto_options[Configuration.collapse_uncollapse][0])) {
                _node_popup_menu_items[i].setEnabled((isCanCollapse() && !node.isExternal()));
            } else if (title.equals(Configuration.clickto_options[Configuration.color_subtree][0])) {
                _node_popup_menu_items[i].setEnabled(isCanColorSubtree());
            } else if (title.equals(Configuration.clickto_options[Configuration.subtree][0])) {
                _node_popup_menu_items[i].setEnabled(isCanSubtree(node));
            } else if (title.equals(Configuration.clickto_options[Configuration.swap][0])) {
                _node_popup_menu_items[i].setEnabled(node.getNumberOfDescendants() == 2);
            } else if (title.equals(Configuration.clickto_options[Configuration.sort_descendents][0])) {
                _node_popup_menu_items[i].setEnabled(node.getNumberOfDescendants() > 1);
            } else if (title.equals(Configuration.clickto_options[Configuration.uncollapse_all][0])) {
                _node_popup_menu_items[i].setEnabled(isCanUncollapseAll(node));
            }
            _node_popup_menu_items[i].addActionListener(this);
            _node_popup_menu.add(_node_popup_menu_items[i]);
        }
    }

    private final void nodeDataAsSB(final PhylogenyNode node, final StringBuilder sb) {
        if (node != null) {
            if (getControlPanel().isShowNodeNames() && (!ForesterUtil.isEmpty(node.getName()))) {
                if (sb.length() > 0) {
                    sb.append(" ");
                }
                // Display-only shortening of an over-long name (e.g. a whole UniProt/NCBI header): the
                // node's actual name is left intact, so export / Find / accession parsing keep the full text.
                sb.append(getControlPanel().isShortenLabels()
                        ? AptxUtil.shortenLabel(node.getName(), AptxConstants.LONG_NODE_NAME_LIMIT)
                        : node.getName());
            }
            if (node.getNodeData().isHasSequence()
                    && (getControlPanel().isShowSeqSymbols()
                    || getControlPanel().isShowGeneNames()
                    || getControlPanel().isShowSeqNames()
                    || getControlPanel().isShowSequenceAcc()
            )) {
                final int s = node.getNodeData().getSequences().size();
                for (int i = 0; i < s; ++i) {
                    final Sequence seq = node.getNodeData().getSequence(i);
                    if (seq != null) {
                        if (s > 1) {
                            if (i > 0) {
                                sb.append("; [");
                            } else {
                                sb.append(" [");
                            }
                            sb.append(i);
                            sb.append("/");
                            sb.append(s - 1);
                            sb.append("]");
                        }
                        if (getControlPanel().isShowSeqSymbols()
                                && (seq.getSymbol().length() > 0)) {
                            if (sb.length() > 0) {
                                sb.append(" ");
                            }
                            sb.append(seq.getSymbol());
                        }
                        if (getControlPanel().isShowGeneNames()
                                && (seq.getGeneName().length() > 0)) {
                            if (sb.length() > 0) {
                                sb.append(" ");
                            }
                            sb.append(seq.getGeneName());
                        }
                        if (getControlPanel().isShowSeqNames()
                                && (seq.getName().length() > 0)) {
                            if (sb.length() > 0) {
                                sb.append(" ");
                            }
                            sb.append(seq.getName());
                        }
                        if (getControlPanel().isShowSequenceAcc()
                                && (seq.getAccession() != null)) {
                            if (sb.length() > 0) {
                                sb.append(" ");
                            }
                            if (!ForesterUtil.isEmpty(seq.getAccession().getSource())) {
                                sb.append(seq.getAccession().getSource());
                                sb.append(":");
                            }
                            sb.append(seq.getAccession().getValue());
                        }
                    }
                }
            }
            if (getControlPanel().isShowProperties() && node.getNodeData().isHasProperties()) {
                if (sb.length() > 0) {
                    sb.append(" ");
                }
                sb.append(propertiesToString(node));
            }
        }
    }

    private final void nodeTaxonomyDataAsSB(final Taxonomy taxonomy, final StringBuilder sb) {
        forEachTaxonomyLabelPart(taxonomy, (text, scientific) -> sb.append(text));
    }

    /** Receives the ordered parts of a taxonomy label; {@code scientific} marks the scientific-name part. */
    private interface TaxonomyLabelPartVisitor {

        void visit(String text, boolean scientific);
    }

    /**
     * Walks a taxonomy label in display order, emitting it one part at a time so callers can render or
     * measure each part differently -- the single place the rank / code / scientific-name / common-name
     * assembly lives. The scientific-name part is flagged so it can be italicized; the rank "[..]", the
     * code, the common name and its parentheses and all separators are non-scientific. The concatenation
     * of the emitted parts is exactly the string {@link #nodeTaxonomyDataAsSB} produces.
     */
    private void forEachTaxonomyLabelPart(final Taxonomy taxonomy, final TaxonomyLabelPartVisitor v) {
        if (_control_panel.isShowTaxonomyRank() && !ForesterUtil.isEmpty(taxonomy.getRank())) {
            v.visit("[" + taxonomy.getRank() + "] ", false);
        }
        if (_control_panel.isShowTaxonomyCode() && !ForesterUtil.isEmpty(taxonomy.getTaxonomyCode())) {
            v.visit(taxonomy.getTaxonomyCode() + " ", false);
        }
        final boolean show_sci = _control_panel.isShowTaxonomyScientificNames();
        final boolean show_common = _control_panel.isShowTaxonomyCommonNames();
        final boolean has_sci = !ForesterUtil.isEmpty(taxonomy.getScientificName());
        final boolean has_common = !ForesterUtil.isEmpty(taxonomy.getCommonName());
        if (show_sci && show_common) {
            if (has_sci && has_common) {
                v.visit(scientificNameForDisplay(taxonomy), true);
                v.visit(" (" + taxonomy.getCommonName() + ") ", false);
            } else if (has_sci) {
                v.visit(scientificNameForDisplay(taxonomy), true);
                v.visit(" ", false);
            } else if (has_common) {
                v.visit(taxonomy.getCommonName() + " ", false);
            }
        } else if (show_sci) {
            if (has_sci) {
                v.visit(scientificNameForDisplay(taxonomy), true);
                v.visit(" ", false);
            }
        } else if (show_common) {
            if (has_common) {
                v.visit(taxonomy.getCommonName() + " ", false);
            }
        }
    }

    /**
     * The scientific name as displayed: abbreviated ("H. sapiens") when that option is on and the name is
     * binomial (has a space), else verbatim.
     */
    private String scientificNameForDisplay(final Taxonomy taxonomy) {
        final String sn = taxonomy.getScientificName();
        if (getOptions().isAbbreviateScientificTaxonNames() && (sn.indexOf(' ') > 0)) {
            return TreePanelUtil.abbreviateScientificName(sn);
        }
        return sn;
    }

    final private void openPdbWeb(final PhylogenyNode node) {
        final List<Accession> pdb_ids = getPdbAccs(node);
        if (ForesterUtil.isEmpty(pdb_ids)) {
            cannotOpenBrowserWarningMessage("PDB");
            return;
        }
        final List<String> uri_strs = TreePanelUtil.createUrisForPdbWeb(node, pdb_ids, getConfiguration(), this);
        if (!ForesterUtil.isEmpty(uri_strs)) {
            for (final String uri_str : uri_strs) {
                try {
                    AptxUtil.launchWebBrowser(new URI(uri_str), "_aptx_seq");
                } catch (final IOException e) {
                    AptxUtil.showErrorMessage(this, e.toString());
                    e.printStackTrace();
                } catch (final URISyntaxException e) {
                    AptxUtil.showErrorMessage(this, e.toString());
                    e.printStackTrace();
                }
            }
        } else {
            cannotOpenBrowserWarningMessage("PDB");
        }
    }

    final private void openSeqWeb(final PhylogenyNode node) {
        if (ForesterUtil.isEmpty(isCanOpenSeqWeb(node))) {
            cannotOpenBrowserWarningMessage("sequence");
            return;
        }
        final String uri_str = TreePanelUtil.createUriForSeqWeb(node, getConfiguration(), this);
        if (!ForesterUtil.isEmpty(uri_str)) {
            try {
                AptxUtil.launchWebBrowser(new URI(uri_str), "_aptx_seq");
            } catch (final IOException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            } catch (final URISyntaxException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            }
        } else {
            cannotOpenBrowserWarningMessage("sequence");
        }
    }

    final private void openSeqWebForSequence(final PhylogenyNode node, final Sequence seq) {
        final String uri_str = TreePanelUtil.createUriForSeqWeb(seq, getConfiguration(), this);
        if (!ForesterUtil.isEmpty(uri_str)) {
            try {
                AptxUtil.launchWebBrowser(new URI(uri_str), "_aptx_seq");
            } catch (final IOException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            } catch (final URISyntaxException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            }
        } else {
            cannotOpenBrowserWarningMessage("sequence");
        }
    }

    final private void openTaxWeb(final PhylogenyNode node) {
        if (!isCanOpenTaxWeb(node)) {
            cannotOpenBrowserWarningMessage("taxonomic");
            return;
        }
        String uri_str = null;
        final Taxonomy tax = node.getNodeData().getTaxonomy();
        if ((tax.getIdentifier() != null) && !ForesterUtil.isEmpty(tax.getIdentifier().getValue())
                && tax.getIdentifier().getValue().startsWith("http://")) {
            try {
                uri_str = new URI(tax.getIdentifier().getValue()).toString();
            } catch (final URISyntaxException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                uri_str = null;
                e.printStackTrace();
            }
        } else if ((tax.getIdentifier() != null) && !ForesterUtil.isEmpty(tax.getIdentifier().getValue())
                && !ForesterUtil.isEmpty(tax.getIdentifier().getProvider())
                && (tax.getIdentifier().getProvider().equalsIgnoreCase("ncbi")
                || tax.getIdentifier().getProvider().equalsIgnoreCase("uniprot"))) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/"
                        + URLEncoder.encode(tax.getIdentifier().getValue(), ForesterConstants.UTF_8);
            } catch (final UnsupportedEncodingException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            }
        } else if (!ForesterUtil.isEmpty(tax.getScientificName())) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/?query="
                        + URLEncoder.encode(tax.getScientificName(), ForesterConstants.UTF_8);
            } catch (final UnsupportedEncodingException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            }
        } else if (!ForesterUtil.isEmpty(tax.getTaxonomyCode())) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/?query="
                        + URLEncoder.encode(tax.getTaxonomyCode(), ForesterConstants.UTF_8);
            } catch (final UnsupportedEncodingException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            }
        } else if (!ForesterUtil.isEmpty(tax.getCommonName())) {
            try {
                uri_str = "http://www.uniprot.org/taxonomy/?query="
                        + URLEncoder.encode(tax.getCommonName(), ForesterConstants.UTF_8);
            } catch (final UnsupportedEncodingException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            }
        }
        if (!ForesterUtil.isEmpty(uri_str)) {
            try {
                AptxUtil.launchWebBrowser(new URI(uri_str), "_aptx_tax");
            } catch (final IOException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            } catch (final URISyntaxException e) {
                AptxUtil.showErrorMessage(this, e.toString());
                e.printStackTrace();
            }
        } else {
            cannotOpenBrowserWarningMessage("taxonomic");
        }
    }

    final private void paintBranchLength(final Graphics2D g,
                                         final PhylogenyNode node,
                                         final boolean to_pdf,
                                         final boolean to_graphics_file) {
        g.setFont(getTreeFontSet().getSmallFont());
        final boolean bl_bw = (to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite();
        // fade a non-hit's branch-length number too (same dimNonMatch wash as the name/taxonomy labels)
        g.setColor(dimNonMatch(inkColor(to_pdf, to_graphics_file, getTreeColorSet().getBranchLengthColor()),
                isInFoundNodes(node), bl_bw));
        if (!node.isRoot()) {
            if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE) {
                TreePanel.drawString(FORMATTER_BRANCH_LENGTH.format(node.getDistanceToParent()),
                        node.getParent().getXcoord() + EURO_D,
                        node.getYcoord() - getTreeFontSet().getSmallMaxDescent(),
                        g);
            } else if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED) {
                TreePanel.drawString(FORMATTER_BRANCH_LENGTH.format(node.getDistanceToParent()),
                        node.getParent().getXcoord() + ROUNDED_D,
                        node.getYcoord() - getTreeFontSet().getSmallMaxDescent(),
                        g);
            } else {
                TreePanel.drawString(FORMATTER_BRANCH_LENGTH.format(node.getDistanceToParent()),
                        node.getParent().getXcoord() + 3,
                        node.getYcoord() - getTreeFontSet().getSmallMaxDescent(),
                        g);
            }
        } else {
            TreePanel.drawString(FORMATTER_BRANCH_LENGTH.format(node.getDistanceToParent()),
                    3,
                    node.getYcoord() - getTreeFontSet().getSmallMaxDescent(),
                    g);
        }
    }

    final private void paintBranchLite(final Graphics2D g,
                                       final float x1,
                                       final float x2,
                                       final float y1,
                                       final float y2,
                                       final PhylogenyNode node) {
        g.setColor(getTreeColorSet().getOvColor());
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR) {
            drawLine(x1, y1, x2, y2, g);
        } else {
            final float x2a = x2;
            final float x1a = x1;
            // draw the vertical line
            if (node.isFirstChildNode() || node.isLastChildNode()) {
                drawLine(x1, y1, x1, y2, g);
            }
            // draw the horizontal line
            drawLine(x1a, y2, x2a, y2, g);
        }
    }

    /**
     * Paint a branch which consists of a vertical and a horizontal bar
     */
    final private void paintBranchRectangular(final Graphics2D g,
                                              final float x1,
                                              final float x2,
                                              final float y1,
                                              final float y2,
                                              final PhylogenyNode node,
                                              final boolean to_pdf,
                                              final boolean to_graphics_file) {
        assignGraphicsForBranchWithColorForParentBranch(node, false, g, to_pdf, to_graphics_file);
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR) {
            drawLine(x1, y1, x2, y2, g);
        } else {
            final float x2a = x2;
            final float x1a = x1;
            // the drawn coords (x1/y1/x2/y2) are LOGICAL (the canvas is rotated by R), so cull against the viewport
            // mapped back into logical space -- correct in a vertical orientation too (was skipped there before)
            final Rectangle2D log_vis = (!to_graphics_file && !to_pdf) ? logicalVisibleRect() : null;
            float y2_r = 0;
            if (node.isFirstChildNode() || node.isLastChildNode()
                    || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE)
                    || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED)) {
                // screen-cull the vertical connector by its y-extent (both endpoints off the same side of the viewport)
                if ((log_vis != null)
                        && (((y2 < (log_vis.getMinY() - 20)) && (y1 < (log_vis.getMinY() - 20)))
                        || ((y2 > (log_vis.getMaxY() + 20)) && (y1 > (log_vis.getMaxY() + 20))))) {
                    // Do nothing.
                } else {
                    if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE) {
                        float x2c = x1 + EURO_D;
                        if (x2c > x2a) {
                            x2c = x2a;
                        }
                        drawLine(x1, y1, x2c, y2, g);
                    } else if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED) {
                        if (y2 > y1) {
                            y2_r = y2 - ROUNDED_D;
                            if (y2_r < y1) {
                                y2_r = y1;
                            }
                            drawLine(x1, y1, x1, y2_r, g);
                        } else {
                            y2_r = y2 + ROUNDED_D;
                            if (y2_r > y1) {
                                y2_r = y1;
                            }
                            drawLine(x1, y1, x1, y2_r, g);
                        }
                    } else {
                        drawLine(x1, y1, x1, y2, g);
                    }
                }
            }
            // draw the horizontal line (cull it when its row is off the top/bottom of the logical-space viewport)
            if ((log_vis != null)
                    && ((y2 < (log_vis.getMinY() - 20)) || (y2 > (log_vis.getMaxY() + 20)))) {
                return;
            }
            float x1_r = 0;
            if (!getControlPanel().isWidthBranches() || (PhylogenyMethods.getBranchWidthValue(node) == 1)) {
                if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED) {
                    x1_r = x1a + ROUNDED_D;
                    if (x1_r < x2a) {
                        drawLine(x1_r, y2, x2a, y2, g);
                    }
                } else if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE) {
                    final float x1c = x1a + EURO_D;
                    if (x1c < x2a) {
                        drawLine(x1c, y2, x2a, y2, g);
                    }
                } else {
                    drawLine(x1a, y2, x2a, y2, g);
                }
            } else {
                final double w = PhylogenyMethods.getBranchWidthValue(node);
                if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED) {
                    x1_r = x1a + ROUNDED_D;
                    if (x1_r < x2a) {
                        drawRectFilled(x1_r, y2 - (w / 2), x2a - x1_r, w, g);
                    }
                } else if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE) {
                    final float x1c = x1a + EURO_D;
                    if (x1c < x2a) {
                        drawRectFilled(x1c, y2 - (w / 2), x2a - x1c, w, g);
                    }
                } else {
                    drawRectFilled(x1a, y2 - (w / 2), x2a - x1a, w, g);
                }
            }
            if ((getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED)) {
                if (x1_r > x2a) {
                    x1_r = x2a;
                }
                if (y2 > y2_r) {
                    final double diff = y2 - y2_r;
                    _arc.setArc(x1, y2_r - diff, 2 * (x1_r - x1), 2 * diff, 180, 90, Arc2D.OPEN);
                } else {
                    _arc.setArc(x1, y2, 2 * (x1_r - x1), 2 * (y2_r - y2), 90, 90, Arc2D.OPEN);
                }
                g.draw(_arc);
            }
        }
        if (node.isExternal()) {
            paintNodeBox(x2, y2, node, g, to_pdf, to_graphics_file);
        }
    }

    final private double paintCirculars(final PhylogenyNode n,
                                        final Phylogeny phy,
                                        final float center_x,
                                        final float center_y,
                                        final double radius,
                                        final boolean radial_labels,
                                        final Graphics2D g,
                                        final boolean to_pdf,
                                        final boolean to_graphics_file) {
        if (n.isExternal() || n.isCollapse()) { //~~circ collapse
            return _urt_nodeid_angle_map.get(n.getId());
        } else {
            final List<PhylogenyNode> descs = n.getDescendants();
            double sum = 0;
            for (final PhylogenyNode desc : descs) {
                sum += paintCirculars(desc,
                        phy,
                        center_x,
                        center_y,
                        radius,
                        radial_labels,
                        g,
                        to_pdf,
                        to_graphics_file);
            }
            double r = 0;
            if (!n.isRoot()) {
                r = 1 - (((double) _circ_max_depth - n.calculateDepth()) / _circ_max_depth);
            }
            final double theta = sum / descs.size();
            n.setXcoord((float) (center_x + (r * radius * Math.cos(theta))));
            n.setYcoord((float) (center_y + (r * radius * Math.sin(theta))));
            _urt_nodeid_angle_map.put(n.getId(), theta);
            for (final PhylogenyNode desc : descs) {
                paintBranchCircular(n, desc, g, radial_labels, to_pdf, to_graphics_file);
            }
            return theta;
        }
    }

    final private void paintCircularsLite(final PhylogenyNode n,
                                          final Phylogeny phy,
                                          final int center_x,
                                          final int center_y,
                                          final int radius,
                                          final Graphics2D g) {
        if (n.isExternal()) {
            return;
        } else {
            final List<PhylogenyNode> descs = n.getDescendants();
            for (final PhylogenyNode desc : descs) {
                paintCircularsLite(desc, phy, center_x, center_y, radius, g);
            }
            float r = 0;
            if (!n.isRoot()) {
                r = 1 - (((float) _circ_max_depth - n.calculateDepth()) / _circ_max_depth);
            }
            final double theta = _urt_nodeid_angle_map.get(n.getId());
            n.setXSecondary((float) (center_x + (radius * r * Math.cos(theta))));
            n.setYSecondary((float) (center_y + (radius * r * Math.sin(theta))));
            for (final PhylogenyNode desc : descs) {
                paintBranchCircularLite(n, desc, g);
            }
        }
    }

    final private void paintCollapsedNode(final Graphics2D g,
                                          final PhylogenyNode node,
                                          final boolean to_graphics_file,
                                          final boolean to_pdf,
                                          final boolean is_in_found_nodes) {
        Color c = null;
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite();
        // A collapsed clade's tips are hidden, so signal any selection/search hits ON the triangle itself:
        // fill it in the found colour when EVERY tip is a hit (reads as "this whole group is selected"), or just
        // outline it (below) when only some are. Selection is clade-wide (collapse is a view state), so a
        // clade-select shows here even though the individual tips are not drawn. Skipped in black-and-white.
        final int[] fc = collapsedCladeFoundCounts(node); // {hits in set 0, hits in set 1, tips hit, total tips}
        Color found_c = null;
        if (!bw && (fc[2] > 0)) {
            found_c = ((fc[0] > 0) && (fc[1] > 0)) ? getTreeColorSet().getFoundColor0and1()
                    : (fc[0] > 0) ? getTreeColorSet().getFoundColor0() : getTreeColorSet().getFoundColor1();
        }
        final boolean all_tips_found = (found_c != null) && (fc[2] == fc[3]);
        if (bw) {
            c = Color.BLACK;
        }
        else if (all_tips_found) {
            c = found_c; // the whole collapsed clade is selected -> paint the triangle in the found colour
        }
        else if (getOptions().isColorLabelsSameAsParentBranch() && getControlPanel().isUseVisualStyles()
                && (PhylogenyMethods.getBranchColorValue(node) != null)) {
            c = PhylogenyMethods.getBranchColorValue(node);
        } else if (to_pdf) {
            g.setColor(getTreeColorSet().getBranchColorForPdf());
        } else {
            c = getTreeColorSet().getCollapseFillColor();
        }
        double d = fc[3];
        float xxx;
        double s = 0;
        if (getControlPanel().isDrawPhylogram()) {
            if (d > 1000) {
                d = 0.75 * _y_distance;
            } else {
                d = 0.25 * Math.log10(d) * _y_distance;
            }
            final float half_box_size = 0.5f * getOptions().getDefaultNodeShapeSize();
            if (d < half_box_size) {
                d = half_box_size;
            }
            _polygon.reset();
            final float xx = node.getXcoord() - (getOptions().getDefaultNodeShapeSize());
            xxx = xx > (node.getParent().getXcoord() + 1) ? xx : node.getParent().getXcoord() + 1;
            _polygon.moveTo(xxx, node.getYcoord() + 0.5);
            _polygon.lineTo(xxx, node.getYcoord() - 0.5);
            s = _options.isCollapsedWithAverageHeigh()
                    ? PhylogenyMethods.calculateAverageTreeHeight(node) * _x_correction_factor
                    : 1;
            _polygon.lineTo(node.getXcoord() + s, node.getYcoord() - d);
            _polygon.lineTo(node.getXcoord() + s, node.getYcoord() + d);
            _polygon.closePath();
        } else {
            if (d > 1000) {
                d = _y_distance;
            } else {
                d = (Math.log10(d) * _y_distance) / 2.5;
            }
            final int box_size = getOptions().getDefaultNodeShapeSize() + 1;
            if (d < box_size) {
                d = box_size;
            }
            final float xx = node.getXcoord() - (2 * box_size);
            xxx = xx > (node.getParent().getXcoord() + 1) ? xx : node.getParent().getXcoord() + 1;
            _polygon.reset();
            _polygon.moveTo(xxx, node.getYcoord());
            _polygon.lineTo(node.getXcoord() + 1, node.getYcoord() - d);
            _polygon.lineTo(node.getXcoord() + 1, node.getYcoord() + d);
            _polygon.closePath();
        }
        if (getOptions().getDefaultNodeFill() == NodeVisualData.NodeFill.SOLID) {
            g.setColor(c);
            g.fill(_polygon);
        } else if (getOptions().getDefaultNodeFill() == NodeVisualData.NodeFill.NONE) {
            g.setColor(getBackground());
            g.fill(_polygon);
            g.setColor(c);
            g.draw(_polygon);
        } else if (getOptions().getDefaultNodeFill() == NodeFill.GRADIENT) {
            g.setPaint(new GradientPaint(xxx,
                    node.getYcoord(),
                    getBackground(),
                    node.getXcoord(),
                    (float) (node.getYcoord() - d),
                    c,
                    false));
            g.fill(_polygon);
            g.setPaint(c);
            g.draw(_polygon);
        }
        if ((found_c != null) && !all_tips_found) {
            // only SOME of the clade's tips are selected/found -> outline the triangle in the found colour as a
            // partial hint (a full fill would over-state it, an omitted mark would hide it -- see also the
            // "[found/total]" count in the collapsed label)
            final Color saved = g.getColor();
            g.setColor(found_c);
            g.draw(_polygon);
            g.setColor(saved);
        }
        if (isVerticalOrientation()) {
            // a collapsed clade reads as a tip: tilt its label 45deg, and (aligned mode) draw its leader as vertical
            // geometry + pivot the label at the aligned column -- like the external-tip labels
            final double sf = s; // capture for the lambda (s is reassigned above)
            if (isAlignedTipLabel(node)) {
                drawConnection(node.getXcoord(), labelTextStartX(node), node.getYcoord(), 5, 20, g);
            }
            withNodeTextFrame(g, labelTextStartX(node), node.getYcoord(), tipLabelAngle(),
                    () -> paintNodeData(g, node, to_graphics_file, to_pdf, is_in_found_nodes, sf));
        } else {
            paintNodeData(g, node, to_graphics_file, to_pdf, is_in_found_nodes, s);
        }
    }

    /**
     * The branch label for a node's confidences: the MAD support (type {@code "MAD"}) first -- only
     * when {@code show_mad} is on, and never hidden by {@code min_confidence} (for MAD, low is good)
     * -- followed by the regular confidences (e.g. bootstrap), each shown only if {@code >=
     * min_confidence} (for those, high is good). Joined by '/', so a MAD-rooted tree reads as e.g.
     * "0/90" (MAD/bootstrap). Non-finite values are skipped. Returns "" when nothing should show.
     */
    static String confidenceLabel(final List<Confidence> confidences, final boolean show_mad,
                                  final double min_confidence, final boolean show_stddev, final int digits) {
        final StringBuilder sb = new StringBuilder();
        boolean not_first = false;
        if (show_mad) {
            for (final Confidence c : confidences) {
                if (PhylogenyMethods.MAD_CONFIDENCE_TYPE.equals(c.getType())) {
                    final double value = c.getValue();
                    if ((value != Confidence.CONFIDENCE_DEFAULT_VALUE) && Double.isFinite(value)) {
                        sb.append(FORMATTER_MAD.format(value)); // own fixed precision, not the "digits" setting
                        not_first = true;
                    }
                }
            }
        }
        for (final Confidence c : confidences) {
            if (PhylogenyMethods.MAD_CONFIDENCE_TYPE.equals(c.getType())) {
                continue; // MAD handled above
            }
            final double value = c.getValue();
            // skip non-finite values, and regular values below the "min. confidence" threshold
            if ((value != Confidence.CONFIDENCE_DEFAULT_VALUE) && Double.isFinite(value) && (value >= min_confidence)) {
                if (not_first) {
                    sb.append("/");
                }
                else {
                    not_first = true;
                }
                sb.append(FORMATTER_CONFIDENCE.format(ForesterUtil.round(value, digits)));
                if (show_stddev && (c.getStandardDeviation() != Confidence.CONFIDENCE_DEFAULT_VALUE)
                        && Double.isFinite(c.getStandardDeviation())) {
                    sb.append("(");
                    sb.append(FORMATTER_CONFIDENCE.format(ForesterUtil.round(c.getStandardDeviation(), digits)));
                    sb.append(")");
                }
            }
        }
        return sb.toString();
    }

    final private void paintConfidenceValues(final Graphics2D g,
                                             final PhylogenyNode node,
                                             final boolean to_pdf,
                                             final boolean to_graphics_file) {
        final List<Confidence> confidences = node.getBranchData().getConfidences();
        Collections.sort(confidences);
        // "Min. confidence shown" is a fraction of the tree's detected support scale (computed once per paint in
        // _confidence_scale_max), so the same setting means the same thing on 0-1 posterior and 0-100 bootstrap.
        // Limitation: the scale is a single per-tree ceiling, so a node carrying confidences on TWO different
        // scales at once (e.g. a 0-1 posterior AND a 0-100 bootstrap) filters both against the one ceiling.
        final double min_confidence = getOptions().getMinConfidenceFraction() * _confidence_scale_max;
        final StringBuilder sb = new StringBuilder(confidenceLabel(confidences, getOptions().isShowMadConfidence(),
                min_confidence, getOptions().isShowConfidenceStddev(),
                getOptions().getNumberOfDigitsAfterCommaForConfidenceValues()));
        if (sb.length() > 0) {
            final float parent_x = node.getParent().getXcoord();
            float x = node.getXcoord();
            g.setFont(getTreeFontSet().getSmallFont());
            if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE) {
                x += EURO_D;
            } else if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED) {
                x += ROUNDED_D;
            }
            final boolean conf_bw = (to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite();
            // fade a non-hit's confidence value too (same dimNonMatch wash as the name/taxonomy labels)
            g.setColor(dimNonMatch(inkColor(to_pdf, to_graphics_file, getTreeColorSet().getConfidenceColor()),
                    isInFoundNodes(node), conf_bw));
            final String conf_str = sb.toString();
            TreePanel.drawString(conf_str,
                    parent_x + ((x - parent_x
                            - getTreeFontSet().getFontMetricsSmall().stringWidth(conf_str)) / 2),
                    (node.getYcoord() + getTreeFontSet().getSmallMaxAscent()) - 1,
                    g);
        }
    }

    // gap (px) between a vertical-orientation branch and its support / branch-length labels drawn beside it
    private static final float VERTICAL_BRANCH_LABEL_PAD = 3f;

    /** Support and/or branch-length for a VERTICAL orientation: centered at the incoming branch's midpoint, to the
     *  RIGHT of the (on-screen vertical) branch, drawn HORIZONTALLY as "support length" -- the support in the
     *  confidence colour, a space, then the length in the branch-length colour (either alone when only one is shown).
     *  Placed explicitly in the base frame; the horizontal-layout positions merely rotated would sit ON the branch. */
    private void paintBranchDataRightVertical(final Graphics2D g, final PhylogenyNode node, final boolean to_pdf,
                                              final boolean to_graphics_file) {
        if (node.isRoot()) {
            return;
        }
        String support = "";
        if (isShowConfidenceValuesForNode(node)) {
            final List<Confidence> confidences = node.getBranchData().getConfidences();
            Collections.sort(confidences);
            final double min_confidence = getOptions().getMinConfidenceFraction() * _confidence_scale_max;
            support = confidenceLabel(confidences, getOptions().isShowMadConfidence(), min_confidence,
                    getOptions().isShowConfidenceStddev(),
                    getOptions().getNumberOfDigitsAfterCommaForConfidenceValues());
        }
        final String length = shouldWriteBranchLength(node)
                ? FORMATTER_BRANCH_LENGTH.format(node.getDistanceToParent()) : "";
        if (support.isEmpty() && length.isEmpty()) {
            return;
        }
        final Point2D.Double mid = new Point2D.Double((node.getParent().getXcoord() + node.getXcoord()) / 2.0,
                node.getYcoord());
        _orientation_R.transform(mid, mid); // midpoint of the (vertically drawn) incoming branch
        final AffineTransform saved = g.getTransform();
        g.setTransform(_orientation_base_transform);
        g.setFont(getTreeFontSet().getSmallFont());
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite();
        final boolean found = isInFoundNodes(node);
        final float baseline = (float) (mid.y + (getTreeFontSet().getSmallMaxAscent() / 2.0) - 1);
        float x = (float) (mid.x + VERTICAL_BRANCH_LABEL_PAD);
        if (!support.isEmpty()) {
            g.setColor(dimNonMatch(inkColor(to_pdf, to_graphics_file, getTreeColorSet().getConfidenceColor()), found, bw));
            TreePanel.drawString(support, x, baseline, g);
            x += getTreeFontSet().getFontMetricsSmall().stringWidth(support + " ");
        }
        if (!length.isEmpty()) {
            g.setColor(dimNonMatch(inkColor(to_pdf, to_graphics_file, getTreeColorSet().getBranchLengthColor()), found,
                    bw));
            TreePanel.drawString(length, x, baseline, g);
        }
        g.setTransform(saved);
    }

    /** Internal-node label for a VERTICAL orientation: horizontal, RIGHT-ALIGNED so it ends just LEFT of the branch,
     *  centered at the incoming branch's midpoint -- mirroring the support/length on the right. Draws the taxonomy
     *  (its own italic/colour) then the node-data string. Tips are skipped (their 45deg/90deg label is drawn
     *  elsewhere); a long label extends leftward toward the neighbouring subtree, as accepted. */
    private void paintInternalLabelLeftVertical(final Graphics2D g, final PhylogenyNode node, final boolean to_pdf,
                                                final boolean to_graphics_file) {
        if (node.isRoot() || !getControlPanel().isShowInternalData()) {
            return;
        }
        final boolean using_visual_font = setFont(g, node);
        final boolean is_in_found_nodes = isInFoundNodes(node);
        final Taxonomy taxonomy = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy() : null;
        final boolean show_tax = (taxonomy != null)
                && (getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyScientificNames()
                        || getControlPanel().isShowTaxonomyCommonNames() || getControlPanel().isShowTaxonomyRank());
        final StringBuilder sb = new StringBuilder();
        nodeDataAsSB(node, sb);
        final String data_str = sb.toString().trim();
        final int tax_w = show_tax ? taxonomyLabel(g, taxonomy, 0, 0, to_pdf, false) : 0;
        final int data_w = data_str.isEmpty() ? 0
                : labelStringWidth(g, data_str, using_visual_font, is_in_found_nodes, to_pdf);
        final int total = tax_w + data_w;
        if (total == 0) {
            return;
        }
        // the taxonomy label always ends with a trailing part-separator space; when it is the RIGHTMOST drawn element
        // (no node-data follows) that space would push the visible label one space-width left of the branch, so
        // right-align on the visible width instead (the invisible space then falls harmlessly into the branch gap)
        final int align_total = TreePanelUtil.internalLabelAlignWidth(tax_w, data_w,
                getFontMetrics(g.getFont()).charWidth(' '), show_tax, data_str.isEmpty());
        final Point2D.Double mid = new Point2D.Double((node.getParent().getXcoord() + node.getXcoord()) / 2.0,
                node.getYcoord());
        _orientation_R.transform(mid, mid);
        final AffineTransform saved = g.getTransform();
        g.setTransform(_orientation_base_transform);
        // right-aligned so it ends just left of the branch; clamp at the canvas left edge so a long label near the
        // breadth edge crowds the branch (accepted) rather than clipping off-canvas at x<0 (unrecoverable)
        final float start_x = Math.max(0f, (float) (mid.x - VERTICAL_BRANCH_LABEL_PAD - align_total));
        final float baseline = (float) (mid.y + (getFontMetrics(g.getFont()).getAscent() / 2.0) - 1);
        if (show_tax) {
            setColor(g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getTaxonomyColor());
            taxonomyLabel(g, taxonomy, start_x, baseline, to_pdf, true);
        }
        if (!data_str.isEmpty()) {
            setColor(g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getSequenceColor());
            TreePanel.drawString(data_str, start_x + tax_w, baseline, g);
        }
        g.setTransform(saved);
    }

    /**
     * Whether the support/confidence value on {@code node}'s incoming branch should be drawn: the toggle is on,
     * the node is an internal, non-root node that carries a confidence, and the layout is one of the horizontal
     * styles that place labels on branches. Shared by the regular-node and collapsed-clade paint paths so a
     * collapsed subtree still shows the support of the branch leading to it (a collapsed clade is still an
     * internal node).
     */
    final private boolean isShowConfidenceValuesForNode(final PhylogenyNode node) {
        return getControlPanel().isShowConfidenceValues() && !node.isExternal() && !node.isRoot()
                && ((getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED)
                        || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR)
                        || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE))
                && node.getBranchData().isHasConfidences();
    }

    final private void paintGainedAndLostCharacters(final Graphics2D g,
                                                    final PhylogenyNode node,
                                                    final String gained,
                                                    final String lost) {
        if (node.getParent() != null) {
            final float parent_x = node.getParent().getXcoord();
            final float x = node.getXcoord();
            g.setFont(getTreeFontSet().getLargeFont());
            g.setColor(getTreeColorSet().getGainedCharactersColor());
            TreePanel.drawString(gained,
                    parent_x + ((x - parent_x
                            - getFontMetricsForLargeDefaultFont().stringWidth(gained)) / 2),
                    (node.getYcoord() - getFontMetricsForLargeDefaultFont().getMaxDescent()),
                    g);
            g.setColor(getTreeColorSet().getLostCharactersColor());
            TreePanel
                    .drawString(lost,
                            parent_x + ((x - parent_x - getFontMetricsForLargeDefaultFont().stringWidth(lost))
                                    / 2),
                            (node.getYcoord() + getFontMetricsForLargeDefaultFont().getMaxAscent()),
                            g);
        }
    }

    private void paintMolecularSequences(final Graphics2D g, final PhylogenyNode node, final boolean to_pdf) {
        final RenderableMsaSequence rs = RenderableMsaSequence.createInstance(node.getNodeData().getSequence()
                .getMolecularSequence(), node.getNodeData().getSequence().getType(), getConfiguration());
        if (rs != null) {
            final int default_height = 8;
            final float y = getYdistance();
            final int h = (y / 2) < default_height ? ForesterUtil.roundToInt(y * 2) : default_height;
            rs.setRenderingHeight(h > 1 ? h : 1);
            if (getControlPanel().isDrawPhylogram()) {
                rs.render((float) ((getMaxDistanceToRoot() * getXcorrectionFactor()) + _length_of_longest_text),
                        node.getYcoord() - (h / 2.0f),
                        g,
                        this,
                        to_pdf);
            } else {
                rs.render(getPhylogeny().getFirstExternalNode().getXcoord()
                        + _length_of_longest_text, node.getYcoord() - (h / 2.0f), g, this, to_pdf);
            }
        }
    }

    /**
     * Draws the branch-support (confidence) symbol on the <i>middle of the branch</i> (support is a
     * branch property, not a node one -- matching where the confidence <i>value</i> text is drawn),
     * per the {@link Options.SUPPORT_VISUALIZATION} mode. Monochrome (branch color) so it does not
     * compete with clade/taxonomy coloring; uses the absolute support scale (1 or 100, detected per
     * paint into {@code _confidence_scale_max}) rather than normalizing to the maximum present, so a
     * given symbol means the same thing across trees. THRESHOLD_MARKS draws a fixed dot only at or
     * above the cutoff; SIZE_SCALED draws a dot whose diameter grows with support.
     */
    final private void paintNodeSupportSymbol(final float x,
                                              final float y,
                                              final PhylogenyNode node,
                                              final Graphics2D g,
                                              final boolean to_pdf,
                                              final boolean to_graphics_file) {
        final Options.SUPPORT_VISUALIZATION mode = getOptions().getSupportVisualization();
        // COLOR_BRANCHES shows support by coloring the branch (assignGraphicsForBranchWithColorForParentBranch),
        // not by a symbol, so it draws no dot here.
        if ((mode == Options.SUPPORT_VISUALIZATION.NONE)
                || (mode == Options.SUPPORT_VISUALIZATION.COLOR_BRANCHES) || node.isExternal()
                || !node.getBranchData().isHasConfidences()) {
            return;
        }
        final double conf = PhylogenyMethods.getConfidenceValue(node);
        if (conf < 0.0) {
            return;
        }
        final float diameter;
        if (mode == Options.SUPPORT_VISUALIZATION.THRESHOLD_MARKS) {
            if (!TreePanelUtil.isSupportAtOrAboveThreshold(conf, _confidence_scale_max,
                    getOptions().getSupportThreshold())) {
                return;
            }
            diameter = AptxConstants.SUPPORT_SYMBOL_MAX_DIAMETER;
        } else {
            diameter = TreePanelUtil.supportSymbolSize(conf,
                    _confidence_scale_max,
                    AptxConstants.SUPPORT_SYMBOL_MIN_DIAMETER,
                    AptxConstants.SUPPORT_SYMBOL_MAX_DIAMETER);
            if (diameter <= 0.0f) {
                return;
            }
        }
        Color c = getGraphicsForNodeBoxWithColorForParentBranch(node);
        if (to_pdf && (c == getTreeColorSet().getBranchColor())) {
            c = getTreeColorSet().getBranchColorForPdf();
        }
        if ((to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite()) {
            c = Color.BLACK;
        }
        g.setColor(c);
        // Center the symbol on the middle of the branch (parent -> node), not at the node on its right
        // end. The root has no branch, so it falls back to the passed-in node position.
        float cx = x;
        float cy = y;
        final PhylogenyNode parent = node.getParent();
        if (parent != null) {
            final boolean radial = (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                    || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR);
            final float[] center = TreePanelUtil.supportSymbolCenter(parent.getXcoord(), node.getXcoord(),
                    parent.getYcoord(), node.getYcoord(), radial);
            cx = center[0];
            cy = center[1];
        }
        final float half = diameter / 2.0f;
        drawOvalFilled(cx - half, cy - half, diameter, diameter, g);
    }

    /**
     * Draw a box at the indicated node.
     *
     * @param x
     * @param y
     * @param node
     * @param g
     */
    final private void paintNodeBox(final float x,
                                    final float y,
                                    final PhylogenyNode node,
                                    final Graphics2D g,
                                    final boolean to_pdf,
                                    final boolean to_graphics_file) {
        if (node.isCollapse()) {
            return;
        }
        paintNodeSupportSymbol(x, y, node, g, to_pdf, to_graphics_file);
        if ((isInFoundNodes(node))
                || (getOptions().isShowDefaultNodeShapesExternal() && node.isExternal())
                || (getOptions().isShowDefaultNodeShapesInternal() && node.isInternal())
                || (getOptions().isShowDefaultNodeShapesForMarkedNodes()
                && (node.getNodeData().getNodeVisualData() != null)
                && (!node.getNodeData().getNodeVisualData().isEmpty()))
                || (getControlPanel().isUseVisualStyles() && ((node.getNodeData().getNodeVisualData() != null)
                && ((node.getNodeData().getNodeVisualData().getNodeColor() != null)
                || (node.getNodeData().getNodeVisualData().getSize() != NodeVisualData.DEFAULT_SIZE)
                || (node.getNodeData().getNodeVisualData().getFillType() != NodeFill.DEFAULT)
                || (node.getNodeData().getNodeVisualData().getShape() != NodeShape.DEFAULT))))
                || (getControlPanel().isEvents() && node.isHasAssignedEvent()
                && (node.getNodeData().getEvent().isDuplication()
                || node.getNodeData().getEvent().isSpeciation()
                || node.getNodeData().getEvent().isSpeciationOrDuplication()))) {
            NodeVisualData vis = null;
            if (getControlPanel().isUseVisualStyles() && (node.getNodeData().getNodeVisualData() != null)
                    && (!node.getNodeData().getNodeVisualData().isEmpty())) {
                vis = node.getNodeData().getNodeVisualData();
            }
            float box_size = getOptions().getDefaultNodeShapeSize();
            if ((vis != null) && (vis.getSize() != NodeVisualData.DEFAULT_SIZE)) {
                box_size = vis.getSize();
            }
            final float half_box_size = box_size / 2.0f;
            Color outline_color = null;
            if ((to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite()) {
                outline_color = Color.BLACK;
            } else if (isInFoundNodes(node)) {
                outline_color = getColorForFoundNode(node);
            } else if (vis != null) {
                if (vis.getNodeColor() != null) {
                    outline_color = vis.getNodeColor();
                } else if (vis.getFontColor() != null) {
                    outline_color = vis.getFontColor();
                }
            } else if (getControlPanel().isEvents() && TreePanelUtil.isHasAssignedEvent(node)) {
                final Event event = node.getNodeData().getEvent();
                if (event.isDuplication()) {
                    outline_color = getTreeColorSet().getDuplicationBoxColor();
                } else if (event.isSpeciation()) {
                    outline_color = getTreeColorSet().getSpecBoxColor();
                } else if (event.isSpeciationOrDuplication()) {
                    outline_color = getTreeColorSet().getDuplicationOrSpeciationColor();
                }
            }
            if (outline_color == null) {
                outline_color = getGraphicsForNodeBoxWithColorForParentBranch(node);
                if (to_pdf && (outline_color == getTreeColorSet().getBranchColor())) {
                    outline_color = getTreeColorSet().getBranchColorForPdf();
                }
            }
            NodeShape shape = null;
            if (vis != null) {
                if (vis.getShape() == NodeShape.CIRCLE) {
                    shape = NodeShape.CIRCLE;
                } else if (vis.getShape() == NodeShape.RECTANGLE) {
                    shape = NodeShape.RECTANGLE;
                }
            }
            if (shape == null) {
                if (getOptions().getDefaultNodeShape() == NodeShape.CIRCLE) {
                    shape = NodeShape.CIRCLE;
                } else if (getOptions().getDefaultNodeShape() == NodeShape.RECTANGLE) {
                    shape = NodeShape.RECTANGLE;
                }
            }
            NodeFill fill = null;
            if (vis != null) {
                if (vis.getFillType() == NodeFill.SOLID) {
                    fill = NodeFill.SOLID;
                } else if (vis.getFillType() == NodeFill.NONE) {
                    fill = NodeFill.NONE;
                } else if (vis.getFillType() == NodeFill.GRADIENT) {
                    fill = NodeFill.GRADIENT;
                }
            }
            if (fill == null) {
                if (getOptions().getDefaultNodeFill() == NodeFill.SOLID) {
                    fill = NodeFill.SOLID;
                } else if (getOptions().getDefaultNodeFill() == NodeFill.NONE) {
                    fill = NodeFill.NONE;
                } else if (getOptions().getDefaultNodeFill() == NodeFill.GRADIENT) {
                    fill = NodeFill.GRADIENT;
                }
            }
            Color vis_fill_color = null;
            if ((vis != null) && (vis.getNodeColor() != null)) {
                vis_fill_color = vis.getNodeColor();
            }
            if (shape == NodeShape.CIRCLE) {
                if (fill == NodeFill.GRADIENT) {
                    drawOvalGradient(x - half_box_size,
                            y - half_box_size,
                            box_size,
                            box_size,
                            g,
                            to_pdf ? Color.WHITE : outline_color,
                            to_pdf ? outline_color : getBackground(),
                            outline_color);
                } else if (fill == NodeFill.NONE) {
                    Color background = getBackground();
                    if (to_pdf) {
                        background = Color.WHITE;
                    }
                    drawOvalGradient(x - half_box_size,
                            y - half_box_size,
                            box_size,
                            box_size,
                            g,
                            background,
                            background,
                            outline_color);
                } else if (fill == NodeVisualData.NodeFill.SOLID) {
                    if (vis_fill_color != null) {
                        g.setColor(vis_fill_color);
                    } else {
                        g.setColor(outline_color);
                    }
                    drawOvalFilled(x - half_box_size, y - half_box_size, box_size, box_size, g);
                }
            } else if (shape == NodeVisualData.NodeShape.RECTANGLE) {
                if (fill == NodeVisualData.NodeFill.GRADIENT) {
                    drawRectGradient(x - half_box_size,
                            y - half_box_size,
                            box_size,
                            box_size,
                            g,
                            to_pdf ? Color.WHITE : outline_color,
                            to_pdf ? outline_color : getBackground(),
                            outline_color);
                } else if (fill == NodeVisualData.NodeFill.NONE) {
                    Color background = getBackground();
                    if (to_pdf) {
                        background = Color.WHITE;
                    }
                    drawRectGradient(x - half_box_size,
                            y - half_box_size,
                            box_size,
                            box_size,
                            g,
                            background,
                            background,
                            outline_color);
                } else if (fill == NodeVisualData.NodeFill.SOLID) {
                    if (vis_fill_color != null) {
                        g.setColor(vis_fill_color);
                    } else {
                        g.setColor(outline_color);
                    }
                    drawRectFilled(x - half_box_size, y - half_box_size, box_size, box_size, g);
                }
            }
        }
    }

    final private int paintNodeData(final Graphics2D g,
                                    final PhylogenyNode node,
                                    final boolean to_graphics_file,
                                    final boolean to_pdf,
                                    final boolean is_in_found_nodes,
                                    final double add) {
        if (isNodeDataInvisible(node) && !to_graphics_file && !to_pdf) {
            return 0;
        }
        // horizontal mode draws the branch-length value here; a vertical orientation draws it (with the support value)
        // horizontally to the RIGHT of the branch instead, in paintNodeRectangular via paintBranchDataRightVertical
        if (!isVerticalOrientation() && shouldWriteBranchLength(node)) {
            paintBranchLength(g, node, to_pdf, to_graphics_file);
        }
        if (!getControlPanel().isShowInternalData() && !node.isExternal() && !node.isCollapse()) {
            return 0;
        }
        if (!getControlPanel().isShowExternalData() && (node.isExternal() || node.isCollapse())) {
            return 0;
        }
        _sb.setLength(0);
        int x = 0;
        if (add > 0) {
            x += add;
        }
        final int half_box_size = getOptions().getDefaultNodeShapeSize() / 2;
        if (((isColorByProperty() && (node.isExternal() || node.isCollapse())))
                || (isSizeByProperty() && node.isExternal())) {
            drawPropertyColorDot(g, node);
        }
        if (usesAboveBranchInternalLabel(node)) {
            return paintInternalLabelAboveBranch(g, node, is_in_found_nodes, to_pdf, to_graphics_file, half_box_size);
        }

        // Whether a taxonomy label was actually drawn. taxonomyLabel no longer leaves the label in _sb
        // (the old paintTaxonomy did, and this flag used to read _sb.length()), so derive it from the
        // painted advance instead -- the collapsed-node label logic below depends on it.
        boolean saw_species = false;
        if ((getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyScientificNames()
                || getControlPanel().isShowTaxonomyCommonNames() || getControlPanel().isShowTaxonomyRank())
                && node.getNodeData().isHasTaxonomy()) {
            final int taxonomy_width = paintTaxonomy(g, node, is_in_found_nodes, to_pdf, to_graphics_file, x);
            x += taxonomy_width;
            saw_species = taxonomy_width > 0;
        }
        setColor(g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getSequenceColor());
        _sb.setLength(0);
        nodeDataAsSB(node, _sb);
        if (node.isCollapse() && ((!node.isRoot() && !node.getParent().isCollapse()) || node.isRoot())) {
            if ((_sb.length() == 0) && !saw_species) {
                if (getOptions().isShowAbbreviatedLabelsForCollapsedNodes()
                        && (getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyScientificNames()
                        || getControlPanel().isShowSeqNames() || getControlPanel().isShowNodeNames())) {
                    final PhylogenyNode first = PhylogenyMethods.getFirstExternalNode(node);
                    final PhylogenyNode last = PhylogenyMethods.getLastExternalNode(node);
                    if (getControlPanel().isShowTaxonomyCode() && first.getNodeData().isHasTaxonomy()
                            && last.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty(first.getNodeData().getTaxonomy().getTaxonomyCode())
                            && !ForesterUtil.isEmpty(last.getNodeData().getTaxonomy().getTaxonomyCode())) {
                        addLabelForCollapsed(first.getNodeData().getTaxonomy().getTaxonomyCode(),
                                last.getNodeData().getTaxonomy().getTaxonomyCode(),
                                node.getAllExternalDescendants().size(),
                                node);
                    } else if (getControlPanel().isShowTaxonomyScientificNames() && first.getNodeData().isHasTaxonomy()
                            && last.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty(first.getNodeData().getTaxonomy().getScientificName())
                            && !ForesterUtil.isEmpty(last.getNodeData().getTaxonomy().getScientificName())) {
                        addLabelForCollapsed(first.getNodeData().getTaxonomy().getScientificName(),
                                last.getNodeData().getTaxonomy().getScientificName(),
                                node.getAllExternalDescendants().size(),
                                node);
                    } else if (getControlPanel().isShowSeqNames() && first.getNodeData().isHasSequence()
                            && last.getNodeData().isHasSequence()
                            && !ForesterUtil.isEmpty(first.getNodeData().getSequence().getName())
                            && !ForesterUtil.isEmpty(last.getNodeData().getSequence().getName())) {
                        addLabelForCollapsed(first.getNodeData().getSequence().getName(),
                                last.getNodeData().getSequence().getName(),
                                node.getAllExternalDescendants().size(),
                                node);
                    } else if (getControlPanel().isShowNodeNames() && !ForesterUtil.isEmpty(first.getName())
                            && !ForesterUtil.isEmpty(last.getName())) {
                        addLabelForCollapsed(first.getName(),
                                last.getName(),
                                node.getAllExternalDescendants().size(),
                                node);
                    }
                }
            } else if ((_sb.length() > 0) || saw_species) {
                _sb.append(" [");
                _sb.append(node.getAllExternalDescendants().size());
                _sb.append("]");
                if ((_found_nodes_0 != null) || (_found_nodes_1 != null)) {
                    final int[] res = calcFoundNodesInSubtree(node);
                    if (res[0] > 0) {
                        _sb.append(" [");
                        _sb.append(res[0]);
                        _sb.append("/");
                        _sb.append(res[1]);
                        _sb.append("]");
                    }
                }
            }
        } else {
            // _sb.setLength( 0 );
        }
        // nodeDataAsSB( node, _sb );
        final boolean using_visual_font = setFont(g, node);
        float down_shift_factor = 3.0f;
        if (!node.isExternal() && (node.getNumberOfDescendants() == 1)) {
            down_shift_factor = 1;
        }
        float pos_x;
        if (tipLabelsBelowColumns() && (node.isExternal() || node.isCollapse())) {
            // clustergram "labels below columns": the label is drawn PAST the tip-aligned columns, at the SAME base
            // (annotationColumnsEndX) that labelTextStartX feeds the rotated-frame pivot -- one shared method so the
            // drawn position and the tilt pivot can never disagree (a mismatch mis-anchors the 45/90 deg label).
            pos_x = labelSegmentStartX(annotationColumnsEndX(), half_box_size, x);
        } else if ((getControlPanel().getTreeDisplayType() == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM)
                && (node.isExternal() || node.isCollapse())) {
            pos_x = labelSegmentStartX((float) ((getMaxDistanceToRoot() * getXcorrectionFactor()) + TreePanel.MOVE
                    + getXdistance()), half_box_size, x);
        } else {
            pos_x = labelSegmentStartX(node.getXcoord(), half_box_size, x);
        }
        float pos_y;
        if (!using_visual_font) {
            pos_y = (node.getYcoord() + (getFontMetricsForLargeDefaultFont().getAscent() / down_shift_factor));
        } else {
            pos_y = (node.getYcoord() + (getFontMetrics(g.getFont()).getAscent() / down_shift_factor));
        }
        if ((getControlPanel().getTreeDisplayType() == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM)
                && (node.isExternal() || node.isCollapse())) {
            // in a vertical orientation the leader is drawn separately as geometry (under R, so it renders vertical)
            // by the caller -- drawing it here would slant it, since this text path runs in the tilted 45deg frame
            if (!isVerticalOrientation()) {
                drawConnection(node.getXcoord(), pos_x - x, node.getYcoord(), 5, 20, g);
            }
            if (node.isCollapse()) {
                pos_x -= add;
            }
        }
        final String sb_str = _sb.toString();
        // GUILHEM_BEG ______________
        if (_control_panel.isShowSequenceRelations() && node.getNodeData().isHasSequence()
                && (_query_sequence != null)) {
            x = paintSequenceRelation(g, node, x, half_box_size, pos_x, pos_y, sb_str);
        }
        // GUILHEM_END _____________
        if (sb_str.length() > 0) {
            TreePanel.drawString(sb_str, pos_x, pos_y, g);
        }
        if (_sb.length() > 0) {
            x += labelStringWidth(g, _sb.toString(), using_visual_font, is_in_found_nodes, to_pdf) + 5;
        }
        if ((getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR)
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE)
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED)) {
            if ((getControlPanel().isShowBinaryCharacters() || getControlPanel().isShowBinaryCharacterCounts())
                    && node.getNodeData().isHasBinaryCharacters()) {
                if ((to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite()) {
                    g.setColor(Color.BLACK);
                } else {
                    g.setColor(getTreeColorSet().getBinaryDomainCombinationsColor());
                }
                if (getControlPanel().isShowBinaryCharacters()) {
                    TreePanel.drawString(
                            node.getNodeData().getBinaryCharacters().getPresentCharactersAsStringBuffer()
                                    .toString(),
                            node.getXcoord() + x + 1 + half_box_size,
                            node.getYcoord() + (getFontMetricsForLargeDefaultFont().getAscent()
                                    / down_shift_factor),
                            g);
                    paintGainedAndLostCharacters(g,
                            node,
                            node.getNodeData().getBinaryCharacters()
                                    .getGainedCharactersAsStringBuffer().toString(),
                            node.getNodeData().getBinaryCharacters()
                                    .getLostCharactersAsStringBuffer().toString());
                } else {
                    TreePanel.drawString(" " + node.getNodeData().getBinaryCharacters().getPresentCount(),
                            node.getXcoord() + x + 4 + half_box_size,
                            node.getYcoord() + (getFontMetricsForLargeDefaultFont().getAscent()
                                    / down_shift_factor),
                            g);
                    paintGainedAndLostCharacters(g,
                            node,
                            "+" + node.getNodeData().getBinaryCharacters().getGainedCount(),
                            "-" + node.getNodeData().getBinaryCharacters().getLostCount());
                }
            }
        }
        return x;
    }

    /**
     * Whether {@code node}'s label should use the publication-style placement: to the LEFT of the node,
     * right-aligned, on top of the incoming branch. Applies to non-root internal nodes (not collapsed)
     * in the rectangular-family layouts when the option is on. Nodes that also show binary characters (a
     * rare legacy feature positioned to the right of the node) keep the old right-placement so nothing
     * is dropped.
     */
    private boolean usesAboveBranchInternalLabel(final PhylogenyNode node) {
        if (!getOptions().isInternalLabelsAboveBranch() || node.isExternal() || node.isCollapse()
                || node.isRoot()) {
            return false;
        }
        final PHYLOGENY_GRAPHICS_TYPE t = getPhylogenyGraphicsType();
        if ((t != PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR) && (t != PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE)
                && (t != PHYLOGENY_GRAPHICS_TYPE.ROUNDED)) {
            return false;
        }
        if ((getControlPanel().isShowBinaryCharacters() || getControlPanel().isShowBinaryCharacterCounts())
                && node.getNodeData().isHasBinaryCharacters()) {
            return false;
        }
        // Sequence relations (the query-sequence highlight / ortholog underline) are drawn to the RIGHT
        // of the node by paintSequenceRelation; keep the old placement for such nodes so that rendering
        // is not dropped -- same rationale as the binary-characters exclusion above.
        if (getControlPanel().isShowSequenceRelations() && (_query_sequence != null)
                && node.getNodeData().isHasSequence()) {
            return false;
        }
        return true;
    }

    /**
     * Draws an internal node's textual label (taxonomy + node-data, e.g. a clade name from "Annotate
     * Clades by Rank") to the LEFT of the node, right-aligned so it ends just left of the node and sits
     * on top of the incoming branch -- the publication-style placement. Confidence values and branch
     * lengths are drawn elsewhere and unaffected. Returns 0: the label grows leftward, consuming no
     * space to the right of the node (so any right-side renderable data starts at the node).
     */
    private int paintInternalLabelAboveBranch(final Graphics2D g,
                                              final PhylogenyNode node,
                                              final boolean is_in_found_nodes,
                                              final boolean to_pdf,
                                              final boolean to_graphics_file,
                                              final int half_box_size) {
        final boolean using_visual_font = setFont(g, node);
        // Each segment is right-aligned to the node, so strip edge whitespace: a measured trailing space
        // (nodeTaxonomyDataAsSB appends one; a property's asText() might) would otherwise push the visible
        // glyphs left of the node. Unlike the classic left-anchored path, trailing space is not harmless.
        String taxo = "";
        if ((getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyScientificNames()
                || getControlPanel().isShowTaxonomyCommonNames() || getControlPanel().isShowTaxonomyRank())
                && node.getNodeData().isHasTaxonomy()) {
            _sb.setLength(0);
            nodeTaxonomyDataAsSB(node.getNodeData().getTaxonomy(), _sb);
            taxo = _sb.toString().trim();
        }
        _sb.setLength(0);
        nodeDataAsSB(node, _sb);
        final String data = _sb.toString().trim();
        if ((taxo.length() == 0) && (data.length() == 0)) {
            return 0;
        }
        // Measure the taxonomy segment via taxonomyLabel (draw=false) -- the same italic-aware, per-part
        // path that draws it below -- so the right-alignment width matches what is painted.
        final int taxo_w = (taxo.length() > 0)
                ? taxonomyLabel(g, node.getNodeData().getTaxonomy(), 0, 0, to_pdf, false) : 0;
        final int data_w = labelStringWidth(g, data, using_visual_font, is_in_found_nodes, to_pdf);
        final int font_descent = using_visual_font ? getFontMetrics(g.getFont()).getMaxDescent()
                : getFontMetricsForLargeDefaultFont().getMaxDescent();
        final float[] layout = TreePanelUtil.internalLabelAboveBranchLayout(node.getXcoord(),
                node.getYcoord(), half_box_size, taxo_w, data_w, INTERNAL_LABEL_SEGMENT_GAP, font_descent,
                INTERNAL_LABEL_MIN_LEFT_MARGIN);
        if (data.length() > 0) {
            setColor(g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getSequenceColor());
            TreePanel.drawString(data, layout[1], layout[2], g);
        }
        if (taxo.length() > 0) {
            setColor(g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getTaxonomyColor());
            taxonomyLabel(g, node.getNodeData().getTaxonomy(), layout[0], layout[2], to_pdf, true);
        }
        return 0;
    }

    /**
     * Advance width of a label string under the same font selection node painting uses everywhere: the
     * fractional-metrics advance for the vector-export path, the large default font otherwise (unless a
     * per-node visual font or a found-node font is in effect). Single source of truth for the width logic
     * in paintTaxonomy, paintNodeData and paintInternalLabelAboveBranch.
     */
    private int labelStringWidth(final Graphics2D g,
                                 final String str,
                                 final boolean using_visual_font,
                                 final boolean is_in_found_nodes,
                                 final boolean to_pdf) {
        if (str.length() == 0) {
            return 0;
        }
        if (to_pdf) {
            return fractionalAdvanceWidth(g, str);
        }
        if (!using_visual_font && !is_in_found_nodes) {
            return getFontMetricsForLargeDefaultFont().stringWidth(str);
        }
        return getFontMetrics(g.getFont()).stringWidth(str);
    }

    private final int paintSequenceRelation(final Graphics2D g,
                                            final PhylogenyNode node,
                                            int x,
                                            final int half_box_size,
                                            final float pos_x,
                                            final float pos_y,
                                            final String sb_str) {
        int nodeTextBoundsWidth = 0;
        if (sb_str.length() > 0) {
            final Rectangle2D node_text_bounds = new TextLayout(sb_str, g.getFont(), _frc).getBounds(); //would like to remove this 'new', but how...
            nodeTextBoundsWidth = (int) node_text_bounds.getWidth();
        }
        if (node.getNodeData().getSequence().equals(_query_sequence)) {
            if (nodeTextBoundsWidth > 0) { // invert font color and background color to show that this is the query sequence
                g.fillRect((int) pos_x - 1, (int) pos_y - 8, nodeTextBoundsWidth + 5, 11);
                g.setColor(getTreeColorSet().getBackgroundColor());
            }
        } else {
            final List<SequenceRelation> seqRelations = node.getNodeData().getSequence().getSequenceRelations();
            for (final SequenceRelation seqRelation : seqRelations) {
                final boolean fGotRelationWithQuery = (seqRelation.getRef0().isEqual(_query_sequence)
                        || seqRelation.getRef1().isEqual(_query_sequence))
                        && seqRelation.getType()
                        .equals(getControlPanel().getSequenceRelationTypeBox().getSelectedItem());
                if (fGotRelationWithQuery) { // we will underline the text to show that this sequence is ortholog to the query
                    final double linePosX = node.getXcoord() + 2 + half_box_size;
                    final String sConfidence = (!getControlPanel().isShowSequenceRelationConfidence()
                            || (seqRelation.getConfidence() == null)) ? null
                            : " (" + seqRelation.getConfidence().getValue() + ")";
                    if (sConfidence != null) {
                        float confidenceX = pos_x;
                        if (sb_str.length() > 0) {
                            confidenceX += new TextLayout(sb_str, g.getFont(), _frc).getBounds().getWidth()
                                    + CONFIDENCE_LEFT_MARGIN;
                        }
                        if (confidenceX > linePosX) { // let's only display confidence value if we are already displaying at least one of Prot/Gene Name and Taxonomy Code
                            final int confidenceWidth = (int) new TextLayout(sConfidence, g.getFont(), _frc)
                                    .getBounds().getWidth();
                            TreePanel.drawString(sConfidence, confidenceX, pos_y, g);
                            x += CONFIDENCE_LEFT_MARGIN + confidenceWidth;
                        }
                    }
                    if ((x + nodeTextBoundsWidth) > 0) /* we only underline if there is something displayed */ {
                        if (nodeTextBoundsWidth == 0) {
                            nodeTextBoundsWidth -= 3; /* the gap between taxonomy code and node name should not be underlined if nothing comes after it */
                        } else {
                            nodeTextBoundsWidth += 2;
                        }
                        g.drawLine((int) linePosX + 1,
                                3 + (int) pos_y,
                                (int) linePosX + x + nodeTextBoundsWidth,
                                3 + (int) pos_y);
                        break;
                    }
                }
            }
        }
        return x;
    }

    private final void drawConnection(final float x1,
                                      final float x2,
                                      final float y,
                                      final int dist_left,
                                      final int dist_right,
                                      final Graphics2D g) {
        if (((x1 + dist_left) < (x2 - dist_right))) {
            final Stroke stroke = g.getStroke();
            final Color saved_color = g.getColor();
            g.setStroke(LEADER_STROKE);
            // The connector is a neutral structural guide linking the (short) branch tip to the lined-up
            // label. It must NOT inherit g's current color, which here is the node's label color: for a
            // search-found node that color is the bright red/green highlight, which would otherwise render
            // the guide as a saturated red/green "branch".
            g.setColor(connectorColor());
            drawLine(x1 + dist_left, y, x2 - dist_right, y, g);
            g.setStroke(stroke);
            g.setColor(saved_color);
        }
    }

    /**
     * The color of the lined-up-data connector: a fixed, neutral guide gray that is deliberately
     * independent of the node's label/highlight color, so the search-found highlight never bleeds
     * into the connector (which exports as a colored "branch").
     */
    static Color connectorColor() {
        return CONNECTOR_GUIDE_COLOR;
    }

    private final void addLabelForCollapsed(final String first,
                                            final String last,
                                            final int size,
                                            final PhylogenyNode node) {
        _sb.append(first.length() < AptxConstants.MAX_LENGTH_FOR_COLLAPSED_NAME ? first
                : first.substring(0, AptxConstants.MAX_LENGTH_FOR_COLLAPSED_NAME - 1));
        _sb.append(" ... ");
        _sb.append(last.length() < AptxConstants.MAX_LENGTH_FOR_COLLAPSED_NAME ? last
                : last.substring(0, AptxConstants.MAX_LENGTH_FOR_COLLAPSED_NAME - 1));
        _sb.append(" (" + size + ")");
        if ((_found_nodes_0 != null) || (_found_nodes_1 != null)) {
            /////
            /////
            final int[] res = calcFoundNodesInSubtree(node);
            if (res[0] > 0) {
                _sb.append(" [");
                _sb.append(res[0]);
                _sb.append("/");
                _sb.append(res[1]);
                _sb.append("]");
            }
        }
    }

    /**
     * How many of a collapsed clade's (hidden) external tips are currently selected/found, so the triangle can
     * signal it: {@code {hits in found set 0, hits in found set 1, tips in either set, total external tips}}.
     * External tips only (selection semantics are tip-based); a tip in both sets counts once toward index 2.
     */
    final int[] collapsedCladeFoundCounts(final PhylogenyNode node) {
        final List<PhylogenyNode> tips = node.getAllExternalDescendants();
        int h0 = 0, h1 = 0, any = 0;
        for (final PhylogenyNode t : tips) {
            final boolean in0 = (_found_nodes_0 != null) && _found_nodes_0.contains(t.getId());
            final boolean in1 = (_found_nodes_1 != null) && _found_nodes_1.contains(t.getId());
            if (in0) {
                ++h0;
            }
            if (in1) {
                ++h1;
            }
            if (in0 || in1) {
                ++any;
            }
        }
        return new int[] { h0, h1, any, tips.size() };
    }

    private final int[] calcFoundNodesInSubtree(final PhylogenyNode node) {
        final List<PhylogenyNode> all_descs = PhylogenyMethods.getAllDescendants(node);
        final int res[] = new int[2];
        int found = 0;
        int total = 0;
        for (final PhylogenyNode desc : all_descs) {
            if (desc.isHasNodeData()) {
                if (((_found_nodes_0 != null) && _found_nodes_0.contains(desc.getId()))
                        || ((_found_nodes_1 != null) && _found_nodes_1.contains(desc.getId()))) {
                    ++found;
                }
                ++total;
            }
        }
        res[0] = found;
        res[1] = total;
        return res;
    }


    final private void paintNodeDataUnrootedCirc(final Graphics2D g,
                                                 final PhylogenyNode node,
                                                 final boolean to_pdf,
                                                 final boolean to_graphics_file,
                                                 final boolean radial_labels,
                                                 final double ur_angle,
                                                 final boolean is_in_found_nodes) {
        if (isNodeDataInvisibleUnrootedCirc(node) && !to_graphics_file && !to_pdf) {
            return;
        }
        _sb.setLength(0);
        _sb.append(" ");
        if (node.getNodeData().isHasTaxonomy()
                && (getControlPanel().isShowTaxonomyCode() || getControlPanel().isShowTaxonomyScientificNames()
                || getControlPanel().isShowTaxonomyCommonNames())) {
            final Taxonomy taxonomy = node.getNodeData().getTaxonomy();
            if (_control_panel.isShowTaxonomyCode() && !ForesterUtil.isEmpty(taxonomy.getTaxonomyCode())) {
                _sb.append(taxonomy.getTaxonomyCode());
                _sb.append(" ");
            }
            if (_control_panel.isShowTaxonomyScientificNames() && _control_panel.isShowTaxonomyCommonNames()) {
                if (!ForesterUtil.isEmpty(taxonomy.getScientificName())
                        && !ForesterUtil.isEmpty(taxonomy.getCommonName())) {
                    _sb.append(taxonomy.getScientificName());
                    _sb.append(" (");
                    _sb.append(taxonomy.getCommonName());
                    _sb.append(") ");
                } else if (!ForesterUtil.isEmpty(taxonomy.getScientificName())) {
                    _sb.append(taxonomy.getScientificName());
                    _sb.append(" ");
                } else if (!ForesterUtil.isEmpty(taxonomy.getCommonName())) {
                    _sb.append(taxonomy.getCommonName());
                    _sb.append(" ");
                }
            } else if (_control_panel.isShowTaxonomyScientificNames()) {
                if (!ForesterUtil.isEmpty(taxonomy.getScientificName())) {
                    _sb.append(taxonomy.getScientificName());
                    _sb.append(" ");
                }
            } else if (_control_panel.isShowTaxonomyCommonNames()) {
                if (!ForesterUtil.isEmpty(taxonomy.getCommonName())) {
                    _sb.append(taxonomy.getCommonName());
                    _sb.append(" ");
                }
            }
        }
        if (node.isCollapse() && ((!node.isRoot() && !node.getParent().isCollapse()) || node.isRoot())) {
            _sb.append(" [");
            _sb.append(node.getAllExternalDescendants().size());
            _sb.append("]");
        }
        if (getControlPanel().isShowNodeNames() && (node.getName().length() > 0)) {
            if (_sb.length() > 0) {
                _sb.append(" ");
            }
            _sb.append(node.getName());
        }
        if (node.getNodeData().isHasSequence()) {
            if (getControlPanel().isShowSequenceAcc()
                    && (node.getNodeData().getSequence().getAccession() != null)) {
                if (_sb.length() > 0) {
                    _sb.append(" ");
                }
                if (!ForesterUtil.isEmpty(node.getNodeData().getSequence().getAccession().getSource())) {
                    _sb.append(node.getNodeData().getSequence().getAccession().getSource());
                    _sb.append(":");
                }
                _sb.append(node.getNodeData().getSequence().getAccession().getValue());
            }
            if (getControlPanel().isShowSeqNames() && (node.getNodeData().getSequence().getName().length() > 0)) {
                if (_sb.length() > 0) {
                    _sb.append(" ");
                }
                _sb.append(node.getNodeData().getSequence().getName());
            }
        }
        //g.setFont( getTreeFontSet().getLargeFont() );
        //if ( is_in_found_nodes ) {
        //    g.setFont( getTreeFontSet().getLargeFont().deriveFont( Font.BOLD ) );
        // }
        if (_sb.length() > 1) {
            setColor(g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getSequenceColor());
            final boolean using_visual_font = setFont(g, node);
            final String sb_str = _sb.toString();
            double m = 0;
            if (_graphics_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) {
                m = _urt_nodeid_angle_map.get(node.getId()) % TWO_PI;
            } else {
                m = (float) (ur_angle % TWO_PI);
            }
            _at = g.getTransform();
            boolean need_to_reset = false;
            final float x_coord = node.getXcoord();
            float y_coord;
            if (!using_visual_font) {
                y_coord = node.getYcoord() + (getFontMetricsForLargeDefaultFont().getAscent() / 3.0f);
            } else {
                y_coord = node.getYcoord() + (getFontMetrics(g.getFont()).getAscent() / 3.0f);
            }
            if (radial_labels) {
                need_to_reset = true;
                boolean left = false;
                if ((m > HALF_PI) && (m < ONEHALF_PI)) {
                    m -= PI;
                    left = true;
                }
                g.rotate(m, x_coord, node.getYcoord());
                if (left) {
                    if (!using_visual_font) {
                        g.translate(-(getFontMetricsForLargeDefaultFont().getStringBounds(sb_str, g).getWidth()),
                                0);
                    } else {
                        g.translate(-(getFontMetrics(g.getFont()).getStringBounds(sb_str, g).getWidth()), 0);
                    }
                }
            } else {
                if ((m > HALF_PI) && (m < ONEHALF_PI)) {
                    need_to_reset = true;
                    if (!using_visual_font) {
                        g.translate(-getFontMetricsForLargeDefaultFont().getStringBounds(sb_str, g).getWidth(), 0);
                    } else {
                        g.translate(-getFontMetrics(g.getFont()).getStringBounds(sb_str, g).getWidth(), 0);
                    }
                }
            }
            TreePanel.drawString(sb_str, x_coord, y_coord, g);
            if (need_to_reset) {
                g.setTransform(_at);
            }
        }
    }

    final private void paintNodeLite(final Graphics2D g, final PhylogenyNode node) {
        if (node.isCollapse()) {
            return;
        }
        if (isInFoundNodes(node)) {
            g.setColor(getColorForFoundNode(node));
            drawRectFilled(node.getXSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                    node.getYSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                    OVERVIEW_FOUND_NODE_BOX_SIZE,
                    OVERVIEW_FOUND_NODE_BOX_SIZE,
                    g);
        }
        float new_x = 0;
        if (!node.isExternal() && !node.isCollapse()) {
            boolean first_child = true;
            float y2 = 0.0f;
            //final int parent_max_branch_to_leaf = getMaxBranchesToLeaf( node );
            // mirror the "Reverse Tip Order" reversal used by paintNodeRectangular, so the overview thumbnail stays
            // consistent with the (flipped) main canvas -- the overview lays out its own YSecondary coords
            final boolean flip = getOptions().isReverseTipOrder();
            final int child_count = node.getNumberOfDescendants();
            for (int i = 0; i < child_count; ++i) {
                final PhylogenyNode child_node = node.getChildNode(flip ? ((child_count - 1) - i) : i);
                final int factor_x = node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes();
                if (first_child) {
                    first_child = false;
                    y2 = node.getYSecondary() - (getOvYDistance()
                            * (node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes()));
                } else {
                    y2 += getOvYDistance() * child_node.getNumberOfExternalNodes();
                }
                final float x2 = calculateOvBranchLengthToParent(child_node, factor_x);
                new_x = x2 + node.getXSecondary();
                final float diff_y = node.getYSecondary() - y2;
                final float diff_x = node.getXSecondary() - new_x;
                if ((diff_y > 2) || (diff_y < -2) || (diff_x > 2) || (diff_x < -2)) {
                    paintBranchLite(g, node.getXSecondary(), new_x, node.getYSecondary(), y2, child_node);
                }
                child_node.setXSecondary(new_x);
                child_node.setYSecondary(y2);
                y2 += getOvYDistance() * child_node.getNumberOfExternalNodes();
            }
        }
    }

    final private void paintNodeRectangular(final Graphics2D g,
                                            final PhylogenyNode node,
                                            final boolean to_pdf,
                                            final boolean dynamically_hide,
                                            final int dynamic_hiding_factor,
                                            final boolean to_graphics_file,
                                            final boolean disallow_shortcutting) {
        final boolean is_in_found_nodes = isInFoundNodes(node);
        final boolean vertical = isVerticalOrientation();
        if (node.isCollapse()) {
            if ((!node.isRoot() && !node.getParent().isCollapse())) {
                paintCollapsedNode(g, node, to_graphics_file, to_pdf, is_in_found_nodes);
                // A collapsed clade is still an internal node: show the support (and length) on its incoming branch.
                if (vertical) {
                    paintBranchDataRightVertical(g, node, to_pdf, to_graphics_file);
                } else if (isShowConfidenceValuesForNode(node)) {
                    paintConfidenceValues(g, node, to_pdf, to_graphics_file);
                }
            }
            return;
        }
        if (node.isExternal()) {
            ++_external_node_index;
        }
        // Support and/or branch-length: in a vertical orientation, horizontal, to the RIGHT of the branch as
        // "support length"; in horizontal mode, confidence here and the branch length inside paintNodeData.
        if (vertical) {
            paintBranchDataRightVertical(g, node, to_pdf, to_graphics_file);
        } else if (isShowConfidenceValuesForNode(node)) {
            paintConfidenceValues(g, node, to_pdf, to_graphics_file);
        }
        // Draw a line to root:
        if (node.isRoot() && _phylogeny.isRooted()) {
            paintRootBranch(g, node.getXcoord(), node.getYcoord(), node, to_pdf, to_graphics_file);
        }
        float new_x = 0;
        float new_x_min = Float.MAX_VALUE;
        float min_dist = 1.5f;
        if (!disallow_shortcutting) {
            if (dynamic_hiding_factor > 4000) {
                min_dist = 4;
            } else if (dynamic_hiding_factor > 1000) {
                min_dist = 3;
            } else if (dynamic_hiding_factor > 100) {
                min_dist = 2;
            }
        }
        if (!node.isExternal() && !node.isCollapse()) {
            boolean first_child = true;
            float y2 = 0.0f;
            // "Reverse Tip Order" is a DISPLAY-ONLY vertical mirror: process the children in REVERSE order so the
            // top-to-bottom tip order reverses, without mutating the tree. The stored y-coords become the real
            // flipped positions, so labels, overlays, hit-testing and every export follow for free (no transform).
            final boolean flip = getOptions().isReverseTipOrder();
            final int child_count = node.getNumberOfDescendants();
            for (int i = 0; i < child_count; ++i) {
                final PhylogenyNode child_node = node.getChildNode(flip ? ((child_count - 1) - i) : i);
                final int factor_x = node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes();
                if (first_child) {
                    first_child = false;
                    y2 = node.getYcoord() - (_y_distance
                            * (node.getNumberOfExternalNodes() - child_node.getNumberOfExternalNodes()));
                } else {
                    y2 += _y_distance * child_node.getNumberOfExternalNodes();
                }
                final float x2 = calculateBranchLengthToParent(child_node, factor_x);
                new_x = x2 + node.getXcoord();
                if (dynamically_hide && (x2 < new_x_min)) {
                    new_x_min = x2;
                }
                final float diff_y = node.getYcoord() - y2;
                final float diff_x = node.getXcoord() - new_x;
                if (disallow_shortcutting || (diff_y > min_dist) || (diff_y < -min_dist) || (diff_x > min_dist)
                        || (diff_x < -min_dist)) {
                    paintBranchRectangular(g,
                            node.getXcoord(),
                            new_x,
                            node.getYcoord(),
                            y2,
                            child_node,
                            to_pdf,
                            to_graphics_file);
                }
                child_node.setXcoord(new_x);
                child_node.setYcoord(y2);
                y2 += _y_distance * child_node.getNumberOfExternalNodes();
            }
            paintNodeBox(node.getXcoord(), node.getYcoord(), node, g, to_pdf, to_graphics_file);
        }
        if (!vertical && getControlPanel().isShowMolSequences() && (node.getNodeData().isHasSequence())
                && (node.getNodeData().getSequence().isMolecularSequenceAligned())
                && (!ForesterUtil.isEmpty(node.getNodeData().getSequence().getMolecularSequence()))) {
            paintMolecularSequences(g, node, to_pdf); // deferred in vertical orientation (increment 1)
        }
        if (dynamically_hide && ((node.isExternal() && ((_external_node_index % dynamic_hiding_factor) != 1))
                || (!node.isExternal() && ((new_x_min < 20)
                || ((_y_distance * node.getNumberOfExternalNodes()) < getFontMetricsForLargeDefaultFont()
                .getHeight()))))) {
            return;
        }
        if (vertical) {
            if (node.isExternal()) {
                // tip label tilts 45deg / 90deg at the tip. Aligned phylogram: draw the tip->label leader as GEOMETRY
                // (it rides R, so it renders vertical) rather than inside the tilted text frame (which would slant it),
                // and pivot the tilted label at the aligned column so it sits on the leader's end.
                if (isAlignedTipLabel(node)) {
                    drawConnection(node.getXcoord(), labelTextStartX(node), node.getYcoord(), 5, 20, g);
                }
                final int[] label_w = { 0 }; // capture the label's pixel width to place the domain track past it
                if (effectiveTipLabelDirection() == Options.TIP_LABEL_DIRECTION.HORIZONTAL) {
                    // UPRIGHT tip label, centred under (root-top) / over (root-bottom) the tip -- the cleanest look for
                    // short names / sparse trees, and the one a rotated-bitmap export can't produce (its text tilts).
                    label_w[0] = paintTipLabelHorizontal(g, node, to_graphics_file, to_pdf, is_in_found_nodes);
                } else {
                    withNodeTextFrame(g, labelTextStartX(node), node.getYcoord(), tipLabelAngle(),
                            () -> label_w[0] = paintNodeData(g, node, to_graphics_file, to_pdf, is_in_found_nodes, 0));
                }
                // renderable domain architecture: a per-tip vertical bar just past the label (boxes ride R, no labels).
                // The other renderable overlays (MSA, binary chars, vectors, sequence relations) stay deferred here.
                paintDomainsVertical(g, node, label_w[0], to_pdf, to_graphics_file);
            } else {
                // internal-node label: horizontal, right-aligned, LEFT of the branch midpoint. This path deliberately
                // does NOT route through paintNodeData, so the other node-data overlays it draws are DEFERRED for
                // internal nodes in a vertical orientation (increment 1): renderable domains, molecular sequences,
                // binary characters / counts, and sequence relations. They each need the rotated-frame re-anchoring
                // the tip labels get; niche for root-top/bottom (focused on internal LABELS, small trees). External
                // nodes still draw them via the tilted paintNodeData path above.
                paintInternalLabelLeftVertical(g, node, to_pdf, to_graphics_file);
            }
        } else {
            final int x = paintNodeData(g, node, to_graphics_file, to_pdf, is_in_found_nodes, 0);
            paintNodeWithRenderableData(x, g, node, to_graphics_file, to_pdf);
        }
    }

    /** Draws an external tip's renderable domain architecture in a VERTICAL orientation: the boxes ride the rotation R
     *  into a thin vertical track just past the tilted tip label (domain-name labels suppressed -- they would collide
     *  with neighbouring tips' tracks). Per tip at its own depth, so the track hangs off the tip. No-op unless domains
     *  are shown, external data is shown, and the tip carries a renderable architecture. */
    private void paintDomainsVertical(final Graphics2D g, final PhylogenyNode node, final int label_w,
                                      final boolean to_pdf, final boolean to_graphics_file) {
        if (!getControlPanel().isShowDomainArchitectures() || !getControlPanel().isShowExternalData()) {
            return;
        }
        if (isNodeDataInvisible(node) && !(to_graphics_file || to_pdf)) {
            return;
        }
        if (!node.getNodeData().isHasSequence() || (node.getNodeData().getSequence().getDomainArchitecture() == null)
                || !(node.getNodeData().getSequence()
                        .getDomainArchitecture() instanceof RenderableDomainArchitecture)) {
            return;
        }
        final RenderableDomainArchitecture rds = (RenderableDomainArchitecture) node.getNodeData().getSequence()
                .getDomainArchitecture();
        final int default_height = 7;
        float yd = getYdistance();
        if (getControlPanel().isDynamicallyHideData()) {
            yd = getTreeFontSet().getFontMetricsLarge().getHeight();
        }
        final int hgt = yd < default_height ? ForesterUtil.roundToInt(yd) : default_height;
        rds.setRenderingHeight(hgt > 1 ? hgt : 2);
        // start just past the label's DEPTH footprint -- the same bounding-box projection depthLabelReserve() reserves
        // (label_w*|sin| + lineH*|cos|), so a tilted 45-degree label's track doesn't overlap the label's lower edge and
        // an upright 0-degree label's track clears the one-line-tall label. Anchor at the label column (labelTextStartX),
        // NOT the tip, so an ALIGNED phylogram's track sits past the lined-up labels. Rides R into a vertical bar.
        final double a = tipLabelAngle();
        final int line_h = getFontMetricsForLargeDefaultFont().getHeight();
        final double depth_reach = (label_w * Math.abs(Math.sin(a))) + (line_h * Math.abs(Math.cos(a)));
        final double start_x = labelTextStartX(node) + depth_reach + VERTICAL_DOMAIN_GAP;
        rds.render((float) start_x, node.getYcoord() - (hgt / 2.0f), g, this, to_pdf, false);
    }

    /** Draws an external tip's label UPRIGHT (0 degrees), centred on the tip along the breadth and placed just past it
     *  along the depth (below the tip for root-top, above for root-bottom) -- the horizontal tip-label direction in a
     *  vertical orientation. Returns the label's pixel width (for the domain-track offset). Chrome-style: it re-anchors
     *  to the upright base frame like {@link #withNodeTextFrame}, but centres instead of pivoting at the tip. */
    private int paintTipLabelHorizontal(final Graphics2D g, final PhylogenyNode node, final boolean to_graphics_file,
                                        final boolean to_pdf, final boolean is_in_found_nodes) {
        final int lw = tipLabelTextWidth(node);
        // anchor at the label START along the depth (past the tip dot, or at the aligned column for an aligned
        // phylogram); its breadth is the tip's breadth, so the same point centres the label AND sets the depth
        final Point2D.Double anchor = screenPoint(labelTextStartX(node), node.getYcoord());
        final FontMetrics fm = getFontMetricsForLargeDefaultFont();
        final boolean root_bottom = getOptions().getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_BOTTOM;
        final double gap = 5.0;
        final double anchor_x = anchor.x - (lw / 2.0); // centre the label on the tip's breadth
        // baseline just past the tip along the depth: below it (root-top) or above it (root-bottom)
        final double anchor_y = root_bottom ? (anchor.y - gap - fm.getDescent()) : (anchor.y + gap + fm.getAscent());
        final AffineTransform saved = g.getTransform();
        g.setTransform(_orientation_base_transform);
        g.translate(anchor_x, anchor_y);
        g.translate(-labelTextStartX(node), -node.getYcoord()); // so paintNodeData's own logical start lands at the anchor
        paintNodeData(g, node, to_graphics_file, to_pdf, is_in_found_nodes, 0);
        g.setTransform(saved);
        return lw;
    }

    /** The pixel width of an external tip's drawn text label (node data + taxonomy) at the large default font -- used
     *  to centre an upright horizontal tip label. Measured the same (non-bold) way as {@code _length_of_longest_text_only}
     *  (the breadth reserve), so the centring and the reserve agree. NOTE: like the reserve, it does not account for a
     *  Bold-Found / visual-style font, so such a label centres a few px off (the documented non-bold-reservation class);
     *  matching the reserve avoids centring a bold label WIDER than the space reserved for it (which would clip). */
    private int tipLabelTextWidth(final PhylogenyNode node) {
        final StringBuilder sb = new StringBuilder();
        nodeDataAsSB(node, sb);
        int w = getFontMetricsForLargeDefaultFont().stringWidth(sb.toString());
        if (node.getNodeData().isHasTaxonomy()) {
            w += taxonomyLabelWidth(node.getNodeData().getTaxonomy(), getTreeFontSet().getLargeFont());
        }
        return w;
    }

    final private void paintNodeWithRenderableData(final int x,
                                                   final Graphics2D g,
                                                   final PhylogenyNode node,
                                                   final boolean to_graphics_file,
                                                   final boolean to_pdf) {
        if (isNodeDataInvisible(node) && !(to_graphics_file || to_pdf)) {
            return;
        }
        if ((!getControlPanel().isShowInternalData() && !node.isExternal())) {
            return;
        }
        if ((!getControlPanel().isShowExternalData() && node.isExternal())) {
            return;
        }
        if (getControlPanel().isShowDomainArchitectures() && node.getNodeData().isHasSequence()
                && (node.getNodeData().getSequence().getDomainArchitecture() != null) && (node.getNodeData()
                .getSequence().getDomainArchitecture() instanceof RenderableDomainArchitecture)) {
            RenderableDomainArchitecture rds = null;
            try {
                rds = (RenderableDomainArchitecture) node.getNodeData().getSequence().getDomainArchitecture();
            } catch (final ClassCastException cce) {
                cce.printStackTrace();
            }
            if (rds != null) {
                final int default_height = 7;
                float y = getYdistance();
                if (getControlPanel().isDynamicallyHideData()) {
                    y = getTreeFontSet().getFontMetricsLarge().getHeight();
                }
                final int h = y < default_height ? ForesterUtil.roundToInt(y) : default_height;
                rds.setRenderingHeight(h > 1 ? h : 2);
                if (getControlPanel().isDrawPhylogram()) {
                    if (getOptions().isLineUpRendarableNodeData()) {
                        if (getOptions().isRightLineUpDomains()) {
                            rds.render((float) ((getMaxDistanceToRoot() * getXcorrectionFactor())
                                            + _length_of_longest_text + _phylogeny.getRoot().getXcoord()
                                            + ((_longest_domain - rds.getTotalLength()) * rds.getRenderingFactorWidth())),
                                    node.getYcoord() - (h / 2.0f),
                                    g,
                                    this,
                                    to_pdf);
                        } else {
                            rds.render((float) ((getMaxDistanceToRoot() * getXcorrectionFactor())
                                            + _length_of_longest_text + _phylogeny.getRoot().getXcoord()),
                                    node.getYcoord() - (h / 2.0f),
                                    g,
                                    this,
                                    to_pdf);
                        }
                    } else {
                        rds.render(node.getXcoord() + x, node.getYcoord() - (h / 2.0f), g, this, to_pdf);
                    }
                } else {
                    if (getOptions().isRightLineUpDomains()) {
                        rds.render(((getPhylogeny().getFirstExternalNode().getXcoord() + _length_of_longest_text)
                                        - 20) + ((_longest_domain - rds.getTotalLength()) * rds.getRenderingFactorWidth()),
                                node.getYcoord() - (h / 2.0f),
                                g,
                                this,
                                to_pdf);
                    } else {
                        rds.render(getPhylogeny().getFirstExternalNode().getXcoord()
                                + _length_of_longest_text, node.getYcoord() - (h / 2.0f), g, this, to_pdf);
                    }
                }
            }
        }
        if (getControlPanel().isShowVectorData() && (node.getNodeData().getVector() != null)
                && (node.getNodeData().getVector().size() > 0) && (getStatisticsForExpressionValues() != null)) {
            final RenderableVector rv = RenderableVector.createInstance(node.getNodeData().getVector(),
                    getStatisticsForExpressionValues(),
                    getConfiguration());
            if (rv != null) {
                double domain_add = 0;
                if (getControlPanel().isShowDomainArchitectures() && node.getNodeData().isHasSequence()
                        && (node.getNodeData().getSequence().getDomainArchitecture() != null)) {
                    domain_add = _domain_structure_width + 10;
                }
                if (getControlPanel().isDrawPhylogram()) {
                    rv.render((float) (node.getXcoord() + x + domain_add), node.getYcoord() - 3, g, this, to_pdf);
                } else {
                    rv.render((float) (getPhylogeny().getFirstExternalNode().getXcoord() + _length_of_longest_text
                            + domain_add), node.getYcoord() - 3, g, this, to_pdf);
                }
            }
        }
        //if ( getControlPanel().isShowMolSequences() && ( node.getNodeData().isHasSequence() )
        //        && ( node.getNodeData().getSequence().isMolecularSequenceAligned() )
        //        && ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) ) {
        //   paintMolecularSequences( g, node, to_pdf );
        //}
    }

    final private void paintOvRectangle(final Graphics2D g) {
        final float w_ratio = ((float) getWidth()) / getVisibleRect().width;
        final float h_ratio = ((float) getHeight()) / getVisibleRect().height;
        final float x_ratio = ((float) getWidth()) / getVisibleRect().x;
        final float y_ratio = ((float) getHeight()) / getVisibleRect().y;
        final float width = getOvMaxWidth() / w_ratio;
        final float height = getOvMaxHeight() / h_ratio;
        final float x = getVisibleRect().x + getOvXPosition() + (getOvMaxWidth() / x_ratio);
        final float y = getVisibleRect().y + getOvYPosition() + (getOvMaxHeight() / y_ratio);
        g.setColor(getTreeColorSet().getFoundColor0());
        getOvRectangle().setRect(x, y, width, height);
        final Stroke s = g.getStroke();
        g.setStroke(STROKE_1);
        if ((width < 6) && (height < 6)) {
            drawRectFilled(x, y, 6, 6, g);
            getOvVirtualRectangle().setRect(x, y, 6, 6);
        } else if (width < 6) {
            drawRectFilled(x, y, 6, height, g);
            getOvVirtualRectangle().setRect(x, y, 6, height);
        } else if (height < 6) {
            drawRectFilled(x, y, width, 6, g);
            getOvVirtualRectangle().setRect(x, y, width, 6);
        } else {
            drawRect(x, y, width, height, g);
            if (isInOvRect()) {
                drawRect(x + 1, y + 1, width - 2, height - 2, g);
            }
            getOvVirtualRectangle().setRect(x, y, width, height);
        }
        g.setStroke(s);
    }

    final private void paintPhylogenyLite(final Graphics2D g) {
        if (_nodes_in_preorder != null) {
            _phylogeny.getRoot().setXSecondary((float) (getVisibleRect().x + getOvXPosition()
                    + (MOVE / (getVisibleRect().width / getOvRectangle().getWidth()))));
            _phylogeny.getRoot().setYSecondary((getVisibleRect().y + getOvYStart()));
            final Stroke s = g.getStroke();
            g.setStroke(STROKE_05);
            for (final PhylogenyNode element : _nodes_in_preorder) {
                paintNodeLite(g, element);
            }
            g.setStroke(s);
            paintOvRectangle(g);
        }
    }

    /**
     * The overview thumbnail for a VERTICAL orientation. Rather than lay out a separate horizontal mini-tree (which
     * would not match the rotated main canvas), it draws the MAIN tree's branches (their logical coords) through one
     * transform: rotate by R into main-screen space, then scale that whole screen extent down into the overview box.
     * Because {@link #paintOvRectangle} maps the viewport into the box with the SAME main-screen->box scaling, the
     * mini-tree and the navigator rectangle stay aligned automatically. The endpoints are transformed by hand (not
     * via {@code g.transform}) so the thin overview stroke is not scaled down to nothing.
     */
    final private void paintPhylogenyLiteVertical(final Graphics2D g) {
        if ((_nodes_in_preorder == null) || (_orientation_R == null)) {
            return;
        }
        final int w_screen = getWidth();
        final int h_screen = getHeight();
        if ((w_screen <= 0) || (h_screen <= 0)) {
            return;
        }
        final double sx = getOvMaxWidth() / (double) w_screen;
        final double sy = getOvMaxHeight() / (double) h_screen;
        final AffineTransform t = AffineTransform.getTranslateInstance(getVisibleRect().x + getOvXPosition(),
                getVisibleRect().y + getOvYPosition());
        t.scale(sx, sy);
        t.concatenate(_orientation_R);
        final Stroke s = g.getStroke();
        g.setStroke(STROKE_05);
        g.setColor(getTreeColorSet().getOvColor());
        final Point2D.Float pa = new Point2D.Float();
        final Point2D.Float pb = new Point2D.Float();
        final Point2D.Float pc = new Point2D.Float();
        for (final PhylogenyNode node : _nodes_in_preorder) {
            // a COLLAPSED node still draws its OWN incoming branch (the elbow to the triangle apex) -- only its hidden
            // descendants (isHiddenUnderCollapse) and the root (no parent) are skipped, matching the horizontal overview
            if (node.isRoot() || isHiddenUnderCollapse(node)) {
                continue;
            }
            final PhylogenyNode parent = node.getParent();
            // rectangular elbow at the node's logical coords: vertical connector at parent x + horizontal leg at node y
            pa.setLocation(parent.getXcoord(), parent.getYcoord());
            pb.setLocation(parent.getXcoord(), node.getYcoord());
            pc.setLocation(node.getXcoord(), node.getYcoord());
            t.transform(pa, pa);
            t.transform(pb, pb);
            t.transform(pc, pc);
            drawLine(pa.x, pa.y, pb.x, pb.y, g);
            drawLine(pb.x, pb.y, pc.x, pc.y, g);
        }
        // mark found/selected nodes, like paintNodeLite
        if ((getFoundNodes0() != null) || (getFoundNodes1() != null)) {
            for (final PhylogenyNode node : _nodes_in_preorder) {
                if (isInFoundNodes(node) && !node.isCollapse() && !isHiddenUnderCollapse(node)) {
                    pa.setLocation(node.getXcoord(), node.getYcoord());
                    t.transform(pa, pa);
                    g.setColor(getColorForFoundNode(node));
                    drawRectFilled(pa.x - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF, pa.y - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                            OVERVIEW_FOUND_NODE_BOX_SIZE, OVERVIEW_FOUND_NODE_BOX_SIZE, g);
                }
            }
        }
        g.setStroke(s);
        paintOvRectangle(g);
    }

    /**
     * Paint the root branch. (Differs from others because it will always be a
     * single horizontal line).
     *
     * @param to_graphics_file
     * @return new x1 value
     */
    final private void paintRootBranch(final Graphics2D g,
                                       final float x1,
                                       final float y1,
                                       final PhylogenyNode root,
                                       final boolean to_pdf,
                                       final boolean to_graphics_file) {
        assignGraphicsForBranchWithColorForParentBranch(root, false, g, to_pdf, to_graphics_file);
        float d = getXdistance();
        if (getControlPanel().isDrawPhylogram() && (root.getDistanceToParent() > 0.0)) {
            d = (float) (getXcorrectionFactor() * root.getDistanceToParent());
        }
        if (d < MIN_ROOT_LENGTH) {
            d = MIN_ROOT_LENGTH;
        }
        if (!getControlPanel().isWidthBranches() || (PhylogenyMethods.getBranchWidthValue(root) == 1)) {
            drawLine(x1 - d, root.getYcoord(), x1, root.getYcoord(), g);
        } else {
            final double w = PhylogenyMethods.getBranchWidthValue(root);
            drawRectFilled(x1 - d, root.getYcoord() - (w / 2), d, w, g);
        }
        paintNodeBox(x1, root.getYcoord(), root, g, to_pdf, to_graphics_file);
    }

    final private void paintScale(final Graphics2D g,
                                  int x1,
                                  int y1,
                                  final boolean to_pdf,
                                  final boolean to_graphics_file) {
        if (isVerticalOrientation()) {
            paintScaleVertical(g, x1, y1, to_pdf, to_graphics_file);
            return;
        }
        x1 += MOVE;
        final double x2 = x1 + (getScaleDistance() * getXcorrectionFactor());
        y1 -= 12;
        final int y2 = y1 - 8;
        final int y3 = y1 - 4;
        g.setFont(getTreeFontSet().getSmallFont());
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        final Stroke s = g.getStroke();
        g.setStroke(STROKE_1);
        drawLine(x1, y1, x1, y2, g);
        drawLine(x2, y1, x2, y2, g);
        drawLine(x1, y3, x2, y3, g);
        if (getScaleLabel() != null) {
            g.drawString(getScaleLabel(), (x1 + 2), y3 - 2);
        }
        g.setStroke(s);
    }

    /** The scale bar for a vertical (root-top/bottom) orientation: the depth axis runs vertically, so the bar is
     *  drawn vertically (its pixel length is still distance * X-correction factor) near the bottom-left, with short
     *  horizontal end ticks and the distance label beside it. */
    final private void paintScaleVertical(final Graphics2D g, int x1, final int y1, final boolean to_pdf,
                                          final boolean to_graphics_file) {
        x1 += MOVE;
        final int y_bottom = y1 - 12;
        final int y_top = (int) Math.round(y_bottom - (getScaleDistance() * getXcorrectionFactor()));
        final int x_tick = x1 + 8;
        final int x_bar = x1 + 4;
        g.setFont(getTreeFontSet().getSmallFont());
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        final Stroke s = g.getStroke();
        g.setStroke(STROKE_1);
        drawLine(x1, y_bottom, x_tick, y_bottom, g);
        drawLine(x1, y_top, x_tick, y_top, g);
        drawLine(x_bar, y_bottom, x_bar, y_top, g);
        if (getScaleLabel() != null) {
            g.drawString(getScaleLabel(), x_tick + 3, (y_bottom + y_top) / 2);
        }
        g.setStroke(s);
    }

    /**
     * Faint vertical reference lines at scale-distance intervals across the tree (phylograms only), drawn
     * BEHIND the branches so branch depths can be read off. Lines use the same distance unit as the scale
     * bar, measured from the root (distance 0).
     */
    final private void paintScaleGrid(final Graphics2D g,
                                      final boolean to_pdf,
                                      final boolean to_graphics_file,
                                      final int graphics_file_y,
                                      final int graphics_file_height) {
        final float origin_x = _phylogeny.getRoot().getXcoord();
        final float spacing = (float) (getScaleDistance() * getXcorrectionFactor());
        final float max_x = (float) (origin_x + (getMaxDistanceToRoot() * getXcorrectionFactor()));
        final float[] xs = TreePanelUtil.scaleGridLineXs(origin_x, spacing, max_x);
        if (xs.length == 0) {
            return;
        }
        // Use the export canvas extent only when it is real; the direct printer path (TreePanel.print)
        // passes to_pdf=true with height 0, so fall back to the panel height there as on screen.
        final boolean use_export_extent = (to_pdf || to_graphics_file) && (graphics_file_height > 0);
        // the grid line spans the full CROSS-tree extent (all tips): the device height in a horizontal layout, but the
        // logical breadth extent in a vertical orientation, where the line is drawn in logical coords and rides R into
        // a horizontal grid line at that depth (spanning the tip breadth).
        final int top = isVerticalOrientation() ? 0 : (use_export_extent ? graphics_file_y : 0);
        final int bottom = isVerticalOrientation() ? treeBreadthExtent()
                : (use_export_extent ? (graphics_file_y + graphics_file_height) : getHeight());
        final Color saved_color = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setColor(scaleGridColor(to_pdf, to_graphics_file));
        g.setStroke(STROKE_05);
        for (final float x : xs) {
            drawLine(x, top, x, bottom, g);
        }
        g.setColor(saved_color);
        g.setStroke(saved_stroke);
    }

    /**
     * A labeled distance axis with tick marks along the bottom (phylograms only): a horizontal line spanning the
     * tree's depth, a tick at each scale-distance interval, and a numeric label under each tick -- so branch lengths
     * can be read off directly. Same distance unit and origin (the root = 0) as the scale bar / grid. Anchored to the
     * TREE bottom (panel height / export extent, like paintScaleGrid), so screen and every export agree (WYSIWYG) and
     * the scale bar + tree name are lifted clear of it. Labels are decimated so they never overlap on a dense tree.
     */
    final private void paintScaleAxis(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file,
                                      final int graphics_file_y, final int graphics_file_height) {
        final double corr = getXcorrectionFactor();
        final float origin_x = _phylogeny.getRoot().getXcoord();
        final double max_dist = getMaxDistanceToRoot();
        final double[] ticks = TreePanelUtil.scaleAxisTickValues(max_dist, getScaleDistance());
        if (ticks.length == 0) {
            return;
        }
        final Font saved_font = g.getFont();
        final Color saved_color = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = g.getFontMetrics();
        final boolean use_export_extent = (to_pdf || to_graphics_file) && (graphics_file_height > 0);
        final int bottom = use_export_extent ? (graphics_file_y + graphics_file_height) : getHeight();
        final int axis_y = bottom - scaleAxisBandHeight();
        final int label_baseline = axis_y + SCALE_AXIS_TICK_LEN + fm.getAscent() + 1;
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        g.setStroke(STROKE_1);
        final int x_end = (int) Math.round(origin_x + (max_dist * corr));
        drawLine(Math.round(origin_x), axis_y, x_end, axis_y, g); // the axis line: root (0) -> deepest tip
        int last_label_right = Integer.MIN_VALUE; // decimate labels (keep all ticks) so they never overlap
        for (final double v : ticks) {
            final int x = (int) Math.round(origin_x + (v * corr));
            drawLine(x, axis_y, x, axis_y + SCALE_AXIS_TICK_LEN, g); // the tick mark (always drawn)
            final String label = TreePanelUtil.formatCompactNumber(v);
            final int half = fm.stringWidth(label) / 2;
            if ((x - half) >= (last_label_right + SCALE_AXIS_LABEL_GAP)) {
                g.drawString(label, x - half, label_baseline);
                last_label_right = x + half;
            }
        }
        // the distance unit once, at the right end -- only if it clears the last drawn tick label
        if (!ForesterUtil.isEmpty(_phylogeny.getDistanceUnit())) {
            final int unit_x = x_end + SCALE_AXIS_UNIT_GAP;
            if (unit_x >= (last_label_right + SCALE_AXIS_LABEL_GAP)) {
                g.drawString("[" + _phylogeny.getDistanceUnit() + "]", unit_x, label_baseline);
            }
        }
        g.setFont(saved_font);
        g.setColor(saved_color);
        g.setStroke(saved_stroke);
    }

    /** The labeled scale axis in a VERTICAL orientation: a ruler down one BREADTH side (the band just past the last
     *  tip that {@link #verticalScaleAxisReserve()} reserves), with tick marks pointing toward the tree and UPRIGHT
     *  numeric labels. Chrome -- drawn after the base frame is restored; each depth position maps to a device point
     *  via {@link #screenPoint(double, double)} (which applies R), so the same distance scale the branches use places
     *  the ticks. The ruler is on the left for ROOT_TOP / on the right for ROOT_BOTTOM (the tree's own handedness). */
    private void paintScaleAxisVertical(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        if (_orientation_R == null) {
            return; // R is built during the geometry pass; before the first vertical paint screenPoint is logical
        }
        final double corr = getXcorrectionFactor();
        if (corr <= 0.0) {
            return;
        }
        final double max_dist = getMaxDistanceToRoot();
        final double[] ticks = TreePanelUtil.scaleAxisTickValues(max_dist, getScaleDistance());
        if (ticks.length == 0) {
            return;
        }
        final float origin_x = _phylogeny.getRoot().getXcoord();
        // the ruler sits at the OUTER edge of the reserved breadth band (2 px inside it). Deriving it from the breadth
        // EXTENT -- not from max(visible tip Ycoord) -- keeps it robust to collapsed clades (a collapsed triangle is an
        // internal row counted by treeBreadthExtent() but absent from visibleExternalTips()) and needs no tip walk.
        final double ruler_ly = treeBreadthExtent() - 2.0;
        final Font saved_font = g.getFont();
        final Color saved_color = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = g.getFontMetrics();
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        g.setStroke(STROKE_1);
        // the ruler runs along the depth (device-y varies) at a constant device-x (the breadth); the ticks/labels
        // point INWARD toward the tree -- derive that device-x sign from a point just inside the ruler (tip side)
        final Point2D.Double r_root = screenPoint(origin_x, ruler_ly);
        final Point2D.Double r_tip = screenPoint(origin_x + (max_dist * corr), ruler_ly);
        final int ruler_x = (int) Math.round(r_root.x);
        final int tip_side_x = (int) Math.round(screenPoint(origin_x, ruler_ly - 16.0).x);
        final int in = (tip_side_x >= ruler_x) ? 1 : -1; // +1 = tree to the right of the ruler (ticks point right)
        drawLine(ruler_x, (int) Math.round(r_root.y), ruler_x, (int) Math.round(r_tip.y), g); // the axis line
        final int center = (fm.getAscent() - fm.getDescent()) / 2; // baseline offset to vertically centre on a tick
        final int min_label_gap = fm.getHeight() + SCALE_AXIS_LABEL_GAP;
        int last_label_y = Integer.MIN_VALUE; // decimate labels (keep every tick) so they never overlap along depth
        for (final double v : ticks) {
            final int ty = (int) Math.round(screenPoint(origin_x + (v * corr), ruler_ly).y);
            drawLine(ruler_x, ty, ruler_x + (in * SCALE_AXIS_TICK_LEN), ty, g); // tick toward the tree
            final String label = TreePanelUtil.formatCompactNumber(v);
            if ((last_label_y == Integer.MIN_VALUE) || (Math.abs(ty - last_label_y) >= min_label_gap)) {
                final int lx = (in > 0) ? (ruler_x + SCALE_AXIS_TICK_LEN + 2)
                        : (ruler_x - SCALE_AXIS_TICK_LEN - 2 - fm.stringWidth(label));
                g.drawString(label, lx, ty + center); // upright, vertically centred on the tick
                last_label_y = ty;
            }
        }
        // (The horizontal axis appends a "[unit]" label at its far end; the vertical ruler deliberately omits it -- the
        // deepest-tip end is where the outermost tip label sits, and the flush side ruler leaves no clear room for it,
        // so a long unit string would overprint the tip label. The numeric ticks convey the scale on their own.)
        g.setFont(saved_font);
        g.setColor(saved_color);
        g.setStroke(saved_stroke);
    }

    /** Vertical space (px) the labeled scale axis occupies below its top (line + ticks + one label row). One source,
     *  used to place the axis AND to lift the scale bar / tree name clear of it. */
    private int scaleAxisBandHeight() {
        return getFontMetrics(getTreeFontSet().getSmallFont()).getHeight() + SCALE_AXIS_TICK_LEN + 4;
    }

    /** True when a labeled scale axis is applicable to the CURRENT layout: a rectangular-family PHYLOGRAM (not a
     *  cladogram, not the radial CIRCULAR/UNROOTED layouts) with a positive scale distance that yields at least one
     *  tick -- i.e. the axis is (or, if the toggle were on, would be) actually drawn and reserves tip-spread space.
     *  Independent of the on/off state, so it also governs whether toggling "Scale Axis" needs a re-fit. */
    boolean scaleAxisAppliesToLayout() {
        if (!getControlPanel().isDrawPhylogram() || (getScaleDistance() <= 0.0)
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            return false;
        }
        return TreePanelUtil.scaleAxisTickValues(getMaxDistanceToRoot(), getScaleDistance()).length > 0;
    }

    /** Breadth (px) reserved at the BOTTOM for the HORIZONTAL scale-axis band (line + ticks + one label row), so the
     *  bottommost tip is not overlapped by the axis on a fit-to-window dense tree. The horizontal-orientation analog of
     *  {@link #verticalScaleAxisReserve()} (which reserves a side band); the two are mutually exclusive by orientation.
     *  Reserved in the breadth budget (calcParametersForPainting) AND the breadth extent (treeBreadthExtent), so paint,
     *  fit, and scroll agree and the ruler sits clear just below the last tip. 0 in a vertical orientation, and 0
     *  wherever the axis is not drawn (cladogram / CIRCULAR / UNROOTED / no ticks) -- see scaleAxisAppliesToLayout. */
    private int scaleAxisBottomReserve() {
        if (isVerticalOrientation() || !getOptions().isShowScaleAxis() || !scaleAxisAppliesToLayout()) {
            return 0;
        }
        return scaleAxisBandHeight();
    }

    /** Breadth (px) reserved for the VERTICAL scale-axis ruler (the tick + the widest tick number + padding), so the
     *  vertical ruler + its upright labels sit clear of the tips just past the last tip. 0 unless the axis is shown in
     *  a vertical orientation on a phylogram with a scale. Reserved in the breadth budget (calcParametersForPainting)
     *  and the logical breadth extent (logicalTreeExtent), so paint, fit, and scroll all agree. */
    private int verticalScaleAxisReserve() {
        if (!isVerticalOrientation() || !getOptions().isShowScaleAxis() || !getControlPanel().isDrawPhylogram()
                || (getScaleDistance() <= 0.0)) {
            return 0;
        }
        final double[] ticks = TreePanelUtil.scaleAxisTickValues(getMaxDistanceToRoot(), getScaleDistance());
        if (ticks.length == 0) {
            return 0; // no ticks -> paintScaleAxisVertical draws nothing, so reserve nothing (gates stay in step)
        }
        final FontMetrics fm = getFontMetrics(getTreeFontSet().getSmallFont());
        int max_label = 0;
        for (final double v : ticks) {
            max_label = Math.max(max_label, fm.stringWidth(TreePanelUtil.formatCompactNumber(v)));
        }
        return SCALE_AXIS_TICK_LEN + max_label + 8;
    }

    /** Ink color shared by the scale bar / scale axis / tree name: black in a B&W export, else the branch-length color. */
    private Color scaleInkColor(final boolean to_pdf, final boolean to_graphics_file) {
        return ((to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite())
                ? Color.BLACK
                : getTreeColorSet().getBranchLengthColor();
    }

    /** Faint, theme-aware color for the scale grid: the background nudged slightly toward the branch-length color. */
    private Color scaleGridColor(final boolean to_pdf, final boolean to_graphics_file) {
        if ((to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite()) {
            return new Color(220, 220, 220);
        }
        return TreePanelUtil.blend(getTreeColorSet().getBackgroundColor(),
                getTreeColorSet().getBranchLengthColor(), SCALE_GRID_BLEND);
    }

    final private void paintTreeName(final Graphics2D g,
                                     final int region_x,
                                     final int region_width,
                                     int y1,
                                     final boolean to_pdf,
                                     final boolean to_graphics_file,
                                     final boolean align_right,
                                     final boolean raise_for_scale_axis) {
        g.setFont(getTreeFontSet().getSmallFont());
        if (raise_for_scale_axis) {
            // the labeled scale axis occupies the bottom band -- lift the name clear above it (same band height the
            // axis reserves, so the two can't drift), never past the top edge
            y1 = Math.max(g.getFontMetrics().getHeight(), y1 - scaleAxisBandHeight());
        }
        y1 -= 12;
        final int y3 = y1 - 4;
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        final Stroke s = g.getStroke();
        g.setStroke(STROKE_1);
        final String name = getPhylogeny().getName();
        // Lower-left by default; when the scale occupies the lower-left, right-align the name into the
        // lower-right so the two never overlap. Whether to draw at all is the caller's decision (the
        // isShowTreeName option) -- the old "getScaleLabel() != null" guard here was a copy-paste leftover from
        // paintScale that made the tree name display depend on the scale-label computation.
        final int left_x = region_x + MOVE + 2;
        // right-align, but never left of the left margin: a name too long for the region -- or a zero-width
        // print region (File > Print passes width 0) -- falls back to the left instead of running off-canvas
        // or to a negative x
        final int x = align_right
                ? Math.max(left_x, region_x + region_width - g.getFontMetrics().stringWidth(name) - MOVE)
                : left_x;
        g.drawString(name, x, y3 - 2);
        g.setStroke(s);
    }

    final private int paintTaxonomy(final Graphics2D g,
                                    final PhylogenyNode node,
                                    final boolean is_in_found_nodes,
                                    final boolean to_pdf,
                                    final boolean to_graphics_file,
                                    final float x_shift) {
        final Taxonomy taxonomy = node.getNodeData().getTaxonomy();
        final boolean using_visual_font = setFont(g, node);
        setColor(g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getTaxonomyColor());
        float start_x = labelSegmentStartX(node.getXcoord(), getOptions().getDefaultNodeShapeSize() / 2, x_shift);
        if ((getControlPanel().getTreeDisplayType() == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM)
                && node.isExternal()) {
            start_x = labelSegmentStartX((float) ((getMaxDistanceToRoot() * getXcorrectionFactor()) + TreePanel.MOVE
                    + getXdistance()), getOptions().getDefaultNodeShapeSize() / 2, x_shift);
        }
        float start_y;
        if (!using_visual_font) {
            start_y = node.getYcoord() + (getFontMetricsForLargeDefaultFont().getAscent()
                    / (node.getNumberOfDescendants() == 1 ? 1 : 3.0f));
        } else {
            start_y = node.getYcoord()
                    + (getFontMetrics(g.getFont()).getAscent() / (node.getNumberOfDescendants() == 1 ? 1 : 3.0f));
        }
        /* GUILHEM_BEG */
        if (_control_panel.isShowSequenceRelations() && (node.getNodeData().isHasSequence())
                && node.getNodeData().getSequence().equals(_query_sequence)) {
            // invert font color and background color to show that this is the query sequence
            final int label_w = taxonomyLabel(g, taxonomy, start_x, start_y, to_pdf, false);
            if (label_w > 0) {
                g.fillRect((int) start_x - 1, (int) start_y - 8, label_w + 4, 11);
                g.setColor(getTreeColorSet().getBackgroundColor());
            }
        }
        /* GUILHEM_END */
        return taxonomyLabel(g, taxonomy, start_x, start_y, to_pdf, true);
    }

    /**
     * Measures and (when {@code draw}) paints a taxonomy label starting at {@code start_x} /
     * {@code baseline_y}, italicizing the scientific-name part when {@link Options#isUseItalicScientificNames()}
     * is on; the rest is drawn in g's current font and color. g's font is restored on return. Returns the
     * total advance width. Each part is measured in the very font it is drawn in, so italic and roman
     * advances stay accurate and adjacent parts (and the node-data drawn after the label) do not overlap.
     */
    private int taxonomyLabel(final Graphics2D g, final Taxonomy taxonomy, final float start_x,
                              final float baseline_y, final boolean to_pdf, final boolean draw) {
        final Font base = g.getFont();
        final Font italic = getOptions().isUseItalicScientificNames() ? italicOf(base) : base;
        final int[] w = { 0 };
        forEachTaxonomyLabelPart(taxonomy, (text, scientific) -> {
            if (text.isEmpty()) {
                return;
            }
            // The scientific-name part uses an italic-derived font; the SVG/EPS backend turns all text
            // (italic included) into glyph outlines so the bundled face is not substituted by the viewer
            // (see OutliningVectorGraphics2D). PDF and screen draw with the real font.
            g.setFont(scientific ? italic : base);
            if (draw) {
                TreePanel.drawString(text, start_x + w[0], baseline_y, g);
            }
            w[0] += to_pdf ? fractionalAdvanceWidth(g, text) : getFontMetrics(g.getFont()).stringWidth(text);
        });
        g.setFont(base);
        return w[0];
    }

    /**
     * Width (px) of a taxonomy label measured part-by-part in the very font each part is drawn in -- the
     * scientific-name part in italics when that option is on -- so layout reservation
     * ({@link #calculateLongestExtNodeInfo}, {@link #calcLengthOfLongestText}) matches what
     * {@link #taxonomyLabel} actually paints. Graphics-free (component-level integer metrics, the same the
     * rest of the layout calc uses); {@code base} is the node's large or visual font.
     */
    private int taxonomyLabelWidth(final Taxonomy taxonomy, final Font base) {
        final Font italic = getOptions().isUseItalicScientificNames() ? italicOf(base) : base;
        final int[] w = { 0 };
        forEachTaxonomyLabelPart(taxonomy, (text, scientific) -> {
            if (!text.isEmpty()) {
                w[0] += getFontMetrics(scientific ? italic : base).stringWidth(text);
            }
        });
        return w[0];
    }

    /** Set by AptxUtil around a PNG export so paintPhylogeny leaves the background transparent. */
    void setExportTransparentBackground(final boolean transparent) {
        _export_transparent_background = transparent;
    }

    /**
     * X at which a node label segment starts: an anchor ({@code base_x} -- the node's x, or the aligned-phylogram
     * right margin) plus half the node-shape size, {@link #LABEL_GAP_AFTER_NODE_SHAPE}, and the total width of any
     * label segment already laid out to its left ({@code prior_width}). The taxonomy label passes the x consumed so
     * far; the node-data label passes that SAME x plus the taxonomy width, so when both segments share the same
     * anchor the node data begins exactly where the taxonomy label ends, not a pixel inside it (which let an italic
     * scientific name's right overhang overlap the following node name). NOTE this abutment holds only when the two
     * segments pass the same base_x: for a collapsed node in aligned-phylogram mode the taxonomy uses the node's x
     * while the node data uses the aligned right margin (pre-existing), so they do not abut there. Shared by the
     * taxonomy and node-data segments; other label bits (e.g. paintSequenceRelation) still position independently.
     */
    static float labelSegmentStartX(final float base_x, final int half_box_size, final float prior_width) {
        return base_x + half_box_size + LABEL_GAP_AFTER_NODE_SHAPE + prior_width;
    }

    /** The italic-derived variant of {@code base}, cached so repeated paints don't re-allocate the Font. */
    private Font italicOf(final Font base) {
        if (base != _italic_base_font) {
            _italic_base_font = base;
            _italic_derived_font = base.deriveFont(base.getStyle() | Font.ITALIC);
        }
        return _italic_derived_font;
    }

    // ADDS bold to the base font's existing style (so a visual-style italic becomes bold-italic, not plain bold),
    // cached like italicOf so a large found-set doesn't allocate a Font per label per repaint.
    private Font boldOf(final Font base) {
        if (base != _bold_base_font) {
            _bold_base_font = base;
            _bold_derived_font = base.deriveFont(base.getStyle() | Font.BOLD);
        }
        return _bold_derived_font;
    }

    final private void paintUnrooted(final PhylogenyNode n,
                                     final double low_angle,
                                     final double high_angle,
                                     final boolean radial_labels,
                                     final Graphics2D g,
                                     final boolean to_pdf,
                                     final boolean to_graphics_file) {
        if (n.isRoot()) {
            n.setXcoord(getWidth() / 2);
            n.setYcoord(getHeight() / 2);
        }
        if (n.isExternal()) {
            paintNodeDataUnrootedCirc(g,
                    n,
                    to_pdf,
                    to_graphics_file,
                    radial_labels,
                    (high_angle + low_angle) / 2,
                    isInFoundNodes(n));
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
        for (int i = 0; i < n.getNumberOfDescendants(); ++i) {
            final PhylogenyNode desc = n.getChildNode(i);
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
            final double arc_size = (desc_num_enclosed / num_enclosed) * (high_angle - low_angle);
            float length;
            if (isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram()) {
                if (desc.getDistanceToParent() < 0) {
                    length = 0;
                } else {
                    length = (float) (desc.getDistanceToParent() * getUrtFactor());
                }
            } else {
                length = getUrtFactor();
            }
            final double mid_angle = current_angle + (arc_size / 2);
            final float new_x = (float) (x + (Math.cos(mid_angle) * length));
            final float new_y = (float) (y + (Math.sin(mid_angle) * length));
            desc.setXcoord(new_x);
            desc.setYcoord(new_y);
            paintUnrooted(desc, current_angle, current_angle + arc_size, radial_labels, g, to_pdf, to_graphics_file);
            current_angle += arc_size;
            assignGraphicsForBranchWithColorForParentBranch(desc, false, g, to_pdf, to_graphics_file);
            drawLine(x, y, new_x, new_y, g);
            paintNodeBox(new_x, new_y, desc, g, to_pdf, to_graphics_file);
        }
        if (n.isRoot()) {
            paintNodeBox(n.getXcoord(), n.getYcoord(), n, g, to_pdf, to_graphics_file);
        }
    }

    final private void paintUnrootedLite(final PhylogenyNode n,
                                         final double low_angle,
                                         final double high_angle,
                                         final Graphics2D g,
                                         final float urt_ov_factor) {
        if (n.isRoot()) {
            final int x_pos = (int) (getVisibleRect().x + getOvXPosition() + (getOvMaxWidth() / 2));
            final int y_pos = (int) (getVisibleRect().y + getOvYPosition() + (getOvMaxHeight() / 2));
            n.setXSecondary(x_pos);
            n.setYSecondary(y_pos);
        }
        if (n.isExternal()) {
            return;
        }
        final float num_enclosed = n.getNumberOfExternalNodes();
        final float x = n.getXSecondary();
        final float y = n.getYSecondary();
        double current_angle = low_angle;
        for (int i = 0; i < n.getNumberOfDescendants(); ++i) {
            final PhylogenyNode desc = n.getChildNode(i);
            final int desc_num_enclosed = desc.getNumberOfExternalNodes();
            final double arc_size = (desc_num_enclosed / num_enclosed) * (high_angle - low_angle);
            float length;
            if (isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram()) {
                if (desc.getDistanceToParent() < 0) {
                    length = 0;
                } else {
                    length = (float) (desc.getDistanceToParent() * urt_ov_factor);
                }
            } else {
                length = urt_ov_factor;
            }
            final double mid_angle = current_angle + (arc_size / 2);
            final float new_x = (float) (x + (Math.cos(mid_angle) * length));
            final float new_y = (float) (y + (Math.sin(mid_angle) * length));
            desc.setXSecondary(new_x);
            desc.setYSecondary(new_y);
            if (isInFoundNodes(desc)) {
                g.setColor(getColorForFoundNode(desc));
                drawRectFilled(desc.getXSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                        desc.getYSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                        OVERVIEW_FOUND_NODE_BOX_SIZE,
                        OVERVIEW_FOUND_NODE_BOX_SIZE,
                        g);
                g.setColor(getTreeColorSet().getOvColor());
            }
            paintUnrootedLite(desc, current_angle, current_angle + arc_size, g, urt_ov_factor);
            current_angle += arc_size;
            drawLine(x, y, new_x, new_y, g);
        }
    }

    final private void pasteSubtree(final PhylogenyNode node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            errorMessageNoCutCopyPasteInUnrootedDisplay();
            return;
        }
        if ((getCutOrCopiedTree() == null) || getCutOrCopiedTree().isEmpty()) {
            JOptionPane.showMessageDialog(this,
                    "No tree in buffer (need to copy or cut a subtree first)",
                    "Attempt to paste with empty buffer",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        final String label = createASimpleTextRepresentationOfANode(getCutOrCopiedTree().getRoot());
        final Object[] options = {"As sibling", "As descendant", "Cancel"};
        final int r = JOptionPane.showOptionDialog(this,
                "How to paste subtree" + label + "?",
                "Paste Subtree",
                JOptionPane.CLOSED_OPTION,
                JOptionPane.QUESTION_MESSAGE,
                null,
                options,
                options[2]);
        boolean paste_as_sibling = true;
        if (r == 1) {
            paste_as_sibling = false;
        } else if (r != 0) {
            return;
        }
        final Phylogeny buffer_phy = getCutOrCopiedTree().copy();
        buffer_phy.setAllNodesToNotCollapse();
        PhylogenyMethods.preOrderReId(buffer_phy);
        buffer_phy.setRooted(true);
        boolean need_to_show_whole = false;
        if (paste_as_sibling) {
            if (node.isRoot()) {
                JOptionPane.showMessageDialog(this,
                        "Cannot paste sibling to root",
                        "Attempt to paste sibling to root",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }
            pushUndoCheckpoint("Paste Subtree");
            buffer_phy.addAsSibling(node);
        } else {
            pushUndoCheckpoint("Paste Subtree");
            if ((node.getNumberOfExternalNodes() == 1) && node.isRoot()) {
                need_to_show_whole = true;
                _phylogeny = buffer_phy;
            } else {
                buffer_phy.addAsChild(node);
            }
        }
        if (getCopiedAndPastedNodes() == null) {
            setCopiedAndPastedNodes(new HashSet<Long>());
        }
        final List<PhylogenyNode> nodes = PhylogenyMethods.obtainAllNodesAsList(buffer_phy);
        final Set<Long> node_ids = new HashSet<>(nodes.size());
        for (final PhylogenyNode n : nodes) {
            node_ids.add(n.getId());
        }
        node_ids.add(node.getId());
        getCopiedAndPastedNodes().addAll(node_ids);
        setNodeInPreorderToNull();
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants(true);
        resetNodeIdToDistToLeafMap();
        setEdited(true);
        if (need_to_show_whole) {
            getControlPanel().showWhole();
        }
        repaint();
    }

    private final StringBuffer propertiesToString(final PhylogenyNode node) {
        // hide internal aptx:* metadata (e.g. the persisted Re-import annotation profile on the root)
        return TreePanelUtil.userVisiblePropertiesText(node.getNodeData().getProperties());
    }

    private void setColor(final Graphics2D g,
                          final PhylogenyNode node,
                          final boolean to_graphics_file,
                          final boolean to_pdf,
                          final boolean is_in_found_nodes,
                          final Color default_color) {
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite();
        Color c;
        if (bw) {
            c = Color.BLACK;
        } else if (is_in_found_nodes) {
            c = getColorForFoundNode(node);
        } else if (isColorByProperty()) {
            c = getPropertyBasedColor(node);
        } else if (getControlPanel().isUseVisualStyles() && (node.getNodeData().getNodeVisualData() != null)
                && (node.getNodeData().getNodeVisualData().getFontColor() != null)) {
            c = node.getNodeData().getNodeVisualData().getFontColor();
        } else if (getOptions().isColorLabelsSameAsParentBranch() && getControlPanel().isUseVisualStyles()
                && (PhylogenyMethods.getBranchColorValue(node) != null)) {
            c = PhylogenyMethods.getBranchColorValue(node);
        } else if (to_pdf) {
            c = Color.BLACK;
        } else {
            c = default_color;
        }
        g.setColor(dimNonMatch(c, is_in_found_nodes, bw));
    }

    // "Dim Non-Matches": while a search/selection with at least one VISIBLE hit is active (_has_visible_found_node,
    // resolved once per paint -- also false when every hit is hidden under a collapse / outside the displayed
    // subtree), fade a NON-hit's label/number toward the background so the hits stand out (WYSIWYG -- screen +
    // exports). Hits and B&W export are never dimmed. Shared by name/taxonomy labels, confidence, and branch length.
    private Color dimNonMatch(final Color c, final boolean is_in_found_nodes, final boolean bw) {
        if (!is_in_found_nodes && !bw && _has_visible_found_node) {
            return TreePanelUtil.blend(c, getTreeColorSet().getBackgroundColor(), DIM_NON_MATCH_FRACTION);
        }
        return c;
    }

    /** The ink color for a text element (branch-length / confidence number): BLACK on any PDF or on a B&W raster
     *  export, else the given scheme color. Centralizes the pdf/black-and-white policy the two number sites share. */
    private Color inkColor(final boolean to_pdf, final boolean to_graphics_file, final Color default_color) {
        return (to_pdf || (to_graphics_file && getOptions().isPrintBlackAndWhite())) ? Color.BLACK : default_color;
    }

    final private void setCopiedAndPastedNodes(final Set<Long> nodeIds) {
        getMainPanel().setCopiedAndPastedNodes(nodeIds);
    }

    final private void setCutOrCopiedTree(final Phylogeny cut_or_copied_tree) {
        getMainPanel().setCutOrCopiedTree(cut_or_copied_tree);
    }

    // Horizontal advance (px) that the vector backends (PDF via OpenPDF glyph outlines, SVG/EPS via
    // VectorGraphics2D <text>) actually use when drawing s in g's current font. They advance glyphs by
    // the font's true *fractional* widths, whereas FontMetrics.stringWidth grid-fits each advance to an
    // integer; at the small leaf-label font that rounding deficit accumulates over a long label into a
    // one-to-two character overlap with the next label (e.g. scientific name running into node name).
    // Raster/screen are unaffected -- AWT grid-fits glyphs to the same integers it measures -- so this
    // is used only on the to_pdf (vector) path, keeping screen and raster output pixel-identical.
    static int fractionalAdvanceWidth(final Graphics2D g, final String s) {
        return (int) Math.ceil(g.getFont().getStringBounds(s, FRC_FRACTIONAL).getWidth());
    }

    private boolean setFont(final Graphics2D g, final PhylogenyNode node) {
        Font visual_font = null;
        if (getControlPanel().isUseVisualStyles() && (node.getNodeData().getNodeVisualData() != null)) {
            visual_font = node.getNodeData().getNodeVisualData().getFont();
            g.setFont(visual_font != null ? visual_font : getTreeFontSet().getLargeFont());
        } else {
            g.setFont(getTreeFontSet().getLargeFont());
        }
        // "Bold Found Labels": found/selected node labels render bold (WYSIWYG; safe in PDF since text is drawn
        // as glyph OUTLINES there, so no stroke-bleed). Off by default -- otherwise highlighting is by color only.
        if (getOptions().isBoldFoundLabels() && isInFoundNodes(node)) {
            g.setFont(boldOf(g.getFont()));
        }
        return visual_font != null;
    }

    final private void setInOv(final boolean in_ov) {
        _in_ov = in_ov;
    }

    final private void setOvMaxHeight(final float ov_max_height) {
        _ov_max_height = ov_max_height;
    }

    final private void setOvMaxWidth(final float ov_max_width) {
        _ov_max_width = ov_max_width;
    }

    final private void setOvXcorrectionFactor(final float f) {
        _ov_x_correction_factor = f;
    }

    final private void setOvXDistance(final float ov_x_distance) {
        _ov_x_distance = ov_x_distance;
    }

    final private void setOvXPosition(final int ov_x_position) {
        _ov_x_position = ov_x_position;
    }

    final private void setOvYDistance(final float ov_y_distance) {
        _ov_y_distance = ov_y_distance;
    }

    final private void setOvYPosition(final int ov_y_position) {
        _ov_y_position = ov_y_position;
    }

    final private void setOvYStart(final int ov_y_start) {
        _ov_y_start = ov_y_start;
    }

    final private void setScaleDistance(final double scale_distance) {
        _scale_distance = scale_distance;
    }

    final private void setScaleLabel(final String scale_label) {
        _scale_label = scale_label;
    }

    private final void setupStroke(final Graphics2D g) {

        float w = _options.getDefaultBranchWidth();


        if (getYdistance() < 0.0001) {
            g.setStroke(STROKE_0025);
        }
        if (getYdistance() < 0.001) {
            g.setStroke(STROKE_005);
        } else if (getYdistance() < 0.01) {
            g.setStroke(STROKE_01);
        } else if (getYdistance() < 0.5) {
            g.setStroke(STROKE_025);
        } else if (getYdistance() < 1) {
            g.setStroke(STROKE_05);
        } else if (getYdistance() < 2) {
            g.setStroke(STROKE_075);
        } else if (_options.getDefaultBranchWidth() > 0) {
            g.setStroke(makeStroke(_options.getDefaultBranchWidth()));
        } else {
            g.setStroke(STROKE_1);
        }
    }

    final private void setUpUrtFactor() {
        final int d = getVisibleRect().width < getVisibleRect().height ? getVisibleRect().width
                : getVisibleRect().height;
        if (isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram()) {
            setUrtFactor((float) (d / (2 * getMaxDistanceToRoot())));
        } else {
            final int max_depth = _circ_max_depth;
            if (max_depth > 0) {
                setUrtFactor(d / (2 * max_depth));
            } else {
                setUrtFactor(d / 2);
            }
        }
        setUrtFactorOv(getUrtFactor());
    }

    final private void setUrtFactor(final float urt_factor) {
        _urt_factor = urt_factor;
    }

    final private void setUrtFactorOv(final float urt_factor_ov) {
        _urt_factor_ov = urt_factor_ov;
    }

    final private void showNodeDataPopup(final MouseEvent e, final PhylogenyNode node) {
        try {
            if ((node.getName().length() > 0)
                    || (node.getNodeData().isHasTaxonomy()
                    && !TreePanelUtil.isTaxonomyEmpty(node.getNodeData().getTaxonomy()))
                    || (node.getNodeData().isHasSequence()
                    && !TreePanelUtil.isSequenceEmpty(node.getNodeData().getSequence()))
                    || (node.getNodeData().isHasDate()) || (node.getNodeData().isHasDistribution())
                    || node.getBranchData().isHasConfidences()) {
                _popup_buffer.setLength(0);
                short lines = 0;
                if (node.getName().length() > 0) {
                    lines++;
                    _popup_buffer.append(node.getName());
                }
                if (node.getNodeData().isHasTaxonomy()
                        && !TreePanelUtil.isTaxonomyEmpty(node.getNodeData().getTaxonomy())) {
                    lines++;
                    boolean enc_data = false;
                    final Taxonomy tax = node.getNodeData().getTaxonomy();
                    if (_popup_buffer.length() > 0) {
                        _popup_buffer.append("\n");
                    }
                    if (!ForesterUtil.isEmpty(tax.getTaxonomyCode())) {
                        _popup_buffer.append("[");
                        _popup_buffer.append(tax.getTaxonomyCode());
                        _popup_buffer.append("]");
                        enc_data = true;
                    }
                    if (!ForesterUtil.isEmpty(tax.getScientificName())) {
                        if (enc_data) {
                            _popup_buffer.append(" ");
                        }
                        _popup_buffer.append(tax.getScientificName());
                        enc_data = true;
                    }
                    if (!ForesterUtil.isEmpty(tax.getCommonName())) {
                        if (enc_data) {
                            _popup_buffer.append(" (");
                        } else {
                            _popup_buffer.append("(");
                        }
                        _popup_buffer.append(tax.getCommonName());
                        _popup_buffer.append(")");
                        enc_data = true;
                    }
                    if (!ForesterUtil.isEmpty(tax.getAuthority())) {
                        if (enc_data) {
                            _popup_buffer.append(" (");
                        } else {
                            _popup_buffer.append("(");
                        }
                        _popup_buffer.append(tax.getAuthority());
                        _popup_buffer.append(")");
                        enc_data = true;
                    }
                    if (!ForesterUtil.isEmpty(tax.getRank())) {
                        if (enc_data) {
                            _popup_buffer.append(" [");
                        } else {
                            _popup_buffer.append("[");
                        }
                        _popup_buffer.append(tax.getRank());
                        _popup_buffer.append("]");
                        enc_data = true;
                    }
                    if (tax.getSynonyms().size() > 0) {
                        if (enc_data) {
                            _popup_buffer.append(" ");
                        }
                        _popup_buffer.append("[");
                        int counter = 1;
                        for (final String syn : tax.getSynonyms()) {
                            if (!ForesterUtil.isEmpty(syn)) {
                                enc_data = true;
                                _popup_buffer.append(syn);
                                if (counter < tax.getSynonyms().size()) {
                                    _popup_buffer.append(", ");
                                }
                            }
                            counter++;
                        }
                        _popup_buffer.append("]");
                    }
                    if (!enc_data) {
                        if ((tax.getIdentifier() != null)
                                && !ForesterUtil.isEmpty(tax.getIdentifier().getValue())) {
                            if (!ForesterUtil.isEmpty(tax.getIdentifier().getProvider())) {
                                _popup_buffer.append("[");
                                _popup_buffer.append(tax.getIdentifier().getProvider());
                                _popup_buffer.append("] ");
                            }
                            _popup_buffer.append(tax.getIdentifier().getValue());
                        }
                    }
                }
                if (node.getNodeData().isHasSequence()
                        && !TreePanelUtil.isSequenceEmpty(node.getNodeData().getSequence())) {
                    lines++;
                    boolean enc_data = false;
                    if (_popup_buffer.length() > 0) {
                        _popup_buffer.append("\n");
                    }
                    final Sequence seq = node.getNodeData().getSequence();
                    if (seq.getAccession() != null) {
                        _popup_buffer.append("[");
                        if (!ForesterUtil.isEmpty(seq.getAccession().getSource())) {
                            _popup_buffer.append(seq.getAccession().getSource());
                            _popup_buffer.append(":");
                        }
                        _popup_buffer.append(seq.getAccession().getValue());
                        _popup_buffer.append("]");
                        enc_data = true;
                    }
                    if (!ForesterUtil.isEmpty(seq.getSymbol())) {
                        if (enc_data) {
                            _popup_buffer.append(" [");
                        } else {
                            _popup_buffer.append("[");
                        }
                        _popup_buffer.append(seq.getSymbol());
                        _popup_buffer.append("]");
                        enc_data = true;
                    }
                    if (!ForesterUtil.isEmpty(seq.getGeneName())) {
                        if (enc_data) {
                            _popup_buffer.append(" [");
                        } else {
                            _popup_buffer.append("[");
                        }
                        _popup_buffer.append(seq.getGeneName());
                        _popup_buffer.append("]");
                        enc_data = true;
                    }
                    if (!ForesterUtil.isEmpty(seq.getName())) {
                        if (enc_data) {
                            _popup_buffer.append(" ");
                        }
                        _popup_buffer.append(seq.getName());
                    }
                }
                if (node.getNodeData().isHasDate()) {
                    lines++;
                    if (_popup_buffer.length() > 0) {
                        _popup_buffer.append("\n");
                    }
                    _popup_buffer.append(node.getNodeData().getDate().asSimpleText());
                }
                if (node.getNodeData().isHasDistribution()) {
                    lines++;
                    if (_popup_buffer.length() > 0) {
                        _popup_buffer.append("\n");
                    }
                    _popup_buffer.append(node.getNodeData().getDistribution().asSimpleText());
                }
                if (node.getBranchData().isHasConfidences()) {
                    final List<Confidence> confs = node.getBranchData().getConfidences();
                    for (final Confidence confidence : confs) {
                        if (!Double.isFinite(confidence.getValue())) {
                            continue; // skip non-finite (NaN/Infinity) confidence values
                        }
                        lines++;
                        if (_popup_buffer.length() > 0) {
                            _popup_buffer.append("\n");
                        }
                        if (!ForesterUtil.isEmpty(confidence.getType())) {
                            _popup_buffer.append("[");
                            _popup_buffer.append(confidence.getType());
                            _popup_buffer.append("] ");
                        }
                        _popup_buffer.append(FORMATTER_CONFIDENCE.format(ForesterUtil
                                .round(confidence.getValue(),
                                        getOptions().getNumberOfDigitsAfterCommaForConfidenceValues())));
                        if ((confidence.getStandardDeviation() != Confidence.CONFIDENCE_DEFAULT_VALUE)
                                && Double.isFinite(confidence.getStandardDeviation())) {
                            _popup_buffer.append(" (sd=");
                            _popup_buffer.append(FORMATTER_CONFIDENCE.format(ForesterUtil
                                    .round(confidence.getStandardDeviation(),
                                            getOptions().getNumberOfDigitsAfterCommaForConfidenceValues())));
                            _popup_buffer.append(")");
                        }
                    }
                }
                if (node.getNodeData().isHasProperties()) {
                    // hide internal aptx:* metadata (e.g. the persisted Re-import annotation profile on the root)
                    final StringBuffer props = TreePanelUtil
                            .userVisiblePropertiesText(node.getNodeData().getProperties());
                    if (props.length() > 0) {
                        if (_popup_buffer.length() > 0) {
                            _popup_buffer.append("\n");
                        }
                        _popup_buffer.append(props);
                    }
                }
                if (_popup_buffer.length() > 0) {
                    // the rollover popup is a canvas overlay, so always match the tree color set
                    // (keeps it consistent with the canvas in both light and dark themes)
                    _rollover_popup.setBorder(BorderFactory.createLineBorder(getTreeColorSet().getBranchColor()));
                    _rollover_popup.setBackground(getTreeColorSet().getBackgroundColor());
                    if (isInFoundNodes0(node) && !isInFoundNodes1(node)) {
                        _rollover_popup.setForeground(getTreeColorSet().getFoundColor0());
                    } else if (!isInFoundNodes0(node) && isInFoundNodes1(node)) {
                        _rollover_popup.setForeground(getTreeColorSet().getFoundColor1());
                    } else if (isInFoundNodes0(node) && isInFoundNodes1(node)) {
                        _rollover_popup.setForeground(getTreeColorSet().getFoundColor0and1());
                    } else {
                        _rollover_popup.setForeground(getTreeColorSet().getSequenceColor());
                    }
                    _rollover_popup.setText(_popup_buffer.toString());
                    _node_desc_popup = PopupFactory.getSharedInstance()
                            .getPopup(null,
                                    _rollover_popup,
                                    e.getLocationOnScreen().x + 10,
                                    e.getLocationOnScreen().y - (lines * 20));
                    _node_desc_popup.show();
                }
            }
        } catch (final Exception ex) {
            // Do nothing.
        }
    }

    final private void showNodeEditFrame(final PhylogenyNode n) {
        if (_node_frame_index < TreePanel.MAX_NODE_FRAMES) {
            // NOTE: node-data edits are NOT yet undoable. NodeEditPanel.writeBack commits fields incrementally
            // (and unconditionally, with no change detection) on every selection change, so a correct checkpoint
            // would have to be a change-gated one inside that legacy commit path; a checkpoint here on mere open
            // would clear the redo stack even when the user only inspects a node. Deferred until writeBack grows
            // change detection. Edits are still safe: writeBack's setEdited(true) clears stale redo (safety net).
            // pop up edit box for single node
            _node_frames[_node_frame_index] = new NodeFrame(n, _phylogeny, this, _node_frame_index, "");
            _node_frame_index++;
        } else {
            JOptionPane.showMessageDialog(this, "too many node windows are open");
        }
    }

    final private void showNodeFrame(final PhylogenyNode n) {
        if (_node_frame_index < TreePanel.MAX_NODE_FRAMES) {
            // pop up edit box for single node
            _node_frames[_node_frame_index] = new NodeFrame(n, _phylogeny, this, _node_frame_index);
            _node_frame_index++;
        } else {
            JOptionPane.showMessageDialog(this, "too many node windows are open");
        }
    }

    final void calcMaxDepth() {
        if (_phylogeny != null) {
            _circ_max_depth = PhylogenyMethods.calculateMaxDepth(_phylogeny);
        }
    }

    /**
     * Set parameters for printing the displayed tree
     */
    final void calcParametersForPainting(final int x_in, final int y_in) {
        // Root-top/bottom: the depth axis (root->tip) is drawn VERTICALLY, so it must fit the viewport HEIGHT, and
        // the breadth axis (tip spread) fits the WIDTH -- the opposite of the horizontal layout. Swapping the two
        // viewport budgets here lets the whole body below (which derives the depth scale from x and the breadth
        // spacing from y) stay unchanged; the paint-time transform R then rotates the logical layout into place.
        final int x = isVerticalOrientation() ? y_in : x_in;
        final int y = isVerticalOrientation() ? x_in : y_in;
        // updateStyle(); not needed?
        if ((_phylogeny != null) && !_phylogeny.isEmpty()) {
            initNodeData();
            calculateLongestExtNodeInfo();
            if ((getLongestExtNodeInfo() > (x * 0.6))
                    && (getTreeFontSet().getLargeFont().getSize() > (2 + TreeFontSet.FONT_SIZE_CHANGE_STEP))) {
                while ((getLongestExtNodeInfo() > (x * 0.7))
                        && (getTreeFontSet().getLargeFont().getSize() > 2)) {
                    getMainPanel().getTreeFontSet().decreaseFontSize(getConfiguration().getMinBaseFontSize(), true);
                    calculateLongestExtNodeInfo();
                }
            } else {
                while ((getLongestExtNodeInfo() < (x * 0.6)) && (getTreeFontSet().getLargeFont()
                        .getSize() <= (getTreeFontSet().getLargeFontMemory().getSize()
                        - TreeFontSet.FONT_SIZE_CHANGE_STEP))) {
                    getMainPanel().getTreeFontSet().increaseFontSize();
                    calculateLongestExtNodeInfo();
                }
            }
            // the overlap auto-fit above may have changed the displayed font size -> reflect it in the slider
            if (getControlPanel() != null) {
                getControlPanel().updateFontSizeSlider();
            }
            _length_of_longest_text = calcLengthOfLongestText(); //~~~
            // resolve AUTO (fit) ONCE for this pass, BEFORE the breadth reserves below read it -- so the reserves and the
            // later paint see the SAME direction (the reserves run pre-setYdistance and would otherwise resolve AUTO
            // against a stale/zero y-distance while the paint resolves it against the fresh one -> clipped edge labels).
            _resolved_auto_tip_dir = (getOptions().getTipLabelDirection() == Options.TIP_LABEL_DIRECTION.AUTO)
                    ? TreePanelUtil.autoTipLabelDirection(2.0 * getYdistance(), _length_of_longest_text_only) : null;
            // the tip-label footprint split by the label tilt: its reach along the depth axis (the x budget) and along
            // the breadth axis (the y budget). In the horizontal orientation these are getLongestExtNodeInfo() and 0,
            // so the fit below is unchanged; in a vertical orientation they let the fit reserve the label on the axis
            // it actually extends along, so the breadth no longer overflows nor the depth over-reserve (see the F fit).
            final int depth_label = depthLabelReserve();
            final int breadth_label = breadthLabelReserve();
            int ext_nodes = _phylogeny.getRoot().getNumberOfExternalNodes();
            final int max_depth = PhylogenyMethods.calculateMaxDepthConsiderCollapsed(_phylogeny) + 1;
            if (ext_nodes == 1) {
                ext_nodes = max_depth;
                if (ext_nodes < 1) {
                    ext_nodes = 1;
                }
            }
            updateOvSizes();
            float xdist = 0;
            float ov_xdist = 0;
            if (!isNonLinedUpCladogram()) {
                xdist = (float) ((x - depth_label - rightMarginExtraWidth() - TreePanel.MOVE)
                        / (ext_nodes + 3.0));
                ov_xdist = (float) (getOvMaxWidth() / (ext_nodes + 3.0));
            } else {
                xdist = ((x - depth_label - rightMarginExtraWidth() - TreePanel.MOVE) / (max_depth + 1));
                ov_xdist = (getOvMaxWidth() / (max_depth + 1));
            }
            // reserve space at the top for the rotated annotation-column headers so they don't overlap the top
            // cells; the tree is compressed into the remaining height and its origin shifted down (see below)
            final int top_reserve = annotationHeaderTopReserve();
            // the labeled scale axis eats into the tip-spread budget so the tips stay clear of it: a side ruler in a
            // vertical orientation (verticalScaleAxisReserve), a bottom band in a horizontal one (scaleAxisBottomReserve)
            // -- the two are mutually exclusive by orientation, and both are 0 when the axis is off. Reserving it here
            // (and in treeBreadthExtent) is what keeps the axis from overlapping the bottommost tip on a dense fit.
            final int axis_reserve = verticalScaleAxisReserve() + scaleAxisBottomReserve();
            float ydist = (float) ((y - TreePanel.MOVE - (2 * verticalBreadthPad()) - top_reserve - breadth_label
                    - axis_reserve) / (ext_nodes * 2.0));
            if (xdist < 0.0) {
                xdist = 0.0f;
            }
            if (ov_xdist < 0.0) {
                ov_xdist = 0.0f;
            }
            if (ydist < 0.0) {
                ydist = 0.0f;
            }
            setXdistance(xdist);
            setYdistance(ydist);
            setOvXDistance(ov_xdist);
            final double height = _phylogeny.calculateHeight(!_options.isCollapsedWithAverageHeigh());
            //final double height = PhylogenyMethods.calculateMaxDepth( _phylogeny );
            if (height > 0) {
                final float corr = (float) ((x - (2.0 * TreePanel.MOVE) - depth_label
                        - rightMarginExtraWidth() - getXdistance()) / height);
                setXcorrectionFactor(corr > 0 ? corr : 0);
                final float ov_corr = (float) ((getOvMaxWidth() - getOvXDistance()) / height);
                setOvXcorrectionFactor(ov_corr > 0 ? ov_corr : 0);
            } else {
                setXcorrectionFactor(0);
                setOvXcorrectionFactor(0);
            }
            _circ_max_depth = max_depth;
            setUpUrtFactor();
            //
            if ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                    && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
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
            //
        }
        invalidateOrientationTransform(); // the layout params (ydistance/xcorrection/...) just changed -> R is stale
    }

    final void calculateLongestExtNodeInfo() {
        if ((_phylogeny == null) || _phylogeny.isEmpty()) {
            return;
        }
        int max_possible_length = ForesterUtil
                .roundToInt((getSize().getWidth() - (2 * MOVE)) * AptxConstants.EXT_NODE_INFO_LENGTH_MAX_RATIO);
        if (max_possible_length < 20) {
            max_possible_length = 20;
        }
        int longest = 30;
        int longest_txt = 0;
        int longest_text_only = 0; // text label (name + taxonomy) width, WITHOUT the domain/vector track widths
        _longest_domain = 0;
        _longest_rendered_domain = 0;
        PhylogenyNode longest_txt_node = _phylogeny.getFirstExternalNode();
        for (final PhylogenyNode node : _phylogeny.getExternalNodes()) {
            int sum = 0;
            if (node.isCollapse()) {
                continue;
            }
            final StringBuilder sb = new StringBuilder();
            nodeDataAsSB(node, sb);
            final boolean has_tax = node.getNodeData().isHasTaxonomy();
            // Combined character count selects the longest-text node; the taxonomy chars are counted here
            // but the taxonomy PIXEL width is measured italic-aware below (not from this roman string).
            int txt = sb.length();
            if (has_tax) {
                final StringBuilder tax_sb = new StringBuilder();
                nodeTaxonomyDataAsSB(node.getNodeData().getTaxonomy(), tax_sb);
                txt += tax_sb.length();
            }
            if (txt > longest_txt) {
                longest_txt = txt;
                longest_txt_node = node;
            }
            boolean use_vis = false;
            final Graphics2D g = (Graphics2D) getGraphics();
            if (getControlPanel().isUseVisualStyles()) {
                use_vis = setFont(g, node);
            }
            final Font base = use_vis ? g.getFont() : getTreeFontSet().getLargeFont();
            sum = (use_vis ? getFontMetrics(g.getFont()) : getFontMetricsForLargeDefaultFont())
                    .stringWidth(sb.toString());
            if (has_tax) {
                sum += taxonomyLabelWidth(node.getNodeData().getTaxonomy(), base);
            }
            if (sum > longest_text_only) {
                longest_text_only = sum; // capture BEFORE the domain/vector/binary track widths are added
            }
            if (getControlPanel().isShowBinaryCharacters() && node.getNodeData().isHasBinaryCharacters()) {
                sum += getFontMetricsForLargeDefaultFont().stringWidth(node.getNodeData().getBinaryCharacters()
                        .getGainedCharactersAsStringBuffer().toString());
            }
            if (getControlPanel().isShowVectorData() && (node.getNodeData().getVector() != null)
                    && (node.getNodeData().getVector().size() > 0)) {
                if (getConfiguration() != null) {
                    sum += getConfiguration().getVectorDataWidth() + 10;
                } else {
                    sum += RenderableVector.VECTOR_DEFAULT_WIDTH + 10;
                }
            }
            if (getControlPanel().isShowDomainArchitectures() && node.getNodeData().isHasSequence()
                    && (node.getNodeData().getSequence().getDomainArchitecture() != null)) {
                // FIXME
                // TODO this might need some clean up
                final DomainArchitecture d = node.getNodeData().getSequence().getDomainArchitecture();
                final double rendered = (_domain_structure_width
                        / ((RenderableDomainArchitecture) d).getOriginalSize().getWidth()) * d.getTotalLength();
                sum += rendered + 10;
                if (d.getTotalLength() > _longest_domain) {
                    _longest_domain = d.getTotalLength();
                }
                if (rendered > _longest_rendered_domain) {
                    _longest_rendered_domain = rendered;
                }
            }
            if (getControlPanel().isShowMolSequences() && (node.getNodeData().isHasSequence())
                    && (node.getNodeData().getSequence().isMolecularSequenceAligned())
                    && (!ForesterUtil.isEmpty(node.getNodeData().getSequence().getMolecularSequence()))) {
                // FIXME
                sum += RenderableMsaSequence.DEFAULT_WIDTH + 30;
            }
            if (sum >= max_possible_length) {
                _longest_ext_node_info = max_possible_length;
                // return; //FIXME why?
            }
            if (sum > longest) {
                longest = sum;
            }
        }
        _ext_node_with_longest_txt_info = longest_txt_node;
        _length_of_longest_text_only = longest_text_only;
        if (longest >= max_possible_length) {
            _longest_ext_node_info = max_possible_length;
        } else {
            _longest_ext_node_info = longest;
        }
        _length_of_longest_text = calcLengthOfLongestText();
    }

    final void calculateScaleDistance() {
        if ((_phylogeny == null) || _phylogeny.isEmpty()) {
            return;
        }
        final double height = getMaxDistanceToRoot();
        if (height > 0) {
            if ((height <= 0.5)) {
                setScaleDistance(0.01);
            } else if (height <= 5.0) {
                setScaleDistance(0.1);
            } else if (height <= 50.0) {
                setScaleDistance(1);
            } else if (height <= 500.0) {
                setScaleDistance(10);
            } else {
                setScaleDistance(100);
            }
        } else {
            setScaleDistance(0.0);
        }
        String scale_label = String.valueOf(getScaleDistance());
        if (!ForesterUtil.isEmpty(_phylogeny.getDistanceUnit())) {
            scale_label += " [" + _phylogeny.getDistanceUnit() + "]";
        }
        setScaleLabel(scale_label);
    }

    final Color calculateTaxonomyBasedColor(final Taxonomy tax) {
        if (getOptions().isColorByTaxonomicGroup()) {
            if (!ForesterUtil.isEmpty(tax.getTaxonomyCode())) {
                boolean ex = false;
                String group = null;
                try {
                    group = TaxonomyUtil.getTaxGroupByTaxCode(tax.getTaxonomyCode());
                } catch (final Exception e) {
                    ex = true;
                }
                if (!ex && !ForesterUtil.isEmpty(group)) {
                    final Color c = ForesterUtil.obtainColorDependingOnTaxonomyGroup(group);
                    if (c != null) {
                        return c;
                    }
                }
            }
            return getTreeColorSet().getTaxonomyColor();
        } else {
            if (ForesterUtil.isEmpty(tax.getTaxonomyCode()) && ForesterUtil.isEmpty(tax.getScientificName())) {
                return getTreeColorSet().getTaxonomyColor();
            }
            Color c = null;
            if (!ForesterUtil.isEmpty(tax.getTaxonomyCode())) {
                c = getControlPanel().getSpeciesColors().get(tax.getTaxonomyCode());
            }
            if ((c == null) && !ForesterUtil.isEmpty(tax.getScientificName())) {
                c = getControlPanel().getSpeciesColors().get(tax.getScientificName());
            }
            if (c == null) {
                if (!ForesterUtil.isEmpty(tax.getTaxonomyCode())) {
                    c = AptxUtil.calculateColorFromString(tax.getTaxonomyCode(), true);
                    getControlPanel().getSpeciesColors().put(tax.getTaxonomyCode(), c);
                } else {
                    c = AptxUtil.calculateColorFromString(tax.getScientificName(), true);
                    getControlPanel().getSpeciesColors().put(tax.getScientificName(), c);
                }
            }
            return c;
        }
    }

    // ---------------------------------------------------------------------------
    // Color leaves on the fly by the value of a chosen phyloXML property.
    // ---------------------------------------------------------------------------
    private PropertyColorScheme _property_color_scheme = null;
    // The selected "Color by" property ref is remembered separately from the scheme so the
    // scheme can be rebuilt for the currently displayed (sub)tree -- and so coloring switches
    // back on when the user returns to a super-tree even if it had no such values in a subtree.
    private String              _color_by_property_ref = null;
    private String              _color_palette_name = PropertyColorScheme.DEFAULT_PALETTE_NAME;
    // "Size by": scale the tip symbol by the value of a chosen numeric phyloXML property (the size counterpart of
    // Color by). Ref remembered separately from the scale so it rebuilds for the displayed (sub)tree, like above.
    private PropertySizeScale   _property_size_scale = null;
    private String              _size_by_property_ref = null;
    // The property-color legend can be dragged: _legend_offset is its top-left relative to the
    // visible area (null = the default top-right corner), so it stays put as the user scrolls.
    // _property_legend_bounds is where it was last drawn on screen (for hit-testing a drag).
    private Point               _legend_offset = null;
    private Rectangle           _property_legend_bounds = null;
    private int                 _legend_grab_dx = 0;
    private int                 _legend_grab_dy = 0;
    // A SECOND, independent legend for "Size by" (own draggable position + last-drawn bounds), so its size key can
    // be shown ALONGSIDE the color/rank legend -- the combined color+size figure needs both keys visible at once.
    private Point               _size_legend_offset = null;
    private Rectangle           _size_legend_bounds = null;
    private boolean             _dragging_size_legend = false; // which legend the active drag moves (shared grab dx/dy)
    // User-assigned per-value colors: ref -> (group key -> color); applied by the color scheme,
    // overriding the automatic palette color, and kept across scheme rebuilds (navigation).
    private final Map<String, Map<String, Color>> _property_color_overrides = new HashMap<>();
    // Layout of the legend's value rows (recorded when drawn) so a click can be mapped to a value.
    private java.util.List<String> _legend_row_labels   = new java.util.ArrayList<>();
    private int                 _legend_rows_top    = 0;
    private int                 _legend_row_height  = 0;
    // In-legend controls (recorded when a categorical legend is drawn; null otherwise): the sort toggle in
    // the title row and the clickable "+N more / show fewer" footer. Hit-tested in handleLegendClick.
    private Rectangle           _legend_sort_toggle_bounds = null;
    private Rectangle           _legend_more_bounds = null;
    // Categorical-legend display prefs (per tab, toggled from the legend itself): sort by count vs A-Z, and
    // how many entries to show before "+N more". Default: most-frequent first, DEFAULT_LEGEND_MAX_ENTRIES cap.
    private boolean             _legend_sort_by_count = true;
    private int                 _legend_max_entries = DEFAULT_LEGEND_MAX_ENTRIES;
    // Memo of the last ordered legend entries so the sort+alloc in orderLegendEntries does not run on every
    // repaint; invalidated (by identity/value comparison) when the source map, counts, sort mode, or cap change.
    private Map<String, Color>   _ordered_legend = null;
    private Map<String, Color>   _ordered_legend_colors_key = null;
    private Map<String, Integer> _ordered_legend_counts_key = null;
    private boolean             _ordered_legend_by_count_key = true;
    private int                 _ordered_legend_max_key = DEFAULT_LEGEND_MAX_ENTRIES;
    // Which legend subject (property ref / rank / column) is currently shown; when it changes, the expand
    // ("show all") state resets so expanding one legend does not leak onto a different, possibly larger one.
    private String              _legend_subject = null;
    // Legend for "Colorize Subtrees via Taxonomic Rank": taxon name -> color (sorted), with a title.
    // Shown (and draggable, via the same machinery as the property legend) when no property-color
    // legend is active. Null when not colorizing by rank.
    private Map<String, Color>  _rank_legend = null;
    private String              _rank_legend_title = null;
    // Per-taxon (visible) tip counts for the rank/clade legend, so its rows show "(N)" and can sort by count
    // like the property legend. Parallel to _rank_legend; may be empty when counts were not computed.
    private Map<String, Integer> _rank_legend_counts = null;
    // Per-rank user color overrides for the rank-colorize / clade-band legend (mirrors
    // _property_color_overrides): rank -> (taxon -> chosen color). Persist across navigation & re-apply.
    private final Map<String, Map<String, Color>> _rank_color_overrides = new HashMap<>();
    // The rank at which branches are currently rank-colorized, or null. Lets a legend recolor re-apply to
    // the branches as well as to clade bands.
    private String              _branch_rank_colorize_rank = null;
    // Clade annotation bands: shaded boxes (behind/over the clade) or right-edge bars by taxonomic rank.
    // See CladeBand / TreePanelUtil.cladeBands; geometry is computed at paint time from the clade's tips.
    enum CLADE_VIS {
        BOXES, BARS, BRACKETS
    }
    private java.util.List<CladeBand> _clade_bands = null;
    private CLADE_VIS                 _clade_bands_mode = CLADE_VIS.BOXES;
    private String                    _clade_bands_rank = null;
    private final static int          CLADE_BOX_ALPHA = 46;
    private final static int          CLADE_BAR_WIDTH = 9;
    private final static int          CLADE_BAR_GAP   = 16;
    // a few px of breathing room past the longest label, so the box/bar clears the text
    private final static int          CLADE_BAND_RIGHT_PAD = 6;
    private final static int          CLADE_BRACKET_TICK = 5;       // end-tick length of the monochrome "]" bracket
    private final static float        CLADE_BRACKET_STROKE = 1.5f;
    // Tip-aligned annotation columns (color strip / heat map / bar / text), drawn right of the labels.
    private java.util.List<AnnotationColumns.ColumnSpec> _annotation_column_specs = null; // the user's selection
    private AnnotationColumns          _annotation_columns = null;                          // built for the current view
    private int[]                     _annotation_col_widths = null;                        // cached per-column pixel widths
    private Font                      _annotation_col_widths_font = null;                   // font they were computed for
    private int                       _annotation_header_top_reserve = 0;                   // cached rotated-header top reserve
    private int                       _focused_annotation_column = -1;                      // header-clicked -> its legend
    private org.forester.archaeopteryx.tools.NodeDataImporter.ImportProfile _last_import_profile = null; // for one-click Re-import Annotations
    private final TreeHistory         _history = new TreeHistory();                         // snapshot-based undo/redo
    private boolean                   _restoring_snapshot = false;                          // guards setEdited during a restore
    private final static int          ANN_COL_GAP  = 5;                                     // gap around each column
    private final static int          ANN_TEXT_MAX = 130;                                   // max width of a text column

    /** The on-screen bounds of the last-drawn property-color legend, or null; used by the drag test. */
    final Rectangle getPropertyLegendBounds() {
        return _property_legend_bounds;
    }

    /** Whether the point of {@code e} is over the (last-drawn) legend box (property-color, rank, or column).
     *  The active-legend test must match the paint dispatch exactly (isColorByProperty, not just a non-null
     *  scheme): a non-null-but-EMPTY scheme -- e.g. after navigating into a subtree where no visible tip has
     *  the property -- draws NO legend, so its stale bounds/controls must not stay clickable. */
    final boolean isOnPropertyLegend(final MouseEvent e) {
        return (isColorByProperty() || hasRankLegend() || hasAnnotationColumnLegend())
                && (_property_legend_bounds != null)
                && _property_legend_bounds.contains(e.getX(), e.getY());
    }

    /** Whether {@code e} is over the (last-drawn) "Size by" legend box -- an independent, separately draggable key. */
    final boolean isOnSizeLegend(final MouseEvent e) {
        return isSizeByProperty() && (_size_legend_bounds != null)
                && _size_legend_bounds.contains(e.getX(), e.getY());
    }

    /** Over either legend (color/rank/column or size) -- used to start a drag / show the move cursor. */
    final boolean isOnAnyLegend(final MouseEvent e) {
        return isOnPropertyLegend(e) || isOnSizeLegend(e);
    }

    final void startLegendDrag(final MouseEvent e) {
        // the size legend is drawn last (on top), so it wins if the two boxes overlap
        _dragging_size_legend = isOnSizeLegend(e);
        final Rectangle b = _dragging_size_legend ? _size_legend_bounds : _property_legend_bounds;
        if (b != null) {
            _legend_grab_dx = e.getX() - b.x;
            _legend_grab_dy = e.getY() - b.y;
            setCursor(MOVE_CURSOR);
        }
    }

    final void dragLegend(final MouseEvent e) {
        final Rectangle b = _dragging_size_legend ? _size_legend_bounds : _property_legend_bounds;
        if (b == null) {
            return;
        }
        final Rectangle vp = getVisibleRect();
        int ox = (e.getX() - _legend_grab_dx) - vp.x;
        int oy = (e.getY() - _legend_grab_dy) - vp.y;
        ox = Math.max(0, Math.min(ox, Math.max(0, vp.width - b.width)));
        oy = Math.max(0, Math.min(oy, Math.max(0, vp.height - b.height)));
        if (_dragging_size_legend) {
            _size_legend_offset = new Point(ox, oy);
        } else {
            _legend_offset = new Point(ox, oy);
        }
        repaint();
    }

    /** A click on the size legend: double-click returns it to its default corner (it has no recolorable rows). */
    final void handleSizeLegendClick(final MouseEvent e) {
        if (e.getClickCount() == 2) {
            _size_legend_offset = null;
            repaint();
        }
    }

    final void endLegendDrag() {
        setCursor(ARROW_CURSOR);
    }

    /** Returns the legend to its default top-right corner. */
    final void resetLegendPosition() {
        _legend_offset = null;
        repaint();
    }

    /** The legend value label under {@code e}, or null (title row, "+N more", or not on the legend). */
    final String legendValueAt(final MouseEvent e) {
        if (!isOnPropertyLegend(e) || _legend_row_labels.isEmpty() || (_legend_row_height <= 0)) {
            return null;
        }
        if (e.getY() < _legend_rows_top) {
            return null; // the title row / padding above the first value row (avoids a negative index truncating to 0)
        }
        final int idx = (e.getY() - _legend_rows_top) / _legend_row_height;
        return (idx < _legend_row_labels.size()) ? _legend_row_labels.get(idx) : null;
    }

    final void handleLegendClick(final MouseEvent e) {
        // in-legend controls take precedence over the value rows and the reset gesture
        if ((_legend_sort_toggle_bounds != null) && _legend_sort_toggle_bounds.contains(e.getX(), e.getY())) {
            _legend_sort_by_count = !_legend_sort_by_count; // toggle between "by count" and "A-Z"
            repaint();
            return;
        }
        if ((_legend_more_bounds != null) && _legend_more_bounds.contains(e.getX(), e.getY())) {
            _legend_max_entries = (_legend_max_entries >= LEGEND_SHOW_ALL) ? DEFAULT_LEGEND_MAX_ENTRIES
                    : LEGEND_SHOW_ALL; // toggle "show all" / "show fewer"
            repaint();
            return;
        }
        final String value = legendValueAt(e);
        if (value != null) {
            if (e.getClickCount() == 1) {
                chooseColorForValue(value); // left-click a value row -> set its color
            }
        } else if (e.getClickCount() == 2) {
            resetLegendPosition(); // double-click off a value row returns the legend to its corner
        }
    }

    /** Test hook: the on-screen bounds of the last-drawn size legend (null if none). */
    Rectangle getSizeLegendBounds() {
        return _size_legend_bounds;
    }

    /** Test hook: draws the "Size by" legend into {@code g} at {@code bounds}. */
    void drawSizeLegendForTest(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        drawSizeLegend(g, bounds, draggable);
    }

    /** Test hook: draws the active legend into {@code g} at {@code bounds}, in on-screen ({@code draggable})
     *  or static-export mode -- so a test can compare the two (exports omit the interactive chips). */
    void drawLegendForTest(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        if (hasAnnotationColumnLegend()) {
            drawAnnotationColumnLegend(g, bounds, draggable);
        } else if (isColorByProperty()) {
            drawPropertyColorLegend(g, bounds, draggable);
        } else if (hasRankLegend()) {
            drawRankLegend(g, bounds, draggable);
        }
    }

    /** Test hooks for the in-legend controls (their last-drawn on-screen bounds, or null). */
    Rectangle legendSortToggleBoundsForTest() {
        return _legend_sort_toggle_bounds;
    }

    Rectangle legendMoreBoundsForTest() {
        return _legend_more_bounds;
    }

    boolean isLegendSortByCount() {
        return _legend_sort_by_count;
    }

    int legendMaxEntries() {
        return _legend_max_entries;
    }

    // Clicking any categorical legend row recolors that value -- property-color, rank-colorize, or
    // clade-band legends all behave the same (the inconsistency was that only the property one did).
    private void chooseColorForValue(final String label) {
        if (hasAnnotationColumnLegend()) {
            return; // a header-focused annotation-column legend is display-only (no per-value recolor)
        }
        if ((_property_color_scheme != null) && (_color_by_property_ref != null)) {
            choosePropertyColorForValue(label);
        } else if (hasRankLegend()) {
            chooseRankColorForValue(label);
        }
    }

    // Opens the shared legend color chooser. Returns 0 = cancel, 1 = set (out[0] = chosen color),
    // 2 = automatic ("Use Automatic Color", enabled only when an override already exists).
    private int showLegendColorDialog(final String label, final Color current, final boolean has_override,
                                      final Color[] out) {
        final javax.swing.JColorChooser chooser = new javax.swing.JColorChooser((current != null) ? current
                : java.awt.Color.GRAY);
        final int[] action = { 0 };
        final javax.swing.JDialog dialog = new javax.swing.JDialog(javax.swing.SwingUtilities.getWindowAncestor(this),
                "Color for \"" + label + "\"", java.awt.Dialog.ModalityType.APPLICATION_MODAL);
        final javax.swing.JButton ok = new javax.swing.JButton("OK");
        final javax.swing.JButton cancel = new javax.swing.JButton("Cancel");
        final javax.swing.JButton auto = new javax.swing.JButton("Use Automatic Color");
        auto.setEnabled(has_override);
        ok.addActionListener(ev -> { action[0] = 1; dialog.dispose(); });
        cancel.addActionListener(ev -> { action[0] = 0; dialog.dispose(); });
        auto.addActionListener(ev -> { action[0] = 2; dialog.dispose(); });
        final javax.swing.JPanel buttons = new javax.swing.JPanel();
        buttons.add(auto);
        buttons.add(cancel);
        buttons.add(ok);
        dialog.getContentPane().add(chooser, java.awt.BorderLayout.CENTER);
        dialog.getContentPane().add(buttons, java.awt.BorderLayout.SOUTH);
        dialog.pack();
        dialog.setLocationRelativeTo(this);
        dialog.setVisible(true);
        out[0] = chooser.getColor();
        return action[0];
    }

    // Property-color legend row: OK assigns the color, "Use Automatic Color" clears a prior override.
    // Overrides are kept per ref so they survive navigation and ref switches.
    private void choosePropertyColorForValue(final String label) {
        final String key = _property_color_scheme.getValueKeys().get(label);
        if (key == null) {
            return;
        }
        final Map<String, Color> per_ref = _property_color_overrides.get(_color_by_property_ref);
        final boolean has_override = (per_ref != null) && per_ref.containsKey(key);
        final Color[] out = new Color[1];
        final int action = showLegendColorDialog(label, _property_color_scheme.getValueColors().get(label),
                has_override, out);
        if (action == 1) {
            setColorOverride(key, out[0]);
        } else if (action == 2) {
            clearColorOverride(key);
        }
    }

    // Rank-colorize / clade-band legend row: same chooser; the override is stored per rank and re-applied
    // to the branches and/or clade bands (see reapplyRankColorization).
    private void chooseRankColorForValue(final String taxon) {
        final String rank = currentRankLegendRank();
        if ((rank == null) || (_rank_legend == null)) {
            return;
        }
        final Map<String, Color> per_rank = _rank_color_overrides.get(rank);
        final boolean has_override = (per_rank != null) && per_rank.containsKey(taxon);
        final Color[] out = new Color[1];
        final int action = showLegendColorDialog(taxon, _rank_legend.get(taxon), has_override, out);
        if (action == 1) {
            setRankColorOverride(rank, taxon, out[0]);
        } else if (action == 2) {
            clearRankColorOverride(rank, taxon);
        }
    }

    /** The rank the current rank/clade legend is for (clade bands take precedence), or null. */
    private String currentRankLegendRank() {
        return hasCladeBands() ? _clade_bands_rank : _branch_rank_colorize_rank;
    }

    private Map<String, Color> rankOverridesFor(final String rank) {
        final Map<String, Color> m = (rank == null) ? null : _rank_color_overrides.get(rank);
        return (m != null) ? m : java.util.Collections.<String, Color> emptyMap();
    }

    void setRankColorOverride(final String rank, final String taxon, final Color color) {
        if (rank == null) {
            return;
        }
        _rank_color_overrides.computeIfAbsent(rank, k -> new HashMap<>()).put(taxon, color);
        reapplyRankColorization();
        setEdited(true);
        repaint();
    }

    void clearRankColorOverride(final String rank, final String taxon) {
        final Map<String, Color> per_rank = (rank == null) ? null : _rank_color_overrides.get(rank);
        if (per_rank != null) {
            per_rank.remove(taxon);
            reapplyRankColorization();
            setEdited(true);
            repaint();
        }
    }

    // Re-derive whatever the rank legend drives, honoring the (just-changed) overrides.
    private void reapplyRankColorization() {
        if (hasCladeBands()) {
            rebuildCladeBands();
        }
        if (_branch_rank_colorize_rank != null) {
            recolorBranchesByRank(_branch_rank_colorize_rank);
        }
    }

    void setColorOverride(final String key, final Color color) {
        Map<String, Color> per_ref = _property_color_overrides.get(_color_by_property_ref);
        if (per_ref == null) {
            per_ref = new HashMap<>();
            _property_color_overrides.put(_color_by_property_ref, per_ref);
        }
        per_ref.put(key, color);
        rebuildPropertyColorScheme();
        setEdited(true);
        repaint();
    }

    void clearColorOverride(final String key) {
        final Map<String, Color> per_ref = _property_color_overrides.get(_color_by_property_ref);
        if (per_ref != null) {
            per_ref.remove(key);
            rebuildPropertyColorScheme();
            setEdited(true);
            repaint();
        }
    }

    // Top-left corner at which to draw the legend box; honors a dragged position on screen AND in
    // exports (PDF/PNG/...), mapping it onto the export target by its viewport fraction.
    private Point legendTopLeft(final Rectangle bounds, final int box_w, final int box_h) {
        return legendTopLeftFor(bounds, getVisibleRect(), _legend_offset, box_w, box_h);
    }

    /**
     * Computes the top-left corner at which to draw the property-color legend. {@code offset} is the
     * user-dragged position within the on-screen {@code viewport} (or {@code null} for the default
     * top-right corner of {@code target}). The dragged position is mapped onto {@code target} -- the
     * viewport on screen, the full page/image on export -- by its fractional position, so a moved
     * legend is honored in PDF/PNG/etc. exports instead of snapping back to the corner. When
     * {@code target} matches the viewport (on-screen rendering and "visible region only" exports)
     * this reproduces the exact dragged pixel position. Tested by {@link LegendPlacementTest}.
     */
    static Point legendTopLeftFor(final Rectangle target, final Rectangle viewport, final Point offset,
                                  final int box_w, final int box_h) {
        if (offset == null) {
            return new Point(Math.max(target.x, (target.x + target.width) - box_w - 10), target.y + 10);
        }
        final double fx = ((viewport.width - box_w) <= 0) ? 0.0
                : clamp01(offset.x / (double) (viewport.width - box_w));
        final double fy = ((viewport.height - box_h) <= 0) ? 0.0
                : clamp01(offset.y / (double) (viewport.height - box_h));
        final int ox = (int) Math.round(fx * Math.max(0, target.width - box_w));
        final int oy = (int) Math.round(fy * Math.max(0, target.height - box_h));
        return new Point(target.x + ox, target.y + oy);
    }

    private static double clamp01(final double v) {
        return (v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v);
    }

    /** Colorize leaves by the given property reference, or turn it off when {@code ref} is empty. */
    void setColorByPropertyRef(final String ref) {
        _color_by_property_ref = ForesterUtil.isEmpty(ref) ? null : ref;
        if (_color_by_property_ref != null) {
            clearRankLegend(); // the property-color legend takes over the (single) legend slot
        }
        rebuildPropertyColorScheme();
    }

    /** Hides the taxonomic-rank legend (after colors are cleared, or another coloring takes over). */
    final void clearRankLegend() {
        _rank_legend = null;
        _rank_legend_counts = null;
        _rank_legend_title = null;
        _branch_rank_colorize_rank = null; // branches are no longer rank-colorized
        repaint();
    }

    /** Whether a taxonomic-rank legend is currently available to show. */
    final boolean hasRankLegend() {
        return (_rank_legend != null) && !_rank_legend.isEmpty();
    }

    /** Test hook: the current legend color for {@code taxon} (rank/clade legend), or null. */
    final Color rankLegendColor(final String taxon) {
        return (_rank_legend == null) ? null : _rank_legend.get(taxon);
    }

    /** Test hook: the tip count recorded for {@code taxon} in the rank/clade legend, or null. */
    final Integer rankLegendCount(final String taxon) {
        return (_rank_legend_counts == null) ? null : _rank_legend_counts.get(taxon);
    }

    String getColorPaletteName() {
        return _color_palette_name;
    }

    /** Resets this panel's per-tab "Color by" AND "Size by" state to a freshly-loaded default -- coloring/sizing OFF,
     *  default palette, and no manual per-value color overrides -- for "Reset to Defaults". Sets the palette field
     *  directly (not via setColorPaletteName) because the shared Options palette has already been reset by the caller. */
    void resetColorStateToDefaults() {
        _color_palette_name = PropertyColorScheme.DEFAULT_PALETTE_NAME;
        _property_color_overrides.clear();
        setColorByPropertyRef( null ); // turns coloring off and rebuilds the scheme (-> null)
        setSizeByPropertyRef( null ); // turns sizing off and rebuilds the scale (-> null), clears the size legend pos
        _legend_offset = null; // also return the color/rank legend to its default corner
        repaint();
    }

    /** Selects the categorical palette (see {@link PropertyColorScheme#paletteNames()}) and recolors. Also records
     *  the choice as the shared default (Options) so new tabs inherit it and it can be persisted. */
    void setColorPaletteName(final String name) {
        if (!ForesterUtil.isEmpty(name) && !name.equals(_color_palette_name)) {
            _color_palette_name = name;
            getOptions().setColorPaletteName(name);
            rebuildPropertyColorScheme();
            repaint();
        }
    }

    /**
     * (Re)builds the property color scheme from the currently displayed (visible) tree for the
     * active "Color by" ref, if any. Called whenever the displayed phylogeny changes -- moving
     * into or out of a subtree, collapsing a clade, or deleting nodes -- so the leaf colors and
     * the legend always describe what is on screen.
     */
    void rebuildPropertyColorScheme() {
        if ((_color_by_property_ref == null) || (_phylogeny == null) || _phylogeny.isEmpty()) {
            _property_color_scheme = null;
        } else {
            _property_color_scheme = new PropertyColorScheme(_phylogeny, _color_by_property_ref,
                    _property_color_overrides.get(_color_by_property_ref), _color_palette_name);
        }
    }

    boolean isColorByProperty() {
        return (_property_color_scheme != null) && !_property_color_scheme.isEmpty();
    }

    PropertyColorScheme getPropertyColorScheme() {
        return _property_color_scheme;
    }

    /** Selects the numeric property to SIZE tip symbols by (or null/empty to turn it off), and rebuilds the scale. */
    void setSizeByPropertyRef(final String ref) {
        _size_by_property_ref = ForesterUtil.isEmpty(ref) ? null : ref;
        if (_size_by_property_ref == null) {
            // turning Size-by OFF: forget the legend's dragged position + stale bounds, so re-enabling later shows
            // it fresh at its default corner (not overlapping a color legend at the old spot)
            _size_legend_offset = null;
            _size_legend_bounds = null;
        }
        rebuildPropertySizeScale();
        repaint();
    }

    /** (Re)builds the size scale from the currently displayed (visible) tree for the active "Size by" ref, if any --
     *  so the size range tracks the on-screen (sub)tree exactly like {@link #rebuildPropertyColorScheme}. */
    void rebuildPropertySizeScale() {
        if ((_size_by_property_ref == null) || (_phylogeny == null) || _phylogeny.isEmpty()) {
            _property_size_scale = null;
        } else {
            _property_size_scale = new PropertySizeScale(_phylogeny, _size_by_property_ref);
        }
    }

    /** Rebuilds BOTH property-driven tip encodings ("Color by" and "Size by") from the currently displayed
     *  (visible) tree. Call this at every site the visible tips change (navigation, prune) so the two encodings
     *  can never drift out of lockstep -- one method means a call site cannot rebuild one but forget the other. */
    void rebuildPropertyDisplays() {
        rebuildPropertyColorScheme();
        rebuildPropertySizeScale();
    }

    boolean isSizeByProperty() {
        return (_property_size_scale != null) && !_property_size_scale.isEmpty();
    }

    PropertySizeScale getPropertySizeScale() {
        return _property_size_scale;
    }

    final Color getPropertyBasedColor(final PhylogenyNode node) {
        final Color c = (_property_color_scheme == null) ? null : _property_color_scheme.colorFor(node);
        return (c != null) ? c : getTreeColorSet().getSequenceColor();
    }

    /** Base tip-symbol diameter (px) that the property dots AND the "Size by" legend's sample dots both scale from --
     *  one source, so the legend's key dots always match the tree's tip dots (see drawPropertyColorDot, drawSizeLegend). */
    float baseDotSize() {
        return getOptions().getDefaultNodeShapeSize() + 3;
    }

    // The single tip symbol drawn for "Color by" and/or "Size by": its COLOR is the "Color by" value (a neutral
    // ink color when only sizing is on) and its DIAMETER is the "Size by" value (the default dot when only coloring
    // is on) -- so the two together make a color+size two-attribute figure on one dot.
    private void drawPropertyColorDot(final Graphics2D g, final PhylogenyNode node) {
        final Color color = (_property_color_scheme == null) ? null : _property_color_scheme.colorFor(node);
        // size the dot by value only when this node actually HAS a size value; a value-less tip in size-only
        // mode draws nothing (like Color-by), so "no data" stays distinct from the smallest value
        final boolean size = isSizeByProperty() && node.isExternal() && _property_size_scale.hasValueFor(node);
        if (!size && (color == null)) {
            return; // nothing to draw
        }
        final float base = baseDotSize();
        final int d = Math.max(1, Math.round(size ? _property_size_scale.diameterFor(node, base) : base));
        final Color saved = g.getColor();
        g.setColor((color != null) ? color : getTreeColorSet().getSequenceColor()); // neutral ink when size-only
        g.fillOval((int) node.getXcoord() - (d / 2), (int) node.getYcoord() - (d / 2), d, d);
        g.setColor(saved);
    }

    // Default number of categorical-legend rows shown before "+N more"; the user can expand (show all) or
    // collapse from the legend itself (see _legend_max_entries / the clickable footer).
    private static final int DEFAULT_LEGEND_MAX_ENTRIES = 20;
    private static final int LEGEND_SHOW_ALL = Integer.MAX_VALUE;
    // The legend is a fixed key, not tree data, so keep it readable even when node-label fonts are
    // set very small: use the small tree font but never below this floor.
    private static final float PROPERTY_LEGEND_MIN_FONT_SIZE = 11f;

    private Font legendFont() {
        final Font small = getTreeFontSet().getSmallFont();
        return small.deriveFont(Math.max(small.getSize2D(), PROPERTY_LEGEND_MIN_FONT_SIZE));
    }

    /**
     * Chooses which categorical-legend entries to show and in what order. Always keeps the {@code max}
     * most-frequent values (so the cap hides the least significant, not an arbitrary alphabetic tail), then
     * orders the survivors for display: by count descending (ties A-Z) when {@code by_count}, else A-Z.
     * {@code counts} may be null/partial (a missing count sorts as 0). Static + package-visible for testing.
     */
    static java.util.LinkedHashMap<String, Color> orderLegendEntries(final Map<String, Color> colors,
            final Map<String, Integer> counts, final int max, final boolean by_count) {
        final java.util.List<String> keys = new java.util.ArrayList<>(colors.keySet());
        keys.sort(new java.util.Comparator<String>() { // rank by count so the cap keeps the most significant

            @Override
            public int compare(final String a, final String b) {
                final int ca = legendCount(counts, a);
                final int cb = legendCount(counts, b);
                return (ca != cb) ? Integer.compare(cb, ca) : String.CASE_INSENSITIVE_ORDER.compare(a, b);
            }
        });
        final java.util.List<String> ordered = new java.util.ArrayList<>((keys.size() > max)
                ? keys.subList(0, Math.max(0, max)) : keys);
        if (!by_count) {
            ordered.sort(String.CASE_INSENSITIVE_ORDER);
        }
        final java.util.LinkedHashMap<String, Color> out = new java.util.LinkedHashMap<>();
        for (final String k : ordered) {
            out.put(k, colors.get(k));
        }
        return out;
    }

    private static int legendCount(final Map<String, Integer> counts, final String key) {
        final Integer c = (counts == null) ? null : counts.get(key);
        return (c == null) ? 0 : c.intValue();
    }

    /**
     * Memoized {@link #orderLegendEntries} for the current sort mode and cap, keyed on the source maps'
     * identity: the legend is repainted every frame (pan/zoom) but its ordering only changes when the scheme
     * is rebuilt (new map) or the user toggles sort / expand -- so the O(n log n) sort + allocation runs then,
     * not on every repaint.
     */
    /**
     * Records the currently-shown legend's subject (a stable id: property ref, rank title, or column header).
     * When it changes -- the user switched to a different legend -- the expand ("show all") cap resets to the
     * default so it does not leak onto the new (possibly much larger) legend. The sort mode stays sticky (a
     * deliberate consistency choice). A no-op across repaints/navigation of the same subject.
     */
    private void noteLegendSubject(final String subject) {
        if ((subject == null) ? (_legend_subject != null) : !subject.equals(_legend_subject)) {
            _legend_max_entries = DEFAULT_LEGEND_MAX_ENTRIES;
            _legend_subject = subject;
        }
    }

    private Map<String, Color> orderedLegend(final Map<String, Color> colors, final Map<String, Integer> counts) {
        if ((_ordered_legend != null) && (colors == _ordered_legend_colors_key)
                && (counts == _ordered_legend_counts_key) && (_ordered_legend_by_count_key == _legend_sort_by_count)
                && (_ordered_legend_max_key == _legend_max_entries)) {
            return _ordered_legend;
        }
        _ordered_legend = orderLegendEntries(colors, counts, _legend_max_entries, _legend_sort_by_count);
        _ordered_legend_colors_key = colors;
        _ordered_legend_counts_key = counts;
        _ordered_legend_by_count_key = _legend_sort_by_count;
        _ordered_legend_max_key = _legend_max_entries;
        return _ordered_legend;
    }

    /**
     * Draws the shared legend chrome -- a constant-1px border around a background-filled box with the title on
     * the first line -- and records the draggable on-screen bounds into the caller's slot ({@code size_legend}
     * selects _size_legend_bounds vs the shared _property_legend_bounds). The stroke is left as STROKE_1 for the
     * caller's own in-box drawing (the caller saves/restores it); returns the title baseline y.
     */
    private int drawLegendBox(final Graphics2D g, final int x, final int y, final int box_w, final int box_h,
                              final int pad, final String title, final FontMetrics fm, final boolean draggable,
                              final boolean size_legend) {
        if (draggable) {
            final Rectangle b = new Rectangle(x, y, box_w, box_h);
            if (size_legend) {
                _size_legend_bounds = b;
            } else {
                _property_legend_bounds = b;
            }
        }
        g.setStroke(STROKE_1);
        final Color fg = getTreeColorSet().getSequenceColor();
        g.setColor(getBackground());
        g.fillRect(x, y, box_w, box_h);
        g.setColor(fg);
        g.drawRect(x, y, box_w, box_h);
        final int baseline = y + pad + fm.getAscent();
        g.drawString(title, x + pad, baseline);
        return baseline;
    }

    /** Clears the categorical-legend row/control hit-test state; used by the gradient/bar legends, which have
     *  no clickable value rows, sort toggle, or expand control. */
    private void clearLegendRowControls() {
        _legend_row_labels = new java.util.ArrayList<>();
        _legend_sort_toggle_bounds = null;
        _legend_more_bounds = null;
    }

    /** Draws a capped value-to-color legend; draggable (recorded for hit-testing) on screen. */
    private void drawPropertyColorLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        if ((_property_color_scheme == null) || _property_color_scheme.isEmpty()) {
            return;
        }
        if (_property_color_scheme.isGradient()) {
            drawPropertyColorGradientLegend(g, bounds, draggable,
                    "Color by: " + PropertyColorScheme.displayName(_property_color_scheme.getRef()),
                    _property_color_scheme);
            return;
        }
        noteLegendSubject("prop:" + _property_color_scheme.getRef());
        final Map<String, Integer> counts = _property_color_scheme.getValueCounts();
        final Map<String, Color> values = orderedLegend(_property_color_scheme.getValueColors(), counts);
        final int more = _property_color_scheme.numberOfValues() - values.size();
        final String title = "Color by: " + PropertyColorScheme.displayName(_property_color_scheme.getRef());
        drawCategoricalLegend(g, bounds, draggable, title, values, counts, more);
    }

    /** Draws the "Colorize Subtrees via Taxonomic Rank" legend (taxon -> color, with tip counts); draggable. */
    private void drawRankLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        if (!hasRankLegend()) {
            return;
        }
        noteLegendSubject("rank:" + _rank_legend_title);
        final Map<String, Color> values = orderedLegend(_rank_legend, _rank_legend_counts);
        drawCategoricalLegend(g, bounds, draggable, _rank_legend_title, values, _rank_legend_counts,
                _rank_legend.size() - values.size());
    }

    /**
     * Draws a boxed title plus swatch/label rows for an ordered value-to-color map. On screen
     * ({@code draggable}) it adds two interactive controls -- a sort toggle ("[by count]"/"[A-Z]") in the
     * title row and a clickable footer that expands ("[show all]") or collapses ("[show fewer]") -- and
     * records their bounds so the shared drag/hit-test machinery can map a click back to a row or control.
     * A static export (not draggable) omits that clickable chrome but keeps the informative "+N more" line so
     * the figure still shows that categories were truncated. Shared by the property-color, taxonomic-rank, and
     * annotation-column legends. {@code counts} may be null (then rows show no count).
     */
    private void drawCategoricalLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable,
                                       final String title, final Map<String, Color> values,
                                       final Map<String, Integer> counts, final int more) {
        final int shown = values.size();
        final int swatch = 10;
        final int gap = 5;
        final int pad = 7;
        final int max_text_px = 200;
        g.setFont(legendFont());
        final FontMetrics fm = g.getFontMetrics();
        final int row_h = fm.getHeight() + 2;
        // The sort toggle and the [show all]/[show fewer] chips are interactive UI, so they are drawn only
        // on screen (draggable); a static export (PDF/SVG/EPS/TIFF/PNG/JPG) keeps the informative "+N more"
        // text but drops the clickable chrome.
        final boolean expanded = _legend_max_entries >= LEGEND_SHOW_ALL;
        final boolean more_footer = more > 0;                            // "+N more" info: on screen AND in exports
        final boolean show_sort = draggable && (shown >= 2);            // sort toggle: on-screen only
        final boolean show_all_chip = draggable && more_footer;         // "[show all]" chip: on-screen only
        final boolean show_fewer = draggable && expanded && !more_footer && (shown > DEFAULT_LEGEND_MAX_ENTRIES);
        final String sort_chip = _legend_sort_by_count ? "[by count]" : "[A-Z]";
        final String more_text = "… +" + more + " more";
        final String show_all_chip_str = "[show all]";
        final String show_fewer_chip = "[show fewer]";
        int text_w = fm.stringWidth(title) + (show_sort ? (gap + fm.stringWidth(sort_chip)) : 0);
        for (final String v : values.keySet()) {
            text_w = Math.max(text_w, swatch + gap + fm.stringWidth(legendRowText(v, counts, fm, max_text_px)));
        }
        if (more_footer) {
            text_w = Math.max(text_w, fm.stringWidth(more_text)
                    + (show_all_chip ? (gap + fm.stringWidth(show_all_chip_str)) : 0));
        } else if (show_fewer) {
            text_w = Math.max(text_w, fm.stringWidth(show_fewer_chip));
        }
        // a few extra px on the right so the longest value clears the border even
        // when PDF font metrics run slightly wider than AWT's stringWidth().
        final int box_w = text_w + (2 * pad) + 4;
        final int box_h = ((1 + shown + ((more_footer || show_fewer) ? 1 : 0)) * row_h) + (2 * pad);
        final Point tl = legendTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = tl.y;
        if (draggable) {
            _legend_row_labels = new java.util.ArrayList<>(values.keySet());
            _legend_rows_top = y + pad + row_h; // first value row starts just below the title row
            _legend_row_height = row_h;
            _legend_sort_toggle_bounds = null; // set below only when actually drawn
            _legend_more_bounds = null;
        }
        // constant 1px stroke for the box outline (drawLegendBox); swatches/chips are fillRect/drawString,
        // so they are already stroke-independent
        final Stroke saved_stroke = g.getStroke();
        final Color fg = getTreeColorSet().getSequenceColor();
        int baseline = drawLegendBox(g, x, y, box_w, box_h, pad, title, fm, draggable, false);
        g.setColor(fg);
        if (show_sort) { // sort toggle, right-aligned in the title row (on-screen only -> draggable)
            final int chip_w = fm.stringWidth(sort_chip);
            final int chip_x = (x + box_w) - pad - chip_w;
            g.drawString(sort_chip, chip_x, baseline);
            _legend_sort_toggle_bounds = new Rectangle(chip_x - 2, y + pad, chip_w + 4, row_h);
        }
        for (final Map.Entry<String, Color> e : values.entrySet()) {
            baseline += row_h;
            g.setColor(e.getValue());
            g.fillRect(x + pad, baseline - fm.getAscent() + ((fm.getAscent() - swatch) / 2) + 1, swatch, swatch);
            g.setColor(fg);
            g.drawString(legendRowText(e.getKey(), counts, fm, max_text_px), x + pad + swatch + gap, baseline);
        }
        if (more_footer) { // "+N more" info (also in exports); the "[show all]" chip is on-screen only
            baseline += row_h;
            g.drawString(more_text, x + pad, baseline);
            if (show_all_chip) {
                final int chip_w = fm.stringWidth(show_all_chip_str);
                final int chip_x = (x + box_w) - pad - chip_w;
                g.drawString(show_all_chip_str, chip_x, baseline);
                _legend_more_bounds = new Rectangle(chip_x - 2, baseline - fm.getAscent(), chip_w + 4, row_h);
            }
        } else if (show_fewer) {
            baseline += row_h;
            g.drawString(show_fewer_chip, x + pad, baseline);
            _legend_more_bounds = new Rectangle(x + pad - 2, baseline - fm.getAscent(),
                    fm.stringWidth(show_fewer_chip) + 4, row_h);
        }
        g.setStroke(saved_stroke);
    }

    /** A legend row: the (clipped) value label followed by its count, e.g. {@code "USA (42)"}. */
    private String legendRowText(final String value, final Map<String, Integer> counts, final FontMetrics fm,
                                 final int max_px) {
        final Integer count = (counts == null) ? null : counts.get(value);
        return clipToWidth(value, fm, max_px) + ((count != null) ? (" (" + count + ")") : "");
    }

    /** Draws a gradient bar legend (low value to high value) for a continuous property. */
    private void drawPropertyColorGradientLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable,
                                                 final String title, final PropertyColorScheme scheme) {
        final int pad = 7;
        final int bar_w = 200;
        final int bar_h = 12;
        g.setFont(legendFont());
        final FontMetrics fm = g.getFontMetrics();
        final String min_lbl = scheme.getGradientMinLabel();
        final String max_lbl = scheme.getGradientMaxLabel();
        final int content_w = Math.max(fm.stringWidth(title), bar_w);
        final int box_w = content_w + (2 * pad) + 4;
        final int box_h = (2 * fm.getHeight()) + bar_h + 6 + (2 * pad);
        final Point tl = legendTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = tl.y;
        if (draggable) {
            clearLegendRowControls(); // a gradient legend has no clickable value rows / sort / expand controls
        }
        // The legend is a fixed UI key, not tree data: draw its borders with a constant 1px stroke
        // rather than inheriting the branch stroke set by setupStroke(), which shrinks to sub-pixel
        // widths when the tree is zoomed out. The gradient bar itself is painted with stroke-
        // independent fillRect columns (see paintGradientBar) so its colors stay saturated at every
        // zoom level. (The legend font is likewise floored independent of the node-label font.)
        final Stroke saved_stroke = g.getStroke();
        final Color fg = getTreeColorSet().getSequenceColor();
        final int baseline = drawLegendBox(g, x, y, box_w, box_h, pad, title, fm, draggable, false);
        final int bar_x = x + pad;
        final int bar_y = baseline + 4;
        paintGradientBar(g, bar_x, bar_y, bar_w, bar_h + 1, t -> scheme.gradientColorAt(t));
        g.setColor(fg);
        g.drawRect(bar_x, bar_y, bar_w - 1, bar_h);
        final int label_baseline = bar_y + bar_h + fm.getAscent() + 2;
        g.drawString(min_lbl, bar_x, label_baseline);
        g.drawString(max_lbl, (bar_x + bar_w) - fm.stringWidth(max_lbl), label_baseline);
        g.setStroke(saved_stroke);
    }

    /**
     * Paints a horizontal left-to-right color bar filling {@code [x, x+w) x [y, y+h)}. Each 1px-wide
     * column is drawn with {@link Graphics2D#fillRect} so its color depends on neither the caller's
     * stroke nor antialiasing -- an axis-aligned, integer-aligned fill always covers whole pixels.
     * That keeps the "Color by:" gradient legend equally saturated at every tree zoom level, instead
     * of washing out when the (zoom-dependent) branch stroke shrinks to sub-pixel widths.
     * {@code colorAt} maps t in [0,1] (left to right) to a column color. Tested by
     * {@link PropertyLegendBarTest}.
     */
    static void paintGradientBar(final Graphics2D g, final int x, final int y, final int w, final int h,
                                 final java.util.function.DoubleFunction<Color> colorAt) {
        for (int i = 0; i < w; ++i) {
            final double t = (w == 1) ? 0.0 : (i / (double) (w - 1));
            g.setColor(colorAt.apply(t));
            g.fillRect(x + i, y, 1, h);
        }
    }

    /**
     * Draws the "Size by" legend/key: a boxed title plus sample dots (min / mid / max value) at their ACTUAL
     * diameters -- the same area-proportional mapping as the tip dots -- each labeled with its value, so a reader
     * (and a published figure) can decode the encoding. The dots are neutral-filled (the size axis is independent of
     * color), so the key reads cleanly ALONGSIDE the "Color by" legend. Own draggable position (recorded on screen).
     */
    private void drawSizeLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        if (!isSizeByProperty()) {
            return;
        }
        final int pad = 7;
        final int gap = 6;        // between the dot column and the value label
        final int row_gap = 4;    // between adjacent value rows
        final int title_gap = 4;  // between the title row and the first value row
        final int max_text = 200; // cap title/label width so a long ref or value can't overrun the figure
        final float base = baseDotSize();
        g.setFont(legendFont());
        final FontMetrics fm = g.getFontMetrics();
        final String title = clipToWidth(
                "Size by: " + PropertyColorScheme.displayName(_property_size_scale.getRef()), fm, max_text);
        final double[] samples = PropertySizeScale.sampleValues(_property_size_scale.getMinValue(),
                _property_size_scale.getMaxValue());
        // one pass: per-sample dot diameter + clipped label, reused for BOTH measuring and drawing (so the box can
        // never be sized for a different dot than the one drawn)
        final int[] diam = new int[samples.length];
        final String[] labels = new String[samples.length];
        int dot_col = 1, label_col = 0;
        for (int i = 0; i < samples.length; ++i) {
            diam[i] = Math.max(1, Math.round(_property_size_scale.diameterForValue(samples[i], base)));
            labels[i] = clipToWidth(PropertySizeScale.formatValue(samples[i]), fm, max_text);
            dot_col = Math.max(dot_col, diam[i]);
            label_col = Math.max(label_col, fm.stringWidth(labels[i]));
        }
        final int content_w = dot_col + gap + label_col;
        final int box_w = Math.max(fm.stringWidth(title), content_w) + (2 * pad);
        int rows_h = 0;
        for (int i = 0; i < diam.length; ++i) {
            rows_h += Math.max(diam[i], fm.getHeight()) + ((i > 0) ? row_gap : 0); // gaps only BETWEEN rows
        }
        final int box_h = (2 * pad) + fm.getHeight() + title_gap + rows_h;
        final Point tl = sizeLegendTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = tl.y;
        // reuse the shared legend chrome (fill/border/title, records _size_legend_bounds when draggable); it leaves
        // STROKE_1 set, so save/restore the branch stroke around it
        final Stroke saved_stroke = g.getStroke();
        drawLegendBox(g, x, y, box_w, box_h, pad, title, fm, draggable, true);
        final Color fg = getTreeColorSet().getSequenceColor();
        // sample dots (smallest at top -> largest at bottom), each labeled with its value
        int row_y = y + pad + fm.getHeight() + title_gap;
        final int dot_center_x = x + pad + (dot_col / 2);
        for (int i = 0; i < samples.length; ++i) {
            final int d = diam[i];
            final int row_h = Math.max(d, fm.getHeight());
            final int cy = row_y + (row_h / 2);
            g.setColor(fg); // neutral fill: the size key is independent of the color-by hue
            g.fillOval(dot_center_x - (d / 2), cy - (d / 2), d, d);
            g.drawString(labels[i], x + pad + dot_col + gap, (cy - (fm.getHeight() / 2)) + fm.getAscent());
            row_y += row_h + row_gap;
        }
        g.setStroke(saved_stroke);
    }

    /** Default position of the size legend: top-right when it is the only legend, else the BOTTOM-right -- clear of
     *  the color/rank legend's top-right default so the two keys never collide by default. Once dragged,
     *  _size_legend_offset maps fractionally like the color legend (so exports honor the moved position). */
    private Point sizeLegendTopLeft(final Rectangle bounds, final int box_w, final int box_h) {
        if (_size_legend_offset != null) {
            return legendTopLeftFor(bounds, getVisibleRect(), _size_legend_offset, box_w, box_h);
        }
        final boolean color_legend_present = isColorByProperty() || hasRankLegend() || hasAnnotationColumnLegend();
        if (!color_legend_present) {
            return legendTopLeftFor(bounds, getVisibleRect(), null, box_w, box_h); // the shared top-right default
        }
        // a color legend already holds the top-right default -> drop to the bottom-right so the two never collide
        return new Point(Math.max(bounds.x, (bounds.x + bounds.width) - box_w - 10),
                Math.max(bounds.y, (bounds.y + bounds.height) - box_h - 10));
    }

    private static String clipToWidth(final String s, final FontMetrics fm, final int max_px) {
        if (fm.stringWidth(s) <= max_px) {
            return s;
        }
        final String ellipsis = "…";
        int len = s.length();
        while ((len > 1) && (fm.stringWidth(s.substring(0, len) + ellipsis) > max_px)) {
            --len;
        }
        return s.substring(0, len) + ellipsis;
    }

    /**
     * Collapse the tree from the given node
     *
     * @param node a PhylogenyNode
     */
    final void collapse(final PhylogenyNode node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            JOptionPane.showMessageDialog(this,
                    "Cannot collapse in unrooted display type",
                    "Attempt to collapse in unrooted display",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }
        if (!node.isExternal() && !node.isRoot()) {
            final boolean collapse = !node.isCollapse();
            TreePanelUtil.collapseSubtree(node, collapse);
            updateSetOfCollapsedExternalNodes();
            _phylogeny.recalculateNumberOfExternalDescendants(true);
            resetNodeIdToDistToLeafMap();
            calculateLongestExtNodeInfo();
            setNodeInPreorderToNull();
            _control_panel.displayedPhylogenyMightHaveChanged(true);
            resetPreferredSize();
            updateOvSizes();
            _main_panel.adjustJScrollPane();
            repaint();
        }
    }

    final void uncollapseAll(final PhylogenyNode node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            JOptionPane.showMessageDialog(this,
                    "Cannot uncollapse in unrooted display type",
                    "Attempt to uncollapse in unrooted display",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }
        if (!node.isExternal()) {
            TreePanelUtil.uncollapseSubtree(node);
            updateSetOfCollapsedExternalNodes();
            _phylogeny.recalculateNumberOfExternalDescendants(true);
            resetNodeIdToDistToLeafMap();
            calculateLongestExtNodeInfo();
            setNodeInPreorderToNull();
            _control_panel.displayedPhylogenyMightHaveChanged(true);
            resetPreferredSize();
            updateOvSizes();
            _main_panel.adjustJScrollPane();
            repaint();
        }
    }


    /**
     * Colorizes subtrees by taxonomic {@code rank} and builds the movable taxon-&gt;color legend
     * (sorted by taxon name). UI-free; returns the number of colorized subtrees. {@link #colorRank}
     * wraps this with the user-facing result dialogs.
     */
    final int colorByRank(final String rank) {
        if ((_phylogeny == null) || (_phylogeny.getNumberOfExternalNodes() < 2)) {
            return 0;
        }
        setWaitCursor();
        // a rank colorization takes over branch coloring, so turn off any active color-by-property
        // (otherwise its legend would be drawn over the rank-colored branches) and sync the dropdown
        setColorByPropertyRef(null);
        _control_panel.populateColorByPropertyBox();
        _rank_legend = null;
        _rank_legend_counts = null;
        _rank_legend_title = null;
        _branch_rank_colorize_rank = null;
        final int colorizations = recolorBranchesByRank(rank); // honors any stored overrides for this rank
        if (colorizations > 0) {
            _control_panel.setColorBranches(true);
            if (_control_panel.getUseVisualStylesCb() != null) {
                _control_panel.getUseVisualStylesCb().setSelected(true);
            }
            _options.setColorLabelsSameAsParentBranch(true);
            if (getMainPanel().getMainFrame()._color_labels_same_as_parent_branch != null) {
                getMainPanel().getMainFrame()._color_labels_same_as_parent_branch.setSelected(true);
            }
            _control_panel.repaint();
        }
        setArrowCursor();
        repaint();
        return colorizations;
    }

    /**
     * Colors the branches by taxon at {@code rank} (honoring this rank's user color overrides) and
     * (re)builds the matching legend. Used both by the menu action ({@link #colorByRank}) and when a
     * legend recolor needs to re-apply. Returns the number of clades colored.
     */
    private int recolorBranchesByRank(final String rank) {
        AptxUtil.removeBranchColors(_phylogeny);
        final Map<String, Color> legend = new HashMap<>();
        final Map<String, Integer> counts = new HashMap<>();
        final int colorizations = TreePanelUtil.colorPhylogenyAccordingToRanks(_phylogeny, rank,
                TreePanelUtil.getDefaultLineageService(), legend, rankOverridesFor(rank), counts);
        if (!legend.isEmpty()) {
            _rank_legend = new java.util.TreeMap<>(legend); // sorted by taxon name
            _rank_legend_counts = counts;
            _rank_legend_title = "Taxonomy: " + rank;
            _branch_rank_colorize_rank = rank;
        }
        return colorizations;
    }

    /**
     * Shows the result dialog after a rank colorization (the colorizing itself is done by
     * {@link #colorByRank}); {@code colorizations} is how many clades were colored. Called on the EDT
     * by {@code MainFrame.colorRank} for the local path and by {@code OnlineTaxonResolver} after
     * a background taxonomy-database resolution.
     */
    final void reportRankColorization(final String rank, final int colorizations) {
        if (colorizations > 0) {
            setEdited(true);
            final String msg = "Taxonomy colorization via " + rank + " completed:\n"
                    + ((colorizations > 1) ? ("colorized " + colorizations + " clades") : "colorized one clade");
            JOptionPane.showMessageDialog(this,
                    msg,
                    "Taxonomy Rank-Colorization Completed (" + rank + ")",
                    JOptionPane.INFORMATION_MESSAGE);
        } else {
            final String msg = "Could not place any tip at rank \"" + rank + "\".\n"
                    + "Try a different rank (e.g. order, family, genus), or check that the tips carry\n"
                    + "taxonomic names the taxonomy database can resolve.";
            JOptionPane
                    .showMessageDialog(this, msg, "Taxonomy Rank-Colorization Failed", JOptionPane.WARNING_MESSAGE);
        }
    }

    // ---- clade-annotation bands (boxes / right-edge bars by taxonomic rank) --------------------

    /** Annotates the tree's clades with boxes or bars at {@code rank}; returns the number of bands. */
    final int setCladeBands(final String rank, final CLADE_VIS mode) {
        _clade_bands_rank = rank;
        _clade_bands_mode = mode;
        rebuildCladeBands();
        repaint();
        return (_clade_bands == null) ? 0 : _clade_bands.size();
    }

    final void clearCladeBands() {
        _clade_bands = null;
        _clade_bands_rank = null;
        // drop the color key with the bands -- unless a branch rank-colorization still owns the legend
        if (_branch_rank_colorize_rank == null) {
            clearRankLegend();
        }
    }

    final boolean hasCladeBands() {
        return (_clade_bands != null) && !_clade_bands.isEmpty();
    }

    CLADE_VIS getCladeBandsMode() {
        return _clade_bands_mode;
    }

    /** Recomputes the bands from the current tree (cache-only); call after navigation swaps the tree. */
    final void rebuildCladeBands() {
        if (ForesterUtil.isEmpty(_clade_bands_rank) || (_phylogeny == null) || _phylogeny.isEmpty()) {
            _clade_bands = null;
            return; // leave any existing (property/rank-colorize) legend untouched
        }
        final Map<String, Integer> counts = new HashMap<>();
        _clade_bands = TreePanelUtil.cladeBands(_phylogeny, _clade_bands_rank,
                TreePanelUtil.getDefaultLineageService(), rankOverridesFor(_clade_bands_rank), counts);
        updateCladeBandLegend(counts);
    }

    /**
     * Mirrors the current clade bands in the shared, draggable taxon-&gt;color legend (with per-taxon tip
     * {@code counts}), so the boxes/bars are labeled even when internal-node data is hidden. Reuses the
     * rank-colorization legend slot (no new control); whichever rank feature was applied last owns the legend.
     */
    private void updateCladeBandLegend(final Map<String, Integer> counts) {
        if (!hasCladeBands()) {
            return;
        }
        if (_clade_bands_mode == CLADE_VIS.BRACKETS) {
            // monochrome brackets carry no color key; drop the clade legend unless a branch
            // rank-colorization still owns the legend slot
            if (_branch_rank_colorize_rank == null) {
                clearRankLegend();
            }
            return;
        }
        final Map<String, Color> legend = new java.util.TreeMap<>(); // sorted by taxon name; dedups repeated taxa
        for (final CladeBand band : _clade_bands) {
            legend.put(band.getTaxon(), band.getColor());
        }
        if (!legend.isEmpty()) {
            _rank_legend = legend;
            _rank_legend_counts = counts;
            _rank_legend_title = "Taxonomy: " + _clade_bands_rank;
        }
    }

    /**
     * Draws the clade annotations over the tree: boxes = translucent color wash; bars = colored right-edge
     * bars with labels (color key in the legend); brackets = monochrome "]" brackets with labels (no color,
     * no legend) for when color already encodes something else.
     */
    // --- tip-aligned annotation columns (color strip / heat map / bar / text) -------------------------------

    final boolean hasAnnotationColumns() {
        return (_annotation_columns != null) && (_annotation_columns.size() > 0);
    }

    java.util.List<AnnotationColumns.ColumnSpec> getAnnotationColumnSpecs() {
        return _annotation_column_specs;
    }

    /** Sets which fields are shown as annotation columns (each spec = a field + render type) and rebuilds. */
    void setAnnotationColumns(final java.util.List<AnnotationColumns.ColumnSpec> specs) {
        _annotation_column_specs = ((specs == null) || specs.isEmpty()) ? null
                : new java.util.ArrayList<AnnotationColumns.ColumnSpec>(specs);
        _focused_annotation_column = -1;
        rebuildAnnotationColumns();
    }

    /** The last annotation import applied to this tree (its source + column mapping), for one-click Re-import, or null. */
    org.forester.archaeopteryx.tools.NodeDataImporter.ImportProfile getLastImportProfile() {
        return _last_import_profile;
    }

    void setLastImportProfile(final org.forester.archaeopteryx.tools.NodeDataImporter.ImportProfile profile) {
        _last_import_profile = profile;
    }

    void clearAnnotationColumns() {
        _annotation_column_specs = null;
        _annotation_columns = null;
        _annotation_col_widths = null;
        _focused_annotation_column = -1;
    }

    /** Rebuilds the column model from the stored specs against the currently displayed tree (visible tips). */
    final void rebuildAnnotationColumns() {
        _annotation_col_widths = null; // invalidate the cached widths -- the model (and its text) may have changed
        if ((_annotation_column_specs == null) || _annotation_column_specs.isEmpty() || (_phylogeny == null)
                || _phylogeny.isEmpty()) {
            _annotation_columns = null;
            _focused_annotation_column = -1;
            return;
        }
        _annotation_columns = new AnnotationColumns(_phylogeny, _annotation_column_specs);
    }

    // --- undo/redo (snapshot-based; see TreeHistory) --------------------------------------------------------

    /**
     * Records the current tree state under {@code label} so the next mutation can be undone. Call this on the
     * EDT immediately BEFORE mutating {@code _phylogeny} (structure, node styles, or node data). Public so the
     * async data tools (Fetch, Infer) in the {@code tools} sub-package can checkpoint at their EDT commit.
     */
    public void pushUndoCheckpoint(final String label) {
        if ((_phylogeny != null) && !_phylogeny.isEmpty()) {
            _history.checkpoint(_phylogeny, label, isEdited());
            notifyEditMenu();
        }
    }

    /**
     * Records an already-captured pre-mutation state as an undo checkpoint — for operations that must run
     * first to know whether they actually changed anything (e.g. import), so they can checkpoint only on a
     * real change instead of leaving a spurious no-op entry. {@code pre_state} must be an independent copy.
     */
    void pushUndoSnapshot(final Phylogeny pre_state, final boolean was_edited, final String label) {
        if ((pre_state != null) && !pre_state.isEmpty()) {
            _history.checkpoint(pre_state, label, was_edited);
            notifyEditMenu();
        }
    }

    boolean canUndo() {
        return _history.canUndo();
    }

    boolean canRedo() {
        return _history.canRedo();
    }

    String undoLabel() {
        return _history.undoLabel();
    }

    String redoLabel() {
        return _history.redoLabel();
    }

    /** Restores the previous tree state; returns true if something was undone. */
    boolean undo() {
        final TreeHistory.Snapshot s = _history.undo(_phylogeny, isEdited());
        if (s == null) {
            return false;
        }
        restoreSnapshot(s);
        return true;
    }

    /** Re-applies the last undone state; returns true if something was redone. */
    boolean redo() {
        final TreeHistory.Snapshot s = _history.redo(_phylogeny, isEdited());
        if (s == null) {
            return false;
        }
        restoreSnapshot(s);
        return true;
    }

    private void restoreSnapshot(final TreeHistory.Snapshot s) {
        _restoring_snapshot = true;
        try {
            final Phylogeny phy = s.getPhylogeny();
            setTree(phy); // also nulls the preorder cache
            setFoundNodes0(null); // the restored copy's search/selection hits from the mutated tree no longer apply
            setFoundNodes1(null);
            phy.externalNodesHaveChanged();
            phy.clearHashIdToNodeMap();
            phy.recalculateNumberOfExternalDescendants(true);
            resetNodeIdToDistToLeafMap();
            updateSetOfCollapsedExternalNodes();
            ensureDomainArchitecturesRenderable(); // Phylogeny.copy() yields plain DomainArchitectures; the
                                                   // layout/paint code requires RenderableDomainArchitecture
            calculateLongestExtNodeInfo();
            recalculateMaxDistanceToRoot();
            resetPreferredSize();
            clearRankLegend(); // the branch rank-colorization legend is UI state, not in the snapshot -- drop it
            rebuildPropertyDisplays(); // color+size schemes summarize the tree -- recompute for the restored one
            rebuildAnnotationColumns();
            rebuildCladeBands();
            setEdited(s.isEdited());
        }
        finally {
            _restoring_snapshot = false;
        }
        notifyEditMenu();
    }

    /**
     * Wraps any plain {@link DomainArchitecture} on an external node's sequence into a
     * {@link RenderableDomainArchitecture}, which the layout ({@link #calculateLongestExtNodeInfo}) and paint
     * code cast to unconditionally. Needed after installing a {@link Phylogeny#copy()} (undo/redo restore),
     * since copy() reproduces the base type, not the renderable subclass. Mirrors the lazy wrapping the paint
     * path already does.
     */
    private void ensureDomainArchitecturesRenderable() {
        if (_phylogeny == null) {
            return;
        }
        for (final PhylogenyNode node : _phylogeny.getExternalNodes()) {
            if (node.getNodeData().isHasSequence()
                    && (node.getNodeData().getSequence().getDomainArchitecture() != null)
                    && !(node.getNodeData().getSequence()
                            .getDomainArchitecture() instanceof RenderableDomainArchitecture)) {
                final RenderableDomainArchitecture rds = SPECIAL_DOMAIN_COLORING
                        ? new RenderableDomainArchitecture(node.getNodeData().getSequence().getDomainArchitecture(),
                                node.getName())
                        : new RenderableDomainArchitecture(node.getNodeData().getSequence().getDomainArchitecture());
                node.getNodeData().getSequence().setDomainArchitecture(rds);
            }
        }
    }

    private void notifyEditMenu() {
        if ((getMainPanel() != null) && (getMainPanel().getMainFrame() != null)) {
            getMainPanel().getMainFrame().updateEditMenu();
        }
    }

    /** Strip / heat map (color key) and bar (length/range key) columns have a legend; text columns do not. */
    private boolean columnHasLegend(final int i) {
        final AnnotationColumns.Type t = _annotation_columns.getColumn(i).getType();
        return (t == AnnotationColumns.Type.COLOR_STRIP) || (t == AnnotationColumns.Type.HEATMAP)
                || (t == AnnotationColumns.Type.MATRIX) || (t == AnnotationColumns.Type.BAR);
    }

    /** Whether the user has explicitly focused (clicked the header of) a legend-bearing annotation column. */
    boolean hasFocusedAnnotationColumn() {
        return (_focused_annotation_column >= 0) && columnLegendReady(_focused_annotation_column);
    }

    /** The annotation column whose legend occupies the shared legend slot: the column the user explicitly focused
     *  (clicked its header) if it has a legend, else -- for the always-on clustergram heat-map legend -- the first
     *  MATRIX column, else -1. So a heat-map matrix shows its shared-scale legend WITHOUT a click, while clicking any
     *  other column's header still overrides to that column's own legend.
     *  <p>The always-on matrix legend DEFERS to an active "Color by" / rank-colorize legend: those share the one
     *  legend slot, and the user explicitly turned them on, so it would be wrong to silently hide them (and their
     *  click-to-recolor) just because a matrix column is present. The matrix legend is still reachable by clicking
     *  its header (an explicit focus, handled above). */
    private int annotationLegendColumn() {
        if (!hasAnnotationColumns()) {
            return -1;
        }
        if ((_focused_annotation_column >= 0) && columnLegendReady(_focused_annotation_column)) {
            return _focused_annotation_column;
        }
        if (isColorByProperty() || hasRankLegend()) {
            return -1; // an active color/rank legend owns the slot; do not usurp it with the always-on matrix legend
        }
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            if ((_annotation_columns.getColumn(i).getType() == AnnotationColumns.Type.MATRIX) && columnLegendReady(i)) {
                return i; // the shared matrix scale -> one always-on legend, drawn from the first matrix column
            }
        }
        return -1;
    }

    /** Whether the shared legend slot should show an annotation-column legend: an explicit focus OR the always-on
     *  matrix legend. */
    boolean hasAnnotationColumnLegend() {
        return annotationLegendColumn() >= 0;
    }

    /** A column has a drawable legend right now: a legend-bearing type whose scheme is present and non-empty (a
     *  column whose visible values vanished after navigating into a subtree has an empty scheme and no legend to
     *  show), and -- for a BAR, whose length key needs a real numeric range -- a gradient rather than a degenerate
     *  single/non-numeric value (else a misleading full-length "0 -> 0" wedge). */
    private boolean columnLegendReady(final int col) {
        if (!hasAnnotationColumns() || (col < 0) || (col >= _annotation_columns.size()) || !columnHasLegend(col)) {
            return false;
        }
        final AnnotationColumns.Column c = _annotation_columns.getColumn(col);
        final PropertyColorScheme s = c.getScheme();
        if ((s == null) || s.isEmpty()) {
            return false;
        }
        if (c.getType() == AnnotationColumns.Type.BAR) {
            return s.isGradient();
        }
        return true;
    }

    /**
     * The annotation column whose header contains ({@code x},{@code y}), or -1. The header band matches where
     * the rotated header is actually drawn (up to {@code anchor_y}), so a click on the visible header toggles
     * its legend even when there is no room above the first cell and the header is clamped down over it.
     */
    /** Device center of an UPRIGHT vertical annotation-column header. Shifted a further half-header-width into the
     *  reserved margin (past the first tip) so the WHOLE header sits before the first cell instead of overlapping it
     *  (the horizontal path does the same by anchoring the rotated header a full text-length up). Shared by the paint,
     *  the hit-test, and the test hook so the drawn header and its click box can't drift. */
    private Point2D.Double verticalHeaderAnchor(final double col_center, final float min_tip_y, final int header_w) {
        return screenPoint(col_center, min_tip_y - getYdistance() - 3.0 - (header_w / 2.0));
    }

    private int annotationHeaderColumnAt(final int x, final int y) {
        if (!hasAnnotationColumns() || (y < 0)) {
            return -1;
        }
        float min_tip_y = Float.MAX_VALUE;
        for (final PhylogenyNode t : visibleExternalTips()) {
            min_tip_y = Math.min(min_tip_y, t.getYcoord());
        }
        if (min_tip_y == Float.MAX_VALUE) {
            return -1;
        }
        final float pad = getYdistance();
        final FontMetrics fm = getFontMetrics(getTreeFontSet().getSmallFont());
        float cx = annotationColumnsStartX();
        if (isVerticalOrientation()) {
            if (_orientation_R == null) {
                return -1; // R is built during paint; before the first vertical paint screenPoint would return logical
            }
            // the header is drawn UPRIGHT in the reserved margin (see verticalHeaderAnchor / paintAnnotationColumnsVertical);
            // hit-test that device box for each column
            final int hh = fm.getHeight();
            for (int i = 0; i < _annotation_columns.size(); ++i) {
                final int w = annotationColumnWidth(i);
                final String header = _annotation_columns.getColumn(i).getHeader();
                if (header.length() > 0) {
                    final int tw = fm.stringWidth(header);
                    final Point2D.Double hp = verticalHeaderAnchor(cx + (w / 2.0), min_tip_y, tw);
                    if ((Math.abs(x - hp.x) <= ((tw / 2.0) + 2)) && (Math.abs(y - hp.y) <= ((hh / 2.0) + 2))) {
                        return i;
                    }
                }
                cx += w + annotationColumnGapAfter(i);
            }
            return -1;
        }
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            final int w = annotationColumnWidth(i);
            if ((x >= cx) && (x <= (cx + w))) {
                // same anchor_y as the drawn header (see paintAnnotationColumns); the header spans y in
                // [anchor_y, anchor_y + textWidth], so a click anywhere from the top down to the header's
                // bottom (just above the first cell) in this column hits it
                final int tw = fm.stringWidth(_annotation_columns.getColumn(i).getHeader());
                final float anchor_y = Math.max(1.0f, (min_tip_y - pad) - 3.0f - tw);
                return (y <= (anchor_y + tw)) ? i : -1;
            }
            cx += w + annotationColumnGapAfter(i);
        }
        return -1;
    }

    /** Toggles the color key (legend) for the annotation column whose header was clicked. */
    final boolean handleAnnotationHeaderClick(final MouseEvent e) {
        final int col = annotationHeaderColumnAt(e.getX(), e.getY());
        if (col < 0) {
            return false;
        }
        setFocusedAnnotationColumn(col);
        return true;
    }

    /** Shows column {@code col}'s legend (toggling it off if it already showed); a text or out-of-range column
     * clears the focus. Strip / heat map (color key) and bar (length/range key) columns have a legend. */
    void setFocusedAnnotationColumn(final int col) {
        if (!hasAnnotationColumns() || (col < 0) || (col >= _annotation_columns.size())
                || !columnHasLegend(col) || (_focused_annotation_column == col)) {
            _focused_annotation_column = -1;
        } else {
            _focused_annotation_column = col;
        }
        repaint();
    }

    /** Test hook: the width (px) reserved for annotation column {@code col}. */
    int annotationColumnWidthForTest(final int col) {
        return annotationColumnWidth(col);
    }

    /** Test hook: the vertical space (px) reserved above the first tip for the rotated column headers. */
    int annotationHeaderTopReserveForTest() {
        return annotationHeaderTopReserve();
    }

    /** A point inside column {@code col}'s clickable header band (for the header-click test), or null. */
    java.awt.Point annotationHeaderMidpointForTest(final int col) {
        if (!hasAnnotationColumns() || (col < 0) || (col >= _annotation_columns.size())) {
            return null;
        }
        float cx = annotationColumnsStartX();
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            final int w = annotationColumnWidth(i);
            if (i == col) {
                if (isVerticalOrientation()) {
                    // the upright header's device center (shared verticalHeaderAnchor -> matches the paint + hit-test)
                    float min_tip_y = Float.MAX_VALUE;
                    for (final PhylogenyNode t : visibleExternalTips()) {
                        min_tip_y = Math.min(min_tip_y, t.getYcoord());
                    }
                    if (min_tip_y == Float.MAX_VALUE) {
                        return null; // no visible tips (all collapsed) -> no header to hit, like the production hit-test
                    }
                    final int tw = getFontMetrics(getTreeFontSet().getSmallFont())
                            .stringWidth(_annotation_columns.getColumn(i).getHeader());
                    final Point2D.Double hp = verticalHeaderAnchor(cx + (w / 2.0), min_tip_y, tw);
                    return new java.awt.Point((int) Math.round(hp.x), (int) Math.round(hp.y));
                }
                return new java.awt.Point(Math.round(cx + (w / 2.0f)), 1); // y=1 is within [0, anchor_y]
            }
            cx += w + annotationColumnGapAfter(i);
        }
        return null;
    }

    private void drawAnnotationColumnLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        final int col_i = annotationLegendColumn();
        if (col_i < 0) {
            return;
        }
        final AnnotationColumns.Column col = _annotation_columns.getColumn(col_i);
        final PropertyColorScheme scheme = col.getScheme();
        if ((scheme == null) || scheme.isEmpty()) {
            return;
        }
        // a MATRIX legend represents the WHOLE shared-scale grid, not one sample column, so title it generically
        // ("Heat map") rather than with the first matrix column's ref (e.g. "s1"); other columns use their header.
        final boolean is_matrix = (col.getType() == AnnotationColumns.Type.MATRIX);
        final String title = is_matrix ? "Heat map" : col.getHeader();
        noteLegendSubject(is_matrix ? "matrix" : ("col:" + col.getHeader()));
        if (col.getType() == AnnotationColumns.Type.BAR) {
            // a bar column encodes magnitude by length, not hue -> show a length wedge with min/max, not a
            // color gradient
            drawAnnotationBarLegend(g, bounds, draggable, title, scheme);
        } else if (scheme.isGradient()) {
            drawPropertyColorGradientLegend(g, bounds, draggable, title, scheme);
        } else {
            final Map<String, Integer> counts = scheme.getValueCounts();
            final Map<String, Color> values = orderedLegend(scheme.getValueColors(), counts);
            drawCategoricalLegend(g, bounds, draggable, title, values, counts,
                    scheme.numberOfValues() - values.size());
        }
    }

    /**
     * Legend for a focused BAR annotation column: a right-growing wedge (short = min value, full = max value)
     * with the numeric range labeled beneath, mirroring how the column itself draws horizontal bars whose
     * length encodes the value. Uses the scheme's gradient min/max as the range endpoints.
     */
    private void drawAnnotationBarLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable,
                                         final String title, final PropertyColorScheme scheme) {
        final int pad = 7;
        final int bar_w = 200;
        final int bar_h = 12;
        g.setFont(legendFont());
        final FontMetrics fm = g.getFontMetrics();
        final String min_lbl = (scheme.getGradientMinLabel() != null) ? scheme.getGradientMinLabel() : "";
        final String max_lbl = (scheme.getGradientMaxLabel() != null) ? scheme.getGradientMaxLabel() : "";
        final int content_w = Math.max(fm.stringWidth(title), bar_w);
        final int box_w = content_w + (2 * pad) + 4;
        final int box_h = (2 * fm.getHeight()) + bar_h + 6 + (2 * pad);
        final Point tl = legendTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = tl.y;
        if (draggable) {
            clearLegendRowControls(); // a bar legend has no clickable value rows / sort / expand controls
        }
        final Stroke saved_stroke = g.getStroke();
        final Color fg = getTreeColorSet().getSequenceColor();
        final int baseline = drawLegendBox(g, x, y, box_w, box_h, pad, title, fm, draggable, false);
        final int bar_x = x + pad;
        final int bar_y = baseline + 4;
        // right-growing wedge: empty (min) at the left, full (max) at the right, echoing the column's bars
        g.fillPolygon(new int[] { bar_x, bar_x + bar_w, bar_x + bar_w },
                new int[] { bar_y + bar_h, bar_y, bar_y + bar_h }, 3);
        g.drawRect(bar_x, bar_y, bar_w - 1, bar_h);
        final int label_baseline = bar_y + bar_h + fm.getAscent() + 2;
        g.drawString(min_lbl, bar_x, label_baseline);
        g.drawString(max_lbl, (bar_x + bar_w) - fm.stringWidth(max_lbl), label_baseline);
        g.setStroke(saved_stroke);
    }

    /** Total horizontal space the annotation columns occupy (0 when none), including the gaps around them.
     * Zero for circular/unrooted layouts, where the columns are not drawn, so no width is reserved. */
    private int annotationColumnsWidth() {
        if (!hasAnnotationColumns() || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            return 0;
        }
        int w = ANN_COL_GAP;
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            w += annotationColumnWidth(i) + annotationColumnGapAfter(i);
        }
        return w;
    }

    /**
     * Vertical space (px) to reserve above the first tip for the rotated annotation-column headers, so they do
     * not overlap the top cells on a tight tree. Because the headers are drawn rotated 90 degrees, the space
     * they need above the tree equals the longest header's horizontal text length. Zero when there are no
     * annotation columns/headers, or for circular/unrooted layouts (which draw no columns).
     */
    private int annotationHeaderTopReserve() {
        if (!hasAnnotationColumns() || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            return 0;
        }
        ensureAnnotationColumnWidths(); // also (re)computes the cached header reserve
        return _annotation_header_top_reserve;
    }

    private int annotationColumnWidth(final int i) {
        ensureAnnotationColumnWidths();
        return _annotation_col_widths[i];
    }

    /** The gap (px) drawn AFTER annotation column {@code i}: 0 between two adjacent heat-map MATRIX columns (so the
     *  matrix reads as one contiguous grid), else the normal {@link #ANN_COL_GAP}. */
    private int annotationColumnGapAfter(final int i) {
        if (((i + 1) < _annotation_columns.size())
                && (_annotation_columns.getColumn(i).getType() == AnnotationColumns.Type.MATRIX)
                && (_annotation_columns.getColumn(i + 1).getType() == AnnotationColumns.Type.MATRIX)) {
            return 0;
        }
        return ANN_COL_GAP;
    }

    /**
     * Populates the per-column width cache (and the rotated-header top reserve) when stale. A TEXT column's
     * width measures every tip's string and the reserve measures every header, so recomputing them on each of
     * the several width queries per repaint is wasteful on large trees; both depend only on the column set and
     * the (zoom-independent) tree font, so they are cached and reused until the model is rebuilt (see
     * {@link #rebuildAnnotationColumns()}) or that font changes (family, size, or style).
     */
    private void ensureAnnotationColumnWidths() {
        final Font large = getFontMetricsForLargeDefaultFont().getFont();
        if ((_annotation_col_widths != null) && (_annotation_col_widths.length == _annotation_columns.size())
                && large.equals(_annotation_col_widths_font)) {
            return;
        }
        final FontMetrics fm = getFontMetricsForLargeDefaultFont();
        final int h = fm.getHeight();
        final int[] widths = new int[_annotation_columns.size()];
        for (int i = 0; i < widths.length; ++i) {
            widths[i] = computeAnnotationColumnWidth(i, fm, h);
        }
        _annotation_col_widths = widths;
        _annotation_col_widths_font = large;
        // the reserve is the longest rotated header (drawn in the small font), computed here so the paint path
        // does not re-measure every header on every repaint (see annotationHeaderTopReserve)
        final FontMetrics sfm = getFontMetrics(getTreeFontSet().getSmallFont());
        int max_header = 0;
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            final String header = _annotation_columns.getColumn(i).getHeader();
            if ((header != null) && (header.length() > 0)) {
                max_header = Math.max(max_header, sfm.stringWidth(header));
            }
        }
        _annotation_header_top_reserve = (max_header > 0) ? (max_header + 6) : 0;
    }

    private int computeAnnotationColumnWidth(final int i, final FontMetrics fm, final int h) {
        switch (_annotation_columns.getColumn(i).getType()) {
            case BAR:
                return Math.max(24, h * 3);
            case TEXT: {
                int max = 0;
                for (final PhylogenyNode t : _phylogeny.getExternalNodes()) {
                    final String v = _annotation_columns.cellText(t, i);
                    if (v.length() > 0) {
                        max = Math.max(max, fm.stringWidth(v));
                    }
                }
                return Math.min(ANN_TEXT_MAX, Math.max(h, max)) + 4;
            }
            default:
                return Math.max(12, h); // color strip / heat map cell
        }
    }

    /** The genuine external nodes currently visible (not hidden under a collapsed clade). */
    private java.util.List<PhylogenyNode> visibleExternalTips() {
        final java.util.List<PhylogenyNode> tips = new java.util.ArrayList<PhylogenyNode>();
        for (final PhylogenyNode t : _phylogeny.getExternalNodes()) {
            if (!isHiddenUnderCollapse(t)) {
                tips.add(t);
            }
        }
        return tips;
    }

    private void paintAnnotationColumns(final Graphics2D g) {
        if (!hasAnnotationColumns()) {
            return;
        }
        final java.util.List<PhylogenyNode> tips = visibleExternalTips();
        final float pad = getYdistance();
        final Font text_font = getTreeFontSet().getSmallFont();
        final Color fg = getTreeColorSet().getSequenceColor();
        final Color saved_color = g.getColor();
        final Font saved_font = g.getFont();
        final AffineTransform saved_tx = g.getTransform();
        float min_tip_y = Float.MAX_VALUE;
        for (final PhylogenyNode t : tips) {
            min_tip_y = Math.min(min_tip_y, t.getYcoord());
        }
        float x = annotationColumnsStartX();
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            final int w = annotationColumnWidth(i);
            final int xi = Math.round(x);
            final AnnotationColumns.Type type = _annotation_columns.getColumn(i).getType();
            for (final PhylogenyNode t : tips) {
                // round each cell's top and bottom to the row boundaries so adjacent cells tile exactly
                // (no 1px seam or overlap) even when getYdistance() is fractional after a fit/zoom
                final int cy = Math.round(t.getYcoord() - pad);
                final int cell_h = Math.max(1, Math.round(t.getYcoord() + pad) - cy);
                switch (type) {
                    case COLOR_STRIP:
                    case HEATMAP:
                    case MATRIX: {
                        final Color c = _annotation_columns.cellColor(t, i);
                        if (c != null) {
                            g.setColor(c);
                            g.fillRect(xi, cy, w, cell_h);
                        }
                        break;
                    }
                    case BAR: {
                        // a present value always draws at least a 1px stub (so the minimum value is visible and
                        // distinct from a missing value, which is NaN and draws nothing)
                        final double f = _annotation_columns.barFraction(t, i);
                        if (!Double.isNaN(f)) {
                            g.setColor(fg);
                            g.fillRect(xi, cy, Math.max(1, (int) Math.round(f * w)), cell_h);
                        }
                        break;
                    }
                    case TEXT: {
                        final String v = _annotation_columns.cellText(t, i);
                        if (v.length() > 0) {
                            g.setFont(text_font);
                            g.setColor(fg);
                            final FontMetrics fm = g.getFontMetrics();
                            g.drawString(fitText(v, w, fm), xi,
                                    Math.round(t.getYcoord() + (fm.getAscent() / 2.0f) - 1));
                        }
                        break;
                    }
                }
            }
            // column header: the field name, rotated to read bottom-to-top. The rotated text extends DOWN from
            // anchor_y (to anchor_y + textWidth), so anchor it a full text-length above the first cell to sit in
            // the reserved band there (see annotationHeaderTopReserve) instead of overlapping the top cells.
            final String header = _annotation_columns.getColumn(i).getHeader();
            if ((header.length() > 0) && (min_tip_y < Float.MAX_VALUE)) {
                g.setFont(text_font);
                g.setColor(fg);
                final FontMetrics hfm = g.getFontMetrics();
                final int tw = hfm.stringWidth(header);
                final float cx = xi + (w / 2.0f);
                final float anchor_y = Math.max(1.0f, (min_tip_y - pad) - 3.0f - tw); // keep the header on-canvas
                g.rotate(-Math.PI / 2.0, cx, anchor_y);
                g.drawString(header, cx - tw, anchor_y + (hfm.getAscent() / 2.0f));
                g.setTransform(saved_tx);
            }
            x += w + annotationColumnGapAfter(i);
        }
        g.setColor(saved_color);
        g.setFont(saved_font);
        g.setTransform(saved_tx);
    }

    /**
     * Vertical-orientation twin of {@link #paintAnnotationColumns}: the tip-aligned columns become horizontal BANDS
     * below the tips. Strip / heat-map / bar CELLS are axis-aligned rectangles, so they ride the rotation R for free
     * (a 90-degree rotation keeps them axis-aligned) -- a logical column (a vertical strip) rotates into a horizontal
     * band. TEXT cells and the column HEADERS are re-anchored to the upright base frame (drawing them under R would
     * render them sideways): a TEXT cell always reads VERTICALLY (90deg) along its band (independent of the "Tip label
     * angle" setting -- a 45deg cell would overrun the fixed band depth), and the header sits upright in the reserved
     * margin just before the first cell. Called while g is still rotated by R (paintPhylogeny).
     */
    private void paintAnnotationColumnsVertical(final Graphics2D g) {
        if (!hasAnnotationColumns()) {
            return;
        }
        final java.util.List<PhylogenyNode> tips = visibleExternalTips();
        if (tips.isEmpty()) {
            return;
        }
        final float pad = getYdistance();
        final Color fg = getTreeColorSet().getSequenceColor();
        final Color saved_color = g.getColor();
        final Font saved_font = g.getFont();
        final AffineTransform withR = g.getTransform();
        final Font text_font = getTreeFontSet().getSmallFont();
        final FontMetrics text_fm = getFontMetrics(text_font);
        // TEXT cells always read VERTICALLY (90deg); the sign follows the orientation so the text is never upside-down
        final double text_angle = (getOptions().getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_BOTTOM)
                ? (-Math.PI / 2.0) : (Math.PI / 2.0);
        float min_tip_y = Float.MAX_VALUE;
        for (final PhylogenyNode t : tips) {
            min_tip_y = Math.min(min_tip_y, t.getYcoord());
        }
        float x = annotationColumnsStartX();
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            final int w = annotationColumnWidth(i);
            final int xi = Math.round(x);
            final float col_center = x + (w / 2.0f);
            final AnnotationColumns.Type type = _annotation_columns.getColumn(i).getType();
            for (final PhylogenyNode t : tips) {
                // round each cell's near/far edge to the row boundaries so adjacent cells tile exactly (as the
                // horizontal path does), then draw the rect in LOGICAL coords -- it rides R into a band cell
                final int cy = Math.round(t.getYcoord() - pad);
                final int cell_h = Math.max(1, Math.round(t.getYcoord() + pad) - cy);
                switch (type) {
                    case COLOR_STRIP:
                    case HEATMAP:
                    case MATRIX: {
                        final Color c = _annotation_columns.cellColor(t, i);
                        if (c != null) {
                            g.setColor(c);
                            g.fillRect(xi, cy, w, cell_h);
                        }
                        break;
                    }
                    case BAR: {
                        final double f = _annotation_columns.barFraction(t, i);
                        if (!Double.isNaN(f)) {
                            g.setColor(fg);
                            g.fillRect(xi, cy, Math.max(1, (int) Math.round(f * w)), cell_h);
                        }
                        break;
                    }
                    case TEXT: {
                        final String v = _annotation_columns.cellText(t, i);
                        if (v.length() > 0) {
                            final String fitted = fitText(v, w, text_fm);
                            final float ty = t.getYcoord();
                            final float tx = x;
                            withNodeTextFrame(g, tx, ty, text_angle, () -> {
                                g.setFont(text_font);
                                g.setColor(fg);
                                g.drawString(fitted, tx, ty + (text_fm.getAscent() / 2.0f) - 1);
                            });
                        }
                        break;
                    }
                }
            }
            // column header: upright, in the reserved margin just BEFORE the first cell (verticalHeaderAnchor shifts
            // it a half-width past the tip so it does not overlap the top cell). min_tip_y is always set (empty tips
            // return early at the top of the method).
            final String header = _annotation_columns.getColumn(i).getHeader();
            if (header.length() > 0) {
                g.setFont(text_font);
                final FontMetrics hfm = g.getFontMetrics();
                final int hw = hfm.stringWidth(header);
                final Point2D.Double hp = verticalHeaderAnchor(col_center, min_tip_y, hw);
                g.setTransform(_orientation_base_transform);
                g.setColor(fg);
                g.drawString(header, (float) (hp.x - (hw / 2.0)), (float) (hp.y + (hfm.getAscent() / 2.0f)));
                g.setTransform(withR);
            }
            x += w + annotationColumnGapAfter(i);
        }
        g.setTransform(withR);
        g.setColor(saved_color);
        g.setFont(saved_font);
    }

    private static String fitText(final String s, final int max_w, final FontMetrics fm) {
        if (fm.stringWidth(s) <= max_w) {
            return s;
        }
        final int ell = fm.stringWidth("…");
        String t = s;
        while ((t.length() > 1) && ((fm.stringWidth(t) + ell) > max_w)) {
            t = t.substring(0, t.length() - 1);
        }
        return t + "…";
    }

    /**
     * "Pulse Found Nodes": a translucent found-color halo disc around each found/selected node that is actually
     * drawn (skipping hits hidden under a collapse). On SCREEN the disc's radius + alpha breathe over
     * PULSE_PERIOD_MS -- an EDT timer ({@link #updatePulseTimer}) repaints just the small halo regions (a CLIPPED
     * repaint: the node walk still runs, but only those regions rasterize). In an EXPORT it renders once at peak
     * phase (a static glow) so the figure still shows the emphasis; suppressed on a black-and-white export (it is
     * inherently colored). Drawn after the node loop (needs the node coords). Rectangular-family layouts only.
     * <p>The screen animation state ({@code _found_halo_bounds} / {@code _has_visible_found_halo}) is cleared at
     * the TOP of {@code paintPhylogeny} and the timer reconciled at its END for EVERY layout -- so an export paint
     * here never clobbers it, and a switch to a circular/unrooted view reliably STOPS the timer (no leak).
     */
    private void paintFoundNodeHalos(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        final boolean to_screen = !to_pdf && !to_graphics_file;
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite();
        if (!getOptions().isPulseFoundNodes() || (_phylogeny == null) || bw || !hasFoundNodes()) {
            return;
        }
        // phase in [0,1]: time-driven breathing on screen (monotonic clock), fixed peak (a static glow) in an export
        final double sin = to_screen
                ? (0.5 + (0.5 * Math.sin((TWO_PI * (System.nanoTime() % PULSE_PERIOD_NS)) / PULSE_PERIOD_NS)))
                : 1.0;
        final float radius = HALO_BASE_RADIUS + (HALO_AMP_RADIUS * (float) sin);
        final int alpha = (int) Math.round(HALO_MIN_ALPHA + ((HALO_MAX_ALPHA - HALO_MIN_ALPHA) * sin));
        final int max_r = Math.round(HALO_BASE_RADIUS + HALO_AMP_RADIUS); // repaint region uses the peak radius
        final Color saved = g.getColor();
        for (final PhylogenyNode node : _nodes_in_preorder) {
            if (!isInFoundNodes(node) || isHiddenUnderCollapse(node)) {
                continue;
            }
            final Color fc = getColorForFoundNode(node);
            g.setColor(new Color(fc.getRed(), fc.getGreen(), fc.getBlue(), alpha));
            final int x = Math.round(node.getXcoord()), y = Math.round(node.getYcoord());
            final int r = Math.round(radius);
            g.fillOval(x - r, y - r, 2 * r, 2 * r);
            if (to_screen) {
                _has_visible_found_halo = true; // screen-only: an export must not touch the animation state
                // the halo is drawn at logical (x,y) but g is rotated by R, so its on-screen centre is R(x,y); record
                // the DEVICE region (the pulse timer repaints these rectangles in device space, not logical)
                final Point2D.Double dev = screenPointFor(node);
                final int dx = Math.round((float) dev.x), dy = Math.round((float) dev.y);
                _found_halo_bounds.add(new Rectangle(dx - max_r - 1, dy - max_r - 1, (2 * max_r) + 2, (2 * max_r) + 2));
            }
        }
        g.setColor(saved);
    }

    /** Starts the pulse timer when there is a visible hit halo to animate on a showing panel, stops it otherwise.
     *  Called at the END of every screen paint (all layouts). The tick also self-terminates, so a hidden tab
     *  (which stops repainting) can't keep it running; and {@link #removeNotify()} stops it on detach. */
    private void updatePulseTimer() {
        if (shouldPulse()) {
            if (_pulse_timer == null) {
                _pulse_timer = new Timer(PULSE_INTERVAL_MS, e -> {
                    if (shouldPulse()) {
                        for (final Rectangle rect : _found_halo_bounds) {
                            repaint(rect.x, rect.y, rect.width, rect.height);
                        }
                    } else {
                        _pulse_timer.stop();
                    }
                });
            }
            if (!_pulse_timer.isRunning()) {
                _pulse_timer.start();
            }
        } else if ((_pulse_timer != null) && _pulse_timer.isRunning()) {
            _pulse_timer.stop();
        }
    }

    /** Whether the found-node pulse should currently be animating: the option is on, at least one hit halo was
     *  drawn on the last screen paint, and the panel is showing. */
    private boolean shouldPulse() {
        return getOptions().isPulseFoundNodes() && _has_visible_found_halo && isShowing();
    }

    @Override
    public void removeNotify() {
        if (_pulse_timer != null) {
            _pulse_timer.stop(); // stop animating once the tab/window is gone (no lingering repaints)
        }
        super.removeNotify();
    }

    /** Test hook: whether the found-node pulse animation timer is currently running. */
    boolean isPulseTimerRunning() {
        return (_pulse_timer != null) && _pulse_timer.isRunning();
    }

    /** Test hook: number of found-node halo repaint regions recorded by the last screen paint (0 = nothing to
     *  pulse, so the timer will stop -- e.g. after switching to a circular/unrooted layout). */
    int getFoundHaloBoundsCountForTest() {
        return _found_halo_bounds.size();
    }

    /** Test hook: does any recorded halo repaint region contain the DEVICE point (x,y)? The pulse timer repaints
     *  these regions in device space, so in a vertical orientation they must be at the node's on-screen position. */
    boolean foundHaloBoundContainsForTest(final int x, final int y) {
        for (final Rectangle r : _found_halo_bounds) {
            if (r.contains(x, y)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Zebra striping: a faint alternating background band behind every other visible tip row, spanning the full
     * width, so a label is easy to track across a wide tree to its annotation columns (the iTOL row-shading aid).
     * Theme-aware and translucent (branches/labels show through); drawn after the node loop (coords are set there).
     * Rectangular layouts only (this is called from the rectangular paint branch). Collapsed-clade triangles count
     * as one row; nodes hidden under a collapse are skipped so the alternation matches the DRAWN rows.
     */
    private void paintZebraStripes(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file,
                                  final int graphics_file_x, final int graphics_file_width) {
        if (!getOptions().isShowZebraStripes() || (_phylogeny == null) || _export_transparent_background) {
            return; // suppressed on a transparent-PNG export: full-width bands would defeat the clean cut-out
        }
        final float row_h = 2f * getYdistance();
        if (!(row_h > 0f)) { // also rejects a NaN pitch (NaN <= 0f is false), not just zero/negative
            return;
        }
        // the stripe spans the full CROSS-tree extent (perpendicular to the tip axis): the device width in a horizontal
        // layout, but the device HEIGHT (the depth axis) in a vertical orientation, where the band rides R -> a vertical
        // stripe over each alternate tip's column. left=0 in vertical (the band starts at the root/logical origin).
        final boolean use_export = (to_pdf || to_graphics_file) && (graphics_file_width > 0);
        final int left = isVerticalOrientation() ? 0 : (use_export ? graphics_file_x : 0);
        final int width = isVerticalOrientation() ? getHeight() : (use_export ? graphics_file_width : getWidth());
        // the color scheme is the single source of truth for light/dark (ExportTheme flips it to LIGHT for a
        // white-background export, so this is correct on every render path -- screen, WYSIWYG, and white-bg export)
        final boolean dark = getTreeColorSet().getCurrentColorScheme() == TreeColorSet.DARK_COLOR_SCHEME;
        final Color saved = g.getColor();
        g.setColor(dark ? ZEBRA_STRIPE_ON_DARK : ZEBRA_STRIPE_ON_LIGHT); // faint translucent, theme-aware
        int row = 0;
        for (final PhylogenyNode node : _nodes_in_preorder) {
            // only the drawn leaf-level rows (tips + collapsed-clade triangles), top to bottom; the cheap shape test
            // first so the O(depth) hidden-under-collapse walk is skipped for the internal-node majority
            if ((!node.isExternal() && !node.isCollapse()) || isHiddenUnderCollapse(node)) {
                continue;
            }
            if ((row % 2) == 1) {
                g.fillRect(left, Math.round(node.getYcoord() - (row_h / 2f)), width, Math.round(row_h));
            }
            row++;
        }
        g.setColor(saved);
    }

    /**
     * Node age (HPD) bars: on a dated phylogram, a translucent horizontal bar at each internal node spanning its age
     * uncertainty -- the FigTree "node bars" standard, showing divergence-time uncertainty. Reads the phyloXML native
     * {@code <date>} (value/min/max, the same model shown in the node popup and editor); the bar is anchored to the
     * node's OWN drawn x plus signed age deltas, so it stays centred on the node even if the tree is not strictly
     * ultrametric or the root has a branch length. Drawn after the node loop (coords are set there), translucent so
     * the branches show through; skips nodes hidden under a collapse. Phylograms only. Assumes the age unit matches
     * the branch-length (time) unit -- the dated-tree use.
     */
    private void paintHpdBars(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        if (!getOptions().isShowHpdBars() || !getControlPanel().isDrawPhylogram() || (_phylogeny == null)) {
            return;
        }
        final double corr = getXcorrectionFactor();
        if (corr <= 0) {
            return; // no branch-length scale -> nothing meaningful to place
        }
        final Color saved = g.getColor();
        g.setColor(((to_pdf || to_graphics_file) && getOptions().isPrintBlackAndWhite()) ? HPD_BAR_COLOR_BW
                : HPD_BAR_COLOR);
        for (final PhylogenyNode node : _nodes_in_preorder) {
            if (node.isExternal() || isHiddenUnderCollapse(node) || !node.getNodeData().isHasDate()) {
                continue;
            }
            final org.forester.phylogeny.data.Date date = node.getNodeData().getDate();
            if ((date.getMin() == null) || (date.getMax() == null)) {
                continue; // need an interval to draw a bar
            }
            final double min = date.getMin().doubleValue();
            final double max = date.getMax().doubleValue();
            final double value = (date.getValue() != null) ? date.getValue().doubleValue() : ((min + max) / 2.0);
            final float[] xr = TreePanelUtil.hpdBarXRange(node.getXcoord(), value, min, max, corr);
            final int left = Math.round(Math.min(xr[0], xr[1])); // robust to swapped/degenerate bounds
            final int w = Math.max(1, Math.round(Math.abs(xr[1] - xr[0])));
            g.fillRect(left, Math.round(node.getYcoord()) - (HPD_BAR_HEIGHT / 2), w, HPD_BAR_HEIGHT);
        }
        g.setColor(saved);
    }

    private void paintCladeBands(final Graphics2D g) {
        if (!hasCladeBands()) {
            return;
        }
        switch (_clade_bands_mode) {
            case BARS:
                drawCladeBars(g);
                break;
            case BRACKETS:
                drawCladeBrackets(g);
                break;
            default:
                drawCladeBoxes(g);
                break;
        }
    }

    /**
     * The x just past the end of the longest currently shown tip label (all text fields plus any aligned
     * domain/vector/MSA columns, via {@link #getLongestExtNodeInfo()}). Mirrors the tree's own label-end
     * logic: in a phylogram labels start at each tip's own x, so the rightmost label ends past the DEEPEST
     * tip; in an aligned/cladogram view all tips share one x. When external labels are hidden ("Show External
     * Data" off) they occupy no width, so annotations sit right after the tips.
     */
    private float labelsRightEdge() {
        // in a vertical orientation the tip labels are TILTED and extend along the DEPTH axis only by their depth
        // component (full width for 90deg, L/sqrt(2) for 45deg), so the columns clear the labels using that reach --
        // not the full horizontal label width, which would leave a large empty gap below the tips
        final float label_w = getControlPanel().isShowExternalData()
                ? (isVerticalOrientation() ? depthLabelReserve() : getLongestExtNodeInfo()) : 0;
        return tipsDepthEdge() + label_w;
    }

    /** The depth (logical x) where the tip branches END -- the anchor for tip labels and tip-aligned columns, WITHOUT
     *  the tip-label reach. (labelsRightEdge() adds the label reach; the clustergram "labels below columns" layout puts
     *  the columns here, right past the tips, and the labels past the columns instead.) */
    private float tipsDepthEdge() {
        if (getControlPanel().isDrawPhylogram()) {
            return (float) ((getMaxDistanceToRoot() * getXcorrectionFactor()) + _phylogeny.getRoot().getXcoord());
        }
        return getPhylogeny().getFirstExternalNode().getXcoord();
    }

    /** Clustergram layout: in a vertical orientation with annotation columns, the tip labels are drawn BELOW the
     *  columns (so the dendrogram sits directly on the tip-aligned grid) rather than between the tree and the columns.
     *  No-op in the horizontal orientation / without annotation columns. */
    boolean tipLabelsBelowColumns() {
        return getOptions().isTipLabelsBelowColumns() && isVerticalOrientation() && hasAnnotationColumns();
    }

    /** The depth (logical x) where the tip-aligned annotation columns START: right past the tips in the clustergram
     *  "labels below" layout, else past the tip labels. */
    private float annotationColumnsStartX() {
        return (tipLabelsBelowColumns() ? tipsDepthEdge() : labelsRightEdge()) + ANN_COL_GAP;
    }

    /** The depth (logical x) where the annotation columns END (one ANN_COL_GAP past the last column). Derived from
     *  {@link #annotationColumnsStartX()} so the two share the one base-selection and can never disagree. */
    private float annotationColumnsEndX() {
        return (annotationColumnsStartX() - ANN_COL_GAP) + annotationColumnsWidth();
    }

    /** The x where clade bands/bars start: past the tip labels and any tip-aligned annotation columns. */
    private float cladeBandRightEdge() {
        return labelsRightEdge() + annotationColumnsWidth() + CLADE_BAND_RIGHT_PAD;
    }

    /**
     * Horizontal space needed to the RIGHT of the label reservation ({@link #getLongestExtNodeInfo()}) so
     * layout / "fit to window" / zoom keep it all on-screen: the tip-aligned annotation columns plus, beyond
     * them, any clade marks (boxes add the small pad; bars/brackets add the gap, the mark, and the rotated
     * label whose horizontal extent is the font height). 0 when neither columns nor bands are shown.
     */
    private int rightMarginExtraWidth() {
        int extra = annotationColumnsWidth();
        if (hasCladeBands()) {
            if (_clade_bands_mode == CLADE_VIS.BOXES) {
                extra += CLADE_BAND_RIGHT_PAD;
            } else {
                final int label_h = getFontMetricsForLargeDefaultFont().getHeight();
                final int mark = (_clade_bands_mode == CLADE_VIS.BARS) ? (CLADE_BAR_WIDTH + 3) : 4;
                extra += CLADE_BAND_RIGHT_PAD + CLADE_BAR_GAP + mark + label_h + 4;
            }
        }
        return extra;
    }

    /** {@code {yTop, yBottom}} of a clade's tips in current paint coordinates, or null if none. */
    private float[] cladeBandYRange(final PhylogenyNode root) {
        float min = Float.MAX_VALUE;
        float max = -Float.MAX_VALUE;
        if (root.isExternal()) {
            min = max = root.getYcoord();
        } else {
            for (final PhylogenyNode t : root.getAllExternalDescendants()) {
                final float y = t.getYcoord();
                if (y < min) {
                    min = y;
                }
                if (y > max) {
                    max = y;
                }
            }
        }
        return (min > max) ? null : new float[] { min, max };
    }

    private void drawCladeBoxes(final Graphics2D g) {
        final float right = cladeBandRightEdge();
        final float pad = getYdistance();
        for (final CladeBand band : _clade_bands) {
            final float[] yr = cladeBandYRange(band.getRoot());
            if (yr == null) {
                continue;
            }
            final float left = band.getRoot().getXcoord();
            final int w = Math.round(right - left);
            final int h = Math.round((yr[1] - yr[0]) + (2 * pad));
            if ((w <= 0) || (h <= 0)) {
                continue;
            }
            final Color c = band.getColor();
            g.setColor(new Color(c.getRed(), c.getGreen(), c.getBlue(), CLADE_BOX_ALPHA));
            g.fillRect(Math.round(left), Math.round(yr[0] - pad), w, h);
        }
    }

    private void drawCladeBars(final Graphics2D g) {
        final float bar_x = cladeBandRightEdge() + CLADE_BAR_GAP;
        final float pad = getYdistance();
        final Font font = getTreeFontSet().getLargeFont();
        for (final CladeBand band : _clade_bands) {
            final float[] yr = cladeBandYRange(band.getRoot());
            if (yr == null) {
                continue;
            }
            final int y = Math.round(yr[0] - pad);
            final int h = Math.max(1, Math.round((yr[1] - yr[0]) + (2 * pad)));
            g.setColor(band.getColor());
            g.fillRect(Math.round(bar_x), y, CLADE_BAR_WIDTH, h); // rides R -> a horizontal bar in a vertical orientation
            // taxon label: rotated 90deg (reading bottom-to-top) beside the vertical bar in a horizontal layout;
            // UPRIGHT (re-anchored to the base frame) beside the horizontal bar in a vertical orientation
            g.setFont(font);
            g.setColor(getTreeColorSet().getSequenceColor());
            final FontMetrics fm = g.getFontMetrics();
            final float label_x = bar_x + CLADE_BAR_WIDTH + 3;
            final float mid_y = (yr[0] + yr[1]) / 2.0f;
            drawCladeBandLabel(g, band.getTaxon(), label_x, mid_y, fm);
        }
    }

    // Monochrome "named clade" annotation: a thin "]" bracket (vertical spine + short end-ticks pointing
    // toward the tips) per clade with the rotated taxon label, all in the foreground color -- no per-clade
    // colors and no legend, for use when color already encodes something else.
    private void drawCladeBrackets(final Graphics2D g) {
        final float spine_x = cladeBandRightEdge() + CLADE_BAR_GAP;
        // less outward pad than the filled bars/boxes, so adjacent brackets keep a clear vertical gap
        // (one "]" per clade) instead of merging into a single continuous line
        final float pad = getYdistance() * 0.3f;
        final Font font = getTreeFontSet().getLargeFont();
        final Stroke saved_stroke = g.getStroke();
        g.setStroke(new BasicStroke(CLADE_BRACKET_STROKE));
        g.setColor(getTreeColorSet().getSequenceColor());
        for (final CladeBand band : _clade_bands) {
            final float[] yr = cladeBandYRange(band.getRoot());
            if (yr == null) {
                continue;
            }
            final int x = Math.round(spine_x);
            final int y0 = Math.round(yr[0] - pad);
            final int y1 = Math.round(yr[1] + pad);
            // "]" : vertical spine plus short end-ticks pointing left, toward the tips
            g.drawLine(x, y0, x, y1);
            g.drawLine(x, y0, x - CLADE_BRACKET_TICK, y0);
            g.drawLine(x, y1, x - CLADE_BRACKET_TICK, y1);
            // taxon label: rotated beside the vertical bracket (horizontal layout), UPRIGHT beside the horizontal
            // bracket (vertical orientation) -- see drawCladeBandLabel
            g.setFont(font);
            final FontMetrics fm = g.getFontMetrics();
            final float label_x = spine_x + 4;
            final float mid_y = (yr[0] + yr[1]) / 2.0f;
            drawCladeBandLabel(g, band.getTaxon(), label_x, mid_y, fm);
        }
        g.setStroke(saved_stroke);
    }

    /** Draws a clade bar/bracket taxon label: rotated 90deg (reading bottom-to-top) beside the vertical mark in the
     *  horizontal layout, or UPRIGHT (re-anchored to the base frame, centered on the clade's breadth) beside the
     *  horizontal mark in a vertical orientation. Restores g's transform (R when vertical) before returning. */
    private void drawCladeBandLabel(final Graphics2D g, final String taxon, final float label_x, final float mid_y,
                                    final FontMetrics fm) {
        final int tw = fm.stringWidth(taxon);
        final AffineTransform saved = g.getTransform();
        if (isVerticalOrientation()) {
            final Point2D.Double lp = screenPoint(label_x, mid_y);
            g.setTransform(_orientation_base_transform);
            g.drawString(taxon, (float) (lp.x - (tw / 2.0)), (float) (lp.y + (fm.getAscent() / 2.0f)));
        }
        else {
            g.rotate(-Math.PI / 2.0, label_x, mid_y);
            g.drawString(taxon, label_x - (tw / 2.0f), mid_y + fm.getAscent());
        }
        g.setTransform(saved);
    }

    final void decreaseDomainStructureEvalueThresholdExp() {
        if (_domain_structure_e_value_thr_exp > -20) {
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
    final PhylogenyNode findNode(final int x, final int y) {
        if ((_phylogeny == null) || _phylogeny.isEmpty()) {
            return null;
        }
        final int half_box_size_plus_wiggle = (getOptions().getDefaultNodeShapeSize() / 2) + WIGGLE;
        // in a vertical orientation the node coords are logical (un-rotated); map the device click back to that space
        final Point2D.Double p = toLogicalPoint(x, y);
        for (final PhylogenyNodeIterator iter = _phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ((_phylogeny.isRooted() || !node.isRoot() || (node.getNumberOfDescendants() > 2))
                    && ((node.getXcoord() - half_box_size_plus_wiggle) <= p.x)
                    && ((node.getXcoord() + half_box_size_plus_wiggle) >= p.x)
                    && ((node.getYcoord() - half_box_size_plus_wiggle) <= p.y)
                    && ((node.getYcoord() + half_box_size_plus_wiggle) >= p.y)
                    // skip nodes hidden under a collapsed ancestor: they keep stale pre-collapse coords, so a
                    // box hit there is a phantom. The collapsed clade-root itself is still returned (it is drawn
                    // as the triangle and IS selectable). Gated behind the cheap box test above (O(depth), rare).
                    && !isHiddenUnderCollapse(node)) {
                return node;
            }
        }
        return null;
    }

    /**
     * The child node whose incoming branch is under ({@code x},{@code y}), or null -- for branch-click clade
     * selection. Picks the nearest branch within a small tolerance (same forgiveness as {@link #findNode}).
     * The rectangular family hit-tests the horizontal leg (at the child's y); the diagonal styles
     * (triangular/convex/curved) hit-test the straight parent-&gt;child line. No-op for circular/unrooted
     * (radial branches). Only branches actually on screen (not hidden under a collapsed ancestor) are tested.
     */
    final PhylogenyNode findBranch(final int x, final int y) {
        final PHYLOGENY_GRAPHICS_TYPE gt = getPhylogenyGraphicsType();
        // Only the layouts whose branch geometry reconstructs EXACTLY from the node coords are supported, so
        // the hit region matches what is drawn: RECTANGULAR (horizontal leg at the child's y + vertical fork
        // connector at the parent's x) and TRIANGULAR (a straight parent->child line). EURO_STYLE/ROUNDED offset
        // the leg near the parent, so branch-click is a no-op there (and no hover cursor is shown) rather than a
        // hit region that disagrees with the painted branch. Ditto the radial layouts (circular/unrooted).
        final boolean diagonal = (gt == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR);
        final boolean rectangular = (gt == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
        if ((_phylogeny == null) || _phylogeny.isEmpty() || !(rectangular || diagonal)) {
            return null;
        }
        final double tol = (getOptions().getDefaultNodeShapeSize() / 2.0) + WIGGLE;
        // in a vertical orientation the node coords are logical (un-rotated); map the device click back to that space
        final Point2D.Double click = toLogicalPoint(x, y);
        PhylogenyNode best = null;
        double best_dist = tol;
        for (final PhylogenyNodeIterator iter = _phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if (n.isCollapse()) {
                continue; // a collapsed clade's tips are hidden (not individually selectable)
            }
            double d = Double.MAX_VALUE;
            if (!n.isRoot()) { // the incoming branch: horizontal leg (at the child's y) or straight diagonal
                final PhylogenyNode p = n.getParent();
                final double a_y = diagonal ? p.getYcoord() : n.getYcoord();
                d = Line2D.ptSegDist(p.getXcoord(), a_y, n.getXcoord(), n.getYcoord(), click.x, click.y);
            }
            if (rectangular && !n.isExternal()) {
                // n's own vertical fork connector (at n's x, spanning its children) selects n's subtree -- this
                // is also how the root becomes selectable (click its fork to select the whole tree's tips)
                d = Math.min(d, Line2D.ptSegDist(n.getXcoord(), n.getFirstChildNode().getYcoord(), n.getXcoord(),
                        n.getLastChildNode().getYcoord(), click.x, click.y));
            }
            // the ancestor-collapse walk is O(depth), so gate it on the cheap distance test (near candidates only)
            if ((d < best_dist) && !isHiddenUnderCollapse(n)) {
                best_dist = d;
                best = n;
            }
        }
        return best;
    }

    /** True if a strict ancestor of {@code n} is collapsed (so {@code n} is not currently drawn). */
    private static boolean isHiddenUnderCollapse(final PhylogenyNode n) {
        PhylogenyNode p = n.getParent();
        while (p != null) {
            if (p.isCollapse()) {
                return true;
            }
            p = p.getParent();
        }
        return false;
    }

    // Branch-hover selection preview (on screen only): the translucent "will be added" alpha (applied to the
    // found-selection color) and the distinct "will be removed" mark (a muted grey, so it reads as
    // "de-selected" and is not confused with the coloured selection whatever the found colour is set to).
    private static final int   BRANCH_HOVER_ALPHA  = 95;
    private static final Color BRANCH_HOVER_REMOVE = new Color(90, 90, 90, 135); // click would DESELECT this clade

    /** Sets the hover-preview target -- a single node ({@code subtree}=false, drawn as one marker on the node)
     *  or a branch ({@code subtree}=true, drawn on the subtree's tips) -- repainting only on an actual change
     *  (so gliding along one branch/over one node does not repaint). */
    private void setHover(final PhylogenyNode node, final boolean subtree) {
        final boolean sub = (node != null) && subtree; // canonical: no subtree flag when nothing is hovered
        if ((node != _hover_node) || (sub != _hover_subtree)) {
            _hover_node = node;
            _hover_subtree = sub;
            repaint();
        }
    }

    /** Applies a hover target, honoring the just-clicked suppression: the target just clicked stays un-previewed
     *  until the pointer moves off it (so it doesn't instantly flip to the "will be removed" grey). */
    private void applyHover(final PhylogenyNode node, final boolean subtree) {
        if (node == _click_suppressed) {
            setHover(null, false);
        } else {
            _click_suppressed = null;
            setHover(node, subtree);
        }
    }

    /** After a click selected {@code target} (a clicked node or branch), suppress its hover preview until the
     *  pointer moves off it, so it doesn't instantly flip to the "will be removed" grey over what was just
     *  selected. */
    private void suppressHoverAfterClick(final PhylogenyNode target) {
        _click_suppressed = target;
        setHover(null, false);
    }

    /** Clears the hover preview (e.g. when the pointer leaves the panel or moves onto the legend/overview). Does
     *  NOT clear the just-clicked suppression -- that is released when the pointer moves onto a different target
     *  (see applyHover) -- so a spurious MOUSE_EXITED (e.g. from a node-data popup) can't re-arm the grey flip. */
    final void clearHoverPreview() {
        setHover(null, false);
    }

    /**
     * On-screen preview of what a click would select/deselect: draws a translucent marker on the affected
     * node(s) so the pending click's effect is visible before committing. For a hovered branch the marks land
     * on the subtree's tips; for a hovered single node (leaf or internal) it's one mark on that node. Direction
     * aware -- targets that would be ADDED get the found-selection colour; when a click would DESELECT (the
     * node, or a fully-selected clade) they get a distinct, larger "remove" grey. No-op off screen (exports),
     * when not hovering, when not in Select-Node(s) mode, or for a node left over from a swapped-out tree.
     */
    private void paintHoverPreview(final Graphics2D g, final boolean to_screen) {
        // _hover_node is always in the currently-displayed tree: it is only set on hover, and any tree swap
        // clears it via setNodeInPreorderToNull (the shared structural-change chokepoint).
        if (!to_screen || (_hover_node == null) || _hover_node.isCollapse()
                || (getControlPanel().getActionWhenNodeClicked() != NodeClickAction.SELECT_NODES)) {
            return; // no hover circle on a collapsed triangle (its own fill/outline already shows selection)
        }
        final Set<Long> found = getFoundNodes0();
        final boolean deselect;                  // whether a click here would DESELECT rather than add
        final java.util.List<PhylogenyNode> marks;
        if (_hover_subtree) {
            final java.util.List<PhylogenyNode> all = _hover_node.getAllExternalDescendants();
            if (all.isEmpty()) {
                return;
            }
            deselect = allTipsSelected(all); // direction over ALL tips (matches what a click toggles)
            marks = new java.util.ArrayList<PhylogenyNode>();
            collectVisibleTips(_hover_node, marks); // only laid-out tips (skip collapsed sub-clades' stale coords)
        } else {
            deselect = (found != null) && found.contains(_hover_node.getId()); // a click toggles just this node
            marks = java.util.Collections.singletonList(_hover_node);
        }
        final Color mark;
        if (deselect) {
            mark = BRANCH_HOVER_REMOVE; // muted grey (its own alpha) -- "will be de-selected"
        } else {
            final Color f = getTreeColorSet().getFoundColor0(); // "will be added" -> preview in the found colour
            mark = new Color(f.getRed(), f.getGreen(), f.getBlue(), BRANCH_HOVER_ALPHA);
        }
        // the remove mark lands on an already-selected node (with its solid found marker) so it is drawn larger
        // to cover it; the add mark sits on a bare node and reads fine smaller
        final int shape = getOptions().getDefaultNodeShapeSize();
        final int d = deselect ? Math.max(10, shape + 4) : Math.max(6, shape);
        final Color saved = g.getColor();
        g.setColor(mark);
        for (final PhylogenyNode t : marks) {
            // a select marks only the not-yet-selected targets; a deselect (or the single hovered node) marks all
            if (!_hover_subtree || deselect || (found == null) || !found.contains(t.getId())) {
                g.fillOval(Math.round(t.getXcoord()) - (d / 2), Math.round(t.getYcoord()) - (d / 2), d, d);
            }
        }
        g.setColor(saved);
    }

    /** Collects the marks to preview for {@code node}'s subtree: its currently-visible external tips only. A
     *  collapsed sub-clade contributes NO mark -- a hover circle over its triangle is ugly and redundant (the
     *  triangle's own fill/outline already shows the clade's selection state); its hidden tips are still toggled
     *  by the click regardless. A leaf adds itself. */
    private static void collectVisibleTips(final PhylogenyNode node, final java.util.List<PhylogenyNode> out) {
        if (node.isExternal()) {
            out.add(node);
            return;
        }
        if (node.isCollapse()) {
            return; // a collapsed clade is drawn as a single triangle; do not circle it in the hover preview
        }
        for (int i = 0; i < node.getNumberOfDescendants(); ++i) {
            collectVisibleTips(node.getChildNode(i), out);
        }
    }

    /** True if every tip in {@code tips} is already selected (found set 0) -- so a branch click would DESELECT
     *  the clade rather than add to it. Shared by the hover preview and the selectSubtreeTips toggle. */
    private boolean allTipsSelected(final java.util.List<PhylogenyNode> tips) {
        final Set<Long> found = getFoundNodes0();
        if (found == null) {
            return false;
        }
        for (final PhylogenyNode t : tips) {
            if (!found.contains(t.getId())) {
                return false;
            }
        }
        return true;
    }

    /** Test hooks for the hover preview: the hovered node (or null), whether it's a subtree (branch) vs single
     *  node hover, and whether a click there would deselect (vs add). */
    PhylogenyNode hoverNodeForTest() {
        return _hover_node;
    }

    boolean hoverIsSubtreeForTest() {
        return _hover_subtree;
    }

    boolean hoverWouldDeselectForTest() {
        if (_hover_node == null) {
            return false;
        }
        if (_hover_subtree) {
            return allTipsSelected(_hover_node.getAllExternalDescendants());
        }
        final Set<Long> found = getFoundNodes0();
        return (found != null) && found.contains(_hover_node.getId());
    }

    /** Test hook: the node currently suppressed from re-previewing after a click (or null). Only released when
     *  the pointer reaches a different target or the tree changes -- NOT by clearHoverPreview / a mouse-exit. */
    PhylogenyNode clickSuppressedForTest() {
        return _click_suppressed;
    }

    /** Test hook: the marks the hover preview would draw for the current hover -- the subtree's visible tips plus
     *  one mark per collapsed sub-clade (or the single hovered node). Empty when not hovering. */
    java.util.List<PhylogenyNode> hoverPreviewMarksForTest() {
        final java.util.List<PhylogenyNode> marks = new java.util.ArrayList<PhylogenyNode>();
        if (_hover_node == null) {
            return marks;
        }
        if (_hover_subtree) {
            collectVisibleTips(_hover_node, marks);
        } else {
            marks.add(_hover_node);
        }
        return marks;
    }

    final Configuration getConfiguration() {
        return _configuration;
    }

    final ControlPanel getControlPanel() {
        return _control_panel;
    }

    final int getDomainStructureEvalueThresholdExp() {
        return _domain_structure_e_value_thr_exp;
    }

    final Set<Long> getFoundNodes0() {
        return _found_nodes_0;
    }

    final Set<Long> getFoundNodes1() {
        return _found_nodes_1;
    }

    List<PhylogenyNode> getFoundNodesAsListOfPhylogenyNodes() {
        final List<PhylogenyNode> additional_nodes = new ArrayList<>();
        if (getFoundNodes0() != null) {
            for (final Long id : getFoundNodes0()) {
                final PhylogenyNode n = _phylogeny.getNode(id);
                if (n != null) {
                    additional_nodes.add(n);
                }
            }
        }
        if (getFoundNodes1() != null) {
            for (final Long id : getFoundNodes1()) {
                if ((getFoundNodes0() == null) || !getFoundNodes0().contains(id)) {
                    final PhylogenyNode n = _phylogeny.getNode(id);
                    if (n != null) {
                        additional_nodes.add(n);
                    }
                }
            }
        }
        return additional_nodes;
    }

    final Color getGraphicsForNodeBoxWithColorForParentBranch(final PhylogenyNode node) {
        if (getControlPanel().isUseVisualStyles() && (PhylogenyMethods.getBranchColorValue(node) != null)) {
            return (PhylogenyMethods.getBranchColorValue(node));
        } else {
            return (getTreeColorSet().getBranchColor());
        }
    }

    final int getLongestExtNodeInfo() {
        return _longest_ext_node_info;
    }

    final Options getOptions() {
        if (_options == null) {
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

    final void increaseDomainStructureEvalueThresholdExp() {
        if (_domain_structure_e_value_thr_exp < 3) {
            _domain_structure_e_value_thr_exp += 1;
        }
    }

    final void initNodeData() {
        if ((_phylogeny == null) || _phylogeny.isEmpty()) {
            return;
        }
        double _max_original_domain_structure_width = 0.0;
        for (final PhylogenyNode node : _phylogeny.getExternalNodes()) {
            if (node.getNodeData().isHasSequence()
                    && (node.getNodeData().getSequence().getDomainArchitecture() != null)) {
                RenderableDomainArchitecture rds = null;
                if (!(node.getNodeData().getSequence()
                        .getDomainArchitecture() instanceof RenderableDomainArchitecture)) {
                    if (SPECIAL_DOMAIN_COLORING) {
                        rds = new RenderableDomainArchitecture(node.getNodeData().getSequence()
                                .getDomainArchitecture(), node.getName());
                    } else {
                        rds = new RenderableDomainArchitecture(node.getNodeData().getSequence()
                                .getDomainArchitecture());
                    }
                    node.getNodeData().getSequence().setDomainArchitecture(rds);
                } else {
                    rds = (RenderableDomainArchitecture) node.getNodeData().getSequence().getDomainArchitecture();
                }
                if (getControlPanel().isShowDomainArchitectures()) {
                    final double dsw = rds.getOriginalSize().getWidth();
                    if (dsw > _max_original_domain_structure_width) {
                        _max_original_domain_structure_width = dsw;
                    }
                }
            }
        }
        if (getControlPanel().isShowDomainArchitectures()) {
            final float ds_factor_width = (float) (_domain_structure_width / _max_original_domain_structure_width);
            for (final PhylogenyNode node : _phylogeny.getExternalNodes()) {
                if (node.getNodeData().isHasSequence()
                        && (node.getNodeData().getSequence().getDomainArchitecture() != null)) {
                    final RenderableDomainArchitecture rds = (RenderableDomainArchitecture) node.getNodeData()
                            .getSequence().getDomainArchitecture();
                    rds.setRenderingFactorWidth(ds_factor_width * 0.9f);
                    rds.setParameter(_domain_structure_e_value_thr_exp);
                }
            }
        }
    }

    final boolean inOv(final MouseEvent e) {
        return ((e.getX() > (getVisibleRect().x + getOvXPosition() + 1))
                && (e.getX() < ((getVisibleRect().x + getOvXPosition() + getOvMaxWidth()) - 1))
                && (e.getY() > (getVisibleRect().y + getOvYPosition() + 1))
                && (e.getY() < ((getVisibleRect().y + getOvYPosition() + getOvMaxHeight()) - 1)));
    }

    final boolean inOvRectangle(final MouseEvent e) {
        return ((e.getX() >= (getOvRectangle().getX() - 1))
                && (e.getX() <= (getOvRectangle().getX() + getOvRectangle().getWidth() + 1))
                && (e.getY() >= (getOvRectangle().getY() - 1))
                && (e.getY() <= (getOvRectangle().getY() + getOvRectangle().getHeight() + 1)));
    }

    final boolean isCanCollapse() {
        return (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED);
    }

    final boolean isCanUncollapseAll(final PhylogenyNode node) {
        if (node.isExternal() || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            return false;
        }
        if (node.isCollapse()) {
            return true;
        }
        final PhylogenyNodeIterator it = new PreorderTreeIterator(node);
        while (it.hasNext()) {
            if (it.next().isCollapse()) {
                return true;
            }
        }
        return false;
    }

    final boolean isCanColorSubtree() {
        return (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED);
    }

    final boolean isCanCopy() {
        return ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && getOptions().isEditable());
    }

    final boolean isCanCut(final PhylogenyNode node) {
        return ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && getOptions().isEditable()
                && !node.isRoot());
    }

    final boolean isCanDelete() {
        return ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && getOptions().isEditable());
    }

    final boolean isCanPaste() {
        return ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && getOptions().isEditable()
                && (getCutOrCopiedTree() != null) && !getCutOrCopiedTree().isEmpty());
    }

    final boolean isCanReroot() {
        return ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && (_subtree_index < 1));
    }

    final boolean isCanSubtree(final PhylogenyNode node) {
        return ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED) && !node.isExternal()
                && (!node.isRoot() || (_subtree_index > 0)));
    }

    final boolean isCurrentTreeIsSubtree() {
        return (_subtree_index > 0);
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
        if ((_phylogeny == null) || (_phylogeny.getNumberOfExternalNodes() < 2)) {
            return;
        }
        if (!_phylogeny.isRerootable()) {
            JOptionPane
                    .showMessageDialog(this, "This is not rerootable", "Not rerootable", JOptionPane.WARNING_MESSAGE);
            return;
        }
        pushUndoCheckpoint("Midpoint-Root");
        setNodeInPreorderToNull();
        setWaitCursor();
        PhylogenyMethods.midpointRoot(_phylogeny);
        PhylogenyMethods.removeMadConfidences(_phylogeny); // a different rooting; MAD support is stale
        resetNodeIdToDistToLeafMap();
        setArrowCursor();
        setEdited(true);
        repaint();
    }

    final void madRoot() {
        if ((_phylogeny == null) || (_phylogeny.getNumberOfExternalNodes() < 2)) {
            return;
        }
        if (!_phylogeny.isRerootable()) {
            JOptionPane
                    .showMessageDialog(this, "This is not rerootable", "Not rerootable", JOptionPane.WARNING_MESSAGE);
            return;
        }
        pushUndoCheckpoint("MAD-Root");
        setNodeInPreorderToNull();
        setWaitCursor();
        PhylogenyMethods.madRoot(_phylogeny);
        resetNodeIdToDistToLeafMap();
        setArrowCursor();
        setEdited(true);
        repaint();
    }

    final void mouseClicked(final MouseEvent e) {
        if (getOptions().isShowOverview() && isOvOn() && isInOv()) {
            final double w_ratio = getVisibleRect().width / getOvRectangle().getWidth();
            final double h_ratio = getVisibleRect().height / getOvRectangle().getHeight();
            double x = (e.getX() - getVisibleRect().x - getOvXPosition() - (getOvRectangle().getWidth() / 2.0))
                    * w_ratio;
            double y = (e.getY() - getVisibleRect().y - getOvYPosition() - (getOvRectangle().getHeight() / 2.0))
                    * h_ratio;
            if (x < 0) {
                x = 0;
            }
            if (y < 0) {
                y = 0;
            }
            final double max_x = getWidth() - getVisibleRect().width;
            final double max_y = getHeight() - getVisibleRect().height;
            if (x > max_x) {
                x = max_x;
            }
            if (y > max_y) {
                y = max_y;
            }
            getMainPanel().getCurrentScrollPane().getViewport()
                    .setViewPosition(new Point(ForesterUtil.roundToInt(x), ForesterUtil.roundToInt(y)));
            setInOvRect(true);
            repaint();
        } else {
            final PhylogenyNode node = findNode(e.getX(), e.getY());
            if (node != null) {
                if (!node.isRoot() && node.getParent().isCollapse()) {
                    return;
                }
                // A collapsed clade is one unit: a selection-gesture click (shift, or plain in Select-Node(s)
                // mode) toggles its WHOLE subtree's tips -- like a branch click -- so the triangle's fill updates
                // and the count reflects the clade. Ctrl/right-click still opens the node popup on the triangle.
                if (node.isCollapse() && ((e.getModifiersEx() & InputEvent.CTRL_DOWN_MASK) == 0)
                        && !SwingUtilities.isRightMouseButton(e)
                        && (((e.getModifiersEx() & InputEvent.SHIFT_DOWN_MASK) != 0)
                        || (_control_panel.getActionWhenNodeClicked() == NodeClickAction.SELECT_NODES))) {
                    selectSubtreeTips(node);
                    suppressHoverAfterClick(node);
                    repaint();
                    return;
                }
                // Check if shift key is down
                if ((e.getModifiersEx() & InputEvent.SHIFT_DOWN_MASK) != 0) {
                    // Yes, so toggle this node in the selection (found set 0), keeping the count label in sync
                    selectNode(node);
                    suppressHoverAfterClick(node); // don't flip the node's preview to "remove" until pointer moves off
                } else if ((e.getModifiersEx() & InputEvent.CTRL_DOWN_MASK) != 0) {
                    // Check if control key is down
                    // Yes, so pop-up menu
                    displayNodePopupMenu(node, e.getX(), e.getY());
                } else {
                    // Handle unadorned click
                    // Check for right mouse button
                    if (SwingUtilities.isRightMouseButton(e)) {
                        displayNodePopupMenu(node, e.getX(), e.getY());
                    } else {
                        // if not in _found_nodes, clear _found_nodes
                        handleClickToAction(_control_panel.getActionWhenNodeClicked(), node);
                        if (_control_panel.getActionWhenNodeClicked() == NodeClickAction.SELECT_NODES) {
                            suppressHoverAfterClick(node); // just selected this node -> suppress its preview
                        }
                    }
                }
            } else {
                // no node under the cursor: a click on a branch selects/deselects that subtree's external tips
                // (extends "Select Node(s)"; also works via shift-click regardless of the click-to mode)
                final boolean shift = (e.getModifiersEx() & InputEvent.SHIFT_DOWN_MASK) != 0;
                final boolean plain_left = SwingUtilities.isLeftMouseButton(e) && !shift
                        && ((e.getModifiersEx() & InputEvent.CTRL_DOWN_MASK) == 0);
                if (shift || (plain_left
                        && (_control_panel.getActionWhenNodeClicked() == NodeClickAction.SELECT_NODES))) {
                    final PhylogenyNode branch = findBranch(e.getX(), e.getY());
                    selectSubtreeTips(branch);
                    suppressHoverAfterClick(branch); // don't flip to the grey "remove" until the pointer moves off
                }
            }
        }
        repaint();
    }

    final void mouseDragInBrowserPanel(final MouseEvent e) {
        setCursor(MOVE_CURSOR);
        final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
        scroll_position.x -= (e.getX() - getLastDragPointX());
        scroll_position.y -= (e.getY() - getLastDragPointY());
        if (scroll_position.x < 0) {
            scroll_position.x = 0;
        } else {
            final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
            if (scroll_position.x > max_x) {
                scroll_position.x = max_x;
            }
        }
        if (scroll_position.y < 0) {
            scroll_position.y = 0;
        } else {
            final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
            if (scroll_position.y > max_y) {
                scroll_position.y = max_y;
            }
        }
        if (isOvOn() || getOptions().isShowScale()) {
            repaint();
        }
        getMainPanel().getCurrentScrollPane().getViewport().setViewPosition(scroll_position);
    }

    final void mouseDragInOvRectangle(final MouseEvent e) {
        setCursor(HAND_CURSOR);
        final double w_ratio = getVisibleRect().width / getOvRectangle().getWidth();
        final double h_ratio = getVisibleRect().height / getOvRectangle().getHeight();
        final Point scroll_position = getMainPanel().getCurrentScrollPane().getViewport().getViewPosition();
        double dx = ((w_ratio * e.getX()) - (w_ratio * getLastDragPointX()));
        double dy = ((h_ratio * e.getY()) - (h_ratio * getLastDragPointY()));
        scroll_position.x = ForesterUtil.roundToInt(scroll_position.x + dx);
        scroll_position.y = ForesterUtil.roundToInt(scroll_position.y + dy);
        if (scroll_position.x <= 0) {
            scroll_position.x = 0;
            dx = 0;
        } else {
            final int max_x = getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getHorizontalScrollBar().getVisibleAmount();
            if (scroll_position.x >= max_x) {
                dx = 0;
                scroll_position.x = max_x;
            }
        }
        if (scroll_position.y <= 0) {
            dy = 0;
            scroll_position.y = 0;
        } else {
            final int max_y = getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getMaximum()
                    - getMainPanel().getCurrentScrollPane().getVerticalScrollBar().getVisibleAmount();
            if (scroll_position.y >= max_y) {
                dy = 0;
                scroll_position.y = max_y;
            }
        }
        repaint();
        getMainPanel().getCurrentScrollPane().getViewport().setViewPosition(scroll_position);
        setLastMouseDragPointX((float) (e.getX() + dx));
        setLastMouseDragPointY((float) (e.getY() + dy));
    }

    final void mouseMoved(final MouseEvent e) {
        requestFocusInWindow();
        if (getControlPanel().isNodeDescPopup()) {
            if (_node_desc_popup != null) {
                _node_desc_popup.hide();
                _node_desc_popup = null;
            }
        }
        if (getOptions().isShowOverview() && isOvOn()) {
            if (inOvVirtualRectangle(e)) {
                if (!isInOvRect()) {
                    setInOvRect(true);
                    repaint();
                }
            } else {
                if (isInOvRect()) {
                    setInOvRect(false);
                    repaint();
                }
            }
        }
        if (inOv(e) && getOptions().isShowOverview() && isOvOn()) {
            if (!isInOv()) {
                setInOv(true);
            }
            setHover(null, false); // over the overview thumbnail -> no hover preview
        } else {
            if (isInOv()) {
                setInOv(false);
            }
            final boolean select_mode = getControlPanel().getActionWhenNodeClicked() == NodeClickAction.SELECT_NODES;
            final PhylogenyNode node = findNode(e.getX(), e.getY());
            if ((node != null) && (node.isRoot() || !node.getParent().isCollapse())) {
                if ((getControlPanel().getActionWhenNodeClicked() == NodeClickAction.CUT_SUBTREE)
                        || (getControlPanel().getActionWhenNodeClicked() == NodeClickAction.COPY_SUBTREE)
                        || (getControlPanel().getActionWhenNodeClicked() == NodeClickAction.PASTE_SUBTREE)
                        || (getControlPanel().getActionWhenNodeClicked() == NodeClickAction.DELETE_NODE_OR_SUBTREE)
                        || (getControlPanel().getActionWhenNodeClicked() == NodeClickAction.REROOT)
                        || (getControlPanel().getActionWhenNodeClicked() == NodeClickAction.ADD_NEW_NODE)) {
                    setCursor(CUT_CURSOR);
                } else {
                    setCursor(HAND_CURSOR);
                    if (getControlPanel().isNodeDescPopup()) {
                        showNodeDataPopup(e, node);
                    }
                }
                // in Select-Node(s) mode, preview the single node a click would toggle; otherwise no preview.
                // A collapsed clade is excluded: a hover circle over its triangle is ugly and redundant -- the
                // hand cursor shows it is clickable and the triangle's own fill/outline shows its selection state.
                applyHover((select_mode && !node.isCollapse()) ? node : null, false);
            } else {
                // not over a node: over a branch (Select-Node(s) mode) -> hand cursor + preview the subtree
                final PhylogenyNode branch = select_mode ? findBranch(e.getX(), e.getY()) : null;
                applyHover(branch, true);
                setCursor((branch != null) ? HAND_CURSOR : ARROW_CURSOR);
            }
        }
    }

    final void mouseReleasedInBrowserPanel(final MouseEvent e) {
        setCursor(ARROW_CURSOR);
    }

    final void multiplyUrtFactor(final float f) {
        _urt_factor *= f;
    }

    final void paintBranchCircular(final PhylogenyNode p,
                                   final PhylogenyNode c,
                                   final Graphics2D g,
                                   final boolean radial_labels,
                                   final boolean to_pdf,
                                   final boolean to_graphics_file) {
        final double angle = _urt_nodeid_angle_map.get(c.getId());
        final double root_x = _root.getXcoord();
        final double root_y = _root.getYcoord();
        final double dx = root_x - p.getXcoord();
        final double dy = root_y - p.getYcoord();
        final double parent_radius = Math.sqrt((dx * dx) + (dy * dy));
        final double arc = (_urt_nodeid_angle_map.get(p.getId())) - angle;
        assignGraphicsForBranchWithColorForParentBranch(c, false, g, to_pdf, to_graphics_file);
        if ((c.isFirstChildNode() || c.isLastChildNode())
                && ((Math.abs(parent_radius * arc) > 1.5) || to_pdf || to_graphics_file)) {
            final double r2 = 2.0 * parent_radius;
            drawArc(root_x - parent_radius, root_y - parent_radius, r2, r2, (-angle - arc), arc, g);
        }
        drawLine(c.getXcoord(),
                c.getYcoord(),
                root_x + (Math.cos(angle) * parent_radius),
                root_y + (Math.sin(angle) * parent_radius),
                g);
        paintNodeBox(c.getXcoord(), c.getYcoord(), c, g, to_pdf, to_graphics_file);
        if (c.isExternal()) {
            final boolean is_in_found_nodes = isInFoundNodes0(c) || isInFoundNodes1(c);
            if ((_dynamic_hiding_factor > 1) && !is_in_found_nodes
                    && ((_urt_nodeid_index_map.get(c.getId()) % _dynamic_hiding_factor) != 1)) {
                return;
            }
            paintNodeDataUnrootedCirc(g, c, to_pdf, to_graphics_file, radial_labels, 0, is_in_found_nodes);
        }
    }

    final void paintBranchCircularLite(final PhylogenyNode p, final PhylogenyNode c, final Graphics2D g) {
        final double angle = _urt_nodeid_angle_map.get(c.getId());
        final double root_x = _root.getXSecondary();
        final double root_y = _root.getYSecondary();
        final double dx = root_x - p.getXSecondary();
        final double dy = root_y - p.getYSecondary();
        final double arc = (_urt_nodeid_angle_map.get(p.getId())) - angle;
        final double parent_radius = Math.sqrt((dx * dx) + (dy * dy));
        g.setColor(getTreeColorSet().getOvColor());
        if ((c.isFirstChildNode() || c.isLastChildNode()) && (Math.abs(arc) > 0.02)) {
            final double r2 = 2.0 * parent_radius;
            drawArc(root_x - parent_radius, root_y - parent_radius, r2, r2, (-angle - arc), arc, g);
        }
        drawLine(c.getXSecondary(),
                c.getYSecondary(),
                root_x + (Math.cos(angle) * parent_radius),
                root_y + (Math.sin(angle) * parent_radius),
                g);
        if (isInFoundNodes(c)) {
            g.setColor(getColorForFoundNode(c));
            drawRectFilled(c.getXSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                    c.getYSecondary() - OVERVIEW_FOUND_NODE_BOX_SIZE_HALF,
                    OVERVIEW_FOUND_NODE_BOX_SIZE,
                    OVERVIEW_FOUND_NODE_BOX_SIZE,
                    g);
        }
    }

    final void paintCircular(final Phylogeny phy,
                             final double starting_angle,
                             final int center_x,
                             final int center_y,
                             final int radius,
                             final Graphics2D g,
                             final boolean to_pdf,
                             final boolean to_graphics_file) {
        final int circ_num_ext_nodes = phy.getNumberOfExternalNodes() - _collapsed_external_nodeid_set.size();
        _root = phy.getRoot();
        _root.setXcoord(center_x);
        _root.setYcoord(center_y);
        final boolean radial_labels = getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL;
        double current_angle = starting_angle;
        int i = 0;
        for (final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if (!n.isCollapse()) {
                n.setXcoord((float) (center_x + (radius * Math.cos(current_angle))));
                n.setYcoord((float) (center_y + (radius * Math.sin(current_angle))));
                _urt_nodeid_angle_map.put(n.getId(), current_angle);
                _urt_nodeid_index_map.put(n.getId(), i++);
                current_angle += (TWO_PI / circ_num_ext_nodes);
            }
        }
        paintCirculars(phy.getRoot(), phy, center_x, center_y, radius, radial_labels, g, to_pdf, to_graphics_file);
        paintNodeBox(_root.getXcoord(), _root.getYcoord(), _root, g, to_pdf, to_graphics_file);
    }

    final void paintCircularLite(final Phylogeny phy,
                                 final double starting_angle,
                                 final int center_x,
                                 final int center_y,
                                 final int radius,
                                 final Graphics2D g) {
        final int circ_num_ext_nodes = phy.getNumberOfExternalNodes();
        _root = phy.getRoot();
        _root.setXSecondary(center_x);
        _root.setYSecondary(center_y);
        double current_angle = starting_angle;
        for (final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            n.setXSecondary((float) (center_x + (radius * Math.cos(current_angle))));
            n.setYSecondary((float) (center_y + (radius * Math.sin(current_angle))));
            _urt_nodeid_angle_map.put(n.getId(), current_angle);
            current_angle += (TWO_PI / circ_num_ext_nodes);
        }
        paintCircularsLite(phy.getRoot(), phy, center_x, center_y, radius, g);
    }

    final void paintPhylogeny(final Graphics2D g,
                              final boolean to_pdf,
                              final boolean to_graphics_file,
                              final int graphics_file_width,
                              final int graphics_file_height,
                              final int graphics_file_x,
                              final int graphics_file_y) {
        if ((_phylogeny == null) || _phylogeny.isEmpty()) {
            return;
        }
        // The support scale ceiling (a single preorder scan) feeds both the support symbols and the "min.
        // confidence shown" label filter (a fraction of this ceiling). Skip the scan -- keeping the cheap 1.0
        // fallback -- when neither consumer is active (the default), so plain repaints on the hot hover/scroll/
        // zoom path don't pay for it.
        _confidence_scale_max = ((getOptions().getSupportVisualization() != Options.SUPPORT_VISUALIZATION.NONE)
                || getControlPanel().isShowConfidenceValues())
                        ? TreePanelUtil.detectConfidenceScaleMax(_phylogeny)
                        : 1.0;
        // "Dim Non-Matches": resolve once per paint whether any hit is actually visible (see anyVisibleFoundNode);
        // the scan is skipped entirely -- keeping the cheap false -- when the option is off (the default).
        _has_visible_found_node = getOptions().isDimNonMatches() && anyVisibleFoundNode();
        // "Pulse Found Nodes": clear the screen halo state each paint (the rectangular branch repopulates it via
        // paintFoundNodeHalos); the timer is reconciled at the END of this method for EVERY layout, so a switch to
        // circular/unrooted (which doesn't repopulate) stops it. Screen-only, so an export never clobbers it.
        if (!to_pdf && !to_graphics_file) {
            _found_halo_bounds.clear();
            _has_visible_found_halo = false;
        }
        if (_control_panel.isShowSequenceRelations()) {
            _query_sequence = _control_panel.getSelectedQuerySequence();
        }
        // Color the background
        if (!to_pdf) {
            final Rectangle r = getVisibleRect();
            g.setColor(getTreeColorSet().getBackgroundColor());
            if (!to_graphics_file) {
                g.fill(r);
            } else if (!_export_transparent_background) {
                if (getOptions().isPrintBlackAndWhite()) {
                    g.setColor(Color.WHITE);
                }
                g.fillRect(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            }
            // else: transparent PNG export -- leave the (ARGB) canvas unfilled
            setupStroke(g);
        } else {
            g.setStroke(new BasicStroke(getOptions().getPrintLineWidth()));
        }
        if ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
            _external_node_index = 0;
            // Position starting X of tree
            if (!_phylogeny.isRooted() /*|| ( _subtree_index > 0 )*/) {
                _phylogeny.getRoot().setXcoord(TreePanel.MOVE);
            } else if ((_phylogeny.getRoot().getDistanceToParent() > 0.0) && getControlPanel().isDrawPhylogram()) {
                _phylogeny.getRoot().setXcoord((float) (TreePanel.MOVE
                        + (_phylogeny.getRoot().getDistanceToParent() * getXcorrectionFactor())));
            } else {
                _phylogeny.getRoot().setXcoord(TreePanel.MOVE + getXdistance());
            }
            // Position starting Y of tree (shifted down by the same header reserve used in calcParametersForPainting)
            _phylogeny.getRoot().setYcoord((getYdistance() * _phylogeny.getRoot().getNumberOfExternalNodes())
                    + (TreePanel.MOVE / 2.0f) + verticalBreadthPad() + annotationHeaderTopReserve());
            final int dynamic_hiding_factor = calcDynamicHidingFactor();
            if (getControlPanel().isDynamicallyHideData()) {
                if (dynamic_hiding_factor > 1) {
                    getControlPanel().setDynamicHidingIsOn(true);
                } else {
                    getControlPanel().setDynamicHidingIsOn(false);
                }
            }
            if (_nodes_in_preorder == null) {
                _nodes_in_preorder = new PhylogenyNode[_phylogeny.getNodeCount()];
                int i = 0;
                for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
                    _nodes_in_preorder[i++] = it.next();
                }
            }
            final boolean disallow_shortcutting = (dynamic_hiding_factor < 40)
                    /* || getControlPanel().isUseVisualStyles() || getOptions().isShowDefaultNodeShapesForMarkedNodes()*/ //TODO check if this is really not needed.
                    || to_graphics_file || to_pdf;
            final boolean vertical = isVerticalOrientation();
            final boolean scale_grid_shown = getOptions().isShowScaleGrid() && getControlPanel().isDrawPhylogram()
                    && (getScaleDistance() > 0.0);
            if (!vertical && scale_grid_shown) {
                paintScaleGrid(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height);
            }
            // Root-top/bottom: the tree is laid out logically (above); now rotate the whole canvas for the geometry
            // pass. Geometry (branches, boxes, triangles, halos) rides R for free; node TEXT is re-anchored upright
            // or at 45deg by withNodeTextFrame inside paintNodeRectangular; the viewport-fixed chrome further below
            // is drawn after the base transform is restored. COMPOSE (g.transform), never setTransform, so an
            // export/print backend's own scale/translate survives.
            final AffineTransform orientation_saved = g.getTransform();
            if (vertical) {
                rebuildOrientationTransform();
                _orientation_base_transform = orientation_saved;
                g.transform(_orientation_R);
                // vertical parity: draw the scale grid INSIDE the R frame (before the node loop, behind the branches)
                // so its logical vertical lines ride R into horizontal grid lines at each depth interval
                if (scale_grid_shown) {
                    paintScaleGrid(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height);
                }
            }
            for (final PhylogenyNode element : _nodes_in_preorder) {
                paintNodeRectangular(g,
                        element,
                        to_pdf,
                        getControlPanel().isDynamicallyHideData() && (dynamic_hiding_factor > 1),
                        dynamic_hiding_factor,
                        to_graphics_file,
                        disallow_shortcutting);
            }
            if (!vertical) {
                // These tree-riding overlays are DEFERRED in the vertical orientation (increment 1): they hardcode a
                // horizontal tip edge, so they are simply not drawn rather than drawn wrong (they will ride R in a
                // later increment). faint alternating row bands first, so annotation columns etc. sit on top.
                paintZebraStripes(g, to_pdf, to_graphics_file, graphics_file_x, graphics_file_width);
                paintHpdBars(g, to_pdf, to_graphics_file); // node-age HPD bars -- node coords set by the loop above
                paintAnnotationColumns(g); // tip-aligned columns (strip/heat map/bar/text), right of the labels
                paintCladeBands(g); // clade boxes/bars over the tree -- node coords set by the loop above
            }
            else {
                // vertical parity: these overlays are drawn while g is rotated by R. Their geometry (zebra bands,
                // annotation cells, clade boxes/bars/brackets, HPD bars) is axis-aligned rects + lines, so it rides R
                // for free; the TEXT (annotation cells/headers, clade labels) is re-anchored upright inside the paint
                // methods. The scale grid rides R inside the R block above; the labeled scale axis is a side ruler
                // drawn as upright chrome after the base frame is restored (paintScaleAxisVertical, further below).
                paintZebraStripes(g, to_pdf, to_graphics_file, graphics_file_x, graphics_file_width); // faint row bands, behind
                paintHpdBars(g, to_pdf, to_graphics_file); // node-age HPD bars: a plain rect at each node -> rides R
                paintAnnotationColumnsVertical(g);
                paintCladeBands(g); // boxes ride R; bars/brackets draw the label upright (isVerticalOrientation branch)
            }
            paintHoverPreview(g, !(to_pdf || to_graphics_file)); // translucent select/deselect hover preview (rides R)
            paintFoundNodeHalos(g, to_pdf, to_graphics_file); // pulsing (screen) / static-glow (export) hit halos
            // restore the upright base frame before the viewport-fixed chrome (scale bar, tree name, overview, legends)
            if (vertical) {
                g.setTransform(orientation_saved);
            }
            final boolean scale_shown = getOptions().isShowScale() && getControlPanel().isDrawPhylogram()
                    && (getScaleDistance() > 0.0);
            final boolean axis_shown = getOptions().isShowScaleAxis() && getControlPanel().isDrawPhylogram()
                    && (getScaleDistance() > 0.0);
            // the horizontal axis owns a reserved bottom band; lift the (viewport-fixed) scale bar clear above it (the
            // tree name is likewise raised, inside paintTreeName) so the three bottom overlays never overprint. Derive
            // both the lift AND the flag from the SAME layout reserve so they stay in lockstep with what is actually
            // drawn/reserved -- 0 in a vertical orientation (side ruler), for a cladogram / circular / unrooted, or an
            // unticked scale (never lift over a band that isn't there).
            final int bottom_reserve = scaleAxisBottomReserve();
            final boolean axis_shown_horizontal = bottom_reserve > 0;
            if (scale_shown) {
                if (!(to_graphics_file || to_pdf)) {
                    paintScale(g,
                            getVisibleRect().x,
                            (getVisibleRect().y + getVisibleRect().height) - bottom_reserve,
                            to_pdf,
                            to_graphics_file);
                } else {
                    paintScale(g, graphics_file_x, (graphics_file_y + graphics_file_height) - bottom_reserve, to_pdf,
                            to_graphics_file);
                }
            }
            if (axis_shown) {
                if (vertical) {
                    paintScaleAxisVertical(g, to_pdf, to_graphics_file); // ruler down the breadth side, upright labels
                } else {
                    paintScaleAxis(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height);
                }
            }
            if (getOptions().isShowTreeName() && !ForesterUtil.isEmpty(getPhylogeny().getName())) {
                // the name sits in the lower-left, but slides to the lower-right when the scale is shown there, and
                // is raised above the scale axis when that occupies the bottom strip -- so it never overlaps either
                if (!(to_graphics_file || to_pdf)) {
                    paintTreeName(g,
                            getVisibleRect().x,
                            getVisibleRect().width,
                            getVisibleRect().y + getVisibleRect().height,
                            to_pdf,
                            to_graphics_file,
                            scale_shown,
                            axis_shown_horizontal);
                } else {
                    paintTreeName(g, graphics_file_x, graphics_file_width, graphics_file_y + graphics_file_height,
                            to_pdf, to_graphics_file, scale_shown, axis_shown_horizontal);
                }
            }
            if (getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf) {
                if (vertical) {
                    paintPhylogenyLiteVertical(g); // rotated thumbnail (scales the main tree into the overview box)
                } else {
                    paintPhylogenyLite(g);
                }
            }
        } else if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            if (getControlPanel().getDynamicallyHideData() != null) {
                getControlPanel().setDynamicHidingIsOn(false);
            }
            final double angle = getStartingAngle();
            final boolean radial_labels = getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL;
            _dynamic_hiding_factor = 0;
            if (getControlPanel().isDynamicallyHideData()) {
                _dynamic_hiding_factor = (int) ((getFontMetricsForLargeDefaultFont().getHeight() * 1.5
                        * getPhylogeny().getNumberOfExternalNodes()) / (TWO_PI * 10));
            }
            if (getControlPanel().getDynamicallyHideData() != null) {
                if (_dynamic_hiding_factor > 1) {
                    getControlPanel().setDynamicHidingIsOn(true);
                } else {
                    getControlPanel().setDynamicHidingIsOn(false);
                }
            }
            paintUnrooted(_phylogeny.getRoot(),
                    angle,
                    (float) (angle + (2 * Math.PI)),
                    radial_labels,
                    g,
                    to_pdf,
                    to_graphics_file);
            if (getOptions().isShowScale()) {
                if (!(to_graphics_file || to_pdf)) {
                    paintScale(g,
                            getVisibleRect().x,
                            getVisibleRect().y + getVisibleRect().height,
                            to_pdf,
                            to_graphics_file);
                } else {
                    paintScale(g, graphics_file_x, graphics_file_y + graphics_file_height, to_pdf, to_graphics_file);
                }
            }
            if (getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf) {
                g.setColor(getTreeColorSet().getOvColor());
                paintUnrootedLite(_phylogeny.getRoot(),
                        angle,
                        angle + (2 * Math.PI),
                        g,
                        (getUrtFactorOv() / (getVisibleRect().width / getOvMaxWidth())));
                paintOvRectangle(g);
            }
        } else if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) {
            final int radius = (int) ((Math.min(getPreferredSize().getWidth(), getPreferredSize().getHeight())
                    / 2) - (MOVE + getLongestExtNodeInfo()));
            final int d = radius + MOVE + getLongestExtNodeInfo();
            _dynamic_hiding_factor = 0;
            if (getControlPanel().isDynamicallyHideData() && (radius > 0)) {
                _dynamic_hiding_factor = (int) ((getFontMetricsForLargeDefaultFont().getHeight() * 1.5
                        * getPhylogeny().getNumberOfExternalNodes()) / (TWO_PI * radius));
            }
            if (getControlPanel().getDynamicallyHideData() != null) {
                if (_dynamic_hiding_factor > 1) {
                    getControlPanel().setDynamicHidingIsOn(true);
                } else {
                    getControlPanel().setDynamicHidingIsOn(false);
                }
            }
            paintCircular(_phylogeny, getStartingAngle(), d, d, radius > 0 ? radius : 0, g, to_pdf, to_graphics_file);
            if (getOptions().isShowOverview() && isOvOn() && !to_graphics_file && !to_pdf) {
                final int radius_ov = (int) (getOvMaxHeight() < getOvMaxWidth() ? getOvMaxHeight() / 2
                        : getOvMaxWidth() / 2);
                double x_scale = 1.0;
                double y_scale = 1.0;
                int x_pos = getVisibleRect().x + getOvXPosition();
                int y_pos = getVisibleRect().y + getOvYPosition();
                if (getWidth() > getHeight()) {
                    x_scale = (double) getHeight() / getWidth();
                    x_pos = ForesterUtil.roundToInt(x_pos / x_scale);
                } else {
                    y_scale = (double) getWidth() / getHeight();
                    y_pos = ForesterUtil.roundToInt(y_pos / y_scale);
                }
                _at = g.getTransform();
                g.scale(x_scale, y_scale);
                paintCircularLite(_phylogeny,
                        getStartingAngle(),
                        x_pos + radius_ov,
                        y_pos + radius_ov,
                        (int) (radius_ov - (getLongestExtNodeInfo()
                                / (getVisibleRect().width / getOvRectangle().getWidth()))),
                        g);
                g.setTransform(_at);
                paintOvRectangle(g);
            }
        }
        if (hasAnnotationColumnLegend() || isColorByProperty() || hasRankLegend()) {
            final boolean to_screen = !(to_pdf || to_graphics_file);
            final Rectangle legend_bounds = to_screen
                    ? getVisibleRect()
                    : new Rectangle(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            // one legend slot; a header-focused annotation-column legend wins (the user just clicked it),
            // else the property-color legend, else the rank legend
            if (hasAnnotationColumnLegend()) {
                drawAnnotationColumnLegend(g, legend_bounds, to_screen);
            } else if (isColorByProperty()) {
                drawPropertyColorLegend(g, legend_bounds, to_screen);
            } else {
                drawRankLegend(g, legend_bounds, to_screen);
            }
        }
        // "Size by" has its OWN legend (a separate, independently placed key), drawn last so it can appear
        // ALONGSIDE the color/rank legend -- the whole point of the combined color+size figure.
        if (isSizeByProperty()) {
            final boolean to_screen = !(to_pdf || to_graphics_file);
            final Rectangle legend_bounds = to_screen ? getVisibleRect()
                    : new Rectangle(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            drawSizeLegend(g, legend_bounds, to_screen);
        }
        // reconcile the "Pulse Found Nodes" animation timer after EVERY screen paint (all layouts): starts it when a
        // hit halo was drawn on a rectangular tree, stops it when none was (option off / no hit / circular-unrooted).
        if (!to_pdf && !to_graphics_file) {
            updatePulseTimer();
        }
    }

    final void recalculateMaxDistanceToRoot() {
        _max_distance_to_root = PhylogenyMethods.calculateMaxDistanceToRoot(getPhylogeny());
        if (getPhylogeny().getRoot().getDistanceToParent() > 0) {
            _max_distance_to_root += getPhylogeny().getRoot().getDistanceToParent();
        }
    }

    /**
     * Remove all edit-node frames
     */
    final void removeAllEditNodeJFrames() {
        for (int i = 0; i <= (TreePanel.MAX_NODE_FRAMES - 1); i++) {
            if (_node_frames[i] != null) {
                _node_frames[i].dispose();
                _node_frames[i] = null;
            }
        }
        _node_frame_index = 0;
    }

    /**
     * Remove a node-edit frame.
     */
    final void removeEditNodeFrame(final int i) {
        _node_frame_index--;
        _node_frames[i] = null;
        if (i < _node_frame_index) {
            for (int j = 0; j < (_node_frame_index - 1); j++) {
                _node_frames[j] = _node_frames[j + 1];
            }
            _node_frames[_node_frame_index] = null;
        }
    }

    final void reRoot(final PhylogenyNode node) {
        if (!getPhylogeny().isRerootable()) {
            JOptionPane
                    .showMessageDialog(this, "This is not rerootable", "Not rerootable", JOptionPane.WARNING_MESSAGE);
            return;
        }
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            JOptionPane.showMessageDialog(this,
                    "Cannot reroot in unrooted display type",
                    "Attempt to reroot tree in unrooted display",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }
        pushUndoCheckpoint("Re-Root");
        if (!node.isRoot()) {
            // a different rooting was chosen manually; any MAD root support is now stale
            PhylogenyMethods.removeMadConfidences(getPhylogeny());
        }
        getPhylogeny().reRoot(node);
        getPhylogeny().recalculateNumberOfExternalDescendants(true);
        resetNodeIdToDistToLeafMap();
        setNodeInPreorderToNull();
        resetPreferredSize();
        getMainPanel().adjustJScrollPane();
        setEdited(true);
        repaint();
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) {
            getControlPanel().showWhole();
        }
    }

    final void resetNodeIdToDistToLeafMap() {
        _nodeid_dist_to_leaf = new HashMap<>();
    }

    final void resetPreferredSize() {
        if ((getPhylogeny() == null) || getPhylogeny().isEmpty()) {
            return;
        }
        final int[] ext = logicalTreeExtent();
        final int w = ext[0]; // logical width  = depth (root->tip along Xcoord) + label column
        final int h = ext[1]; // logical height = breadth (tip spread along Ycoord)
        if (isVerticalOrientation()) {
            // root-top/bottom: depth (w) becomes the vertical extent, breadth (h) the horizontal. w already includes the
            // downward tip-label reach; the SIDEWAYS reach (breadthLabelReserve, the outermost tilted label extending
            // past its tip) is added to the width here so the label is not clipped -- the fit reserved the same amount.
            setPreferredSize(new Dimension(h + breadthLabelReserve(), w));
        } else {
            setPreferredSize(new Dimension(w, h));
        }
        invalidateOrientationTransform(); // the logical extents may have changed -> R must be rebuilt on the next paint
    }

    /** The tip-label footprint reserved along the DEPTH axis (root-&gt;tip). In a vertical orientation the tip labels
     *  tilt, so the depth reach is the label bounding box's projection onto the depth: {@code L*|sin| + lineH*|cos|}
     *  -- the full width L for a 90-degree label, L/sqrt(2)+ for a 45-degree one, and one line height for an upright
     *  (0-degree) label sitting below/above the tip. In the horizontal orientation the label lies flat, so the full
     *  width is reserved (unchanged behavior). */
    private int depthLabelReserve() {
        if (!isVerticalOrientation()) {
            return getLongestExtNodeInfo();
        }
        final double a = tipLabelAngle();
        final int line_h = getFontMetricsForLargeDefaultFont().getHeight();
        // the TEXT label projected onto the depth (L_text*|sin| + lineH*|cos|) + the axis-aligned domain track past it.
        // Uses the TEXT-ONLY longest (not getLongestExtNodeInfo, which already folds in the domain width -- using that
        // here would count the domain twice, over-compressing the depth axis when domains are shown).
        return (int) Math.ceil((_length_of_longest_text_only * Math.abs(Math.sin(a))) + (line_h * Math.abs(Math.cos(a))))
                + verticalDomainReserve();
    }

    /** Extra depth (px) reserved past the tilted tip labels for the axis-aligned domain-architecture track that a
     *  vertical orientation draws as a per-tip vertical bar (0 unless domains are shown). The label reserve above tilts
     *  by sin(angle), but the domain boxes are axis-aligned (their FULL rendered width runs along the depth), so add
     *  the widest rendered track + the label gap + the render's internal 20px lead-in. */
    private int verticalDomainReserve() {
        if (!isVerticalOrientation() || !getControlPanel().isShowDomainArchitectures()) {
            return 0;
        }
        return (int) Math.ceil(_longest_rendered_domain) + VERTICAL_DOMAIN_GAP + 20;
    }

    /** The tip-label footprint reserved along the BREADTH axis (tip spread), which exists only in a vertical
     *  orientation: the labels' sideways component -- ~a line height for a 90-degree label, L/sqrt(2) (plus a line
     *  height) for a 45-degree one -- so the edge tips' labels do not clip and the fit does not overflow the width.
     *  Zero in the horizontal orientation, where tip labels do not extend along the breadth axis. */
    private int breadthLabelReserve() {
        if (!isVerticalOrientation()) {
            return 0;
        }
        if (effectiveTipLabelDirection() == Options.TIP_LABEL_DIRECTION.HORIZONTAL) {
            return 0; // an upright label is CENTRED on its tip, so its L/2-per-side reach is reserved symmetrically by
                      // verticalBreadthPad() (both edge tips fit), not one-sided here
        }
        return (int) Math.ceil((_length_of_longest_text_only * Math.abs(Math.cos(tipLabelAngle())))
                + getFontMetricsForLargeDefaultFont().getHeight());
    }

    /** The logical breadth (tip-spread) extent -- the {@code h} element of {@link #logicalTreeExtent()}, but computed
     *  WITHOUT the depth (so no {@code calculateHeight} tree walk). The vertical scale grid + scale-axis ruler paint
     *  need only the breadth span, so they call this per repaint instead of the full extent. Grows by the vertical
     *  scale-axis reserve so the ruler band is inside the canvas; the tree is shifted down by MOVE + the header reserve
     *  (see the root Ycoord in paintPhylogeny), so the scroll/paint extents agree. */
    private int treeBreadthExtent() {
        return TreePanel.MOVE + (2 * verticalBreadthPad()) + annotationHeaderTopReserve() + verticalScaleAxisReserve()
                + scaleAxisBottomReserve() // horizontal-orientation bottom axis band (0 in vertical; mutually exclusive)
                + ForesterUtil.roundToInt(getYdistance() * getPhylogeny().getRoot().getNumberOfExternalNodes() * 2);
    }

    /** The per-side breadth reserve for a fitted vertical (root-top/bottom) tree: a small aesthetic margin so the tree
     *  is not flush against the left/right edges, PLUS -- for an upright (0-degree) tip label, which is CENTRED on its
     *  tip and reaches L/2 to EACH breadth side -- that half-label, so both edge tips fit (tilted labels reach one side
     *  and are handled by breadthLabelReserve). 0 in the horizontal orientation. Reserved symmetrically: {@code 2x} in
     *  the breadth extent + the fit budget, and the tree origin shifted by {@code 1x} -- keep the three call sites in
     *  step. */
    private int verticalBreadthPad() {
        if (!isVerticalOrientation()) {
            return 0;
        }
        final int centred_half = (effectiveTipLabelDirection() == Options.TIP_LABEL_DIRECTION.HORIZONTAL)
                ? ((_length_of_longest_text_only + 1) / 2) : 0;
        return VERTICAL_BREADTH_PAD + centred_half;
    }

    /** The tree's LOGICAL (root-on-left) extent {width, height}: width = depth + label column, height = tip
     *  spread. A single source of truth for both {@link #resetPreferredSize()} and the orientation transform R
     *  (so the scroll extent and R's output box agree). In a vertical orientation the tip-label footprint is split
     *  by the label tilt: its depth component grows the width, its breadth component grows the height. */
    private int[] logicalTreeExtent() {
        // include the annotation-header top reserve: paintPhylogeny shifts the whole tree down by it (see the
        // root Ycoord), so the scrollable height must grow by the same amount or the bottom tips/cells clip
        // h is R's breadth translate (the tip-spread extent, WITHOUT the sideways tip-label reach): the tilted labels
        // extend PAST the outermost tip, so their breadth reach is added to the preferred WIDTH / fitHeight budget (see
        // resetPreferredSize / logicalBreadthExtent), not here -- otherwise R would push the outermost label off-canvas.
        final int h = treeBreadthExtent();
        int w;
        if (getControlPanel().isDrawPhylogram()) {
            w = TreePanel.MOVE + depthLabelReserve() + rightMarginExtraWidth()
                    + ForesterUtil.roundToInt((getXcorrectionFactor()
                    * getPhylogeny().calculateHeight(!_options.isCollapsedWithAverageHeigh()))
                    + getXdistance());
        } else if (!isNonLinedUpCladogram()) {
            w = TreePanel.MOVE + depthLabelReserve() + rightMarginExtraWidth() + ForesterUtil
                    .roundToInt(getXdistance() * (getPhylogeny().getRoot().getNumberOfExternalNodes() + 2));
        } else {
            w = TreePanel.MOVE + depthLabelReserve() + rightMarginExtraWidth() + ForesterUtil
                    .roundToInt(getXdistance() * (PhylogenyMethods.calculateMaxDepth(getPhylogeny()) + 1));
        }
        return new int[] { w, h };
    }

    /** The tree's LOGICAL breadth (tip-spread) extent in px -- the logical height, which equals the vertical
     *  orientation's preferred WIDTH (breadth axis). fitHeight keeps this as its breadth budget so the horizontal
     *  (breadth) zoom doesn't drift each press: it includes the same breadth tip-label reserve the fit subtracts,
     *  so feeding it back reproduces the current y-distance exactly (the fitHeight-idempotence property). */
    final int logicalBreadthExtent() {
        return logicalTreeExtent()[1] + breadthLabelReserve();
    }

    /** True when the tree is drawn in a vertical (root-top / root-bottom) orientation. Always false for the
     *  radial CIRCULAR/UNROOTED layouts (orientation is a rectangular-family concept). */
    final boolean isVerticalOrientation() {
        final Options.TREE_ORIENTATION o = getOptions().getTreeOrientation();
        return ((o == Options.TREE_ORIENTATION.ROOT_TOP) || (o == Options.TREE_ORIENTATION.ROOT_BOTTOM))
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED);
    }

    /** (Re)builds the logical->screen rotation R (and its inverse) for the current vertical orientation, from the
     *  logical extents. Pure rotations (determinant +1, no mirror); the translate keeps the tree in the positive
     *  quadrant. ROOT_TOP turns the page 90deg clockwise: (x,y)->(H-y, x); ROOT_BOTTOM 90deg CCW: (x,y)->(y, W-x). */
    final void rebuildOrientationTransform() {
        final Options.TREE_ORIENTATION current = getOptions().getTreeOrientation();
        if (!_orientation_transform_dirty && (_orientation_R != null) && (_orientation_R_built_for == current)) {
            return; // cached R is still valid -- no layout/structure change since it was last built
        }
        final int[] ext = logicalTreeExtent();
        final AffineTransform r = TreePanelUtil.orientationTransformFor(current, ext[0], ext[1]);
        _orientation_R = r;
        try {
            _orientation_R_inverse = r.createInverse();
        } catch (final java.awt.geom.NoninvertibleTransformException e) {
            _orientation_R_inverse = new AffineTransform(); // a pure rotation is always invertible; identity fallback
        }
        _orientation_R_built_for = current;
        _orientation_transform_dirty = false;
    }

    /** Marks the cached orientation transform R stale; called from the layout chokepoints (never during a repaint). */
    private void invalidateOrientationTransform() {
        _orientation_transform_dirty = true;
    }

    /** Test hook: the current (possibly cached) orientation transform R -- object identity reveals a rebuild vs reuse. */
    AffineTransform getOrientationRForTest() {
        return _orientation_R;
    }

    /** Test hook: whether this node's data is screen-culled at the current viewport + orientation. */
    boolean isNodeDataInvisibleForTest(final PhylogenyNode node) {
        return isNodeDataInvisible(node);
    }

    /** Test hook: the resolved tip-label tilt (radians) for the current "Tip label angle" setting + orientation. */
    double tipLabelAngleForTest() {
        return tipLabelAngle();
    }

    /** Test hook: a node's position in the (rotated) overview thumbnail, in VIEWPORT-relative coords (matching a
     *  {@code JViewport.printAll} image) -- the paint's {@code getVisibleRect()} translate and printAll's own
     *  {@code -visibleRect} translate cancel, so only the overview corner offset remains. */
    Point2D.Double overviewPointForTest(final PhylogenyNode node) {
        if (_orientation_R == null) {
            return null;
        }
        final double sx = getOvMaxWidth() / (double) getWidth();
        final double sy = getOvMaxHeight() / (double) getHeight();
        final AffineTransform t = AffineTransform.getTranslateInstance(getOvXPosition(), getOvYPosition());
        t.scale(sx, sy);
        t.concatenate(_orientation_R);
        final Point2D.Double p = new Point2D.Double(node.getXcoord(), node.getYcoord());
        t.transform(p, p);
        return p;
    }

    /** The tip-label tilt (radians) in a vertical orientation, from the user's "Tip label angle" setting: VERTICAL
     *  (90deg, never crosses neighboring branches -- best for unaligned phylograms with varying tip depths) or ANGLED
     *  (45deg, reads better on aligned phylograms / cladograms / ultrametric trees whose tips line up on a baseline).
     *  The sign follows the orientation so the label always extends AWAY from the tree (down for root-top, up for
     *  root-bottom). */
    private double tipLabelAngle() {
        final Options.TIP_LABEL_DIRECTION dir = effectiveTipLabelDirection();
        if (dir == Options.TIP_LABEL_DIRECTION.HORIZONTAL) {
            return 0.0; // upright labels, centred under/over each tip (see paintTipLabelHorizontal)
        }
        final double base = (dir == Options.TIP_LABEL_DIRECTION.VERTICAL) ? (Math.PI / 2.0) : (Math.PI / 4.0);
        return (getOptions().getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_BOTTOM) ? -base : base;
    }

    /** The concrete tip-label direction, resolving AUTO (fit) from the current tip spacing and longest tip label:
     *  upright when labels fit between tips, else 45°, else 90° (see {@link TreePanelUtil#autoTipLabelDirection}). */
    private Options.TIP_LABEL_DIRECTION effectiveTipLabelDirection() {
        final Options.TIP_LABEL_DIRECTION dir = getOptions().getTipLabelDirection();
        if (dir != Options.TIP_LABEL_DIRECTION.AUTO) {
            return dir;
        }
        // use the value cached by the current calcParametersForPainting pass so a reserve and the paint agree; fall back
        // to a live resolve only when no layout pass has run yet (e.g. a paint before the first layout).
        return (_resolved_auto_tip_dir != null) ? _resolved_auto_tip_dir
                : TreePanelUtil.autoTipLabelDirection(2.0 * getYdistance(), _length_of_longest_text_only);
    }


    /** True for a tip whose label is lined up in the far-right aligned column (aligned-phylogram mode). */
    private boolean isAlignedTipLabel(final PhylogenyNode node) {
        return (getControlPanel().getTreeDisplayType() == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM)
                && (node.isExternal() || node.isCollapse());
    }

    /** The logical X of the far-right aligned label column (aligned-phylogram mode), where all tip labels line up. */
    private float alignedLabelColumnX() {
        return (float) ((getMaxDistanceToRoot() * getXcorrectionFactor()) + TreePanel.MOVE + getXdistance());
    }

    /** The logical X where a node's label text begins: the aligned column for an aligned tip, else just right of the
     *  node. Serves as the vertical-mode tilt pivot AND the aligned leader's far endpoint, so the tilted label sits
     *  on the end of its (vertical) leader instead of being pushed diagonally off it. Matches the label anchors in
     *  {@link #paintNodeData} / {@link #paintTaxonomy}. */
    private float labelTextStartX(final PhylogenyNode node) {
        final int half_box = getOptions().getDefaultNodeShapeSize() / 2;
        // clustergram "labels below columns": tip/collapsed labels are drawn past the tip-aligned columns (aligned at
        // the far edge), so the dendrogram sits directly on the grid and the sample labels run along the bottom
        if (tipLabelsBelowColumns() && (node.isExternal() || node.isCollapse())) {
            return labelSegmentStartX(annotationColumnsEndX(), half_box, 0);
        }
        return labelSegmentStartX(isAlignedTipLabel(node) ? alignedLabelColumnX() : node.getXcoord(), half_box, 0);
    }

    /** Whether a branch-length value should be drawn for {@code node} -- the same gate {@link #paintNodeData} used
     *  inline, factored out so the vertical-orientation path (which draws it separately) shares one condition. */
    private boolean shouldWriteBranchLength(final PhylogenyNode node) {
        return getControlPanel().isWriteBranchLengthValues()
                && ((getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR)
                        || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED)
                        || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE))
                && !node.isRoot() && (node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT);
    }

    /** The device (screen) position of a node's logical coords under the current orientation -- the forward of
     *  {@link #toLogicalPoint}. Returns the raw coords in the horizontal orientation. (Used by tests and any code
     *  needing a node's on-screen position in a vertical orientation.) */
    Point2D.Double screenPointFor(final PhylogenyNode node) {
        return screenPoint(node.getXcoord(), node.getYcoord());
    }

    /** The device (screen) position of an arbitrary LOGICAL point under the current orientation. Raw in horizontal. */
    Point2D.Double screenPoint(final double logical_x, final double logical_y) {
        final Point2D.Double p = new Point2D.Double(logical_x, logical_y);
        if (isVerticalOrientation() && (_orientation_R != null)) {
            _orientation_R.transform(p, p);
        }
        return p;
    }

    /** Maps a device (mouse) point back to LOGICAL node-coordinate space via R inverse, for hit-testing in a
     *  vertical orientation. A pass-through (returns the point as-is) in the horizontal orientation. */
    private Point2D.Double toLogicalPoint(final int x, final int y) {
        if (isVerticalOrientation() && (_orientation_R_inverse != null)) {
            final Point2D.Double p = new Point2D.Double(x, y);
            _orientation_R_inverse.transform(p, p);
            return p;
        }
        return new Point2D.Double(x, y);
    }

    /** The viewport's visible rectangle expressed in LOGICAL (un-rotated) node-coordinate space, so the coord-based
     *  screen-culls (which compare against node getXcoord()/getYcoord()) stay correct in EVERY orientation. In a
     *  vertical orientation the device visible rect's corners are mapped back through R-inverse and their bounding box
     *  returned (a 90-degree rotation keeps the box axis-aligned); a pass-through of getVisibleRect() in horizontal. */
    private Rectangle2D.Double logicalVisibleRect() {
        final Rectangle v = getVisibleRect();
        if ((v.width <= 0) || (v.height <= 0)) {
            return null; // no meaningful viewport (offscreen render / File>Print / not-yet-realized) -> never cull
        }
        if (!isVerticalOrientation() || (_orientation_R_inverse == null)) {
            return new Rectangle2D.Double(v.x, v.y, v.width, v.height);
        }
        final Point2D.Double[] corners = { toLogicalPoint(v.x, v.y), toLogicalPoint(v.x + v.width, v.y),
                toLogicalPoint(v.x, v.y + v.height), toLogicalPoint(v.x + v.width, v.y + v.height) };
        double minx = corners[0].x, maxx = corners[0].x, miny = corners[0].y, maxy = corners[0].y;
        for (final Point2D.Double p : corners) {
            minx = Math.min(minx, p.x);
            maxx = Math.max(maxx, p.x);
            miny = Math.min(miny, p.y);
            maxy = Math.max(maxy, p.y);
        }
        return new Rectangle2D.Double(minx, miny, maxx - minx, maxy - miny);
    }

    /** Draws node TEXT so it does NOT ride the tree-rotation transform R (which would render it sideways): it
     *  re-anchors g to the upright base frame at the node's SCREEN position, rotates by {@code angle} (45deg for
     *  tip labels, 0 for upright internal text), and translates back so the existing coord-based label painters
     *  draw correctly with zero internal changes. A no-op (runs the painter directly) in horizontal orientation. */
    private void withNodeTextFrame(final Graphics2D g, final PhylogenyNode node, final double angle,
                                   final Runnable painter) {
        withNodeTextFrame(g, node.getXcoord(), node.getYcoord(), angle, painter);
    }

    /** As above but pivoting the tilted frame at an explicit logical point (e.g. the aligned label column) rather
     *  than the node -- so an aligned tip label lands on the end of its leader, not offset diagonally from the node. */
    private void withNodeTextFrame(final Graphics2D g, final double pivot_x, final double pivot_y, final double angle,
                                   final Runnable painter) {
        if (!isVerticalOrientation() || (_orientation_R == null)) {
            painter.run();
            return;
        }
        final AffineTransform saved = g.getTransform();
        final Point2D.Double screen = new Point2D.Double(pivot_x, pivot_y);
        _orientation_R.transform(screen, screen);
        g.setTransform(_orientation_base_transform);
        g.translate(screen.x, screen.y);
        g.rotate(angle);
        g.translate(-pivot_x, -pivot_y);
        painter.run();
        g.setTransform(saved);
    }

    final void selectNode(final PhylogenyNode node) {
        if ((getFoundNodes0() != null) && getFoundNodes0().contains(node.getId())) {
            getFoundNodes0().remove(node.getId());
        } else {
            ensureFoundNodes0Visible();
            getFoundNodes0().add(node.getId());
        }
        refreshFoundNodes0Bookkeeping();
    }

    /**
     * Adds or removes an entire subtree's external tips from the manual selection (found set 0), reusing the
     * "Select Node(s)" machinery -- for branch-click clade selection. All-or-nothing toggle: if every tip of
     * {@code node}'s subtree is already selected, deselect them all; otherwise select the ones still missing.
     */
    final void selectSubtreeTips(final PhylogenyNode node) {
        if (node == null) {
            return;
        }
        final java.util.List<PhylogenyNode> tips = node.getAllExternalDescendants(); // the node itself if a leaf
        if ((tips == null) || tips.isEmpty()) {
            return;
        }
        if (allTipsSelected(tips)) { // all-or-nothing: fully selected -> clear it, otherwise add the missing tips
            for (final PhylogenyNode t : tips) {
                getFoundNodes0().remove(t.getId());
            }
        } else {
            ensureFoundNodes0Visible();
            for (final PhylogenyNode t : tips) {
                getFoundNodes0().add(t.getId());
            }
        }
        refreshFoundNodes0Bookkeeping();
    }

    /** Makes the found-set-0 search UI (count label + reset button) visible and ensures the set exists. */
    private void ensureFoundNodes0Visible() {
        getControlPanel().getSearchFoundCountsLabel0().setVisible(true);
        getControlPanel().getSearchResetButton0().setEnabled(true);
        getControlPanel().getSearchResetButton0().setVisible(true);
        if (getFoundNodes0() == null) {
            setFoundNodes0(new HashSet<Long>());
        }
    }

    /** After changing found set 0: refresh its count label, and reset the search-0 state when it became empty
     *  (searchReset0 only nulls the set, so the "Found: N" label must be set to 0 here first). */
    private void refreshFoundNodes0Bookkeeping() {
        final Set<Long> found = getFoundNodes0();
        final int n = (found == null) ? 0 : found.size();
        getControlPanel().setSearchFoundCountsOnLabel0(n);
        _search_hit_index = -1; // a click that changed the manual selection restarts the step-through
        if (n < 1) {
            getControlPanel().searchReset0(); // nulls the set -> setFoundNodes0 refreshes the navigator
        }
        else {
            refreshSearchHitNavigation(); // in-place add/remove bypasses setFoundNodes0, so refresh here
        }
    }

    final void setArrowCursor() {
        setCursor(ARROW_CURSOR);
        repaint();
    }

    final void setControlPanel(final ControlPanel atv_control) {
        _control_panel = atv_control;
    }

    final void setFoundNodes0(final Set<Long> found_nodes) {
        // Only a CHANGED hit set restarts the step-through: re-running the same search (e.g. the search box's
        // keyReleased re-fires while the user presses Cmd-G) must NOT snap the position back to the first hit.
        if (!Objects.equals(_found_nodes_0, found_nodes)) {
            _search_hit_index = -1;
        }
        _found_nodes_0 = found_nodes;
        refreshSearchHitNavigation(); // single chokepoint for every found-set-0 REPLACEMENT (search, reset, undo, ...)
    }

    final void setFoundNodes1(final Set<Long> found_nodes) {
        if (!Objects.equals(_found_nodes_1, found_nodes)) {
            _search_hit_index = -1;
        }
        _found_nodes_1 = found_nodes;
        refreshSearchHitNavigation();
    }

    /** Refreshes the left panel's "k / N" step-through navigator; null-safe during construction/teardown. */
    private void refreshSearchHitNavigation() {
        final ControlPanel cp = getControlPanel();
        if (cp != null) {
            cp.updateSearchHitNavigation();
        }
    }

    final void setInOvRect(final boolean in_ov_rect) {
        _in_ov_rect = in_ov_rect;
    }

    final void setLastMouseDragPointX(final float x) {
        _last_drag_point_x = x;
    }

    final void setLastMouseDragPointY(final float y) {
        _last_drag_point_y = y;
    }

    final void setNodeInPreorderToNull() {
        _nodes_in_preorder = null;
        // a tree-structure change (navigation, collapse, delete, ...) invalidates a branch-hover preview whose
        // node may now be detached/relaid-out; this is the shared chokepoint for those changes
        _hover_node = null;
        _click_suppressed = null;
    }

    final void setOvOn(final boolean ov_on) {
        _ov_on = ov_on;
    }

    final void setPhylogenyGraphicsType(final PHYLOGENY_GRAPHICS_TYPE graphics_type) {
        _graphics_type = graphics_type;
        setTextAntialias();
    }

    final void setStartingAngle(final double starting_angle) {
        _urt_starting_angle = starting_angle;
    }

    void setStatisticsForExpressionValues(final DescriptiveStatistics statistics_for_expression_values) {
        _statistics_for_vector_data = statistics_for_expression_values;
    }

    final void setTextAntialias() {
        if ((_phylogeny != null) && !_phylogeny.isEmpty()) {
            if (_phylogeny.getNumberOfExternalNodes() <= LIMIT_FOR_HQ_RENDERING) {
                _rendering_hints.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
            } else {
                _rendering_hints.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED);
            }
        }
        // Screen antialiasing is always on: it is free on modern hardware and essential for legible
        // labels (the old on/off toggle was a relic of when antialiasing slowed redraws too much).
        _rendering_hints.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        _rendering_hints.put(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_LCD_HRGB);
    }

    final void setTreeFile(final File treefile) {
        _treefile = treefile;
    }

    final void setXcorrectionFactor(final float f) {
        _x_correction_factor = f;
    }

    final void setXdistance(final float x) {
        _x_distance = x;
    }

    final void setYdistance(final float y) {
        _y_distance = y;
    }

    final void sortDescendants(final PhylogenyNode node) {
        if (!node.isExternal()) {
            pushUndoCheckpoint("Sort Descendants");
            DESCENDANT_SORT_PRIORITY pri = DESCENDANT_SORT_PRIORITY.NODE_NAME;
            if (getControlPanel().isShowTaxonomyScientificNames() || getControlPanel().isShowTaxonomyCode()) {
                pri = DESCENDANT_SORT_PRIORITY.TAXONOMY;
            } else if (getControlPanel().isShowSeqNames() || getControlPanel().isShowSeqSymbols()
                    || getControlPanel().isShowGeneNames()) {
                pri = DESCENDANT_SORT_PRIORITY.SEQUENCE;
            }
            PhylogenyMethods.sortNodeDescendents(node, pri);
            setNodeInPreorderToNull();
            _phylogeny.externalNodesHaveChanged();
            _phylogeny.clearHashIdToNodeMap();
            _phylogeny.recalculateNumberOfExternalDescendants(true);
            resetNodeIdToDistToLeafMap();
            setEdited(true);
        }
        repaint();
    }

    final void subTree(final PhylogenyNode clicked_node) {
        if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED) {
            JOptionPane.showMessageDialog(this,
                    "Cannot get a sub/super tree in unrooted display",
                    "Attempt to get sub/super tree in unrooted display",
                    JOptionPane.WARNING_MESSAGE);
            return;
        }
        // A leaf cannot root a sub-tree, so show the smallest enclosing clade -- the leaf's
        // parent -- instead. (This used to pop up a modal warning, which could freeze the app.)
        final boolean redirected_from_leaf = clicked_node.isExternal();
        final PhylogenyNode node = redirected_from_leaf ? clicked_node.getParent() : clicked_node;
        if (node == null) {
            return; // a single-node tree: nothing to show
        }
        setNodeInPreorderToNull();
        if (!node.isExternal() && !node.isRoot() && (_subtree_index <= (TreePanel.MAX_SUBTREES - 1))) {
            _sub_phylogenies[_subtree_index] = _phylogeny;
            _sub_phylogenies_temp_roots[_subtree_index] = node;
            ++_subtree_index;
            _phylogeny = TreePanelUtil.subTree(node, _phylogeny);
            if (_phylogeny.getRoot().isCollapse()) {
                _phylogeny.getRoot().setCollapse(false);
            }
            _phylogeny.externalNodesHaveChanged();
            _phylogeny.clearHashIdToNodeMap();
            _phylogeny.recalculateNumberOfExternalDescendants(true);
            rebuildPropertyDisplays(); // color+size schemes summarize the visible tree -- recompute for the subtree
            rebuildCladeBands(); // band roots referenced the old tree -- recompute for the subtree
            rebuildAnnotationColumns(); // rescale gradients / regroup categories to the subtree's visible tips
            // undo history is per displayed (sub)tree: clear it at a navigation boundary so an undo can never
            // restore a snapshot that belongs to a different view (which would desync the sub-tree stack)
            _history.clear();
            notifyEditMenu();
            updateSubSuperTreeButton();
            getMainPanel().getControlPanel().search0();
            getMainPanel().getControlPanel().search1();
            getMainPanel().getControlPanel().updateDomainStructureEvaluethresholdDisplay();
        } else if (node.isRoot() && isCurrentTreeIsSubtree() && !redirected_from_leaf) {
            // an explicit click on the displayed (sub-tree) root climbs one branch toward the
            // root -- the same gesture as the "R1" button -- rather than jumping all the way
            // back; a leaf whose parent is that same root does nothing (already shown)
            superTreeOneLevel();
        }
        _main_panel.getControlPanel().showWhole();
        repaint();
    }

    final void superTree() {
        setNodeInPreorderToNull();
        final PhylogenyNode temp_root = _sub_phylogenies_temp_roots[_subtree_index - 1];
        for (final PhylogenyNode n : temp_root.getDescendants()) {
            n.setParent(temp_root);
        }
        _sub_phylogenies[_subtree_index] = null;
        _sub_phylogenies_temp_roots[_subtree_index] = null;
        _phylogeny = _sub_phylogenies[--_subtree_index];
        _phylogeny.externalNodesHaveChanged();
        _phylogeny.clearHashIdToNodeMap();
        _phylogeny.recalculateNumberOfExternalDescendants(true);
        rebuildPropertyDisplays(); // color+size schemes summarize the visible tree -- recompute for the restored tree
        rebuildCladeBands(); // band roots referenced the old (sub)tree -- recompute for the restored tree
        rebuildAnnotationColumns(); // recompute the columns' schemes for the restored tree
        _history.clear(); // navigation boundary -- see subTree()
        notifyEditMenu();
        getMainPanel().getControlPanel().search0();
        getMainPanel().getControlPanel().search1();
        getMainPanel().getControlPanel().updateDomainStructureEvaluethresholdDisplay();
        updateSubSuperTreeButton();
    }

    /**
     * Move the displayed sub-tree up by exactly one branch in the tree's topology: show the
     * sub-tree rooted at the <i>parent</i> of the current sub-tree's root. This differs from
     * {@link #superTree()}, which pops one frame off the navigation stack -- when the user
     * descended by clicking a leaf (a single stack frame that can span many topological levels),
     * popping the frame jumps far, whereas this always climbs one branch. A no-op on the whole
     * tree. Implemented as "pop the current frame, then re-descend into the parent clade".
     */
    final void superTreeOneLevel() {
        if (!isCurrentTreeIsSubtree()) {
            return;
        }
        final PhylogenyNode current_root = _sub_phylogenies_temp_roots[_subtree_index - 1];
        final PhylogenyNode parent = current_root.getParent();
        superTree(); // back to the phylogeny we descended from (restores current_root's children)
        if ((parent != null) && !parent.isRoot()) {
            subTree(parent); // re-descend one branch up; subTree() fits and repaints
        } else {
            // the parent is the root of that phylogeny, which is now exactly what is displayed
            getMainPanel().getControlPanel().showWhole();
            repaint();
        }
    }

    final void orderSubtree(final PhylogenyNode node) {
        if (node.isExternal()) {
            return;
        }
        DESCENDANT_SORT_PRIORITY pri = DESCENDANT_SORT_PRIORITY.NODE_NAME;
        if (getControlPanel().isShowTaxonomyScientificNames() || getControlPanel().isShowTaxonomyCode()) {
            pri = DESCENDANT_SORT_PRIORITY.TAXONOMY;
        } else if (getControlPanel().isShowSeqNames() || getControlPanel().isShowSeqSymbols()
                || getControlPanel().isShowGeneNames()) {
            pri = DESCENDANT_SORT_PRIORITY.SEQUENCE;
        }
        pushUndoCheckpoint("Order Subtree");
        PhylogenyMethods.orderAppearanceX(node, true, pri);
        setNodeInPreorderToNull();
        getPhylogeny().externalNodesHaveChanged();
        getPhylogeny().clearHashIdToNodeMap();
        getPhylogeny().recalculateNumberOfExternalDescendants(true);
        resetNodeIdToDistToLeafMap();
        setEdited(true);
        getControlPanel().displayedPhylogenyMightHaveChanged(true);
        repaint();
    }

    final void swap(final PhylogenyNode node) {
        if (node.isExternal() || (node.getNumberOfDescendants() < 2)) {
            return;
        }
        if (node.getNumberOfDescendants() > 2) {
            JOptionPane.showMessageDialog(this,
                    "Cannot swap descendants of nodes with more than 2 descendants",
                    "Cannot swap descendants",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        if (!node.isExternal()) {
            pushUndoCheckpoint("Swap Descendants");
            node.swapChildren();
            setNodeInPreorderToNull();
            _phylogeny.externalNodesHaveChanged();
            _phylogeny.clearHashIdToNodeMap();
            _phylogeny.recalculateNumberOfExternalDescendants(true);
            resetNodeIdToDistToLeafMap();
            setEdited(true);
        }
        repaint();
    }


    final void updateOvSettings() {
        switch (getOptions().getOvPlacement()) {
            case LOWER_LEFT:
                setOvXPosition(OV_BORDER);
                setOvYPosition(ForesterUtil.roundToInt(getVisibleRect().height - OV_BORDER - getOvMaxHeight()));
                setOvYStart(ForesterUtil.roundToInt(getOvYPosition() + (getOvMaxHeight() / 2)));
                break;
            case LOWER_RIGHT:
                setOvXPosition(ForesterUtil.roundToInt(getVisibleRect().width - OV_BORDER - getOvMaxWidth()));
                setOvYPosition(ForesterUtil.roundToInt(getVisibleRect().height - OV_BORDER - getOvMaxHeight()));
                setOvYStart(ForesterUtil.roundToInt(getOvYPosition() + (getOvMaxHeight() / 2)));
                break;
            case UPPER_RIGHT:
                setOvXPosition(ForesterUtil.roundToInt(getVisibleRect().width - OV_BORDER - getOvMaxWidth()));
                setOvYPosition(OV_BORDER);
                setOvYStart(ForesterUtil.roundToInt(OV_BORDER + (getOvMaxHeight() / 2)));
                break;
            default:
                setOvXPosition(OV_BORDER);
                setOvYPosition(OV_BORDER);
                setOvYStart(ForesterUtil.roundToInt(OV_BORDER + (getOvMaxHeight() / 2)));
                break;
        }
    }

    final void updateOvSizes() {
        if ((getWidth() > (1.05 * getVisibleRect().width))
                || (getHeight() > (1.05 * getVisibleRect().height))) {
            setOvOn(true);
            float l = getLongestExtNodeInfo();
            final float w_ratio = getOvMaxWidth() / getWidth();
            l *= w_ratio;
            final int ext_nodes = _phylogeny.getRoot().getNumberOfExternalNodes();
            setOvYDistance(getOvMaxHeight() / (2 * ext_nodes));
            float ov_xdist = 0;
            if (!isNonLinedUpCladogram()) {
                ov_xdist = ((getOvMaxWidth() - l) / (ext_nodes));
            } else {
                ov_xdist = ((getOvMaxWidth() - l) / (PhylogenyMethods.calculateMaxDepth(_phylogeny)));
            }
            float ydist = (float) ((getOvMaxWidth() / (ext_nodes * 2.0)));
            if (ov_xdist < 0.0) {
                ov_xdist = 0.0f;
            }
            if (ydist < 0.0) {
                ydist = 0.0f;
            }
            setOvXDistance(ov_xdist);
            final double height = _phylogeny.calculateHeight(!_options.isCollapsedWithAverageHeigh());
            if (height > 0) {
                final float ov_corr = (float) (((getOvMaxWidth() - l) - getOvXDistance()) / height);
                setOvXcorrectionFactor(ov_corr > 0 ? ov_corr : 0);
            } else {
                setOvXcorrectionFactor(0);
            }
        } else {
            setOvOn(false);
        }
    }

    void updateSetOfCollapsedExternalNodes() {
        final Phylogeny phy = getPhylogeny();
        _collapsed_external_nodeid_set.clear();
        if (phy != null) {
            E:
            for (final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
                final PhylogenyNode ext_node = it.next();
                PhylogenyNode n = ext_node;
                while (!n.isRoot()) {
                    if (n.isCollapse()) {
                        _collapsed_external_nodeid_set.add(ext_node.getId());
                        ext_node.setCollapse(true);
                        continue E;
                    }
                    n = n.getParent();
                }
            }
        }
    }

    final void updateSubSuperTreeButton() {
        if (_subtree_index < 1) {
            getControlPanel().deactivateButtonsToReturnToSuperTree();
        } else {
            getControlPanel().activateButtonsToReturnToSuperTree();
        }
    }

    final void updateButtonToUncollapseAll() {
        if (PhylogenyMethods.isHasCollapsedNodes(_phylogeny)) {
            getControlPanel().activateButtonToUncollapseAll();
        } else {
            getControlPanel().deactivateButtonToUncollapseAll();
        }
    }

    final void zoomInDomainStructure() {
        if (_domain_structure_width < 2000) {
            _domain_structure_width *= 1.2;
        }
    }

    final void zoomOutDomainStructure() {
        if (_domain_structure_width > 20) {
            _domain_structure_width *= 0.8;
        }
    }

    private final static void colorizeNodesHelper(final Color c, final PhylogenyNode node) {
        if (node.getNodeData().getNodeVisualData() == null) {
            node.getNodeData().setNodeVisualData(new NodeVisualData());
        }
        node.getNodeData().getNodeVisualData().setFontColor(new Color(c.getRed(), c.getGreen(), c.getBlue()));
    }

    final private static void drawString(final String str, final float x, final float y, final Graphics2D g) {
        g.drawString(str, x, y);
    }

    final private static boolean plusPressed(final int key_code) {
        return ((key_code == KeyEvent.VK_ADD) || (key_code == KeyEvent.VK_PLUS)
                || (key_code == KeyEvent.VK_EQUALS) || (key_code == KeyEvent.VK_SEMICOLON)
                || (key_code == KeyEvent.VK_1));
    }

}
