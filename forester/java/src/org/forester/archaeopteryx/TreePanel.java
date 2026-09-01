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
import java.awt.geom.AffineTransform;
import java.awt.geom.Arc2D;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
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
import java.util.LinkedHashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

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
import javax.swing.UIManager;

import org.forester.archaeopteryx.ControlPanel.NodeClickAction;
import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.io.parsers.json.AuspiceJsonParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyMethods.DESCENDANT_SORT_PRIORITY;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

public final class TreePanel extends JPanel implements ActionListener, MouseWheelListener {

    final private class SubtreeColorizationActionListener implements ActionListener {

        JColorChooser _chooser = null;
        PhylogenyNode _node = null;

        SubtreeColorizationActionListener(final JColorChooser chooser, final PhylogenyNode node) {
            _chooser = chooser;
            _node = node;
        }

        @Override
        public void actionPerformed(final ActionEvent e) {
            final Color c = _chooser.getColor();
            if (c != null) {
                colorizeSubtree(c, _node);
            }
        }
    }

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
    // Gap (px) between the node shape and the first label glyph. The taxonomy label and the node-data label that
    // follows it MUST use the same value (via labelSegmentStartX) so the node data begins exactly where the
    // taxonomy label ends. When the two differed (taxonomy 3, node data 2), the node data started a pixel inside
    // the taxonomy box and an italic scientific name's right overhang overlapped the following node name.
    private final static int LABEL_GAP_AFTER_NODE_SHAPE = 2;
    // tip images: a fixed-width slot per imaged tip (= size * max-aspect) so the label offset is deterministic
    // (known without waiting for the async image load); an image is scaled to fit the slot, aspect preserved.
    private final static float TIP_IMAGE_MAX_ASPECT      = 1.8f;
    private final static int   TIP_IMAGE_GAP             = 3;
    private final static int   MSA_TRACK_GAP             = 8;                                  // gap before the MSA track
    private final static float MSA_LETTER_MIN_WIDTH      = 7f;                                 // draw residue letters at/above this cell width
    private final static double MSA_MAX_VIEWPORT_FRACTION = 0.6;                               // the MSA window takes at most this share of the viewport (the tree keeps the rest)
    private final static int   MSA_MIN_BAND_PX           = 120;                                // but never narrower than this
    private final static int   MSA_RULER_TICK_LEN        = 4;                                  // column-position ruler tick length
    // The domain-architecture box height tracks the tip-row spacing (getYdistance()), CLAMPED to [MIN, MAX]: it grows
    // as the tree is expanded vertically (responds to Y+/Y-) instead of staying at a fixed size; MIN keeps a very
    // zoomed-out tree's boxes visible, MAX keeps a very zoomed-in tree's boxes as bars, not tall blocks (and << the
    // row pitch 2*yDistance, so no overlap). Tuned with the user against real domain trees.
    private final static int DOMAIN_STRUCTURE_HEIGHT_MIN = 6;
    private final static int DOMAIN_STRUCTURE_HEIGHT_MAX = 16;
    // Radial gap (px) between a tip's label column and the start of its domain architecture in circular/unrooted.
    private final static int DOMAIN_RADIAL_GAP = 4;
    // Gap (px) between a tip and its UPRIGHT (horizontal) label along the DEPTH axis in a vertical orientation
    // (see paintTipLabelHorizontal). depthLabelReserve() must reserve this too, else the outermost tip's label
    // pokes past the depth edge and clips (the top in root-bottom, the bottom in root-top).
    private final static int TIP_LABEL_DEPTH_GAP = 5;
    // Aesthetic breathing room (px) reserved past the tip labels at the FAR depth edge in a vertical orientation, so the
    // outermost upright label doesn't sit flush against the window/canvas edge. The 0.11.55 depthLabelReserve() only
    // guaranteed no-CLIP (margin ~0 = flush, which reads too tight); this backs the labels off the edge a little.
    private final static int TIP_LABEL_DEPTH_EDGE_PAD = 8;
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
    // In a radial (circular/unrooted) layout the tip labels radiate OUTWARD, so their reach is reserved from the
    // radius. Cap that reservation at this fraction of the available radius, so a tree whose labels are long relative
    // to the canvas still draws a real circle/fan (the labels then extend past the canvas edge -- zoom out or use
    // "Shorten Labels" to see them) instead of collapsing the whole tree onto the centre point.
    private final static double RADIAL_LABEL_MAX_RATIO = 0.5;
    // Floor the circular tip-ring at this fraction of the available radius so a wide domain track (reserved fully)
    // can't collapse the ring to nothing; only a pathologically huge domain track then clips.
    private final static double RADIAL_MIN_TREE_RATIO = 0.2;
    // Cap the domain track at this fraction of the canvas HALF-radius in a radial layout: the rectangular default
    // (~0.25*viewport) is ~half the circular radius, far too thick -- this keeps it a legible ring that fits.
    private final static double RADIAL_DOMAIN_MAX_FRACTION = 0.2;
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
    // radial (circular/unrooted) square-canvas side bounds (px): a floor so it never vanishes, a cap so an extreme
    // zoom can't blow up the preferred size, and a fallback for when the viewport isn't measurable yet.
    private static final int MIN_RADIAL_DIAMETER = 80;
    private static final int MAX_RADIAL_DIAMETER = 20000;
    private static final int DEFAULT_RADIAL_DIAMETER = 600;
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
    // "Break Long Branches": a branch longer than this multiple of the median (strictly-positive) branch length is
    // drawn capped, with an axis-break glyph. 8x is conservative -- a clock-like / well-behaved tree has no branch
    // that long, so it shows no breaks; only a genuine outlier (an outgroup, a fast lineage) is truncated.
    final static double LONG_BRANCH_BREAK_MULTIPLIER = 8.0;
    // Test seam: set false to render the CAPPED layout WITHOUT the break glyph, so a test can isolate the glyph as the
    // pixel difference against the same layout (production leaves it true; see BreakLongBranchRenderTest).
    static boolean PAINT_BREAK_GLYPH = true;
    // Break-glyph geometry (tree-coordinate units, so it scales with zoom like the node/support marks).
    private static final float BRANCH_BREAK_GLYPH_HALF_HEIGHT = 5.0f;
    private static final float BRANCH_BREAK_GLYPH_SLANT = 2.0f;
    private static final float BRANCH_BREAK_GLYPH_GAP = 3.0f;
    // Position of the break glyph along the branch (0 = parent end, 1 = node end). Off the midpoint so it clears the
    // support/branch-length values, which are centered on the branch midpoint (both are, in a vertical orientation).
    private static final float BRANCH_BREAK_GLYPH_FRACTION = 0.72f;
    // How far the (faint) scale-grid color is blended from the background toward the branch-length color.
    private static final double SCALE_GRID_BLEND = 0.18;
    private static final int    GEOLOGIC_RING_ALPHA = 70; // translucency of the circular geologic-band annuli
    private static final int    GEOLOGIC_AXIS_EDGE_GAP = 6; // px between the rectangular geologic axis and the edge
    private static final int    CALENDAR_AXIS_EDGE_GAP = 12; // px between the calendar tick ruler and the window edge
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
    // Node-age SPINDLE (an alternative shape to the flat HPD bar): a tapered lens peaking at the point estimate. Its
    // MAX half-thickness (a touch fatter than the bar's half-height so the tapered shape reads distinctly), and the
    // sampling step (px) at which the smooth outline is walked.
    private static final double SPINDLE_HALF_HEIGHT = 5.0;
    private static final double SPINDLE_SAMPLE_STEP = 2.0;
    // Fossil stratigraphic-range (FAD/LAD) bars: a solid-ish sepia bar spanning a TIP's observed first->last
    // appearance datum, with short end-caps so it reads as a bracketed interval (the strap "|--|" range convention).
    // A more OPAQUE earth tone than the translucent HPD bar so the two overlays stay visually distinct on a tree that
    // carries both node-age intervals AND fossil tip ranges.
    private static final int    FOSSIL_BAR_HEIGHT = 5;
    private static final int    FOSSIL_BAR_CAP = 3; // half-length (px) of the end-cap ticks
    private static final Color  FOSSIL_BAR_COLOR = new Color(150, 100, 55, 220);  // opaque-ish sepia
    private static final Color  FOSSIL_BAR_COLOR_BW = new Color(60, 60, 60, 220);  // near-solid gray for B&W export
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
    // PER-TAB view state, held here rather than on the (per-WINDOW) Options/ControlPanel so two tabs can be
    // looked at differently -- one tree root-left with its internal labels on, another a root-top clustergram
    // with them off. The shared controls are re-seeded from whichever tab is current (ControlPanel.tabChanged),
    // exactly as the P/A/C buttons and the tree style already are.
    // A new tab is SEEDED from the current shared state: the orientation from Options (which is also what
    // persists across restarts), the two Display-Data toggles from the checkboxes as they stand (they have no
    // Options field and so do not persist -- they are a per-figure decluttering choice, not a preference).
    private Options.TREE_ORIENTATION _tree_orientation = Options.TREE_ORIENTATION.ROOT_LEFT;
    private boolean                  _show_internal_data = true;
    private boolean                  _show_external_data = true;
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
    private final Rectangle2D _rectangle = new Rectangle2D.Float();
    private final Path2D.Float _diamond = new Path2D.Float();
    // "Time tree" badge state: the expensive DATED detection is cached per tree; a user-confirmed ULTRAMETRIC tree
    // is remembered by the tree OBJECT it was confirmed for, so replacing the panel's tree (undo/paste/sub-tree)
    // cleanly drops a stale confirmation without any per-path reset.
    private Phylogeny _time_tree_cached_for    = null;
    private boolean   _time_tree_dated         = false;
    private String    _time_tree_unit          = null;
    private Phylogeny _confirmed_time_tree_for = null;
    // TIME AXIS -- per-tree (per-panel) state. The TYPE + the two refinement toggles are panel-scoped; the type is null
    // == "follow auto-derive" (AptxUtil.deriveTimeAxisType, cached below by tree identity). A saved aptx:time_axis
    // property (TimeAxisConfig) restored on load sets these explicitly and wins over auto-derive.
    private Options.TIME_AXIS_TYPE _time_axis_type = null; // null == auto-derive
    private boolean   _time_axis_grid           = false;
    private boolean   _time_axis_ages           = false;
    private Options.TIME_AXIS_TYPE _derived_time_axis_type = Options.TIME_AXIS_TYPE.NONE;
    private Phylogeny _derived_time_axis_for    = null;   // identity key for the derived-type cache
    // time-AXIS calibration: an explicit user root-age (Ma) override wins; else derive from the oldest <date>, cached
    private double    _time_axis_root_age       = 0;
    private Phylogeny _time_axis_root_age_for    = null;
    private double    _time_axis_present_date    = 0;
    private Phylogeny _time_axis_present_date_for = null;
    private Phylogeny _time_axis_age_cached_for = null;
    private double    _time_axis_root_age_cache = 0;
    // Auspice/Nextstrain: a per-panel, reversible time<->divergence branch-length display mode. Both metrics are
    // retained on the tree (the <date> values + the nextstrain:div properties), so switching just rewrites the branch
    // lengths from the retained metadata -- a pure display mode (no setEdited / no undo). Default TIME = the loaded default.
    private NEXTSTRAIN_BRANCH_MODE _nextstrain_branch_mode = NEXTSTRAIN_BRANCH_MODE.TIME;
    // applicability (carries both a date AND a nextstrain:div) is a function of the tree's data, so cache it by tree
    // identity (like the derived-time-axis cache); a mode toggle doesn't change it, only a tree replacement does
    private Phylogeny _nextstrain_applicable_for = null;
    private boolean   _nextstrain_applicable     = false;
    private final RenderingHints _rendering_hints = new RenderingHints(RenderingHints.KEY_RENDERING,
            RenderingHints.VALUE_RENDER_DEFAULT);
    private JTextArea _rollover_popup;
    private PhylogenyNode _root;
    private final StringBuilder _sb = new StringBuilder();
    // Set transiently by AptxUtil around a PNG export so paintPhylogeny leaves the background unfilled
    // (transparent). Off for screen and every other export format.
    private boolean _export_transparent_background = false;
    // "Break Long Branches": the model-length cap, cached by tree IDENTITY (like maxNodeDateValue) -- it depends only
    // on the tree's branch-length distribution. <=0 = no positive branch length (capping inactive).
    private Phylogeny _break_cap_for = null;
    private double _break_cap = 0;
    private double _break_capped_height = 0;
    // the RADIAL normalizer: max capped distance-to-root over the tips (root branch EXCLUDED, collapse-unaware) -- the
    // capped analogue of getMaxDistanceToRoot, so a positive root branch does not under-fill the ring (differs from
    // _break_capped_height, which is root-INCLUDED + collapse-aware for the rectangular depth scale)
    private double _break_capped_radial_max = 0;
    // the capped height depends on the display-collapse state too (a collapsed clade is measured only to its root),
    // which changes WITHOUT replacing _phylogeny -- so the cache is also keyed by the collapsed-tip-set size
    private int _break_capped_height_collapse_sig = -1;
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
    private final Phylogeny[] _sub_phylogenies = new Phylogeny[TreePanel.MAX_SUBTREES];
    private final PhylogenyNode[] _sub_phylogenies_temp_roots = new PhylogenyNode[TreePanel.MAX_SUBTREES];
    private int _subtree_index = 0;
    private File _treefile = null;
    private TipImageCache _tip_image_cache = null; // lazily created; loads/caches the tip images (local + URL)
    private int _msa_col_offset = 0;               // first alignment column shown in the MSA track's window
    private javax.swing.JScrollBar _msa_scrollbar = null; // this tab's dedicated horizontal scroller for the MSA window
    /** Per-column conservation over the tips CURRENTLY ON SCREEN. Rebuilt from the navigation choke point
     *  (rebuildPropertyDisplays) -- the same hook the property colour scheme rebuilds on -- never per paint:
     *  measured 25 ms at 1000 tips x 5000 columns and 324 ms at 5000 x 10000, which is invisible once per
     *  navigation and would be intolerable on the hover/scroll repaint path. */
    private MsaConservation _msa_conservation = null;
    private float _urt_factor = 1;
    private float _urt_factor_ov = 1;
    // the radial (circular/unrooted) layout is drawn in a SQUARE canvas of this side (px); it is the single knob for
    // radial size -- set to fit the viewport by showWhole, scaled by radial zoom -- decoupled from the rectangular
    // x/y-distance machinery. 0 = not yet initialised (lazy-fit on first use). See resetPreferredSize/setUpUrtFactor.
    private int _radial_diameter = 0;
    // The ring centre + radius the LAST circular paint used (screen: the padded panel centre + the zoom-diameter
    // radius; export: the export-canvas centre + radius). circularLabelAnchor / the ring hit-test read these so they
    // ALWAYS match the drawn tree (the radius is decoupled from the padded preferred size).
    private int _circular_center_x = 0;
    private int _circular_center_y = 0;
    private int _circular_radius = 0;
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
    private int _length_of_longest_text_only; // longest tip TEXT label (name+taxonomy) only, excluding the domain track
    // AUTO tip-label direction resolved ONCE per calcParametersForPainting pass (null = option is not AUTO). Caching it
    // keeps the breadth RESERVES and the PAINT in agreement within a pass -- both read this instead of re-resolving AUTO
    // against a y-distance that the reserves see stale (pre-setYdistance) but the paint sees fresh, which clipped the
    // outermost upright labels when the fresh spacing flipped the resolved angle. See effectiveTipLabelDirection().
    private Options.TIP_LABEL_DIRECTION _resolved_auto_tip_dir = null;
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
        _tree_orientation = tjp.getOptions().getTreeOrientation(); // the default a new tab inherits
        // ... and likewise the two Display-Data toggles, so a second tree opens looking like the one you are on.
        // They are read from the shared checkboxes rather than from Options because they have no Options field
        // (and so, unlike the orientation, they do NOT persist across restarts -- deliberately: they are a
        // per-figure decluttering choice, not a standing preference).
        if (tjp.getControlPanel() != null) {
            // a BRAND-NEW tab is seeded from the shared widgets (there is no tab of our own to read yet)
            _show_internal_data = tjp.getControlPanel().isShowInternalData();
            _show_external_data = tjp.getControlPanel().isShowExternalData();
        }
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
        // Register for tooltips so getToolTipText(MouseEvent) is consulted -- that is the alignment residue
        // readout. It returns null everywhere except over the alignment, so the canvas shows no other tooltip.
        javax.swing.ToolTipManager.sharedInstance().registerComponent(this);
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
                        rotateRadial(true);
                    }
                } else {
                    for (int i = 0; i < notches; ++i) {
                        rotateRadial(false);
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
            // a radial layout has ONE uniform zoom (the square canvas), so a plain wheel notch scales it ONCE via
            // zoomInX/zoomOutX -- calling zoomInY/zoomOutY too would double it; the rectangular layout zooms both axes.
            final boolean radial = (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                    || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED);
            if (notches < 0) {
                for (int i = 0; i < (-notches); ++i) {
                    getControlPanel().zoomInX(AptxConstants.WHEEL_ZOOM_IN_FACTOR,
                            AptxConstants.WHEEL_ZOOM_IN_X_CORRECTION_FACTOR);
                    if (!radial) {
                        getControlPanel().zoomInY(AptxConstants.WHEEL_ZOOM_IN_FACTOR);
                    }
                    getControlPanel().displayedPhylogenyMightHaveChanged(false);
                }
            } else {
                for (int i = 0; i < notches; ++i) {
                    if (!radial) {
                        getControlPanel().zoomOutY(AptxConstants.WHEEL_ZOOM_OUT_FACTOR);
                    }
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
        if ((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) {
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
        } else if (shows(DisplayOption.USE_STYLE) && (PhylogenyMethods.getBranchColorValue(node) != null)) {
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

    private final int calcDynamicHidingFactor() {
        return (int) (0.5 + (getFontMetricsForLargeDefaultFont().getHeight() / (1.5 * getYdistance())));
    }

    /** Whether the CURRENT layout hides tip labels via "Dyna Hide" (the option is on AND the labels are too dense
     *  to fit) -- used to warn on a graphics export at a small size that some labels were dropped. Reads the live
     *  layout (getYdistance), so it is meaningful only right after a layout/paint. Rectangular family only: the
     *  radial layouts have their own hiding logic, so return false there rather than mis-warn. */
    final boolean labelsDynamicallyHidden() {
        return !isRadialLayout() && (getControlPanel() != null) && shows(DisplayOption.DYNAMICALLY_HIDE_DATA)
                && (calcDynamicHidingFactor() > 1);
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
            double dtp = node.getDistanceToParent();
            if (dtp < 0.0) {
                return 0.0f;
            }
            if (breakLongBranchesActive() && (dtp > breakLongBranchCap())) {
                dtp = breakLongBranchCap(); // "Break Long Branches": draw an outlier capped (glyph in paintBranchRectangular)
            }
            return (float) (getXcorrectionFactor() * dtp);
        } else {
            if ((factor == 0) || isNonLinedUpCladogram()) {
                return getXdistance();
            }
            return getXdistance() * factor;
        }
    }

    final private float calculateOvBranchLengthToParent(final PhylogenyNode node, final int factor) {
        if (getControlPanel().isDrawPhylogram()) {
            double dtp = node.getDistanceToParent();
            if (dtp < 0.0) {
                return 0.0f;
            }
            if (breakLongBranchesActive() && (dtp > breakLongBranchCap())) {
                dtp = breakLongBranchCap(); // keep the overview thumbnail consistent with the capped main view
            }
            return (float) (getOvXcorrectionFactor() * dtp);
        } else {
            if ((factor == 0) || isNonLinedUpCladogram()) {
                return getOvXDistance();
            }
            return getOvXDistance() * factor;
        }
    }

    final boolean isBreakLongBranches() {
        return getOptions().isBreakLongBranches();
    }

    /** Compute + cache the break cap AND the resulting capped tree height, by tree identity (both depend only on the
     *  tree's branch-length distribution). NOTE: keyed by tree IDENTITY, so an in-place branch-length edit (NodeEditPanel,
     *  which mutates _phylogeny without replacing it) leaves this stale until the tree object is swapped (nav / undo /
     *  paste) -- the same accepted cache-invalidation class as {@link #maxNodeDateValue()} and the color-by caches. */
    private void ensureBreakCap() {
        final int collapse_sig = _collapsed_external_nodeid_set.size();
        if ((_break_cap_for != _phylogeny) || (_break_capped_height_collapse_sig != collapse_sig)) {
            _break_cap = TreePanelUtil.longBranchBreakCap(_phylogeny, LONG_BRANCH_BREAK_MULTIPLIER);
            // collapse-aware, matching the collapse-aware calculateHeight the non-break depth scale uses
            _break_capped_height = (_break_cap > 0) ? TreePanelUtil.cappedTreeHeight(_phylogeny, _break_cap,
                    !_options.isCollapsedWithAverageHeigh()) : 0;
            _break_capped_radial_max = (_break_cap > 0)
                    ? TreePanelUtil.cappedMaxDistanceToRoot(_phylogeny, _break_cap) : 0;
            _break_cap_for = _phylogeny;
            _break_capped_height_collapse_sig = collapse_sig;
        }
    }

    /** The model-length cap above which a branch is drawn broken. <=0 when the tree has no positive branch length. */
    final double breakLongBranchCap() {
        ensureBreakCap();
        return _break_cap;
    }

    /** The deepest root-to-tip path length after capping (the drawn depth extent when capping is active). */
    private double breakCappedHeight() {
        ensureBreakCap();
        return _break_capped_height;
    }

    /** The RADIAL normalizer while capping: the largest capped distance-from-root over the tips (root branch EXCLUDED),
     *  so the deepest capped tip lands on the outer ring/diameter -- the capped analogue of {@link #getMaxDistanceToRoot()}
     *  (which the radial layouts use when not capping). Distinct from {@link #breakCappedHeight()}, which is
     *  root-INCLUDED for the rectangular depth scale. See circularRadiusFraction / setUpUrtFactor. */
    private double breakCappedRadialMax() {
        ensureBreakCap();
        return _break_capped_radial_max;
    }

    /** The max distance from the root that is actually DRAWN: the capped tree height while Break Long Branches is
     *  active, else the true {@link #getMaxDistanceToRoot()}. Anchors the aligned tip column / domain lineup to the
     *  drawn tips instead of off-canvas at the uncapped outlier distance -- so "A" (aligned phylogram) is capped too. */
    private double displayedMaxDistanceToRoot() {
        return breakLongBranchesActive() ? breakCappedHeight() : getMaxDistanceToRoot();
    }

    /** The tree HEIGHT used for the DRAWN depth extent (preferred size / overview): the capped height while Break Long
     *  Branches is active -- matching the corr computed in calcParametersForPainting -- else calculateHeight. Without
     *  this the extent would be corr * the UNCAPPED height, ballooning the preferred size (a capped internal branch
     *  makes corr large while the uncapped height stays huge -> a hugely oversized scroll extent -> clipping). */
    private double displayedTreeHeight() {
        double h = breakLongBranchesActive() ? breakCappedHeight()
                : getPhylogeny().calculateHeight(!_options.isCollapsedWithAverageHeigh());
        // a subtree's root branch is drawn as a fixed stub (see displayedRootBranchLength), so its length is NOT part
        // of the depth -- exclude it (BOTH calculateHeight AND breakCappedHeight fold in the (capped) root branch) or
        // the extent over-reserves depth by a stub-sized empty margin
        if (isCurrentTreeIsSubtree()) {
            h -= cappedRootBranchLength();
        }
        return h;
    }

    /** The root's OWN branch length as DRAWN to scale (the "root edge" offset) in the current phylogram layout, else 0.
     *  In a rooted phylogram the root is positioned at {@code MOVE + rootBranch * xcorr} (paintPhylogeny); its callers
     *  apply the break-cap. Returns 0 for a SUBTREE: its inherited branch to the (hidden) former parent is meaningless
     *  in isolation, so it is drawn as a fixed short stub (MOVE + xdist) rather than to scale -- the same behaviour a
     *  normal root_dtp=0 tree already gets, which is why a subtree then lays out (and its domains fit) consistently. */
    private double displayedRootBranchLength() {
        if ((_phylogeny == null) || (getControlPanel() == null) || !getControlPanel().isDrawPhylogram()
                || isCurrentTreeIsSubtree()) {
            return 0;
        }
        final double dtp = _phylogeny.getRoot().getDistanceToParent();
        return dtp > 0 ? dtp : 0;
    }

    /** The root's branch length capped exactly as the DRAWN tree caps it under "Break Long Branches" (an outlier root
     *  branch is truncated to the cap). {@code calculateHeight()} and {@code breakCappedHeight()} both fold this in, so
     *  a subtree (drawn with a fixed root stub, not to scale) subtracts it from those heights. */
    private double cappedRootBranchLength() {
        if (_phylogeny == null) {
            return 0;
        }
        final double dtp = _phylogeny.getRoot().getDistanceToParent();
        if (dtp <= 0) {
            return 0;
        }
        return (breakLongBranchesActive() && (dtp > breakLongBranchCap())) ? breakLongBranchCap() : dtp;
    }

    /** The root's OWN branch length as FOLDED INTO {@link #displayedMaxDistanceToRoot()} in the current phylogram
     *  layout: under Break Long Branches the depth is {@code breakCappedHeight()}, which includes the CAPPED root
     *  branch; otherwise {@code getMaxDistanceToRoot()}, which (recalculateMaxDistanceToRoot) includes the UNCAPPED
     *  root branch only for a FULL tree -- a subtree excludes it there. A tip-aligned column anchored at the root's
     *  Xcoord subtracts THIS (not the drawn offset displayedRootBranchLength) to land on the deepest tip without
     *  double-counting the root edge, in every combination of subtree / full-tree-root-edge / break-long-branches. */
    private double rootBranchInMaxDistance() {
        if ((_phylogeny == null) || (getControlPanel() == null) || !getControlPanel().isDrawPhylogram()) {
            return 0;
        }
        final double dtp = _phylogeny.getRoot().getDistanceToParent();
        if (dtp <= 0) {
            return 0;
        }
        if (breakLongBranchesActive()) {
            return (dtp > breakLongBranchCap()) ? breakLongBranchCap() : dtp; // breakCappedHeight folds in the capped root branch
        }
        return isCurrentTreeIsSubtree() ? 0 : dtp; // recalculateMaxDistanceToRoot adds it only for a full tree
    }

    /** The x of the common right-edge domain-architecture column in a phylogram: the deepest tip (rootX +
     *  root-to-tip distance*xcorr) plus the longest tip label. The root's Xcoord already carries the root-edge
     *  offset and displayedMaxDistanceToRoot() ALSO folds in the root branch, so subtract rootBranchInMaxDistance()
     *  once -- else the column double-counts the root edge and, in a subtree (whose root branch is nonzero) or a tree
     *  with a genuine root edge, shifts right and clips the domains off the near edge. Shared by the draw and the
     *  fit test so the two can never disagree. Valid only after the paint that assigns the root's Xcoord. */
    final float alignedPhylogramDomainColumnX() {
        return (float) (((displayedMaxDistanceToRoot() - rootBranchInMaxDistance()) * getXcorrectionFactor())
                + _length_of_longest_text + _phylogeny.getRoot().getXcoord());
    }

    /** Test hook: the right edge of the widest aligned domain track in a phylogram (column x + the widest drawn
     *  architecture, which fills {@code effectiveDomainStructureWidth * 0.9}); must stay within the preferred width. */
    final float alignedDomainColumnRightEdgeForTest() {
        return alignedPhylogramDomainColumnX() + (float) (effectiveDomainStructureWidth() * 0.9);
    }

    /** Test hook: the depth cache getMaxDistanceToRoot() (includes the root branch only when drawn to scale). */
    final double getMaxDistanceToRootForTest() {
        return getMaxDistanceToRoot();
    }

    /** Test hook: the root's drawn x (the fixed stub start MOVE+xdist for a subtree, else MOVE+rootBranch*xcorr). */
    final float rootXcoordForTest() {
        return _phylogeny.getRoot().getXcoord();
    }

    /** Test hook: the fixed-stub root x (MOVE + xdist) a subtree root should be drawn at. */
    final float stubRootXForTest() {
        return TreePanel.MOVE + getXdistance();
    }

    /** Test hook: the depth height feeding the preferred-width extent (excludes a subtree's stubbed root branch). */
    final double displayedTreeHeightForTest() {
        return displayedTreeHeight();
    }

    /** The distance the numeric scale axis spans from the root (origin_x) to the deepest tip. origin_x (rootX) already
     *  carries the root-edge offset that getMaxDistanceToRoot() also folds in, so subtract rootBranchInMaxDistance()
     *  -- else the axis line/ticks overshoot the deepest tip by the root branch and disagree with the scale grid.
     *  Shared by paintScaleAxis and the fit test so the two can never disagree. */
    final double numericScaleAxisMaxDist() {
        return getMaxDistanceToRoot() - rootBranchInMaxDistance();
    }

    /** Whether the "Break Long Branches" layout COULD apply here, independent of the option being on: a rectangular-
     *  family phylogram (UNALIGNED "P" or ALIGNED "A") with branch lengths. The aligned tip column / domain lineup are
     *  anchored to the capped extent via {@link #displayedMaxDistanceToRoot()}, so aligned caps too; the radial layouts
     *  derive radius from the distance directly (not corr), so they are not capped in v1. Drives the re-fit on toggle
     *  (ON or OFF) so the depth scale is recomputed either way. */
    final boolean breakLongBranchesRelevantToLayout() {
        return (getControlPanel() != null) && getControlPanel().isDrawPhylogram() && !isRadialLayout()
                && isPhyHasBranchLengths();
    }

    /** Whether long branches are actively being capped in the current view: the option is on AND the layout is one
     *  that caps AND there is a positive branch length to reference. Gates the capping in calculateBranchLengthToParent
     *  + the depth scale, the break glyph, and the suppression of the numeric distance scale (a broken branch is not to
     *  scale, so a distance ruler/grid over it would be misleading). */
    final boolean breakLongBranchesActive() {
        return isBreakLongBranches() && breakLongBranchesRelevantToLayout() && (breakLongBranchCap() > 0);
    }

    /** Whether Break Long Branches is capping in the CIRCULAR phylogram (radius encodes distance-from-root = the
     *  capped distance): a branch over the cap is drawn as a shorter radial leg + a break glyph, pulling its clade
     *  inward so the informative rings decompress. See circularRadiusFraction. */
    final boolean breakLongBranchesActiveCircular() {
        return isBreakLongBranches() && isCircularPhylogram() && (breakLongBranchCap() > 0);
    }

    /** Whether Break Long Branches is capping in the UNROOTED layout (each branch's radial spoke length is capped;
     *  the urt factor is derived from the capped height so the informative part fans out). See setUpUrtFactor /
     *  paintUnrooted. */
    final boolean breakLongBranchesActiveUnrooted() {
        return isBreakLongBranches() && (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                && getControlPanel().isDrawPhylogram() && isPhyHasBranchLengths() && (breakLongBranchCap() > 0);
    }

    /** Option-INDEPENDENT: whether a RADIAL layout (circular / unrooted phylogram with branch lengths) is one the
     *  feature caps -- drives the re-fit on toggle (ON to apply, OFF to restore), like breakLongBranchesRelevantToLayout
     *  does for the rectangular family. (Unrooted derives its urt factor in the fit pass, so a re-fit is needed.) */
    final boolean breakLongBranchesRelevantToRadialLayout() {
        return (getControlPanel() != null) && isRadialLayout() && getControlPanel().isDrawPhylogram()
                && isPhyHasBranchLengths();
    }

    /** Draw an axis-break glyph ("//" with a small gap in the branch) centred at (mx, my) on a capped branch drawn
     *  along {@code branch_angle} (0 = a horizontal rectangular segment; the radial layouts pass the spoke angle so the
     *  gap rides the branch and the slashes cross it). */
    private void paintBranchBreakGlyph(final Graphics2D g,
                                       final float mx,
                                       final float my,
                                       final double branch_angle,
                                       final boolean to_graphics_file) {
        if (!PAINT_BREAK_GLYPH) {
            return;
        }
        final Color ink = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        final java.awt.geom.AffineTransform saved_tx = g.getTransform();
        if (branch_angle != 0) {
            g.rotate(branch_angle, mx, my); // align the glyph with the branch direction (radial spoke)
        }
        final float h = BRANCH_BREAK_GLYPH_HALF_HEIGHT;
        final float run = BRANCH_BREAK_GLYPH_SLANT;
        final float gap = BRANCH_BREAK_GLYPH_GAP;
        // erase a short slice of the branch behind the glyph so it reads as a real break -- but NOT on a transparent
        // export, where a background-coloured cover would paint a solid blob onto the intended cut-out.
        if (!(to_graphics_file && _export_transparent_background)) {
            g.setColor(getTreeColorSet().getBackgroundColor());
            g.setStroke(new BasicStroke(3.0f));
            g.draw(new Line2D.Float(mx - gap - run, my, mx + gap + run, my));
        }
        // two parallel "/" strokes in the branch ink (bottom-left to top-right; screen y grows downward)
        g.setColor(ink);
        g.setStroke(new BasicStroke(1.4f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        g.draw(new Line2D.Float(mx - gap - run, my + h, mx - gap + run, my - h));
        g.draw(new Line2D.Float(mx + gap - run, my + h, mx + gap + run, my - h));
        g.setColor(ink);
        g.setStroke(saved_stroke);
        g.setTransform(saved_tx);
    }

    final private void cannotOpenBrowserWarningMessage(final String type_type) {
        JOptionPane.showMessageDialog(this,
                "Cannot launch web browser for " + type_type + " data of this node",
                "Cannot launch web browser",
                JOptionPane.WARNING_MESSAGE);
    }

    final private void colorizeSubtree(final Color c, final PhylogenyNode node) {
        _control_panel.setColorBranches(true);
        if (_control_panel.getUseVisualStylesCb() != null) {
            _control_panel.getUseVisualStylesCb().setSelected(true);
        }
        if (node != null) {
            for (final PreorderTreeIterator it = new PreorderTreeIterator(node); it.hasNext(); ) {
                it.next().getBranchData().setBranchColor(new BranchColor(c));
            }
        }
        repaint();
    }

    /** The "Node Style" click-to action: opens the {@link NodeStyleDialog} for the single clicked node so the user
     *  can change its font (style/size/colour) and node mark (shape/fill/size/colour). */
    private void editNodeStyle(final PhylogenyNode node) {
        final List<PhylogenyNode> targets = new ArrayList<>();
        targets.add(node);
        new NodeStyleDialog(this, targets, true).setVisible(true);
    }

    /**
     * Applies a per-node visual-style edit (from {@link NodeStyleDialog}) to {@code targets}: checkpoints undo,
     * writes only the spec's ticked attributes ({@link NodeStyleEditor}), turns on "Use Visual Styles" so the edit
     * is visible, marks the tree edited and repaints. A tree-data MUTATION (undoable, saved), but -- being a pure
     * visual style -- it does NOT write a provenance sentence to the description (per the "font/color changes don't
     * write provenance" rule). Returns the number of nodes changed.
     */
    final int applyNodeStyleEdit(final List<PhylogenyNode> targets, final NodeStyleEditor.Spec spec) {
        if ((targets == null) || targets.isEmpty() || (spec == null) || spec.isEmpty()) {
            return 0;
        }
        boolean has_target = false;
        for (final PhylogenyNode t : targets) {
            if (t != null) {
                has_target = true;
                break;
            }
        }
        if (!has_target) {
            return 0; // nothing to change -> don't checkpoint (which would clear redo) for a no-op
        }
        pushUndoCheckpoint("Node Style");
        // use the BASE (user-chosen) font, not getLargeFont() -- that is the transient, auto-shrunk DISPLAYED font,
        // so a font edit while zoomed out would otherwise pin the node to a tiny size
        final Font tree_font = getTreeFontSet().getBaseFont();
        final int n = NodeStyleEditor.apply(targets, spec, tree_font.getFamily(), tree_font.getSize());
        if (n > 0) {
            // make the edit visible: "Use Visual Styles" is the render gate (isUseVisualStyles() reads the checkbox
            // when present, else the color-branches flag), so set both
            getControlPanel().setColorBranches(true);
            if (getControlPanel().getUseVisualStylesCb() != null) {
                getControlPanel().getUseVisualStylesCb().setSelected(true);
            }
            setEdited(true);
            repaint();
        }
        return n;
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
        final SubtreeColorizationActionListener al = new SubtreeColorizationActionListener(_color_chooser, node);
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

    /**
     * Everything the panel must recompute after the tree's STRUCTURE changed (a node or subtree deleted, tips
     * pruned): the traversal caches, and every display model that either summarizes the tips or HOLDS NODE
     * REFERENCES of its own.
     * <p>
     * Centralised because it was duplicated across the delete paths and drifted apart: the subtree /
     * return-to-whole-tree paths rebuilt the clade bands, the Tools prune rebuilt the property displays and the
     * annotation columns but not the bands, and the click-to delete rebuilt nothing at all. A band keeps the node
     * its mark spans, so a deleted band root was still painted -- and walking a detached node's external
     * descendants threw, over and over, inside the EDT paint loop.
     */
    void afterTreeStructureChanged() {
        setNodeInPreorderToNull();
        if ((_phylogeny != null) && !_phylogeny.isEmpty()) {
            _phylogeny.externalNodesHaveChanged();
            _phylogeny.clearHashIdToNodeMap();
            _phylogeny.recalculateNumberOfExternalDescendants(true);
        }
        resetNodeIdToDistToLeafMap();
        rebuildPropertyDisplays();
        rebuildAnnotationColumns();
        rebuildCladeBands(); // band roots hold NODE REFERENCES: a deleted one must not survive into the next paint
        setHover(null, false); // the pointer's node may be the one just deleted; the focus glow walks its tips
    }

    /** The delete itself, past the confirmation dialog -- package-visible so the behaviour can be tested without
     *  driving a modal dialog. */
    void deleteNodeOrSubtreeConfirmed(final PhylogenyNode node, final boolean node_only) {
        pushUndoCheckpoint(node_only ? "Delete Node" : "Delete Subtree");
        if (node_only) {
            PhylogenyMethods.removeNode(node, _phylogeny);
        } else {
            _phylogeny.deleteSubtree(node, true);
        }
        afterTreeStructureChanged();
        setEdited(true);
        repaint();
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
        deleteNodeOrSubtreeConfirmed(node, node_only);
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

    /** A diamond (rhombus) inscribed in the box [x,x+width] x [y,y+height]: vertices at the four edge midpoints. */
    private void setDiamond(final double x, final double y, final double width, final double heigth) {
        final float cx = (float) (x + (width / 2.0));
        final float cy = (float) (y + (heigth / 2.0));
        _diamond.reset();
        _diamond.moveTo(cx, (float) y);
        _diamond.lineTo((float) (x + width), cy);
        _diamond.lineTo(cx, (float) (y + heigth));
        _diamond.lineTo((float) x, cy);
        _diamond.closePath();
    }

    final private void drawDiamondFilled(final double x,
                                         final double y,
                                         final double width,
                                         final double heigth,
                                         final Graphics2D g) {
        setDiamond(x, y, width, heigth);
        g.fill(_diamond);
    }

    final private void drawDiamondGradient(final float x,
                                           final float y,
                                           final float width,
                                           final float heigth,
                                           final Graphics2D g,
                                           final Color color_1,
                                           final Color color_2,
                                           final Color color_border) {
        setDiamond(x, y, width, heigth);
        g.setPaint(new GradientPaint(x, y, color_1, (x + width), (y + heigth), color_2, false));
        g.fill(_diamond);
        if (color_border != null) {
            g.setPaint(color_border);
            g.draw(_diamond);
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
            case NODE_STYLE:
                editNodeStyle(node);
                break;
            case OPEN_SEQ_WEB:
                openSeqWeb(node);
                break;
            case OPEN_TAX_WEB:
                openTaxWeb(node);
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

    /** Tree-validated breakdown of the highlighted nodes for the "Found / Selected" menu-bar counter tooltip:
     *  {@code [inA, inB, inBoth]}, counted over the nodes actually present in the current (sub)tree -- so it stays
     *  consistent with {@link #getSearchHitCount()} (which is inA + inB - inBoth) even after a prune leaves stale ids
     *  in the found sets. */
    final int[] foundSelectedBreakdown() {
        int a = 0, b = 0, both = 0;
        if ((_phylogeny != null) && hasFoundNodes()) {
            for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                final boolean in0 = isInFoundNodes0(n);
                final boolean in1 = isInFoundNodes1(n);
                if (in0) {
                    a++;
                }
                if (in1) {
                    b++;
                }
                if (in0 && in1) {
                    both++;
                }
            }
        }
        return new int[] { a, b, both };
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

    /** Whether a radial LABEL anchored at (x,y) is off-screen (device coords, 20px slop). Takes the point rather than
     *  the node so the aligned circular phylogram can cull on the tip's LABEL position (the outer ring), not the node's
     *  branch-length position -- else a short-branch tip whose node is off-screen but whose ring label is visible would
     *  be wrongly culled (and its leader would dangle to nothing). */
    private boolean isRadialLabelPointInvisible(final double x, final double y) {
        return ((y < (getVisibleRect().getMinY() - 20)) || (y > (getVisibleRect().getMaxY() + 20))
                || (x < (getVisibleRect().getMinX() - 20)) || (x > (getVisibleRect().getMaxX() + 20)));
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
                if (isRadialLayout()) {
                    // "W" is the node-label-direction flip ("L") in a radial layout -- keep Alt+W matching the button.
                    getControlPanel().toggleNodeLabelDirection();
                } else if (isVerticalOrientation()) {
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
                // a radial layout has ONE uniform zoom, so scale it ONCE (zoomOutX) -- calling zoomOutY too would
                // double it (same reason the wheel plain-zoom is gated); the rectangular layout zooms both axes
                if (!isRadialLayout()) {
                    getMainPanel().getControlPanel().zoomOutY(AptxConstants.WHEEL_ZOOM_OUT_FACTOR);
                }
                getMainPanel().getControlPanel().zoomOutX(AptxConstants.WHEEL_ZOOM_OUT_FACTOR,
                        AptxConstants.WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR);
                getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged(false);
                handled = true;
            } else if (plusPressed(e.getKeyCode())) {
                getMainPanel().getControlPanel().zoomInX(AptxConstants.WHEEL_ZOOM_IN_FACTOR,
                        AptxConstants.WHEEL_ZOOM_IN_FACTOR);
                if (!isRadialLayout()) {
                    getMainPanel().getControlPanel().zoomInY(AptxConstants.WHEEL_ZOOM_IN_FACTOR);
                }
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
                if (isRadialLayout()) {
                    rotateRadial(true);
                    handled = true;
                }
            } else if (e.getKeyCode() == KeyEvent.VK_A) {
                if (isRadialLayout()) {
                    rotateRadial(false);
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
            if (title.equals(ClickToOption.OPEN_SEQ_WEB.title())) {
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
            } else if (title.equals(ClickToOption.OPEN_TAX_WEB.title())) {
                _node_popup_menu_items[i].setEnabled(isCanOpenTaxWeb(node));
            } else if (title.equals(ClickToOption.DELETE_SUBTREE_OR_NODE.title())) {
                if (!getOptions().isEditable()) {
                    continue;
                }
                _node_popup_menu_items[i].setEnabled(isCanDelete());
            } else if (title.equals(ClickToOption.CUT_SUBTREE.title())) {
                if (!getOptions().isEditable()) {
                    continue;
                }
                _node_popup_menu_items[i].setEnabled(isCanCut(node));
            } else if (title.equals(ClickToOption.COPY_SUBTREE.title())) {
                if (!getOptions().isEditable()) {
                    continue;
                }
                _node_popup_menu_items[i].setEnabled(isCanCopy());
            } else if (title.equals(ClickToOption.PASTE_SUBTREE.title())) {
                if (!getOptions().isEditable()) {
                    continue;
                }
                _node_popup_menu_items[i].setEnabled(isCanPaste());
            } else if (title.equals(ClickToOption.EDIT_NODE_DATA.title())) {
                if (!getOptions().isEditable()) {
                    continue;
                }
            } else if (title.equals(ClickToOption.ADD_NEW_NODE.title())) {
                if (!getOptions().isEditable()) {
                    continue;
                }
            } else if (title.equals(ClickToOption.REROOT.title())) {
                _node_popup_menu_items[i].setEnabled(isCanReroot());
            } else if (title.equals(ClickToOption.COLLAPSE_UNCOLLAPSE.title())) {
                _node_popup_menu_items[i].setEnabled((isCanCollapse() && !node.isExternal()));
            } else if (title.equals(ClickToOption.COLOR_SUBTREE.title())) {
                _node_popup_menu_items[i].setEnabled(isCanColorSubtree());
            } else if (title.equals(ClickToOption.SUBTREE.title())) {
                _node_popup_menu_items[i].setEnabled(isCanSubtree(node));
            } else if (title.equals(ClickToOption.SWAP.title())) {
                _node_popup_menu_items[i].setEnabled(node.getNumberOfDescendants() == 2);
            } else if (title.equals(ClickToOption.UNCOLLAPSE_ALL.title())) {
                _node_popup_menu_items[i].setEnabled(isCanUncollapseAll(node));
            }
            _node_popup_menu_items[i].addActionListener(this);
            _node_popup_menu.add(_node_popup_menu_items[i]);
        }
    }

    private final void nodeDataAsSB(final PhylogenyNode node, final StringBuilder sb) {
        if (node != null) {
            if (shows(DisplayOption.SHOW_NODE_NAMES) && (!ForesterUtil.isEmpty(node.getName()))
                    && !nameDuplicatesShownTaxonomy(node)) {
                if (sb.length() > 0) {
                    sb.append(" ");
                }
                // Display-only shortening of an over-long name (e.g. a whole UniProt/NCBI header): the
                // node's actual name is left intact, so export / Find / accession parsing keep the full text.
                sb.append(shows(DisplayOption.SHORTEN_LABELS)
                        ? AptxUtil.shortenLabel(node.getName(), AptxConstants.LONG_NODE_NAME_LIMIT)
                        : node.getName());
            }
            if (node.getNodeData().isHasSequence()
                    && (shows(DisplayOption.SHOW_SEQ_SYMBOLS)
                    || shows(DisplayOption.SHOW_GENE_NAMES)
                    || shows(DisplayOption.SHOW_SEQ_NAMES)
                    || shows(DisplayOption.SHOW_SEQUENCE_ACC)
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
                        if (shows(DisplayOption.SHOW_SEQ_SYMBOLS)
                                && (seq.getSymbol().length() > 0)) {
                            if (sb.length() > 0) {
                                sb.append(" ");
                            }
                            sb.append(seq.getSymbol());
                        }
                        if (shows(DisplayOption.SHOW_GENE_NAMES)
                                && (seq.getGeneName().length() > 0)) {
                            if (sb.length() > 0) {
                                sb.append(" ");
                            }
                            sb.append(seq.getGeneName());
                        }
                        if (shows(DisplayOption.SHOW_SEQ_NAMES)
                                && (seq.getName().length() > 0)) {
                            if (sb.length() > 0) {
                                sb.append(" ");
                            }
                            sb.append(seq.getName());
                        }
                        if (shows(DisplayOption.SHOW_SEQUENCE_ACC)
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
            if (shows(DisplayOption.SHOW_PROPERTIES) && node.getNodeData().isHasProperties()) {
                // may be empty (every field deselected, or only internal metadata) -- do not leave a trailing space
                final String props = propertiesToString(node);
                if (props.length() > 0) {
                    if (sb.length() > 0) {
                        sb.append(" ");
                    }
                    sb.append(props);
                }
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
                v.visit(scientificNameForDisplay(taxonomy), isItalicName(taxonomy));
                v.visit(" (" + taxonomy.getCommonName() + ") ", false);
            } else if (has_sci) {
                v.visit(scientificNameForDisplay(taxonomy), isItalicName(taxonomy));
                v.visit(" ", false);
            } else if (has_common) {
                v.visit(taxonomy.getCommonName() + " ", false);
            }
        } else if (show_sci) {
            if (has_sci) {
                v.visit(scientificNameForDisplay(taxonomy), isItalicName(taxonomy));
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

    /**
     * Whether a taxonomy's scientific name should be drawn in italics -- gated on the display option, then true only
     * for genus/species-level names (see {@link TreePanelUtil#scientificNameIsItalic}), so family/order/higher render
     * upright per the taxonomic convention.
     */
    private boolean isItalicName(final Taxonomy taxonomy) {
        return getOptions().isUseItalicScientificNames()
                && TreePanelUtil.scientificNameIsItalic(taxonomy.getRank(), taxonomy.getScientificName());
    }

    /**
     * True when the node's own name merely repeats a taxonomy label already being drawn (e.g. a node named
     * "Filoviridae" that also carries a {@code <taxonomy>} whose scientific/common name is "Filoviridae"), so the
     * same word is not printed twice. Only suppresses when that taxonomy field is actually shown.
     */
    private boolean nameDuplicatesShownTaxonomy(final PhylogenyNode node) {
        if (!node.getNodeData().isHasTaxonomy()) {
            return false;
        }
        final Taxonomy tax = node.getNodeData().getTaxonomy();
        return TreePanelUtil.nodeNameDuplicatesTaxonomy(node.getName(), tax.getScientificName(), tax.getCommonName(),
                shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES), shows(DisplayOption.SHOW_TAXONOMY_COMMON_NAMES));
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
        final boolean bl_bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
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
            // "Break Long Branches": a capped triangular chord gets the break glyph at its midpoint too, so the branch
            // is never silently shortened without a marker.
            if (breakLongBranchesActive() && (node.getDistanceToParent() > breakLongBranchCap())) {
                paintBranchBreakGlyph(g, x1 + ((x2 - x1) * BRANCH_BREAK_GLYPH_FRACTION),
                        y1 + ((y2 - y1) * BRANCH_BREAK_GLYPH_FRACTION), 0, to_graphics_file);
            }
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
            if (!shows(DisplayOption.WIDTH_BRANCHES) || (PhylogenyMethods.getBranchWidthValue(node) == 1)) {
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
            // "Break Long Branches": mark a capped branch with an axis-break glyph across the middle of its horizontal
            // segment (the segment is drawn shortened; the true length is still shown by the branch-length label).
            if (breakLongBranchesActive() && (node.getDistanceToParent() > breakLongBranchCap())) {
                paintBranchBreakGlyph(g, x1a + ((x2a - x1a) * BRANCH_BREAK_GLYPH_FRACTION), y2, 0, to_graphics_file);
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
            final double r = circularRadiusFraction(n); // depth-based (cladogram) or distance-based (phylogram)
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
        }
        if (n.isCollapse()) {
            // do NOT recurse into a collapsed clade's hidden descendants: hidden INTERNAL nodes have no
            // _urt_nodeid_angle_map entry, so recursing autounboxed a null angle at the theta lookup below -> NPE on a
            // circular tree with the overview shown (the overview counterpart of the 0.11.3 main-path collapse fix).
            // Position the collapsed root at its ring slot from its (main-paint) angle so its branch still draws.
            // NOTE: this stops the crash; the overview tip loop (paintCircularLite) still spaces the VISIBLE tips by
            // the full external count, so a collapsed thumbnail is slightly out of scale -- a cosmetic, screen-only
            // follow-up (the whole *Lite overview path doesn't honor collapse for layout).
            final Double a = _urt_nodeid_angle_map.get(n.getId());
            if (a != null) {
                final double r = circularRadiusFraction(n);
                n.setXSecondary((float) (center_x + (radius * r * Math.cos(a))));
                n.setYSecondary((float) (center_y + (radius * r * Math.sin(a))));
            }
            return;
        }
        final List<PhylogenyNode> descs = n.getDescendants();
        for (final PhylogenyNode desc : descs) {
            paintCircularsLite(desc, phy, center_x, center_y, radius, g);
        }
        final double r = circularRadiusFraction(n);
        final double theta = _urt_nodeid_angle_map.get(n.getId());
        n.setXSecondary((float) (center_x + (radius * r * Math.cos(theta))));
        n.setYSecondary((float) (center_y + (radius * r * Math.sin(theta))));
        for (final PhylogenyNode desc : descs) {
            paintBranchCircularLite(n, desc, g);
        }
    }

    final private void paintCollapsedNode(final Graphics2D g,
                                          final PhylogenyNode node,
                                          final boolean to_graphics_file,
                                          final boolean to_pdf,
                                          final boolean is_in_found_nodes) {
        Color c = null;
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
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
        else if (getOptions().isColorLabelsSameAsParentBranch() && shows(DisplayOption.USE_STYLE)
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
            final boolean conf_bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
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
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
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

    /**
     * Support (confidence) + branch-length NUMBERS for a RADIAL (circular/unrooted) layout: drawn at the incoming
     * branch's midpoint, ROTATED to ride the branch (like the radial labels), centred on the midpoint and offset just
     * off the line -- the support in the confidence colour, then a space, then the length in the branch-length colour.
     * The rectangular layout draws these in paintConfidenceValues / paintBranchLength (which are horizontal-branch
     * only, gated by isShowConfidenceValuesForNode / shouldWriteBranchLength); this is the radial equivalent, with the
     * gate inlined (no layout term -- this method is only called from the radial paths). {@code branch_angle} is the
     * direction of the node's incoming branch; {@code (mid_x,mid_y)} its device midpoint.
     */
    private void paintBranchDataRadial(final Graphics2D g, final PhylogenyNode node, final double mid_x,
                                       final double mid_y, final double branch_angle, final boolean to_pdf,
                                       final boolean to_graphics_file) {
        if (node.isRoot()) {
            return;
        }
        String support = "";
        if (shows(DisplayOption.WRITE_CONFIDENCE_VALUES) && !node.isExternal()
                && node.getBranchData().isHasConfidences()) {
            final List<Confidence> confidences = node.getBranchData().getConfidences();
            Collections.sort(confidences);
            final double min_confidence = getOptions().getMinConfidenceFraction() * _confidence_scale_max;
            support = confidenceLabel(confidences, getOptions().isShowMadConfidence(), min_confidence,
                    getOptions().isShowConfidenceStddev(),
                    getOptions().getNumberOfDigitsAfterCommaForConfidenceValues());
        }
        final String length = (shows(DisplayOption.WRITE_BRANCH_LENGTH_VALUES)
                && (node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT))
                        ? FORMATTER_BRANCH_LENGTH.format(node.getDistanceToParent()) : "";
        if (support.isEmpty() && length.isEmpty()) {
            return;
        }
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = getTreeFontSet().getFontMetricsSmall();
        final int gap_w = (support.isEmpty() || length.isEmpty()) ? 0 : fm.stringWidth(" ");
        final int total_w = fm.stringWidth(support) + gap_w + fm.stringWidth(length);
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
        final boolean found = isInFoundNodes(node);
        double m = branch_angle % TWO_PI;
        if (m < 0) {
            m += TWO_PI;
        }
        boolean flipped = false;
        if ((m > HALF_PI) && (m < ONEHALF_PI)) {
            m -= PI; // keep the numbers upright on the far half of the fan
            flipped = true;
        }
        final AffineTransform saved = g.getTransform();
        g.rotate(m, mid_x, mid_y);
        float x = (float) (mid_x - (total_w / 2.0)); // centre the "support length" on the branch midpoint
        // sit just off the branch line, not on it; flip the perpendicular offset when the frame was flipped so the
        // numbers stay on the SAME visual side of the branch across both halves of the fan
        final float baseline = (float) (mid_y + (flipped ? 2.0 : -2.0));
        if (!support.isEmpty()) {
            g.setColor(dimNonMatch(inkColor(to_pdf, to_graphics_file, getTreeColorSet().getConfidenceColor()), found,
                    bw));
            TreePanel.drawString(support, x, baseline, g);
            x += fm.stringWidth(support) + gap_w;
        }
        if (!length.isEmpty()) {
            g.setColor(dimNonMatch(inkColor(to_pdf, to_graphics_file, getTreeColorSet().getBranchLengthColor()), found,
                    bw));
            TreePanel.drawString(length, x, baseline, g);
        }
        g.setTransform(saved);
    }

    /** Draws a collapse marker for a collapsed clade-root in a radial (CIRCULAR/UNROOTED) layout: a filled triangle
     *  whose apex is AT the node and which opens OUTWARD along {@code out_angle} (the radial/branch direction, away
     *  from the root), plus the hidden-tip count "(N)" riding the branch just beyond it -- the radial analogue of the
     *  rectangular collapsed-clade triangle ({@link #paintCollapsedNode}). The node box is skipped for collapsed nodes
     *  ({@link #paintNodeBox}), so this IS the marker. Fill + found-state coloring mirror paintCollapsedNode: filled in
     *  the found colour when EVERY hidden tip is selected, outlined when only some are. The triangle's size grows with
     *  the tip count on a log scale, like the rectangular triangle. */
    private void paintRadialCollapsedMarker(final Graphics2D g, final PhylogenyNode node, final double out_angle,
                                            final boolean to_pdf, final boolean to_graphics_file) {
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
        final int[] fc = collapsedCladeFoundCounts(node); // {hits set 0, hits set 1, tips hit, total tips}
        Color found_c = null;
        if (!bw && (fc[2] > 0)) {
            found_c = ((fc[0] > 0) && (fc[1] > 0)) ? getTreeColorSet().getFoundColor0and1()
                    : (fc[0] > 0) ? getTreeColorSet().getFoundColor0() : getTreeColorSet().getFoundColor1();
        }
        final boolean all_tips_found = (found_c != null) && (fc[2] == fc[3]);
        final Color fill;
        if (bw) {
            fill = Color.BLACK;
        } else if (all_tips_found) {
            fill = found_c; // the whole collapsed clade is selected -> paint the triangle in the found colour
        } else if (getOptions().isColorLabelsSameAsParentBranch() && shows(DisplayOption.USE_STYLE)
                && (PhylogenyMethods.getBranchColorValue(node) != null)) {
            fill = PhylogenyMethods.getBranchColorValue(node);
        } else if (to_pdf) {
            fill = getTreeColorSet().getBranchColorForPdf();
        } else {
            fill = getTreeColorSet().getCollapseFillColor();
        }
        // triangle: apex AT the node, opening OUTWARD; depth (and base half-width) grow gently with the hidden-tip
        // count on a log scale (as in paintCollapsedNode), with a visible floor and a cap so a huge clade stays sane
        final int tips = Math.max(fc[3], 2);
        final float base = getOptions().getDefaultNodeShapeSize() + 2f;
        double depth = base * (2.2 + Math.log10(tips));
        final double cap = base * 5.0;
        if (depth > cap) {
            depth = cap;
        }
        final double half_w = depth * 0.45;
        final double cos = Math.cos(out_angle), sin = Math.sin(out_angle);
        final double px = node.getXcoord(), py = node.getYcoord();
        final double bx = px + (cos * depth), by = py + (sin * depth); // base centre, outward along the branch
        final double perp_x = -sin, perp_y = cos; // unit perpendicular to the branch
        _polygon.reset();
        _polygon.moveTo((float) px, (float) py); // apex at the node (toward the root)
        _polygon.lineTo((float) (bx + (perp_x * half_w)), (float) (by + (perp_y * half_w)));
        _polygon.lineTo((float) (bx - (perp_x * half_w)), (float) (by - (perp_y * half_w)));
        _polygon.closePath();
        final Color saved = g.getColor();
        // fill respecting the node-fill mode, mirroring paintCollapsedNode: SOLID = fill; NONE = fill the background
        // first (cut out any branches showing through), then outline; GRADIENT = background->fill along the branch
        final NodeVisualData.NodeFill node_fill = getOptions().getDefaultNodeFill();
        if (node_fill == NodeVisualData.NodeFill.NONE) {
            g.setColor(getBackground());
            g.fill(_polygon);
            g.setColor(fill);
            g.draw(_polygon);
        } else if (node_fill == NodeVisualData.NodeFill.GRADIENT) {
            g.setPaint(new GradientPaint((float) px, (float) py, getBackground(), (float) bx, (float) by, fill, false));
            g.fill(_polygon);
            g.setPaint(fill);
            g.draw(_polygon);
        } else {
            g.setColor(fill);
            g.fill(_polygon);
        }
        if ((found_c != null) && !all_tips_found) {
            // only SOME tips selected -> outline in the found colour (partial hint; matches paintCollapsedNode)
            g.setColor(found_c);
            g.draw(_polygon);
        }
        // the "(N)" hidden-tip count, centred on a point just beyond the base and riding the branch (rotated upright
        // on both halves of the fan -- centring on the anchor keeps the position flip-independent)
        final String label = "(" + fc[3] + ")";
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = getTreeFontSet().getFontMetricsSmall();
        final double anchor = depth + 3 + (fm.stringWidth(label) / 2.0); // near text edge clears the base
        final double lx = px + (cos * anchor), ly = py + (sin * anchor);
        double m = out_angle % TWO_PI;
        if (m < 0) {
            m += TWO_PI;
        }
        if ((m > HALF_PI) && (m < ONEHALF_PI)) {
            m -= PI; // keep the count upright on the far half of the fan
        }
        final AffineTransform saved_t = g.getTransform();
        g.rotate(m, lx, ly);
        // dim the count only when the clade holds NO hit -- a hit-containing collapsed clade keeps its count vivid
        // (tip-based, matching the triangle's found fill; the collapsed ROOT is never itself in the found set)
        g.setColor(dimNonMatch(inkColor(to_pdf, to_graphics_file, getTreeColorSet().getSequenceColor()),
                fc[2] > 0, bw));
        TreePanel.drawString(label, (float) (lx - (fm.stringWidth(label) / 2.0)),
                (float) (ly + (fm.getAscent() / 2.0) - (fm.getDescent() / 2.0)), g);
        g.setTransform(saved_t);
        g.setColor(saved);
    }

    /** Internal-node label for a VERTICAL orientation: horizontal, RIGHT-ALIGNED so it ends just LEFT of the branch,
     *  centered at the incoming branch's midpoint -- mirroring the support/length on the right. Draws the taxonomy
     *  (its own italic/colour) then the node-data string. Tips are skipped (their 45deg/90deg label is drawn
     *  elsewhere); a long label extends leftward toward the neighbouring subtree, as accepted. */
    private void paintInternalLabelLeftVertical(final Graphics2D g, final PhylogenyNode node, final boolean to_pdf,
                                                final boolean to_graphics_file) {
        if (node.isRoot() || !isShowInternalDataForThisTab()) {
            return;
        }
        final boolean using_visual_font = setFont(g, node);
        final boolean is_in_found_nodes = isInFoundNodes(node);
        final Taxonomy taxonomy = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy() : null;
        final boolean show_tax = (taxonomy != null)
                && (shows(DisplayOption.SHOW_TAX_CODE) || shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES)
                        || shows(DisplayOption.SHOW_TAXONOMY_COMMON_NAMES) || shows(DisplayOption.SHOW_TAX_RANK))
                && !TreePanelUtil.isDuplicateOfAncestorTaxon(node, this::internalTaxonomyLabelText);
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
        return shows(DisplayOption.WRITE_CONFIDENCE_VALUES) && !node.isExternal() && !node.isRoot()
                && ((getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.ROUNDED)
                        || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR)
                        || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE))
                && node.getBranchData().isHasConfidences();
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
        if ((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) {
            c = Color.BLACK;
        }
        g.setColor(c);
        // Center the symbol on the middle of the branch (parent -> node), not at the node on its right
        // end. The root has no branch, so it falls back to the passed-in node position.
        float cx = x;
        float cy = y;
        final PhylogenyNode parent = node.getParent();
        if (parent != null) {
            final float[] center;
            if (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) {
                // the branch is a radial leg along the node's spoke, NOT a straight line to the parent -> place the
                // symbol on that leg (the ring centre = the root's coords in the circular layout)
                center = TreePanelUtil.circularSupportSymbolCenter(_phylogeny.getRoot().getXcoord(),
                        _phylogeny.getRoot().getYcoord(), node.getXcoord(), node.getYcoord(), parent.getXcoord(),
                        parent.getYcoord());
            }
            else {
                // UNROOTED draws a straight parent->node line, so the Cartesian midpoint IS on the branch;
                // rectangular uses the branch-mid x at the node's y.
                final boolean radial = getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED;
                center = TreePanelUtil.supportSymbolCenter(parent.getXcoord(), node.getXcoord(), parent.getYcoord(),
                        node.getYcoord(), radial);
            }
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
                || (shows(DisplayOption.USE_STYLE) && ((node.getNodeData().getNodeVisualData() != null)
                && ((node.getNodeData().getNodeVisualData().getNodeColor() != null)
                || (node.getNodeData().getNodeVisualData().getSize() != NodeVisualData.DEFAULT_SIZE)
                || (node.getNodeData().getNodeVisualData().getFillType() != NodeFill.DEFAULT)
                || (node.getNodeData().getNodeVisualData().getShape() != NodeShape.DEFAULT))))
                || (shows(DisplayOption.WRITE_EVENTS) && node.isHasAssignedEvent()
                && (node.getNodeData().getEvent().isDuplication()
                || node.getNodeData().getEvent().isSpeciation()
                || node.getNodeData().getEvent().isSpeciationOrDuplication()))) {
            NodeVisualData vis = null;
            if (shows(DisplayOption.USE_STYLE) && (node.getNodeData().getNodeVisualData() != null)
                    && (!node.getNodeData().getNodeVisualData().isEmpty())) {
                vis = node.getNodeData().getNodeVisualData();
            }
            float box_size = getOptions().getDefaultNodeShapeSize();
            if ((vis != null) && (vis.getSize() != NodeVisualData.DEFAULT_SIZE)) {
                box_size = vis.getSize();
            }
            final float half_box_size = box_size / 2.0f;
            Color outline_color = null;
            if ((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) {
                outline_color = Color.BLACK;
            } else if (isInFoundNodes(node)) {
                outline_color = getColorForFoundNode(node);
            } else if (vis != null) {
                if (vis.getNodeColor() != null) {
                    outline_color = vis.getNodeColor();
                } else if (vis.getFontColor() != null) {
                    outline_color = vis.getFontColor();
                }
            } else if (shows(DisplayOption.WRITE_EVENTS) && TreePanelUtil.isHasAssignedEvent(node)) {
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
            if ((vis != null) && (vis.getShape() != NodeShape.DEFAULT)) {
                shape = vis.getShape(); // CIRCLE / RECTANGLE / DIAMOND
            }
            if ((shape == null) && (getOptions().getDefaultNodeShape() != NodeShape.DEFAULT)) {
                shape = getOptions().getDefaultNodeShape();
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
            } else if (shape == NodeShape.DIAMOND) {
                if (fill == NodeVisualData.NodeFill.GRADIENT) {
                    drawDiamondGradient(x - half_box_size,
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
                    drawDiamondGradient(x - half_box_size,
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
                    drawDiamondFilled(x - half_box_size, y - half_box_size, box_size, box_size, g);
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
        if (!isShowInternalDataForThisTab() && !node.isExternal() && !node.isCollapse()) {
            return 0;
        }
        if (!isShowExternalDataForThisTab() && (node.isExternal() || node.isCollapse())) {
            return 0;
        }
        _sb.setLength(0);
        int x = 0;
        if (add > 0) {
            x += add;
        }
        final int half_box_size = effectiveNodeHalfBoxSize(node);
        x += tipImageAdvance(node); // shift this tip's label past its image (0 for a non-imaged/aligned tip)
        final boolean want_dot = (isColorByProperty() && (node.isExternal() || node.isCollapse()))
                || (isSizeByProperty() && node.isExternal());
        // at a pie node the pie IS the marker, so suppress the plain color/size dot (which the later-drawn pie
        // would otherwise sit on top of -- and a large Size-by dot could peek out from behind the smaller pie).
        // Reads the per-node distribution CACHE (built in rebuildAncestralPieColors), so this is an O(1) lookup.
        if (want_dot && !(isShowAncestralPies() && (_ancestral_pie_dist != null)
                && _ancestral_pie_dist.containsKey(node))) {
            drawPropertyColorDot(g, node);
        }
        if (usesAboveBranchInternalLabel(node)) {
            return paintInternalLabelAboveBranch(g, node, is_in_found_nodes, to_pdf, to_graphics_file, half_box_size);
        }

        // Whether a taxonomy label was actually drawn. taxonomyLabel no longer leaves the label in _sb
        // (the old paintTaxonomy did, and this flag used to read _sb.length()), so derive it from the
        // painted advance instead -- the collapsed-node label logic below depends on it.
        boolean saw_species = false;
        if ((shows(DisplayOption.SHOW_TAX_CODE) || shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES)
                || shows(DisplayOption.SHOW_TAXONOMY_COMMON_NAMES) || shows(DisplayOption.SHOW_TAX_RANK))
                && node.getNodeData().isHasTaxonomy()
                && !TreePanelUtil.isDuplicateOfAncestorTaxon(node, this::internalTaxonomyLabelText)) {
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
                        && (shows(DisplayOption.SHOW_TAX_CODE) || shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES)
                        || shows(DisplayOption.SHOW_SEQ_NAMES) || shows(DisplayOption.SHOW_NODE_NAMES))) {
                    // prefer the clade's COMMON taxon (deepest shared among its tips) when taxonomy is shown -- more
                    // informative than the boundary tip names; falls back to first...last when none is derivable
                    final String common = (shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES)
                            || shows(DisplayOption.SHOW_TAX_CODE)) ? collapsedCommonTaxon(node) : "";
                    final PhylogenyNode first = PhylogenyMethods.getFirstExternalNode(node);
                    final PhylogenyNode last = PhylogenyMethods.getLastExternalNode(node);
                    if (!ForesterUtil.isEmpty(common)) {
                        addLabelForCollapsedCommonTaxon(common, node.getAllExternalDescendants().size(), node);
                    } else if (shows(DisplayOption.SHOW_TAX_CODE) && first.getNodeData().isHasTaxonomy()
                            && last.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty(first.getNodeData().getTaxonomy().getTaxonomyCode())
                            && !ForesterUtil.isEmpty(last.getNodeData().getTaxonomy().getTaxonomyCode())) {
                        addLabelForCollapsed(first.getNodeData().getTaxonomy().getTaxonomyCode(),
                                last.getNodeData().getTaxonomy().getTaxonomyCode(),
                                node.getAllExternalDescendants().size(),
                                node);
                    } else if (shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES) && first.getNodeData().isHasTaxonomy()
                            && last.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty(first.getNodeData().getTaxonomy().getScientificName())
                            && !ForesterUtil.isEmpty(last.getNodeData().getTaxonomy().getScientificName())) {
                        addLabelForCollapsed(first.getNodeData().getTaxonomy().getScientificName(),
                                last.getNodeData().getTaxonomy().getScientificName(),
                                node.getAllExternalDescendants().size(),
                                node);
                    } else if (shows(DisplayOption.SHOW_SEQ_NAMES) && first.getNodeData().isHasSequence()
                            && last.getNodeData().isHasSequence()
                            && !ForesterUtil.isEmpty(first.getNodeData().getSequence().getName())
                            && !ForesterUtil.isEmpty(last.getNodeData().getSequence().getName())) {
                        addLabelForCollapsed(first.getNodeData().getSequence().getName(),
                                last.getNodeData().getSequence().getName(),
                                node.getAllExternalDescendants().size(),
                                node);
                    } else if (shows(DisplayOption.SHOW_NODE_NAMES) && !ForesterUtil.isEmpty(first.getName())
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
            pos_x = labelSegmentStartX((float) ((displayedMaxDistanceToRoot() * getXcorrectionFactor()) + TreePanel.MOVE
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
        if (sb_str.length() > 0) {
            TreePanel.drawString(sb_str, pos_x, pos_y, g);
        }
        if (_sb.length() > 0) {
            x += labelStringWidth(g, _sb.toString(), using_visual_font, is_in_found_nodes, to_pdf) + 5;
        }
        return x;
    }

    /**
     * Whether {@code node}'s label should use the publication-style placement: to the LEFT of the node,
     * right-aligned, on top of the incoming branch. Applies to non-root internal nodes (not collapsed)
     * in the rectangular-family layouts when the option is on.
     */
    /** The exact taxonomy text the paint path would draw for {@code node} under the active taxonomy checkboxes
     *  (rank/code/scientific/common), trimmed -- the label basis for {@link TreePanelUtil#isDuplicateOfAncestorTaxon}
     *  so the redundancy test matches what is actually shown. Empty when the node has no taxonomy or none of its
     *  parts are displayed. */
    private String internalTaxonomyLabelText(final PhylogenyNode node) {
        if (!node.getNodeData().isHasTaxonomy()) {
            return "";
        }
        final StringBuilder sb = new StringBuilder();
        nodeTaxonomyDataAsSB(node.getNodeData().getTaxonomy(), sb);
        return sb.toString().trim();
    }

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
        if ((shows(DisplayOption.SHOW_TAX_CODE) || shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES)
                || shows(DisplayOption.SHOW_TAXONOMY_COMMON_NAMES) || shows(DisplayOption.SHOW_TAX_RANK))
                && node.getNodeData().isHasTaxonomy()
                && !TreePanelUtil.isDuplicateOfAncestorTaxon(node, this::internalTaxonomyLabelText)) {
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

    // collapsed-clade common-taxon labels: the deepest taxon shared by a collapsed clade's tips (from their cached /
    // stored lineages). _tip_lineages is the whole-tree tip->lineage map (cache-only, built lazily); both caches are
    // rebuilt whenever _phylogeny is REPLACED (identity change -- covers setTree/subTree/superTree/paste/restore in
    // one place, no per-call-site clearing). Keyed by node identity.
    private Map<PhylogenyNode, org.forester.ws.seqdb.TaxonLineage> _tip_lineages          = null;
    private Phylogeny                                              _tip_lineages_for      = null;
    private final Map<PhylogenyNode, String>                       _collapsed_taxon_cache = new java.util.IdentityHashMap<>();

    /**
     * The clade's COMMON taxon (the deepest taxon shared by ALL its tips, from their cached / stored lineages) for a
     * collapsed-node label, or "" when none is derivable. Cached per node; the whole-tree tip lineages are built
     * lazily and cache-only (no network at paint) and rebuilt when the displayed tree is replaced. A collapse /
     * uncollapse does not change a node's tips, so the per-node cache stays valid across those.
     */
    private String collapsedCommonTaxon(final PhylogenyNode node) {
        if ((_tip_lineages == null) || (_tip_lineages_for != _phylogeny)) {
            _tip_lineages = TreePanelUtil.tipLineages(_phylogeny, TreePanelUtil.getDefaultLineageService());
            _tip_lineages_for = _phylogeny;
            _collapsed_taxon_cache.clear();
        }
        final String cached = _collapsed_taxon_cache.get(node);
        if (cached != null) {
            return cached;
        }
        final org.forester.ws.seqdb.TaxonLineage.Ancestor a = org.forester.analysis.AncestralTaxonomyInference
                .commonTaxonOf(node, _tip_lineages);
        final String label = (a == null) ? "" : a.getName();
        _collapsed_taxon_cache.put(node, label);
        return label;
    }

    /** Test hook: the computed common-taxon label for a collapsed clade rooted at {@code node} (see
     *  {@link #collapsedCommonTaxon}), "" when none is derivable. */
    String collapsedCommonTaxonForTest(final PhylogenyNode node) {
        return collapsedCommonTaxon(node);
    }

    /** Appends a collapsed-clade label built from its COMMON taxon: e.g. "Carnivora (23)" + the found-count suffix. */
    private void addLabelForCollapsedCommonTaxon(final String taxon, final int size, final PhylogenyNode node) {
        _sb.append(taxon);
        _sb.append(" (" + size + ")");
        if ((_found_nodes_0 != null) || (_found_nodes_1 != null)) {
            final int[] res = calcFoundNodesInSubtree(node);
            if (res[0] > 0) {
                _sb.append(" [" + res[0] + "/" + res[1] + "]");
            }
        }
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
        // Where this node's label is anchored: its own point, EXCEPT an external tip in the aligned circular phylogram,
        // whose label is pinned to the common outer ring (circularLabelAnchor). Cull on the ANCHOR (the label position),
        // and reuse it for the label placement + the aligned leader below.
        final Point2D.Double anchor = circularLabelAnchor(node);
        if (isRadialLabelPointInvisible(anchor.x, anchor.y) && !to_graphics_file && !to_pdf) {
            return;
        }
        // Leave a collapsed clade-ROOT label-less radially (as before internal labels were enabled). Its coords ARE
        // valid now (circular positions it in assignCircularDisplayedTipAngles; unrooted via the recursion), but a
        // radial collapse indicator (marker + count) is a deliberate follow-up, so we don't draw a bare label here.
        if (node.isCollapse()) {
            return;
        }
        // show-internal / show-external gates (mirror the rectangular paintNodeData) -- previously moot because this
        // was called for external nodes only; now internal-node labels ride the branch radially too
        if (!isShowInternalDataForThisTab() && !node.isExternal()) {
            return; // (collapsed nodes already returned above)
        }
        if (!isShowExternalDataForThisTab() && node.isExternal()) {
            return;
        }
        // A tip image (if any) is drawn UPRIGHT just outside the tip along its spoke -- in DEVICE space, before the
        // label's rotate/flip below -- so a frog looks like a frog. Offset by HALF the image's spoke footprint so its
        // near edge just clears the node box (not the branch), and the label is pushed out past its far edge via
        // `img_footprint` in `gap` below. Drawn here (before the text early-return) so an image-only tip still shows.
        final double img_theta = hasTipImage(node)
                ? ((_graphics_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) ? _urt_nodeid_angle_map.get(node.getId())
                        : ur_angle)
                : 0;
        final double img_footprint = hasTipImage(node) ? tipImageRadialFootprint(node, img_theta) : 0;
        if (hasTipImage(node)) {
            final double r_off = effectiveNodeHalfBoxSize(node) + 3.0 + (img_footprint / 2.0);
            drawTipImageUpright(g, node, (float) (anchor.x + (r_off * Math.cos(img_theta))),
                    (float) (anchor.y + (r_off * Math.sin(img_theta))), getOptions().getTipImageSize());
        }
        // The TAXONOMY label is drawn via the shared taxonomyLabel part-walker (so the scientific-name part is italic,
        // the rank is included, and the "Abbreviate Scientific Names" option applies) -- matching the rectangular
        // layout. The trailing node name / sequence text is built into _sb and drawn after it.
        final boolean show_tax = node.getNodeData().isHasTaxonomy()
                && (shows(DisplayOption.SHOW_TAX_CODE) || shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES)
                        || shows(DisplayOption.SHOW_TAXONOMY_COMMON_NAMES) || shows(DisplayOption.SHOW_TAX_RANK));
        _sb.setLength(0);
        if (shows(DisplayOption.SHOW_NODE_NAMES) && (node.getName().length() > 0)) {
            _sb.append(" ");
            _sb.append(node.getName());
        }
        if (node.getNodeData().isHasSequence()) {
            final Sequence seq = node.getNodeData().getSequence();
            if (shows(DisplayOption.SHOW_SEQUENCE_ACC) && (seq.getAccession() != null)) {
                _sb.append(" ");
                if (!ForesterUtil.isEmpty(seq.getAccession().getSource())) {
                    _sb.append(seq.getAccession().getSource());
                    _sb.append(":");
                }
                _sb.append(seq.getAccession().getValue());
            }
            // gene symbol + gene name (were missing radially -- rectangular parity, see paintNodeData)
            if (shows(DisplayOption.SHOW_SEQ_SYMBOLS) && (seq.getSymbol().length() > 0)) {
                _sb.append(" ");
                _sb.append(seq.getSymbol());
            }
            if (shows(DisplayOption.SHOW_GENE_NAMES) && (seq.getGeneName().length() > 0)) {
                _sb.append(" ");
                _sb.append(seq.getGeneName());
            }
            if (shows(DisplayOption.SHOW_SEQ_NAMES) && (seq.getName().length() > 0)) {
                _sb.append(" ");
                _sb.append(seq.getName());
            }
        }
        // node PROPERTIES, in the same position the rectangular label puts them (nodeDataAsSB) -- the radial label
        // used to stop at the sequence data, so the "Properties" display option and the Annotation Fields label
        // selection silently did nothing in the circular and unrooted layouts
        if (shows(DisplayOption.SHOW_PROPERTIES) && node.getNodeData().isHasProperties()) {
            final String props = propertiesToString(node);
            if (props.length() > 0) {
                _sb.append(" ");
                _sb.append(props);
            }
        }
        String rest = _sb.toString();
        if (!show_tax && (rest.length() < 1)) {
            return; // nothing to draw
        }
        setColor(g, node, to_graphics_file, to_pdf, is_in_found_nodes, getTreeColorSet().getSequenceColor());
        setFont(g, node);
        final Font base_font = g.getFont();
        // start the label just off the node box (per-node size aware), pushed out past a tip image's spoke footprint
        final float gap = effectiveNodeHalfBoxSize(node) + 3f
                + (hasTipImage(node) ? (float) (img_footprint + TIP_IMAGE_GAP) : 0);
        final int tax_w = show_tax ? taxonomyLabelWidth(node.getNodeData().getTaxonomy(), base_font) : 0;
        // Radial layouts: a full label in a circle runs ~twice the radius (off the canvas), so truncate the node-name/
        // sequence part with an ellipsis to the label budget left after the taxonomy, keeping tree + labels + domains
        // on-canvas. (No cap in a rectangular layout -- radialMaxLabelWidth returns MAX_VALUE there.)
        if ((rest.length() > 0) && isRadialLayout()) {
            rest = TreePanelUtil.truncateToPixelWidth(getFontMetrics(base_font), rest,
                    Math.max(0, radialMaxLabelWidth() - tax_w));
        }
        final int rest_w = getFontMetrics(base_font).stringWidth(rest);
        final double total_w = gap + tax_w + rest_w; // full extent from the node, for the left-half flip
        double m = (_graphics_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                ? (_urt_nodeid_angle_map.get(node.getId()) % TWO_PI) : (ur_angle % TWO_PI);
        _at = g.getTransform();
        boolean need_to_reset = false;
        // The label is anchored at `anchor` (computed at the top): the node, or the outer ring in the aligned circular
        // phylogram. All the rotation/flip logic below is reused verbatim -- only the pivot/start point moves to the ring.
        final float x_coord = (float) anchor.x;
        final float y_coord = (float) anchor.y + (getFontMetrics(base_font).getAscent() / 3.0f);
        // On the left half of the fan the label is flipped so it still reads left-to-right; translating by
        // total_w + gap (rather than just total_w) keeps the SAME small clearance from the node on both halves
        // (otherwise the flipped label butts right up against the node box).
        if (radial_labels) {
            need_to_reset = true;
            boolean left = false;
            if ((m > HALF_PI) && (m < ONEHALF_PI)) {
                m -= PI;
                left = true;
            }
            g.rotate(m, x_coord, (float) anchor.y);
            if (left) {
                g.translate(-(total_w + gap), 0);
            }
        } else if ((m > HALF_PI) && (m < ONEHALF_PI)) {
            need_to_reset = true;
            g.translate(-(total_w + gap), 0);
        }
        float x = x_coord + gap;
        if (show_tax) {
            x += taxonomyLabel(g, node.getNodeData().getTaxonomy(), x, y_coord, to_pdf, true);
        }
        if (rest.length() > 0) {
            g.setFont(base_font);
            TreePanel.drawString(rest, x, y_coord, g);
        }
        if (need_to_reset) {
            g.setTransform(_at);
        }
        // Aligned circular phylogram: a faint dotted radial leader from the tip (branch-length radius) OUT to its label
        // on the ring -- drawn HERE (transform restored to device space), so it appears EXACTLY when a label was drawn
        // (past the collapse/content/dynamic-hiding/off-screen returns above) and its endpoint IS the label anchor, so it
        // can never dangle to an empty ring point. Reuses LEADER_STROKE + connectorColor like the rectangular aligned
        // tip->label connector. Skipped for a deepest tip already at the ring (leader length < 1px).
        if (node.isExternal() && isAlignedCircularPhylogram()) {
            final double dx = anchor.x - node.getXcoord();
            final double dy = anchor.y - node.getYcoord();
            if (((dx * dx) + (dy * dy)) >= 1.0) {
                final Stroke leader_saved_stroke = g.getStroke();
                final Color leader_saved_color = g.getColor();
                g.setStroke(LEADER_STROKE);
                g.setColor(connectorColor());
                drawLine(node.getXcoord(), node.getYcoord(), anchor.x, anchor.y, g);
                g.setStroke(leader_saved_stroke);
                g.setColor(leader_saved_color);
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
                if (effectiveTipLabelDirection() == Options.TIP_LABEL_DIRECTION.HORIZONTAL) {
                    // UPRIGHT tip label, centred under (root-top) / over (root-bottom) the tip -- the cleanest look for
                    // short names / sparse trees, and the one a rotated-bitmap export can't produce (its text tilts).
                    paintTipLabelHorizontal(g, node, to_graphics_file, to_pdf, is_in_found_nodes);
                } else {
                    withNodeTextFrame(g, labelTextStartX(node), node.getYcoord(), tipLabelAngle(),
                            () -> paintNodeData(g, node, to_graphics_file, to_pdf, is_in_found_nodes, 0));
                }
                // renderable domain architecture: a vertical bar in the COMMON aligned column past the labels (boxes
                // ride R, no labels) -- lined up like the horizontal layout, matching what depthLabelReserve() reserves.
                paintDomainsVertical(g, node, to_pdf, to_graphics_file);
            } else {
                // internal-node label: horizontal, right-aligned, LEFT of the branch midpoint. This path deliberately
                // does NOT route through paintNodeData, so the renderable domain overlay it draws is DEFERRED for
                // internal nodes in a vertical orientation (increment 1): it needs the rotated-frame re-anchoring the
                // tip labels get; niche for root-top/bottom (focused on internal LABELS, small trees). External nodes
                // still draw it via the tilted paintNodeData path above.
                paintInternalLabelLeftVertical(g, node, to_pdf, to_graphics_file);
            }
        } else {
            paintNodeData(g, node, to_graphics_file, to_pdf, is_in_found_nodes, 0);
            paintNodeWithRenderableData(g, node, to_graphics_file, to_pdf);
        }
    }

    /** Draws an external tip's renderable domain architecture in a VERTICAL orientation: the boxes ride the rotation R
     *  into a thin vertical track just past the tilted tip label (domain-name labels suppressed -- they would collide
     *  with neighbouring tips' tracks). Per tip at its own depth, so the track hangs off the tip. No-op unless domains
     *  are shown, external data is shown, and the tip carries a renderable architecture. */
    private void paintDomainsVertical(final Graphics2D g, final PhylogenyNode node,
                                      final boolean to_pdf, final boolean to_graphics_file) {
        if (!shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES) || !isShowExternalDataForThisTab()) {
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
        // Height tracks the ACTUAL tip-row spacing (getYdistance()) -- NOT the font height even under dynamic hiding
        // (which pinned it to a fixed size, so it neither responded to vertical zoom NOR shrank when the rows pack
        // tight, crowding the boxes). Clamped into [MIN, MAX] so it stays readable zoomed out and bars-not-blocks
        // zoomed in (shared with the other domain-height site via TreePanelUtil.domainBoxHeight).
        final float yd = getYdistance();
        final int hgt = TreePanelUtil.domainBoxHeight(yd, DOMAIN_STRUCTURE_HEIGHT_MIN, DOMAIN_STRUCTURE_HEIGHT_MAX);
        rds.setRenderingHeight(hgt);
        // draw in the COMMON aligned column (past the deepest tip + the longest label's depth footprint) so every
        // tip's track lines up -- the alignment the horizontal layout gives via alignedPhylogramDomainColumnX(), and
        // exactly the depth depthLabelReserve() reserves. Rides R into a vertical bar.
        rds.render(verticalDomainColumnStart(), node.getYcoord() - (hgt / 2.0f), g, this, to_pdf, false);
    }

    /** The COMMON depth (logical x) where the domain track starts in a VERTICAL orientation, so all tips' tracks line
     *  up in one band: past the deepest tip + the longest tip label's tilt-projected DEPTH footprint. Mirrors the
     *  {@link #depthLabelReserve()} label-end (which reserves exactly this), so the aligned tracks never clip. */
    private float verticalDomainColumnStart() {
        final double a = tipLabelAngle();
        final int line_h = getFontMetricsForLargeDefaultFont().getHeight();
        final int anchor_offset = maxTipEffectiveHalfBoxSize() + LABEL_GAP_AFTER_NODE_SHAPE;
        final int upright_gap = (effectiveTipLabelDirection() == Options.TIP_LABEL_DIRECTION.HORIZONTAL)
                ? TIP_LABEL_DEPTH_GAP : 0;
        final double text_depth = (_length_of_longest_text_only * Math.abs(Math.sin(a)))
                + (line_h * Math.abs(Math.cos(a)));
        final float deepest_tip_x = getControlPanel().isDrawPhylogram()
                ? (float) (((displayedMaxDistanceToRoot() - rootBranchInMaxDistance()) * getXcorrectionFactor())
                        + _phylogeny.getRoot().getXcoord())
                : getPhylogeny().getFirstExternalNode().getXcoord();
        return (float) (deepest_tip_x + anchor_offset + text_depth + upright_gap + VERTICAL_DOMAIN_GAP);
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
        final boolean root_bottom = getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_BOTTOM;
        final double gap = TIP_LABEL_DEPTH_GAP;
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

    final private void paintNodeWithRenderableData(final Graphics2D g,
                                                   final PhylogenyNode node,
                                                   final boolean to_graphics_file,
                                                   final boolean to_pdf) {
        if (isNodeDataInvisible(node) && !(to_graphics_file || to_pdf)) {
            return;
        }
        if ((!isShowInternalDataForThisTab() && !node.isExternal())) {
            return;
        }
        if ((!isShowExternalDataForThisTab() && node.isExternal())) {
            return;
        }
        if (shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES) && node.getNodeData().isHasSequence()
                && (node.getNodeData().getSequence().getDomainArchitecture() != null) && (node.getNodeData()
                .getSequence().getDomainArchitecture() instanceof RenderableDomainArchitecture)) {
            RenderableDomainArchitecture rds = null;
            try {
                rds = (RenderableDomainArchitecture) node.getNodeData().getSequence().getDomainArchitecture();
            } catch (final ClassCastException cce) {
                cce.printStackTrace();
            }
            if (rds != null) {
                final float y = getYdistance(); // track the actual row spacing (see the horizontal path above)
                final int h = TreePanelUtil.domainBoxHeight(y, DOMAIN_STRUCTURE_HEIGHT_MIN, DOMAIN_STRUCTURE_HEIGHT_MAX);
                rds.setRenderingHeight(h);
                // Domain architectures always line up in a common right-edge column (the standard,
                // iTOL-style comparable layout); the phylogram column is past the deepest tip + longest
                // label, the cladogram column just past the aligned tips.
                if (getControlPanel().isDrawPhylogram()) {
                    rds.render(alignedPhylogramDomainColumnX(),
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
        // a subtree draws a fixed short stub (d stays getXdistance); a full-tree phylogram root edge is drawn to scale
        if (displayedRootBranchLength() > 0.0) {
            double root_dtp = displayedRootBranchLength();
            if (breakLongBranchesActive() && (root_dtp > breakLongBranchCap())) {
                root_dtp = breakLongBranchCap();
            }
            d = (float) (getXcorrectionFactor() * root_dtp);
        }
        if (d < MIN_ROOT_LENGTH) {
            d = MIN_ROOT_LENGTH;
        }
        if (!shows(DisplayOption.WIDTH_BRANCHES) || (PhylogenyMethods.getBranchWidthValue(root) == 1)) {
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
        final double x2 = x1 + (displayScaleDistance() * getXcorrectionFactor());
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
        if (displayScaleLabel() != null) {
            g.drawString(displayScaleLabel(), (x1 + 2), y3 - 2);
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
        final int y_top = (int) Math.round(y_bottom - (displayScaleDistance() * getXcorrectionFactor()));
        final int x_tick = x1 + 8;
        final int x_bar = x1 + 4;
        g.setFont(getTreeFontSet().getSmallFont());
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        final Stroke s = g.getStroke();
        g.setStroke(STROKE_1);
        drawLine(x1, y_bottom, x_tick, y_bottom, g);
        drawLine(x1, y_top, x_tick, y_top, g);
        drawLine(x_bar, y_bottom, x_bar, y_top, g);
        if (displayScaleLabel() != null) {
            g.drawString(displayScaleLabel(), x_tick + 3, (y_bottom + y_top) / 2);
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
        // grid spans root (origin_x) to the deepest tip; origin_x already carries the root-edge offset that
        // displayedMaxDistanceToRoot() also includes, so subtract rootBranchInMaxDistance() to avoid extending
        // the grid a root-branch past the tips (same double-count the domain column above avoids)
        final float max_x = (float) (origin_x
                + ((displayedMaxDistanceToRoot() - rootBranchInMaxDistance()) * getXcorrectionFactor()));
        final float[] xs = TreePanelUtil.scaleGridLineXs(origin_x, spacing, max_x);
        if (xs.length == 0) {
            return;
        }
        // Use the export canvas extent only when it is real; the on-screen paint passes to_pdf=false with
        // height 0, so fall back to the panel height there.
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

    /** Optional (off by default) faint reference lines across the tree at the active time axis's tick positions -- the
     *  FINE geologic band's interval boundaries (e.g. the Early/Middle/Late Triassic and the Triassic/Permian
     *  boundaries) OR the CALENDAR year ticks. Drawn BEHIND the tree, mirroring {@link #paintScaleGrid} (vertical lines
     *  in a horizontal layout; logical lines that ride R into horizontal grid lines in a vertical orientation; the
     *  CIRCULAR analogue is the geologic boundary rings / the calendar year rings). */
    private void paintTimeAxisGridLines(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file,
                                        final int graphics_file_y, final int graphics_file_height) {
        if (!isTimeAxisGrid()) {
            return;
        }
        final double corr = getXcorrectionFactor();
        if (corr <= 0) {
            return;
        }
        final float origin_x = _phylogeny.getRoot().getXcoord();
        // the grid-line x positions come from whichever time axis is active: the fine geologic band's boundaries
        // (branch-aligned via ageToX, bounded by the youngest tip), or the calendar year ticks. Skip the tips/root edge.
        final java.util.List<Integer> xs = new java.util.ArrayList<>();
        if (geologicAxisApplies()) {
            final double root_age = timeAxisRootAgeMa();
            if (root_age <= 0) {
                return;
            }
            final GeologicTimeScale.Rank fine = GeologicTimeScale.bandRanks(root_age)[1];
            final double young_bound = timeAxisYoungestAgeMa(); // >0 for a fossil-only tree: no grid lines past the tips
            for (final GeologicTimeScale.Interval iv : GeologicTimeScale.overlapping(fine, young_bound, root_age)) {
                final double b = Math.min(root_age, iv.oldMa());
                if ((b > young_bound) && (b < root_age)) {
                    xs.add((int) Math.round(ageToX(b, root_age, origin_x, corr)));
                }
            }
        }
        else if (calendarAxisApplies()) {
            final double max_dist = getMaxDistanceToRoot();
            final float tip_x = (float) (origin_x + (max_dist * corr));
            final double present = timeAxisPresentDate();
            if ((present <= 0) || (max_dist <= 0)) {
                return;
            }
            final double root_year = present - max_dist;
            for (final double year : TreePanelUtil.calendarTickYears(root_year, present)) {
                final int x = (int) Math.round(origin_x + ((year - root_year) * corr));
                if ((x > Math.round(origin_x)) && (x < tip_x)) {
                    xs.add(x);
                }
            }
        }
        else {
            return;
        }
        final boolean use_export_extent = (to_pdf || to_graphics_file) && (graphics_file_height > 0);
        final int top = isVerticalOrientation() ? 0 : (use_export_extent ? graphics_file_y : 0);
        final int bottom = isVerticalOrientation() ? treeBreadthExtent()
                : (use_export_extent ? (graphics_file_y + graphics_file_height) : getHeight());
        final Color saved_color = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setColor(scaleGridColor(to_pdf, to_graphics_file));
        g.setStroke(STROKE_05);
        for (final int x : xs) {
            drawLine(x, top, x, bottom, g);
        }
        g.setColor(saved_color);
        g.setStroke(saved_stroke);
    }

    /** Screen-only guidance (never drawn in exports -- it is UI help, not figure content): a dated tree carries a time
     *  axis (Geologic / Calendar), but a time axis needs a PHYLOGRAM to draw. When the tree is shown as a CLADOGRAM the
     *  axis silently steps aside; a faint one-line hint at the bottom of the view says why, and how to get it back. */
    private void paintTimeAxisHint(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        if (to_pdf || to_graphics_file || (_phylogeny == null) || (getControlPanel() == null)
                || getControlPanel().isDrawPhylogram()) {
            return;
        }
        final Options.TIME_AXIS_TYPE t = effectiveTimeAxisType();
        if ((t != Options.TIME_AXIS_TYPE.GEOLOGIC) && (t != Options.TIME_AXIS_TYPE.CALENDAR)) {
            return;
        }
        final java.awt.Rectangle vr = getVisibleRect();
        final Font saved_font = g.getFont();
        final Color saved_color = g.getColor();
        g.setFont(saved_font.deriveFont(Font.ITALIC));
        g.setColor(scaleInkColor(false, false));
        g.drawString("Time axis hidden: display as a phylogram (P) to show it", vr.x + 8, (vr.y + vr.height) - 8);
        g.setFont(saved_font);
        g.setColor(saved_color);
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
        final double max_dist = numericScaleAxisMaxDist();
        final double[] ticks = TreePanelUtil.scaleAxisTickValues(max_dist, getScaleDistance());
        if (ticks.length == 0) {
            return;
        }
        final Font saved_font = g.getFont();
        final Color saved_color = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = g.getFontMetrics();
        // Screen: FLOAT the axis at the viewport bottom so it never scrolls out of view when zoomed in (PearTree-
        // style), exactly like the viewport-fixed scale bar. A file export stays anchored to the tree/export bottom
        // so figures remain WYSIWYG.
        final Rectangle vr = getVisibleRect();
        final int bottom = TreePanelUtil.scaleAxisFloatingBottom(to_pdf, to_graphics_file, graphics_file_y,
                graphics_file_height, getHeight(), vr.y + vr.height);
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

    /** The CALENDAR (absolute-date) time axis along the bottom: a horizontal ruler with a tick + calendar-year label at
     *  each "nice" year (like {@link #paintScaleAxis}, but the tick values are calendar years and the position maps the
     *  branch-length depth to time). Calibrated by {@link #timeAxisPresentDate()} (the most-recent tip = the present) so
     *  {@code root_year = present - maxDist}; a node at distance-from-root {@code d} sits at year {@code root_year + d}. */
    private void paintCalendarAxis(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file,
                                   final int graphics_file_y, final int graphics_file_height) {
        final double present = timeAxisPresentDate();
        final double corr = getXcorrectionFactor();
        final double max_dist = getMaxDistanceToRoot();
        if ((present <= 0) || (corr <= 0) || (max_dist <= 0)) {
            return;
        }
        final double root_year = present - max_dist;
        final double[] years = TreePanelUtil.calendarTickYears(root_year, present);
        if (years.length == 0) {
            return;
        }
        final float origin_x = _phylogeny.getRoot().getXcoord();
        final Font saved_font = g.getFont();
        final Color saved_color = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = g.getFontMetrics();
        final Rectangle vr = getVisibleRect();
        final int bottom = TreePanelUtil.scaleAxisFloatingBottom(to_pdf, to_graphics_file, graphics_file_y,
                graphics_file_height, getHeight(), vr.y + vr.height);
        // draw one CALENDAR_AXIS_EDGE_GAP up from the edge so the ruler + labels don't sit flush against the border
        final int axis_y = bottom - scaleAxisBandHeight() - CALENDAR_AXIS_EDGE_GAP;
        final int label_baseline = axis_y + SCALE_AXIS_TICK_LEN + fm.getAscent() + 1;
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        g.setStroke(STROKE_1);
        final int x_end = (int) Math.round(origin_x + (max_dist * corr));
        drawLine(Math.round(origin_x), axis_y, x_end, axis_y, g); // root (root_year) -> deepest tip (present)
        int last_label_right = Integer.MIN_VALUE;
        for (final double year : years) {
            final int x = (int) Math.round(origin_x + ((year - root_year) * corr));
            drawLine(x, axis_y, x, axis_y + SCALE_AXIS_TICK_LEN, g);
            final String label = calendarYearLabel(year);
            final int half = fm.stringWidth(label) / 2;
            if ((x - half) >= (last_label_right + SCALE_AXIS_LABEL_GAP)) {
                g.drawString(label, x - half, label_baseline);
                last_label_right = x + half;
            }
        }
        g.setFont(saved_font);
        g.setColor(saved_color);
        g.setStroke(saved_stroke);
    }

    /** A calendar-year tick label: the whole year (e.g. "2021") -- ticks are whole-year multiples. */
    private static String calendarYearLabel(final double year) {
        return String.valueOf((int) Math.round(year));
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
        final int anchored_x = (int) Math.round(r_root.x); // the tree-anchored breadth position (used for exports)
        final int tip_side_x = (int) Math.round(screenPoint(origin_x, ruler_ly - 16.0).x);
        final int in = (tip_side_x >= anchored_x) ? 1 : -1; // +1 = tree to the right of the ruler (ticks point right)
        // Screen: FLOAT the ruler to the viewport breadth EDGE on its own side (away from the tree) so it stays
        // visible when zoomed/scrolled along the breadth (PearTree-style); every export keeps the tree-anchored
        // breadth position (treeBreadthExtent) so figures remain WYSIWYG. The tick DEPTH positions (device-y from
        // screenPoint) are unchanged, so they stay aligned with the branches.
        final Rectangle vr = getVisibleRect();
        final int ruler_x = TreePanelUtil.scaleAxisRulerX(to_pdf, to_graphics_file, anchored_x, in, vr.x, vr.width);
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

    /** The CALENDAR time axis in a VERTICAL orientation: a ruler down one BREADTH side with a tick + upright
     *  calendar-year label at each nice year, exactly like {@link #paintScaleAxisVertical} but with calendar-year
     *  ticks. Floats to the viewport breadth edge on screen; tree-anchored in exports. */
    private void paintCalendarAxisVertical(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        if (_orientation_R == null) {
            return;
        }
        final double corr = getXcorrectionFactor();
        final double max_dist = getMaxDistanceToRoot();
        final double present = timeAxisPresentDate();
        if ((corr <= 0.0) || (max_dist <= 0.0) || (present <= 0.0)) {
            return;
        }
        final double root_year = present - max_dist;
        final double[] years = TreePanelUtil.calendarTickYears(root_year, present);
        if (years.length == 0) {
            return;
        }
        final float origin_x = _phylogeny.getRoot().getXcoord();
        final double ruler_ly = treeBreadthExtent() - 2.0;
        final Font saved_font = g.getFont();
        final Color saved_color = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = g.getFontMetrics();
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        g.setStroke(STROKE_1);
        final Point2D.Double r_root = screenPoint(origin_x, ruler_ly);
        final Point2D.Double r_tip = screenPoint(origin_x + (max_dist * corr), ruler_ly);
        final int anchored_x = (int) Math.round(r_root.x);
        final int tip_side_x = (int) Math.round(screenPoint(origin_x, ruler_ly - 16.0).x);
        final int in = (tip_side_x >= anchored_x) ? 1 : -1;
        final Rectangle vr = getVisibleRect();
        // float to the viewport breadth edge, then step one CALENDAR_AXIS_EDGE_GAP INWARD (toward the tree, i.e. by
        // +in) so the ruler + its year labels don't sit flush against the window border
        final int ruler_x = TreePanelUtil.scaleAxisRulerX(to_pdf, to_graphics_file, anchored_x, in, vr.x, vr.width)
                + (in * CALENDAR_AXIS_EDGE_GAP);
        drawLine(ruler_x, (int) Math.round(r_root.y), ruler_x, (int) Math.round(r_tip.y), g);
        final int center = (fm.getAscent() - fm.getDescent()) / 2;
        final int min_label_gap = fm.getHeight() + SCALE_AXIS_LABEL_GAP;
        int last_label_y = Integer.MIN_VALUE;
        for (final double year : years) {
            final int ty = (int) Math.round(screenPoint(origin_x + ((year - root_year) * corr), ruler_ly).y);
            drawLine(ruler_x, ty, ruler_x + (in * SCALE_AXIS_TICK_LEN), ty, g);
            final String label = calendarYearLabel(year);
            if ((last_label_y == Integer.MIN_VALUE) || (Math.abs(ty - last_label_y) >= min_label_gap)) {
                final int lx = (in > 0) ? (ruler_x + SCALE_AXIS_TICK_LEN + 2)
                        : (ruler_x - SCALE_AXIS_TICK_LEN - 2 - fm.stringWidth(label));
                g.drawString(label, lx, ty + center);
                last_label_y = ty;
            }
        }
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
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                || breakLongBranchesActive()) { // a capped tree has no single linear distance scale (see below)
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
        if (isVerticalOrientation()) {
            return 0;
        }
        if (geologicAxisApplies()) {
            return geologicAxisReserve(); // the geologic axis takes the bottom strip (+ an edge gap) vs the numeric axis
        }
        if (calendarAxisApplies()) {
            // the calendar tick axis reserves one label row (+ an edge gap from the window border), like the numeric
            // scale axis -- but only when it actually yields year ticks (a sub-year span within one calendar year has
            // none, so it draws nothing -> reserve nothing)
            final double present = timeAxisPresentDate();
            if (TreePanelUtil.calendarTickYears(present - getMaxDistanceToRoot(), present).length > 0) {
                return scaleAxisBandHeight() + CALENDAR_AXIS_EDGE_GAP;
            }
            return 0;
        }
        if (!getOptions().isShowScaleAxis() || !scaleAxisAppliesToLayout()) {
            return 0;
        }
        return scaleAxisBandHeight();
    }

    /** Breadth (px) reserved for the VERTICAL scale-axis ruler (the tick + the widest tick number + padding), so the
     *  vertical ruler + its upright labels sit clear of the tips just past the last tip. 0 unless the axis is shown in
     *  a vertical orientation on a phylogram with a scale. Reserved in the breadth budget (calcParametersForPainting)
     *  and the logical breadth extent (logicalTreeExtent), so paint, fit, and scroll all agree. */
    int verticalScaleAxisReserve() {
        if (!isVerticalOrientation()) {
            return 0;
        }
        if (geologicAxisApplies()) {
            return geologicAxisReserve(); // the geologic two-band axis takes the breadth side band (+ edge gap) in vertical
        }
        if (calendarAxisApplies()) {
            final double present = timeAxisPresentDate();
            final double[] years = TreePanelUtil.calendarTickYears(present - getMaxDistanceToRoot(), present);
            if (years.length == 0) {
                return 0;
            }
            final FontMetrics fmc = getFontMetrics(getTreeFontSet().getSmallFont());
            int max_year_label = 0;
            for (final double y : years) {
                max_year_label = Math.max(max_year_label, fmc.stringWidth(calendarYearLabel(y)));
            }
            return SCALE_AXIS_TICK_LEN + max_year_label + 8 + CALENDAR_AXIS_EDGE_GAP;
        }
        if (!getOptions().isShowScaleAxis() || !getControlPanel().isDrawPhylogram() || (getScaleDistance() <= 0.0)
                || breakLongBranchesActive()) { // capping suppresses the numeric axis (see scaleAxisAppliesToLayout)
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
        return ((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite())
                ? Color.BLACK
                : getTreeColorSet().getBranchLengthColor();
    }

    /** Faint, theme-aware color for the scale grid: the background nudged slightly toward the branch-length color. */
    private Color scaleGridColor(final boolean to_pdf, final boolean to_graphics_file) {
        if ((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) {
            return new Color(220, 220, 220);
        }
        return TreePanelUtil.blend(getTreeColorSet().getBackgroundColor(),
                getTreeColorSet().getBranchLengthColor(), SCALE_GRID_BLEND);
    }

    /** Concentric distance RINGS for a circular PHYLOGRAM: faint circles at scale-distance intervals from the ring
     *  centre, each labelled with its distance value at 12 o'clock, so distance-from-root reads off at any angle (the
     *  radial analogue of the rectangular scale grid). Drawn BEHIND the tree; a no-op unless the circular layout is a
     *  phylogram. Labels are decimated so they never stack when the rings are dense. */
    private void paintCircularScaleRings(final Graphics2D g, final int cx, final int cy, final int radius,
                                         final boolean to_pdf, final boolean to_graphics_file) {
        // the rings ARE the circular scale, so they follow the same "Scale" option that draws the bottom bar in the
        // rectangular/unrooted layouts (the circular phylogram POSITIONS by branch length regardless; only the scale
        // overlay is gated). Off unless it is a circular phylogram with the scale shown.
        if (!isCircularPhylogram() || !getOptions().isShowScale() || (radius <= 0)) {
            return;
        }
        final double spacing = getScaleDistance();
        final double max = getMaxDistanceToRoot();
        if ((spacing <= 0) || (max <= 0)) {
            return;
        }
        final Color saved = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        final Font saved_font = g.getFont();
        g.setStroke(STROKE_05);
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = getTreeFontSet().getFontMetricsSmall();
        final int line_h = fm.getHeight();
        final Color ring_c = scaleGridColor(to_pdf, to_graphics_file);
        final Color label_c = scaleInkColor(to_pdf, to_graphics_file);
        int last_label_y = Integer.MAX_VALUE; // rings go inner->outer, so labels move UP the top spoke
        int k = 1;
        for (double d = spacing; (d <= (max + 1e-9)) && (k <= 1000); d += spacing, ++k) {
            final int rr = (int) Math.round((d / max) * radius);
            if (rr <= 0) {
                continue;
            }
            g.setColor(ring_c);
            g.drawOval(cx - rr, cy - rr, 2 * rr, 2 * rr);
            final int ly = (cy - rr) + fm.getAscent() + 1; // just inside the ring at the top
            if ((last_label_y - ly) >= line_h) { // enough vertical gap since the last drawn label
                final String label = TreePanelUtil.formatCompactNumber(d);
                g.setColor(label_c);
                g.drawString(label, cx - (fm.stringWidth(label) / 2f), ly);
                last_label_y = ly;
            }
        }
        g.setColor(saved);
        g.setStroke(saved_stroke);
        g.setFont(saved_font);
    }

    /** Concentric CALENDAR-year rings for a circular PHYLOGRAM: the radial axis is time, so a faint ring at each nice
     *  year (radius = distance-from-root scaled) labelled with the year up the top spoke -- the polar analogue of the
     *  rectangular calendar ruler. A no-op unless {@link #calendarRingsApplyCircular()}. */
    private void paintCalendarRingsCircular(final Graphics2D g, final int cx, final int cy, final int radius,
                                            final boolean to_pdf, final boolean to_graphics_file) {
        if (!calendarRingsApplyCircular() || (radius <= 0)) {
            return;
        }
        final double present = timeAxisPresentDate();
        final double max = getMaxDistanceToRoot();
        if ((present <= 0) || (max <= 0)) {
            return;
        }
        final double root_year = present - max;
        final double[] years = TreePanelUtil.calendarTickYears(root_year, present);
        if (years.length == 0) {
            return;
        }
        final Color saved = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        final Font saved_font = g.getFont();
        g.setStroke(STROKE_05);
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = getTreeFontSet().getFontMetricsSmall();
        final int line_h = fm.getHeight();
        final Color ring_c = scaleGridColor(to_pdf, to_graphics_file);
        final Color label_c = scaleInkColor(to_pdf, to_graphics_file);
        int last_label_y = Integer.MAX_VALUE; // rings inner->outer (old->recent), labels move UP the top spoke
        for (final double year : years) {
            final int rr = (int) Math.round(((year - root_year) / max) * radius);
            if ((rr <= 0) || (rr > radius)) {
                continue;
            }
            g.setColor(ring_c);
            g.drawOval(cx - rr, cy - rr, 2 * rr, 2 * rr);
            final int ly = (cy - rr) + fm.getAscent() + 1;
            if ((last_label_y - ly) >= line_h) {
                final String label = calendarYearLabel(year);
                g.setColor(label_c);
                g.drawString(label, cx - (fm.stringWidth(label) / 2f), ly);
                last_label_y = ly;
            }
        }
        g.setColor(saved);
        g.setStroke(saved_stroke);
        g.setFont(saved_font);
    }

    /** age (Ma) -> radius (px) for the circular geologic axis, aligned to the TREE's OWN radial scale (radius per
     *  distance/time unit == {@code radius / maxDistanceToRoot}): the root (oldest age = {@code root_age}) at the centre
     *  (r=0), each Ma one distance unit outward, so a fossil-only tree's rings line up with the branch radii instead of
     *  pinning age 0 to the outer tip ring. For a clock tree (maxDist == root_age) this equals the older
     *  {@code (root_age-age)/root_age * radius} mapping. Callers are gated on {@link #geologicRingsApplyCircular()},
     *  which guarantees maxDistanceToRoot &gt; 0. */
    private int geologicRadiusPx(final double age, final double root_age, final int radius) {
        return (int) Math.round(((root_age - age) / getMaxDistanceToRoot()) * radius);
    }

    /** Concentric ICS geologic bands for a circular PHYLOGRAM: the radial axis IS time (radius = distance-from-root), so
     *  each geologic PERIOD fills a translucent coloured annulus from its old-boundary radius (inner) to its
     *  young-boundary radius (outer), drawn BEHIND the tree, with faint EPOCH boundary rings for the finer subdivision
     *  and the period name labelled up the top (12 o'clock) spoke. The polar analogue of the rectangular two-band axis;
     *  a no-op unless {@link #geologicRingsApplyCircular()}. */
    private void paintGeologicRingsCircular(final Graphics2D g, final int cx, final int cy, final int radius,
                                            final boolean to_pdf, final boolean to_graphics_file) {
        if (!geologicRingsApplyCircular() || (radius <= 0)) {
            return;
        }
        final double root_age = timeAxisRootAgeMa();
        if (root_age <= 0) {
            return;
        }
        final double young_bound = timeAxisYoungestAgeMa(); // >0 for a fossil-only tree: the outer ring is youngest, not 0
        final Color saved = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        // the coarse+fine rank pair adapts to the tree's depth (Period/Epoch -> Era/Period -> Eon/Era); the coarse
        // rank fills the coloured annuli, the fine rank draws the boundary rings
        final GeologicTimeScale.Rank[] ranks = GeologicTimeScale.bandRanks(root_age);
        // translucent coarse-rank annuli (age -> radius: age root_age at the centre, the youngest tip at the outer ring)
        for (final GeologicTimeScale.Interval iv : GeologicTimeScale.overlapping(ranks[0], young_bound, root_age)) {
            final double young = Math.max(young_bound, iv.youngMa());
            final double old = Math.min(root_age, iv.oldMa());
            final int r_outer = geologicRadiusPx(young, root_age, radius); // younger -> larger radius
            final int r_inner = geologicRadiusPx(old, root_age, radius);
            if ((r_outer - r_inner) <= 0) {
                continue;
            }
            final Area ring = new Area(new Ellipse2D.Double(cx - r_outer, cy - r_outer, 2 * r_outer, 2 * r_outer));
            ring.subtract(new Area(new Ellipse2D.Double(cx - r_inner, cy - r_inner, 2 * r_inner, 2 * r_inner)));
            g.setColor(geologicRingFill(iv.color(), to_pdf, to_graphics_file));
            g.fill(ring);
        }
        // faint fine-rank boundary rings (the finer subdivision, as ring outlines within the coarse-rank colours) --
        // the circular analogue of the geologic grid lines, so they follow the same optional (off by default) toggle
        if (isTimeAxisGrid()) {
            g.setStroke(STROKE_05);
            g.setColor(scaleGridColor(to_pdf, to_graphics_file));
            for (final GeologicTimeScale.Interval iv : GeologicTimeScale.overlapping(ranks[1], young_bound, root_age)) {
                final int rr = geologicRadiusPx(Math.max(young_bound, iv.youngMa()), root_age, radius);
                if ((rr > 0) && (rr < radius)) {
                    g.drawOval(cx - rr, cy - rr, 2 * rr, 2 * rr);
                }
            }
        }
        g.setColor(saved);
        g.setStroke(saved_stroke);
    }

    /** The coarse-rank band names (Period, or Era/Eon for a deep tree) for the circular geologic bands, drawn UP the
     *  top (12 o'clock) spoke ON TOP of the tree
     *  (so a branch at 12 o'clock can't hide them), each at its annulus mid-radius, decimated so they never stack. Split
     *  from {@link #paintGeologicRingsCircular} (which fills the annuli behind the tree) so the labels stay legible. */
    private void paintGeologicRingLabelsCircular(final Graphics2D g, final int cx, final int cy, final int radius,
                                                 final boolean to_pdf, final boolean to_graphics_file) {
        if (!geologicRingsApplyCircular() || (radius <= 0)) {
            return;
        }
        final double root_age = timeAxisRootAgeMa();
        if (root_age <= 0) {
            return;
        }
        final Font saved_font = g.getFont();
        final Color saved = g.getColor();
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = getTreeFontSet().getFontMetricsSmall();
        final int line_h = fm.getHeight();
        final double young_bound = timeAxisYoungestAgeMa(); // >0 for a fossil-only tree: no annuli past the tips
        // collect each coarse-rank band's label baseline y up the top spoke (at its annulus mid-radius), then draw
        // INNER->OUTER (ascending radius) greedily keeping a >=line_h gap -- order-independent of what overlapping() gives
        final java.util.List<GeologicTimeScale.Interval> periods = new java.util.ArrayList<>(
                GeologicTimeScale.overlapping(GeologicTimeScale.bandRanks(root_age)[0], young_bound, root_age));
        periods.sort((x, y) -> Double.compare(x.youngMa(), y.youngMa())); // young first = outer first (larger radius)
        int last_label_y = Integer.MIN_VALUE / 2; // outer->inner, ly increases; keep a >=line_h gap (half-min: no overflow)
        for (final GeologicTimeScale.Interval iv : periods) {
            final int r_outer = geologicRadiusPx(Math.max(young_bound, iv.youngMa()), root_age, radius);
            final int r_inner = geologicRadiusPx(Math.min(root_age, iv.oldMa()), root_age, radius);
            if ((r_outer - r_inner) <= 0) {
                continue;
            }
            final int ly = (cy - ((r_inner + r_outer) / 2)) + (fm.getAscent() / 2);
            if ((ly - last_label_y) < line_h) {
                continue;
            }
            final int lx = cx - (fm.stringWidth(iv.name()) / 2);
            final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
            if (!bw) {
                // an opaque ICS-coloured chip so the name reads over both the pale band and any branch at 12 o'clock;
                // ink contrast-picked against that chip colour
                g.setColor(iv.color());
                g.fillRect(lx - 2, ly - fm.getAscent(), fm.stringWidth(iv.name()) + 4, line_h);
                g.setColor(labelInkOn(iv.color()));
            }
            else {
                // B&W export: no chip; the annulus behind is the light grey geologicRingFill paints, so use black ink
                // (NOT labelInkOn(iv.color()), which would pick white for a dark period and vanish on the grey)
                g.setColor(Color.BLACK);
            }
            g.drawString(iv.name(), lx, ly);
            last_label_y = ly;
        }
        g.setFont(saved_font);
        g.setColor(saved);
    }

    /** Optional Ma age labels at the COARSE band's boundary RADII, up the top (12 o'clock) spoke, for the circular
     *  geologic axis -- the circular analogue of the rectangular boundary-age row. Drawn ON TOP with a faint background
     *  so they read over the band + any branch, decimated so they never stack. Gated on "Geologic Boundary Ages". */
    private void paintGeologicBoundaryAgesCircular(final Graphics2D g, final int cx, final int cy, final int radius,
                                                   final boolean to_pdf, final boolean to_graphics_file) {
        if (!geologicRingsApplyCircular() || (radius <= 0) || !isTimeAxisAges()) {
            return;
        }
        final double root_age = timeAxisRootAgeMa();
        if (root_age <= 0) {
            return;
        }
        final Font saved_font = g.getFont();
        final Color saved = g.getColor();
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = getTreeFontSet().getFontMetricsSmall();
        final int line_h = fm.getHeight();
        final Color bg = getTreeColorSet().getBackgroundColor();
        final Color ink = scaleInkColor(to_pdf, to_graphics_file);
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
        // each coarse interval's OLDER edge is a boundary; youngest first = outer first (smallest ly), decimate down.
        // Skip boundaries in the crowded OUTER region (near the tip ring, where the thin outer annuli + their name
        // chips collide) -- scale-adaptive, so a wider zoom shows more ages; the deep, well-spaced boundaries stay.
        // skip the crowded outer region, but never go negative (on a tiny disc that would drop ALL ages, not just the
        // outer ones) -- then only the age-0 edge is excluded and the deep, well-spaced boundaries still show
        final int outer_cutoff = (radius > (2 * line_h)) ? (radius - (2 * line_h)) : radius;
        final double young_bound = timeAxisYoungestAgeMa(); // >0 for a fossil-only tree: no boundaries past the tips
        final java.util.List<Double> ages = new java.util.ArrayList<>();
        for (final GeologicTimeScale.Interval iv : GeologicTimeScale.overlapping(GeologicTimeScale.bandRanks(root_age)[0],
                young_bound, root_age)) {
            final double b = iv.oldMa();
            if ((b > young_bound) && (b < root_age)) {
                ages.add(b);
            }
        }
        ages.sort((a, c) -> Double.compare(a, c)); // ascending age = younger boundary first = larger radius (outer)
        int last_label_y = Integer.MIN_VALUE / 2;
        for (final double b : ages) {
            final int rr = geologicRadiusPx(b, root_age, radius);
            if ((rr <= 0) || (rr >= outer_cutoff)) {
                continue;
            }
            final int ly = (cy - rr) + (fm.getAscent() / 2);
            if ((ly - last_label_y) < line_h) {
                continue;
            }
            final String label = TreePanelUtil.formatCompactNumber(b);
            final int w = fm.stringWidth(label);
            final int lx = cx - (w / 2);
            if (!bw) {
                g.setColor(new Color(bg.getRed(), bg.getGreen(), bg.getBlue(), 225));
                g.fillRect(lx - 1, ly - fm.getAscent(), w + 2, line_h);
            }
            g.setColor(ink);
            g.drawString(label, lx, ly);
            last_label_y = ly;
        }
        g.setFont(saved_font);
        g.setColor(saved);
    }

    /** Translucent fill for a circular geologic band: the official ICS colour softened with alpha so the tree + rings
     *  read over it (a light grey in a B&W export -- the band is a colour key, muted when colour is unavailable). */
    private Color geologicRingFill(final Color c, final boolean to_pdf, final boolean to_graphics_file) {
        if ((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) {
            return new Color(235, 235, 235);
        }
        return new Color(c.getRed(), c.getGreen(), c.getBlue(), GEOLOGIC_RING_ALPHA);
    }

    /** The "Time tree" badge text for the current tree, or null when it is not a time tree. A DATED tree (node
     *  {@code <date>}s) is auto-labeled -- with its declared unit if any; an ULTRAMETRIC tree is labeled only after
     *  the user confirmed the load-time offer. The expensive DATED detection is cached per tree. */
    String timeTreeBadgeLabel() {
        if (_time_tree_cached_for != _phylogeny) {
            _time_tree_cached_for = _phylogeny;
            _time_tree_dated = (AptxUtil.detectTimeTree(_phylogeny) == AptxUtil.TIME_TREE_KIND.DATED);
            _time_tree_unit = _time_tree_dated ? AptxUtil.timeTreeUnit(_phylogeny) : null;
        }
        if (_time_tree_dated) {
            return ForesterUtil.isEmpty(_time_tree_unit) ? "Time tree" : ("Time tree · " + _time_tree_unit);
        }
        if (_confirmed_time_tree_for == _phylogeny) { // confirmation is bound to THIS tree object
            return "Time tree";
        }
        return null;
    }

    /** Marks the panel's CURRENT (ultrametric) tree as a time tree after the user accepted the load-time offer, then
     *  repaints. The confirmation is tied to the tree object, so it lapses automatically if the tree is replaced. */
    void setConfirmedTimeTree(final boolean confirmed) {
        _confirmed_time_tree_for = confirmed ? _phylogeny : null;
        repaint();
    }

    boolean isConfirmedTimeTreeForTest() {
        return _confirmed_time_tree_for == _phylogeny;
    }

    // ---- time axis (geologic) ----------------------------------------------------------------------------------

    /** The absolute age (Ma before present) at the ROOT, used to calibrate the time axis: an explicit user
     *  calibration if set (&gt; 0), else the oldest node {@code <date>} value (a dated tree). 0 when the tree carries
     *  no absolute time information (an uncalibrated ultrametric tree). Cached per tree. */
    double timeAxisRootAgeMa() {
        // an explicit override applies only to the tree it was set on -- navigating into a subtree / undo / paste
        // replaces _phylogeny, and a stale override must not calibrate the now-different displayed tree
        if ((_time_axis_root_age > 0) && (_time_axis_root_age_for == _phylogeny)) {
            return _time_axis_root_age;
        }
        return maxNodeDateValue();
    }

    /** The maximum node {@code <date>} value in the tree, cached per tree. Its MEANING depends on the date convention:
     *  for a geologic tree (ages before present) it is the OLDEST age (the root age); for a calendar tip-dated tree
     *  (calendar-year dates) it is the MOST-RECENT tip's date (the present). 0 when the tree carries no dates. */
    private double maxNodeDateValue() {
        if (_time_axis_age_cached_for != _phylogeny) {
            _time_axis_age_cached_for = _phylogeny;
            double max = 0;
            if (_phylogeny != null) {
                for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if (n.getNodeData().isHasDate() && (n.getNodeData().getDate().getValue() != null)) {
                        max = Math.max(max, n.getNodeData().getDate().getValue().doubleValue());
                    }
                }
            }
            _time_axis_root_age_cache = max;
        }
        return _time_axis_root_age_cache;
    }

    /** The calendar date (decimal year) at the MOST-RECENT tip (the "present"), calibrating the calendar time axis: an
     *  explicit user override if set, else the largest node {@code <date>} value (a tip-dated tree whose dates are
     *  calendar years). 0 when the tree carries no absolute dates. */
    double timeAxisPresentDate() {
        if ((_time_axis_present_date > 0) && (_time_axis_present_date_for == _phylogeny)) {
            return _time_axis_present_date;
        }
        return maxNodeDateValue();
    }

    /** Sets an explicit most-recent-tip (present) calendar date for the calendar axis; 0 clears it. */
    void setTimeAxisPresentDate(final double year) {
        _time_axis_present_date = Math.max(0, year);
        _time_axis_present_date_for = _phylogeny; // scope the override to the tree it was set on
        repaint();
    }

    /** Sets an explicit root-age calibration (Ma) for a tree without dates (the "set root age" dialog); 0 clears it. */
    void setTimeAxisRootAge(final double ma) {
        _time_axis_root_age = Math.max(0, ma);
        _time_axis_root_age_for = _phylogeny; // scope the override to the tree it was set on (see timeAxisRootAgeMa)
        repaint();
    }

    /** The Time-Axis type in effect for THIS tree (tab): an explicit per-panel choice if one was made, else the type
     *  auto-derived from the tree's own {@code <date>} data. Per-tree, so a Dinosaur (geologic) tab and a SARS-CoV-2
     *  (calendar) tab can show different axes at the same time. */
    Options.TIME_AXIS_TYPE effectiveTimeAxisType() {
        return (_time_axis_type != null) ? _time_axis_type : derivedTimeAxisType();
    }

    /** The EXPLICIT per-tab Time-Axis type override, or {@code null} when the type follows auto-derive ("Auto"). Unlike
     *  {@link #effectiveTimeAxisType()} (which resolves null to the derived type) this returns the raw choice, so the
     *  Settings dropdown can show "Auto" vs an explicit Off/Geologic/Calendar. */
    Options.TIME_AXIS_TYPE getTimeAxisTypeOverride() {
        return _time_axis_type;
    }

    /** The auto-derived Time-Axis type for the current tree ({@link AptxUtil#deriveTimeAxisType}), cached by tree
     *  identity (recomputed when {@code _phylogeny} is replaced by a subtree / undo / paste). */
    private Options.TIME_AXIS_TYPE derivedTimeAxisType() {
        if (_derived_time_axis_for != _phylogeny) {
            _derived_time_axis_for = _phylogeny;
            _derived_time_axis_type = AptxUtil.deriveTimeAxisType(_phylogeny);
        }
        return _derived_time_axis_type;
    }

    /** Sets an explicit per-tab Time-Axis type override (Off / Geologic / Calendar); wins over auto-derive and survives
     *  a tree replacement within this panel (a tab-level view choice). */
    void setTimeAxisType(final Options.TIME_AXIS_TYPE type) {
        _time_axis_type = type;
        repaint();
    }

    boolean isTimeAxisGrid() {
        return _time_axis_grid;
    }

    void setTimeAxisGrid(final boolean b) {
        _time_axis_grid = b;
        repaint();
    }

    boolean isTimeAxisAges() {
        return _time_axis_ages;
    }

    void setTimeAxisAges(final boolean b) {
        _time_axis_ages = b;
        // the boundary-age row grows the geologic-axis reserve, so toggling it must re-fit (like the scale axis)
        if (geologicAxisApplies()) {
            getControlPanel().showWhole();
        }
        else {
            repaint();
        }
    }

    /** Restore a per-tree Time-Axis config read from a saved tree ({@code aptx:time_axis}); it wins over auto-derive.
     *  A {@code null} cfg (no / unparsable saved property) leaves this panel on auto-derive. */
    void applyTimeAxisConfig(final TimeAxisConfig cfg) {
        if (cfg == null) {
            return;
        }
        _time_axis_type = cfg.getType();
        _time_axis_grid = cfg.isGrid();
        _time_axis_ages = cfg.isAges();
        if (cfg.getRootAgeOverride() > 0) {
            setTimeAxisRootAge(cfg.getRootAgeOverride());
        }
        if (cfg.getPresentDateOverride() > 0) {
            setTimeAxisPresentDate(cfg.getPresentDateOverride());
        }
    }

    /** A snapshot of this panel's live effective Time-Axis config, for persisting into the tree on save. The two
     *  calibration overrides are included only when set AND still scoped to the current tree. */
    TimeAxisConfig currentTimeAxisConfig() {
        final double root_override = ((_time_axis_root_age > 0) && (_time_axis_root_age_for == _phylogeny))
                ? _time_axis_root_age : 0;
        final double present_override = ((_time_axis_present_date > 0) && (_time_axis_present_date_for == _phylogeny))
                ? _time_axis_present_date : 0;
        // snapshot the RAW type field (null == follow auto-derive), NOT the resolved effectiveTimeAxisType() -- so a
        // tree still on auto-derive persists no type, and a refinement-only (grid/ages) deviation keeps type auto
        return new TimeAxisConfig(_time_axis_type, root_override, present_override, _time_axis_grid, _time_axis_ages);
    }

    /** Reset the per-tab Time-Axis state to auto-derive (Reset to Defaults): clears the type override, both refinement
     *  toggles, and the calibration overrides. */
    void resetTimeAxisToAutoDerive() {
        _time_axis_type = null;
        _time_axis_grid = false;
        _time_axis_ages = false;
        _time_axis_root_age = 0;
        _time_axis_present_date = 0;
        repaint();
    }

    /** Persist this panel's Time-Axis config into the tree (an {@code aptx:time_axis} root property) so it survives a
     *  save/reload -- but ONLY when it DEVIATES from what auto-derive would produce (else strip any stale property, to
     *  keep files clean). Called at the save choke points. Pure view state -> does NOT {@code setEdited}. */
    void syncTimeAxisConfigToTree() {
        final TimeAxisConfig live = currentTimeAxisConfig();
        if (live.isDefault()) { // pure auto-derive (no type override, no toggles, no calibration override)
            TimeAxisConfig.writeToTree(_phylogeny, null); // clean default -> no property
        }
        else {
            TimeAxisConfig.writeToTree(_phylogeny, live);
        }
    }

    // ---- Auspice/Nextstrain time<->divergence display toggle (Increment 2) --------------------------------------

    /** True when the current tree carries BOTH a date (time signal) AND a {@code nextstrain:div} property (divergence
     *  signal), so the time&harr;divergence branch-length toggle is meaningful (an Auspice/Nextstrain tree). Cached by
     *  tree identity (recomputed on a tree replacement; an in-place node-data edit that adds/removes the metric is the
     *  same accepted cache-staleness class as the derived-time-axis / color-by caches). */
    boolean isNextstrainTimeDivergenceApplicable() {
        if ( _nextstrain_applicable_for != _phylogeny ) {
            _nextstrain_applicable = AuspiceJsonParser.hasTimeAndDivergence( _phylogeny );
            _nextstrain_applicable_for = _phylogeny;
        }
        return _nextstrain_applicable;
    }

    NEXTSTRAIN_BRANCH_MODE getNextstrainBranchMode() {
        return _nextstrain_branch_mode;
    }

    /** Switch the branch lengths (and the time axis) between the TIME view (num_date deltas + auto calendar axis) and
     *  the DIVERGENCE view (nextstrain:div deltas + no calendar axis -- the numeric substitutions/site scale). A pure,
     *  reversible display mode: both metrics stay retained on the tree, so it never {@code setEdited}s or checkpoints
     *  undo. Re-fits to the viewport (the depth axis changes meaning, so a fit -- not just a repaint -- is needed). */
    void setNextstrainBranchMode( final NEXTSTRAIN_BRANCH_MODE mode ) {
        if ( ( mode == null ) || ( mode == _nextstrain_branch_mode ) || !isNextstrainTimeDivergenceApplicable() ) {
            return;
        }
        _nextstrain_branch_mode = mode;
        if ( mode == NEXTSTRAIN_BRANCH_MODE.DIVERGENCE ) {
            AuspiceJsonParser.applyDivergenceBranchLengths( _phylogeny );
            _phylogeny.setDistanceUnit( "subs/site" );
            _time_axis_type = Options.TIME_AXIS_TYPE.NONE; // a calendar axis is meaningless for divergence lengths
        }
        else {
            AuspiceJsonParser.applyTimeBranchLengths( _phylogeny );
            _phylogeny.setDistanceUnit( "year" );
            _time_axis_type = null; // back to auto -> re-derives CALENDAR from the retained <date> unit
        }
        recalculateMaxDistanceToRoot(); // the branch-length change invalidates the depth cache
        if ( getControlPanel() != null ) {
            // showWhole() recomputes the layout (displayedPhylogenyMightHaveChanged) then fits to the VIEWPORT --
            // drift-free (the depth scale changes drastically between year deltas and div deltas).
            getControlPanel().showWhole();
        }
        else {
            repaint();
        }
        // this flipped the time axis (CALENDAR<->off); re-seed an open modeless Settings dialog so its axis combo isn't
        // stale (guarded internally by isShowing(), so a no-op when the dialog is closed)
        if ( ( getMainPanel() != null ) && ( getMainPanel().getMainFrame() != null ) ) {
            getMainPanel().getMainFrame().refreshOpenSettingsDialog();
        }
    }

    /** Reset the branch-length view to the TIME default (used by Reset to Defaults). Rewrites only THIS panel's model
     *  + invalidates its depth cache (NO cross-panel re-fit, so it is safe in a batch per-tab reset loop -- the layout
     *  refreshes when the tab is next shown). No-op if already TIME or the tree is not applicable. */
    void resetNextstrainBranchModeToDefault() {
        if ( _nextstrain_branch_mode == NEXTSTRAIN_BRANCH_MODE.TIME ) {
            return;
        }
        _nextstrain_branch_mode = NEXTSTRAIN_BRANCH_MODE.TIME;
        if ( isNextstrainTimeDivergenceApplicable() ) {
            AuspiceJsonParser.applyTimeBranchLengths( _phylogeny );
            _phylogeny.setDistanceUnit( "year" );
            _time_axis_type = null;
            recalculateMaxDistanceToRoot();
        }
        repaint();
    }

    /** Whether the branch lengths currently represent TIME, so a node's {@code <date>} interval (in time units) maps to
     *  an x / radius distance. FALSE only for an Auspice/Nextstrain tree switched to the DIVERGENCE view, where the
     *  branch lengths are substitutions/site: a date-based node-age (HPD) bar/spindle or fossil-range bar would then be
     *  scaled by the (huge) divergence corr and paint absurdly long (covering the whole canvas). Gates those overlays. */
    boolean isBranchLengthTimeCalibrated() {
        return _nextstrain_branch_mode != NEXTSTRAIN_BRANCH_MODE.DIVERGENCE;
    }

    /** Invalidate the identity-keyed time-axis caches (derived type + max node-date) after an IN-PLACE change to the
     *  node {@code <date>}s on the SAME tree object (e.g. Extract Dates from Labels) -- the caches are keyed by tree
     *  identity, which does not change on an in-place edit, so the Calendar axis would otherwise stay stale. */
    void invalidateTimeAxisDerivation() {
        _derived_time_axis_for = null;
        _time_axis_age_cached_for = null;
    }

    /** Whether the two-band geologic time axis is ON and drawable in a RECTANGULAR-family layout: mode GEOLOGIC, a
     *  phylogram in any of the three rectangular orientations (root-left / root-top / root-bottom), and an absolute
     *  root-age calibration. The circular analogue is {@link #geologicRingsApplyCircular()}; UNROOTED is N/A (an
     *  approved biological exception -- a distance-from-root radial has no single time axis to band). */
    boolean geologicAxisApplies() {
        return (effectiveTimeAxisType() == Options.TIME_AXIS_TYPE.GEOLOGIC)
                && getControlPanel().isDrawPhylogram()
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                && !breakLongBranchesActive() // a capped tree's age->x no longer lines up with the bands (see below)
                && (timeAxisRootAgeMa() > 0);
    }

    /** Whether the geologic time scale is ON and drawable as concentric coloured RINGS in the CIRCULAR layout: mode
     *  GEOLOGIC, a circular phylogram (radius encodes distance-from-root = time), and an absolute root-age calibration.
     *  The polar analogue of {@link #geologicAxisApplies()} -- the radial time axis is banded by the ICS periods. */
    boolean geologicRingsApplyCircular() {
        return (effectiveTimeAxisType() == Options.TIME_AXIS_TYPE.GEOLOGIC)
                && isCircularPhylogram()
                && !breakLongBranchesActiveCircular() // capped radius no longer lines up with the age rings
                && (timeAxisRootAgeMa() > 0);
    }

    /** Whether the CALENDAR (absolute-date) time axis is ON and drawable in a RECTANGULAR-family layout: mode CALENDAR,
     *  a phylogram in any of the three rectangular orientations, branch lengths present, and a most-recent-tip date
     *  calibration. A labeled year/decade ruler (like the numeric scale axis) rather than coloured bands. The polar
     *  analogue is {@link #calendarRingsApplyCircular()}; UNROOTED is N/A (an approved biological exception). */
    boolean calendarAxisApplies() {
        return (effectiveTimeAxisType() == Options.TIME_AXIS_TYPE.CALENDAR)
                && getControlPanel().isDrawPhylogram()
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                && !breakLongBranchesActive() // a capped tree's date->x no longer lines up with the year ticks
                && isPhyHasBranchLengths()
                && (getMaxDistanceToRoot() > 0)
                && (timeAxisPresentDate() > 0);
    }

    /** Whether the CALENDAR time axis is ON and drawable as concentric year RINGS in the CIRCULAR layout. */
    boolean calendarRingsApplyCircular() {
        return (effectiveTimeAxisType() == Options.TIME_AXIS_TYPE.CALENDAR)
                && isCircularPhylogram()
                && !breakLongBranchesActiveCircular() // capped radius no longer lines up with the year rings
                && (timeAxisPresentDate() > 0);
    }

    private int geologicAxisRowHeight() {
        return getFontMetrics(getTreeFontSet().getSmallFont()).getHeight() + 2;
    }

    /** Height of the two drawn geologic rows (coarse + fine band), WITHOUT the edge gap or the optional age-label row. */
    private int geologicAxisBandHeight() {
        return 2 * geologicAxisRowHeight();
    }

    /** Height of the optional boundary-age label row (drawn between the band and the edge when "Geologic Boundary Ages"
     *  is on), 0 otherwise. */
    private int geologicAgeRowHeight() {
        return isTimeAxisAges() ? geologicAxisRowHeight() : 0;
    }

    /** Height of the numeric "Ma before present" age ruler drawn at the outer edge of the geologic axis (baseline +
     *  ticks + one label row). Always present when the geologic axis applies -- it is the numeric axis itself. */
    private int geologicAgeRulerHeight() {
        return scaleAxisBandHeight();
    }

    /** The drawn content of the rectangular geologic axis: the two band rows + the optional age-label row + the numeric
     *  Ma ruler, WITHOUT the edge gap -- i.e. from the band's tree-side edge to its outermost drawn edge. */
    private int geologicAxisContentHeight() {
        return geologicAxisBandHeight() + geologicAgeRowHeight() + geologicAgeRulerHeight();
    }

    /** The full breadth reserved for the rectangular geologic axis: the drawn content ({@link #geologicAxisContentHeight})
     *  PLUS a few pixels of {@link #GEOLOGIC_AXIS_EDGE_GAP} between the axis and the canvas/window edge, so the axis
     *  doesn't sit flush against the border. Used by both the horizontal-bottom and vertical-side reserves and to place
     *  the band; the band is drawn at (edge - reserve), so the age row lands between the band and the edge. */
    private int geologicAxisReserve() {
        return geologicAxisContentHeight() + GEOLOGIC_AXIS_EDGE_GAP;
    }

    /** age (Ma) -> device x, aligned to the TREE's OWN branch scale: the root (oldest age = {@code root_age}) sits at
     *  {@code origin_x}, and each Ma spans exactly one branch-length unit ({@code corr} px), so the geologic band cell
     *  boundaries and the Ma ruler ticks line up with the branch nodes -- even when the tree has NO extant (age-0) tip.
     *  x = origin_x + (root_age - age)*corr. For a clock tree whose deepest tip is the present, root_age ==
     *  maxDistanceToRoot and this equals the older tip_x-anchored mapping; for a fossil-only tree (maxDist < root_age)
     *  it no longer wrongly pins age 0 to the youngest tip, so the bands align to the branches. */
    private static double ageToX(final double age, final double root_age, final float origin_x, final double corr) {
        return origin_x + ((root_age - age) * corr);
    }

    /** The age (Ma) of the YOUNGEST drawn tip = {@code root_age - maxDistanceToRoot}: 0 for a clock tree with an extant
     *  (present-day) tip, but &gt; 0 for a fossil-only tree whose most-recent taxon is still older than the present. The
     *  geologic axis / ruler / rings are drawn over {@code [youngestAge, root_age]} (mapped to [tip edge, root]) so they
     *  stay aligned to the branches and don't extend past the tips toward the present. */
    private double timeAxisYoungestAgeMa() {
        return Math.max(0, timeAxisRootAgeMa() - getMaxDistanceToRoot());
    }

    /** Test hook: the geologic axis's device x for a given age, using the current tree's root-age calibration and branch
     *  scale. The alignment invariant is that a consistent time tree's node sits at {@code ageToX(nodeAge)} == its
     *  drawn x; a fossil-only tree (maxDist &lt; rootAge) only satisfies it with the branch-aligned mapping. */
    double geologicAgeToXForTest(final double age) {
        return ageToX(age, timeAxisRootAgeMa(), _phylogeny.getRoot().getXcoord(), getXcorrectionFactor());
    }

    /** Test hook: {@link #timeAxisYoungestAgeMa()} (0 for a clock tree, &gt; 0 for a fossil-only tree). */
    double timeAxisYoungestAgeMaForTest() {
        return timeAxisYoungestAgeMa();
    }

    /** Draws the two-band colored geologic (ICS) time axis in the reserved bottom strip: the coarse rank over the fine
     *  rank (adaptive to depth -- Period/Epoch, Era/Period, or Eon/Era; see {@link GeologicTimeScale#bandRanks}), each cell
     *  filled with its official ICS colour and (where it fits) labelled, boundaries at the true ages. */
    private void paintGeologicTimeAxis(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file,
                                       final int graphics_file_y, final int graphics_file_height) {
        final double root_age = timeAxisRootAgeMa();
        final double corr = getXcorrectionFactor();
        if ((root_age <= 0) || (corr <= 0)) {
            return;
        }
        final float origin_x = _phylogeny.getRoot().getXcoord();
        final Rectangle vr = getVisibleRect();
        final int bottom = TreePanelUtil.scaleAxisFloatingBottom(to_pdf, to_graphics_file, graphics_file_y,
                graphics_file_height, getHeight(), vr.y + vr.height);
        final int row_h = geologicAxisRowHeight();
        final int top_y = bottom - geologicAxisReserve(); // band ends one GEOLOGIC_AXIS_EDGE_GAP short of the edge
        final Font saved_font = g.getFont();
        final Color saved_color = g.getColor();
        // the coarse+fine rank pair adapts to the tree's depth (Period/Epoch -> Era/Period -> Eon/Era)
        final GeologicTimeScale.Rank[] ranks = GeologicTimeScale.bandRanks(root_age);
        paintGeologicBand(g, ranks[0], root_age, origin_x, corr, top_y, row_h, to_pdf, to_graphics_file);
        paintGeologicBand(g, ranks[1], root_age, origin_x, corr, top_y + row_h, row_h, to_pdf, to_graphics_file);
        // optional Ma age labels at the coarse-band boundaries, in the reserved row between the band and the edge
        paintGeologicBoundaryAges(g, ranks[0], root_age, origin_x, corr, top_y + geologicAxisBandHeight(), to_pdf,
                to_graphics_file);
        // the numeric "Ma before present" ruler at the outer edge (the axis itself) -- ages increase toward the root
        paintGeologicAgeRuler(g, root_age, origin_x, corr,
                top_y + geologicAxisBandHeight() + geologicAgeRowHeight(), to_pdf, to_graphics_file);
        g.setFont(saved_font);
        g.setColor(saved_color);
    }

    /** The two-band geologic axis in a VERTICAL orientation (root-top / root-bottom): the SAME coloured band cells
     *  as the horizontal axis ({@link #paintGeologicBand}), but drawn INSIDE the R frame in LOGICAL coordinates,
     *  so the axis-aligned band cells AND their labels ride R into a side band down the breadth edge, in the reserve
     *  {@link #verticalScaleAxisReserve()} sets aside just past the last tip. Called while g is rotated by R (like the
     *  other vertical tree overlays), before the upright base frame is restored. */
    private void paintGeologicTimeAxisVertical(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        final double root_age = timeAxisRootAgeMa();
        final double corr = getXcorrectionFactor();
        if ((root_age <= 0) || (corr <= 0)) {
            return;
        }
        final float origin_x = _phylogeny.getRoot().getXcoord();
        final int row_h = geologicAxisRowHeight();
        // the reserved side band, just past the last tip; band ends one GEOLOGIC_AXIS_EDGE_GAP short of the breadth edge.
        // On screen it FLOATS to the viewport breadth edge (like the numeric vertical scale ruler) so it stays visible
        // when a zoomed tree is scrolled along the breadth; exports keep the tree-anchored position (WYSIWYG).
        final int band_top = floatVerticalGeologicBandTop(treeBreadthExtent() - geologicAxisReserve(), origin_x, to_pdf,
                to_graphics_file);
        final Font saved_font = g.getFont();
        final Color saved_color = g.getColor();
        final GeologicTimeScale.Rank[] ranks = GeologicTimeScale.bandRanks(root_age); // Period/Epoch -> Era/Period -> Eon/Era
        paintGeologicBand(g, ranks[0], root_age, origin_x, corr, band_top, row_h, to_pdf, to_graphics_file);
        paintGeologicBand(g, ranks[1], root_age, origin_x, corr, band_top + row_h, row_h, to_pdf, to_graphics_file);
        // optional Ma age labels at the coarse-band boundaries; ride R (logical coords) like the band cell labels
        paintGeologicBoundaryAges(g, ranks[0], root_age, origin_x, corr, band_top + geologicAxisBandHeight(), to_pdf,
                to_graphics_file);
        // the numeric "Ma before present" ruler at the outer edge; rides R (logical coords) like the band labels
        paintGeologicAgeRuler(g, root_age, origin_x, corr,
                band_top + geologicAxisBandHeight() + geologicAgeRowHeight(), to_pdf, to_graphics_file);
        g.setFont(saved_font);
        g.setColor(saved_color);
    }

    /** Draws the optional Ma age labels at the COARSE band's interval boundaries (e.g. "201.4" at the base of the
     *  Jurassic), on a baseline just past the band toward the edge, decimated so they never overlap. Reused by the
     *  horizontal axis (device coords) and the vertical axis (logical coords, riding R into rotated labels). */
    private void paintGeologicBoundaryAges(final Graphics2D g, final GeologicTimeScale.Rank coarse, final double root_age,
                                           final float origin_x, final double corr, final int band_bottom_y,
                                           final boolean to_pdf, final boolean to_graphics_file) {
        if (!isTimeAxisAges()) {
            return;
        }
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = g.getFontMetrics();
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        final int baseline_y = band_bottom_y + fm.getAscent() + 1;
        final double young_bound = timeAxisYoungestAgeMa(); // >0 for a fossil-only tree: no boundaries past the tips
        // each coarse interval's OLDER edge is a boundary (its younger edge is the older edge of the next interval);
        // collect the in-range boundaries as (x, age), then place left->right with a min-gap decimation
        final java.util.List<double[]> pts = new java.util.ArrayList<>();
        for (final GeologicTimeScale.Interval iv : GeologicTimeScale.overlapping(coarse, young_bound, root_age)) {
            final double b = iv.oldMa();
            if ((b > young_bound) && (b < root_age)) {
                pts.add(new double[] { ageToX(b, root_age, origin_x, corr), b });
            }
        }
        pts.sort((a, c) -> Double.compare(a[0], c[0]));
        int last_right = Integer.MIN_VALUE;
        for (final double[] pt : pts) {
            final String label = TreePanelUtil.formatCompactNumber(pt[1]);
            final int w = fm.stringWidth(label);
            final int left = (int) Math.round(pt[0]) - (w / 2);
            if (left >= (last_right + SCALE_AXIS_LABEL_GAP)) {
                g.drawString(label, left, baseline_y);
                last_right = left + w;
            }
        }
    }

    /** Draws the numeric "Ma before present" age ruler for the geologic axis: a baseline from the tips (age 0) to the
     *  root (the oldest age), a tick mark + Ma label at each "nice" round age (0, 50, 100 ... via
     *  {@link TreePanelUtil#maAxisTickValues}), the unit "Ma" at the age-0 end, decimated so labels never overlap. Ages
     *  INCREASE toward the root -- the standard age-before-present axis. Because the ruler is anchored to the tree's own
     *  root-age calibration ({@link #ageToX}), not a manual offset/reverse, it can't show a wrong root age (unlike the
     *  FigTree reverse-axis footgun). Reused by the horizontal axis (device coords) and the vertical axis (logical
     *  coords, riding R into a rotated ruler), exactly like {@link #paintGeologicBoundaryAges}. */
    private void paintGeologicAgeRuler(final Graphics2D g, final double root_age, final float origin_x,
                                       final double corr, final int ruler_y, final boolean to_pdf,
                                       final boolean to_graphics_file) {
        final double[] ticks = TreePanelUtil.maAxisTickValues(root_age);
        if (ticks.length == 0) {
            return;
        }
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = g.getFontMetrics();
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        final Stroke saved_stroke = g.getStroke();
        g.setStroke(STROKE_1);
        final double young_bound = timeAxisYoungestAgeMa(); // >0 for a fossil-only tree: the tip end is youngest, not 0
        // the tip end of the axis line: youngestAge maps here (== origin_x + maxDist*corr); for a clock tree that is age 0
        final int tip_x = (int) Math.round(ageToX(young_bound, root_age, origin_x, corr));
        final int baseline_y = ruler_y + SCALE_AXIS_TICK_LEN + fm.getAscent() + 1;
        // the axis line from the root (oldest age, origin_x) to the deepest tip (youngest age, tip_x)
        drawLine((int) Math.round(origin_x), ruler_y, tip_x, ruler_y, g);
        // ticks run right (youngest) to left (root age); collect (x, age) then place LEFT->RIGHT with a min-gap decimation.
        // a fossil-only tree's youngest tip is > 0 Ma, so skip ticks younger than it (they would fall past the tip end)
        final java.util.List<double[]> pts = new java.util.ArrayList<>();
        for (final double age : ticks) {
            if (age < (young_bound - 1e-9)) {
                continue;
            }
            pts.add(new double[] { ageToX(age, root_age, origin_x, corr), age });
        }
        pts.sort((a, c) -> Double.compare(a[0], c[0]));
        // a fossil-only tree's tip end sits at the youngest age (e.g. the K-Pg 66 Ma), which is usually not a round
        // tick -- so reserve its label span at the tip end and draw it explicitly, the feature's headline value. For a
        // clock tree (young_bound == 0) age 0 IS a round tick, so this block is skipped and the ruler is unchanged.
        final boolean label_tip = young_bound > 1e-9;
        final int tip_half = label_tip ? (fm.stringWidth(TreePanelUtil.formatCompactNumber(young_bound)) / 2) : 0;
        final int tip_reserve_left = label_tip ? ((tip_x - tip_half) - SCALE_AXIS_LABEL_GAP) : Integer.MAX_VALUE;
        int last_right = Integer.MIN_VALUE;
        for (final double[] pt : pts) {
            final int x = (int) Math.round(pt[0]);
            drawLine(x, ruler_y, x, ruler_y + SCALE_AXIS_TICK_LEN, g); // the tick mark (always drawn)
            final String label = TreePanelUtil.formatCompactNumber(pt[1]);
            final int half = fm.stringWidth(label) / 2;
            if (((x - half) >= (last_right + SCALE_AXIS_LABEL_GAP)) && ((x + half) <= tip_reserve_left)) {
                g.drawString(label, x - half, baseline_y);
                last_right = x + half;
            }
        }
        if (label_tip) { // the youngest-tip age at the tip end, always labelled (the round ticks reserved space for it)
            drawLine(tip_x, ruler_y, tip_x, ruler_y + SCALE_AXIS_TICK_LEN, g);
            g.drawString(TreePanelUtil.formatCompactNumber(young_bound), tip_x - tip_half, baseline_y);
            last_right = tip_x + tip_half;
        }
        // the unit "Ma" just past the tip end, if it clears the last drawn tick label
        final int unit_x = last_right + SCALE_AXIS_UNIT_GAP;
        if (unit_x >= (last_right + SCALE_AXIS_LABEL_GAP)) {
            g.drawString("Ma", unit_x, baseline_y);
        }
        g.setStroke(saved_stroke);
    }

    /** Shifts the vertical geologic band's logical breadth so its OUTER edge (away from the tips) floats to the viewport
     *  breadth edge on screen -- so the two-band axis stays visible when a zoomed vertical tree is scrolled along the
     *  breadth, exactly like the numeric vertical scale ruler. Exports/print keep the tree-anchored position (WYSIWYG).
     *  Works by mapping the desired device-x edge back to a logical breadth via the R inverse (depth positions, hence
     *  the tree-alignment along the time axis, are untouched). */
    private int floatVerticalGeologicBandTop(final int anchored_band_top, final float origin_x, final boolean to_pdf,
                                             final boolean to_graphics_file) {
        if (to_pdf || to_graphics_file || (_orientation_R_inverse == null)) {
            return anchored_band_top;
        }
        final Rectangle vr = getVisibleRect();
        if ((vr.width <= 0) || (vr.height <= 0)) {
            return anchored_band_top;
        }
        // the OUTERMOST drawn edge (past the band + the optional age-label row) is what floats to the viewport edge
        final int outer_ly = anchored_band_top + geologicAxisContentHeight(); // outermost drawn edge (larger breadth)
        final Point2D.Double outer_dev = screenPoint(origin_x, outer_ly);
        final double toward_tree_x = screenPoint(origin_x, outer_ly - 16.0).x; // 16 breadth toward the tips
        // the outer edge floats to the viewport edge on the side AWAY from the tree (a small gap in from the border)
        final boolean tree_to_right = toward_tree_x > outer_dev.x;
        final int target_x = tree_to_right ? (vr.x + GEOLOGIC_AXIS_EDGE_GAP)
                : ((vr.x + vr.width) - GEOLOGIC_AXIS_EDGE_GAP);
        final Point2D.Double floated = toLogicalPoint(target_x, (int) Math.round(outer_dev.y));
        return (int) Math.round(floated.y) - geologicAxisContentHeight();
    }

    private void paintGeologicBand(final Graphics2D g, final GeologicTimeScale.Rank rank, final double root_age,
                                   final float origin_x, final double corr, final int y, final int h,
                                   final boolean to_pdf, final boolean to_graphics_file) {
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = g.getFontMetrics();
        final Color ink = scaleInkColor(to_pdf, to_graphics_file);
        final double young_bound = timeAxisYoungestAgeMa(); // >0 for a fossil-only tree: don't draw past the tips
        for (final GeologicTimeScale.Interval iv : GeologicTimeScale.overlapping(rank, young_bound, root_age)) {
            final double young = Math.max(young_bound, iv.youngMa());
            final double old = Math.min(root_age, iv.oldMa());
            final int x_young = (int) Math.round(ageToX(young, root_age, origin_x, corr));
            final int x_old = (int) Math.round(ageToX(old, root_age, origin_x, corr));
            final int left = Math.min(x_old, x_young);
            final int w = Math.abs(x_young - x_old);
            if (w <= 0) {
                continue;
            }
            g.setColor(iv.color()); // the official ICS colour (kept even in B&W export -- the timescale IS a colour key)
            g.fillRect(left, y, w, h);
            g.setColor(ink);
            g.drawRect(left, y, w, h);
            if ((fm.stringWidth(iv.name()) + 4) <= w) {
                g.setColor(labelInkOn(iv.color()));
                g.drawString(iv.name(), left + ((w - fm.stringWidth(iv.name())) / 2), y + fm.getAscent() + 1);
            }
        }
    }

    /** Black or white label ink, whichever reads better on the given band colour (per-cell contrast). */
    private static Color labelInkOn(final Color c) {
        final double luminance = ((0.299 * c.getRed()) + (0.587 * c.getGreen()) + (0.114 * c.getBlue())) / 255.0;
        return (luminance > 0.55) ? Color.BLACK : Color.WHITE;
    }

    /** Draws the small "Time tree" badge at the top-right of the drawing region (viewport-fixed on screen, the export
     *  extent for a file), when the current tree is a time tree. WYSIWYG. */
    final private void paintTimeTreeBadge(final Graphics2D g,
                                          final int region_x,
                                          final int region_width,
                                          final int region_y_top,
                                          final boolean to_pdf,
                                          final boolean to_graphics_file) {
        final String label = timeTreeBadgeLabel();
        if (label == null) {
            return;
        }
        final Font saved = g.getFont();
        g.setFont(getTreeFontSet().getSmallFont());
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        final int tw = g.getFontMetrics().stringWidth(label);
        // top-right, but never left of the left margin (a zero-width region falls back to the left)
        final int x = Math.max(region_x + MOVE, (region_x + region_width) - tw - MOVE);
        final int y = region_y_top + g.getFontMetrics().getAscent() + 4;
        g.drawString(label, x, y);
        g.setFont(saved);
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
            y1 = Math.max(g.getFontMetrics().getHeight(), y1 - scaleAxisBottomReserve());
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
        float start_x = labelSegmentStartX(node.getXcoord(), effectiveNodeHalfBoxSize(node), x_shift);
        if ((getControlPanel().getTreeDisplayType() == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM)
                && node.isExternal()) {
            start_x = labelSegmentStartX((float) ((displayedMaxDistanceToRoot() * getXcorrectionFactor()) + TreePanel.MOVE
                    + getXdistance()), effectiveNodeHalfBoxSize(node), x_shift);
        }
        float start_y;
        if (!using_visual_font) {
            start_y = node.getYcoord() + (getFontMetricsForLargeDefaultFont().getAscent()
                    / (node.getNumberOfDescendants() == 1 ? 1 : 3.0f));
        } else {
            start_y = node.getYcoord()
                    + (getFontMetrics(g.getFont()).getAscent() / (node.getNumberOfDescendants() == 1 ? 1 : 3.0f));
        }
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
     * taxonomy and node-data segments.
     */
    static float labelSegmentStartX(final float base_x, final int half_box_size, final float prior_width) {
        return base_x + half_box_size + LABEL_GAP_AFTER_NODE_SHAPE + prior_width;
    }

    /** Half the node MARK size actually DRAWN for {@code node}: the per-node visual size when "Use Visual Styles" is on
     *  and the node carries one, else the global default node-shape size. Label positions offset by this so a large
     *  per-node node shape doesn't overlap its label (a plain default node uses the default, unchanged). */
    int effectiveNodeHalfBoxSize(final PhylogenyNode node) {
        float box_size = getOptions().getDefaultNodeShapeSize();
        if (shows(DisplayOption.USE_STYLE) && (node != null)) {
            final NodeVisualData vis = node.getNodeData().getNodeVisualData();
            if ((vis != null) && (vis.getSize() != NodeVisualData.DEFAULT_SIZE)) {
                box_size = vis.getSize();
            }
        }
        return (int) (box_size / 2.0f); // truncate (floor) to match the previous default int-division exactly
    }

    /** The largest node-mark half-box drawn among the external tips: the global default unless "Use Visual Styles" is
     *  on AND some tip carries a larger custom size. {@link #depthLabelReserve()} offsets the vertical tip-label reserve
     *  by this (not the bare default) so a large custom tip node can't push its upright label past the depth edge. The
     *  O(tips) walk runs only when visual styles are on; otherwise every tip resolves to the (constant) default. */
    private int maxTipEffectiveHalfBoxSize() {
        final int def = (int) (getOptions().getDefaultNodeShapeSize() / 2.0);
        if ((_phylogeny == null) || (getControlPanel() == null) || !shows(DisplayOption.USE_STYLE)) {
            return def;
        }
        int max = def;
        for (final PhylogenyNode t : _phylogeny.getExternalNodes()) {
            max = Math.max(max, effectiveNodeHalfBoxSize(t));
        }
        return max;
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

    /** The X of the radial (circular/unrooted) canvas centre. For a FULL export the drawing canvas is the passed
     *  graphics-file width -- a FIXED-SIZE export lays the tree into a canvas whose size differs from the on-screen
     *  panel, so centring at getWidth()/2 would push the ring off the export canvas. On screen, and for a visible-only
     *  export (whose Graphics is translated to crop), the drawing space is the panel's own width, so getWidth() is
     *  correct. (For the ordinary full export the caller passes getWidth() as the canvas, so the two agree.) */
    private int radialCanvasCenterX(final int graphics_file_width, final boolean to_pdf,
                                    final boolean to_graphics_file) {
        final boolean full_export = (to_pdf || to_graphics_file) && !getOptions().isGraphicsExportVisibleOnly();
        return (full_export ? graphics_file_width : getWidth()) / 2;
    }

    private int radialCanvasCenterY(final int graphics_file_height, final boolean to_pdf,
                                    final boolean to_graphics_file) {
        final boolean full_export = (to_pdf || to_graphics_file) && !getOptions().isGraphicsExportVisibleOnly();
        return (full_export ? graphics_file_height : getHeight()) / 2;
    }

    final private void paintUnrooted(final PhylogenyNode n,
                                     final double low_angle,
                                     final double high_angle,
                                     final boolean radial_labels,
                                     final Graphics2D g,
                                     final boolean to_pdf,
                                     final boolean to_graphics_file,
                                     final int graphics_file_width,
                                     final int graphics_file_height) {
        if (n.isRoot()) {
            n.setXcoord(radialCanvasCenterX(graphics_file_width, to_pdf, to_graphics_file));
            n.setYcoord(radialCanvasCenterY(graphics_file_height, to_pdf, to_graphics_file));
        }
        if (n.isExternal()) {
            paintNodeDataUnrootedCirc(g,
                    n,
                    to_pdf,
                    to_graphics_file,
                    radial_labels,
                    (high_angle + low_angle) / 2,
                    isInFoundNodes(n));
            // the tip's domain architecture rides its spoke, extending outward past the label (like circular) -- only
            // with RADIAL labels; under horizontal labels it would clash with the upright labels (see domainBoxesDrawn)
            if (radial_labels && (getControlPanel() != null) && shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES)) {
                final int num_ext = _phylogeny.getNumberOfExternalNodes();
                // tips sit near the periphery (radius ~ radialDiameter/2); estimate the arc between adjacent tips there
                final double spacing = (Math.PI * radialDiameter()) / Math.max(1, num_ext);
                final int height = TreePanelUtil.domainBoxHeight((float) spacing, DOMAIN_STRUCTURE_HEIGHT_MIN,
                        DOMAIN_STRUCTURE_HEIGHT_MAX);
                paintDomainArchitectureRadial(g, n, n.getXcoord(), n.getYcoord(), (high_angle + low_angle) / 2,
                        _length_of_longest_text_only + DOMAIN_RADIAL_GAP, height, to_pdf);
            }
            return;
        }
        // honor collapse: a collapsed clade is a single stub here -- its incoming branch + node box are drawn by the
        // PARENT's iteration (coords were set there), so we just stop, NOT recursing into the hidden subtree. Without
        // this, unrooted drew collapsed clades fully EXPANDED (circular already honors collapse), which disagreed with
        // hit-testing / halos / hover that treat the subtree as hidden. A bare stub for now (no label -- see
        // paintNodeDataUnrootedCirc); a radial collapse marker + [N] count is a planned follow-up.
        if (n.isCollapse()) {
            return;
        }
        // internal-node label (clade names, node/seq names) rides the branch radially, rotated to this node's own
        // wedge midpoint (= the direction of its incoming branch). The root sits at the canvas centre, too cramped
        // for a label, so it is skipped. Gated on "Show Internal Data" inside the method (not dynamic-hiding-culled,
        // matching the rectangular layout).
        if (!n.isRoot()) {
            paintNodeDataUnrootedCirc(g, n, to_pdf, to_graphics_file, radial_labels, (high_angle + low_angle) / 2,
                    isInFoundNodes(n));
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
                double dtp = desc.getDistanceToParent();
                if (dtp < 0) {
                    dtp = 0;
                } else if (breakLongBranchesActiveUnrooted() && (dtp > breakLongBranchCap())) {
                    dtp = breakLongBranchCap(); // "Break Long Branches": cap the spoke length (glyph drawn below)
                }
                length = (float) (dtp * getUrtFactor());
            } else {
                length = getUrtFactor();
            }
            final double mid_angle = current_angle + (arc_size / 2);
            final float new_x = (float) (x + (Math.cos(mid_angle) * length));
            final float new_y = (float) (y + (Math.sin(mid_angle) * length));
            desc.setXcoord(new_x);
            desc.setYcoord(new_y);
            paintUnrooted(desc, current_angle, current_angle + arc_size, radial_labels, g, to_pdf, to_graphics_file,
                    graphics_file_width, graphics_file_height);
            current_angle += arc_size;
            assignGraphicsForBranchWithColorForParentBranch(desc, false, g, to_pdf, to_graphics_file);
            drawLine(x, y, new_x, new_y, g);
            // "Break Long Branches": mark a capped spoke with a break glyph, rotated to the branch direction, at 0.72
            // along it (clear of the support/length numbers at the midpoint)
            if (breakLongBranchesActiveUnrooted() && (desc.getDistanceToParent() > breakLongBranchCap())) {
                paintBranchBreakGlyph(g, (float) (x + ((new_x - x) * BRANCH_BREAK_GLYPH_FRACTION)),
                        (float) (y + ((new_y - y) * BRANCH_BREAK_GLYPH_FRACTION)), mid_angle, to_graphics_file);
            }
            // support + branch-length numbers ride the middle of this branch, rotated to its direction (mid_angle)
            paintBranchDataRadial(g, desc, (x + new_x) / 2.0, (y + new_y) / 2.0, mid_angle, to_pdf, to_graphics_file);
            paintNodeBox(new_x, new_y, desc, g, to_pdf, to_graphics_file);
            if (desc.isCollapse()) {
                // collapsed clade-root stub (paintUnrooted returned early for it -> no subtree): draw its collapse
                // marker (triangle + count) opening outward along this branch's direction (mid_angle)
                paintRadialCollapsedMarker(g, desc, mid_angle, to_pdf, to_graphics_file);
            }
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
        if (n.isCollapse()) {
            return; // honor collapse in the thumbnail too: draw the clade as a stub, matching the main paintUnrooted
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
                double dtp = desc.getDistanceToParent();
                if (dtp < 0) {
                    dtp = 0;
                } else if (breakLongBranchesActiveUnrooted() && (dtp > breakLongBranchCap())) {
                    dtp = breakLongBranchCap(); // "Break Long Branches": cap the spoke here too, matching paintUnrooted
                }
                length = (float) (dtp * urt_ov_factor);
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

    /** The node's properties as one-line label text: the chosen fields' VALUES, comma-joined (internal aptx:*
     *  metadata, such as the persisted Re-import annotation profile, is never shown). May be empty. */
    private final String propertiesToString(final PhylogenyNode node) {
        return TreePanelUtil.labelPropertiesText(node.getNodeData().getProperties(), _label_property_refs,
                _annotation_column_refs);
    }

    private void setColor(final Graphics2D g,
                          final PhylogenyNode node,
                          final boolean to_graphics_file,
                          final boolean to_pdf,
                          final boolean is_in_found_nodes,
                          final Color default_color) {
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
        Color c;
        if (bw) {
            c = Color.BLACK;
        } else if (is_in_found_nodes) {
            c = getColorForFoundNode(node);
        } else if (isColorByProperty()) {
            c = getPropertyBasedColor(node);
        } else if (shows(DisplayOption.USE_STYLE) && (node.getNodeData().getNodeVisualData() != null)
                && (node.getNodeData().getNodeVisualData().getFontColor() != null)) {
            c = node.getNodeData().getNodeVisualData().getFontColor();
        } else if (getOptions().isColorLabelsSameAsParentBranch() && shows(DisplayOption.USE_STYLE)
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
        return (to_pdf || (to_graphics_file && getOptions().isExportBlackAndWhite())) ? Color.BLACK : default_color;
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
        if (shows(DisplayOption.USE_STYLE) && (node.getNodeData().getNodeVisualData() != null)) {
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
        // radial: scale the unrooted spread by the radial-zoom diameter (reserving the label margin on both sides so
        // the outermost labels fit the square), so it grows/shrinks in lockstep with the circular radius. Non-radial:
        // the urt factor is unused, keep the old viewport-based value.
        final int d;
        if (isRadialLayout()) {
            // reserve the label margin on both sides, but CAP it (like the circular radius) so a tree with long
            // labels relative to the canvas still fans out instead of collapsing to a tiny region. When domains are
            // drawn, reserve the FULL reach (labels + the bounded domain track) -- capped so the fan can't collapse --
            // so the domains fit rather than overflowing the canvas.
            final double half = radialDiameter() / 2.0;
            final double cap = domainBoxesDrawnInCurrentLayout() ? (1 - RADIAL_MIN_TREE_RATIO) : RADIAL_LABEL_MAX_RATIO;
            final int per_side = (int) Math.min(MOVE + getLongestExtNodeInfo(), half * cap);
            d = Math.max(MIN_RADIAL_DIAMETER, radialDiameter() - (2 * per_side));
        } else {
            d = getVisibleRect().width < getVisibleRect().height ? getVisibleRect().width : getVisibleRect().height;
        }
        if (isPhyHasBranchLengths() && getControlPanel().isDrawPhylogram()) {
            // "Break Long Branches": scale by the CAPPED radial max so the informative part fans out (matches
            // paintUnrooted, which caps each branch's spoke length; root-EXCLUDED, like getMaxDistanceToRoot)
            final double max = breakLongBranchesActiveUnrooted() ? breakCappedRadialMax() : getMaxDistanceToRoot();
            setUrtFactor((float) (d / (2 * max)));
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

    final void showNodeEditFrame(final PhylogenyNode n) { // package-visible so a test can open one the way the UI does
        if (_node_frame_index < TreePanel.MAX_NODE_FRAMES) {
            // Node-data edits ARE undoable ("Edit Node Data"), but the checkpoint is NOT taken here: writeBack
            // commits fields on every selection change and on close, so a checkpoint on mere open would push a
            // no-op undo -- and clear the redo stack -- for someone who only inspects a node. NodeEditPanel
            // instead snapshots on the first write that FOLLOWS a committed cell edit, so one editor visit is
            // exactly one undo step, and a visit that changes nothing leaves the history untouched.
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
     * Set parameters for painting/rendering the displayed tree
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
            final int axis_reserve = verticalScaleAxisReserve() + scaleAxisBottomReserve() + msaRulerReserve()
                    + msaConservationReserve();
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
            double height = _phylogeny.calculateHeight(!_options.isCollapsedWithAverageHeigh());
            //final double height = PhylogenyMethods.calculateMaxDepth( _phylogeny );
            // "Break Long Branches": derive the depth scale from the CAPPED height so the informative part reclaims the
            // width a broken outlier branch would otherwise consume (the branch itself is drawn capped, both here for
            // the scale and in calculateBranchLengthToParent for its x). With no branch over the cap this equals the
            // ordinary height, so a well-behaved tree is unchanged. Off / cladogram / aligned / radial: full height.
            if (breakLongBranchesActive() && (breakCappedHeight() > 0)) {
                height = breakCappedHeight();
            }
            // a subtree's root branch is drawn as a fixed stub (displayedRootBranchLength), not to scale, so it is NOT
            // part of the depth -- exclude it (BOTH calculateHeight AND breakCappedHeight fold in the (capped) root
            // branch) so the depth scale fills the width instead of leaving a stub-sized empty margin (mirrors
            // displayedTreeHeight; ov_corr below reads the same height)
            if (isCurrentTreeIsSubtree()) {
                height -= cappedRootBranchLength();
            }
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
        updateMsaScrollBar(); // the MSA window / scrollability may have changed with the new layout (or orientation)
    }

    final void calculateLongestExtNodeInfo() {
        if ((_phylogeny == null) || _phylogeny.isEmpty()) {
            return;
        }
        // When "Show External Data" is off, no external-node labels/domains are drawn (see the
        // isShowExternalData guards in paintNodeData), so reserve ZERO width for them -- otherwise fit-to-window ("F")
        // and fit-width ("W") leave a large empty right margin where the (undrawn) labels would have gone. This is the
        // single source of the label-reach reservation (depthLabelReserve/breadthLabelReserve, the radial label reach,
        // annotation-column x), so zeroing it here fixes the wasted space in EVERY display type. Node shapes at the
        // tips are covered by the MOVE margin the fit always keeps.
        if ((getControlPanel() != null) && !isShowExternalDataForThisTab()) {
            _longest_ext_node_info = 0;
            _length_of_longest_text = 0;
            _length_of_longest_text_only = 0;
            _longest_rendered_domain = 0;
            _ext_node_with_longest_txt_info = _phylogeny.getFirstExternalNode();
            return;
        }
        int max_possible_length = ForesterUtil
                .roundToInt((getSize().getWidth() - (2 * MOVE)) * AptxConstants.EXT_NODE_INFO_LENGTH_MAX_RATIO);
        if (max_possible_length < 20) {
            max_possible_length = 20;
        }
        int longest = 30;
        int longest_txt = 0;
        int longest_text_only = 0; // text label (name + taxonomy) width, WITHOUT the domain track width
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
            // null when the panel is not SHOWING (an off-screen render, a component not yet realized). Per-node
            // visual fonts cannot be measured then, so fall back to the default font rather than throwing.
            final Graphics2D g = (Graphics2D) getGraphics();
            if ((g != null) && shows(DisplayOption.USE_STYLE)) {
                use_vis = setFont(g, node);
            }
            final Font base = use_vis ? g.getFont() : getTreeFontSet().getLargeFont();
            sum = (use_vis ? getFontMetrics(g.getFont()) : getFontMetricsForLargeDefaultFont())
                    .stringWidth(sb.toString());
            if (has_tax) {
                sum += taxonomyLabelWidth(node.getNodeData().getTaxonomy(), base);
            }
            if (isRadialLayout()) {
                // radial labels are truncated to this reach (see paintNodeDataUnrootedCirc) -- so the reservation, the
                // domain start (radius + longest_text) and the drawn labels all agree, keeping everything on-canvas
                sum = Math.min(sum, radialMaxLabelWidth());
            }
            if (sum > longest_text_only) {
                longest_text_only = sum; // capture BEFORE the domain track width is added
            }
            if (!isRadialLayout()) {
                // reserve the tip-image slot (the label is shifted right/along-depth by tipImageAdvance) so annotation
                // columns and clade bands anchor PAST the image instead of overprinting it. Radial reserves the image's
                // spoke footprint separately (see paintNodeDataUnrootedCirc); 0 for a non-imaged/aligned/off tip.
                sum += tipImageAdvance(node);
            }
            if (shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES) && node.getNodeData().isHasSequence()
                    && (node.getNodeData().getSequence().getDomainArchitecture() != null)) {
                // FIXME
                // TODO this might need some clean up
                final DomainArchitecture d = node.getNodeData().getSequence().getDomainArchitecture();
                final double rendered = (effectiveDomainStructureWidth()
                        / ((RenderableDomainArchitecture) d).getOriginalSize().getWidth()) * d.getTotalLength();
                sum += rendered + 10;
                if (rendered > _longest_rendered_domain) {
                    _longest_rendered_domain = rendered;
                }
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
        setScaleDistance(TreePanelUtil.niceScaleBarDistance(getMaxDistanceToRoot()));
        setScaleLabel(scaleBarLabel(getScaleDistance()));
    }

    /** The scale-bar label for a given distance: the number plus the tree's distance unit in brackets, if any. */
    private String scaleBarLabel(final double distance) {
        String label = String.valueOf(distance);
        if (!ForesterUtil.isEmpty(_phylogeny.getDistanceUnit())) {
            label += " [" + _phylogeny.getDistanceUnit() + "]";
        }
        return label;
    }

    /** The scale-bar distance to DRAW: while Break Long Branches caps the tree, size the bar from the DRAWN (capped)
     *  extent so it reflects the un-broken (ingroup) scale, not the outlier-inflated one; otherwise the cached
     *  {@link #getScaleDistance()}. (The bar reads correctly for the un-broken tree; the broken branch is off-scale,
     *  marked by the break glyph -- the numeric scale AXIS / grid, which would span the whole width, stay suppressed.) */
    private double displayScaleDistance() {
        return breakLongBranchesActive() ? TreePanelUtil.niceScaleBarDistance(displayedMaxDistanceToRoot())
                : getScaleDistance();
    }

    private String displayScaleLabel() {
        return breakLongBranchesActive() ? scaleBarLabel(displayScaleDistance()) : getScaleLabel();
    }

    final Color calculateTaxonomyBasedColor(final Taxonomy tax) {
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
    // Ancestral-state pie charts: the selected discrete/geographic trait (null = off) and the stable state->color
    // map (all distinct states in the displayed tree) shared by every pie AND the legend, so a state keeps its
    // color across nodes. See TreePanelUtil.ancestralStateTraits / stateDistribution / collectAncestralStates.
    private String              _ancestral_pie_trait = null;
    private Map<String, Color>  _ancestral_pie_colors = null;
    // Per-node distribution cache (node identity -> its {state,prob} list), rebuilt with the color map on every
    // navigation/prune, so the per-frame paint + the dot-suppression guard do a map lookup instead of re-parsing
    // the BEAST brace strings on every repaint.
    private Map<PhylogenyNode, List<TreePanelUtil.StateProbability>> _ancestral_pie_dist = null;
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
    // A THIRD, independent legend for the ancestral-state pie charts (own draggable position + last-drawn bounds):
    // the state->color key for the pie wedges, shown alongside any color/rank/size legend.
    private Point               _ancestral_pie_legend_offset = null;
    private Rectangle           _ancestral_pie_legend_bounds = null;
    // A FOURTH, independent legend: the "Internal Taxonomy Key" -- a text key of the distinct internal-node taxa by
    // rank (no colors), shown alongside any other legend when isShowInternalTaxonomyKey() and the tree carries them.
    private Point               _internal_taxa_key_offset = null;
    private Rectangle           _internal_taxa_key_bounds = null;
    // A FIFTH, independent legend: the protein-domain legend (domain name -> its box colour), shown when the domain
    // label mode is LEGEND. E-value-cutoff aware -- it lists exactly the domains that pass the current threshold.
    private Point               _domain_legend_offset     = null;
    private Rectangle           _domain_legend_bounds     = null;
    // Which legend the active drag moves (they share _legend_grab_dx/dy). Generalized from a single "size?" boolean
    // once a third (pie) draggable legend was added.
    enum DRAGGED_LEGEND {
        PROPERTY, SIZE, ANCESTRAL_PIE, INTERNAL_TAXA, DOMAIN
    }
    // The two branch-length metrics an Auspice/Nextstrain tree can be laid out by (Increment 2 display toggle).
    enum NEXTSTRAIN_BRANCH_MODE {
        TIME( "Time" ), DIVERGENCE( "Divergence" );

        private final String _label;

        NEXTSTRAIN_BRANCH_MODE( final String label ) {
            _label = label;
        }

        String label() {
            return _label;
        }
    }
    private DRAGGED_LEGEND      _dragged_legend = DRAGGED_LEGEND.PROPERTY;
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
    // When SEVERAL ranks are annotated at once the legend is split into a titled block per rank -- otherwise the
    // reader gets one alphabetical heap and cannot tell a family row from a genus row. taxon -> its rank, plus the
    // rank order to print the blocks in (broadest first, the way taxonomy is written). Null for a single level,
    // which keeps the plain flat legend exactly as it was.
    private Map<String, String>  _rank_legend_sections = null;
    private java.util.List<String> _rank_legend_section_order = null;
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
    // Up to CladeLevel.MAX_LEVELS ranks can be annotated at once (genus inside family inside order), drawn as
    // nested bar/bracket columns. Split like the annotation columns: the SPECS are the user's choice and survive
    // navigation, the LEVELS are rebuilt from them for whatever tree is displayed now. Both are held in DRAWING
    // order (finest rank first = nearest the tips), decided once by CladeLevel.order.
    private java.util.List<CladeLevel.Spec> _clade_level_specs = null;
    private java.util.List<CladeLevel>      _clade_levels = null;
    private CLADE_VIS                 _clade_bands_mode = CLADE_VIS.BOXES;
    private boolean                   _clade_bands_skip_singletons = true; // BARS/BRACKETS: skip a single-tip clade

    /** On-screen angle of a clade bar/bracket taxon label in the root-left layout (root-top/bottom labels are always
     *  upright, circular labels ride the spoke). VERTICAL (90deg, reads up -- the compact default) can overlap when
     *  small clades sit close together; DIAGONAL (45deg) and HORIZONTAL (0deg, reads right) trade horizontal space
     *  for less vertical overlap. */
    enum CLADE_LABEL_ANGLE {
        VERTICAL, DIAGONAL, HORIZONTAL
    }
    private final static int          CLADE_BOX_ALPHA = 46;
    private final static int          CLADE_BAR_WIDTH = 9;
    private final static int          CLADE_BAR_GAP   = 16;
    // a few px of breathing room past the longest label, so the box/bar clears the text
    private final static int          CLADE_BAND_RIGHT_PAD = 6;
    private final static int          CLADE_BRACKET_TICK = 5;       // end-tick length of the monochrome "]" bracket
    private final static float        CLADE_BRACKET_STROKE = 1.5f;
    // Tip-aligned annotation columns (color strip / heat map / bar / text), drawn right of the labels.
    private java.util.List<AnnotationColumns.ColumnSpec> _annotation_column_specs = null; // the user's selection
    // which node properties the tip label shows, in order; null = every user-visible property that is not already
    // drawn as a column (see _annotation_column_refs) -- the default
    private java.util.List<String>     _label_property_refs = null;
    private java.util.Set<String>      _annotation_column_refs = java.util.Collections.emptySet();
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

    /** Over the ancestral-state pie legend -- guarded by isShowAncestralPies (pies + their legend draw in every
     *  layout, so no layout gate is needed). */
    final boolean isOnAncestralPieLegend(final MouseEvent e) {
        return isShowAncestralPies() && (_ancestral_pie_legend_bounds != null)
                && _ancestral_pie_legend_bounds.contains(e.getX(), e.getY());
    }

    /** Over the internal-taxonomy key -- guarded by the option; the recorded bounds mean it was drawn (has content). */
    final boolean isOnInternalTaxaKey(final MouseEvent e) {
        return getOptions().isShowInternalTaxonomyKey() && (_internal_taxa_key_bounds != null)
                && _internal_taxa_key_bounds.contains(e.getX(), e.getY());
    }

    /** Over the protein-domain legend -- guarded by the LEGEND label mode AND that domain boxes are actually shown; the
     *  recorded bounds mean it was drawn. */
    final boolean isOnDomainLegend(final MouseEvent e) {
        return (getOptions().getDomainLabelMode() == Options.DOMAIN_LABEL_MODE.LEGEND) && (getControlPanel() != null)
                && shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES) && (_domain_legend_bounds != null)
                && _domain_legend_bounds.contains(e.getX(), e.getY());
    }

    /** Over any legend (color/rank/column, size, ancestral pie, internal-taxonomy key, or domain legend). */
    final boolean isOnAnyLegend(final MouseEvent e) {
        return isOnPropertyLegend(e) || isOnSizeLegend(e) || isOnAncestralPieLegend(e) || isOnInternalTaxaKey(e)
                || isOnDomainLegend(e);
    }

    private Rectangle draggedLegendBounds() {
        switch (_dragged_legend) {
            case DOMAIN:
                return _domain_legend_bounds;
            case INTERNAL_TAXA:
                return _internal_taxa_key_bounds;
            case ANCESTRAL_PIE:
                return _ancestral_pie_legend_bounds;
            case SIZE:
                return _size_legend_bounds;
            default:
                return _property_legend_bounds;
        }
    }

    final void startLegendDrag(final MouseEvent e) {
        // the legends are drawn color/rank -> size -> pie -> internal-taxa (last on top), so an overlap grab must
        // hit-test in reverse draw order and give the top-most legend priority
        if (isOnDomainLegend(e)) {
            _dragged_legend = DRAGGED_LEGEND.DOMAIN;
        } else if (isOnInternalTaxaKey(e)) {
            _dragged_legend = DRAGGED_LEGEND.INTERNAL_TAXA;
        } else if (isOnAncestralPieLegend(e)) {
            _dragged_legend = DRAGGED_LEGEND.ANCESTRAL_PIE;
        } else if (isOnSizeLegend(e)) {
            _dragged_legend = DRAGGED_LEGEND.SIZE;
        } else {
            _dragged_legend = DRAGGED_LEGEND.PROPERTY;
        }
        final Rectangle b = draggedLegendBounds();
        if (b != null) {
            _legend_grab_dx = e.getX() - b.x;
            _legend_grab_dy = e.getY() - b.y;
            setCursor(MOVE_CURSOR);
        }
    }

    final void dragLegend(final MouseEvent e) {
        final Rectangle b = draggedLegendBounds();
        if (b == null) {
            return;
        }
        final Rectangle vp = getVisibleRect();
        int ox = (e.getX() - _legend_grab_dx) - vp.x;
        int oy = (e.getY() - _legend_grab_dy) - vp.y;
        ox = Math.max(0, Math.min(ox, Math.max(0, vp.width - b.width)));
        oy = Math.max(0, Math.min(oy, Math.max(0, vp.height - b.height)));
        switch (_dragged_legend) {
            case DOMAIN:
                _domain_legend_offset = new Point(ox, oy);
                break;
            case INTERNAL_TAXA:
                _internal_taxa_key_offset = new Point(ox, oy);
                break;
            case ANCESTRAL_PIE:
                _ancestral_pie_legend_offset = new Point(ox, oy);
                break;
            case SIZE:
                _size_legend_offset = new Point(ox, oy);
                break;
            default:
                _legend_offset = new Point(ox, oy);
                break;
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

    /** A click on the ancestral-pie legend: double-click returns it to its default corner (no recolorable rows). */
    final void handleAncestralPieLegendClick(final MouseEvent e) {
        if (e.getClickCount() == 2) {
            _ancestral_pie_legend_offset = null;
            repaint();
        }
    }

    /** Test hook: the last-drawn internal-taxonomy key bounds (null when it was not drawn / had no content). */
    Rectangle getInternalTaxaKeyBoundsForTest() {
        return _internal_taxa_key_bounds;
    }

    /** A click on the internal-taxonomy key: double-click returns it to its default corner (no recolorable rows). */
    final void handleInternalTaxaKeyClick(final MouseEvent e) {
        if (e.getClickCount() == 2) {
            _internal_taxa_key_offset = null;
            repaint();
        }
    }

    /** A click on the domain legend: double-click returns it to its default corner (no recolorable rows). */
    final void handleDomainLegendClick(final MouseEvent e) {
        if (e.getClickCount() == 2) {
            _domain_legend_offset = null;
            repaint();
        }
    }

    /** Test hook: the last-drawn domain-legend bounds (null when it was not drawn / had no content). */
    Rectangle getDomainLegendBoundsForTest() {
        return _domain_legend_bounds;
    }

    /** Test hook: draw the domain legend directly (records the bounds when {@code draggable}). */
    void drawDomainLegendForTest(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        drawDomainLegend(g, bounds, draggable);
    }

    /** Test hook: set the domain E-value threshold exponent (normally driven by the ControlPanel +/- controls). */
    void setDomainEvalueThresholdExpForTest(final int exp) {
        _domain_structure_e_value_thr_exp = exp;
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
        final String rank = rankLegendRankFor(taxon);
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

    /** The rank the current rank/clade legend is for (clade bands take precedence), or null. With several levels
     *  annotated the finest one owns the legend rows, which is where most of the taxa are. */
    private String currentRankLegendRank() {
        if (hasCladeBands()) {
            final java.util.List<String> ranks = cladeLevelRanks();
            return ranks.isEmpty() ? null : ranks.get(0);
        }
        return _branch_rank_colorize_rank;
    }

    /** The rank a legend ROW belongs to. With several levels each row has its own rank, so a colour chosen on a
     *  family row must be stored against "family" -- not against whichever rank happens to own the legend. */
    private String rankLegendRankFor(final String taxon) {
        if ((_rank_legend_sections != null) && _rank_legend_sections.containsKey(taxon)) {
            return _rank_legend_sections.get(taxon);
        }
        return currentRankLegendRank();
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
        if ((_legend_offset == null) && defaultLegendGoesLeft()) {
            // The never-dragged default is the top-RIGHT corner, which is exactly where the clade bar/bracket
            // columns are -- an annotated figure would open with its legend sitting on top of the marks it
            // describes (and, with several nested columns, hiding whole ranks of them). Start on the left
            // instead, above the root, where a root-left tree has room. Dragging still overrides this.
            return new Point(bounds.x + 10, bounds.y + 10);
        }
        return legendTopLeftFor(bounds, getVisibleRect(), _legend_offset, box_w, box_h);
    }

    /** Whether the default legend corner must move off the right edge: only when clade bar/bracket columns are
     *  actually drawn there, i.e. in the rectangular root-left layout (the vertical orientations put their marks
     *  along the top/bottom, and the radial layouts ring the tree, so neither collides with the corner). */
    private boolean defaultLegendGoesLeft() {
        return hasCladeBands() && (_clade_bands_mode != CLADE_VIS.BOXES) && !isRadialLayout()
                && !isVerticalOrientation() && (drawnCladeBandCount() > 0);
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
        _rank_legend_sections = null;
        _rank_legend_section_order = null;
        _branch_rank_colorize_rank = null; // branches are no longer rank-colorized
        repaint();
    }

    /** Whether a taxonomic-rank legend is currently available to show. */
    /** For tests: the rank legend's title ("Taxonomy: &lt;rank&gt;"), or null when there is no rank legend. Which
     *  rank the on-screen key actually describes is otherwise invisible from outside. */
    String rankLegendTitleForTest() {
        return _rank_legend_title;
    }

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
        setAncestralPieTrait( null ); // turns ancestral-state pies off, clears the pie legend pos
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

    /** Rebuilds the property-driven displays ("Color by", "Size by", ancestral-pie colors) from the currently
     *  displayed (visible) tree. Call this at every site the visible tips change (navigation, prune) so the
     *  encodings can never drift out of lockstep -- one method means a call site cannot rebuild one but forget
     *  another. */
    void rebuildPropertyDisplays() {
        rebuildPropertyColorScheme();
        rebuildPropertySizeScale();
        rebuildAncestralPieColors();
        _msa_conservation = null; // recomputed for the tips now on screen (entering a subtree, collapsing, pruning)
    }

    /** Show ancestral-state pie charts for the given discrete/geographic trait, or turn them off when
     *  {@code trait} is empty. Rebuilds the stable state->color map over the displayed tree. */
    void setAncestralPieTrait(final String trait) {
        _ancestral_pie_trait = ForesterUtil.isEmpty(trait) ? null : trait;
        if (_ancestral_pie_trait == null) {
            // turning pies OFF: forget the legend's dragged position + stale bounds so re-enabling shows it fresh
            _ancestral_pie_legend_offset = null;
            _ancestral_pie_legend_bounds = null;
        }
        rebuildAncestralPieColors();
        repaint();
    }

    String getAncestralPieTrait() {
        return _ancestral_pie_trait;
    }

    /** Recompute the ancestral-pie state->color map AND the per-node distribution cache for the currently displayed
     *  (visible) tree, so a state keeps a distinct color across every node; recomputed on navigation like the other
     *  property-driven encodings. */
    void rebuildAncestralPieColors() {
        if ((_ancestral_pie_trait == null) || (_phylogeny == null) || _phylogeny.isEmpty()) {
            _ancestral_pie_colors = null;
            _ancestral_pie_dist = null;
            return;
        }
        _ancestral_pie_colors = AptxUtil
                .assignDistinctColors(TreePanelUtil.collectAncestralStates(_phylogeny, _ancestral_pie_trait));
        // cache each node's distribution once, so the hot paint path does not re-parse the brace strings per frame
        final Map<PhylogenyNode, List<TreePanelUtil.StateProbability>> dist = new java.util.IdentityHashMap<>();
        for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext();) {
            final PhylogenyNode n = it.next();
            final List<TreePanelUtil.StateProbability> d = TreePanelUtil.stateDistribution(n, _ancestral_pie_trait);
            if (!d.isEmpty()) {
                dist.put(n, d);
            }
        }
        _ancestral_pie_dist = dist;
    }

    /** Whether ancestral-state pie charts are currently shown (a trait is selected and it yields >=1 state). */
    boolean isShowAncestralPies() {
        return (_ancestral_pie_trait != null) && (_ancestral_pie_colors != null) && !_ancestral_pie_colors.isEmpty();
    }

    /** Whether the current layout is a radial one (circular or unrooted), as opposed to the rectangular family.
     *  Used to suppress tip-mark legends (Color-by / Size-by dots, annotation columns) whose marks are not drawn
     *  radially. */
    boolean isRadialLayout() {
        return (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR);
    }

    /** Rotates the radial (circular/unrooted) layout by one step: clockwise increments the starting angle (screen y
     *  points down, so a larger angle sweeps the fan clockwise), counter-clockwise decrements it (wrapped into
     *  [0,2pi)). The single rotate seam behind the S/A keys, the Shift+mouse-wheel, and the X-/X+ zoom-cluster
     *  buttons (which become rotate controls in a radial layout, since X/Y zoom both drive the one radial diameter).
     *  A no-op outside a radial layout. */
    void rotateRadial(final boolean clockwise) {
        if (!isRadialLayout()) {
            return;
        }
        if (clockwise) {
            setStartingAngle((getStartingAngle() % TWO_PI) + ANGLE_ROTATION_UNIT);
        }
        else {
            setStartingAngle((getStartingAngle() % TWO_PI) - ANGLE_ROTATION_UNIT);
            if (getStartingAngle() < 0) {
                setStartingAngle(TWO_PI + getStartingAngle());
            }
        }
        getControlPanel().displayedPhylogenyMightHaveChanged(false);
    }

    /** Test hook: the assigned color for a pie state, or null when pies are off / the state is unknown. */
    Color ancestralPieColor(final String state) {
        return (_ancestral_pie_colors == null) ? null : _ancestral_pie_colors.get(state);
    }

    /** Test hook: the last-drawn ancestral-pie legend bounds (for drag/click hit-testing), or null. */
    Rectangle getAncestralPieLegendBounds() {
        return _ancestral_pie_legend_bounds;
    }

    /** Test hook: draws the ancestral-pie legend into {@code g} at {@code bounds} (draggable records its bounds; bw
     *  grays the swatches). */
    void drawAncestralPieLegendForTest(final Graphics2D g, final Rectangle bounds, final boolean draggable,
                                       final boolean bw) {
        drawAncestralPieLegend(g, bounds, draggable, bw);
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
        if (_rank_legend_sections != null) { // several ranks at once -> one titled block per rank
            drawSectionedRankLegend(g, bounds, draggable);
            return;
        }
        final Map<String, Color> values = orderedLegend(_rank_legend, _rank_legend_counts);
        drawCategoricalLegend(g, bounds, draggable, _rank_legend_title, values, _rank_legend_counts,
                _rank_legend.size() - values.size());
    }

    /**
     * The clade legend when SEVERAL ranks are annotated at once: one titled block per rank, broadest first (the
     * order taxonomy is written in), each block listing that rank's taxa. Without the blocks the three ranks
     * arrive as a single alphabetical heap in which a family row is indistinguishable from a genus row -- which
     * defeats the point of drawing the ranks together.
     * <p>
     * The per-rank cap keeps every block visible: with one shared cap a big rank would crowd the others out
     * entirely, so each gets its own share and the footer reports everything hidden across all of them. Header
     * rows are recorded as {@code null} labels in the row list, so the shared click-to-recolour hit test maps a
     * click on a header to "no row" instead of to the wrong taxon.
     */
    private void drawSectionedRankLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        final int swatch = 10, gap = 5, pad = 7, max_text_px = 200;
        g.setFont(legendFont());
        final FontMetrics fm = g.getFontMetrics();
        final int row_h = fm.getHeight() + 2;
        final int sections = Math.max(1, _rank_legend_section_order.size());
        final int per_section = (_legend_max_entries >= LEGEND_SHOW_ALL) ? LEGEND_SHOW_ALL
                : Math.max(3, _legend_max_entries / sections);
        // rows in draw order: a header (null label) then that rank's taxa, ordered by the shared sort toggle
        final java.util.List<String> labels = new ArrayList<>();
        final java.util.List<String> headers = new ArrayList<>();
        int hidden = 0;
        for (final String rank : _rank_legend_section_order) {
            final Map<String, Color> in_rank = new java.util.TreeMap<>();
            for (final Map.Entry<String, Color> e : _rank_legend.entrySet()) {
                if (rank.equals(_rank_legend_sections.get(e.getKey()))) {
                    in_rank.put(e.getKey(), e.getValue());
                }
            }
            if (in_rank.isEmpty()) {
                continue;
            }
            final Map<String, Color> shown = orderLegendEntries(in_rank, _rank_legend_counts, per_section,
                    _legend_sort_by_count);
            hidden += in_rank.size() - shown.size();
            headers.add(rank);
            labels.add(null); // the header row itself is not a clickable taxon
            labels.addAll(shown.keySet());
        }
        final boolean more_footer = hidden > 0;
        final String more_text = "… +" + hidden + " more";
        // Same rule as the flat legend: the expand/collapse control goes in the TITLE row. It matters even more
        // here -- a three-rank legend expanded to "show all" is long, and this renderer previously drew NO control
        // once everything was shown, so there was no way back to the short list at all.
        final int shown_rows = labels.size() - headings( labels );
        final String expand_chip = draggable ? legendExpandChip( shown_rows, hidden ) : null;
        // the sort order still drives the rows inside every block, so the control that changes it has to be here
        // too -- without it, annotating a second rank silently froze the legend on whatever order it was left in
        final boolean show_sort = draggable && (shown_rows >= 2);
        final String sort_chip = _legend_sort_by_count ? "[by count]" : "[A-Z]";
        int text_w = fm.stringWidth(_rank_legend_title)
                + ((expand_chip != null) ? (gap + fm.stringWidth(expand_chip)) : 0)
                + (show_sort ? (gap + fm.stringWidth(sort_chip)) : 0);
        for (final String h : headers) {
            text_w = Math.max(text_w, fm.stringWidth(h.toUpperCase()));
        }
        for (final String v : labels) {
            if (v != null) {
                text_w = Math.max(text_w, swatch + gap + fm.stringWidth(legendRowText(v, _rank_legend_counts, fm,
                        max_text_px)));
            }
        }
        if (more_footer) {
            text_w = Math.max(text_w, fm.stringWidth(more_text));
        }
        final int box_w = text_w + (2 * pad) + 4;
        final int box_h = ((1 + labels.size() + (more_footer ? 1 : 0)) * row_h) + (2 * pad);
        final Point tl = legendTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = clampLegendTop(bounds, tl.y, row_h, pad);
        if (draggable) {
            _legend_row_labels = labels; // nulls at the header rows -> a click there resolves to no taxon
            _legend_rows_top = y + pad + row_h;
            _legend_row_height = row_h;
            _legend_sort_toggle_bounds = null; // set below only when actually drawn
            _legend_more_bounds = null;
        }
        final Color fg = getTreeColorSet().getSequenceColor();
        int baseline = drawLegendBox(g, x, y, box_w, box_h, pad, _rank_legend_title, fm, draggable, false);
        g.setColor(fg);
        int chip_right = (x + box_w) - pad; // title-row chips are laid out right-to-left from the box edge
        if (expand_chip != null) { // in the title row, so it stays reachable however long the legend gets
            final int chip_w = fm.stringWidth(expand_chip);
            final int chip_x = chip_right - chip_w;
            g.drawString(expand_chip, chip_x, baseline);
            _legend_more_bounds = new Rectangle(chip_x - 2, y + pad, chip_w + 4, row_h);
            chip_right = chip_x - gap;
        }
        if (show_sort) {
            final int chip_w = fm.stringWidth(sort_chip);
            final int chip_x = chip_right - chip_w;
            g.drawString(sort_chip, chip_x, baseline);
            _legend_sort_toggle_bounds = new Rectangle(chip_x - 2, y + pad, chip_w + 4, row_h);
        }
        int header_i = 0;
        for (final String label : labels) {
            baseline += row_h;
            if (label == null) { // a rank heading, set apart from its rows by being bold and un-swatched
                g.setColor(fg);
                final Font plain = g.getFont();
                g.setFont(plain.deriveFont(Font.BOLD));
                g.drawString(headers.get(header_i++).toUpperCase(), x + pad, baseline);
                g.setFont(plain);
                continue;
            }
            g.setColor(_rank_legend.get(label));
            g.fillRect(x + pad + gap, baseline - fm.getAscent() + ((fm.getAscent() - swatch) / 2) + 1, swatch,
                    swatch);
            g.setColor(fg);
            g.drawString(legendRowText(label, _rank_legend_counts, fm, max_text_px), x + pad + gap + swatch + gap,
                    baseline);
        }
        if (more_footer) { // informative only; its control is in the title row
            baseline += row_h;
            g.setColor(fg);
            g.drawString(more_text, x + pad, baseline);
        }
    }

    /** How many of a sectioned legend's rows are rank HEADINGS rather than taxa (headings are null labels). */
    private static int headings(final java.util.List<String> labels) {
        int n = 0;
        for (final String label : labels) {
            if (label == null) {
                ++n;
            }
        }
        return n;
    }

    /** Test hook: whether the default (never-dragged) legend corner is moved off the right edge to clear the
     *  clade bar/bracket columns. */
    final boolean legendDefaultGoesLeftForTest() {
        return defaultLegendGoesLeft();
    }

    /** Test hook: the legend rows as drawn -- {@code null} where a rank heading sits. */
    final java.util.List<String> legendRowLabelsForTest() {
        return _legend_row_labels;
    }

    /**
     * The label of the legend's expand/collapse control, or null when neither applies: "[show all]" while entries
     * are hidden, "[show fewer]" once everything is shown and the list is longer than the default.
     */
    private String legendExpandChip(final int shown, final int hidden) {
        if (hidden > 0) {
            return "[show all]";
        }
        return ((_legend_max_entries >= LEGEND_SHOW_ALL) && (shown > DEFAULT_LEGEND_MAX_ENTRIES)) ? "[show fewer]"
                : null;
    }

    /**
     * Keeps the legend's TITLE row inside the drawing area however tall the box has grown. The expand/collapse
     * control lives in that row, and an expanded legend is routinely taller than the window -- with the control at
     * the BOTTOM (where it used to be) "show all" was a one-way trip on any long legend, because the way back was
     * scrolled off the screen.
     */
    private int clampLegendTop(final Rectangle bounds, final int y, final int row_h, final int pad) {
        final int lowest = (bounds.y + bounds.height) - row_h - (2 * pad);
        return Math.max(bounds.y, Math.min(y, Math.max(bounds.y, lowest)));
    }

    /**
     * Draws a boxed title plus swatch/label rows for an ordered value-to-color map. On screen
     * ({@code draggable}) the TITLE row carries both interactive controls -- a sort toggle
     * ("[by count]"/"[A-Z]") and the expand/collapse chip ("[show all]"/"[show fewer]") -- and their bounds are
     * recorded so the shared drag/hit-test machinery can map a click back to a row or control. Both sit in the
     * title row rather than the footer so they stay reachable when an expanded legend grows taller than the
     * window. A static export (not draggable) omits that clickable chrome but keeps the informative "+N more"
     * line so the figure still shows that categories were truncated. Shared by the property-color,
     * taxonomic-rank, and annotation-column legends. {@code counts} may be null (then rows show no count).
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
        final boolean more_footer = more > 0;                            // "+N more" info: on screen AND in exports
        final boolean show_sort = draggable && (shown >= 2);            // sort toggle: on-screen only
        // The expand/collapse control sits in the TITLE row, not the footer: an expanded legend is often taller
        // than the window, and a control at the bottom of it cannot be reached (see clampLegendTop).
        final String expand_chip = draggable ? legendExpandChip(shown, more) : null;
        final String sort_chip = _legend_sort_by_count ? "[by count]" : "[A-Z]";
        final String more_text = "… +" + more + " more";
        int text_w = fm.stringWidth(title) + (show_sort ? (gap + fm.stringWidth(sort_chip)) : 0)
                + ((expand_chip != null) ? (gap + fm.stringWidth(expand_chip)) : 0);
        for (final String v : values.keySet()) {
            text_w = Math.max(text_w, swatch + gap + fm.stringWidth(legendRowText(v, counts, fm, max_text_px)));
        }
        if (more_footer) {
            text_w = Math.max(text_w, fm.stringWidth(more_text));
        }
        // a few extra px on the right so the longest value clears the border even
        // when PDF font metrics run slightly wider than AWT's stringWidth().
        final int box_w = text_w + (2 * pad) + 4;
        final int box_h = ((1 + shown + (more_footer ? 1 : 0)) * row_h) + (2 * pad);
        final Point tl = legendTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = clampLegendTop(bounds, tl.y, row_h, pad);
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
        int chip_right = (x + box_w) - pad; // title-row chips are laid out right-to-left from the box edge
        if (expand_chip != null) { // expand/collapse, always reachable because the title row is always on screen
            final int chip_w = fm.stringWidth(expand_chip);
            final int chip_x = chip_right - chip_w;
            g.drawString(expand_chip, chip_x, baseline);
            _legend_more_bounds = new Rectangle(chip_x - 2, y + pad, chip_w + 4, row_h);
            chip_right = chip_x - gap;
        }
        if (show_sort) { // sort toggle, right-aligned in the title row (on-screen only -> draggable)
            final int chip_w = fm.stringWidth(sort_chip);
            final int chip_x = chip_right - chip_w;
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
        if (more_footer) { // "+N more": informative only now (also in exports); its control is in the title row
            baseline += row_h;
            g.drawString(more_text, x + pad, baseline);
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

    /**
     * Draws the ancestral-state pie legend: a boxed "Ancestral state: &lt;trait&gt;" title over one swatch+label row
     * per state, colored from the same stable {@link #_ancestral_pie_colors} map the wedges use (grayed by the same
     * {@link TreePanelUtil#grayShade} ramp when {@code bw}, so the key matches a black-and-white figure's wedges).
     * Modeled on the size legend (own draggable position, reuses the shared {@link #drawLegendBox} chrome); its rows
     * are not recolorable (auto-assigned palette), so it needs no sort/expand controls -- but it caps at
     * {@link #DEFAULT_LEGEND_MAX_ENTRIES} rows with a "+N more" line so a many-state trait can't run off the canvas.
     */
    private void drawAncestralPieLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable,
                                        final boolean bw) {
        if (!isShowAncestralPies()) {
            return;
        }
        final int pad = 7;
        final int swatch = 10;
        final int gap = 5;
        final int max_text = 200;
        g.setFont(legendFont());
        final FontMetrics fm = g.getFontMetrics();
        final int row_h = fm.getHeight() + 2;
        final int n = _ancestral_pie_colors.size();
        final int shown = Math.min(n, DEFAULT_LEGEND_MAX_ENTRIES);
        final int more = n - shown;
        final String title = clipToWidth("Ancestral state: " + _ancestral_pie_trait, fm, max_text);
        final String more_text = "… +" + more + " more";
        int text_w = fm.stringWidth(title);
        for (final String state : _ancestral_pie_colors.keySet()) {
            text_w = Math.max(text_w, swatch + gap + fm.stringWidth(clipToWidth(state, fm, max_text)));
        }
        if (more > 0) {
            text_w = Math.max(text_w, fm.stringWidth(more_text));
        }
        final int box_w = text_w + (2 * pad) + 4;
        final int box_h = ((1 + shown + ((more > 0) ? 1 : 0)) * row_h) + (2 * pad);
        final Point tl = ancestralPieLegendTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = tl.y;
        if (draggable) {
            _ancestral_pie_legend_bounds = new Rectangle(x, y, box_w, box_h);
        }
        final Stroke saved_stroke = g.getStroke();
        final Color fg = getTreeColorSet().getSequenceColor();
        // draggable=false here: this legend records its OWN bounds (above), not the property/size slot drawLegendBox
        // would otherwise write
        int baseline = drawLegendBox(g, x, y, box_w, box_h, pad, title, fm, false, false);
        int i = 0;
        for (final Map.Entry<String, Color> e : _ancestral_pie_colors.entrySet()) {
            if (i >= shown) {
                break;
            }
            baseline += row_h;
            // swatch by GLOBAL index so it matches the wedge (color, or the gray ramp in B&W)
            g.setColor(bw ? TreePanelUtil.grayShade(i, n) : e.getValue());
            g.fillRect(x + pad, baseline - fm.getAscent() + ((fm.getAscent() - swatch) / 2) + 1, swatch, swatch);
            g.setColor(fg);
            g.drawString(clipToWidth(e.getKey(), fm, max_text), x + pad + swatch + gap, baseline);
            ++i;
        }
        if (more > 0) {
            baseline += row_h;
            g.setColor(fg);
            g.drawString(more_text, x + pad, baseline);
        }
        g.setStroke(saved_stroke);
    }

    /** Default position of the ancestral-pie legend: TOP-RIGHT (the shared primary-legend corner) -- clear of the
     *  overview thumbnail (top-left by default) and the tree name + scale (bottom-left). If a color/rank/annotation
     *  legend already holds the top-right slot, drop to the BOTTOM-LEFT so the two don't collide. Once dragged,
     *  _ancestral_pie_legend_offset maps fractionally like the other legends (so exports honor the moved position). */
    private Point ancestralPieLegendTopLeft(final Rectangle bounds, final int box_w, final int box_h) {
        if (_ancestral_pie_legend_offset != null) {
            return legendTopLeftFor(bounds, getVisibleRect(), _ancestral_pie_legend_offset, box_w, box_h);
        }
        // which legend actually holds the top-right slot: the color/annotation legends are SUPPRESSED radially (see
        // the legend tail), so only the rank legend (which draws in every layout) occupies it there
        final boolean top_right_taken = hasRankLegend()
                || (!isRadialLayout() && isColorByProperty()) || annotationLegendVisible();
        if (top_right_taken) {
            return new Point(Math.max(bounds.x, bounds.x + 10),
                    Math.max(bounds.y, (bounds.y + bounds.height) - box_h - 10)); // bottom-left
        }
        return legendTopLeftFor(bounds, getVisibleRect(), null, box_w, box_h); // top-right
    }

    /**
     * Draws the "Internal Taxonomy Key": a boxed title over one row per RANK listing the distinct internal-node taxa
     * at that rank with counts, e.g. {@code order: Carnivora (3)  Rodentia (2)  Primates (1)}. A TEXT key (no colors,
     * so it does not duplicate the rank color legend); own draggable position, reusing {@link #drawLegendBox} chrome.
     * A no-op when the option is off or the tree carries no internal taxa (bounds nulled so there's no stale hit box).
     */
    private void drawInternalTaxonomyKey(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        if (!getOptions().isShowInternalTaxonomyKey() || (_phylogeny == null)) {
            return;
        }
        final java.util.LinkedHashMap<String, java.util.LinkedHashMap<String, Integer>> by_rank = TreePanelUtil
                .internalTaxaByRank(_phylogeny);
        if (by_rank.isEmpty()) {
            _internal_taxa_key_bounds = null;
            return;
        }
        final int pad = 7;
        final int max_text = 340;
        g.setFont(legendFont());
        final FontMetrics fm = g.getFontMetrics();
        final int row_h = fm.getHeight() + 2;
        final String title = "Internal taxa";
        final java.util.List<String> rows = new java.util.ArrayList<>();
        for (final Map.Entry<String, java.util.LinkedHashMap<String, Integer>> e : by_rank.entrySet()) {
            final StringBuilder sb = new StringBuilder(e.getKey()).append(": ");
            boolean first = true;
            for (final Map.Entry<String, Integer> t : e.getValue().entrySet()) {
                if (!first) {
                    sb.append("  ");
                }
                sb.append(t.getKey()).append(" (").append(t.getValue()).append(")");
                first = false;
            }
            rows.add(clipToWidth(sb.toString(), fm, max_text));
        }
        int text_w = fm.stringWidth(title);
        for (final String r : rows) {
            text_w = Math.max(text_w, fm.stringWidth(r));
        }
        final int box_w = text_w + (2 * pad) + 4;
        final int box_h = ((1 + rows.size()) * row_h) + (2 * pad);
        final Point tl = internalTaxaKeyTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = tl.y;
        if (draggable) {
            _internal_taxa_key_bounds = new Rectangle(x, y, box_w, box_h);
        }
        final Stroke saved_stroke = g.getStroke();
        int baseline = drawLegendBox(g, x, y, box_w, box_h, pad, title, fm, false, false);
        g.setColor(getTreeColorSet().getSequenceColor());
        for (final String r : rows) {
            baseline += row_h;
            g.drawString(r, x + pad, baseline);
        }
        g.setStroke(saved_stroke);
    }

    /** Default position of the internal-taxonomy key: BOTTOM-RIGHT (clear of the top-right primary legends, the
     *  top-left overview, and the bottom-left tree name / scale); once dragged, maps fractionally like the others. */
    private Point internalTaxaKeyTopLeft(final Rectangle bounds, final int box_w, final int box_h) {
        if (_internal_taxa_key_offset != null) {
            return legendTopLeftFor(bounds, getVisibleRect(), _internal_taxa_key_offset, box_w, box_h);
        }
        return new Point(Math.max(bounds.x, (bounds.x + bounds.width) - box_w - 10),
                Math.max(bounds.y, (bounds.y + bounds.height) - box_h - 10));
    }

    /**
     * Draws the protein-domain legend (domain name -> its box colour). Shown only in the LEGEND label mode, and
     * E-value-cutoff AWARE: it lists exactly the domain names that pass the current threshold across the displayed
     * tips -- so it stays in sync as the threshold moves. Its OWN draggable slot (default bottom-right); once dragged
     * it maps fractionally like the other legends. Painted last, so an overlap grab gives it priority.
     */
    private void drawDomainLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable) {
        if ((getOptions().getDomainLabelMode() != Options.DOMAIN_LABEL_MODE.LEGEND) || (_phylogeny == null)
                || !domainBoxesDrawnInCurrentLayout()) {
            // no boxes drawn in this layout -> no orphan legend (rectangular always; radial only with RADIAL labels)
            _domain_legend_bounds = null;
            return;
        }
        final LinkedHashMap<String, Color> values = new LinkedHashMap<String, Color>();
        final Map<String, Integer> counts = new HashMap<String, Integer>();
        collectDisplayedDomains(values, counts);
        if (values.isEmpty()) {
            _domain_legend_bounds = null;
            return;
        }
        final int pad = 7;
        final int swatch = 10;
        final int gap = 5;
        final int max_text = 240;
        g.setFont(legendFont());
        final FontMetrics fm = g.getFontMetrics();
        final int row_h = fm.getHeight() + 2;
        final String title = "Protein domains (E ≤ 1e" + _domain_structure_e_value_thr_exp + ")";
        int text_w = fm.stringWidth(title);
        for (final String name : values.keySet()) {
            text_w = Math.max(text_w, swatch + gap + fm.stringWidth(domainLegendRow(name, counts, fm, max_text)));
        }
        final int box_w = text_w + (2 * pad) + 4;
        final int box_h = ((1 + values.size()) * row_h) + (2 * pad);
        final Point tl = domainLegendTopLeft(bounds, box_w, box_h);
        final int x = tl.x;
        final int y = tl.y;
        if (draggable) {
            _domain_legend_bounds = new Rectangle(x, y, box_w, box_h);
        }
        final Stroke saved_stroke = g.getStroke();
        int baseline = drawLegendBox(g, x, y, box_w, box_h, pad, title, fm, false, false);
        final Color fg = getTreeColorSet().getSequenceColor();
        for (final Map.Entry<String, Color> en : values.entrySet()) {
            baseline += row_h;
            g.setColor(en.getValue());
            g.fillRect(x + pad, (baseline - fm.getAscent()) + ((fm.getAscent() - swatch) / 2) + 1, swatch, swatch);
            g.setColor(fg);
            g.drawString(domainLegendRow(en.getKey(), counts, fm, max_text), x + pad + swatch + gap, baseline);
        }
        g.setStroke(saved_stroke);
    }

    private String domainLegendRow(final String name, final Map<String, Integer> counts, final FontMetrics fm,
                                   final int max_text) {
        final Integer c = counts.get(name);
        return clipToWidth(name + ((c != null) ? (" (" + c + ")") : ""), fm, max_text);
    }

    /** Collect the DISPLAYED domain names passing the E-value threshold -> box colour, ordered by first appearance,
     *  plus per-name instance counts (the number of drawn boxes). */
    private void collectDisplayedDomains(final LinkedHashMap<String, Color> values, final Map<String, Integer> counts) {
        final double thr = Math.pow(10, _domain_structure_e_value_thr_exp);
        for (final PhylogenyNode n : _phylogeny.getExternalNodes()) {
            if (isHiddenUnderCollapse(n)) {
                continue; // a tip hidden under a collapsed clade draws no boxes, so it is not in the legend/counts
            }
            if (n.getNodeData().isHasSequence() && (n.getNodeData().getSequence().getDomainArchitecture() != null)) {
                final DomainArchitecture da = n.getNodeData().getSequence().getDomainArchitecture();
                for (final ProteinDomain d : da.getDomains().values()) {
                    if ((d.getName() != null) && (d.getConfidence() <= thr)) {
                        values.putIfAbsent(d.getName(), RenderableDomainArchitecture.colorFor(d.getName()));
                        counts.merge(d.getName(), 1, Integer::sum);
                    }
                }
            }
        }
    }

    /** Default position of the domain legend: BOTTOM-RIGHT (clear of the top-right primary legends and the top-left
     *  overview); once dragged, maps fractionally like the others. */
    private Point domainLegendTopLeft(final Rectangle bounds, final int box_w, final int box_h) {
        if (_domain_legend_offset != null) {
            return legendTopLeftFor(bounds, getVisibleRect(), _domain_legend_offset, box_w, box_h);
        }
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
            // Collapsing is display state, but it LIVES ON THE TREE (PhylogenyNode._collapse, which
            // copyNodeData() carries), so the snapshot history restores it like any other mutation -- and
            // restoreSnapshot already refreshes the collapse-derived caches.
            pushUndoCheckpoint(collapse ? "Collapse Clade" : "Uncollapse Clade");
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

    /** Whether {@code node} or anything under it is collapsed -- i.e. whether an "uncollapse" would change
     *  anything at all. Guards the undo checkpoint, so an uncollapse-all over an already-expanded tree does not
     *  push a snapshot that restores an identical tree (and does not clear the redo stack). */
    static boolean hasCollapsedNodeIn(final PhylogenyNode node) {
        if (node == null) {
            return false;
        }
        if (node.isCollapse()) {
            return true;
        }
        for (final PhylogenyNodeIterator it = new PreorderTreeIterator(node); it.hasNext();) {
            if (it.next().isCollapse()) {
                return true;
            }
        }
        return false;
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
            // Guard only the CHECKPOINT -- the refresh below (re-fit, scrollpane, repaint) must still run, as it
            // always has; it is only the undo entry that would restore an identical tree.
            if (hasCollapsedNodeIn(node)) {
                pushUndoCheckpoint("Uncollapse All");
            }
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
        _rank_legend_sections = null;
        _rank_legend_section_order = null;
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
            _rank_legend_sections = null; // one rank -> the plain flat legend
            _rank_legend_section_order = null;
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

    /** Annotates the tree's clades with boxes or bars at one {@code rank}; returns the number of bands. */
    final int setCladeBands(final String rank, final CLADE_VIS mode) {
        return setCladeBands(rank, mode, _clade_bands_skip_singletons, CLADE_LABEL_ANGLE.VERTICAL);
    }

    final int setCladeBands(final String rank, final CLADE_VIS mode, final boolean skip_singletons) {
        return setCladeBands(rank, mode, skip_singletons, CLADE_LABEL_ANGLE.VERTICAL);
    }

    /** One level, at the given label angle -- the single-rank form of {@link #setCladeLevels}. */
    final int setCladeBands(final String rank, final CLADE_VIS mode, final boolean skip_singletons,
                            final CLADE_LABEL_ANGLE label_angle) {
        return setCladeLevels(java.util.Collections.singletonList(new CladeLevel.Spec(rank, label_angle)), mode,
                              skip_singletons);
    }

    /**
     * Annotates the tree's clades at ONE TO {@link CladeLevel#MAX_LEVELS} ranks at once -- nested columns of bars
     * or brackets (genus inside family inside order). The specs are re-ordered finest-rank-first by
     * {@link CladeLevel#order}, so the caller never has to get the placement right; blanks and duplicate ranks
     * are dropped there too. Returns the number of marks that will actually be DRAWN across all levels (see
     * {@link #drawnCladeBandCount}).
     */
    final int setCladeLevels(final java.util.List<CladeLevel.Spec> specs, final CLADE_VIS mode,
                             final boolean skip_singletons) {
        _clade_level_specs = CladeLevel.order(specs);
        _clade_bands_mode = mode;
        _clade_bands_skip_singletons = skip_singletons;
        rebuildCladeBands();
        repaint();
        return drawnCladeBandCount();
    }

    /** The levels currently drawn, finest rank first (nearest the tips); never null. */
    private java.util.List<CladeLevel> cladeLevels() {
        return (_clade_levels == null) ? java.util.Collections.<CladeLevel> emptyList() : _clade_levels;
    }

    /** How many ranks are annotated at once. */
    final int cladeLevelCount() {
        return cladeLevels().size();
    }

    /** The ranks actually DRAWN, finest first -- for the report and the legend titles. A rank that placed no
     *  clade is not listed, so the report cannot claim a level the reader cannot see. */
    final java.util.List<String> cladeLevelRanks() {
        final java.util.List<String> out = new ArrayList<>();
        for (final CladeLevel level : drawableCladeLevels()) {
            out.add(level.getRank());
        }
        return out;
    }

    /**
     * BOXES stay SINGLE-level by design: the wash is alpha-composited, so a clade inside another clade would
     * paint as the product of two washes -- darker than either, and no longer the colour its legend row claims.
     * Bars and brackets are opaque marks in their own columns, so they nest cleanly. When boxes are somehow asked
     * for several levels, only the finest is drawn.
     */
    private java.util.List<CladeLevel> drawableCladeLevels() {
        // A level whose rank placed nothing must be dropped, not merely drawn empty: it would still claim a full
        // column of right margin (gap + mark + a label line, ~45 px) that nothing occupies, and would still be
        // counted in the "drew N marks across M levels" report. That happens for real -- ask for a rank the tree
        // cannot resolve alongside one it can, and only the resolvable level has bands.
        final java.util.List<CladeLevel> levels = new ArrayList<>();
        for (final CladeLevel level : cladeLevels()) {
            if (!level.getBands().isEmpty()) {
                levels.add(level);
            }
        }
        return ((_clade_bands_mode == CLADE_VIS.BOXES) && (levels.size() > 1)) ? levels.subList(0, 1) : levels;
    }

    /** A clade represented by a single tip -- a degenerate one-row bar/bracket. Counts the STRUCTURAL external
     *  descendants (collapse-independent), NOT getNumberOfExternalNodes(): the latter returns 1 for a DISPLAY-collapsed
     *  node, which would wrongly skip a real multi-tip clade's bar the moment the user collapses it (and would drift
     *  from the drawn-count computed at annotate time). getAllExternalDescendants() is a pure tree walk, so a collapsed
     *  clade still counts its real tips. */
    private boolean isSingleMemberClade(final CladeBand band) {
        // A band root that is no longer in the tree is skipped rather than walked: walking a detached node THROWS,
        // and this is reached from the paint (see the note on cladeBandYRange). The rebuild on every structural
        // change means it should never be stale -- this only makes the failure a missing bar, not a dead window.
        if (!isInCurrentTree(band.getRoot())) {
            return true;
        }
        return band.getRoot().getAllExternalDescendants().size() <= 1;
    }

    /** True when {@code band}'s bar/bracket is suppressed as a single-member clade. Only BARS and BRACKETS honor
     *  the option (a BOXES single-row wash is a harmless soft highlight, so boxes always draw). */
    private boolean skipCladeBand(final CladeBand band) {
        return _clade_bands_skip_singletons && (_clade_bands_mode != CLADE_VIS.BOXES) && isSingleMemberClade(band);
    }

    /** How many clade marks are actually DRAWN in the current mode -- the count of bands the render does NOT skip
     *  (single-member bars/brackets are skipped when the option is on), for the "drew N bar(s)/bracket(s)" report.
     *  Uses the SAME {@link #skipCladeBand} predicate the render loops do, so the count can't drift from what draws.
     *  The legend still lists every taxon. */
    private int drawnCladeBandCount() {
        int n = 0;
        for (final CladeLevel level : drawableCladeLevels()) {
            for (final CladeBand band : level.getBands()) {
                if (!skipCladeBand(band)) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Total clade bands built for the current rank -- every placed taxon, including single-member clades whose
     *  bar/bracket is skipped (they still get a legend row). Compare with {@link #drawnCladeBandCount()} (drawn marks)
     *  to tell "no tip could be placed" (0 total) apart from "all placed clades are single-member and were skipped"
     *  (total &gt; 0 but drawn 0). */
    final int cladeBandCount() {
        int n = 0;
        for (final CladeLevel level : cladeLevels()) {
            n += level.getBands().size();
        }
        return n;
    }

    /**
     * WRITES each clade's taxon at {@code rank} into the tree's internal {@code <taxonomy>} (the persistent
     * counterpart of {@link #setCladeBands}; see {@link TreePanelUtil#writeCladeTaxonomies}), appends a provenance
     * sentence, and marks the tree edited when anything was written. A tree-data MUTATION: the caller must
     * {@link #pushUndoCheckpoint} first (mirrors {@link #colorByRank}). Returns the number of nodes written.
     */
    final int writeCladeTaxonomiesByRank(final String rank, final boolean overwrite) {
        final int n = TreePanelUtil.writeCladeTaxonomies(_phylogeny, rank, TreePanelUtil.getDefaultLineageService(),
                overwrite);
        if (n > 0) {
            final String sentence = TreePanelUtil.cladeTaxonomyProvenance(_phylogeny, rank, n, overwrite);
            final String existing = _phylogeny.getDescription();
            _phylogeny.setDescription(ForesterUtil.isEmpty(existing) ? sentence : existing + " " + sentence);
            setEdited(true);
        }
        return n;
    }

    // --- accessors the figure spec captures from (see FigureSpec) ---------------------------------------------

    String getColorByPropertyRef() {
        return _color_by_property_ref;
    }

    String getSizeByPropertyRef() {
        return _size_by_property_ref;
    }

    boolean isCladeBandsSkipSingletons() {
        return _clade_bands_skip_singletons;
    }

    java.util.List<CladeLevel.Spec> getCladeLevelSpecs() {
        return _clade_level_specs;
    }

    /** The phylogram / aligned-phylogram / cladogram choice for THIS tab (already per-tab; see ControlPanel).
     *  Keyed on THIS panel's own index, not the current tab's -- "Save All" captures every tab's figure while only
     *  one of them is in front, and the current tab's choice is not this tab's. */
    Options.PHYLOGENY_DISPLAY_TYPE getTreeDisplayTypeForThisTab() {
        if ((getControlPanel() == null) || (getMainPanel() == null)) {
            return null;
        }
        int index = getMainPanel().getTreePanels().indexOf(this);
        if (index < 0) {
            index = getMainPanel().getCurrentTabIndex(); // not in a tab yet (being constructed)
        }
        return getControlPanel().treeDisplayTypeAt(index);
    }

    /** Stamps THIS tab's figure onto its own tree, so a saved file reproduces what is on screen. Called from the
     *  save paths beside {@link #syncTimeAxisConfigToTree()}, which is the same kind of embedded per-tree state. */
    void syncFigureToTree() {
        if (_phylogeny != null) {
            FigureSpec.writeToTree(_phylogeny, FigureSpec.capture(this));
        }
    }

    /** For tests: how many clade LEVELS are configured, regardless of whether any band could actually be placed.
     *  The specs are stored even when the tree resolves nothing at that rank, so this is the honest way to assert
     *  that a reset dropped the configuration on a fixture whose taxonomy is not resolvable offline. */
    int cladeLevelSpecCountForTest() {
        return (_clade_level_specs == null) ? 0 : _clade_level_specs.size();
    }

    final void clearCladeBands() {
        _clade_levels = null;
        _clade_level_specs = null;
        if (_branch_rank_colorize_rank == null) {
            clearRankLegend(); // nothing else owns the color key -> it goes with the bands
        }
        else {
            // A branch rank-colorization still owns the legend, but the CLADE bands overwrote its contents
            // (updateCladeBandLegend replaces the rows and the title). Leaving those rows behind would show a key
            // for marks that are gone, at a rank the branches are not colored by -- and a color picked on one of
            // them would be stored against the BRANCH rank, because currentRankLegendRank() falls back to it once
            // the bands are gone. Re-derive the branch legend instead.
            recolorBranchesByRank(_branch_rank_colorize_rank);
        }
    }

    final boolean hasCladeBands() {
        for (final CladeLevel level : cladeLevels()) {
            if (!level.getBands().isEmpty()) {
                return true;
            }
        }
        return false;
    }

    CLADE_VIS getCladeBandsMode() {
        return _clade_bands_mode;
    }

    /** Recomputes the bands from the current tree (cache-only); call after navigation swaps the tree. */
    final void rebuildCladeBands() {
        if ((_clade_level_specs == null) || _clade_level_specs.isEmpty() || (_phylogeny == null)
                || _phylogeny.isEmpty()) {
            _clade_levels = null;
            return; // leave any existing (property/rank-colorize) legend untouched
        }
        // Counts are gathered PER LEVEL and merged first-wins (the specs run finest-first), not accumulated into
        // one shared map: a taxon name that is assigned at two annotated ranks would otherwise have its two tip
        // counts SUMMED, inflating its "(N)" in the legend and skewing the by-count sort.
        final Map<String, Integer> counts = new HashMap<>();
        final java.util.List<java.util.List<CladeBand>> per_level = new ArrayList<>();
        for (final CladeLevel.Spec spec : _clade_level_specs) {
            final Map<String, Integer> level_counts = new HashMap<>();
            // each level runs its OWN maximal-monophyletic assignment at its rank; nesting is therefore whatever
            // the taxonomy actually says, not something imposed -- a clade whose broader rank is unresolved simply
            // has no mark outside it, which is honest about the gap rather than inventing a parent
            per_level.add(TreePanelUtil.cladeBands(_phylogeny, spec.getRank(),
                    TreePanelUtil.getDefaultLineageService(), rankOverridesFor(spec.getRank()), level_counts));
            for (final Map.Entry<String, Integer> e : level_counts.entrySet()) {
                counts.putIfAbsent(e.getKey(), e.getValue()); // finest level first -> the finer count wins
            }
        }
        // Several levels: re-colour them so each finer clade's colour sits inside its containing clade's hue slice
        // (a single level has no hierarchy to express and keeps the plain distinct-colour palette). Runs BEFORE the
        // per-taxon overrides are re-applied below, so a colour the user chose by hand still wins.
        final java.util.List<java.util.List<CladeBand>> colored = CladeHuePalette.hueBand(per_level);
        final java.util.List<CladeLevel> levels = new ArrayList<>();
        for (int i = 0; i < _clade_level_specs.size(); ++i) {
            final CladeLevel.Spec spec = _clade_level_specs.get(i);
            levels.add(new CladeLevel(spec, applyCladeColorOverrides(colored.get(i), spec.getRank())));
        }
        _clade_levels = levels;
        updateCladeBandLegend(counts);
    }

    /** Re-applies this rank's per-taxon user colour overrides after hue-banding, so a colour picked by hand on a
     *  legend row is never silently replaced by a computed one. */
    private java.util.List<CladeBand> applyCladeColorOverrides(final java.util.List<CladeBand> bands,
                                                               final String rank) {
        final Map<String, Color> overrides = rankOverridesFor(rank);
        if (overrides.isEmpty()) {
            return bands;
        }
        final java.util.List<CladeBand> out = new ArrayList<>();
        for (final CladeBand band : bands) {
            final Color chosen = overrides.get(band.getTaxon());
            out.add((chosen == null) ? band : new CladeBand(band.getTaxon(), chosen, band.getRoot()));
        }
        return out;
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
        final java.util.List<CladeLevel> levels = drawableCladeLevels();
        final Map<String, String> sections = new HashMap<>();
        final java.util.List<String> section_order = new ArrayList<>();
        for (final CladeLevel level : levels) {
            for (final CladeBand band : level.getBands()) {
                legend.putIfAbsent(band.getTaxon(), band.getColor());
                // levels run FINEST first, so first-wins really does file a taxon under the finer of two ranks
                // (a plain put would have let the broader, later level overwrite it -- and then a colour chosen
                // on that legend row would have been stored against the wrong rank)
                sections.putIfAbsent(band.getTaxon(), level.getRank());
            }
            section_order.add(0, level.getRank()); // levels are finest-first; print the blocks broadest-first
        }
        if (levels.size() > 1) {
            _rank_legend_sections = sections;
            _rank_legend_section_order = section_order;
        }
        else {
            _rank_legend_sections = null;
            _rank_legend_section_order = null;
        }
        if (!legend.isEmpty()) {
            _rank_legend = legend;
            _rank_legend_counts = counts;
            _rank_legend_title = "Taxonomy: " + String.join(" / ", cladeLevelRanks());
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
        rebuildAnnotationColumnRefs();
        rebuildAnnotationColumns();
        labelTextChanged(); // a field that became a column leaves the DEFAULT label set -> the labels got shorter
    }

    /** Re-derives the set of refs drawn as columns, which the default label selection excludes (one role per field). */
    private void rebuildAnnotationColumnRefs() {
        if (_annotation_column_specs == null) {
            _annotation_column_refs = java.util.Collections.emptySet();
            return;
        }
        final java.util.Set<String> refs = new java.util.HashSet<String>();
        for (final AnnotationColumns.ColumnSpec s : _annotation_column_specs) {
            if (s._type != AnnotationColumns.Type.LABEL) {
                refs.add(s._ref);
            }
        }
        _annotation_column_refs = refs;
    }

    /** The last annotation import applied to this tree (its source + column mapping), for one-click Re-import, or null. */
    org.forester.archaeopteryx.tools.NodeDataImporter.ImportProfile getLastImportProfile() {
        return _last_import_profile;
    }

    void setLastImportProfile(final org.forester.archaeopteryx.tools.NodeDataImporter.ImportProfile profile) {
        _last_import_profile = profile;
    }

    /** Which node properties the tip label shows, in display order, or null for the default (all of them). */
    java.util.List<String> getLabelPropertyRefs() {
        return _label_property_refs;
    }

    /**
     * Sets which node properties the tip label shows and in what order (the fields set to "In tip label" in
     * Tools > Annotation Fields). An EMPTY list means "no properties in the label" -- distinct from null, which
     * restores the default of showing every user-visible property. Either way the "Properties" display checkbox
     * still has to be on for any of them to be drawn.
     */
    void setLabelPropertyRefs(final java.util.List<String> refs) {
        _label_property_refs = (refs == null) ? null : new java.util.ArrayList<String>(refs);
        labelTextChanged();
    }

    /**
     * The tip-label TEXT changed, so the cached longest-label width has to be RECOMPUTED. Not merely zeroed: 0 is a
     * legitimate value rather than a sentinel (an all-external-data-hidden tab reserves exactly 0), so a caller that
     * only repaints would reserve no label width at all -- and that one number feeds the fit, the annotation-column
     * x, the clade-band right edge and the radial radius, so the columns would be drawn over the labels.
     */
    private void labelTextChanged() {
        if ((_phylogeny != null) && !_phylogeny.isEmpty()) {
            calculateLongestExtNodeInfo();
        }
    }

    /** For tests: the label text the renderer would draw for {@code node}, assembled by the very method the paint
     *  path uses -- so the "Properties" display gate, the chosen label fields and their order are all exercised. */
    String nodeLabelTextForTest(final PhylogenyNode node) {
        final StringBuilder sb = new StringBuilder();
        nodeDataAsSB(node, sb);
        return sb.toString();
    }

    void clearAnnotationColumns() {
        _annotation_column_specs = null;
        _annotation_columns = null;
        _annotation_col_widths = null;
        _focused_annotation_column = -1;
        // the same chooser assigns the label fields, so a reset drops those too (back to "all properties")
        _label_property_refs = null;
        _annotation_column_refs = java.util.Collections.emptySet();
        labelTextChanged(); // both halves of the reset change the label text
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
            // Any open node window points at a node of the tree we are ABOUT to replace. Left open it would edit a
            // detached node: the user's changes would vanish with no feedback while still marking the file dirty.
            closeAllNodeFrames();
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
            if (getControlPanel() != null) {
                // the ladderize toggle icon is view state, not in the snapshot -- re-derive it from the restored tree
                getControlPanel().syncOrderButtonIconToTree(this);
            }
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
                final RenderableDomainArchitecture rds = new RenderableDomainArchitecture(
                        node.getNodeData().getSequence().getDomainArchitecture());
                node.getNodeData().getSequence().setDomainArchitecture(rds);
            }
        }
    }

    private void notifyEditMenu() {
        if ((getMainPanel() != null) && (getMainPanel().getMainFrame() != null)) {
            getMainPanel().getMainFrame().updateEditMenu();
        }
    }

    /** Strip / symbol (categorical color key), heat map (color key) and bar (length/range key) columns have a
     *  legend; text columns do not. */
    private boolean columnHasLegend(final int i) {
        final AnnotationColumns.Type t = _annotation_columns.getColumn(i).getType();
        return (t == AnnotationColumns.Type.COLOR_STRIP) || (t == AnnotationColumns.Type.SYMBOL)
                || (t == AnnotationColumns.Type.HEATMAP) || (t == AnnotationColumns.Type.MATRIX)
                || (t == AnnotationColumns.Type.BAR) || AnnotationColumns.isMergedType(t);
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

    /** Whether the annotation-column legend is actually drawn in the current layout: the columns render as tip-aligned
     *  columns (rectangular family) OR concentric rings (circular), so the legend is shown in both; it is suppressed in
     *  the UNROOTED layout, where the columns/rings aren't drawn. */
    private boolean annotationLegendVisible() {
        return hasAnnotationColumnLegend()
                && (!isRadialLayout() || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR));
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
        if (AnnotationColumns.isMergedType(c.getType())) {
            // a merged (stacked-bar / pie) column has no PropertyColorScheme -- its legend is the series colour key,
            // ready whenever it has at least one series
            return !_annotation_columns.stackColors(col).isEmpty();
        }
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

    /** Toggles the color key (legend) for the annotation column whose header (rectangular) or ring (circular) was
     *  clicked. */
    final boolean handleAnnotationHeaderClick(final MouseEvent e) {
        final int col = (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                ? circularAnnotationRingAt(e.getX(), e.getY())
                : annotationHeaderColumnAt(e.getX(), e.getY());
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

    /** Test hook: the number of RESOLVED annotation columns (a STACKED_BAR group of several fields counts as one). */
    int annotationColumnCountForTest() {
        return hasAnnotationColumns() ? _annotation_columns.size() : 0;
    }

    /** Test hook: the horizontal label advance (px) a tip image reserves for {@code node} (0 for a non-imaged tip). */
    int tipImageAdvanceForTest(final PhylogenyNode node) {
        return tipImageAdvance(node);
    }

    /** Test hook: synchronously wait (off the EDT) until every imaged tip's image has loaded or failed, so a
     *  subsequent render is deterministic (the loads are otherwise asynchronous). */
    void preloadTipImagesForTest(final long timeout_ms) {
        if (_phylogeny == null) {
            return;
        }
        final File base = imageBaseDir();
        final long end = System.currentTimeMillis() + timeout_ms;
        for (final PhylogenyNode t : _phylogeny.getExternalNodes()) {
            final String ref = TipImages.imageRefFor(t);
            if (ref == null) {
                continue;
            }
            while ((tipImageCache().get(ref, base) == null) && !tipImageCache().isFailed(ref, base)
                    && (System.currentTimeMillis() < end)) {
                try {
                    Thread.sleep(20);
                }
                catch (final InterruptedException e) {
                    Thread.currentThread().interrupt();
                    break;
                }
            }
        }
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
        if (AnnotationColumns.isMergedType(col.getType())) {
            // a merged (stacked-bar / pie) column has no scheme -- its key is the series colours (name + swatch)
            noteLegendSubject("series:" + col.getHeader());
            drawSeriesLegend(g, bounds, draggable, col.getHeader(), col_i);
            return;
        }
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

    /**
     * Legend for a focused MERGED (stacked-bar / pie) annotation column: the series colour key -- one swatch + field
     * name per series, in series order, reusing the categorical-legend renderer (no per-series counts). A stacked bar's
     * absolute-vs-normalized mode is a property of the whole column, not of the key, so it is not shown here.
     */
    private void drawSeriesLegend(final Graphics2D g, final Rectangle bounds, final boolean draggable,
                                  final String title, final int col) {
        final java.util.List<Color> colors = _annotation_columns.stackColors(col);
        final java.util.List<String> headers = _annotation_columns.stackHeaders(col);
        final Map<String, Color> values = new java.util.LinkedHashMap<String, Color>(); // series order preserved
        for (int k = 0; k < Math.min(colors.size(), headers.size()); ++k) {
            // guarantee a UNIQUE legend key so two series whose display names collide (refs differing only by
            // namespace, e.g. data:count / sample:count -> "Count") each keep their own row + colour instead of one
            // overwriting the other in the map
            String label = headers.get(k);
            if (values.containsKey(label)) {
                label = label + " (" + (k + 1) + ")";
            }
            values.put(label, colors.get(k)); // no per-series count -> null
        }
        drawCategoricalLegend(g, bounds, draggable, title, values, null, 0);
    }

    /** Total horizontal space the annotation columns occupy (0 when none), including the gaps around them.
     * Zero for circular/unrooted layouts, where the columns are not drawn, so no width is reserved. */
    private int annotationColumnsWidth() {
        if (!hasAnnotationColumns() || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                || (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.UNROOTED)) {
            return 0;
        }
        return sumAnnotationColumnWidths();
    }

    /** The total depth (px) of the annotation-column stack: a leading gap plus each column's width and trailing gap.
     *  Layout-agnostic (callers gate on layout); shared by the rectangular {@link #annotationColumnsWidth()} and the
     *  circular {@link #circularAnnotationRingsReserve()} so the two can never disagree on the stack size. */
    private int sumAnnotationColumnWidths() {
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
            case STACKED_BAR:
                return Math.max(36, h * 4); // a touch wider than a single bar so the segments are legible
            case PIE:
                return Math.max(20, h * 2); // a square-ish cell big enough for legible wedges
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
                    case SYMBOL: {
                        // a centered shape glyph: FILLED when present, a hollow OUTLINE for an explicitly
                        // false/absent value, nothing for a missing value -- colored like a COLOR_STRIP cell
                        final AnnotationColumns.Fill fill = _annotation_columns.symbolFill(t, i);
                        if (fill != AnnotationColumns.Fill.NONE) {
                            final Color c = _annotation_columns.cellColor(t, i);
                            if (c != null) {
                                drawSymbolGlyph(g, c, xi + (w / 2.0f), t.getYcoord(), Math.min(w, cell_h) - 2,
                                        fill == AnnotationColumns.Fill.FILLED, _annotation_columns.symbolShape(i));
                            }
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
                    case STACKED_BAR:
                        paintStackedBarRow(g, _annotation_columns.stackFractions(t, i),
                                _annotation_columns.stackColors(i), xi, cy, w, cell_h);
                        break;
                    case PIE:
                        drawPieGlyph(g, _annotation_columns.stackFractions(t, i), _annotation_columns.stackColors(i),
                                xi + (w / 2.0f), t.getYcoord(), Math.min(w, cell_h) - 2, fg);
                        break;
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
        final double text_angle = (getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_BOTTOM)
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
                    case SYMBOL: {
                        // draw the glyph UPRIGHT in the base frame (NOT under R): a circle/square/diamond is
                        // symmetric under R's 90 deg, but a TRIANGLE would ride the rotation and point sideways.
                        // So place it at the R-mapped device centre and draw it un-rotated -- parity with the
                        // rectangular/circular glyph (apex up in every orientation), like the TEXT cells are anchored
                        // to the upright frame.
                        final AnnotationColumns.Fill fill = _annotation_columns.symbolFill(t, i);
                        if (fill != AnnotationColumns.Fill.NONE) {
                            final Color c = _annotation_columns.cellColor(t, i);
                            if (c != null) {
                                final Point2D.Double gp = screenPoint(xi + (w / 2.0), t.getYcoord());
                                g.setTransform(_orientation_base_transform);
                                drawSymbolGlyph(g, c, (float) gp.x, (float) gp.y, Math.min(w, cell_h) - 2,
                                        fill == AnnotationColumns.Fill.FILLED, _annotation_columns.symbolShape(i));
                                g.setTransform(withR);
                            }
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
                    case STACKED_BAR:
                        // axis-aligned segment rects ride R for free (a 90deg rotation keeps them axis-aligned),
                        // just like the BAR / COLOR_STRIP cells -- so a logical stacked column becomes a horizontal band
                        paintStackedBarRow(g, _annotation_columns.stackFractions(t, i),
                                _annotation_columns.stackColors(i), xi, cy, w, cell_h);
                        break;
                    case PIE: {
                        // draw the pie UPRIGHT in the base frame (its wedges would rotate under R) at the R-mapped
                        // device centre -- parity with the SYMBOL glyph and the rectangular/circular pie
                        final Point2D.Double gp = screenPoint(xi + (w / 2.0), t.getYcoord());
                        g.setTransform(_orientation_base_transform);
                        drawPieGlyph(g, _annotation_columns.stackFractions(t, i), _annotation_columns.stackColors(i),
                                (float) gp.x, (float) gp.y, Math.min(w, cell_h) - 2, fg);
                        g.setTransform(withR);
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

    /**
     * Draws a SYMBOL-column glyph: a shape (circle / square / diamond / triangle) centered at (center_x,
     * center_y) inscribed in a {@code diameter}-side box (floored so it never vanishes), either FILLED (solid) or
     * a hollow OUTLINE, in {@code color}. Shared by the rectangular, vertical, and circular annotation-column
     * paint paths -- each passes DEVICE coordinates and draws the glyph UPRIGHT (the vertical caller maps the
     * logical centre through R and draws in the base frame, so a triangle points up in every orientation).
     * Self-contained: the outline branch saves and restores g's stroke.
     */
    private void drawSymbolGlyph(final Graphics2D g, final Color color, final float center_x, final float center_y,
                                 final float diameter, final boolean filled,
                                 final AnnotationColumns.SymbolShape shape) {
        final float d = (float) Math.max(3.0, diameter);
        final float x = center_x - (d / 2.0f);
        final float y = center_y - (d / 2.0f);
        g.setColor(color);
        final Stroke saved_stroke = filled ? null : g.getStroke();
        if (!filled) {
            g.setStroke(STROKE_1);
        }
        switch (shape) {
            case SQUARE:
                _rectangle.setFrame(x, y, d, d);
                if (filled) { g.fill(_rectangle); } else { g.draw(_rectangle); }
                break;
            case DIAMOND:
                setDiamond(x, y, d, d);
                if (filled) { g.fill(_diamond); } else { g.draw(_diamond); }
                break;
            case TRIANGLE:
                _polygon.reset();
                _polygon.moveTo(x + (d / 2.0f), y); // apex, top-centre
                _polygon.lineTo(x + d, y + d);      // base-right
                _polygon.lineTo(x, y + d);          // base-left
                _polygon.closePath();
                if (filled) { g.fill(_polygon); } else { g.draw(_polygon); }
                break;
            case CIRCLE:
            default:
                if (filled) { drawOvalFilled(x, y, d, d, g); } else { drawOval(x, y, d, d, g); }
                break;
        }
        if (!filled) {
            g.setStroke(saved_stroke);
        }
    }

    /**
     * Draws one tip's STACKED_BAR row: the segments ({@code fractions}, parallel to {@code colors}) tiled
     * left-to-right across the column, each filled in its series colour, with cumulative rounding at the boundaries
     * so adjacent segments abut without a seam or overlap. A zero-length segment (missing / non-positive series
     * value) is skipped. Shared by the rectangular and vertical annotation-column paths -- both fill axis-aligned
     * rects (device coords, and logical coords riding R, respectively), so a stacked column becomes a horizontal
     * band in the vertical orientations for free.
     */
    private void paintStackedBarRow(final Graphics2D g, final double[] fractions, final java.util.List<Color> colors,
                                    final int xi, final int cy, final int w, final int cell_h) {
        double off = 0;
        for (int k = 0; k < fractions.length; ++k) {
            final double f = fractions[k];
            if (f > 0) {
                final int x0 = xi + (int) Math.round(off * w);
                final int x1 = xi + (int) Math.round((off + f) * w);
                // a sub-pixel segment rounds to zero width and is simply invisible (a composition slice, not a
                // single BAR whose minimum must show a stub); NO forced 1px, so it can't overpaint the next segment
                // or, at the column edge, poke past into the gap
                if (x1 > x0) {
                    g.setColor(colors.get(k));
                    g.fillRect(x0, cy, x1 - x0, cell_h);
                }
            }
            off += f;
        }
    }

    /**
     * Draws one tip's PIE glyph: the series ({@code fractions}, parallel to {@code colors}) as wedges of a disc of
     * diameter {@code diameter} centred at ({@code cx},{@code cy}), starting at 12 o'clock and sweeping clockwise, then
     * a thin {@code outline} circle around it (mirroring the ancestral-state pies via the shared {@code _arc}/
     * {@code _ellipse}). A tip with no data (all fractions 0) draws NOTHING (no empty circle). The fractions are the
     * tip's own proportions (they sum to 1 for a tip with data), so the wedge angles are {@code fraction * 360}. Shared
     * by the rectangular, vertical, and circular paint paths -- each passes DEVICE coordinates and draws upright.
     */
    private void drawPieGlyph(final Graphics2D g, final double[] fractions, final java.util.List<Color> colors,
                              final float cx, final float cy, final float diameter, final Color outline) {
        final float d = Math.max(3.0f, diameter);
        final double x = cx - (d / 2.0), y = cy - (d / 2.0);
        double start = 90.0; // begin at 12 o'clock and sweep clockwise (negative arc angle), like the ancestral pies
        boolean any = false;
        for (int k = 0; k < fractions.length; ++k) {
            final double f = fractions[k];
            if (f > 0) {
                final double sweep = 360.0 * f;
                g.setColor(colors.get(k));
                _arc.setArc(x, y, d, d, start, -sweep, Arc2D.PIE);
                g.fill(_arc);
                start -= sweep;
                any = true;
            }
        }
        if (any) { // outline only a real pie -- an empty circle would look like a tip with (blank) data
            final Stroke saved_stroke = g.getStroke();
            g.setStroke(STROKE_1);
            g.setColor(outline);
            _ellipse.setFrame(x, y, d, d);
            g.draw(_ellipse);
            g.setStroke(saved_stroke);
        }
    }

    /** The lazily-created tip-image cache; on a completed background load it invalidates the label-width reservation
     *  and repaints, so a freshly-loaded image (and its label shift) appears on the next frame. */
    private TipImageCache tipImageCache() {
        if (_tip_image_cache == null) {
            _tip_image_cache = new TipImageCache();
            _tip_image_cache.setRepaintCallback(() -> {
                calculateLongestExtNodeInfo(); // the tips' label reach may change once images load
                repaint();
            });
        }
        return _tip_image_cache;
    }

    /** The directory local (relative) image paths resolve against: the loaded tree file's directory, or null (then a
     *  relative path resolves against the working directory). Keep the tree, the annotation table, and the images
     *  together. */
    private File imageBaseDir() {
        return (_treefile != null) ? _treefile.getParentFile() : null;
    }

    /** Whether {@code node} is an external tip that carries an image reference and tip images are shown. Cheap
     *  (no image loading) -- used by the label offset and the paint loop. */
    private boolean hasTipImage(final PhylogenyNode node) {
        return getOptions().isShowTipImages() && (node != null) && node.isExternal()
                && (TipImages.imageRefFor(node) != null);
    }

    /** The fixed slot width (px) reserved for a tip image = the target height times the max aspect. */
    private int tipImageSlotWidth() {
        return Math.round(getOptions().getTipImageSize() * TIP_IMAGE_MAX_ASPECT);
    }

    /**
     * The extra horizontal advance (px) the tip label is shifted by to sit AFTER the tip image, so the image occupies
     * the branch end (its true position on a time tree) and the name follows. Zero unless the node is an imaged tip
     * whose label starts at the branch end -- so an aligned-phylogram / clustergram label (drawn at a common far
     * column, past all branch ends, where an at-branch-end image doesn't collide) is not shifted.
     */
    private int tipImageAdvance(final PhylogenyNode node) {
        if (!hasTipImage(node) || isAlignedTipLabel(node) || tipLabelsBelowColumns()) {
            return 0;
        }
        // the image's extent along the DEPTH axis (which the label is pushed past): the image is drawn UPRIGHT, so in
        // a vertical orientation the depth axis is the image HEIGHT (the target size), in root-left it is the WIDTH
        // (up to the slot width).
        return (isVerticalOrientation() ? getOptions().getTipImageSize() : tipImageSlotWidth()) + TIP_IMAGE_GAP;
    }

    /**
     * The tip image's footprint (px) along the SPOKE at device angle {@code theta} -- for a radial (circular /
     * unrooted) layout, where the image is drawn UPRIGHT but the branch and label run along the spoke. An upright
     * {@code w x h} rectangle's extent along a direction is {@code w*|cos| + h*|sin|}, so a WIDE photo at an oblique
     * spoke reaches further along the branch than its height would suggest -- using this (not just the height) for the
     * image offset AND the label push keeps the image clear of both the branch and the label. Uses the loaded image's
     * drawn size when available, else a conservative default until it loads.
     */
    private double tipImageRadialFootprint(final PhylogenyNode node, final double theta) {
        final int size = getOptions().getTipImageSize();
        final java.awt.image.BufferedImage img = tipImageCache().get(TipImages.imageRefFor(node), imageBaseDir());
        int dw, dh;
        if (img != null) {
            final int[] wh = TipImages.scaledSize(img.getWidth(), img.getHeight(), size, tipImageSlotWidth());
            dw = wh[0];
            dh = wh[1];
        }
        else {
            dw = tipImageSlotWidth(); // not loaded / broken marker: reserve the widest possible slot
            dh = size;
        }
        return (Math.abs(dw * Math.cos(theta)) + Math.abs(dh * Math.sin(theta)));
    }

    /**
     * Draws each visible external tip's image at the branch end (rectangular root-left family), scaled to the target
     * height (aspect preserved, clamped to the slot width), centred on the tip row and starting just past the node
     * mark -- in the gap the label was shifted right to open ({@link #tipImageAdvance}). A not-yet-loaded image draws
     * nothing this frame (the async load repaints when ready); a broken reference draws nothing. Only in the
     * rectangular non-vertical layout for now.
     */
    private void paintTipImages(final Graphics2D g) {
        if (!getOptions().isShowTipImages()
                || ((getControlPanel() != null) && !isShowExternalDataForThisTab())) {
            return; // tip images are external data -> suppressed with "Show External Data" off, as in the radial layout
        }
        final int target_h = getOptions().getTipImageSize();
        final int slot_w = tipImageSlotWidth();
        final File base = imageBaseDir();
        final Object hint = g.getRenderingHint(RenderingHints.KEY_INTERPOLATION);
        g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
        for (final PhylogenyNode node : _phylogeny.getExternalNodes()) {
            if (isHiddenUnderCollapse(node) || !hasTipImage(node)) {
                continue;
            }
            final String ref = TipImages.imageRefFor(node);
            final java.awt.image.BufferedImage img = tipImageCache().get(ref, base);
            if (img == null) {
                if ((target_h > 0) && tipImageCache().isFailed(ref, base)) {
                    // the reference was recognized but the image could not be loaded (missing file, or a URL that
                    // isn't a direct image -- e.g. a web PAGE link rather than the image file) -> a faint broken-image
                    // marker, so the failure is visible rather than silent (and the reserved space is not just blank)
                    drawBrokenImagePlaceholder(g, node, target_h);
                }
                continue; // else still loading -> nothing this frame (the async load repaints when ready)
            }
            final int[] wh = TipImages.scaledSize(img.getWidth(), img.getHeight(), target_h, slot_w);
            if ((wh[0] <= 0) || (wh[1] <= 0)) {
                continue;
            }
            final int x = Math.round(node.getXcoord()) + effectiveNodeHalfBoxSize(node) + LABEL_GAP_AFTER_NODE_SHAPE;
            final int y = Math.round(node.getYcoord() - (wh[1] / 2.0f));
            g.drawImage(img, x, y, wh[0], wh[1], null);
        }
        restoreInterpolationHint(g, hint);
    }

    /** Restores the interpolation hint we set to BILINEAR for image scaling; when it was unset (null, the default),
     *  reset it to the Java2D default (nearest-neighbour) rather than leaking BILINEAR onto whatever paints next. */
    private static void restoreInterpolationHint(final Graphics2D g, final Object previous) {
        g.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                (previous != null) ? previous : RenderingHints.VALUE_INTERPOLATION_NEAREST_NEIGHBOR);
    }

    /** A faint broken-image marker at a tip whose image failed to load (missing file / a non-image URL such as a web
     *  page link), with its LEFT edge at the branch end (rectangular root-left) -- so a broken reference is VISIBLE. */
    private void drawBrokenImagePlaceholder(final Graphics2D g, final PhylogenyNode node, final int size) {
        final int w = Math.max(8, Math.round(size * 0.85f));
        final int left = Math.round(node.getXcoord()) + effectiveNodeHalfBoxSize(node) + LABEL_GAP_AFTER_NODE_SHAPE;
        drawBrokenImagePlaceholderAt(g, left + (w / 2.0f), node.getYcoord(), size);
    }

    /** A faint outlined box with an X, centred at ({@code cx},{@code cy}) -- the shared broken-image marker (drawn in
     *  device coords, so the vertical/circular callers pass an R-mapped, upright centre). */
    private void drawBrokenImagePlaceholderAt(final Graphics2D g, final float cx, final float cy, final int size) {
        final int w = Math.max(8, Math.round(size * 0.85f));
        final int h = Math.max(8, size);
        final int x = Math.round(cx - (w / 2.0f));
        final int y = Math.round(cy - (h / 2.0f));
        final Color fg = getTreeColorSet().getSequenceColor();
        final Stroke saved = g.getStroke();
        g.setStroke(STROKE_1);
        g.setColor(new Color(fg.getRed(), fg.getGreen(), fg.getBlue(), 95)); // faint
        g.drawRect(x, y, w, h);
        g.drawLine(x, y, x + w, y + h);
        g.drawLine(x + w, y, x, y + h);
        g.setStroke(saved);
        g.setColor(fg);
    }

    /**
     * Vertical-orientation (root-top/bottom) twin of {@link #paintTipImages}: each imaged tip's picture drawn UPRIGHT
     * (never rotated under R -- a frog looks like a frog in every layout) at the slot centre just past the tip along
     * the depth axis, mapped to device coords. Called while g is rotated by R (paintPhylogeny), so it maps the centre
     * through screenPoint and draws in the base frame, like the SYMBOL glyph and TEXT cells.
     */
    private void paintTipImagesVertical(final Graphics2D g) {
        if (!getOptions().isShowTipImages()
                || ((getControlPanel() != null) && !isShowExternalDataForThisTab())) {
            return; // tip images are external data -> suppressed with "Show External Data" off, as in the radial layout
        }
        final int size = getOptions().getTipImageSize();
        final AffineTransform withR = g.getTransform();
        final Object hint = g.getRenderingHint(RenderingHints.KEY_INTERPOLATION);
        g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
        for (final PhylogenyNode node : _phylogeny.getExternalNodes()) {
            if (isHiddenUnderCollapse(node) || !hasTipImage(node)) {
                continue;
            }
            final double lx = node.getXcoord() + effectiveNodeHalfBoxSize(node) + LABEL_GAP_AFTER_NODE_SHAPE
                    + (size / 2.0);
            final Point2D.Double gp = screenPoint(lx, node.getYcoord());
            g.setTransform(_orientation_base_transform);
            drawTipImageUpright(g, node, (float) gp.x, (float) gp.y, size);
            g.setTransform(withR);
        }
        restoreInterpolationHint(g, hint);
    }

    /**
     * Draws {@code node}'s tip image (or its broken-image marker) UPRIGHT, centred at device ({@code cx},{@code cy}),
     * scaled to height {@code size} (aspect preserved, width clamped to the slot). Shared by the vertical and radial
     * (circular / unrooted) paths, which draw the image un-rotated at an R-mapped / spoke-offset centre.
     */
    private void drawTipImageUpright(final Graphics2D g, final PhylogenyNode node, final float cx, final float cy,
                                     final int size) {
        final File base = imageBaseDir();
        final String ref = TipImages.imageRefFor(node);
        final java.awt.image.BufferedImage img = tipImageCache().get(ref, base);
        if (img != null) {
            final int[] wh = TipImages.scaledSize(img.getWidth(), img.getHeight(), size, tipImageSlotWidth());
            if ((wh[0] > 0) && (wh[1] > 0)) {
                g.drawImage(img, Math.round(cx - (wh[0] / 2.0f)), Math.round(cy - (wh[1] / 2.0f)), wh[0], wh[1], null);
            }
        }
        else if (tipImageCache().isFailed(ref, base)) {
            drawBrokenImagePlaceholderAt(g, cx, cy, size);
        }
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
     * inherently colored). Drawn after the node loop (needs the node coords). Dispatched in EVERY layout
     * (rectangular family + circular/unrooted) -- the halo rides the node's device coords, which are set radially too.
     * <p>The screen animation state ({@code _found_halo_bounds} / {@code _has_visible_found_halo}) is cleared at
     * the TOP of {@code paintPhylogeny} and the timer reconciled at its END for EVERY layout -- so an export paint
     * here never clobbers it, and turning the option off reliably STOPS the timer (no leak).
     */
    private void paintFoundNodeHalos(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        final boolean to_screen = !to_pdf && !to_graphics_file;
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
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
        // iterate the phylogeny directly (not _nodes_in_preorder, which is built ONLY in the rectangular branch and
        // is null for a tree opened straight into a radial layout) so the halos work in every layout
        for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
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

    /** Zebra striping for the CIRCULAR layout -- the polar analogue of {@link #paintZebraStripes}: a faint translucent
     *  angular WEDGE (a full pie slice, centre out past the tips + labels + rings) behind every other displayed tip's
     *  angular slice, so a tip is easy to track outward to its annotation ring. Iterates the DISPLAYED leaf-level rows
     *  in tree (= angle) order via iteratorPreorder (NOT _nodes_in_preorder, which is built only in the rectangular
     *  branch), collapsed-clade roots counting as one row; drawn over the tree but under the rings (translucent). */
    private void paintZebraStripesCircular(final Graphics2D g, final int cx, final int cy, final int radius) {
        if (!getOptions().isShowZebraStripes() || (_phylogeny == null) || _export_transparent_background
                || (radius <= 0)) {
            return; // suppressed on a transparent-PNG export, like the rectangular path
        }
        final int displayed = countCircularDisplayedTips(_phylogeny.getRoot());
        if (displayed <= 0) {
            return;
        }
        final double half_step = Math.PI / displayed; // half a tip's angular slice
        final double r_outer = radius + getLongestExtNodeInfo() + circularAnnotationRingsReserve() + CLADE_BAND_RIGHT_PAD;
        final boolean dark = getTreeColorSet().getCurrentColorScheme() == TreeColorSet.DARK_COLOR_SCHEME;
        final Color saved = g.getColor();
        g.setColor(dark ? ZEBRA_STRIPE_ON_DARK : ZEBRA_STRIPE_ON_LIGHT); // faint translucent, theme-aware
        int row = 0;
        for (final java.util.Iterator<PhylogenyNode> it = _phylogeny.iteratorPreorder(); it.hasNext();) {
            final PhylogenyNode node = it.next();
            if ((!node.isExternal() && !node.isCollapse()) || isHiddenUnderCollapse(node)) {
                continue; // only the DRAWN leaf-level rows (tips + collapsed-clade triangles), in angle order
            }
            if ((row % 2) == 1) {
                final Double a = _urt_nodeid_angle_map.get(node.getId());
                if (a != null) {
                    g.fill(annularSector(cx, cy, 0, r_outer, a - half_step, a + half_step)); // full pie wedge
                }
            }
            row++;
        }
        g.setColor(saved);
    }

    /** Whether domain-architecture BOXES are actually drawn in the current layout: always in the rectangular family;
     *  in circular/unrooted ONLY when labels are RADIAL (ride the spoke). Under horizontal ("Radial Labels" off)
     *  labels a spoke-riding domain bar would clash with the upright labels and has no clean radial track, so domains
     *  -- and their legend -- are suppressed there. */
    boolean domainBoxesDrawnInCurrentLayout() {
        return (getControlPanel() != null) && shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES)
                && (!isRadialLayout() || (getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL));
    }

    /** The domain-structure target width to use for the CURRENT layout: the user's {@code _domain_structure_width}
     *  in a rectangular layout, but capped in a radial layout to a fraction of the canvas half-radius (the rectangular
     *  ~0.25*viewport width extends ~half the circular radius and clips). Used consistently by the reservation
     *  ({@code calculateLongestExtNodeInfo}) and the draw factor ({@code initNodeData}) so they agree. */
    private double effectiveDomainStructureWidth() {
        if (!isRadialLayout()) {
            return _domain_structure_width;
        }
        // the radial canvas is the SQUARE radialDiameter (what circularRadius / the preferred size use), NOT getSize()
        // (the panel's actual size, which after a fit is the possibly-large preferred size) -> cap to its half-radius
        final double half = radialDiameter() / 2.0;
        return Math.min(_domain_structure_width, half * RADIAL_DOMAIN_MAX_FRACTION);
    }

    /** The maximum pixel reach a tip label may extend outward in a radial (circular/unrooted) layout: a fraction of
     *  the canvas half-radius, so that tree ring + labels + domains all FIT on fit-to-window (labels longer than this
     *  are truncated with an ellipsis -- full labels in a circle would run ~twice the radius, off the canvas). No cap
     *  in a rectangular layout. */
    private int radialMaxLabelWidth() {
        if (!isRadialLayout()) {
            return Integer.MAX_VALUE;
        }
        return (int) ((radialDiameter() / 2.0) * RADIAL_LABEL_MAX_RATIO);
    }


    /** Protein-domain architectures for the CIRCULAR layout: each external tip's architecture rides its spoke, drawn
     *  as a bar extending radially OUTWARD from a common start radius just past the tip labels -- the iTOL
     *  circular-domains look (a clean CONCENTRIC ring in both aligned and unaligned phylograms). circularRadius already
     *  reserves the concentric worst-case reach, so this only DRAWS. On-box domain-name labels are suppressed in radial
     *  (they would rotate illegibly on a spoke and collide with neighbours); the draggable, E-value-aware domain LEGEND
     *  names them instead -- the same approved exception as the vertical (root-top/bottom) orientation. */
    private void paintDomainsCircular(final Graphics2D g, final int cx, final int cy, final int radius,
                                      final boolean to_pdf, final boolean to_graphics_file) {
        if (!domainBoxesDrawnInCurrentLayout() || (_phylogeny == null) || (radius <= 0)) {
            return;
        }
        final int displayed = countCircularDisplayedTips(_phylogeny.getRoot());
        if (displayed <= 0) {
            return;
        }
        // ALL architectures start at a common radius just past the longest tip label -> a clean CONCENTRIC domain
        // ring (the iTOL look) in BOTH aligned and unaligned phylograms (radius + longest label clears every tip's
        // label, since even a shallow tip's label reaches at most radius + longest_text).
        final double start_r = radius + _length_of_longest_text_only + DOMAIN_RADIAL_GAP;
        // clamp the box thickness to the arc between adjacent spokes at that radius, so dense trees don't overlap
        final int height = TreePanelUtil.domainBoxHeight((float) (start_r * (TWO_PI / displayed)),
                DOMAIN_STRUCTURE_HEIGHT_MIN, DOMAIN_STRUCTURE_HEIGHT_MAX);
        for (final java.util.Iterator<PhylogenyNode> it = _phylogeny.iteratorPreorder(); it.hasNext();) {
            final PhylogenyNode node = it.next();
            if (!node.isExternal() || isHiddenUnderCollapse(node)) {
                continue;
            }
            final Double a = _urt_nodeid_angle_map.get(node.getId());
            if (a != null) {
                paintDomainArchitectureRadial(g, node, cx, cy, a, start_r, height, to_pdf);
            }
        }
    }

    /** Draw {@code node}'s domain architecture as a bar riding the spoke at {@code angle}, extending outward from
     *  {@code (pivot_x, pivot_y)} starting {@code start_dist} px out, {@code height} px thick (perpendicular to the
     *  spoke). Rotating the graphics turns the renderer's horizontal bar into a radial one; on-box labels are
     *  suppressed (radial exception). Shared by the circular pass and the unrooted per-tip draw. */
    private void paintDomainArchitectureRadial(final Graphics2D g, final PhylogenyNode node, final double pivot_x,
                                               final double pivot_y, final double angle, final double start_dist,
                                               final int height, final boolean to_pdf) {
        final RenderableDomainArchitecture rds = renderableDomainArchitectureOf(node);
        if (rds == null) {
            return;
        }
        rds.setRenderingHeight(height);
        final java.awt.geom.AffineTransform saved = g.getTransform();
        g.rotate(angle, pivot_x, pivot_y); // the spoke at `angle` becomes the local +x axis
        rds.render((float) (pivot_x + start_dist), (float) (pivot_y - (height / 2.0)), g, this, to_pdf, false);
        g.setTransform(saved);
    }

    /** The node's domain architecture as a RenderableDomainArchitecture, or null (no sequence / no architecture / a
     *  plain DomainArchitecture from an undo copy not yet re-wrapped). */
    private RenderableDomainArchitecture renderableDomainArchitectureOf(final PhylogenyNode node) {
        if (!node.getNodeData().isHasSequence()) {
            return null;
        }
        final org.forester.phylogeny.data.DomainArchitecture da = node.getNodeData().getSequence()
                .getDomainArchitecture();
        // a plain DomainArchitecture (e.g. from an undo copy not yet re-wrapped) is not renderable -> null
        return (da instanceof RenderableDomainArchitecture) ? (RenderableDomainArchitecture) da : null;
    }

    /** Node age (HPD) bars for the CIRCULAR PHYLOGRAM -- the polar analogue of {@link #paintHpdBars}: a translucent
     *  RADIAL segment along each dated internal node's spoke, spanning its age interval (the radius encodes
     *  distance-from-root = time, so an age range is a radial range). The older bound sits at a SMALLER radius (toward
     *  the root), the younger at a larger one; the segment is anchored to the node's OWN drawn radius plus signed age
     *  deltas (like the rectangular bar's x-anchor), scaled by radius/maxDistanceToRoot -- the same scale the phylogram
     *  draws the node radius with. A thick round-capped stroke gives a constant pixel width at any radius. Circular
     *  phylograms only (the radius must encode distance); skips nodes hidden under a collapse. */
    private void paintHpdBarsCircular(final Graphics2D g, final int cx, final int cy, final int radius,
                                      final boolean to_pdf, final boolean to_graphics_file) {
        if (!getOptions().isShowHpdBars() || !isCircularPhylogram() || (_phylogeny == null) || (radius <= 0)
                || !isBranchLengthTimeCalibrated()) {
            return;
        }
        // while capping, use the CAPPED radial normalizer so the interval scale matches the (capped) r_node -- like the
        // rectangular paintHpdBars, whose corr becomes the capped corr; else the bar is scaled off the uncapped max
        final double max_dist = breakLongBranchesActiveCircular() ? breakCappedRadialMax() : getMaxDistanceToRoot();
        if (max_dist <= 0) {
            return;
        }
        final double radial_corr0 = radius / max_dist; // px per distance/time unit along the spoke (== the phylogram scale)
        // a CALENDAR tree's date increases toward the tips (larger radius); negate so the earlier bound sits at a
        // smaller radius (toward the centre) and the later bound at a larger radius (toward the tips)
        final double radial_corr = (effectiveTimeAxisType() == Options.TIME_AXIS_TYPE.CALENDAR) ? -radial_corr0
                : radial_corr0;
        final Color saved = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setColor(((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) ? HPD_BAR_COLOR_BW
                : HPD_BAR_COLOR);
        g.setStroke(new BasicStroke(HPD_BAR_HEIGHT, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        for (final java.util.Iterator<PhylogenyNode> it = _phylogeny.iteratorPreorder(); it.hasNext();) {
            final PhylogenyNode node = it.next();
            if (node.isExternal() || isHiddenUnderCollapse(node) || !node.getNodeData().isHasDate()) {
                continue;
            }
            final org.forester.phylogeny.data.Date date = node.getNodeData().getDate();
            if ((date.getMin() == null) || (date.getMax() == null)) {
                continue; // need an interval to draw a bar
            }
            final Double ang = _urt_nodeid_angle_map.get(node.getId());
            if (ang == null) {
                continue; // no circular angle (e.g. hidden) -> nothing to place
            }
            final double min = date.getMin().doubleValue();
            final double max = date.getMax().doubleValue();
            final double value = (date.getValue() != null) ? date.getValue().doubleValue() : ((min + max) / 2.0);
            final double r_node = circularRadiusFraction(node) * radius;
            double r_low = r_node - ((max - value) * radial_corr);  // older bound -> smaller radius (toward the root)
            double r_high = r_node + ((value - min) * radial_corr); // younger bound -> larger radius (toward the tips)
            r_low = Math.max(0, Math.min(r_low, r_high)); // robust to swapped/degenerate bounds; never past the centre
            r_high = Math.max(r_low + 1, r_high); // >= 1px floor so a dated node always shows a mark (rectangular parity)
            if (getOptions().getNodeAgeShape() == Options.NODE_AGE_SHAPE.SPINDLE) {
                fillRadialNodeAgeSpindle(g, cx, cy, ang, r_low, r_high, r_node); // a radial lens peaking at the node's radius
            }
            else {
                final double cos = Math.cos(ang), sin = Math.sin(ang);
                drawLine(cx + (r_low * cos), cy + (r_low * sin), cx + (r_high * cos), cy + (r_high * sin), g);
            }
        }
        g.setColor(saved);
        g.setStroke(saved_stroke);
    }

    /** Fills a radial node-age SPINDLE (circular twin of {@link #fillNodeAgeSpindle}): a lens along the spoke at angle
     *  {@code ang}, from the older bound radius {@code r_low} to the younger {@code r_high}, its perpendicular
     *  half-thickness peaking at the node's own radius {@code r_peak}. */
    private void fillRadialNodeAgeSpindle(final Graphics2D g, final int cx, final int cy, final double ang,
                                          final double r_low, final double r_high, final double r_peak) {
        final double cos = Math.cos(ang), sin = Math.sin(ang);
        if ((r_high - r_low) <= SPINDLE_SAMPLE_STEP) { // too narrow to sample a lens: a short radial mark so the node shows
            drawLine(cx + (r_low * cos), cy + (r_low * sin), cx + (r_high * cos), cy + (r_high * sin), g);
            return;
        }
        final double px = -sin, py = cos; // unit vector perpendicular to the spoke
        final java.awt.geom.Path2D.Double path = new java.awt.geom.Path2D.Double();
        path.moveTo(cx + (r_low * cos), cy + (r_low * sin)); // inner tip
        for (double r = r_low + SPINDLE_SAMPLE_STEP; r < r_high; r += SPINDLE_SAMPLE_STEP) {
            final double h = TreePanelUtil.spindleHalfHeightAt(r, r_low, r_high, r_peak, SPINDLE_HALF_HEIGHT);
            path.lineTo(cx + (r * cos) + (h * px), cy + (r * sin) + (h * py));
        }
        path.lineTo(cx + (r_high * cos), cy + (r_high * sin)); // outer tip
        for (double r = r_high - SPINDLE_SAMPLE_STEP; r > r_low; r -= SPINDLE_SAMPLE_STEP) {
            final double h = TreePanelUtil.spindleHalfHeightAt(r, r_low, r_high, r_peak, SPINDLE_HALF_HEIGHT);
            path.lineTo(cx + (r * cos) - (h * px), cy + (r * sin) - (h * py));
        }
        path.closePath();
        g.fill(path);
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
        if (!getOptions().isShowHpdBars() || !getControlPanel().isDrawPhylogram() || (_phylogeny == null)
                || !isBranchLengthTimeCalibrated()) {
            return;
        }
        final double corr = getXcorrectionFactor();
        if (corr <= 0) {
            return; // no branch-length scale -> nothing meaningful to place
        }
        // a CALENDAR tree's date INCREASES toward the tips (opposite of geologic age, which decreases toward the tips);
        // negate corr so hpdBarXRange places the earlier bound to the left and the later bound to the right
        final double signed_corr = (effectiveTimeAxisType() == Options.TIME_AXIS_TYPE.CALENDAR) ? -corr : corr;
        final Color saved = g.getColor();
        g.setColor(((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) ? HPD_BAR_COLOR_BW
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
            final float[] xr = TreePanelUtil.hpdBarXRange(node.getXcoord(), value, min, max, signed_corr);
            final double y = node.getYcoord();
            if (getOptions().getNodeAgeShape() == Options.NODE_AGE_SHAPE.SPINDLE) {
                // a tapered lens from the older bound to the younger, peaking at the node's own (point-estimate) x
                fillNodeAgeSpindle(g, Math.min(xr[0], xr[1]), Math.max(xr[0], xr[1]), node.getXcoord(), y);
            }
            else {
                final int left = Math.round(Math.min(xr[0], xr[1])); // robust to swapped/degenerate bounds
                final int w = Math.max(1, Math.round(Math.abs(xr[1] - xr[0])));
                g.fillRect(left, (int) Math.round(y) - (HPD_BAR_HEIGHT / 2), w, HPD_BAR_HEIGHT);
            }
        }
        g.setColor(saved);
    }

    /** Fills a horizontal node-age SPINDLE (the tapered alternative to the flat HPD bar): a smooth lens from the older
     *  bound {@code x_left} to the younger {@code x_right}, peaking at the point-estimate x {@code x_peak}, centred on
     *  {@code y}. The half-thickness at each x comes from {@link TreePanelUtil#spindleHalfHeightAt}. */
    private void fillNodeAgeSpindle(final Graphics2D g, final double x_left, final double x_right, final double x_peak,
                                    final double y) {
        if ((x_right - x_left) <= SPINDLE_SAMPLE_STEP) { // too narrow to sample a lens: a tiny mark so the node still shows
            g.fillRect((int) Math.round(x_left) - 1, (int) Math.round(y) - (HPD_BAR_HEIGHT / 2), 2, HPD_BAR_HEIGHT);
            return;
        }
        final java.awt.geom.Path2D.Double path = new java.awt.geom.Path2D.Double();
        path.moveTo(x_left, y); // left tip (half-height 0)
        for (double x = x_left + SPINDLE_SAMPLE_STEP; x < x_right; x += SPINDLE_SAMPLE_STEP) {
            path.lineTo(x, y - TreePanelUtil.spindleHalfHeightAt(x, x_left, x_right, x_peak, SPINDLE_HALF_HEIGHT));
        }
        path.lineTo(x_right, y); // right tip
        for (double x = x_right - SPINDLE_SAMPLE_STEP; x > x_left; x -= SPINDLE_SAMPLE_STEP) {
            path.lineTo(x, y + TreePanelUtil.spindleHalfHeightAt(x, x_left, x_right, x_peak, SPINDLE_HALF_HEIGHT));
        }
        path.closePath();
        g.fill(path);
    }

    /**
     * Fossil stratigraphic-range (FAD/LAD) bars: on a dated phylogram, a solid horizontal bar at each TIP that carries
     * a node-age {@code <date>} interval, spanning its observed range -- the First Appearance Datum (oldest, {@code max})
     * to the Last Appearance Datum (youngest, {@code min}) -- with short end-caps so it reads as a bracketed
     * stratigraphic range. This is the TIP analogue of {@link #paintHpdBars} (which draws INTERNAL-node age
     * uncertainty): here it shows how long an extinct taxon is known to have existed, the strap/FAD-LAD convention no
     * interactive tree viewer draws. Reuses {@link TreePanelUtil#hpdBarXRange} for the geometry (anchored to the tip's
     * OWN drawn x plus signed age deltas -- FAD to the left/root, LAD to the right/tips), so it stays put on the tip
     * even when the tree is not strictly ultrametric. Rectangular family (phylograms only); a circular twin is
     * {@link #paintFossilRangeBarsCircular}. Assumes the age unit matches the branch-length (time) unit.
     */
    private void paintFossilRangeBars(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        if (!getOptions().isShowFossilRangeBars() || !getControlPanel().isDrawPhylogram() || (_phylogeny == null)
                || !isBranchLengthTimeCalibrated()) {
            return;
        }
        final double corr = getXcorrectionFactor();
        if (corr <= 0) {
            return; // no branch-length scale -> nothing meaningful to place
        }
        final Color saved = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setColor(((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) ? FOSSIL_BAR_COLOR_BW
                : FOSSIL_BAR_COLOR);
        g.setStroke(STROKE_1); // thin, deterministic end-cap ticks (independent of the ambient branch stroke)
        for (final PhylogenyNode node : _nodes_in_preorder) {
            if (!node.isExternal() || isHiddenUnderCollapse(node) || !node.getNodeData().isHasDate()) {
                continue; // fossil RANGE bars are for TIPS (external nodes)
            }
            final org.forester.phylogeny.data.Date date = node.getNodeData().getDate();
            if ((date.getMin() == null) || (date.getMax() == null)) {
                continue; // need a FAD/LAD range to draw a bar
            }
            final double min = date.getMin().doubleValue();
            final double max = date.getMax().doubleValue();
            final double value = (date.getValue() != null) ? date.getValue().doubleValue() : ((min + max) / 2.0);
            final float[] xr = TreePanelUtil.hpdBarXRange(node.getXcoord(), value, min, max, corr);
            final int left = Math.round(Math.min(xr[0], xr[1])); // robust to swapped/degenerate bounds
            final int right = Math.round(Math.max(xr[0], xr[1]));
            final int y = Math.round(node.getYcoord());
            g.fillRect(left, y - (FOSSIL_BAR_HEIGHT / 2), Math.max(1, right - left), FOSSIL_BAR_HEIGHT);
            // FAD/LAD end-caps: short vertical ticks so the range reads as a bracketed interval
            drawLine(left, y - FOSSIL_BAR_CAP, left, y + FOSSIL_BAR_CAP, g);
            drawLine(right, y - FOSSIL_BAR_CAP, right, y + FOSSIL_BAR_CAP, g);
        }
        g.setColor(saved);
        g.setStroke(saved_stroke);
    }

    /**
     * Circular twin of {@link #paintFossilRangeBars}: on a circular PHYLOGRAM the radius encodes distance-from-root =
     * time, so a fossil tip's stratigraphic range is a RADIAL segment along its spoke (FAD toward the root, LAD toward
     * the outer ring). Mirrors {@link #paintHpdBarsCircular} but for external tips, using the fossil colour; no
     * end-caps in the fan (a perpendicular tick would foul neighbouring spokes). Unrooted is N/A (no distance ring).
     */
    private void paintFossilRangeBarsCircular(final Graphics2D g, final int cx, final int cy, final int radius,
                                              final boolean to_pdf, final boolean to_graphics_file) {
        if (!getOptions().isShowFossilRangeBars() || !isCircularPhylogram() || (_phylogeny == null) || (radius <= 0)
                || !isBranchLengthTimeCalibrated()) {
            return;
        }
        // while capping, use the CAPPED radial normalizer so the range-bar scale matches the (capped) r_node (parity
        // with the rectangular fossil bars, which ride the capped corr)
        final double max_dist = breakLongBranchesActiveCircular() ? breakCappedRadialMax() : getMaxDistanceToRoot();
        if (max_dist <= 0) {
            return;
        }
        final double radial_corr = radius / max_dist; // px per distance/time unit along the spoke (== the phylogram scale)
        final Color saved = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setColor(((to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite()) ? FOSSIL_BAR_COLOR_BW
                : FOSSIL_BAR_COLOR);
        g.setStroke(new BasicStroke(FOSSIL_BAR_HEIGHT, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
        for (final java.util.Iterator<PhylogenyNode> it = _phylogeny.iteratorPreorder(); it.hasNext();) {
            final PhylogenyNode node = it.next();
            if (!node.isExternal() || isHiddenUnderCollapse(node) || !node.getNodeData().isHasDate()) {
                continue; // fossil RANGE bars are for TIPS (external nodes)
            }
            final org.forester.phylogeny.data.Date date = node.getNodeData().getDate();
            if ((date.getMin() == null) || (date.getMax() == null)) {
                continue; // need a FAD/LAD range to draw a bar
            }
            final Double ang = _urt_nodeid_angle_map.get(node.getId());
            if (ang == null) {
                continue; // no circular angle (e.g. hidden) -> nothing to place
            }
            final double min = date.getMin().doubleValue();
            final double max = date.getMax().doubleValue();
            final double value = (date.getValue() != null) ? date.getValue().doubleValue() : ((min + max) / 2.0);
            final double r_node = circularRadiusFraction(node) * radius;
            double r_low = r_node - ((max - value) * radial_corr);  // FAD (older) -> smaller radius (toward the root)
            double r_high = r_node + ((value - min) * radial_corr); // LAD (younger) -> larger radius (toward the tips)
            r_low = Math.max(0, Math.min(r_low, r_high)); // robust to swapped/degenerate bounds; never past the centre
            r_high = Math.max(r_low + 1, r_high); // >= 1px floor so a dated tip always shows a mark
            final double cos = Math.cos(ang), sin = Math.sin(ang);
            drawLine(cx + (r_low * cos), cy + (r_low * sin), cx + (r_high * cos), cy + (r_high * sin), g);
        }
        g.setColor(saved);
        g.setStroke(saved_stroke);
    }

    /** Diameter (px, tree/device space) of an ancestral-state pie: ~2.5x the node dot, floored at 10px so the
     *  wedges stay legible even at a tiny node-shape size. */
    private double ancestralPieDiameter() {
        return Math.max(10.0, 2.5 * getOptions().getDefaultNodeShapeSize());
    }

    /**
     * Draws an ancestral-state pie at each drawn node: wedges sized by the node's posterior state probabilities and
     * colored from the stable state->color map ({@link #_ancestral_pie_colors}). Tips with a single observed state
     * render as one full-circle wedge (a solid state-colored disc). Drawn AFTER the node-paint loop (which sets the
     * coords), alongside the other tree-riding overlays; the pie replaces the plain node dot at pie nodes (see the
     * suppression guard at the drawPropertyColorDot call site). Rectangular family; in the root-top/bottom
     * orientations it rides the rotated frame (the disc stays a disc; the wedge orientation rotates with the tree).
     */
    private void paintAncestralPies(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        if (!isShowAncestralPies() || (_ancestral_pie_dist == null)) {
            return;
        }
        final boolean bw = (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite();
        final int n = _ancestral_pie_colors.size();
        final double d = ancestralPieDiameter();
        final double r = d / 2.0;
        final boolean radial = isRadialLayout();
        final Color saved = g.getColor();
        final Stroke saved_stroke = g.getStroke();
        g.setStroke(STROKE_05);
        final Color outline = bw ? Color.BLACK : getTreeColorSet().getBranchColor();
        // Walk the tree in DETERMINISTIC preorder (not the per-node cache's IdentityHashMap address-order, which
        // varies run to run) and look up each node's cached distribution -- so overlapping-pie z-order is reproducible
        // in exports. Works in every layout: each node's device coords are set by whichever layout ran (the cache is
        // used only to avoid re-parsing, NOT for iteration; _nodes_in_preorder is built only in the rectangular branch).
        for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext();) {
            final PhylogenyNode node = it.next();
            final List<TreePanelUtil.StateProbability> dist = _ancestral_pie_dist.get(node);
            if (dist == null) {
                continue; // no distribution for this node
            }
            // skip nodes whose device coords are not validly set: any node hidden under a collapse, plus a collapsed
            // clade-ROOT in a radial layout -- paintCirculars/paintUnrooted don't position collapsed nodes (unlike the
            // rectangular layout, which does), so it would otherwise draw a phantom pie at stale coords
            if (isHiddenUnderCollapse(node) || (radial && node.isCollapse())) {
                continue;
            }
            drawAncestralPie(g, node.getXcoord(), node.getYcoord(), r, d, dist, n, outline, bw);
        }
        g.setStroke(saved_stroke);
        g.setColor(saved);
    }

    /** Fills the wedges of one pie in GLOBAL state order (iterating {@code _ancestral_pie_colors}, so a state occupies
     *  a comparable arc + color across every node), then strokes a thin outline circle. In B&W export the wedges map
     *  to the same evenly-spaced gray ramp the legend uses (by global index). Allocation-free: the small per-node
     *  {@code dist} list is scanned directly (no per-node map). */
    private void drawAncestralPie(final Graphics2D g, final double cx, final double cy, final double r,
                                  final double d, final List<TreePanelUtil.StateProbability> dist, final int n,
                                  final Color outline, final boolean bw) {
        double start = 90.0; // begin at 12 o'clock and sweep clockwise (negative arc angle)
        int i = 0;
        for (final Map.Entry<String, Color> e : _ancestral_pie_colors.entrySet()) {
            double prob = 0.0;
            for (final TreePanelUtil.StateProbability sp : dist) { // dist is tiny (2-5 states); sum defends a dup state
                if (e.getKey().equals(sp.getState())) {
                    prob += sp.getProbability();
                }
            }
            if (prob > 0.0) {
                final double sweep = 360.0 * prob;
                g.setColor(bw ? TreePanelUtil.grayShade(i, n) : e.getValue());
                _arc.setArc(cx - r, cy - r, d, d, start, -sweep, Arc2D.PIE);
                g.fill(_arc);
                start -= sweep;
            }
            ++i;
        }
        g.setColor(outline);
        _ellipse.setFrame(cx - r, cy - r, d, d);
        g.draw(_ellipse);
    }

    /** Draws the Color-by / Size-by tip dot at each node in a RADIAL (circular/unrooted) layout -- the rectangular
     *  path draws it inside paintNodeData, which the radial paths don't use. Mirrors the rectangular gate exactly:
     *  Color-by covers external + collapsed nodes, Size-by covers external nodes; both require "Show External Data";
     *  and a pie node is skipped (the pie IS the marker there). Dispatched after the radial recursion has set coords. */
    private void paintRadialPropertyDots(final Graphics2D g) {
        if ((!isColorByProperty() && !isSizeByProperty()) || !isShowExternalDataForThisTab()) {
            return;
        }
        for (final PhylogenyNodeIterator it = _phylogeny.iteratorPreorder(); it.hasNext();) {
            final PhylogenyNode node = it.next();
            // skip nodes hidden under a collapse AND collapsed clade-roots themselves: radial layouts don't position
            // collapsed nodes, so a dot there would be a phantom at stale coords
            if (isHiddenUnderCollapse(node) || node.isCollapse()) {
                continue;
            }
            final boolean want_dot = (isColorByProperty() && node.isExternal())
                    || (isSizeByProperty() && node.isExternal());
            if (want_dot && !(isShowAncestralPies() && (_ancestral_pie_dist != null)
                    && _ancestral_pie_dist.containsKey(node))) {
                drawPropertyColorDot(g, node);
            }
        }
    }

    /** The post-layout overlays for the radial (CIRCULAR/UNROOTED) branches, dispatched identically by both after the
     *  tree geometry has set every node's device coords: Color-by/Size-by tip dots (the node loop draws them for the
     *  rectangular family), ancestral-state pies, the translucent hover preview, and the pulsing found-node halos.
     *  Kept as one method so a future radial overlay is added once, not in two places that can drift apart. */
    private void paintRadialOverlays(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file) {
        paintRadialPropertyDots(g); // Color-by / Size-by tip dots (rectangular draws them in the node loop)
        paintAncestralPies(g, to_pdf, to_graphics_file); // per-node state pies
        paintHoverPreview(g, !(to_pdf || to_graphics_file)); // translucent select/deselect hover preview
        paintFoundNodeHalos(g, to_pdf, to_graphics_file); // pulsing (screen) / static-glow (export) hit halos
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

    /** Radial space (px) reserved OUTSIDE the tip labels for the circular clade BARS/BRACKETS (the arc itself; a long
     *  taxon label riding the spoke past it may still overflow -- a known limitation, like the rectangular label reach).
     *  0 for the BOXES mode (annular sectors sit within the labels) or when no clade bands are shown. */
    private int circularCladeBandReserve() {
        if ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) || !hasCladeBands()
                || (_clade_bands_mode == CLADE_VIS.BOXES)) {
            return 0;
        }
        // Every level's ring must be inside the reserve, and so must every label EXCEPT the outermost one's --
        // that last one riding past the edge is the pre-existing, documented limitation above, and reserving it
        // now would shrink every single-level circular tree for no reason.
        final java.util.List<CladeLevel> levels = drawableCladeLevels();
        double r = 0;
        for (int i = 0; i < levels.size(); ++i) {
            r += circularCladeLevelRadialExtent(levels.get(i));
            if (i == (levels.size() - 1)) {
                r -= maxCladeLabelWidth(levels.get(i)); // the outermost label may still overflow, as before
            }
        }
        return (int) Math.ceil(r);
    }

    /** The angular extent {@code {a0, a1}} (radians, a0<=a1 in tree order = increasing angle) of a clade's DISPLAYED
     *  tips in the circular layout, or null if none. Walks the DISPLAYED tips -- a visible external, OR a collapsed
     *  clade-root (one angular slot; the walk stops there) -- exactly like assignCircularDisplayedTipAngles, so a band
     *  containing a collapsed sub-clade still covers that sub-clade's slot. */
    private double[] circularCladeAngleRange(final PhylogenyNode clade_root) {
        final Double[] first_last = new Double[2];
        collectCircularDisplayedTipAngles(clade_root, first_last);
        return (first_last[0] == null) ? null : new double[] { first_last[0], first_last[1] };
    }

    /** Records the FIRST and LAST displayed-tip angle under {@code node} (tree order) into {@code first_last}. Stops at
     *  a visible external or a collapsed clade-root (its hidden descendants are not displayed). */
    private void collectCircularDisplayedTipAngles(final PhylogenyNode node, final Double[] first_last) {
        if (node.isExternal() || node.isCollapse()) {
            final Double a = _urt_nodeid_angle_map.get(node.getId());
            if (a != null) {
                if (first_last[0] == null) {
                    first_last[0] = a;
                }
                first_last[1] = a;
            }
            return;
        }
        for (int i = 0; i < node.getNumberOfDescendants(); ++i) {
            collectCircularDisplayedTipAngles(node.getChildNode(i), first_last);
        }
    }

    /** Clade bands for the CIRCULAR layout -- the polar analogue of paintCladeBands. BOXES become translucent ANNULAR
     *  SECTORS spanning each clade's angular extent (clade-root radius out past the tips + labels); BARS a solid colour
     *  arc at that outer radius with a radial taxon label; BRACKETS a thin arc + radial end-ticks + label. cx/cy = the
     *  ring centre, radius = the tip ring radius. */
    private void paintCladeBandsCircular(final Graphics2D g, final int cx, final int cy, final int radius) {
        if (!hasCladeBands() || (radius <= 0)) {
            return;
        }
        // just past the tips + labels AND any annotation rings (so the clade marks sit outside the ring stack)
        final double base = radius + getLongestExtNodeInfo() + circularAnnotationRingsReserve() + CLADE_BAND_RIGHT_PAD;
        final int displayed = countCircularDisplayedTips(_phylogeny.getRoot());
        final double half_step = (displayed > 0) ? (Math.PI / displayed) : 0; // half a tip's angular slice
        double r_outer = base; // the radius cursor: each level's ring sits outside the one before it
        for (final CladeLevel level : drawableCladeLevels()) {
            for (final CladeBand band : level.getBands()) {
                if (skipCladeBand(band)) { // single-member clades: no bar/bracket (boxes still draw)
                    continue;
                }
                final double[] ar = circularCladeAngleRange(band.getRoot());
                if (ar == null) {
                    continue;
                }
                final double a0 = ar[0] - half_step, a1 = ar[1] + half_step; // cover the outer tips fully
                switch (_clade_bands_mode) {
                    case BARS:
                        drawCircularCladeArc(g, cx, cy, r_outer, a0, a1, band);
                        break;
                    case BRACKETS:
                        drawCircularCladeBracket(g, cx, cy, r_outer, a0, a1, band);
                        break;
                    default:
                        drawCircularCladeSector(g, cx, cy, band, a0, a1, r_outer);
                        break;
                }
            }
            r_outer += circularCladeLevelRadialExtent(level);
        }
    }

    /** The radial thickness one circular bar/bracket level occupies: the gap, the mark, and the taxon label riding
     *  the spoke outward from it (so the NEXT level's ring clears this level's names). */
    private double circularCladeLevelRadialExtent(final CladeLevel level) {
        final int mark = (_clade_bands_mode == CLADE_VIS.BARS) ? (CLADE_BAR_WIDTH + 3) : (CLADE_BRACKET_TICK + 4);
        return CLADE_BAR_GAP + mark + maxCladeLabelWidth(level) + 4;
    }

    /** A translucent annular SECTOR (the circular "clade box"): from the clade root's radius out to r_outer over the
     *  clade's angular span, built as the outer PIE minus the inner PIE so the tree shows through. */
    private void drawCircularCladeSector(final Graphics2D g, final int cx, final int cy, final CladeBand band,
                                         final double a0, final double a1, final double r_outer) {
        final PhylogenyNode root = band.getRoot();
        final double r_inner = Math.hypot(root.getXcoord() - cx, root.getYcoord() - cy);
        if ((r_outer - r_inner) < 1) {
            return;
        }
        final Area sector = annularSector(cx, cy, r_inner, r_outer, a0, a1);
        final Color c = band.getColor();
        final Color saved = g.getColor();
        g.setColor(new Color(c.getRed(), c.getGreen(), c.getBlue(), CLADE_BOX_ALPHA));
        g.fill(sector);
        g.setColor(saved);
    }

    /** A solid colour ARC (the circular "clade bar") just past r_outer, a CLADE_BAR_WIDTH-thick annular sector over the
     *  clade's angular span, plus the taxon label riding the spoke at the mid-angle beyond it. */
    private void drawCircularCladeArc(final Graphics2D g, final int cx, final int cy, final double r_outer,
                                      final double a0, final double a1, final CladeBand band) {
        final double r0 = r_outer + CLADE_BAR_GAP;
        final Area arc = annularSector(cx, cy, r0, r0 + CLADE_BAR_WIDTH, a0, a1);
        final Color saved = g.getColor();
        g.setColor(band.getColor());
        g.fill(arc);
        g.setColor(getTreeColorSet().getSequenceColor());
        drawRadialTextCentered(g, band.getTaxon(), cx, cy, r0 + CLADE_BAR_WIDTH + 3, (a0 + a1) / 2.0);
        g.setColor(saved);
    }

    /** A thin arc BRACKET (the circular monochrome "]" clade mark) just past r_outer with short radial end-ticks + label. */
    private void drawCircularCladeBracket(final Graphics2D g, final int cx, final int cy, final double r_outer,
                                          final double a0, final double a1, final CladeBand band) {
        final double r = r_outer + CLADE_BAR_GAP;
        final Stroke saved_stroke = g.getStroke();
        final Color saved = g.getColor();
        g.setStroke(new BasicStroke(CLADE_BRACKET_STROKE));
        g.setColor(getTreeColorSet().getSequenceColor());
        drawArc(cx - r, cy - r, 2 * r, 2 * r, -a0, -(a1 - a0), g); // the spine arc (drawArc takes radians)
        for (final double a : new double[] { a0, a1 }) { // short radial end-ticks pointing inward toward the tips
            final double cos = Math.cos(a), sin = Math.sin(a);
            drawLine(cx + (r * cos), cy + (r * sin), cx + ((r - CLADE_BRACKET_TICK) * cos),
                    cy + ((r - CLADE_BRACKET_TICK) * sin), g);
        }
        drawRadialTextCentered(g, band.getTaxon(), cx, cy, r + 4, (a0 + a1) / 2.0);
        g.setStroke(saved_stroke);
        g.setColor(saved);
    }

    /** An annular-sector {@link Area} (outer PIE minus inner PIE) spanning device angles a0..a1 (radians, y-down).
     *  Arc2D angles are CCW in a y-up frame; device space is y-down, so a device angle theta maps to Arc2D angle -theta. */
    private Area annularSector(final int cx, final int cy, final double r_inner, final double r_outer, final double a0,
                               final double a1) {
        final double start_deg = -Math.toDegrees(a0);
        final double extent_deg = -Math.toDegrees(a1 - a0);
        final Area outer = new Area(new Arc2D.Double(cx - r_outer, cy - r_outer, 2 * r_outer, 2 * r_outer, start_deg,
                extent_deg, Arc2D.PIE));
        outer.subtract(new Area(new Arc2D.Double(cx - r_inner, cy - r_inner, 2 * r_inner, 2 * r_inner, start_deg,
                extent_deg, Arc2D.PIE)));
        return outer;
    }

    /** Draws {@code text} centred on the point at radius {@code r} and device angle {@code angle} from (cx,cy), rotated
     *  to ride the spoke and kept upright on the far half of the fan (the radial-label convention). */
    private void drawRadialTextCentered(final Graphics2D g, final String text, final int cx, final int cy,
                                        final double r, final double angle) {
        if (ForesterUtil.isEmpty(text)) {
            return;
        }
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = getTreeFontSet().getFontMetricsSmall();
        double m = angle % TWO_PI;
        if (m < 0) {
            m += TWO_PI;
        }
        if ((m > HALF_PI) && (m < ONEHALF_PI)) {
            m -= PI; // upright on the far half of the fan
        }
        final double ar = r + (fm.stringWidth(text) / 2.0); // centre it a half-width out so its NEAR edge sits at r
        final double lx = cx + (ar * Math.cos(angle)), ly = cy + (ar * Math.sin(angle));
        final AffineTransform saved = g.getTransform();
        g.rotate(m, lx, ly);
        TreePanel.drawString(text, (float) (lx - (fm.stringWidth(text) / 2.0)),
                (float) (ly + (fm.getAscent() / 2.0) - (fm.getDescent() / 2.0)), g);
        g.setTransform(saved);
    }

    /** Radial space (px) reserved OUTSIDE the tip labels for the concentric annotation-column RINGS in the circular
     *  layout: the leading gap plus each ring's width and inter-ring gap (the polar analogue of the horizontal column
     *  stack the rectangular layout reserves). 0 when not circular or no annotation columns are shown. */
    private int circularAnnotationRingsReserve() {
        if ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) || !hasAnnotationColumns()) {
            return 0;
        }
        return sumAnnotationColumnWidths();
    }

    /**
     * The circular tip-ring radius for a canvas of {@code pref_w} x {@code pref_h}: the available half-canvas minus
     * the fixed margins (MOVE + clade-band + annotation-ring reserves) minus the tip-label reach -- but the label
     * reach is CAPPED at {@link #RADIAL_LABEL_MAX_RATIO} of the available radius, so a tree whose labels are long
     * relative to the canvas still draws a real circle (the labels then extend past the canvas edge) instead of the
     * radius going negative and collapsing every tip onto the centre. Never negative. Shared by the paint + the
     * ring-click hit-test so they agree.
     */
    private int circularRadius(final double side) {
        final double avail = Math.max(0, (side / 2.0) - MOVE - circularCladeBandReserve()
                - circularAnnotationRingsReserve());
        if (domainBoxesDrawnInCurrentLayout()) {
            // Domains draw as a CONCENTRIC ring: every architecture starts at radius + _length_of_longest_text_only
            // (+ gap + the renderer's 20px lead-in) and extends by up to _longest_rendered_domain. The worst-case
            // reach is the SUM of those two SEPARATE maxima -- the longest-label tip and the widest-domain tip may be
            // DIFFERENT tips -- NOT getLongestExtNodeInfo() (a per-tip max of text+domain, which under-reserves and
            // clips when they differ). Reserve the sum so both fit on fit-to-window; floor the ring so it can't
            // collapse (this SHRINKS the ring rather than clipping -- the user wants domains + labels to fit).
            final double reach = _length_of_longest_text_only + _longest_rendered_domain + DOMAIN_RADIAL_GAP + 20;
            return (int) Math.max(avail - reach, avail * RADIAL_MIN_TREE_RATIO);
        }
        final int label = (int) Math.min(getLongestExtNodeInfo(), avail * RADIAL_LABEL_MAX_RATIO);
        return (int) (avail - label);
    }

    /** The radius (px from the ring centre) where the annotation-column rings START: just past the tips + labels. */
    private double circularAnnotationRingStart(final int radius) {
        return radius + getLongestExtNodeInfo() + ANN_COL_GAP;
    }

    /** Annotation columns for the CIRCULAR layout -- the polar analogue of {@link #paintAnnotationColumns}: each column
     *  becomes a concentric RING outside the tips + labels, each tip's cell an arc segment over that tip's angular slice
     *  (colour strip / heat-map / matrix = a filled sector; bar = a sector grown radially by the value; text = a radial
     *  label). Each ring carries its field-name header as a tangential label at the fan seam. cx/cy = the ring centre,
     *  radius = the tip ring radius. */
    private void paintAnnotationColumnsCircular(final Graphics2D g, final int cx, final int cy, final int radius) {
        if (!hasAnnotationColumns() || (radius <= 0)) {
            return;
        }
        final java.util.List<PhylogenyNode> tips = visibleExternalTips();
        final int displayed = countCircularDisplayedTips(_phylogeny.getRoot());
        final double half_step = (displayed > 0) ? (Math.PI / displayed) : 0; // half a tip's angular slice
        final Color fg = getTreeColorSet().getSequenceColor();
        final Color saved_color = g.getColor();
        final Font saved_font = g.getFont();
        double r = circularAnnotationRingStart(radius);
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            final int w = annotationColumnWidth(i);
            final double r0 = r, r1 = r + w;
            final AnnotationColumns.Type type = _annotation_columns.getColumn(i).getType();
            for (final PhylogenyNode t : tips) {
                final Double a = _urt_nodeid_angle_map.get(t.getId());
                if (a == null) {
                    continue; // a tip with no circular angle (hidden under a collapse) has no cell
                }
                final double a0 = a - half_step, a1 = a + half_step;
                switch (type) {
                    case COLOR_STRIP:
                    case HEATMAP:
                    case MATRIX: {
                        final Color c = _annotation_columns.cellColor(t, i);
                        if (c != null) {
                            g.setColor(c);
                            g.fill(annularSector(cx, cy, r0, r1, a0, a1));
                        }
                        break;
                    }
                    case SYMBOL: {
                        // a shape glyph centred in this tip's ring cell (mid-radius, on its spoke), drawn UPRIGHT
                        // (no spoke rotation) so every shape -- incl. a triangle -- reads the same as in the other layouts
                        final AnnotationColumns.Fill fill = _annotation_columns.symbolFill(t, i);
                        if (fill != AnnotationColumns.Fill.NONE) {
                            final Color c = _annotation_columns.cellColor(t, i);
                            if (c != null) {
                                final double rmid = (r0 + r1) / 2.0;
                                final float px = (float) (cx + (rmid * Math.cos(a)));
                                final float py = (float) (cy + (rmid * Math.sin(a)));
                                final double arc = rmid * (a1 - a0); // ~ the tip's arc width at this radius
                                drawSymbolGlyph(g, c, px, py, (float) (Math.min(w, arc) - 2),
                                        fill == AnnotationColumns.Fill.FILLED, _annotation_columns.symbolShape(i));
                            }
                        }
                        break;
                    }
                    case BAR: {
                        // a present value always fills at least a 1px-thick stub (so the minimum value is visible and
                        // distinct from a missing value, which is NaN and fills nothing); the bar grows RADIALLY outward
                        final double f = _annotation_columns.barFraction(t, i);
                        if (!Double.isNaN(f)) {
                            g.setColor(fg);
                            g.fill(annularSector(cx, cy, r0, Math.max(r0 + 1, r0 + (f * w)), a0, a1));
                        }
                        break;
                    }
                    case STACKED_BAR: {
                        // the segments stack RADIALLY outward within this ring's [r0, r1] band, each a coloured annulus
                        final double[] fr = _annotation_columns.stackFractions(t, i);
                        final java.util.List<Color> cols = _annotation_columns.stackColors(i);
                        double off = 0;
                        for (int k = 0; k < fr.length; ++k) {
                            final double f = fr[k];
                            if (f > 0) {
                                // f > 0 -> outer radius > inner, so the sector is non-empty; a sub-pixel slice draws
                                // as ~nothing (no forced 1px that would overshoot the ring's outer edge)
                                final double rs0 = r0 + (off * w);
                                g.setColor(cols.get(k));
                                g.fill(annularSector(cx, cy, rs0, r0 + ((off + f) * w), a0, a1));
                            }
                            off += f;
                        }
                        break;
                    }
                    case PIE: {
                        // a pie centred in this tip's ring cell (mid-radius, on its spoke), drawn UPRIGHT so the wedge
                        // order reads the same as in the other layouts
                        final double rmid = (r0 + r1) / 2.0;
                        final float px = (float) (cx + (rmid * Math.cos(a)));
                        final float py = (float) (cy + (rmid * Math.sin(a)));
                        final double arc = rmid * (a1 - a0);
                        drawPieGlyph(g, _annotation_columns.stackFractions(t, i), _annotation_columns.stackColors(i),
                                px, py, (float) (Math.min(w, arc) - 2), fg);
                        break;
                    }
                    case TEXT: {
                        final String v = _annotation_columns.cellText(t, i);
                        if (v.length() > 0) {
                            g.setColor(fg);
                            drawRadialTextCentered(g, v, cx, cy, r0, a); // rides the spoke outward from the ring's inner edge
                        }
                        break;
                    }
                }
            }
            // the field-name header: a horizontal label centred over the tree at the TOP (12 o'clock), at the ring's
            // own radius so the stack reads top-to-bottom by ring -- legible and within its band (unlike a tangential
            // label at the seam, whose chord crosses neighbouring thin rings)
            final String header = _annotation_columns.getColumn(i).getHeader();
            if (header.length() > 0) {
                g.setColor(fg);
                drawTangentialText(g, header, cx, cy, (r0 + r1) / 2.0, -HALF_PI);
            }
            r = r1 + annotationColumnGapAfter(i);
        }
        g.setColor(saved_color);
        g.setFont(saved_font);
    }

    /** The annotation column whose RING band contains the click (circular layout), or -1. Lets a click anywhere on a
     *  ring focus that column's colour key -- the polar analogue of clicking a rectangular column header. */
    private int circularAnnotationRingAt(final int x, final int y) {
        if (!hasAnnotationColumns() || (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
            return -1;
        }
        final int radius = _circular_radius; // the exact radius the last paint used, so clicks hit the drawn rings
        if (radius <= 0) {
            return -1;
        }
        final double rr = Math.hypot(x - _circular_center_x, y - _circular_center_y);
        double r = circularAnnotationRingStart(radius);
        for (int i = 0; i < _annotation_columns.size(); ++i) {
            final int w = annotationColumnWidth(i);
            if ((rr >= r) && (rr <= (r + w))) {
                return i;
            }
            r += w + annotationColumnGapAfter(i);
        }
        return -1;
    }

    /** Draws {@code text} centred on the point at radius {@code r} / device angle {@code angle}, rotated TANGENTIALLY
     *  (reading along the ring, perpendicular to the spoke) and kept upright. Used for the circular ring headers. */
    private void drawTangentialText(final Graphics2D g, final String text, final int cx, final int cy, final double r,
                                    final double angle) {
        if (ForesterUtil.isEmpty(text)) {
            return;
        }
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics fm = getTreeFontSet().getFontMetricsSmall();
        double m = (angle + HALF_PI) % TWO_PI; // the tangent direction
        if (m < 0) {
            m += TWO_PI;
        }
        if ((m > HALF_PI) && (m < ONEHALF_PI)) {
            m -= PI; // keep upright
        }
        final double px = cx + (r * Math.cos(angle)), py = cy + (r * Math.sin(angle));
        final AffineTransform saved = g.getTransform();
        g.rotate(m, px, py);
        TreePanel.drawString(text, (float) (px - (fm.stringWidth(text) / 2.0)),
                (float) (py + (fm.getAscent() / 2.0) - (fm.getDescent() / 2.0)), g);
        g.setTransform(saved);
    }

    /**
     * The x just past the end of the longest currently shown tip label (all text fields plus any aligned
     * domain track, via {@link #getLongestExtNodeInfo()}). Mirrors the tree's own label-end
     * logic: in a phylogram labels start at each tip's own x, so the rightmost label ends past the DEEPEST
     * tip; in an aligned/cladogram view all tips share one x. When external labels are hidden ("Show External
     * Data" off) they occupy no width, so annotations sit right after the tips.
     */
    private float labelsRightEdge() {
        // in a vertical orientation the tip labels are TILTED and extend along the DEPTH axis only by their depth
        // component (full width for 90deg, L/sqrt(2) for 45deg), so the columns clear the labels using that reach --
        // not the full horizontal label width, which would leave a large empty gap below the tips
        final float label_w = isShowExternalDataForThisTab()
                ? (isVerticalOrientation() ? depthLabelReserve() : getLongestExtNodeInfo()) : 0;
        return tipsDepthEdge() + label_w;
    }

    /** The depth (logical x) where the tip branches END -- the anchor for tip labels and tip-aligned columns, WITHOUT
     *  the tip-label reach. (labelsRightEdge() adds the label reach; the clustergram "labels below columns" layout puts
     *  the columns here, right past the tips, and the labels past the columns instead.) */
    private float tipsDepthEdge() {
        if (getControlPanel().isDrawPhylogram()) {
            return (float) ((displayedMaxDistanceToRoot() * getXcorrectionFactor()) + _phylogeny.getRoot().getXcoord());
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

    /** The x where clade bands/bars start: past the tip labels, any tip-aligned annotation columns, and the MSA track. */
    private float cladeBandRightEdge() {
        return labelsRightEdge() + annotationColumnsWidth() + msaTrackWidth() + CLADE_BAND_RIGHT_PAD;
    }

    // ------------------------------------------------------------------------------------------------------------------
    // Sequence-alignment (MSA) track -- colored residue cells beside the tips, rectangular root-left only.
    // ------------------------------------------------------------------------------------------------------------------

    /** Whether the alignment track is drawn now: option on, a root-left rectangular layout, and the tree carries an
     *  aligned sequence (so vertical/circular reserve nothing and draw nothing -- the approved per-item exception). */
    private boolean msaShown() {
        return getOptions().isShowMsa() && !isVerticalOrientation() && !isRadialLayout() && (alignmentLength() > 0);
    }

    /** The number of alignment columns = the longest aligned tip sequence (uncached; O(tips), called a few times per
     *  layout/paint). Tips with a shorter/absent sequence simply run out of cells. */
    private int alignmentLength() {
        int max = 0;
        if ((_phylogeny == null) || _phylogeny.isEmpty()) {
            return 0;
        }
        for (final PhylogenyNode ext : _phylogeny.getExternalNodes()) {
            if (ext.getNodeData().isHasSequence()) {
                final Sequence s = ext.getNodeData().getSequence();
                if (s.isMolecularSequenceAligned()) {
                    final String mol = s.getMolecularSequence();
                    if ((mol != null) && (mol.length() > max)) {
                        max = mol.length();
                    }
                }
            }
        }
        return max;
    }

    /** Whether the point is inside the OVERVIEW box, which floats over the viewport and -- in both RIGHT
     *  placements -- sits exactly where the alignment window is drawn. Something painted on top of the alignment
     *  hides it, so a readout for the residue underneath would describe a cell the user cannot even see. */
    private boolean isOverFloatingOverlay(final int x, final int y) {
        if (!isOvOn() || (getOvMaxWidth() <= 0) || (getOvMaxHeight() <= 0)) {
            return false;
        }
        final Rectangle vis = getVisibleRect();
        // the overview is positioned in VIEWPORT coordinates, so compare against the scrolled viewport origin
        final int ox = vis.x + getOvXPosition();
        final int oy = vis.y + getOvYPosition();
        return (x >= ox) && (x < (ox + getOvMaxWidth())) && (y >= oy) && (y < (oy + getOvMaxHeight()));
    }

    /** For tests: {x, y, w, h} of the overview box in panel coordinates, or null when it is not drawn. */
    int[] floatingOverlayRectForTest() {
        if (!isOvOn() || (getOvMaxWidth() <= 0) || (getOvMaxHeight() <= 0)) {
            return null;
        }
        final Rectangle vis = getVisibleRect();
        return new int[] { vis.x + getOvXPosition(), vis.y + getOvYPosition(), Math.round(getOvMaxWidth()),
                           Math.round(getOvMaxHeight()) };
    }

    /** One cell of the alignment display: which tip's row, and which alignment column (0-based). */
    static final class MsaCell {

        final PhylogenyNode _tip;
        final int           _column;

        MsaCell(final PhylogenyNode tip, final int column) {
            _tip = tip;
            _column = column;
        }
    }

    /**
     * The alignment cell under a device point, or null when the point is not over a drawn residue.
     * <p>
     * Deliberately mirrors {@code paintMsaTrack}'s geometry term for term -- same origin
     * ({@code annotationColumnsEndX() + MSA_TRACK_GAP}), same column width, same window offset, same per-tip row
     * band -- because a hit-test that drifts from what is painted is a readout that lies. {@code MsaHitTestTest}
     * pins the two together by hit-testing the centre of a painted cell.
     * <p>
     * Cheap by construction: the column is arithmetic, and the row is the same order of scan the panel already runs
     * on every mouse move in {@link #findNode}.
     */
    MsaCell msaCellAt(final int x, final int y) {
        if (!msaShown() || isOverFloatingOverlay(x, y)) {
            return null;
        }
        final float cw = getOptions().getMsaColumnWidth();
        if (cw <= 0) {
            return null;
        }
        final float origin_x = annotationColumnsEndX() + MSA_TRACK_GAP;
        final int visible = msaVisibleColumns();
        // The paint places cell i at Math.round(origin_x + i*cw); origin_x is fractional, so a plain division does
        // NOT invert it -- rounding shifts a boundary by up to a pixel, and that pixel column would resolve to the
        // neighbouring cell (naming the wrong residue). Start from the division, then settle on the cell whose
        // ROUNDED span actually contains x.
        int i = -1;
        final int guess = (int) Math.floor((x - origin_x) / cw);
        for (int cand = guess - 1; cand <= (guess + 1); ++cand) {
            if ((cand < 0) || (cand >= visible)) {
                continue;
            }
            final int left = Math.round(origin_x + (cand * cw));
            final int right = Math.round(origin_x + ((cand + 1) * cw));
            if ((x >= left) && (x < right)) {
                i = cand;
                break;
            }
        }
        if (i < 0) {
            return null;
        }
        final int column = msaColumnOffset() + i;
        if ((column < 0) || (column >= alignmentLength())) {
            return null;
        }
        final float pad = getYdistance();
        for (final PhylogenyNode tip : visibleExternalTips()) {
            if (!tip.getNodeData().isHasSequence()
                    || !tip.getNodeData().getSequence().isMolecularSequenceAligned()) {
                continue;
            }
            final int cy = Math.round(tip.getYcoord() - pad);
            final int cell_h = Math.max(1, Math.round(tip.getYcoord() + pad) - cy);
            if ((y >= cy) && (y < (cy + cell_h))) {
                final String mol = tip.getNodeData().getSequence().getMolecularSequence();
                return ForesterUtil.isEmpty(mol) || (column >= mol.length()) ? null : new MsaCell(tip, column);
            }
        }
        return null;
    }

    /**
     * The alignment display's residue readout: hovering a cell names the alignment column, the residue's own
     * ungapped number within that sequence, and what the residue IS (letter, full name, colour class, and for an
     * amino acid its Kyte-Doolittle hydropathy).
     * <p>
     * A Swing tooltip rather than a painted overlay ON PURPOSE. Swing asks for this text only when it is about to
     * show or update a tooltip, and draws it in its own popup -- so the canvas is never repainted. A hover repaint
     * costs 2-3.6 ms (measured for the focus glow); assembling this text costs microseconds.
     */
    @Override
    public String getToolTipText(final java.awt.event.MouseEvent event) {
        if (event == null) {
            return null;
        }
        final MsaCell cell = msaCellAt(event.getX(), event.getY());
        if (cell == null) {
            return null; // not over the alignment -> no tooltip at all (the canvas has no other tooltip)
        }
        return ResidueInfo.describeCell(cell._tip.getName(),
                                        cell._tip.getNodeData().getSequence().getMolecularSequence(),
                                        cell._column,
                                        msaIsNucleotide());
    }

    /** The MSA window's total pixel width RESERVED to the right of the annotation columns (0 when not shown). The
     *  window is bounded to a share of the viewport so a LONG alignment never drags the tree off-screen -- a short
     *  alignment shows in full, a longer one shows a scrollable window (its own {@link #_msa_scrollbar}). */
    private int msaTrackWidth() {
        if (!msaShown()) {
            return 0;
        }
        final int full_px = Math.round(alignmentLength() * getOptions().getMsaColumnWidth());
        return MSA_TRACK_GAP + Math.min(full_px, msaBandBudgetPx());
    }

    /** The maximum pixel width the MSA window may occupy: a share of the scroll-pane viewport (the tree keeps the rest).
     *  Viewport-derived, so it is stable across a layout pass -- independent of the tree's own depth reservation. */
    private int msaBandBudgetPx() {
        int viewport_w = 0;
        if ((getMainPanel() != null) && (getMainPanel().getSizeOfViewport() != null)) {
            viewport_w = getMainPanel().getSizeOfViewport().width;
        }
        if (viewport_w <= 0) {
            viewport_w = (getWidth() > 0) ? getWidth() : 800; // fallback before the viewport is realized
        }
        return Math.max(MSA_MIN_BAND_PX, (int) Math.round(viewport_w * MSA_MAX_VIEWPORT_FRACTION));
    }

    /** How many alignment columns fit in the drawn window (the reserved band minus the leading gap). */
    private int msaVisibleColumns() {
        if (!msaShown()) {
            return 0;
        }
        final float cw = getOptions().getMsaColumnWidth();
        if (cw <= 0) {
            return alignmentLength();
        }
        return Math.max(1, Math.min(alignmentLength(), (int) Math.floor((msaTrackWidth() - MSA_TRACK_GAP) / cw)));
    }

    /** True when the alignment is wider than its window, so the dedicated MSA scroller is needed. */
    private boolean isMsaScrollable() {
        return msaShown() && (alignmentLength() > msaVisibleColumns());
    }

    /** Height (px) of the column-position ruler band (a line + ticks + one label row), reserved just below the tips. */
    private int msaRulerBandHeight() {
        return getFontMetrics(getTreeFontSet().getSmallFont()).getHeight() + MSA_RULER_TICK_LEN + 4;
    }

    /** Breadth (px) reserved at the bottom for the MSA column ruler, so it sits clear just below the last tip (added to
     *  the breadth budget + extent alongside {@link #scaleAxisBottomReserve()}). 0 unless the alignment is shown. */
    private int msaRulerReserve() {
        return msaShown() ? msaRulerBandHeight() : 0;
    }

    /** Whether the conservation/consensus track is drawn: the alignment is shown AND the option is on. */
    private boolean msaConservationShown() {
        return msaShown() && getOptions().isShowMsaConservation();
    }

    /** The consensus letters are drawn by the same rule the residue letters use -- only when a column is wide
     *  enough to read one. The band height follows, so no space is reserved for letters that are not drawn. */
    private boolean msaConsensusLettersDrawn() {
        return getOptions().getMsaColumnWidth() >= MSA_LETTER_MIN_WIDTH;
    }

    /** Height of the bar area alone. ONE definition, used by the reserve below AND by the painter: were the two
     *  to drift, the band would either overlap the last alignment row or float clear of the ruler. */
    private int msaConservationBarHeight() {
        return Math.max(AptxConstants.MSA_CONSERVATION_BAR_HEIGHT_MIN,
                getFontMetrics(getTreeFontSet().getSmallFont()).getHeight() * 2);
    }

    /** Height of the conservation band: the bar area plus, when they are drawn, a row for the consensus letters. */
    private int msaConservationBandHeight() {
        final int fh = getFontMetrics(getTreeFontSet().getSmallFont()).getHeight();
        return msaConservationBarHeight() + (msaConsensusLettersDrawn() ? fh : 0)
                + AptxConstants.MSA_CONSERVATION_TOP_GAP;
    }

    private int msaConservationReserve() {
        return msaConservationShown() ? msaConservationBandHeight() : 0;
    }

    /** The per-column conservation over the tips on screen, computed on first use after a navigation. */
    private MsaConservation msaConservation() {
        if (_msa_conservation == null) {
            final java.util.List<String> rows = new java.util.ArrayList<String>();
            for (final PhylogenyNode tip : visibleExternalTips()) {
                if (tip.getNodeData().isHasSequence()) {
                    final Sequence s = tip.getNodeData().getSequence();
                    if (s.isMolecularSequenceAligned() && !ForesterUtil.isEmpty(s.getMolecularSequence())) {
                        rows.add(s.getMolecularSequence());
                    }
                }
            }
            _msa_conservation = MsaConservation.compute(rows, alignmentLength(), msaIsNucleotide());
        }
        return _msa_conservation;
    }

    /** For tests: the score the track would draw for a 0-based alignment column. */
    double msaConservationScoreForTest(final int col) {
        return msaConservation().scoreAt(col, getOptions().getMsaConservationMeasure());
    }

    /** For tests: the consensus residue the track would draw for a 0-based alignment column (0 when none). */
    char msaConsensusForTest(final int col) {
        return msaConservation().consensusAt(col);
    }

    /** For tests: the vertical space the conservation band takes from the tree (0 when it is not drawn). */
    int msaConservationReserveForTest() {
        return msaConservationReserve();
    }

    /** A "nice" tick step (1/2/5 x 10^k) for a column ruler spanning {@code span} columns (~5-8 ticks). */
    private static int niceColumnStep(final int span) {
        final double raw = Math.max(1, span) / 6.0;
        if (raw <= 1) {
            return 1;
        }
        final double mag = Math.pow(10, Math.floor(Math.log10(raw)));
        final double norm = raw / mag;
        final double nice = (norm <= 1) ? 1 : (norm <= 2) ? 2 : (norm <= 5) ? 5 : 10;
        return Math.max(1, (int) Math.round(nice * mag));
    }

    /** The first alignment column shown in the window, clamped so the window never runs past the last column. */
    private int msaColumnOffset() {
        final int max_off = Math.max(0, alignmentLength() - msaVisibleColumns());
        if (_msa_col_offset > max_off) {
            _msa_col_offset = max_off;
        }
        if (_msa_col_offset < 0) {
            _msa_col_offset = 0;
        }
        return _msa_col_offset;
    }

    /** Package hook for the scroller listener (and tests): set the first shown column. */
    void setMsaColumnOffset(final int offset) {
        _msa_col_offset = offset;
    }

    int getMsaColumnOffset() {
        return msaColumnOffset();
    }

    boolean isMsaScrollableForTest() {
        return isMsaScrollable();
    }

    int getMsaVisibleColumnsForTest() {
        return msaVisibleColumns();
    }

    /** Whether the real START-of-alignment boundary line is drawn now (column 0 is in the window). */
    boolean msaStartBoundaryVisibleForTest() {
        return msaShown() && (msaColumnOffset() <= 0);
    }

    /** Whether the real END-of-alignment boundary line is drawn now (the last column is in the window). */
    boolean msaEndBoundaryVisibleForTest() {
        return msaShown() && ((msaColumnOffset() + msaVisibleColumns()) >= alignmentLength());
    }

    int msaRulerReserveForTest() {
        return msaRulerReserve();
    }

    static int niceColumnStepForTest(final int span) {
        return niceColumnStep(span);
    }

    /** Wires this tab's dedicated horizontal MSA scroller (created by MainPanel, shown under the tree canvas). */
    void setMsaScrollBar(final javax.swing.JScrollBar bar) {
        _msa_scrollbar = bar;
    }

    /** Syncs the MSA scroller to the current window (value = first column, extent = visible columns, max = total
     *  columns) and shows it only when the alignment is scrollable. A {@code setValues} with an unchanged value fires
     *  no adjustment event, so this can run from either the layout or the paint without a feedback loop. */
    void updateMsaScrollBar() {
        if (_msa_scrollbar == null) {
            return;
        }
        if (isMsaScrollable()) {
            final int total = alignmentLength();
            final int extent = msaVisibleColumns();
            final int value = msaColumnOffset();
            if ((_msa_scrollbar.getValue() != value) || (_msa_scrollbar.getVisibleAmount() != extent)
                    || (_msa_scrollbar.getMaximum() != total)) {
                _msa_scrollbar.setValues(value, extent, 0, total);
                _msa_scrollbar.setUnitIncrement(1);
                _msa_scrollbar.setBlockIncrement(Math.max(1, extent - 1));
            }
            if (!_msa_scrollbar.isVisible()) {
                _msa_scrollbar.setVisible(true);
            }
        }
        else if (_msa_scrollbar.isVisible()) {
            _msa_scrollbar.setVisible(false);
        }
    }

    /** True iff the alignment reads as nucleotide (samples the first aligned tip); drives the residue color scheme. */
    private boolean msaIsNucleotide() {
        for (final PhylogenyNode ext : _phylogeny.getExternalNodes()) {
            if (ext.getNodeData().isHasSequence()) {
                final Sequence s = ext.getNodeData().getSequence();
                if (s.isMolecularSequenceAligned() && !ForesterUtil.isEmpty(s.getMolecularSequence())) {
                    return MsaColors.isNucleotide(s.getMolecularSequence());
                }
            }
        }
        return false;
    }

    /**
     * Draws the alignment as colored residue cells, one row per visible tip, starting just right of the annotation
     * columns. Rectangular root-left only (dispatched in the {@code !vertical} block and gated by {@link #msaShown()}).
     * Only the columns within the current clip are drawn, so a very wide alignment stays cheap. Residue letters are
     * overlaid once the column is wide enough to read. Rides every export path via {@code paintPhylogeny} (WYSIWYG).
     */
    private void paintMsaTrack(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file,
                               final int graphics_file_y, final int graphics_file_height) {
        updateMsaScrollBar(); // keep the dedicated scroller in sync (also hides it when off / short)
        if (!msaShown()) {
            return;
        }
        final int total = alignmentLength();
        final int offset = msaColumnOffset();   // the first column shown in the window
        final int visible = msaVisibleColumns(); // columns the window can show
        final float cw = getOptions().getMsaColumnWidth();
        final float origin_x = annotationColumnsEndX() + MSA_TRACK_GAP;
        // draw only window positions [first_i, last_i] (0-based within the window), clipped to the graphics clip too
        int first_i = 0;
        int last_i = visible - 1;
        final Rectangle clip = g.getClipBounds();
        if ((clip != null) && (cw > 0)) {
            first_i = Math.max(0, (int) Math.floor((clip.x - origin_x) / cw));
            last_i = Math.min(visible - 1, (int) Math.ceil(((clip.x + clip.width) - origin_x) / cw));
        }
        if (first_i > last_i) {
            return;
        }
        final boolean nucleotide = msaIsNucleotide();
        final float row_h = 2f * getYdistance();
        final boolean draw_letters = cw >= MSA_LETTER_MIN_WIDTH;
        Font letter_font = null;
        FontMetrics fm = null;
        if (draw_letters) {
            final int size = Math.max(6, Math.min(13, Math.round(Math.min(cw, row_h) * 0.8f)));
            letter_font = new Font(Font.MONOSPACED, Font.PLAIN, size);
            fm = getFontMetrics(letter_font);
        }
        final Color saved_color = g.getColor();
        final Font saved_font = g.getFont();
        final Stroke saved_stroke = g.getStroke();
        // A gap is drawn as a faint horizontal line (not blank) so the alignment's extent -- where each row's aligned
        // sequence begins/ends -- is visible; consecutive gaps join into a continuous dash.
        final Color fg = getTreeColorSet().getSequenceColor();
        final Color gap_line_color = new Color(fg.getRed(), fg.getGreen(), fg.getBlue(), 38);
        final Color boundary_color = new Color(fg.getRed(), fg.getGreen(), fg.getBlue(), 165);
        g.setStroke(STROKE_1);
        int y_top = Integer.MAX_VALUE; // vertical span of the aligned rows (for the start/end boundary lines)
        int y_bottom = Integer.MIN_VALUE;
        for (final PhylogenyNode tip : visibleExternalTips()) {
            if (!tip.getNodeData().isHasSequence()) {
                continue;
            }
            final Sequence s = tip.getNodeData().getSequence();
            if (!s.isMolecularSequenceAligned()) {
                continue;
            }
            final String mol = s.getMolecularSequence();
            if (ForesterUtil.isEmpty(mol)) {
                continue;
            }
            final float pad = getYdistance();
            final int cy = Math.round(tip.getYcoord() - pad);
            final int cell_h = Math.max(1, Math.round(tip.getYcoord() + pad) - cy);
            y_top = Math.min(y_top, cy);
            y_bottom = Math.max(y_bottom, cy + cell_h);
            for (int i = first_i; i <= last_i; ++i) {
                final int c = offset + i; // the actual alignment column at window position i
                if ((c >= total) || (c >= mol.length())) {
                    break; // ran out of columns / this tip's sequence
                }
                final char res = mol.charAt(c);
                final int cx = Math.round(origin_x + (i * cw)); // window position i, so the window stays put
                final int cell_w = Math.max(1, Math.round(origin_x + ((i + 1) * cw)) - cx);
                final Color col = MsaColors.colorFor(res, nucleotide);
                if (col == null) {
                    final int mid = cy + (cell_h / 2); // a gap -> faint horizontal line spanning the cell
                    g.setColor(gap_line_color);
                    g.drawLine(cx, mid, cx + cell_w, mid);
                    continue;
                }
                g.setColor(col);
                g.fillRect(cx, cy, cell_w, cell_h);
                if (draw_letters && (cell_h >= (fm.getAscent() + fm.getDescent()))) {
                    final String ch = String.valueOf(Character.toUpperCase(res));
                    final int lw = fm.stringWidth(ch);
                    g.setFont(letter_font);
                    g.setColor(MsaColors.letterColor(col));
                    g.drawString(ch, cx + ((cell_w - lw) / 2), cy + (((cell_h + fm.getAscent()) - fm.getDescent()) / 2));
                }
            }
        }
        // Boundary lines at the REAL start/end of the alignment -- drawn only when that edge is actually in the window
        // (offset 0 shows the start; the window reaching the last column shows the end), so the user can tell a true
        // edge from a scroll cutoff. A middle-scrolled window shows neither.
        if (y_bottom > y_top) {
            g.setColor(boundary_color);
            if (offset <= 0) {
                final int sx = Math.round(origin_x);
                g.drawLine(sx, y_top, sx, y_bottom);
            }
            if ((offset + visible) >= total) {
                final int ex = Math.round(origin_x + ((total - offset) * cw));
                g.drawLine(ex, y_top, ex, y_bottom);
            }
        }
        paintMsaConservation(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height, origin_x, cw, offset,
                visible, total);
        paintMsaRuler(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height, origin_x, cw, offset, visible,
                total);
        g.setColor(saved_color);
        g.setFont(saved_font);
        g.setStroke(saved_stroke);
    }

    /**
     * The conservation/consensus track: one bar per alignment column in the reserved band between the alignment and
     * the column ruler, with the consensus residue under it once the columns are wide enough to read one.
     * <p>
     * Scored over the tips CURRENTLY ON SCREEN (see {@link #msaConservation()}), so entering a subtree or collapsing
     * a clade re-scores the profile for what is actually being looked at -- the same policy the "Color by" scheme and
     * the annotation columns already follow. Which of the two measures the bar shows is a Setting; the measures
     * themselves, and what gaps do to them, are documented on {@link MsaConservation}.
     * <p>
     * Floats with the ruler at the viewport bottom on screen (the alignment rows are the TIPS, so a band anchored
     * under the last one would scroll out of view) and is tree-anchored in an export, exactly like the ruler.
     */
    private void paintMsaConservation(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file,
                                      final int graphics_file_y, final int graphics_file_height, final float origin_x,
                                      final float cw, final int offset, final int visible, final int total) {
        if (!msaConservationShown()) {
            return;
        }
        final MsaConservation cons = msaConservation();
        if (cons.rows() < 1) {
            return;
        }
        final Rectangle vr = getVisibleRect();
        final int canvas_bottom = TreePanelUtil.scaleAxisFloatingBottom(to_pdf, to_graphics_file, graphics_file_y,
                graphics_file_height, getHeight(), vr.y + vr.height);
        // sits directly on top of the ruler band, which itself sits on top of the distance-scale band
        final int band_bottom = (canvas_bottom - scaleAxisBottomReserve() - msaRulerBandHeight()) - 1;
        final Font small = getTreeFontSet().getSmallFont();
        final FontMetrics fm = getFontMetrics(small);
        final boolean letters = msaConsensusLettersDrawn();
        final int bar_bottom = band_bottom - (letters ? fm.getHeight() : 0);
        final int bar_h = msaConservationBarHeight(); // the same value the layout reserve is built from
        final Color saved_color = g.getColor();
        final Font saved_font = g.getFont();
        final Stroke saved_stroke = g.getStroke();
        // The axis ink (already black on a black-and-white export), softened so the bars read as a chart rather
        // than as text -- greyscale chart furniture, like the scale grid.
        final Color ink = scaleInkColor(to_pdf, to_graphics_file);
        final Color bar_color = new Color(ink.getRed(), ink.getGreen(), ink.getBlue(), 170);
        g.setStroke(STROKE_1);
        final int shown = Math.min(visible, total - offset);
        // A very faint full-height track behind the bars, so the chart has a visible scale: on a well-conserved
        // alignment nearly every bar is near the top, and without the empty part showing they read as one solid
        // block rather than as "0.95 and 0.83".
        final int track_x0 = Math.round(origin_x);
        final int track_w = Math.round(origin_x + (shown * cw)) - track_x0;
        // ...but not on a transparent-background export, where a wide translucent fill would defeat the cut-out
        // (the same call the zebra stripes make). The bars themselves still draw; they just lose their backdrop.
        if ((track_w > 0) && !(to_graphics_file && _export_transparent_background)) {
            g.setColor(new Color(ink.getRed(), ink.getGreen(), ink.getBlue(), 26));
            g.fillRect(track_x0, bar_bottom - bar_h, track_w, bar_h);
        }
        for (int i = 0; i < shown; ++i) {
            final int c = offset + i;
            final double score = cons.scoreAt(c, getOptions().getMsaConservationMeasure());
            final int cx = Math.round(origin_x + (i * cw));
            final int cell_w = Math.max(1, Math.round(origin_x + ((i + 1) * cw)) - cx);
            final int h = (int) Math.round(score * bar_h);
            if (h > 0) {
                g.setColor(bar_color);
                g.fillRect(cx, bar_bottom - h, cell_w, h);
            }
            if (letters) {
                final char res = cons.consensusAt(c);
                if (res != 0) {
                    g.setFont(small);
                    g.setColor(ink);
                    final String ch = String.valueOf(res);
                    g.drawString(ch, cx + ((cell_w - fm.stringWidth(ch)) / 2), bar_bottom + fm.getAscent());
                }
            }
        }
        // a baseline under the bars, so a column of zero conservation still reads as a measured zero, not a gap in
        // the chart -- and so the bar heights have something to be read against
        g.setColor(ink);
        drawLine(Math.round(origin_x), bar_bottom, Math.round(origin_x + (shown * cw)), bar_bottom, g);
        g.setColor(saved_color);
        g.setFont(saved_font);
        g.setStroke(saved_stroke);
    }

    /**
     * A column-position ruler in the reserved band just below the alignment: a baseline spanning the visible window with
     * ticks + 1-based column numbers at "nice" intervals, plus the first (1) and last (total) columns whenever they are
     * in view -- so the user can read WHICH columns the scrolled window is showing. Floats at the viewport bottom on
     * screen (always visible, like the distance scale axis), anchored to the tree/export bottom in an export (WYSIWYG),
     * and sits above the distance-scale band if that is also shown.
     */
    private void paintMsaRuler(final Graphics2D g, final boolean to_pdf, final boolean to_graphics_file,
                               final int graphics_file_y, final int graphics_file_height, final float origin_x,
                               final float cw, final int offset, final int visible, final int total) {
        final Rectangle vr = getVisibleRect();
        final int canvas_bottom = TreePanelUtil.scaleAxisFloatingBottom(to_pdf, to_graphics_file, graphics_file_y,
                graphics_file_height, getHeight(), vr.y + vr.height);
        final int ruler_y = (canvas_bottom - scaleAxisBottomReserve() - msaRulerBandHeight()) + 2;
        g.setFont(getTreeFontSet().getSmallFont());
        final FontMetrics rfm = g.getFontMetrics();
        g.setColor(scaleInkColor(to_pdf, to_graphics_file));
        g.setStroke(STROKE_1);
        final int shown = Math.min(visible, total - offset); // columns actually in the window
        final int base_x0 = Math.round(origin_x);
        final int base_x1 = Math.round(origin_x + (shown * cw));
        drawLine(base_x0, ruler_y, base_x1, ruler_y, g); // the ruler baseline, spanning the visible columns
        final int label_baseline = ruler_y + MSA_RULER_TICK_LEN + rfm.getAscent() + 1;
        final int step = niceColumnStep(shown);
        final int first_1b = offset + 1;
        final int last_1b = offset + shown;
        final java.util.TreeSet<Integer> marks = new java.util.TreeSet<Integer>();
        if (offset <= 0) {
            marks.add(1); // the true start
        }
        if ((offset + visible) >= total) {
            marks.add(total); // the true end
        }
        for (int p = ((first_1b + step - 1) / step) * step; p <= last_1b; p += step) {
            if (p >= 1) {
                marks.add(p);
            }
        }
        int last_label_right = Integer.MIN_VALUE; // decimate labels (keep all ticks) so numbers never overlap
        for (final int p : marks) {
            final int i = (p - 1) - offset; // window position of the 1-based column p
            final int cx = Math.round(origin_x + (i * cw) + (cw / 2f)); // the column centre
            drawLine(cx, ruler_y, cx, ruler_y + MSA_RULER_TICK_LEN, g);
            final String label = Integer.toString(p);
            final int half = rfm.stringWidth(label) / 2;
            if ((cx - half) >= (last_label_right + 4)) {
                g.drawString(label, cx - half, label_baseline);
                last_label_right = cx + half;
            }
        }
    }

    /**
     * Horizontal space needed to the RIGHT of the label reservation ({@link #getLongestExtNodeInfo()}) so
     * layout / "fit to window" / zoom keep it all on-screen: the tip-aligned annotation columns plus, beyond
     * them, any clade marks (boxes add the small pad; bars/brackets add the gap, the mark, and the rotated
     * label whose horizontal extent is the font height). 0 when neither columns nor bands are shown.
     */
    private int rightMarginExtraWidth() {
        int extra = annotationColumnsWidth() + msaTrackWidth();
        if (hasCladeBands()) {
            if (_clade_bands_mode == CLADE_VIS.BOXES) {
                extra += CLADE_BAND_RIGHT_PAD;
            } else if (drawnCladeBandCount() > 0) { // no bar/bracket reserve if every band is a skipped single-member clade
                final int label_h = getFontMetricsForLargeDefaultFont().getHeight();
                extra += CLADE_BAND_RIGHT_PAD;
                for (final CladeLevel level : drawableCladeLevels()) { // every nested column needs its own room
                    extra += cladeLevelDepthExtent(level, label_h);
                }
            }
        }
        return extra;
    }

    /** The depth-axis (horizontal, in the root-left layout) footprint of a clade bar/bracket taxon label, given its
     *  angle: HORIZONTAL reserves the full label width, DIAGONAL ~width/sqrt2, VERTICAL just one line height.
     *  Root-top/bottom labels are always upright, so their depth footprint is one line height. */
    private int cladeLabelDepthExtent(final CladeLevel level, final int label_h) {
        if (isVerticalOrientation()) {
            return label_h;
        }
        switch (level.getLabelAngle()) {
            case HORIZONTAL:
                return Math.max(label_h, maxCladeLabelWidth(level));
            case DIAGONAL:
                return Math.max(label_h, Math.round((maxCladeLabelWidth(level) + label_h) / 1.41f));
            default: // VERTICAL
                return label_h;
        }
    }

    /**
     * Where each bar/bracket level's mark COLUMN starts, finest first (nearest the tips) -- the single cursor the
     * nested-column geometry is defined by. Both draw loops walk it and {@link #cladeLevelStartsForTest} reports
     * it, so what is drawn and what is asserted cannot drift apart.
     */
    private float[] cladeLevelStarts() {
        final java.util.List<CladeLevel> levels = drawableCladeLevels();
        final int label_h = getFontMetricsForLargeDefaultFont().getHeight();
        final float[] xs = new float[levels.size()];
        float x = cladeBandRightEdge();
        for (int i = 0; i < levels.size(); ++i) {
            xs[i] = x + CLADE_BAR_GAP;
            x += cladeLevelDepthExtent(levels.get(i), label_h);
        }
        return xs;
    }

    /** Test hook: the depth (x) each nested clade column starts at, finest first. */
    final float[] cladeLevelStartsForTest() {
        return cladeLevelStarts();
    }

    /**
     * The FULL depth (horizontal, in the root-left layout) that one bar/bracket level occupies: the gap in front
     * of it, the mark itself, its label, and a little air before the next level. The draw loops advance their
     * cursor by exactly this and {@link #rightMarginExtraWidth} reserves exactly the sum of it, so a level can
     * never be drawn into space the layout did not set aside (which would clip the outermost column on a
     * fit-to-window).
     */
    private int cladeLevelDepthExtent(final CladeLevel level, final int label_h) {
        final int mark = (_clade_bands_mode == CLADE_VIS.BARS) ? (CLADE_BAR_WIDTH + 3) : 4;
        return CLADE_BAR_GAP + mark + cladeLabelDepthExtent(level, label_h) + 4;
    }

    /** Test hook: the horizontal reserve past the labels (annotation columns + any clade marks and their labels). */
    int cladeBandRightReserveForTest() {
        return rightMarginExtraWidth();
    }

    /** The widest DRAWN clade taxon label at ONE level (in the large font), for reserving that level's
     *  horizontal/diagonal label space; 0 if none. Skipped single-member clades don't draw, so they don't reserve. */
    private int maxCladeLabelWidth(final CladeLevel level) {
        final FontMetrics fm = getFontMetricsForLargeDefaultFont();
        int max = 0;
        for (final CladeBand band : level.getBands()) {
            if (!skipCladeBand(band)) {
                max = Math.max(max, fm.stringWidth(band.getTaxon()));
            }
        }
        return max;
    }

    /** Whether {@code node} still hangs off the displayed tree's root. A node deleted from the tree keeps its own
     *  children, so it looks perfectly normal in isolation -- only walking UP reveals that it is detached. */
    private boolean isInCurrentTree(final PhylogenyNode node) {
        if ((_phylogeny == null) || _phylogeny.isEmpty() || (node == null)) {
            return false;
        }
        final PhylogenyNode root = _phylogeny.getRoot();
        PhylogenyNode n = node;
        while (n != null) {
            if (n == root) {
                return true;
            }
            n = n.isRoot() ? null : n.getParent();
        }
        return false;
    }

    /** {@code {yTop, yBottom}} of a clade's tips in current paint coordinates, or null if none.
     *  <p>
     *  Returns null for a band root that is no longer part of the tree. That should not happen -- every structural
     *  change goes through {@link #afterTreeStructureChanged()}, which rebuilds the bands -- but walking a detached
     *  node's external descendants THROWS, and this runs inside the paint loop, where a throw repeats forever and
     *  the tree stops drawing entirely. A missing band is a far better failure than a dead window. */
    private float[] cladeBandYRange(final PhylogenyNode root) {
        if (!isInCurrentTree(root)) {
            return null;
        }
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
        final java.util.List<CladeLevel> levels = drawableCladeLevels(); // boxes: at most one (see drawableCladeLevels)
        if (levels.isEmpty()) {
            return;
        }
        final float right = cladeBandRightEdge();
        final float pad = getYdistance();
        for (final CladeBand band : levels.get(0).getBands()) {
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

    /**
     * The nested bar columns: one column per annotated rank, finest nearest the tips, each a solid colour bar per
     * clade plus its taxon label. The x cursor walks outward level by level, so the broader ranks sit outside the
     * finer ones they contain -- the nested-bracket look of a published figure. Every level's width contribution
     * is {@link #cladeLevelDepthExtent}, the same function the layout reserve uses, so what is drawn and what is
     * reserved can never disagree.
     */
    private void drawCladeBars(final Graphics2D g) {
        final float pad = getYdistance();
        final Font font = getTreeFontSet().getLargeFont();
        final java.util.List<CladeLevel> levels = drawableCladeLevels();
        final float[] starts = cladeLevelStarts();
        for (int li = 0; li < levels.size(); ++li) {
            final CladeLevel level = levels.get(li);
            final float bar_x = starts[li];
            for (final CladeBand band : level.getBands()) {
                final float[] yr = cladeBandYRange(band.getRoot());
                if (yr == null) {
                    continue;
                }
                if (skipCladeBand(band)) { // "Skip single-member clades": a single-tip clade's bar is a degenerate stub
                    continue;
                }
                final int y = Math.round(yr[0] - pad);
                final int h = Math.max(1, Math.round((yr[1] - yr[0]) + (2 * pad)));
                g.setColor(band.getColor());
                g.fillRect(Math.round(bar_x), y, CLADE_BAR_WIDTH, h); // rides R -> horizontal bar in a vertical orientation
                // taxon label: at this level's own angle in a horizontal layout (an outer rank has few, long names
                // worth reading straight; inner ones stay narrow), UPRIGHT in a vertical orientation
                g.setFont(font);
                g.setColor(getTreeColorSet().getSequenceColor());
                final FontMetrics fm = g.getFontMetrics();
                final float label_x = bar_x + CLADE_BAR_WIDTH + 3;
                final float mid_y = (yr[0] + yr[1]) / 2.0f;
                drawCladeBandLabel(g, band.getTaxon(), label_x, mid_y, fm, level.getLabelAngle());
            }
        }
    }

    // Monochrome "named clade" annotation: a thin "]" bracket (vertical spine + short end-ticks pointing
    // toward the tips) per clade with the rotated taxon label, all in the foreground color -- no per-clade
    // colors and no legend, for use when color already encodes something else.
    /** The nested bracket columns -- the monochrome twin of {@link #drawCladeBars}: one "]" column per annotated
     *  rank, finest nearest the tips, the x cursor walking outward by the same {@link #cladeLevelDepthExtent}. */
    private void drawCladeBrackets(final Graphics2D g) {
        // less outward pad than the filled bars/boxes, so adjacent brackets keep a clear vertical gap
        // (one "]" per clade) instead of merging into a single continuous line
        final float pad = getYdistance() * 0.3f;
        final Font font = getTreeFontSet().getLargeFont();
        final Stroke saved_stroke = g.getStroke();
        g.setStroke(new BasicStroke(CLADE_BRACKET_STROKE));
        g.setColor(getTreeColorSet().getSequenceColor());
        final java.util.List<CladeLevel> levels = drawableCladeLevels();
        final float[] starts = cladeLevelStarts();
        for (int li = 0; li < levels.size(); ++li) {
            final CladeLevel level = levels.get(li);
            final float spine_x = starts[li];
            for (final CladeBand band : level.getBands()) {
                final float[] yr = cladeBandYRange(band.getRoot());
                if (yr == null) {
                    continue;
                }
                if (skipCladeBand(band)) { // "Skip single-member clades": a single-tip clade's bracket is a stub
                    continue;
                }
                final int x = Math.round(spine_x);
                final int y0 = Math.round(yr[0] - pad);
                final int y1 = Math.round(yr[1] + pad);
                // "]" : vertical spine plus short end-ticks pointing left, toward the tips
                g.drawLine(x, y0, x, y1);
                g.drawLine(x, y0, x - CLADE_BRACKET_TICK, y0);
                g.drawLine(x, y1, x - CLADE_BRACKET_TICK, y1);
                // taxon label: at this level's own angle (horizontal layout), UPRIGHT in a vertical orientation
                g.setFont(font);
                final FontMetrics fm = g.getFontMetrics();
                final float label_x = spine_x + 4;
                final float mid_y = (yr[0] + yr[1]) / 2.0f;
                drawCladeBandLabel(g, band.getTaxon(), label_x, mid_y, fm, level.getLabelAngle());
            }
        }
        g.setStroke(saved_stroke);
    }

    /** Draws a clade bar/bracket taxon label: rotated 90deg (reading bottom-to-top) beside the vertical mark in the
     *  horizontal layout, or UPRIGHT (re-anchored to the base frame, centered on the clade's breadth) beside the
     *  horizontal mark in a vertical orientation. Restores g's transform (R when vertical) before returning. */
    /**
     * The baseline for an UPRIGHT clade-band label in a vertical orientation, anchored at screen y
     * {@code anchor_y} -- the point just past the mark along the depth. The glyph box
     * {@code [baseline - ascent, baseline + descent]} must fall wholly on the far side of that anchor: BELOW it
     * with the root at the top, ABOVE it with the root at the bottom. Centring the text on the anchor instead
     * (what this replaced) left half of every name printed inside the coloured bar it labels.
     * <p>
     * Pure, so the rule can be checked without rendering.
     */
    static double uprightCladeLabelBaseline(final double anchor_y, final int ascent, final int descent,
                                            final boolean root_bottom) {
        return root_bottom ? (anchor_y - descent) : (anchor_y + ascent);
    }

    private void drawCladeBandLabel(final Graphics2D g, final String taxon, final float label_x, final float mid_y,
                                    final FontMetrics fm, final CLADE_LABEL_ANGLE angle) {
        final int tw = fm.stringWidth(taxon);
        final AffineTransform saved = g.getTransform();
        if (isVerticalOrientation()) {
            // root-top/bottom: always upright (the "horizontal" reading), re-anchored to the base frame.
            // label_x is the point just PAST the mark along the depth, so the text must be drawn BEYOND it --
            // below it in root-top, above it in root-bottom -- exactly like an upright tip label
            // (paintTipLabelHorizontal). Centring the text vertically on that anchor instead made every label
            // straddle the very bar it names, so half of "Ctenophora" sat inside the coloured band.
            final Point2D.Double lp = screenPoint(label_x, mid_y);
            final boolean root_bottom = getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_BOTTOM;
            g.setTransform(_orientation_base_transform);
            g.drawString(taxon, (float) (lp.x - (tw / 2.0)),
                    (float) uprightCladeLabelBaseline(lp.y, fm.getAscent(), fm.getDescent(), root_bottom));
        }
        else {
            switch (angle) {
                case HORIZONTAL: // reads left-to-right, beside the mark -- least vertical overlap, most horizontal space
                    g.drawString(taxon, label_x, mid_y + ((fm.getAscent() - fm.getDescent()) / 2.0f));
                    break;
                case DIAGONAL: // reads up-and-to-the-right from the mark (45 deg)
                    g.rotate(-Math.PI / 4.0, label_x, mid_y);
                    g.drawString(taxon, label_x + 2, mid_y + ((fm.getAscent() - fm.getDescent()) / 2.0f));
                    break;
                default: // VERTICAL (90 deg): reads bottom-to-top, centered on the clade -- the compact default
                    g.rotate(-Math.PI / 2.0, label_x, mid_y);
                    g.drawString(taxon, label_x - (tw / 2.0f), mid_y + fm.getAscent());
                    break;
            }
        }
        g.setTransform(saved);
    }

    final void decreaseDomainStructureEvalueThresholdExp() {
        if (_domain_structure_e_value_thr_exp > -20) {
            _domain_structure_e_value_thr_exp -= 1;
            AptxUtil.assignDomainPalette(null, getMainPanel()); // the drawn domain set changed -> recolour it
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
        // in a vertical orientation the node coords are logical (un-rotated); map the device click back to that space
        final Point2D.Double p = toLogicalPoint(x, y);
        for (final PhylogenyNodeIterator iter = _phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            // size the click target to the mark ACTUALLY drawn (per-node when styled) so a large custom node shape
            // is clickable across its full extent, not just the default-sized box
            final int half_box_size_plus_wiggle = effectiveNodeHalfBoxSize(node) + WIGGLE;
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
     * selection. Picks the nearest branch within a small tolerance (same forgiveness as {@link #findNode}). Each
     * layout hit-tests the branch geometry it actually DRAWS: RECTANGULAR the horizontal leg (at the child's y) plus
     * the vertical fork connector; TRIANGULAR / UNROOTED the straight parent-&gt;child line; CIRCULAR the radial leg
     * from the child to the point on the parent's radius at the child's angle. EURO_STYLE/ROUNDED offset the leg near
     * the parent, so branch-click stays a no-op there rather than a hit region that disagrees with the painted branch.
     * Only branches actually on screen (not hidden under a collapsed ancestor) are tested.
     */
    final PhylogenyNode findBranch(final int x, final int y) {
        final PHYLOGENY_GRAPHICS_TYPE gt = getPhylogenyGraphicsType();
        final boolean unrooted = (gt == PHYLOGENY_GRAPHICS_TYPE.UNROOTED);
        final boolean circular = (gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR);
        // "diagonal" = a straight parent->child line: TRIANGULAR and UNROOTED draw exactly that
        final boolean diagonal = (gt == PHYLOGENY_GRAPHICS_TYPE.TRIANGULAR) || unrooted;
        final boolean rectangular = (gt == PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR);
        if ((_phylogeny == null) || _phylogeny.isEmpty() || !(rectangular || diagonal || circular)) {
            return null;
        }
        final double tol = (getOptions().getDefaultNodeShapeSize() / 2.0) + WIGGLE;
        // in a vertical orientation the node coords are logical (un-rotated); map the device click back to that space
        final Point2D.Double click = toLogicalPoint(x, y);
        // circular: the root is at the ring centre; each node's angle is derivable from its coords (atan2), so the
        // radial leg reconstructs from coords WITHOUT the fragile _urt_nodeid_angle_map global
        final double root_x = circular ? _phylogeny.getRoot().getXcoord() : 0;
        final double root_y = circular ? _phylogeny.getRoot().getYcoord() : 0;
        PhylogenyNode best = null;
        double best_dist = tol;
        for (final PhylogenyNodeIterator iter = _phylogeny.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode n = iter.next();
            if (n.isCollapse()) {
                continue; // a collapsed clade's tips are hidden (not individually selectable)
            }
            double d = Double.MAX_VALUE;
            if (!n.isRoot()) {
                final PhylogenyNode p = n.getParent();
                if (circular) {
                    // radial leg: from n out to the point on the PARENT's radius at n's angle (what paintBranchCircular draws)
                    final double angle = Math.atan2(n.getYcoord() - root_y, n.getXcoord() - root_x);
                    final double pdx = p.getXcoord() - root_x, pdy = p.getYcoord() - root_y;
                    final double parent_radius = Math.sqrt((pdx * pdx) + (pdy * pdy));
                    d = Line2D.ptSegDist(n.getXcoord(), n.getYcoord(), root_x + (Math.cos(angle) * parent_radius),
                            root_y + (Math.sin(angle) * parent_radius), click.x, click.y);
                }
                else { // horizontal leg (at the child's y) or straight parent->child diagonal
                    final double a_y = diagonal ? p.getYcoord() : n.getYcoord();
                    d = Line2D.ptSegDist(p.getXcoord(), a_y, n.getXcoord(), n.getYcoord(), click.x, click.y);
                }
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
    // The FOCUS GLOW: concentric translucent discs, largest/faintest first, so the overlap accumulates toward the
    // centre and reads as a halo rather than as flat rings. Static on purpose -- hover is transient, so an animation
    // ramp would lag the pointer and barely play, and a pulsing hover would compete with the breathing halo that
    // already marks SEARCH HITS (paintFoundNodeHalos).
    private static final float[] HOVER_GLOW_RADII   = { 1.65f, 1.15f, 0.75f };
    private static final int[]   HOVER_GLOW_ALPHAS  = { 34, 44, 58 };
    private static final int     HOVER_GLOW_MIN_DIA = 18;                        // readable on a tiny default node
    private static final Color   HOVER_GLOW_ACCENT_FALLBACK = new Color(0x26, 0x75, 0xbf);

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

    /** Test/probe hook: drive the focus glow directly, without synthesizing mouse motion. */
    void setHoverForTest(final PhylogenyNode node, final boolean subtree) {
        setHover(node, subtree);
    }

    /** Applies a hover target, honoring the just-clicked suppression: the target just clicked stays un-glowed
     *  until the pointer moves off it (so it doesn't instantly flip to the "will be removed" grey). */
    private void applyHover(final PhylogenyNode node, final boolean subtree) {
        // The suppression exists ONLY to stop a just-selected node flipping to the grey "will be removed" under a
        // stationary pointer -- a selection-semantics problem. Outside Select-Node(s) mode the glow is a neutral
        // focus mark with nothing to flip, so honouring the suppression there would just blank the focus ring for
        // no reason (shift-click arms it in EVERY mode, not only Select).
        if ((node == _click_suppressed)
                && (getControlPanel().getActionWhenNodeClicked() == NodeClickAction.SELECT_NODES)) {
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
     * The FOCUS GLOW: a soft halo around the node under the pointer, drawn in EVERY click-to mode.
     * <p>
     * The hand cursor says "something here is clickable"; it never says WHICH node, which on a dense tree is the
     * thing you actually need -- most of all in the modes that reroot or delete. Before this, on-canvas hover
     * feedback existed in Select-Node(s) mode alone.
     * <p>
     * There is exactly ONE circle on the hovered node in every mode. In Select-Node(s) the glow CARRIES the
     * meaning the old flat preview disc used to (found colour = a click adds, grey = a click removes) instead of
     * drawing a second disc inside it; elsewhere it is a neutral accent, so a focus ring can never be mistaken for
     * a selection state. Hovering a BRANCH additionally keeps the flat per-tip marks, because one glow cannot say
     * "these forty tips".
     * <p>
     * Screen only -- hover is transient state and must never reach an export.
     */
    private void paintHoverPreview(final Graphics2D g, final boolean to_screen) {
        // _hover_node is always in the currently-displayed tree: it is only set on hover, and any tree swap
        // clears it via setNodeInPreorderToNull (the shared structural-change chokepoint).
        if (!to_screen || (_hover_node == null)) {
            return;
        }
        final boolean select_mode = getControlPanel().getActionWhenNodeClicked() == NodeClickAction.SELECT_NODES;
        final Set<Long> found = getFoundNodes0();
        boolean deselect = false;                // whether a click here would DESELECT rather than add
        java.util.List<PhylogenyNode> marks = null;
        if (select_mode) {
            if (_hover_subtree) {
                final java.util.List<PhylogenyNode> all = _hover_node.getAllExternalDescendants();
                if (all.isEmpty()) {
                    return;
                }
                deselect = allTipsSelected(all); // direction over ALL tips (matches what a click toggles)
                marks = new java.util.ArrayList<PhylogenyNode>();
                collectVisibleTips(_hover_node, marks); // only laid-out tips (skip collapsed sub-clades' stale coords)
            }
            else if (_hover_node.isCollapse()) {
                // A collapsed clade is one unit: mouseClicked routes a selection click on the triangle to
                // selectSubtreeTips, which toggles its HIDDEN TIPS -- the clade root's own id is never in the
                // found set. Reading the root here would show "a click adds" over a fully-selected clade that
                // the click is about to clear.
                deselect = allTipsSelected(_hover_node.getAllExternalDescendants());
            }
            else {
                deselect = (found != null) && found.contains(_hover_node.getId()); // a click toggles just this node
            }
        }
        final Color saved = g.getColor();
        paintFocusGlow(g, _hover_node, hoverGlowColor(select_mode, deselect));
        if (marks != null) {
            // the clade's tips: the glow marks the clade ROOT, these say which tips the click will take
            final Color f = getTreeColorSet().getFoundColor0();
            g.setColor(deselect ? BRANCH_HOVER_REMOVE
                    : new Color(f.getRed(), f.getGreen(), f.getBlue(), BRANCH_HOVER_ALPHA));
            final int shape = getOptions().getDefaultNodeShapeSize();
            // the remove mark lands on an already-selected node (with its solid found marker) so it is drawn larger
            // to cover it; the add mark sits on a bare node and reads fine smaller
            final int d = deselect ? Math.max(10, shape + 4) : Math.max(6, shape);
            for (final PhylogenyNode t : marks) {
                // a select marks only the not-yet-selected targets; a deselect marks all
                if (deselect || (found == null) || !found.contains(t.getId())) {
                    g.fillOval(Math.round(t.getXcoord()) - (d / 2), Math.round(t.getYcoord()) - (d / 2), d, d);
                }
            }
        }
        g.setColor(saved);
    }

    /** The colour the focus glow takes -- see {@link #paintHoverPreview}. */
    private Color hoverGlowColor(final boolean select_mode, final boolean deselect) {
        if (select_mode) {
            if (deselect) {
                return BRANCH_HOVER_REMOVE; // "a click here will de-select"
            }
            return getTreeColorSet().getFoundColor0(); // "a click here will select"
        }
        final Color accent = UIManager.getColor("Component.accentColor");
        return (accent != null) ? accent : HOVER_GLOW_ACCENT_FALLBACK;
    }

    /**
     * Draws the halo itself: concentric translucent discs, largest and faintest first, so their overlap builds up
     * toward the centre. The node keeps showing through -- this marks a node, it must not hide it.
     */
    private void paintFocusGlow(final Graphics2D g, final PhylogenyNode node, final Color base) {
        final int dia = Math.max(HOVER_GLOW_MIN_DIA, getOptions().getDefaultNodeShapeSize() * 3);
        final int x = Math.round(node.getXcoord()), y = Math.round(node.getYcoord());
        for (int i = 0; i < HOVER_GLOW_RADII.length; i++) {
            g.setColor(new Color(base.getRed(), base.getGreen(), base.getBlue(), HOVER_GLOW_ALPHAS[i]));
            final int d = Math.round(dia * HOVER_GLOW_RADII[i]);
            g.fillOval(x - (d / 2), y - (d / 2), d, d);
        }
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
        if (_hover_subtree || _hover_node.isCollapse()) { // a collapsed clade toggles its tips -- see paintHoverPreview
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
        if (shows(DisplayOption.USE_STYLE) && (PhylogenyMethods.getBranchColorValue(node) != null)) {
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
            AptxUtil.assignDomainPalette(null, getMainPanel()); // the drawn domain set changed -> recolour it
        }
    }

    /** Set the initial domain-architecture horizontal scale so the WIDEST architecture spans roughly {@code fraction}
     *  of {@code viewport_width} pixels -- called when the domain display is first fitted on load, so domains open at a
     *  useful, screen-proportional size instead of a fixed ~90px. (The final drawn width is this times the ~0.9 fit
     *  headroom applied in {@link #initNodeData}.) The user's later domain zoom (+/-) then scales from here. */
    void fitDomainWidthToScreen(final double fraction, final int viewport_width) {
        if (viewport_width > 0) {
            _domain_structure_width = fraction * viewport_width;
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
                    rds = new RenderableDomainArchitecture(node.getNodeData().getSequence().getDomainArchitecture());
                    node.getNodeData().getSequence().setDomainArchitecture(rds);
                } else {
                    rds = (RenderableDomainArchitecture) node.getNodeData().getSequence().getDomainArchitecture();
                }
                if (shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES)) {
                    final double dsw = rds.getOriginalSize().getWidth();
                    if (dsw > _max_original_domain_structure_width) {
                        _max_original_domain_structure_width = dsw;
                    }
                }
            }
        }
        if (shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES)) {
            final float ds_factor_width = (float) (effectiveDomainStructureWidth() / _max_original_domain_structure_width);
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
        if (shows(DisplayOption.NODE_DATA_POPUP)) {
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
                    if (shows(DisplayOption.NODE_DATA_POPUP)) {
                        showNodeDataPopup(e, node);
                    }
                }
                // The focus glow marks the node under the pointer in EVERY mode -- the cursor says "clickable",
                // the glow says WHICH. Collapsed clades included: they are clickable everywhere too, and a soft
                // halo behind the triangle reads cleanly where the old flat disc over it did not.
                applyHover(node, false);
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

    /** The radial square-canvas side (px) that fits the current viewport at zoom 1: the smaller viewport dimension,
     *  with a floor so a tiny/again-hidden panel still yields a sane value. */
    private int radialFitDiameter() {
        final int d = Math.min(getVisibleRect().width, getVisibleRect().height);
        return (d >= MIN_RADIAL_DIAMETER) ? d : DEFAULT_RADIAL_DIAMETER;
    }

    /** The current radial canvas side, lazily initialised to the fit diameter (so a freshly-shown radial tree fits). */
    final int radialDiameter() {
        if (_radial_diameter <= 0) {
            _radial_diameter = radialFitDiameter();
        }
        return _radial_diameter;
    }

    /** Re-fit the radial layout to an explicit viewport size (used by showWhole, where getVisibleRect may lag). */
    final void fitRadialTo(final int viewport_w, final int viewport_h) {
        final int d = Math.min(viewport_w, viewport_h);
        _radial_diameter = (d >= MIN_RADIAL_DIAMETER) ? d : DEFAULT_RADIAL_DIAMETER;
    }

    /** Scale the radial canvas -- the single radial zoom -- clamped so it can neither vanish nor blow up the preferred
     *  size. Both layouts follow: the circular radius is half the (square) preferred size, the unrooted spread is
     *  scaled by an urt-factor derived from it (see setUpUrtFactor). */
    final void multiplyRadialDiameter(final float f) {
        long d = Math.round(radialDiameter() * (double) f);
        if (d < MIN_RADIAL_DIAMETER) {
            d = MIN_RADIAL_DIAMETER;
        }
        if (d > MAX_RADIAL_DIAMETER) {
            d = MAX_RADIAL_DIAMETER;
        }
        _radial_diameter = (int) d;
    }

    /** Mark the radial canvas for a lazy re-fit on the next paint (used when switching TO a radial layout). */
    final void invalidateRadialDiameter() {
        _radial_diameter = 0;
    }

    /**
     * Lays the tree out to FILL a fixed-size export canvas of {@code w} x {@code h} points (see
     * {@link AptxUtil#writePhylogenyToGraphicsFileAtSize} and {@link ExportSizeSpec}), for all five display types:
     * the rectangular family via the depth/breadth machinery in {@link #calcParametersForPainting}, and
     * circular/unrooted via the radial diameter. Returns an opaque token capturing the prior on-screen layout (the
     * component size, and the radial zoom) to hand to {@link #restoreLayoutAfterExport}. MUST be paired with a
     * restore (in a finally), so a fixed-size export can never leave the on-screen tree laid out at the export size.
     */
    final int[] layoutForExportSize(final int w, final int h) {
        final int[] prior = { getWidth(), getHeight(), _radial_diameter };
        if (isRadialLayout()) {
            // the radial canvas is a SQUARE driven by _radial_diameter; fit it to the export frame BEFORE
            // calcParametersForPainting so setUpUrtFactor (invoked from there) reads the fitted diameter
            fitRadialTo(w, h);
        }
        calcParametersForPainting(w, h);
        return prior;
    }

    /** Restores the on-screen layout captured by {@link #layoutForExportSize} and repaints: the radial zoom from the
     *  captured diameter, and the layout by re-fitting at the component size. A fixed-size (usually narrower) export
     *  can auto-shrink the shared displayed font, and the auto-fit is hysteretic (it shrinks and grows at different
     *  thresholds), so re-fitting from the shrunk size could settle on a different font than the on-screen one; reset
     *  the font to the user/base size FIRST so the re-fit approaches from the base -- the canonical on-screen state. */
    final void restoreLayoutAfterExport(final int[] prior) {
        if (prior == null) {
            return;
        }
        if (isRadialLayout()) {
            _radial_diameter = prior[2];
        }
        final TreeFontSet fonts = getMainPanel().getTreeFontSet();
        fonts.setUserFontSize(fonts.getUserFontSize()); // reset displayed->base + clear the export's auto-shrink
        calcParametersForPainting(prior[0], prior[1]);
        repaint();
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
        final double inward_x = root_x + (Math.cos(angle) * parent_radius);
        final double inward_y = root_y + (Math.sin(angle) * parent_radius);
        drawLine(c.getXcoord(), c.getYcoord(), inward_x, inward_y, g);
        // "Break Long Branches": mark a capped radial leg with a break glyph, rotated to the spoke, at 0.72 along the
        // leg (clear of the support/length numbers centred at the midpoint)
        if (breakLongBranchesActiveCircular() && (c.getDistanceToParent() > breakLongBranchCap())) {
            paintBranchBreakGlyph(g,
                    (float) (inward_x + ((c.getXcoord() - inward_x) * BRANCH_BREAK_GLYPH_FRACTION)),
                    (float) (inward_y + ((c.getYcoord() - inward_y) * BRANCH_BREAK_GLYPH_FRACTION)), angle,
                    to_graphics_file);
        }
        // support + branch-length numbers ride the middle of this radial leg (rectangular draws them via
        // paintConfidenceValues/paintBranchLength, which are horizontal-branch only)
        paintBranchDataRadial(g, c, (c.getXcoord() + inward_x) / 2.0, (c.getYcoord() + inward_y) / 2.0, angle, to_pdf,
                to_graphics_file);
        paintNodeBox(c.getXcoord(), c.getYcoord(), c, g, to_pdf, to_graphics_file);
        if (c.isCollapse()) {
            // a collapsed clade-root is a stub here (no box/subtree); draw its collapse marker (triangle + count),
            // opening outward along the node's ring angle. paintNodeDataUnrootedCirc below returns early for it.
            paintRadialCollapsedMarker(g, c, angle, to_pdf, to_graphics_file);
        }
        final boolean is_in_found_nodes = isInFoundNodes0(c) || isInFoundNodes1(c);
        if (c.isExternal()) {
            if ((_dynamic_hiding_factor > 1) && !is_in_found_nodes
                    && ((_urt_nodeid_index_map.get(c.getId()) % _dynamic_hiding_factor) != 1)) {
                return;
            }
            paintNodeDataUnrootedCirc(g, c, to_pdf, to_graphics_file, radial_labels, 0, is_in_found_nodes);
        } else {
            // internal-node label (clade names from rank annotation, node/seq names) rides the branch radially; the
            // node's angle is read from _urt_nodeid_angle_map inside. Gated on "Show Internal Data" inside the method.
            // Not dynamic-hiding-culled -- same as the rectangular layout, which also draws every internal label (a
            // shared, deferred perf/clutter concern on very large trees; zoom to declutter).
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
        _root = phy.getRoot();
        _root.setXcoord(center_x);
        _root.setYcoord(center_y);
        final boolean radial_labels = getOptions().getNodeLabelDirection() == NODE_LABEL_DIRECTION.RADIAL;
        // Assign each DISPLAYED tip an angle around the ring: a visible external OR a collapsed clade-ROOT (one
        // pseudo-tip each). A collapsed root MUST get an explicit angle + coords here because it is INTERNAL and so
        // is skipped by the paintCirculars pass below -- without it, reading its angle later NPE'd (the long-standing
        // circular-collapse crash). Hidden externals under a collapse are never reached (the walk stops at the
        // collapsed root), so they get no angle, as intended.
        final int num_tips = countCircularDisplayedTips(_root);
        final double angle_increment = (num_tips > 0) ? (TWO_PI / num_tips) : TWO_PI;
        assignCircularDisplayedTipAngles(_root, center_x, center_y, radius, angle_increment,
                new double[] { starting_angle }, new int[] { 0 });
        paintCirculars(phy.getRoot(), phy, center_x, center_y, radius, radial_labels, g, to_pdf, to_graphics_file);
        paintNodeBox(_root.getXcoord(), _root.getYcoord(), _root, g, to_pdf, to_graphics_file);
    }

    /** Number of DISPLAYED tips in a circular layout: a visible external, or a collapsed clade-root (a pseudo-tip;
     *  the walk stops there, so its hidden descendants are not counted). */
    private int countCircularDisplayedTips(final PhylogenyNode node) {
        if (node.isCollapse() || node.isExternal()) {
            return 1;
        }
        int n = 0;
        for (int i = 0; i < node.getNumberOfDescendants(); ++i) {
            n += countCircularDisplayedTips(node.getChildNode(i));
        }
        return n;
    }

    /** Whether the CIRCULAR layout is drawn as a PHYLOGRAM: a node's ring RADIUS then encodes its distance-from-root
     *  (branch lengths), like the unrooted layout, instead of topological depth. Requires "Draw Phylogram" on, branch
     *  lengths present, and a positive tree height. */
    private boolean isCircularPhylogram() {
        return (getPhylogenyGraphicsType() == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) && getControlPanel().isDrawPhylogram()
                && isPhyHasBranchLengths() && (getMaxDistanceToRoot() > 0);
    }

    /** The circular ALIGNED phylogram (the "A" tree-shape button in circular): a phylogram whose branches end at each
     *  tip's branch-length radius, but whose external tip LABELS are all pinned to the common OUTER ring (radius) with a
     *  dotted radial leader bridging the gap -- the iTOL aligned-tips signature look. The polar twin of the rectangular
     *  ALIGNED_PHYLOGRAM (labels at a common right column + leader). UNALIGNED ("P") keeps labels at each tip's radius. */
    private boolean isAlignedCircularPhylogram() {
        return isCircularPhylogram()
                && (getControlPanel().getTreeDisplayType() == Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM);
    }

    /** The device point where {@code node}'s circular tip label is anchored: the tip's OWN position, EXCEPT in the
     *  aligned circular phylogram where an external tip's label is pinned to the common outer ring at the tip's angle
     *  (radius = the tree ring), so all labels line up on a circle. The single source of the anchor shared by the label
     *  paint + leader ({@link #paintNodeDataUnrootedCirc}) and the render test, so they can never drift. The centre +
     *  radius are derived from {@link #getPreferredSize()} + {@link #circularRadius} -- IDENTICAL to how the enclosing
     *  paintCircular set the node coords in the same pass, so anchor and drawn tree agree on screen and in exports. */
    private Point2D.Double circularLabelAnchor(final PhylogenyNode node) {
        if (node.isExternal() && (_graphics_type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR) && isAlignedCircularPhylogram()
                && (_circular_radius > 0)) {
            final Double a = _urt_nodeid_angle_map.get(node.getId());
            if (a != null) {
                final double m = a % TWO_PI;
                // the ring centre + radius the enclosing paintCircular used this pass (screen: padded-panel centre +
                // zoom radius; export: export-canvas centre + radius), so the anchor matches the drawn tree exactly
                return new Point2D.Double(_circular_center_x + (_circular_radius * Math.cos(m)),
                        _circular_center_y + (_circular_radius * Math.sin(m)));
            }
        }
        return new Point2D.Double(node.getXcoord(), node.getYcoord());
    }

    /** Test hook: the circular tip-label anchor point (see {@link #circularLabelAnchor}). */
    Point2D.Double circularLabelAnchorForTest(final PhylogenyNode node) {
        return circularLabelAnchor(node);
    }

    /** A node's fraction [0,1] of the ring RADIUS in the circular layout: the root at the centre (0); in a PHYLOGRAM
     *  its distance-to-root over the tree height (branch-length scaled); in a CLADOGRAM the tips on the outer ring (1)
     *  and internal/collapsed nodes by topological depth. Used by the main paint AND the overview so both agree. */
    private double circularRadiusFraction(final PhylogenyNode node) {
        if (node.isRoot()) {
            return 0;
        }
        if (isCircularPhylogram()) {
            // "Break Long Branches": the radius encodes the CAPPED distance-from-root, normalised by the capped height,
            // so an outlier branch is drawn as a shorter radial leg and the informative rings reclaim the ring.
            final boolean cap = breakLongBranchesActiveCircular();
            final double dist = cap ? TreePanelUtil.cappedDistanceToRoot(node, breakLongBranchCap())
                    : node.calculateDistanceToRoot();
            // root-EXCLUDED capped normalizer (matches cappedDistanceToRoot), so the deepest tip fills the ring exactly
            final double max = cap ? breakCappedRadialMax() : getMaxDistanceToRoot();
            final double f = (max > 0) ? (dist / max) : 0;
            return (f < 0) ? 0 : ((f > 1) ? 1 : f); // clamp a root-branch / rounding overshoot onto the ring
        }
        if (node.isExternal()) {
            return 1.0; // cladogram: all tips on the outer ring
        }
        return 1 - (((double) _circ_max_depth - node.calculateDepth()) / _circ_max_depth);
    }

    /** Assigns each displayed tip its ring angle (advancing {@code angle[0]} by {@code angle_inc}), in tree order:
     *  a visible external goes on the outer ring (or its distance radius in a phylogram); a collapsed clade-root is
     *  positioned at its own depth/distance radius (its incoming branch ends there, where its collapse marker is
     *  drawn -- see {@link #paintRadialCollapsedMarker}). {@code index[0]} is the external ordinal used for
     *  dynamic-hiding. */
    private void assignCircularDisplayedTipAngles(final PhylogenyNode node, final int cx, final int cy,
                                                  final int radius, final double angle_inc, final double[] angle,
                                                  final int[] index) {
        if (node.isCollapse()) {
            final double m = angle[0];
            final double r = circularRadiusFraction(node); // collapsed clade-root at its depth/distance radius
            node.setXcoord((float) (cx + (r * radius * Math.cos(m))));
            node.setYcoord((float) (cy + (r * radius * Math.sin(m))));
            _urt_nodeid_angle_map.put(node.getId(), m);
            angle[0] += angle_inc;
        }
        else if (node.isExternal()) {
            final double m = angle[0];
            final double r = circularRadiusFraction(node); // cladogram: outer ring (1); phylogram: distance-to-root
            node.setXcoord((float) (cx + (r * radius * Math.cos(m))));
            node.setYcoord((float) (cy + (r * radius * Math.sin(m))));
            _urt_nodeid_angle_map.put(node.getId(), m);
            _urt_nodeid_index_map.put(node.getId(), index[0]++);
            angle[0] += angle_inc;
        }
        else {
            for (int i = 0; i < node.getNumberOfDescendants(); ++i) {
                assignCircularDisplayedTipAngles(node.getChildNode(i), cx, cy, radius, angle_inc, angle, index);
            }
        }
    }

    final void paintCircularLite(final Phylogeny phy,
                                 final int center_x,
                                 final int center_y,
                                 final int radius,
                                 final Graphics2D g) {
        _root = phy.getRoot();
        _root.setXSecondary(center_x);
        _root.setYSecondary(center_y);
        // reuse the main paint's DISPLAYED-tip angles (collapse-aware; _urt_nodeid_angle_map was just populated by
        // paintCircular, which also folds in the starting angle) and skip tips hidden under a collapse -- so the
        // thumbnail's tip spacing MATCHES the main tree, instead of the old full-external-count spacing (which spread
        // every structural tip, including the ones hidden under a collapse, leaving a collapsed thumbnail out of scale).
        for (final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if (isHiddenUnderCollapse(n)) {
                continue;
            }
            final Double angle = _urt_nodeid_angle_map.get(n.getId());
            if (angle == null) {
                continue;
            }
            final double r = circularRadiusFraction(n); // outer ring (cladogram) or distance-to-root (phylogram)
            n.setXSecondary((float) (center_x + (r * radius * Math.cos(angle))));
            n.setYSecondary((float) (center_y + (r * radius * Math.sin(angle))));
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
                || shows(DisplayOption.WRITE_CONFIDENCE_VALUES))
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
        // Color the background
        if (!to_pdf) {
            final Rectangle r = getVisibleRect();
            g.setColor(getTreeColorSet().getBackgroundColor());
            if (!to_graphics_file) {
                g.fill(r);
            } else if (!_export_transparent_background) {
                if (getOptions().isExportBlackAndWhite()) {
                    g.setColor(Color.WHITE);
                }
                g.fillRect(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            }
            // else: transparent PNG export -- leave the (ARGB) canvas unfilled
            setupStroke(g);
        } else {
            g.setStroke(new BasicStroke(getOptions().getPdfLineWidth()));
        }
        if ((getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED)
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)) {
            _external_node_index = 0;
            // Position starting X of tree
            if (!_phylogeny.isRooted() /*|| ( _subtree_index > 0 )*/) {
                _phylogeny.getRoot().setXcoord(TreePanel.MOVE);
            } else if (displayedRootBranchLength() > 0.0) {
                // draw the root edge to scale (only for a full-tree phylogram root; a subtree draws a fixed stub via
                // the else branch -- see displayedRootBranchLength)
                double root_dtp = displayedRootBranchLength();
                if (breakLongBranchesActive() && (root_dtp > breakLongBranchCap())) {
                    root_dtp = breakLongBranchCap(); // cap a pathological root branch too (rare -- root usually has none)
                }
                _phylogeny.getRoot().setXcoord((float) (TreePanel.MOVE + (root_dtp * getXcorrectionFactor())));
            } else {
                _phylogeny.getRoot().setXcoord(TreePanel.MOVE + getXdistance());
            }
            // Position starting Y of tree (shifted down by the same header reserve used in calcParametersForPainting)
            _phylogeny.getRoot().setYcoord((getYdistance() * _phylogeny.getRoot().getNumberOfExternalNodes())
                    + (TreePanel.MOVE / 2.0f) + verticalBreadthPad() + annotationHeaderTopReserve());
            final int dynamic_hiding_factor = calcDynamicHidingFactor();
            if (shows(DisplayOption.DYNAMICALLY_HIDE_DATA)) {
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
                    /* || shows(DisplayOption.USE_STYLE) || getOptions().isShowDefaultNodeShapesForMarkedNodes()*/ //TODO check if this is really not needed.
                    || to_graphics_file || to_pdf;
            final boolean vertical = isVerticalOrientation();
            // the geologic axis is an alternative time-scale representation; suppress the numeric grid lines when it is
            // on (like the numeric scale bar + axis), so the two differently-spaced tick systems don't clash
            final boolean scale_grid_shown = getOptions().isShowScaleGrid() && getControlPanel().isDrawPhylogram()
                    && (getScaleDistance() > 0.0) && !geologicAxisApplies() && !calendarAxisApplies()
                    && !breakLongBranchesActive();
            if (!vertical && scale_grid_shown) {
                paintScaleGrid(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height);
            }
            if (!vertical) {
                // optional geologic grid lines at the fine-band boundaries (self-gated on the option + geologic axis)
                paintTimeAxisGridLines(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height);
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
                // optional geologic grid lines ride R the same way (logical lines -> horizontal grid lines by depth)
                paintTimeAxisGridLines(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height);
            }
            for (final PhylogenyNode element : _nodes_in_preorder) {
                paintNodeRectangular(g,
                        element,
                        to_pdf,
                        shows(DisplayOption.DYNAMICALLY_HIDE_DATA) && (dynamic_hiding_factor > 1),
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
                paintFossilRangeBars(g, to_pdf, to_graphics_file); // FAD/LAD stratigraphic-range bars on fossil tips
                paintAnnotationColumns(g); // tip-aligned columns (strip/heat map/bar/text), right of the labels
                paintMsaTrack(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height); // MSA cells + ruler
                paintCladeBands(g); // clade boxes/bars over the tree -- node coords set by the loop above
                paintAncestralPies(g, to_pdf, to_graphics_file); // per-node state pies, on top -- coords set above
                paintTipImages(g); // tip images at the branch end (the label was shifted right to make room)
            }
            else {
                // vertical parity: these overlays are drawn while g is rotated by R. Their geometry (zebra bands,
                // annotation cells, clade boxes/bars/brackets, HPD bars) is axis-aligned rects + lines, so it rides R
                // for free; the TEXT (annotation cells/headers, clade labels) is re-anchored upright inside the paint
                // methods. The scale grid rides R inside the R block above; the labeled scale axis is a side ruler
                // drawn as upright chrome after the base frame is restored (paintScaleAxisVertical, further below).
                paintZebraStripes(g, to_pdf, to_graphics_file, graphics_file_x, graphics_file_width); // faint row bands, behind
                paintHpdBars(g, to_pdf, to_graphics_file); // node-age HPD bars: a plain rect at each node -> rides R
                paintFossilRangeBars(g, to_pdf, to_graphics_file); // FAD/LAD tip range bars: axis-aligned rects -> ride R
                paintAnnotationColumnsVertical(g);
                paintCladeBands(g); // boxes ride R; bars/brackets draw the label upright (isVerticalOrientation branch)
                paintAncestralPies(g, to_pdf, to_graphics_file); // pies ride R: the disc stays a disc, wedges rotate
                paintTipImagesVertical(g); // tip images drawn UPRIGHT (not rotated under R) at the branch end
                if (geologicAxisApplies()) {
                    // the two-band geologic axis rides R into a side band down the breadth edge (bands + labels rotate)
                    paintGeologicTimeAxisVertical(g, to_pdf, to_graphics_file);
                }
            }
            paintHoverPreview(g, !(to_pdf || to_graphics_file)); // translucent select/deselect hover preview (rides R)
            paintFoundNodeHalos(g, to_pdf, to_graphics_file); // pulsing (screen) / static-glow (export) hit halos
            // restore the upright base frame before the viewport-fixed chrome (scale bar, tree name, overview, legends)
            if (vertical) {
                g.setTransform(orientation_saved);
            }
            // the geologic time axis takes over the bottom strip; suppress the numeric scale bar + axis when it is on
            final boolean geo_axis = geologicAxisApplies();
            final boolean calendar_axis = calendarAxisApplies();
            final boolean time_axis = geo_axis || calendar_axis; // a geologic/calendar time axis replaces the numeric one
            // a capped tree ("Break Long Branches") has no single linear distance scale across its whole width, so
            // suppress the full-width scale AXIS + GRID (the break glyph marks the discontinuity). The small scale BAR
            // stays -- it reads correctly for the un-broken (ingroup) part and is sized from the drawn extent (see
            // displayScaleDistance); use displayScaleDistance() for the >0 gate so a capped tree still shows it.
            final boolean break_active = breakLongBranchesActive();
            final boolean scale_shown = !time_axis && getOptions().isShowScale()
                    && getControlPanel().isDrawPhylogram() && (displayScaleDistance() > 0.0);
            final boolean axis_shown = !time_axis && !break_active && getOptions().isShowScaleAxis()
                    && getControlPanel().isDrawPhylogram() && (getScaleDistance() > 0.0);
            // the horizontal axis owns a reserved bottom band; lift the (viewport-fixed) scale bar clear above it (the
            // tree name is likewise raised, inside paintTreeName) so the three bottom overlays never overprint. Derive
            // both the lift AND the flag from the SAME layout reserve so they stay in lockstep with what is actually
            // drawn/reserved -- 0 in a vertical orientation (side ruler), for a cladogram / circular / unrooted, or an
            // unticked scale (never lift over a band that isn't there).
            final int bottom_reserve = scaleAxisBottomReserve();
            final boolean axis_shown_horizontal = bottom_reserve > 0;
            paintTimeAxisHint(g, to_pdf, to_graphics_file); // a dated CLADOGRAM: say why the time axis isn't showing
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
            if (geo_axis && !vertical) {
                // horizontal (root-left) geologic axis: chrome floating at the viewport bottom (the vertical variant
                // rode R inside the frame above, before the base transform was restored)
                paintGeologicTimeAxis(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height);
            }
            if (calendar_axis) {
                // the calendar (absolute-date) axis is a tick ruler like the numeric scale axis -- chrome (upright
                // labels), a bottom axis in a horizontal layout and a side ruler in a vertical one
                if (vertical) {
                    paintCalendarAxisVertical(g, to_pdf, to_graphics_file);
                } else {
                    paintCalendarAxis(g, to_pdf, to_graphics_file, graphics_file_y, graphics_file_height);
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
            if (!(to_graphics_file || to_pdf)) {
                paintTimeTreeBadge(g, getVisibleRect().x, getVisibleRect().width, getVisibleRect().y, to_pdf,
                        to_graphics_file);
            } else {
                paintTimeTreeBadge(g, graphics_file_x, graphics_file_width, graphics_file_y, to_pdf, to_graphics_file);
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
            if (shows(DisplayOption.DYNAMICALLY_HIDE_DATA)) {
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
                    to_graphics_file,
                    graphics_file_width,
                    graphics_file_height);
            paintRadialOverlays(g, to_pdf, to_graphics_file); // dots + pies + hover preview + halos (coords set above)
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
            // centre the ring in the drawing area (the viewport on screen, the export canvas on export -- both are the
            // tp's own width/height) so the circle is CENTRED and fills it, like the unrooted layout already does. The
            // old code placed the centre at min(w,h)/2, i.e. top-left-biased whenever the canvas is not square (the
            // radial preferred size is the rectangular tip-spread extent), which pushed the circle into a corner.
            // The ring is sized to the ZOOM diameter (WYSIWYG on screen AND export -- the domain-width cap + label
            // truncation are radialDiameter-based too, so all three agree) and CENTRED at the panel/canvas centre
            // getWidth()/getHeight(): on screen resetPreferredSize pads the panel to at least the viewport, so when the
            // tree is zoomed out BELOW fit the panel still fills the window and the ring stays CENTRED (was pinned to
            // the top-left corner, shrinking away). On export getWidth()/getHeight() is the export canvas (a
            // visible-only export TRANSLATES g to crop, so drawing at the full-panel centre yields the correct crop --
            // using graphics_file_width/height here would instead re-centre the whole ring into the cropped canvas).
            final int radius = circularRadius(radialDiameter());
            // centre in the drawing canvas: getWidth()/getHeight() on screen (and for a visible-only crop export),
            // the passed graphics-file size for a FULL export -- so a fixed-size export (canvas != panel) is centred,
            // not pushed off-canvas. For an ordinary full export the caller passes getWidth(), so the two agree.
            final int center_x = radialCanvasCenterX(graphics_file_width, to_pdf, to_graphics_file);
            final int center_y = radialCanvasCenterY(graphics_file_height, to_pdf, to_graphics_file);
            _circular_center_x = center_x; // so circularLabelAnchor / the ring hit-test match this exact geometry
            _circular_center_y = center_y;
            _circular_radius = radius;
            _dynamic_hiding_factor = 0;
            if (shows(DisplayOption.DYNAMICALLY_HIDE_DATA) && (radius > 0)) {
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
            // concentric ICS geologic bands behind the tree (a no-op unless GEOLOGIC + a circular phylogram); it is an
            // alternative radial time-scale representation, so it suppresses the numeric distance rings when on
            paintGeologicRingsCircular(g, center_x, center_y, radius > 0 ? radius : 0, to_pdf, to_graphics_file);
            // concentric CALENDAR-year rings behind the tree (a no-op unless CALENDAR + a circular phylogram)
            paintCalendarRingsCircular(g, center_x, center_y, radius > 0 ? radius : 0, to_pdf, to_graphics_file);
            if (!geologicRingsApplyCircular() && !calendarRingsApplyCircular() && !breakLongBranchesActiveCircular()) {
                // concentric distance rings behind the tree (a no-op unless this is a circular PHYLOGRAM); suppressed
                // while capping -- the radius is the CAPPED distance, so linear distance rings would misalign
                paintCircularScaleRings(g, center_x, center_y, radius > 0 ? radius : 0, to_pdf, to_graphics_file);
            }
            paintCircular(_phylogeny, getStartingAngle(), center_x, center_y, radius > 0 ? radius : 0, g, to_pdf,
                    to_graphics_file);
            // (aligned circular phylogram: the tip->ring leaders are drawn per-tip inside paintNodeDataUnrootedCirc,
            // right where the label is, so a leader appears iff its label does)
            // faint alternating angular wedges (row-tracking aid), over the tree but UNDER the rings/clade bands
            paintZebraStripesCircular(g, center_x, center_y, radius > 0 ? radius : 0);
            // node-age (HPD) bars as radial age-range segments (circular phylogram only), over the tree
            paintHpdBarsCircular(g, center_x, center_y, radius > 0 ? radius : 0, to_pdf, to_graphics_file);
            // fossil stratigraphic-range (FAD/LAD) bars as radial segments on the tips (circular phylogram only)
            paintFossilRangeBarsCircular(g, center_x, center_y, radius > 0 ? radius : 0, to_pdf, to_graphics_file);
            // tip-aligned annotation columns as concentric rings (strip/heat-map/bar/text), just past the tips + labels
            paintAnnotationColumnsCircular(g, center_x, center_y, radius > 0 ? radius : 0);
            // clade bands as polar sectors/arcs, over the tree (coords set above), like the rectangular wash
            paintCladeBandsCircular(g, center_x, center_y, radius > 0 ? radius : 0);
            // the geologic band names, up the top spoke ON TOP of the tree (the annuli are drawn behind, above)
            paintGeologicRingLabelsCircular(g, center_x, center_y, radius > 0 ? radius : 0, to_pdf, to_graphics_file);
            // optional Ma age labels at the coarse-band boundary radii, up the spoke (gated on "Geologic Boundary Ages")
            paintGeologicBoundaryAgesCircular(g, center_x, center_y, radius > 0 ? radius : 0, to_pdf, to_graphics_file);
            // protein-domain architectures riding each tip's spoke, just past the labels (iTOL circular-domains look)
            paintDomainsCircular(g, center_x, center_y, radius > 0 ? radius : 0, to_pdf, to_graphics_file);
            paintRadialOverlays(g, to_pdf, to_graphics_file); // dots + pies + hover preview + halos (coords set above)
            paintTimeAxisHint(g, to_pdf, to_graphics_file); // a dated circular CLADOGRAM: say why the ring axis isn't showing
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
                        x_pos + radius_ov,
                        y_pos + radius_ov,
                        // cap the (thumbnail-scaled) label reach like the main circle, so a long-label tree's overview
                        // thumbnail still draws a real circle instead of collapsing
                        (int) (radius_ov - Math.min(
                                getLongestExtNodeInfo() / (getVisibleRect().width / getOvRectangle().getWidth()),
                                radius_ov * RADIAL_LABEL_MAX_RATIO)),
                        g);
                g.setTransform(_at);
                paintOvRectangle(g);
            }
        }
        // The Color-by / Size-by tip DOTS render in every layout (paintRadialPropertyDots), so their legends are keyed
        // and shown everywhere. The annotation COLUMNS render as tip-aligned columns (rectangular) or concentric RINGS
        // (circular), so their legend is shown there too and suppressed ONLY in the UNROOTED layout, where the columns
        // aren't drawn (its bounds nulled so a stale hit region from a prior paint isn't clickable). The RANK legend
        // keys BRANCH colors, which render in every layout.
        final boolean draw_annotation_legend = annotationLegendVisible();
        final boolean draw_color_legend = isColorByProperty();
        if (draw_annotation_legend || draw_color_legend || hasRankLegend()) {
            final boolean to_screen = !(to_pdf || to_graphics_file);
            final Rectangle legend_bounds = to_screen
                    ? getVisibleRect()
                    : new Rectangle(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            // one legend slot; a header-focused annotation-column legend wins (the user just clicked it),
            // else the property-color legend, else the rank legend
            if (draw_annotation_legend) {
                drawAnnotationColumnLegend(g, legend_bounds, to_screen);
            } else if (draw_color_legend) {
                drawPropertyColorLegend(g, legend_bounds, to_screen);
            } else {
                drawRankLegend(g, legend_bounds, to_screen);
            }
        } else {
            _property_legend_bounds = null; // nothing in the shared slot -> no stale hit region
        }
        // "Size by" has its OWN legend (a separate, independently placed key), drawn last so it can appear
        // ALONGSIDE the color/rank legend -- the whole point of the combined color+size figure. Its size dots now
        // render in every layout, so the legend shows in every layout too.
        if (isSizeByProperty()) {
            final boolean to_screen = !(to_pdf || to_graphics_file);
            final Rectangle legend_bounds = to_screen ? getVisibleRect()
                    : new Rectangle(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            drawSizeLegend(g, legend_bounds, to_screen);
        }
        // ancestral-state pies have their OWN key (state -> color), drawn last so it can appear alongside the others.
        // Pies draw in EVERY layout (rectangular family + circular/unrooted), so the legend is never orphaned.
        if (isShowAncestralPies()) {
            final boolean to_screen = !(to_pdf || to_graphics_file);
            final Rectangle legend_bounds = to_screen ? getVisibleRect()
                    : new Rectangle(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            drawAncestralPieLegend(g, legend_bounds, to_screen,
                    (to_pdf || to_graphics_file) && getOptions().isExportBlackAndWhite());
        }
        // the internal-taxonomy key has its OWN text key (taxa by rank), drawn last (on top) so it can appear
        // alongside the others; a no-op unless the option is on AND the tree carries internal taxa.
        if (getOptions().isShowInternalTaxonomyKey()) {
            final boolean to_screen = !(to_pdf || to_graphics_file);
            final Rectangle legend_bounds = to_screen ? getVisibleRect()
                    : new Rectangle(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            drawInternalTaxonomyKey(g, legend_bounds, to_screen);
        }
        // the protein-domain legend (LEGEND label mode): its own draggable slot, drawn last so an overlap grab wins;
        // E-value-cutoff aware. Layout-agnostic (rides every layout + export), like the other legends.
        {
            final boolean to_screen = !(to_pdf || to_graphics_file);
            final Rectangle legend_bounds = to_screen ? getVisibleRect()
                    : new Rectangle(graphics_file_x, graphics_file_y, graphics_file_width, graphics_file_height);
            drawDomainLegend(g, legend_bounds, to_screen);
        }
        // reconcile the "Pulse Found Nodes" animation timer after EVERY screen paint (all layouts): starts it when a
        // hit halo was drawn (rectangular OR radial), stops it when none was (option off / no hit).
        if (!to_pdf && !to_graphics_file) {
            updatePulseTimer();
        }
    }

    final void recalculateMaxDistanceToRoot() {
        _max_distance_to_root = PhylogenyMethods.calculateMaxDistanceToRoot(getPhylogeny());
        // include the root's own branch length only when it is actually drawn to scale (the root-edge offset at
        // MOVE + rootBranch*xcorr) -- NOT for a subtree, whose inherited root branch is drawn as a fixed stub, so its
        // getMaxDistanceToRoot must exclude it (else the aligned tip-label / domain column shifts right and clips)
        if (!isCurrentTreeIsSubtree() && (getPhylogeny().getRoot().getDistanceToParent() > 0)) {
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
    /** Closes every open node window. Called before an undo/redo installs a different tree, since a node frame
     *  holds a direct reference to a node of the tree being replaced. */
    /** For tests: how many node windows this panel currently has open. */
    int openNodeFrameCountForTest() {
        return _node_frame_index;
    }

    private void closeAllNodeFrames() {
        // Disposed and cleared DIRECTLY rather than through NodeFrame.close(), which calls back into
        // removeEditNodeFrame and compacts the array underneath the loop (its stored _index is not updated by the
        // compaction, so a callback-driven loop can close a frame twice).
        for (int i = 0; i < _node_frames.length; ++i) {
            if (_node_frames[i] != null) {
                _node_frames[i].dispose();
                _node_frames[i] = null;
            }
        }
        _node_frame_index = 0;
    }

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
        if (isRadialLayout()) {
            // radial layouts use a SQUARE canvas of _radial_diameter (the single radial-zoom knob), decoupled from the
            // rectangular x/y-distance extent -- the circle is sized to it (see the CIRCULAR paint block), and the
            // unrooted spread is scaled by an urt-factor derived from it. PAD the panel to at least the viewport so a
            // tree zoomed out BELOW fit still fills the window (the scroll pane keeps it CENTRED) instead of shrinking
            // into the top-left corner; zoomed IN past the viewport the panel is larger -> scrollbars + panning.
            final int d = radialDiameter();
            int vw = d, vh = d;
            final java.awt.Container parent = getParent();
            if (parent instanceof javax.swing.JViewport) {
                final Dimension ext = ((javax.swing.JViewport) parent).getExtentSize();
                if (ext.width > 0) {
                    vw = ext.width;
                }
                if (ext.height > 0) {
                    vh = ext.height;
                }
            }
            setPreferredSize(new Dimension(Math.max(d, vw), Math.max(d, vh)));
            invalidateOrientationTransform();
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
        // every vertical tip label is anchored at labelTextStartX -- the tip's own effectiveNodeHalfBoxSize +
        // LABEL_GAP_AFTER_NODE_SHAPE PAST the node along the depth -- so that offset eats into the reserve; an UPRIGHT
        // (0-degree) label is drawn one extra TIP_LABEL_DEPTH_GAP past the tip (paintTipLabelHorizontal). Reserve both,
        // else after a fit the OUTERMOST tip's label pokes past the near depth edge (the top in root-bottom, the bottom
        // in root-top) and clips (was ~9px in root-bottom). Uses the MAX tip half-box (not the global default) so a
        // large custom node mark on a tip can't push its label past the edge.
        final int anchor_offset = maxTipEffectiveHalfBoxSize() + LABEL_GAP_AFTER_NODE_SHAPE;
        final int upright_gap = (effectiveTipLabelDirection() == Options.TIP_LABEL_DIRECTION.HORIZONTAL)
                ? TIP_LABEL_DEPTH_GAP : 0;
        // the TEXT label projected onto the depth (L_text*|sin| + lineH*|cos|) + the anchor offset + the upright gap +
        // the axis-aligned domain track past it. Uses the TEXT-ONLY longest (not getLongestExtNodeInfo, which already
        // folds in the domain width -- that would count the domain twice, over-compressing the depth when domains show).
        // + a little breathing room at the far depth edge (scaled up a touch for large/HiDPI fonts), so the outermost
        // tip label doesn't sit flush against the window edge (the 0.11.55 reserve only guaranteed no-clip = flush)
        final int edge_pad = Math.max(TIP_LABEL_DEPTH_EDGE_PAD, line_h / 3);
        return (int) Math.ceil((_length_of_longest_text_only * Math.abs(Math.sin(a))) + (line_h * Math.abs(Math.cos(a))))
                + anchor_offset + upright_gap + verticalDomainReserve() + edge_pad;
    }

    /** Extra depth (px) reserved past the tilted tip labels for the axis-aligned domain-architecture track that a
     *  vertical orientation draws as a per-tip vertical bar (0 unless domains are shown). The label reserve above tilts
     *  by sin(angle), but the domain boxes are axis-aligned (their FULL rendered width runs along the depth), so add
     *  the widest rendered track + the label gap + the render's internal 20px lead-in. */
    private int verticalDomainReserve() {
        if (!isVerticalOrientation() || !shows(DisplayOption.SHOW_DOMAIN_ARCHITECTURES)) {
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
                + msaRulerReserve() // the MSA column ruler's bottom band (root-left only; 0 otherwise)
                + msaConservationReserve() // ...and the conservation/consensus band just above it
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
                    + ForesterUtil.roundToInt((getXcorrectionFactor() * displayedTreeHeight()) + getXdistance());
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
        final Options.TREE_ORIENTATION o = getTreeOrientation();
        return ((o == Options.TREE_ORIENTATION.ROOT_TOP) || (o == Options.TREE_ORIENTATION.ROOT_BOTTOM))
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR)
                && (getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED);
    }

    /** (Re)builds the logical->screen rotation R (and its inverse) for the current vertical orientation, from the
     *  logical extents. Pure rotations (determinant +1, no mirror); the translate keeps the tree in the positive
     *  quadrant. ROOT_TOP turns the page 90deg clockwise: (x,y)->(H-y, x); ROOT_BOTTOM 90deg CCW: (x,y)->(y, W-x). */
    final void rebuildOrientationTransform() {
        final Options.TREE_ORIENTATION current = getTreeOrientation();
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

    /** Test hook: for the current vertical-orientation layout with UPRIGHT (HORIZONTAL) tip labels, the smallest
     *  signed margin (px) between an external tip's upright label FAR edge and the near canvas DEPTH edge -- the top
     *  edge in ROOT_BOTTOM, the bottom edge in ROOT_TOP. Negative means the outermost label is clipped off the
     *  canvas (the defect {@link #depthLabelReserve()} must prevent). Mirrors {@link #paintTipLabelHorizontal}'s
     *  anchor math exactly, so it measures the ACTUAL drawn label position, not a re-derivation of the reserve. */
    double minUprightLabelDepthMarginForTest() {
        final FontMetrics fm = getFontMetricsForLargeDefaultFont();
        final boolean root_bottom = getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_BOTTOM;
        double min = Double.MAX_VALUE;
        for (final PhylogenyNode t : _phylogeny.getExternalNodes()) {
            final Point2D.Double a = screenPoint(labelTextStartX(t), t.getYcoord());
            // far edge of the upright label along the depth (baseline +/- gap +/- the ascent/descent it spans)
            final double far = root_bottom ? (a.y - TIP_LABEL_DEPTH_GAP - fm.getAscent() - fm.getDescent())
                    : (a.y + TIP_LABEL_DEPTH_GAP + fm.getAscent() + fm.getDescent());
            min = Math.min(min, root_bottom ? far : (getHeight() - far));
        }
        return min;
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
        return (getTreeOrientation() == Options.TREE_ORIENTATION.ROOT_BOTTOM) ? -base : base;
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
    float alignedLabelColumnX() {
        return (float) ((displayedMaxDistanceToRoot() * getXcorrectionFactor()) + TreePanel.MOVE + getXdistance());
    }

    /** The logical X where a node's label text begins: the aligned column for an aligned tip, else just right of the
     *  node. Serves as the vertical-mode tilt pivot AND the aligned leader's far endpoint, so the tilted label sits
     *  on the end of its (vertical) leader instead of being pushed diagonally off it. Matches the label anchors in
     *  {@link #paintNodeData} / {@link #paintTaxonomy}. */
    private float labelTextStartX(final PhylogenyNode node) {
        final int half_box = effectiveNodeHalfBoxSize(node);
        // clustergram "labels below columns": tip/collapsed labels are drawn past the tip-aligned columns (aligned at
        // the far edge), so the dendrogram sits directly on the grid and the sample labels run along the bottom
        if (tipLabelsBelowColumns() && (node.isExternal() || node.isCollapse())) {
            return labelSegmentStartX(annotationColumnsEndX(), half_box, 0);
        }
        // an imaged tip's label is shifted right by the image slot so the image occupies the branch end (tipImageAdvance
        // is 0 for a non-imaged tip / aligned label, so nothing else moves)
        return labelSegmentStartX(isAlignedTipLabel(node) ? alignedLabelColumnX() : node.getXcoord(), half_box,
                tipImageAdvance(node));
    }

    /** Whether a branch-length value should be drawn for {@code node} -- the same gate {@link #paintNodeData} used
     *  inline, factored out so the vertical-orientation path (which draws it separately) shares one condition. */
    private boolean shouldWriteBranchLength(final PhylogenyNode node) {
        return shows(DisplayOption.WRITE_BRANCH_LENGTH_VALUES)
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
            return null; // no meaningful viewport (offscreen render / not-yet-realized) -> never cull
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

    /** THIS tab's orientation. Options.getTreeOrientation() is only the default new tabs are seeded with. */
    final Options.TREE_ORIENTATION getTreeOrientation() {
        return _tree_orientation;
    }

    /** Sets THIS tab's orientation. The caller re-fits: orientation swaps the layout's width and height, so a
     *  plain repaint would leave the old scroll extent. */
    final void setTreeOrientation(final Options.TREE_ORIENTATION orientation) {
        if (orientation != null) {
            _tree_orientation = orientation;
        }
    }

    /**
     * Every {@link DisplayOption}'s state FOR THIS TAB. The checkboxes in the control panel are shared by the whole
     * window and show whichever tab is in front; this is the tab's own copy, which is what a figure is made of and
     * what a saved figure has to restore. Continues the migration that already moved the orientation and the two
     * Display-Data toggles off the shared {@code Options} (see {@code PerTabViewStateTest}).
     * <p>
     * An option with no entry falls back to the shared widget, so a tab that has never been seeded behaves exactly
     * as it did before.
     */
    private final java.util.EnumMap<DisplayOption, Boolean> _display_state =
            new java.util.EnumMap<DisplayOption, Boolean>(DisplayOption.class);

    /** Whether {@code which} is on FOR THIS TAB. A tab with no opinion yet falls back to what the paint would
     *  have read before this state existed ({@link ControlPanel#currentValueOf}), so an unseeded tab behaves
     *  exactly as it always did -- including the readers that are not plain checkbox reads. */
    final boolean shows(final DisplayOption which) {
        final Boolean b = _display_state.get(which);
        if (b != null) {
            return b.booleanValue();
        }
        return (getControlPanel() != null) && getControlPanel().currentValueOf(which);
    }

    final void setShows(final DisplayOption which, final boolean on) {
        _display_state.put(which, Boolean.valueOf(on));
    }

    /** For tests: whether this tab has its own opinion on {@code which} yet. */
    final boolean hasOwnDisplayState(final DisplayOption which) {
        return _display_state.containsKey(which);
    }

    /** Whether THIS tab shows internal-node labels (the shared "Show Internal Data" checkbox reflects the
     *  current tab). */
    final boolean isShowInternalDataForThisTab() {
        return _show_internal_data;
    }

    final void setShowInternalDataForThisTab(final boolean show) {
        _show_internal_data = show;
    }

    final boolean isShowExternalDataForThisTab() {
        return _show_external_data;
    }

    final void setShowExternalDataForThisTab(final boolean show) {
        _show_external_data = show;
    }

    final void setPhylogenyGraphicsType(final PHYLOGENY_GRAPHICS_TYPE graphics_type) {
        _graphics_type = graphics_type;
        // a hover-preview target from the previous layout would otherwise draw a spurious select preview in the new
        // one (the pointer isn't over that branch there); it self-corrects on the next mouse move, but clear it now
        clearHoverPreview();
        if (isRadialLayout()) {
            invalidateRadialDiameter(); // switching TO a radial layout re-fits the square canvas to the viewport
        }
        setTextAntialias();
    }

    final void setStartingAngle(final double starting_angle) {
        _urt_starting_angle = starting_angle;
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
        if (shows(DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES) || shows(DisplayOption.SHOW_TAX_CODE)) {
            pri = DESCENDANT_SORT_PRIORITY.TAXONOMY;
        } else if (shows(DisplayOption.SHOW_SEQ_NAMES) || shows(DisplayOption.SHOW_SEQ_SYMBOLS)
                || shows(DisplayOption.SHOW_GENE_NAMES)) {
            pri = DESCENDANT_SORT_PRIORITY.SEQUENCE;
        }
        pushUndoCheckpoint("Ladderize Subtree");
        PhylogenyMethods.orderAppearanceX(node, true, pri);
        final String prov = TreePanelUtil.ladderizeProvenanceSentence(false, null, getPhylogeny().getName(),
                getPhylogeny().getNumberOfExternalNodes());
        final String existing = getPhylogeny().getDescription();
        getPhylogeny().setDescription(ForesterUtil.isEmpty(existing) ? prov : existing + " " + prov);
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
            double height = _phylogeny.calculateHeight(!_options.isCollapsedWithAverageHeigh());
            // keep the overview x-scale consistent with the capped main view (see calcParametersForPainting)
            if (breakLongBranchesActive() && (breakCappedHeight() > 0)) {
                height = breakCappedHeight();
            }
            // a subtree draws its root as a fixed stub, not to scale -> exclude it here too, so the overview mini-tree
            // fills its box at the same scale as the main view (else it renders too small)
            if (isCurrentTreeIsSubtree()) {
                height -= cappedRootBranchLength();
            }
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

    final private static void drawString(final String str, final float x, final float y, final Graphics2D g) {
        g.drawString(str, x, y);
    }

    final private static boolean plusPressed(final int key_code) {
        return ((key_code == KeyEvent.VK_ADD) || (key_code == KeyEvent.VK_PLUS)
                || (key_code == KeyEvent.VK_EQUALS) || (key_code == KeyEvent.VK_SEMICOLON)
                || (key_code == KeyEvent.VK_1));
    }

}
