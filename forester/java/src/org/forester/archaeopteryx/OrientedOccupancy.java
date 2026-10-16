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

/**
 * "What has already been drawn here" for marks that are NOT axis-aligned: the tip labels of the radial layouts,
 * which ride their spokes at every angle. The same job as {@link LabelOccupancy}, which keeps the rectangular
 * layout's numbers off each other with axis-aligned boxes -- but a rotated label's axis-aligned bounds are up to
 * twice its area (a 100 x 14 px label at 45 degrees bounds 80 x 80), and reserving those would hide neighbours that
 * do not touch. So a box here is ORIENTED -- centre, half extents, angle -- and two boxes overlap when no edge of
 * either separates them (the separating-axis test, four axes).
 * <p>
 * <b>Why.</b> The unrooted layout never thinned its tip labels at all, and circular hid every k-th by index, a
 * proxy that hides labels which do not overlap and keeps ones that do. Measured on the bat phylogeny at 1100 x 850:
 * 66 overlapping pairs among 34 horizontal labels unrooted, 21 circular; with labels along the spoke 40 and 0. A
 * label is drawn only if nothing already drawn is in its way, first come, and a found node's label always -- the
 * rule the numbers follow, and the one the manual states.
 * <p>
 * <b>Cost.</b> A uniform cell grid in primitive arrays, one intrusive list node per (box, cell); the oriented test
 * runs only against the boxes threaded through the cells a query touches. A box is threaded through the cells it
 * COVERS -- column by column, the rows between its lowest and highest point in that column -- not through every
 * cell of its axis-aligned bounds: a long label along a 45-degree spoke covers about two cells per column, where
 * its bounds hold a square of them (a 10,000 x 14 px box at 16 px cells: ~1,300 entries against 196,000; until
 * 0.11.181 the bounds were used, see {@link BranchObstacles} for the paint that ran out of heap). A pass that
 * needed more than {@link #SHRINK_ABOVE} slots is abnormal; what it took is given back at the next {@link #reset}.
 */
final class OrientedOccupancy {

    private final static float MIN_CELL  = 8.0f;
    private final static long  EMPTY_KEY = Long.MIN_VALUE;
    /** A slot table larger than this after a pass is an abnormal pass's; {@link #reset} gives it back. */
    final static int           SHRINK_ABOVE = 1 << 18;

    private long[]  _slot_key;
    private int[]   _slot_head;
    private int[]   _node_box;
    private int[]   _node_next;
    private int     _node_count;
    private float[] _cx;
    private float[] _cy;
    private float[] _hw;
    private float[] _hh;
    private float[] _cos;
    private float[] _sin;
    private int     _box_count;
    private int     _slot_used;
    private float   _cell = 32.0f;
    // scratch for the box being threaded or looked up: its four corners, and the rows it reaches in one column
    private final double[] _kx = new double[4];
    private final double[] _ky = new double[4];
    private int     _row_lo;
    private int     _row_hi;

    OrientedOccupancy() {
        _slot_key = new long[1 << 10];
        _slot_head = new int[1 << 10];
        java.util.Arrays.fill(_slot_key, EMPTY_KEY);
        _node_box = new int[1024];
        _node_next = new int[1024];
        _cx = new float[256];
        _cy = new float[256];
        _hw = new float[256];
        _hh = new float[256];
        _cos = new float[256];
        _sin = new float[256];
    }

    /** Starts a new pass. {@code cell_size} should be about the height of the marks being placed. */
    void reset(final float cell_size) {
        if (_slot_key.length > SHRINK_ABOVE) {
            // an abnormal pass grew the tables; start the next one from the small ones rather than keep that memory
            _slot_key = new long[1 << 10];
            _slot_head = new int[1 << 10];
            _node_box = new int[1024];
            _node_next = new int[1024];
            _cx = new float[256];
            _cy = new float[256];
            _hw = new float[256];
            _hh = new float[256];
            _cos = new float[256];
            _sin = new float[256];
        }
        java.util.Arrays.fill(_slot_key, EMPTY_KEY);
        _box_count = 0;
        _node_count = 0;
        _slot_used = 0;
        _cell = Math.max(MIN_CELL, cell_size);
    }

    /**
     * Whether a box -- centre {@code (cx, cy)}, half extents {@code hw} along its own direction {@code theta} and
     * {@code hh} across it -- overlaps any box already claimed this pass. A read: nothing is recorded.
     */
    boolean overlaps(final double cx, final double cy, final double hw, final double hh, final double theta) {
        return find(cx, cy, hw, hh, theta) >= 0;
    }

    /** Records a box: for a mark the caller has decided to draw -- after {@link #overlaps} said the way was clear,
     *  or regardless of it (a found node's label) -- so that what comes after keeps clear of it. Asking and placing
     *  are separate calls so that a caller can ask for SEVERAL boxes (a name and its image) before committing any;
     *  a single ask-and-record {@code claim} existed and nothing in the program used it (a review find,
     *  2026-09-27). A refused mark is simply never placed, so it never blocks a later one. */
    void place(final double cx, final double cy, final double hw, final double hh, final double theta) {
        if ((hw <= 0) || (hh <= 0)) {
            return;
        }
        record(cx, cy, hw, hh, theta);
    }

    /** The index of a recorded box overlapping the query, or -1. */
    private int find(final double cx, final double cy, final double hw, final double hh, final double theta) {
        final double cos = Math.cos(theta), sin = Math.sin(theta);
        final double ex = Math.abs(hw * cos) + Math.abs(hh * sin); // half the axis-aligned bounds, along x
        final int col0 = (int) Math.floor((cx - ex) / _cell);
        final int col1 = (int) Math.floor((cx + ex) / _cell);
        corners(cx, cy, hw, hh, cos, sin);
        for (int c = col0; c <= col1; ++c) {
            if (!stripRows(Math.max(cx - ex, c * (double) _cell), Math.min(cx + ex, (c + 1) * (double) _cell))) {
                continue;
            }
            for (int r = _row_lo; r <= _row_hi; ++r) {
                for (int n = head(key(c, r)); n >= 0; n = _node_next[n]) {
                    final int i = _node_box[n];
                    if (overlap(cx, cy, hw, hh, cos, sin, _cx[i], _cy[i], _hw[i], _hh[i], _cos[i], _sin[i])) {
                        return i;
                    }
                }
            }
        }
        return -1;
    }

    /**
     * Separating-axis test for two oriented rectangles: they overlap unless one of the four edge directions
     * separates their projections. Touching edges do not overlap (strict), as in {@link LabelOccupancy}.
     */
    static boolean overlap(final double ax, final double ay, final double ahw, final double ahh, final double acos,
                           final double asin, final double bx, final double by, final double bhw, final double bhh,
                           final double bcos, final double bsin) {
        final double dx = bx - ax, dy = by - ay;
        // the four candidate axes: a's direction, a's normal, b's direction, b's normal
        return !separated(dx, dy, acos, asin, ahw, ahh, acos, asin, bhw, bhh, bcos, bsin)
                && !separated(dx, dy, -asin, acos, ahw, ahh, acos, asin, bhw, bhh, bcos, bsin)
                && !separated(dx, dy, bcos, bsin, ahw, ahh, acos, asin, bhw, bhh, bcos, bsin)
                && !separated(dx, dy, -bsin, bcos, ahw, ahh, acos, asin, bhw, bhh, bcos, bsin);
    }

    /** Whether the axis {@code (ux, uy)} separates the two boxes: the distance between their centres along it
     *  exceeds the sum of their half-widths along it. */
    private static boolean separated(final double dx, final double dy, final double ux, final double uy,
                                     final double ahw, final double ahh, final double acos, final double asin,
                                     final double bhw, final double bhh, final double bcos, final double bsin) {
        final double ra = (ahw * Math.abs((acos * ux) + (asin * uy))) + (ahh * Math.abs((-asin * ux) + (acos * uy)));
        final double rb = (bhw * Math.abs((bcos * ux) + (bsin * uy))) + (bhh * Math.abs((-bsin * ux) + (bcos * uy)));
        return Math.abs((dx * ux) + (dy * uy)) >= (ra + rb);
    }

    private void record(final double cx, final double cy, final double hw, final double hh, final double theta) {
        if (_box_count == _cx.length) {
            final int n = _box_count * 2;
            _cx = java.util.Arrays.copyOf(_cx, n);
            _cy = java.util.Arrays.copyOf(_cy, n);
            _hw = java.util.Arrays.copyOf(_hw, n);
            _hh = java.util.Arrays.copyOf(_hh, n);
            _cos = java.util.Arrays.copyOf(_cos, n);
            _sin = java.util.Arrays.copyOf(_sin, n);
        }
        final int box = _box_count++;
        _cx[box] = (float) cx;
        _cy[box] = (float) cy;
        _hw[box] = (float) hw;
        _hh[box] = (float) hh;
        _cos[box] = (float) Math.cos(theta);
        _sin[box] = (float) Math.sin(theta);
        final double ex = Math.abs(hw * _cos[box]) + Math.abs(hh * _sin[box]);
        final int col0 = (int) Math.floor((cx - ex) / _cell);
        final int col1 = (int) Math.floor((cx + ex) / _cell);
        corners(cx, cy, hw, hh, _cos[box], _sin[box]);
        for (int c = col0; c <= col1; ++c) {
            if (!stripRows(Math.max(cx - ex, c * (double) _cell), Math.min(cx + ex, (c + 1) * (double) _cell))) {
                continue;
            }
            for (int r = _row_lo; r <= _row_hi; ++r) {
                link(key(c, r), box);
            }
        }
    }

    /** The four corners of a box into the scratch arrays: centre +/- hw along (cos, sin), +/- hh along the normal. */
    private void corners(final double cx, final double cy, final double hw, final double hh, final double cos,
                         final double sin) {
        final double ux = hw * cos, uy = hw * sin;   // half the box along its direction
        final double vx = -hh * sin, vy = hh * cos;  // half the box across it
        _kx[0] = cx + ux + vx;
        _ky[0] = cy + uy + vy;
        _kx[1] = cx - ux + vx;
        _ky[1] = cy - uy + vy;
        _kx[2] = cx - ux - vx;
        _ky[2] = cy - uy - vy;
        _kx[3] = cx + ux - vx;
        _ky[3] = cy + uy - vy;
    }

    /**
     * The rows the scratch box reaches inside the column strip {@code x in [xa, xb]}, into {@code _row_lo..hi}:
     * the y-range of its four edges clipped to the strip, which for a convex box is the y-range of the box itself
     * there (its extreme points in a strip lie on its boundary). False when no edge reaches the strip (a strip
     * inside the bounds always has one). No slack against rounding, for the reason given at
     * {@code BranchObstacles.linkAlong}: an extreme is computed inexactly only on a column edge, where it is the
     * single point of the box in the cell it could be filed out of.
     */
    private boolean stripRows(final double xa, final double xb) {
        double lo = Double.POSITIVE_INFINITY, hi = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < 4; ++i) {
            final int j = (i + 1) & 3;
            final double ax = _kx[i], ay = _ky[i], bx = _kx[j], by = _ky[j];
            final double dx = bx - ax;
            final double ya, yb;
            if (dx == 0) {
                if ((ax < xa) || (ax > xb)) {
                    continue; // a vertical edge outside the strip
                }
                ya = ay;
                yb = by;
            } else {
                final double ta = (xa - ax) / dx, tb = (xb - ax) / dx;
                final double t0 = Math.max(0, Math.min(ta, tb)), t1 = Math.min(1, Math.max(ta, tb));
                if (t0 > t1) {
                    continue; // the edge never enters the strip
                }
                ya = ay + (t0 * (by - ay));
                yb = ay + (t1 * (by - ay));
            }
            lo = Math.min(lo, Math.min(ya, yb));
            hi = Math.max(hi, Math.max(ya, yb));
        }
        if (lo > hi) {
            return false;
        }
        _row_lo = (int) Math.floor(lo / _cell);
        _row_hi = (int) Math.floor(hi / _cell);
        return true;
    }

    /** For tests: how many marks have been placed this pass. */
    int claimedCountForTest() {
        return _box_count;
    }

    /** For tests: how many (box, cell) entries this pass holds -- the cost of threading. */
    int cellEntryCountForTest() {
        return _node_count;
    }

    /** For tests: the slot table's size -- what a pass grew it to, or what a reset gave back. */
    int slotTableSizeForTest() {
        return _slot_key.length;
    }

    /** For the error log's out-of-memory entry: the sizes of what this grid holds. */
    String describeForLog() {
        return "boxes " + _box_count + ", cell entries " + _node_count + " of " + _node_box.length + ", cells "
                + _slot_used + " in a table of " + _slot_key.length;
    }

    private int head(final long cell_key) {
        final int slot = findSlot(cell_key);
        return (_slot_key[slot] == EMPTY_KEY) ? -1 : _slot_head[slot];
    }

    private void link(final long cell_key, final int box) {
        int slot = findSlot(cell_key);
        if (_slot_key[slot] == EMPTY_KEY) {
            if (((_slot_used + 1) * 2) >= _slot_key.length) {
                growSlots();
                slot = findSlot(cell_key);
            }
            _slot_key[slot] = cell_key;
            _slot_head[slot] = -1;
            ++_slot_used;
        }
        if (_node_count == _node_box.length) {
            _node_box = java.util.Arrays.copyOf(_node_box, _node_count * 2);
            _node_next = java.util.Arrays.copyOf(_node_next, _node_count * 2);
        }
        _node_box[_node_count] = box;
        _node_next[_node_count] = _slot_head[slot];
        _slot_head[slot] = _node_count++;
    }

    /** Linear probing; the table is kept under half full, so a free slot always exists. */
    private int findSlot(final long cell_key) {
        final int mask = _slot_key.length - 1;
        int i = ((int) (cell_key ^ (cell_key >>> 32)) * 0x9E3779B1) & mask;
        while ((_slot_key[i] != EMPTY_KEY) && (_slot_key[i] != cell_key)) {
            i = (i + 1) & mask;
        }
        return i;
    }

    private void growSlots() {
        final long[] old_keys = _slot_key;
        final int[] old_heads = _slot_head;
        _slot_key = new long[old_keys.length * 2];
        _slot_head = new int[old_keys.length * 2];
        java.util.Arrays.fill(_slot_key, EMPTY_KEY);
        for (int i = 0; i < old_keys.length; ++i) {
            if (old_keys[i] != EMPTY_KEY) {
                final int slot = findSlot(old_keys[i]);
                _slot_key[slot] = old_keys[i];
                _slot_head[slot] = old_heads[i];
            }
        }
    }

    private static long key(final int col, final int row) {
        return (((long) col) << 32) ^ (row & 0xffffffffL);
    }
}
