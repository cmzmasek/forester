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
 * <b>Cost.</b> A uniform cell grid in primitive arrays, one intrusive list node per (box, cell), keyed by each box's
 * axis-aligned bounds; the oriented test runs only against the boxes threaded through the cells a query touches.
 */
final class OrientedOccupancy {

    private final static float MIN_CELL  = 8.0f;
    private final static long  EMPTY_KEY = Long.MIN_VALUE;

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
        final double ex = Math.abs(hw * cos) + Math.abs(hh * sin); // half the axis-aligned bounds
        final double ey = Math.abs(hw * sin) + Math.abs(hh * cos);
        final int col0 = (int) Math.floor((cx - ex) / _cell);
        final int col1 = (int) Math.floor((cx + ex) / _cell);
        final int row0 = (int) Math.floor((cy - ey) / _cell);
        final int row1 = (int) Math.floor((cy + ey) / _cell);
        for (int c = col0; c <= col1; ++c) {
            for (int r = row0; r <= row1; ++r) {
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
        final double ey = Math.abs(hw * _sin[box]) + Math.abs(hh * _cos[box]);
        final int col0 = (int) Math.floor((cx - ex) / _cell);
        final int col1 = (int) Math.floor((cx + ex) / _cell);
        final int row0 = (int) Math.floor((cy - ey) / _cell);
        final int row1 = (int) Math.floor((cy + ey) / _cell);
        for (int c = col0; c <= col1; ++c) {
            for (int r = row0; r <= row1; ++r) {
                link(key(c, r), box);
            }
        }
    }

    /** For tests: how many marks have been placed this pass. */
    int claimedCountForTest() {
        return _box_count;
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
