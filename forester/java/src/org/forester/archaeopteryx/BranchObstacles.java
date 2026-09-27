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
 * The BRANCH LINES drawn in one radial (circular / unrooted) paint pass, kept so that a branch NUMBER can ask
 * whether its ink would cross a branch other than the one it labels.
 * <p>
 * <b>Why.</b> {@link LabelOccupancy} keeps numbers off each other; nothing kept them off the tree. In the unrooted
 * layout sibling branches leave one node a few degrees apart, and a 24 px number centred on a 16 px branch reaches
 * past both ends of it into that fan, so it lay across its sibling's branch or back over the parent's -- measured
 * on real trees at 3 of 16 drawn numbers (bat phylogeny), 4 of 25 (animal tree of life), 8 of 38 and 12 of 53 in
 * the circular layout, where a number on a short leg reaches the arc or the children's legs instead. A number
 * across a foreign branch is worse than a missing one: it reads as that branch's value. So a branch is an obstacle
 * exactly as another number is, under the same "Auto-hide Crowded Data" switch, and nothing is hidden unless a line
 * really runs through the ink -- the test is the number's ROTATED ink box against every segment, not a bound.
 * <p>
 * <b>What is registered.</b> Every leg the two radial painters draw and, in the circular layout, every arc (as a
 * polyline whose chords sag under half a pixel). A zero-length branch is drawn as nothing and registered as
 * nothing, so its number is never blocked by a line that is not there ("zero is a value, not an absence").
 * <p>
 * <b>Cost.</b> The same shape as {@link LabelOccupancy}: a uniform cell grid in primitive arrays, no allocation in
 * the steady state, one intrusive list node per (segment, cell). A query touches the cells under the number's
 * bounds and tests only the segments threaded there. A segment threaded through several of those cells is tested
 * once per cell: a per-query stamp to test it once instead was built and measured at nothing at all (a 2,339-clade
 * unrooted paint, 18.7 ms best with it and without), because the legs and chords registered here are short and a
 * number's bounds cover a cell or two, so it was removed rather than kept on the strength of a story.
 */
final class BranchObstacles {

    /** Floor on the cell size: cells far smaller than a number put every segment in many cells for nothing. */
    private final static float MIN_CELL  = 8.0f;
    /** A key no real cell can produce, so an untouched slot is recognisable without a separate occupancy array. */
    private final static long  EMPTY_KEY = Long.MIN_VALUE;
    /** Largest sagitta a chord of a registered arc may leave, in px. */
    private final static double ARC_SAG  = 0.5;
    /** The owner of a line that is nobody's own -- an arc. A number sits flush against its own LEG, which is why a
     *  leg is excluded for its owner; the fork arc at that leg's inner end is a line the number may reach past the
     *  node into, and a number lying across it is struck through whoever drew it. */
    final static long NO_OWNER = Long.MIN_VALUE;

    private long[]  _slot_key;  // cell key per open-addressing slot (EMPTY_KEY when free)
    private int[]   _slot_head; // index of the first list node in that cell, or -1
    private int[]   _node_seg;  // one node per (segment, cell) pair
    private int[]   _node_next;
    private int     _node_count;
    private float[] _x0;
    private float[] _y0;
    private float[] _x1;
    private float[] _y1;
    private long[]  _owner;     // the node whose branch the segment is part of
    private int     _seg_count;
    private int     _slot_used;
    private float   _cell = 32.0f;

    BranchObstacles() {
        _slot_key = new long[1 << 10];
        _slot_head = new int[1 << 10];
        java.util.Arrays.fill(_slot_key, EMPTY_KEY);
        _node_seg = new int[1024];
        _node_next = new int[1024];
        _x0 = new float[256];
        _y0 = new float[256];
        _x1 = new float[256];
        _y1 = new float[256];
        _owner = new long[256];
    }

    /** Starts a new pass. {@code cell_size} should be about the size of the numbers that will be tested. */
    void reset(final float cell_size) {
        java.util.Arrays.fill(_slot_key, EMPTY_KEY);
        _seg_count = 0;
        _node_count = 0;
        _slot_used = 0;
        _cell = Math.max(MIN_CELL, cell_size);
    }

    /**
     * Registers one drawn line of {@code owner}'s branch. A zero-length line is skipped: the painter draws nothing
     * for it, so there is nothing a number could cross.
     */
    void add(final long owner, final float x0, final float y0, final float x1, final float y1) {
        if ((x0 == x1) && (y0 == y1)) {
            return;
        }
        if (_seg_count == _x0.length) {
            final int n = _seg_count * 2;
            _x0 = java.util.Arrays.copyOf(_x0, n);
            _y0 = java.util.Arrays.copyOf(_y0, n);
            _x1 = java.util.Arrays.copyOf(_x1, n);
            _y1 = java.util.Arrays.copyOf(_y1, n);
            _owner = java.util.Arrays.copyOf(_owner, n);
        }
        final int seg = _seg_count++;
        _x0[seg] = x0;
        _y0[seg] = y0;
        _x1[seg] = x1;
        _y1[seg] = y1;
        _owner[seg] = owner;
        final int col0 = (int) Math.floor(Math.min(x0, x1) / _cell);
        final int col1 = (int) Math.floor(Math.max(x0, x1) / _cell);
        final int row0 = (int) Math.floor(Math.min(y0, y1) / _cell);
        final int row1 = (int) Math.floor(Math.max(y0, y1) / _cell);
        for (int c = col0; c <= col1; ++c) {
            for (int r = row0; r <= row1; ++r) {
                link(key(c, r), seg);
            }
        }
    }

    /**
     * Registers a fork arc -- centre {@code (cx, cy)}, radius {@code r}, from angle {@code a0} to {@code a1} in
     * radians in the painter's own convention (a point at angle {@code a} sits at {@code (cx + r cos a, cy + r sin
     * a)}) -- as chords that never sag more than half a pixel from it. Under {@link #NO_OWNER}: no number is excused
     * from an arc, the branch that drew it included.
     */
    void addArc(final double cx, final double cy, final double r, final double a0, final double a1) {
        if (r <= 0) {
            return;
        }
        final double sweep = a1 - a0;
        // sagitta of a chord over angle t is r (1 - cos(t/2)); keep it under ARC_SAG
        final double max_piece = (r > ARC_SAG) ? (2.0 * Math.acos(1.0 - (ARC_SAG / r))) : Math.PI;
        final int pieces = Math.max(1, (int) Math.ceil(Math.abs(sweep) / max_piece));
        float px = (float) (cx + (r * Math.cos(a0)));
        float py = (float) (cy + (r * Math.sin(a0)));
        for (int i = 1; i <= pieces; ++i) {
            final double a = a0 + ((sweep * i) / pieces);
            final float nx = (float) (cx + (r * Math.cos(a)));
            final float ny = (float) (cy + (r * Math.sin(a)));
            add(NO_OWNER, px, py, nx, ny);
            px = nx;
            py = ny;
        }
    }

    /**
     * Whether any registered line of a branch OTHER than {@code owner}'s runs through the box that is {@code half_w}
     * either side of {@code (mid_x, mid_y)} along the direction {@code m} and from {@code top} to {@code bottom}
     * across it (local y, positive to the right of the direction of travel on a y-down canvas), grown by
     * {@code margin} on every side. A line touching the grown box counts: {@code margin} is what turns "the ink
     * touches the line" into a crossing, so callers pass half the line's width plus its antialiasing.
     */
    boolean crosses(final long owner, final double mid_x, final double mid_y, final double m, final double half_w,
                    final double top, final double bottom, final double margin) {
        final double hw = half_w + margin;
        final double t = top - margin;
        final double b = bottom + margin;
        final double cos = Math.cos(m), sin = Math.sin(m);
        // the box's device bounds, to know which cells to look in
        final double rot_w = (Math.abs(2 * hw * cos) + Math.abs((b - t) * sin)) / 2.0;
        final double rot_h = (Math.abs(2 * hw * sin) + Math.abs((b - t) * cos)) / 2.0;
        final double cy = (t + b) / 2.0; // the box's centre in the local frame is off the line
        final double cx_dev = mid_x - (cy * sin);
        final double cy_dev = mid_y + (cy * cos);
        final int col0 = (int) Math.floor((cx_dev - rot_w) / _cell);
        final int col1 = (int) Math.floor((cx_dev + rot_w) / _cell);
        final int row0 = (int) Math.floor((cy_dev - rot_h) / _cell);
        final int row1 = (int) Math.floor((cy_dev + rot_h) / _cell);
        for (int c = col0; c <= col1; ++c) {
            for (int r = row0; r <= row1; ++r) {
                for (int n = head(key(c, r)); n >= 0; n = _node_next[n]) {
                    final int seg = _node_seg[n];
                    if (_owner[seg] == owner) {
                        continue;
                    }
                    // the segment in the box's own frame, then a clip against the axis-aligned box there
                    final double dx0 = _x0[seg] - mid_x, dy0 = _y0[seg] - mid_y;
                    final double dx1 = _x1[seg] - mid_x, dy1 = _y1[seg] - mid_y;
                    final double u0 = (dx0 * cos) + (dy0 * sin), v0 = (dy0 * cos) - (dx0 * sin);
                    final double u1 = (dx1 * cos) + (dy1 * sin), v1 = (dy1 * cos) - (dx1 * sin);
                    if (segmentMeetsBox(u0, v0, u1, v1, -hw, t, hw, b)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /** Liang-Barsky: does the segment (a, b) meet the box [x0, x1] x [y0, y1], touching included? */
    static boolean segmentMeetsBox(final double ax, final double ay, final double bx, final double by,
                                   final double x0, final double y0, final double x1, final double y1) {
        final double dx = bx - ax, dy = by - ay;
        double t_in = 0, t_out = 1;
        for (int side = 0; side < 4; ++side) {
            final double p, q;
            switch (side) {
                case 0:  p = -dx; q = ax - x0; break;
                case 1:  p = dx;  q = x1 - ax; break;
                case 2:  p = -dy; q = ay - y0; break;
                default: p = dy;  q = y1 - ay; break;
            }
            if (p == 0) {
                if (q < 0) {
                    return false; // parallel to this side and wholly outside it
                }
            } else {
                final double r = q / p;
                if (p < 0) {
                    t_in = Math.max(t_in, r);
                } else {
                    t_out = Math.min(t_out, r);
                }
            }
        }
        return t_in <= t_out;
    }

    /** For tests: how many line pieces are registered this pass (an arc counts once per chord). */
    int segmentCountForTest() {
        return _seg_count;
    }

    private int head(final long cell_key) {
        final int slot = findSlot(cell_key);
        return (_slot_key[slot] == EMPTY_KEY) ? -1 : _slot_head[slot];
    }

    private void link(final long cell_key, final int seg) {
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
        if (_node_count == _node_seg.length) {
            _node_seg = java.util.Arrays.copyOf(_node_seg, _node_count * 2);
            _node_next = java.util.Arrays.copyOf(_node_next, _node_count * 2);
        }
        _node_seg[_node_count] = seg;
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
