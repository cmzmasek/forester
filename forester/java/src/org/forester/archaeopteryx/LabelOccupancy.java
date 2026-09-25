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
 * "What has already been drawn here" for ONE paint pass -- the mechanism behind auto-hiding crowded branch data.
 * A mark {@link #claim(float, float, float, float)}s the rectangle it is about to occupy; it is drawn only if
 * nothing already drawn overlaps it.
 * <p>
 * <b>Why not a branch-length threshold.</b> The obvious rule -- drop a support value when its branch is drawn too
 * short to carry the text -- was built and measured first, and it is a poor proxy for what actually makes a number
 * unreadable. A number becomes unreadable when it overlaps ANOTHER number, and that depends on how close the two
 * are in <i>y</i> as much as in <i>x</i>; branch length sees only <i>x</i>. On a real 825-tip tree at a comfortable
 * zoom, a "branch under 1/3 of the text" rule hid 285 numbers of which only 51 actually collided -- 234 that read
 * perfectly well -- while still leaving 22 pile-ups on screen. No threshold value did better: every setting both
 * over-hid and under-fixed. Asking the question directly costs about the same and is exact: nothing is hidden
 * unless something really is in the way, so a lone zero-length branch keeps its number ("zero is a value, not an
 * absence") with no special case for it.
 * <p>
 * <b>Cost.</b> A uniform grid, cell-sized to the marks themselves, so a claim touches a handful of cells and the
 * whole pass is linear in the number of marks rather than quadratic.
 * <p>
 * <b>Determinism.</b> First claim wins, and the paint walks the tree in preorder, so the mark nearer the root
 * keeps its place and the same tree at the same size always drops the same marks -- on screen and in every export.
 */
final class LabelOccupancy {

    /** Grid cell size floor: a cell far smaller than the marks would put every mark in many cells for nothing. */
    private final static float MIN_CELL = 4.0f;

    /**
     * The grid, held in primitive arrays rather than a {@code Map<Long, List<float[]>>}.
     * <p>
     * This runs in the paint hot path: a 30x12 number box over a ~14 px cell touches 4-6 cells, and the map form
     * boxed a {@link Long} for every one of them TWICE (probe, then insert) plus an {@code ArrayList} per cell --
     * on an 825-tip tree carrying confidences and branch lengths that is tens of thousands of short-lived objects
     * per frame, in exactly the paint the FPS readout was added to measure. Here a cell is an open-addressed slot
     * holding the head of an intrusive list threaded through {@code _next}, and the boxes live in four flat float
     * arrays. No allocation at all in the steady state: {@link #reset} keeps the arrays and only clears the slots.
     */
    private long[]  _slot_key;   // cell key per open-addressing slot (EMPTY_KEY when free)
    private int[]   _slot_head;  // index of the first LIST NODE in that cell, or -1
    // One node per (box, cell) pair -- NOT per box. A box wider than a cell belongs to several cells at once, so
    // threading the list through a per-box "next" cannot represent it: each extra cell overwrote the previous
    // cell's link and silently truncated its chain, so a second box in a cell stopped being seen. Caught by
    // running this grid and the Map<Long,List<float[]>> it replaced over the same 20k claims: 12000 placed
    // against 12070. Nodes are cheap (two ints) and reused across passes.
    private int[]   _node_box;
    private int[]   _node_next;
    private int     _node_count;
    private float[] _bx;
    private float[] _by;
    private float[] _bw;
    private float[] _bh;
    private int     _box_count;
    private int     _slot_used;
    private float   _cell = 16.0f;

    /** A key no real cell can produce, so an untouched slot is recognisable without a separate occupancy array. */
    private final static long EMPTY_KEY = Long.MIN_VALUE;

    LabelOccupancy() {
        allocate(1 << 10, 256);
    }

    private void allocate(final int slots, final int boxes) {
        _slot_key = new long[slots];
        _slot_head = new int[slots];
        java.util.Arrays.fill(_slot_key, EMPTY_KEY);
        _node_box = new int[boxes * 4];
        _node_next = new int[boxes * 4];
        _bx = new float[boxes];
        _by = new float[boxes];
        _bw = new float[boxes];
        _bh = new float[boxes];
    }

    /** Starts a new pass. {@code cell_size} should be about the size of the marks being placed. */
    void reset(final float cell_size) {
        java.util.Arrays.fill(_slot_key, EMPTY_KEY);
        _box_count = 0;
        _node_count = 0;
        _slot_used = 0;
        _cell = Math.max(MIN_CELL, cell_size);
    }

    /**
     * Claims the rectangle for a mark about to be drawn.
     *
     * @return true when nothing already claimed overlaps it (and it is now recorded), false when the space is
     *         taken -- in which case nothing is recorded, so a rejected mark never blocks a later one
     */
    boolean claim(final float x, final float y, final float w, final float h) {
        if ((w <= 0) || (h <= 0)) {
            // Not a behavioural guard: under the strict comparisons below a zero-size box can never overlap
            // anything, so it would be granted anyway. What this avoids is RECORDING marks that reserve no
            // space -- grid entries that can only ever cost memory and comparisons. (Measured by
            // LabelOccupancyTest, which asserts the map stays empty.)
            return true;
        }
        final int col0 = (int) Math.floor(x / _cell);
        final int col1 = (int) Math.floor((x + w) / _cell);
        final int row0 = (int) Math.floor(y / _cell);
        final int row1 = (int) Math.floor((y + h) / _cell);
        for (int c = col0; c <= col1; ++c) {
            for (int r = row0; r <= row1; ++r) {
                for (int n = head(key(c, r)); n >= 0; n = _node_next[n]) {
                    final int i = _node_box[n];
                    if ((x < (_bx[i] + _bw[i])) && ((x + w) > _bx[i]) && (y < (_by[i] + _bh[i]))
                            && ((y + h) > _by[i])) {
                        return false;
                    }
                }
            }
        }
        final int box = addBox(x, y, w, h);
        for (int c = col0; c <= col1; ++c) {
            for (int r = row0; r <= row1; ++r) {
                link(key(c, r), box);
            }
        }
        return true;
    }

    /** The first box recorded in {@code cell_key}, or -1. */
    private int head(final long cell_key) {
        final int slot = findSlot(cell_key);
        return (_slot_key[slot] == EMPTY_KEY) ? -1 : _slot_head[slot];
    }

    /** Threads a new node for {@code box} onto {@code cell_key}'s list, creating the slot if this cell is new. */
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

    private int addBox(final float x, final float y, final float w, final float h) {
        if (_box_count == _bx.length) {
            final int n = _box_count * 2;
            _bx = java.util.Arrays.copyOf(_bx, n);
            _by = java.util.Arrays.copyOf(_by, n);
            _bw = java.util.Arrays.copyOf(_bw, n);
            _bh = java.util.Arrays.copyOf(_bh, n);
        }
        _bx[_box_count] = x;
        _by[_box_count] = y;
        _bw[_box_count] = w;
        _bh[_box_count] = h;
        return _box_count++;
    }

    /** For tests: how many marks have been placed this pass (a box spanning several cells still counts once). */
    int claimedCountForTest() {
        return _box_count;
    }

    private static long key(final int col, final int row) {
        return (((long) col) << 32) ^ (row & 0xffffffffL);
    }
}
