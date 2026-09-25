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

import java.util.EnumMap;
import java.util.Map;

/**
 * Divides the width of a rectangular tree view between the TREE itself and everything drawn beside it -- tip labels,
 * the domain track, the sequence alignment, the annotation columns, the legend column and the clade bands.
 * <p>
 * <b>Why this exists.</b> Each of those used to reserve its width on its own, capped against the WHOLE panel width and
 * blind to the others: labels took up to 95% (70% before the font auto-fit gave up), the alignment 60% of the viewport,
 * the legend 40%. Those caps sum to well past 100%, and the tree was never allocated anything at all -- it simply got
 * whatever was left over. On a deep tree carrying domains and an alignment the leftover went NEGATIVE
 * (measured on a 825-tip tree at 1400 px: labels+domains 525, alignment 877, legend 249 -- the tree −271 px), the depth
 * scale clamped to zero, and the tree collapsed to a single vertical line. "Where is the tree?! There is only a line!"
 * <p>
 * <b>The rule.</b> The tree is allocated FIRST -- it is the point of the figure, so it is never the remainder. It gets
 * {@code tree_share} of the width, and more when the side components do not want it all. That share is a TARGET, not
 * an absolute floor: the parts that cannot be shrunk are honoured first -- the tip labels once the font auto-fit has
 * reached its readability floor, the clade bands, and the legend column, which is what keeps the legend off the
 * tracks -- so a crowded figure can still leave the tree under its share (measured at ~32% on an 825-tip tree
 * carrying long labels, domains, an alignment AND a legend). What it can no longer do is leave the tree NOTHING,
 * which is exactly what the independent per-component caps did. The side components
 * then share what is left: if their requests fit, everyone gets what they asked for; if not, each is first given its
 * MINIMUM and the remaining room is split between them in proportion to how much each asked for ABOVE that minimum.
 * That is deterministic, order-independent and monotone -- asking for more never gets you less.
 * <p>
 * Pure arithmetic: no Swing, no panel state, so it can be tested on its own.
 */
final class LayoutWidthBudget {

    /** The things that compete for the width beside the tree. */
    enum Part {
        /** Tip labels: the text, any tip-image slot, AND the protein-domain track, which rides the same
         *  reservation ({@code calculateLongestExtNodeInfo} folds it in) rather than being requested separately.
         *  Shrinks by reducing the font, not by clipping. */
        LABELS,
        /** The sequence-alignment window. Already scrolls, so a smaller grant simply shows fewer columns. */
        MSA,
        /** The stack of tip-aligned annotation columns. */
        ANNOTATIONS,
        /** The reserved column that keeps a right-edge legend off the other tracks. */
        LEGEND,
        /** Clade bars / brackets and their labels. */
        CLADE_BANDS
    }

    private final Map<Part, Integer> _granted = new EnumMap<Part, Integer>(Part.class);
    private final int                _tree_width;

    private LayoutWidthBudget(final Map<Part, Integer> granted, final int tree_width) {
        _granted.putAll(granted);
        _tree_width = tree_width;
    }

    /** The width granted to {@code part} (0 when it asked for nothing, or was never requested). */
    int granted(final Part part) {
        final Integer g = _granted.get(part);
        return (g == null) ? 0 : g.intValue();
    }

    /** The width left for the tree's own depth axis -- at least {@code tree_share} of the total whenever the
     *  minimums allow it, and more when the side components ask for less than their share of the room. */
    int treeWidth() {
        return _tree_width;
    }

    /** The sum of every grant. {@code treeWidth() + sideWidth() + fixed == total} holds whenever the requests
     *  fit or were squeezed to fit; it does NOT when even the minimums do not fit (the one branch that honours
     *  the minimums anyway and lets the tree take what is left), where the sum can exceed the usable width. */
    int sideWidth() {
        int sum = 0;
        for (final Integer g : _granted.values()) {
            sum += g.intValue();
        }
        return sum;
    }

    /** Collects the requests, then allocates. A part asks for {@code want} and will accept no less than {@code min}. */
    static final class Builder {

        private final Map<Part, int[]> _requests = new EnumMap<Part, int[]>(Part.class);

        /**
         * @param want what this part would draw at if nothing else competed (never negative)
         * @param min  the least it can usefully be drawn at; clamped into [0, want], so a part that wants nothing
         *             is granted nothing however large its stated minimum
         */
        Builder request(final Part part, final int want, final int min) {
            final int w = Math.max(0, want);
            _requests.put(part, new int[] { w, Math.max(0, Math.min(min, w)) });
            return this;
        }

        /**
         * @param total      the whole width being divided (a panel or viewport width)
         * @param fixed      width that is not negotiable and belongs to nobody (the fixed outer margin)
         * @param tree_share the share of {@code total} the tree is guaranteed, clamped into
         *                   [{@link AptxConstants#TREE_WIDTH_SHARE_MIN}, {@link AptxConstants#TREE_WIDTH_SHARE_MAX}]
         */
        LayoutWidthBudget allocate(final int total, final int fixed, final double tree_share) {
            final Map<Part, Integer> granted = new EnumMap<Part, Integer>(Part.class);
            final int usable = Math.max(0, total - Math.max(0, fixed));
            final double share = Math.min(AptxConstants.TREE_WIDTH_SHARE_MAX,
                    Math.max(AptxConstants.TREE_WIDTH_SHARE_MIN, tree_share));
            // The tree is allocated FIRST, so what the sides may spend is what remains after its guaranteed share.
            final int tree_floor = (int) Math.round(usable * share);
            final int side_budget = Math.max(0, usable - tree_floor);
            int want_total = 0;
            int min_total = 0;
            for (final int[] r : _requests.values()) {
                want_total += r[0];
                min_total += r[1];
            }
            if (want_total <= side_budget) {
                // Everything fits: nobody is squeezed, and the tree keeps the whole remainder -- which is MORE than
                // its share. Reserving exactly the share here instead would waste width on a tree with no side data.
                for (final Map.Entry<Part, int[]> e : _requests.entrySet()) {
                    granted.put(e.getKey(), Integer.valueOf(e.getValue()[0]));
                }
                return new LayoutWidthBudget(granted, Math.max(0, usable - want_total));
            }
            if (min_total >= side_budget) {
                // Even the minimums do not fit. Honour them anyway (a track drawn below its minimum is unreadable,
                // so there is nothing to be gained by shaving it further) and let the tree take what is left. This
                // is the one case where the tree can fall under its share, and it needs a panel narrower than the
                // minimums combined -- roughly a 300 px window with every track switched on.
                for (final Map.Entry<Part, int[]> e : _requests.entrySet()) {
                    granted.put(e.getKey(), Integer.valueOf(e.getValue()[1]));
                }
                return new LayoutWidthBudget(granted, Math.max(0, usable - min_total));
            }
            // The squeeze: everyone gets their minimum, and the room left over is split in proportion to how much
            // each asked for ABOVE it. Floor every share first so the total can only be under budget, then hand the
            // rounding remainder out one pixel at a time in enum order -- deterministic, and it never overspends.
            final int headroom = side_budget - min_total;
            final int want_extra_total = want_total - min_total; // > 0: want_total > side_budget >= min_total
            int spent = 0;
            for (final Part part : Part.values()) {
                final int[] r = _requests.get(part);
                if (r == null) {
                    continue;
                }
                final int extra = (int) (((long) headroom * (r[0] - r[1])) / want_extra_total);
                granted.put(part, Integer.valueOf(r[1] + extra));
                spent += r[1] + extra;
            }
            for (final Part part : Part.values()) {
                if (spent >= side_budget) {
                    break;
                }
                final int[] r = _requests.get(part);
                if ((r == null) || (granted.get(part).intValue() >= r[0])) {
                    continue; // already at what it asked for -- a leftover pixel would be over-granting
                }
                granted.put(part, Integer.valueOf(granted.get(part).intValue() + 1));
                ++spent;
            }
            return new LayoutWidthBudget(granted, Math.max(0, usable - spent));
        }
    }
}
