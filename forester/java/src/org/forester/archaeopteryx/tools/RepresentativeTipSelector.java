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

package org.forester.archaeopteryx.tools;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Tree-based dereplication: reduce redundancy by grouping the external nodes (tips) of a phylogeny into
 * clades whose members are all close in evolutionary distance, and keeping one representative per group.
 * This is the tree-native analogue of sequence-clustering tools such as cd-hit — but instead of realigning
 * sequences, it uses the tree the user already has: the patristic distance (sum of branch lengths along the
 * path between two tips) is the distance model.
 *
 * <p>Clustering is <em>complete-linkage on clades</em> (the "max diameter" method of TreeCluster): a group is
 * a maximal clade whose diameter — the largest patristic distance between any two of its tips — is at most the
 * cutoff. So the guarantee is honest: no two tips placed in the same group are farther apart than the cutoff.
 * Because a clade's diameter can only grow as you move rootward, the maximal qualifying clades are well defined
 * and the whole clustering is a single O(n) tree traversal.
 *
 * <p>Two entry points: {@link #selectByCutoff} (the user supplies the distance cutoff) and
 * {@link #selectByTargetCount} (the user supplies roughly how many representatives to keep; the cutoff that
 * comes closest to that count is found by exploiting the fact that the group count is monotonic in the cutoff).
 *
 * <p>If the tree carries no branch lengths, patristic distance is undefined, so the class falls back to
 * <em>topological</em> distance (one unit per edge) and flags it on the {@link SelectionResult}.
 *
 * <p>Pure and side-effect-free: it reads the tree and returns a {@link SelectionResult}; it never mutates the
 * phylogeny. The GUI glue highlights the representatives and can extract them into a new tree.
 */
public final class RepresentativeTipSelector {

    private static final double EPS = 1e-9;

    public enum RepresentativePick {
        /** The most central tip: the one with the smallest total patristic distance to the group's other tips. */
        MEDOID,
        /** The most divergent tip: the one with the longest terminal (external) branch. */
        LONGEST_BRANCH
    }

    private RepresentativeTipSelector() {
        // utility class
    }

    /**
     * Groups tips into maximal clades whose diameter is at most {@code cutoff} and keeps one representative per
     * group.
     *
     * @param phy    the phylogeny (its tips are the elements to dereplicate); not mutated
     * @param cutoff the maximum allowed distance between any two tips in a group (in the tree's branch-length
     *               unit, or edge count when the tree has no branch lengths); must be finite and &gt;= 0
     * @param pick   how to pick the representative within each group
     */
    public static SelectionResult selectByCutoff(final Phylogeny phy, final double cutoff,
                                                 final RepresentativePick pick) {
        return selectByCutoff(phy, cutoff, pick, Collections.emptySet());
    }

    /**
     * As {@link #selectByCutoff(Phylogeny, double, RepresentativePick)}, but the tips whose ids are in
     * {@code protected_tip_ids} are never dropped: a protected tip is always kept, and it stands in as its
     * group's representative (a group with several protected tips keeps them all). Useful for preserving
     * reference/type sequences the user has selected.
     */
    public static SelectionResult selectByCutoff(final Phylogeny phy, final double cutoff,
                                                 final RepresentativePick pick, final Set<Long> protected_tip_ids) {
        validate(phy);
        if (Double.isNaN(cutoff) || Double.isInfinite(cutoff) || (cutoff < 0.0)) {
            throw new IllegalArgumentException("cutoff must be a finite, non-negative number");
        }
        final boolean topological = !hasUsableBranchLengths(phy);
        final Map<Long, Double> diameter = computeDiameters(phy, topological);
        final List<PhylogenyNode> roots = clusterRoots(phy, diameter, cutoff);
        return buildResult(phy, roots, pick, topological, cutoff, -1, protected_tip_ids);
    }

    /**
     * Groups tips so that the number of representatives is as close as possible to {@code target}, and keeps one
     * representative per group. Because the achievable group counts are discrete (they change only at internal
     * clade diameters), the exact target may not be reachable; the count actually produced is reported on the
     * result. Ties in closeness favor the larger count (keeping more representatives).
     *
     * @param phy    the phylogeny; not mutated
     * @param target the desired number of representatives; must be &gt;= 1
     * @param pick   how to pick the representative within each group
     */
    public static SelectionResult selectByTargetCount(final Phylogeny phy, final int target,
                                                      final RepresentativePick pick) {
        return selectByTargetCount(phy, target, pick, Collections.emptySet());
    }

    /**
     * As {@link #selectByTargetCount(Phylogeny, int, RepresentativePick)}, but the tips whose ids are in
     * {@code protected_tip_ids} are never dropped. The clustering still aims for {@code target} groups; because
     * protected tips are always kept, the number of kept tips can exceed the target (reported on the result).
     */
    public static SelectionResult selectByTargetCount(final Phylogeny phy, final int target,
                                                      final RepresentativePick pick,
                                                      final Set<Long> protected_tip_ids) {
        validate(phy);
        if (target < 1) {
            throw new IllegalArgumentException("target number of representatives must be at least 1");
        }
        final boolean topological = !hasUsableBranchLengths(phy);
        final Map<Long, Double> diameter = computeDiameters(phy, topological);
        final int num_ext = phy.getNumberOfExternalNodes();
        final List<PhylogenyNode> roots;
        final double effective_cutoff;
        if (target >= num_ext) {
            // keep every tip: each is its own representative (a cutoff below every clade diameter)
            roots = externalNodes(phy);
            effective_cutoff = 0.0;
        }
        else if (target <= 1) {
            // one representative: the whole tree is a single group
            roots = new ArrayList<>();
            roots.add(phy.getRoot());
            effective_cutoff = diameter.get(phy.getRoot().getId());
        }
        else {
            final double cutoff = cutoffForTargetCount(phy, diameter, target, num_ext);
            roots = clusterRoots(phy, diameter, cutoff);
            // cutoffForTargetCount may return the -1 all-singletons sentinel; report 0 (each kept tip is its
            // own group, so the within-group distance is 0), never a nonsensical negative distance
            effective_cutoff = Math.max(0.0, cutoff);
        }
        return buildResult(phy, roots, pick, topological, effective_cutoff, target, protected_tip_ids);
    }

    /**
     * Finds the cutoff whose resulting group count is closest to {@code target}. The candidate cutoffs are the
     * distinct internal-clade diameters (the only values at which the count changes); the count is monotonically
     * non-increasing in the cutoff, so a binary search locates the boundary in O(n log n).
     */
    private static double cutoffForTargetCount(final Phylogeny phy, final Map<Long, Double> diameter,
                                               final int target, final int num_ext) {
        final double[] cand = sortedDistinctInternalDiameters(phy, diameter);
        // Find the smallest candidate whose count is <= target (such a candidate always exists: the root's
        // diameter yields a single group).
        int lo = 0;
        int hi = cand.length - 1;
        int boundary = cand.length - 1;
        while (lo <= hi) {
            final int mid = (lo + hi) >>> 1;
            if (countGroups(phy, diameter, cand[mid]) <= target) {
                boundary = mid;
                hi = mid - 1;
            }
            else {
                lo = mid + 1;
            }
        }
        final double t_high = cand[boundary];
        final int c_high = countGroups(phy, diameter, t_high); // <= target
        // The neighboring option that yields MORE groups than t_high: the previous candidate, or (below the
        // smallest internal diameter) every tip as its own group.
        final double t_low;
        final int c_low;
        if (boundary > 0) {
            t_low = cand[boundary - 1];
            c_low = countGroups(phy, diameter, t_low);
        }
        else {
            t_low = -1.0; // below every clade diameter -> all singletons
            c_low = num_ext;
        }
        // Pick the closer count; on a tie keep more representatives (the larger count = smaller cutoff).
        if (Math.abs(c_low - target) <= Math.abs(c_high - target)) {
            return t_low;
        }
        return t_high;
    }

    // --- clustering primitives -------------------------------------------------------------------------------

    /**
     * One post-order pass computing, for every node, the diameter of its subtree = the largest patristic
     * distance between any two tips below it. Also tracks each subtree's height (max branch-distance down to a
     * tip) as the cross-pair term.
     */
    private static Map<Long, Double> computeDiameters(final Phylogeny phy, final boolean topological) {
        final Map<Long, Double> height = new HashMap<>();
        final Map<Long, Double> diameter = new HashMap<>();
        for (final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext();) {
            final PhylogenyNode n = it.next();
            if (n.isExternal()) {
                height.put(n.getId(), 0.0);
                diameter.put(n.getId(), 0.0);
                continue;
            }
            double best1 = 0.0; // largest (child height + edge)
            double best2 = 0.0; // second largest
            double max_child_diameter = 0.0;
            int child_count = 0;
            for (int i = 0; i < n.getNumberOfDescendants(); ++i) {
                final PhylogenyNode c = n.getChildNode(i);
                final double reach = height.get(c.getId()) + edgeLength(c, topological);
                if (reach >= best1) {
                    best2 = best1;
                    best1 = reach;
                }
                else if (reach > best2) {
                    best2 = reach;
                }
                max_child_diameter = Math.max(max_child_diameter, diameter.get(c.getId()));
                ++child_count;
            }
            height.put(n.getId(), best1);
            // a tip pair crossing this node exists only with >= 2 children
            final double cross_pair = (child_count >= 2) ? (best1 + best2) : 0.0;
            diameter.put(n.getId(), Math.max(max_child_diameter, cross_pair));
        }
        return diameter;
    }

    /**
     * Top-down cut: the group roots are the highest nodes whose subtree diameter is at most {@code cutoff}
     * (a tip is always a group root, since it cannot be subdivided).
     */
    private static List<PhylogenyNode> clusterRoots(final Phylogeny phy, final Map<Long, Double> diameter,
                                                    final double cutoff) {
        final List<PhylogenyNode> roots = new ArrayList<>();
        final Deque<PhylogenyNode> stack = new ArrayDeque<>();
        stack.push(phy.getRoot());
        while (!stack.isEmpty()) {
            final PhylogenyNode n = stack.pop();
            if (n.isExternal() || (diameter.get(n.getId()) <= (cutoff + EPS))) {
                roots.add(n);
            }
            else {
                for (int i = 0; i < n.getNumberOfDescendants(); ++i) {
                    stack.push(n.getChildNode(i));
                }
            }
        }
        return roots;
    }

    private static int countGroups(final Phylogeny phy, final Map<Long, Double> diameter, final double cutoff) {
        int count = 0;
        final Deque<PhylogenyNode> stack = new ArrayDeque<>();
        stack.push(phy.getRoot());
        while (!stack.isEmpty()) {
            final PhylogenyNode n = stack.pop();
            if (n.isExternal() || (diameter.get(n.getId()) <= (cutoff + EPS))) {
                ++count;
            }
            else {
                for (int i = 0; i < n.getNumberOfDescendants(); ++i) {
                    stack.push(n.getChildNode(i));
                }
            }
        }
        return count;
    }

    private static double[] sortedDistinctInternalDiameters(final Phylogeny phy, final Map<Long, Double> diameter) {
        final List<Double> vals = new ArrayList<>();
        for (final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext();) {
            final PhylogenyNode n = it.next();
            if (!n.isExternal()) {
                vals.add(diameter.get(n.getId()));
            }
        }
        vals.sort(null);
        // de-duplicate
        final double[] tmp = new double[vals.size()];
        int m = 0;
        for (final double v : vals) {
            if ((m == 0) || (v > (tmp[m - 1] + EPS))) {
                tmp[m++] = v;
            }
        }
        final double[] out = new double[m];
        System.arraycopy(tmp, 0, out, 0, m);
        return out;
    }

    // --- result assembly -------------------------------------------------------------------------------------

    private static SelectionResult buildResult(final Phylogeny phy, final List<PhylogenyNode> roots,
                                               final RepresentativePick pick, final boolean topological,
                                               final double effective_cutoff, final int requested_target,
                                               final Set<Long> protected_tip_ids) {
        final List<Cluster> clusters = new ArrayList<>();
        int protected_kept = 0;
        for (final PhylogenyNode root : roots) {
            final List<PhylogenyNode> pre = preorderSubtree(root);
            final List<PhylogenyNode> members = new ArrayList<>();
            for (final PhylogenyNode n : pre) {
                if (n.isExternal()) {
                    members.add(n);
                }
            }
            final List<PhylogenyNode> protectd = new ArrayList<>();
            if (!protected_tip_ids.isEmpty()) {
                for (final PhylogenyNode m : members) {
                    if (protected_tip_ids.contains(m.getId())) {
                        protectd.add(m);
                    }
                }
            }
            final List<PhylogenyNode> kept;
            if (!protectd.isEmpty()) {
                // a protected tip is never dropped and stands in as its group's representative
                kept = protectd;
                protected_kept += protectd.size();
            }
            else {
                kept = List.of(chooseRepresentative(pre, members, root, pick, topological));
            }
            clusters.add(new Cluster(root, members, kept));
        }
        // present groups in a stable order (by representative node id, which for a freshly loaded tree
        // increases in tip order) so the result is deterministic and readable
        clusters.sort((a, b) -> Long.compare(a.getRepresentative().getId(), b.getRepresentative().getId()));
        return new SelectionResult(phy.getNumberOfExternalNodes(), clusters, effective_cutoff, topological, pick,
                requested_target, protected_kept);
    }

    private static PhylogenyNode chooseRepresentative(final List<PhylogenyNode> pre,
                                                      final List<PhylogenyNode> members,
                                                      final PhylogenyNode root, final RepresentativePick pick,
                                                      final boolean topological) {
        if (members.size() == 1) {
            return members.get(0);
        }
        if (pick == RepresentativePick.LONGEST_BRANCH) {
            return longestBranchMember(members, topological);
        }
        return medoidMember(pre, members, root, topological);
    }

    private static PhylogenyNode longestBranchMember(final List<PhylogenyNode> members, final boolean topological) {
        // members are in tree order, so the first-encountered wins ties -> deterministic across a tree copy
        PhylogenyNode best = members.get(0);
        double best_len = edgeLength(best, topological);
        for (int i = 1; i < members.size(); ++i) {
            final double len = edgeLength(members.get(i), topological);
            if (len > (best_len + EPS)) {
                best = members.get(i);
                best_len = len;
            }
        }
        return best;
    }

    /**
     * The medoid = the tip with the smallest total patristic distance to the other tips in its group, computed
     * with an O(group size) sum-of-distances rerooting DP (never the O(group size^2) all-pairs sum), so a huge
     * group produced by an aggressive cutoff / tiny target does not stall.
     */
    private static PhylogenyNode medoidMember(final List<PhylogenyNode> pre, final List<PhylogenyNode> members,
                                              final PhylogenyNode root, final boolean topological) {
        final Map<Long, Integer> subtree_tips = new HashMap<>();
        final Map<Long, Double> down = new HashMap<>(); // sum of distances from the node down to tips below it
        // post-order (reverse of pre-order)
        for (int i = pre.size() - 1; i >= 0; --i) {
            final PhylogenyNode n = pre.get(i);
            if (n.isExternal()) {
                subtree_tips.put(n.getId(), 1);
                down.put(n.getId(), 0.0);
                continue;
            }
            int s = 0;
            double d = 0.0;
            for (int c = 0; c < n.getNumberOfDescendants(); ++c) {
                final PhylogenyNode child = n.getChildNode(c);
                final double e = edgeLength(child, topological);
                final int cs = subtree_tips.get(child.getId());
                s += cs;
                d += down.get(child.getId()) + (e * cs);
            }
            subtree_tips.put(n.getId(), s);
            down.put(n.getId(), d);
        }
        final int total = subtree_tips.get(root.getId());
        final Map<Long, Double> up = new HashMap<>(); // sum of distances from the node up to tips NOT below it
        up.put(root.getId(), 0.0);
        // pre-order
        for (final PhylogenyNode n : pre) {
            final double un = up.get(n.getId());
            final double dn = down.get(n.getId());
            for (int c = 0; c < n.getNumberOfDescendants(); ++c) {
                final PhylogenyNode child = n.getChildNode(c);
                final double e = edgeLength(child, topological);
                final int cs = subtree_tips.get(child.getId());
                final double uc = un + dn - down.get(child.getId()) - (e * cs) + (e * (total - cs));
                up.put(child.getId(), uc);
            }
        }
        // a tip's total distance to all group-mates is down(0) + up = up
        PhylogenyNode best = members.get(0);
        double best_total = up.get(best.getId());
        for (int i = 1; i < members.size(); ++i) {
            final double t = up.get(members.get(i).getId());
            if (t < (best_total - EPS)) {
                best = members.get(i);
                best_total = t;
            }
        }
        return best;
    }

    // --- small helpers ---------------------------------------------------------------------------------------

    private static List<PhylogenyNode> preorderSubtree(final PhylogenyNode root) {
        final List<PhylogenyNode> out = new ArrayList<>();
        final Deque<PhylogenyNode> stack = new ArrayDeque<>();
        stack.push(root);
        while (!stack.isEmpty()) {
            final PhylogenyNode n = stack.pop();
            out.add(n);
            for (int i = n.getNumberOfDescendants() - 1; i >= 0; --i) {
                stack.push(n.getChildNode(i));
            }
        }
        return out;
    }

    private static List<PhylogenyNode> externalNodes(final Phylogeny phy) {
        return new ArrayList<>(phy.getExternalNodes());
    }

    /** Branch length to the parent, clamped to be non-negative; 1 unit per edge when using topological distance. */
    private static double edgeLength(final PhylogenyNode n, final boolean topological) {
        if (topological) {
            return 1.0;
        }
        final double d = n.getDistanceToParent();
        return (d > 0.0) ? d : 0.0;
    }

    /**
     * True if the tree carries at least one real branch length. When false, distance-cutoff selection is
     * meaningless (the GUI hides it) and only target-count selection — which clusters by topological distance —
     * makes sense.
     */
    public static boolean hasUsableBranchLengths(final Phylogeny phy) {
        for (final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext();) {
            final PhylogenyNode n = it.next();
            if (!n.isRoot() && (n.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT)) {
                return true;
            }
        }
        return false;
    }

    private static void validate(final Phylogeny phy) {
        if ((phy == null) || phy.isEmpty()) {
            throw new IllegalArgumentException("phylogeny is null or empty");
        }
        if (phy.getNumberOfExternalNodes() < 1) {
            throw new IllegalArgumentException("phylogeny has no external nodes");
        }
    }

    static String formatDistance(final double d) {
        final double r = Math.round(d * 1e5) / 1e5;
        if (r == Math.rint(r)) {
            return String.valueOf((long) r);
        }
        return String.valueOf(r);
    }

    // --- value types -----------------------------------------------------------------------------------------

    /** One group of tips and the tip(s) kept to represent them. */
    public static final class Cluster {

        private final PhylogenyNode _clade_root;
        private final List<PhylogenyNode> _members;
        private final List<PhylogenyNode> _kept_members;

        Cluster(final PhylogenyNode clade_root, final List<PhylogenyNode> members,
                final List<PhylogenyNode> kept_members) {
            _clade_root = clade_root;
            _members = List.copyOf(members);
            _kept_members = List.copyOf(kept_members);
        }

        /** The clade whose subtree contains exactly this group's tips. */
        public PhylogenyNode getCladeRoot() {
            return _clade_root;
        }

        /** The tips in this group (&gt;= 1), in tree order. */
        public List<PhylogenyNode> getMembers() {
            return _members;
        }

        /**
         * The tips kept from this group: normally the single computed representative, or — when the group
         * contains protected tips — all of those protected tips.
         */
        public List<PhylogenyNode> getKeptMembers() {
            return _kept_members;
        }

        /** The primary kept tip (the representative, or the first protected tip). */
        public PhylogenyNode getRepresentative() {
            return _kept_members.get(0);
        }

        public int size() {
            return _members.size();
        }
    }

    /** The outcome of a representative-selection run. Immutable. */
    public static final class SelectionResult {

        private final int _total_tips;
        private final List<Cluster> _clusters;
        private final double _effective_cutoff;
        private final boolean _used_topological_distance;
        private final RepresentativePick _pick;
        private final int _requested_target; // -1 for cutoff mode
        private final int _protected_kept;

        SelectionResult(final int total_tips, final List<Cluster> clusters, final double effective_cutoff,
                        final boolean used_topological_distance, final RepresentativePick pick,
                        final int requested_target, final int protected_kept) {
            _total_tips = total_tips;
            _clusters = List.copyOf(clusters);
            _effective_cutoff = effective_cutoff;
            _used_topological_distance = used_topological_distance;
            _pick = pick;
            _requested_target = requested_target;
            _protected_kept = protected_kept;
        }

        public int getTotalTips() {
            return _total_tips;
        }

        public List<Cluster> getClusters() {
            return _clusters;
        }

        /** Number of groups the tips were clustered into. */
        public int getClusterCount() {
            return _clusters.size();
        }

        /**
         * Total number of tips kept: one representative per group, plus any extra protected tips (a group with
         * two or more protected tips keeps them all). Equals {@link #getClusterCount()} when nothing is
         * protected.
         */
        public int getKeptCount() {
            int n = 0;
            for (final Cluster c : _clusters) {
                n += c.getKeptMembers().size();
            }
            return n;
        }

        /** Number of kept tips that were protected from removal. */
        public int getProtectedKeptCount() {
            return _protected_kept;
        }

        /** The diameter threshold that produced these groups (0 when every tip is kept). */
        public double getEffectiveCutoff() {
            return _effective_cutoff;
        }

        public boolean usedTopologicalDistance() {
            return _used_topological_distance;
        }

        public RepresentativePick getPick() {
            return _pick;
        }

        /** The primary representative tip of each group, one per group, in tree order. */
        public List<PhylogenyNode> getRepresentatives() {
            final List<PhylogenyNode> reps = new ArrayList<>(_clusters.size());
            for (final Cluster c : _clusters) {
                reps.add(c.getRepresentative());
            }
            return reps;
        }

        /** All tips kept, in tree order: the representatives plus any protected tips. */
        public List<PhylogenyNode> getKeptTips() {
            final List<PhylogenyNode> kept = new ArrayList<>();
            for (final Cluster c : _clusters) {
                kept.addAll(c.getKeptMembers());
            }
            return kept;
        }

        /** Ids of all kept tips (representatives plus protected tips), valid on the phylogeny computed from. */
        public Set<Long> representativeIds() {
            final Set<Long> ids = new HashSet<>();
            for (final Cluster c : _clusters) {
                for (final PhylogenyNode k : c.getKeptMembers()) {
                    ids.add(k.getId());
                }
            }
            return ids;
        }

        public String summary() {
            final StringBuilder sb = new StringBuilder();
            sb.append("Grouped ").append(_total_tips).append(_total_tips == 1 ? " tip" : " tips")
                    .append(" into ").append(getClusterCount())
                    .append(getClusterCount() == 1 ? " group." : " groups.");
            if ((_requested_target > 0) && (getClusterCount() != _requested_target)) {
                sb.append(" (requested ").append(_requested_target).append(")");
            }
            sb.append("\nEach group's tips are within a distance of ")
                    .append(formatDistance(_effective_cutoff)).append(" of each other");
            if (_used_topological_distance) {
                sb.append(" (topological distance — the tree has no branch lengths)");
            }
            sb.append(".");
            if (_protected_kept > 0) {
                sb.append("\nKeeping ").append(getKeptCount()).append(getKeptCount() == 1 ? " tip, including "
                        : " tips, including ").append(_protected_kept)
                        .append(_protected_kept == 1 ? " selected tip protected from removal."
                                : " selected tips protected from removal.");
            }
            sb.append("\nRepresentative per group: ")
                    .append(_pick == RepresentativePick.LONGEST_BRANCH ? "most divergent (longest branch)"
                            : "most central (medoid)")
                    .append(".");
            return sb.toString();
        }
    }
}
