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

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.forester.archaeopteryx.tools.RepresentativeTipSelector.RepresentativePick;
import org.forester.archaeopteryx.tools.RepresentativeTipSelector.SelectionResult;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;

public class RepresentativeTipSelectorTest {

    public static boolean test() {
        return cutoffClustering() && monotonicity() && completeLinkageGuarantee() && targetCount() && medoid()
                && longestBranch() && topologicalFallback() && extractionDeterminism()
                && copyPreservesExternalOrder() && protection() && accessors() && edgeCases();
    }

    public static void main(final String[] args) {
        System.out.println(test() ? "OK" : "FAILED");
    }

    /** Two balanced cherries far apart: root d=0.82, each cherry d=0.02. */
    private static Phylogeny twoCherries() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode x = child(root, null, 0.4);
        final PhylogenyNode y = child(root, null, 0.4);
        child(x, "A", 0.01);
        child(x, "B", 0.01);
        child(y, "C", 0.01);
        child(y, "D", 0.01);
        return finish(root);
    }

    /**
     * Caterpillar with asymmetric branches so the medoid is unique: leaf A hangs far out on a short branch,
     * and C sits on the shortest terminal branch nearest the tree's center, making C the sole minimum-total tip
     * (totals A=1.4, B=1.1, C=1.0, D=1.1).
     */
    private static Phylogeny caterpillar() {
        final PhylogenyNode root = new PhylogenyNode();
        child(root, "A", 0.05);
        final PhylogenyNode m = child(root, null, 0.2);
        child(m, "B", 0.1);
        final PhylogenyNode n = child(m, null, 0.2);
        child(n, "C", 0.05);
        child(n, "D", 0.1);
        return finish(root);
    }

    /** Nested unbalanced tree whose group counts are {4, 3, 2, 1} (used for target-count). */
    private static Phylogeny nested() {
        final PhylogenyNode root = new PhylogenyNode();
        child(root, "A", 0.01);
        final PhylogenyNode p = child(root, null, 0.5);
        child(p, "B", 0.01);
        final PhylogenyNode q = child(p, null, 0.5);
        child(q, "C", 0.01);
        child(q, "D", 0.01);
        return finish(root);
    }

    private static boolean cutoffClustering() {
        final Phylogeny phy = twoCherries();
        final SelectionResult r = RepresentativeTipSelector.selectByCutoff(phy, 0.02, RepresentativePick.MEDOID);
        if (r.getClusterCount() != 2) {
            return fail("cutoff 0.02 should give 2 groups, got " + r.getClusterCount());
        }
        if (r.usedTopologicalDistance()) {
            return fail("tree has branch lengths; should not be topological");
        }
        for (final RepresentativeTipSelector.Cluster c : r.getClusters()) {
            if (c.size() != 2) {
                return fail("each group should have 2 members, got " + c.size());
            }
        }
        if (r.getTotalTips() != 4) {
            return fail("total tips should be 4");
        }
        // one representative per group, all distinct
        if (r.representativeIds().size() != 2) {
            return fail("expected 2 distinct representatives");
        }
        // whole tree in one group at/above the root diameter
        if (RepresentativeTipSelector.selectByCutoff(phy, 1.0, RepresentativePick.MEDOID).getClusterCount() != 1) {
            return fail("cutoff 1.0 should give 1 group");
        }
        // every tip its own group at cutoff 0 (all internal diameters are positive)
        if (RepresentativeTipSelector.selectByCutoff(phy, 0.0, RepresentativePick.MEDOID).getClusterCount() != 4) {
            return fail("cutoff 0 should give 4 singleton groups");
        }
        return true;
    }

    private static boolean monotonicity() {
        final Phylogeny phy = nested();
        int prev = Integer.MAX_VALUE;
        // group count must be non-increasing as the cutoff grows
        for (final double t : new double[] { 0.0, 0.01, 0.02, 0.3, 0.52, 0.9, 1.02, 5.0 }) {
            final int c = RepresentativeTipSelector.selectByCutoff(phy, t, RepresentativePick.MEDOID)
                    .getClusterCount();
            if (c > prev) {
                return fail("group count increased with cutoff at t=" + t + " (" + c + " > " + prev + ")");
            }
            prev = c;
        }
        return true;
    }

    private static boolean completeLinkageGuarantee() {
        final Phylogeny phy = nested();
        final double cutoff = 0.52;
        final SelectionResult r = RepresentativeTipSelector.selectByCutoff(phy, cutoff, RepresentativePick.MEDOID);
        for (final RepresentativeTipSelector.Cluster c : r.getClusters()) {
            final double d = maxPairwise(c.getMembers());
            if (d > (cutoff + 1e-9)) {
                return fail("group violates the cutoff: max pairwise " + d + " > " + cutoff);
            }
            // the representative must be a member of its group
            if (!c.getMembers().contains(c.getRepresentative())) {
                return fail("representative is not a member of its group");
            }
        }
        return true;
    }

    private static boolean targetCount() {
        final Phylogeny phy = nested(); // achievable counts: 4, 3, 2, 1
        if (count(phy, 3) != 3) {
            return fail("target 3 should be reachable exactly, got " + count(phy, 3));
        }
        if (count(phy, 2) != 2) {
            return fail("target 2 should give 2, got " + count(phy, 2));
        }
        if (count(phy, 1) != 1) {
            return fail("target 1 should give 1, got " + count(phy, 1));
        }
        if (count(phy, 4) != 4) {
            return fail("target 4 (== #tips) should keep all 4");
        }
        if (count(phy, 9) != 4) {
            return fail("target beyond #tips should keep all 4");
        }
        // tie in closeness favors keeping MORE representatives: on twoCherries the counts are {4,2,1};
        // target 3 is equidistant from 4 and 2, so it should choose 4
        if (count(twoCherries(), 3) != 4) {
            return fail("tie should favor more representatives (expected 4), got " + count(twoCherries(), 3));
        }
        // requested-vs-produced is reported
        final SelectionResult r = RepresentativeTipSelector.selectByTargetCount(twoCherries(), 3,
                RepresentativePick.MEDOID);
        if (!r.summary().contains("requested 3")) {
            return fail("summary should note the requested target when it differs");
        }
        return true;
    }

    private static boolean medoid() {
        // one group of all four caterpillar leaves; the medoid must minimize total distance
        final Phylogeny phy = caterpillar();
        final SelectionResult r = RepresentativeTipSelector.selectByCutoff(phy, 10.0, RepresentativePick.MEDOID);
        if (r.getClusterCount() != 1) {
            return fail("cutoff 10 should give 1 group");
        }
        final RepresentativeTipSelector.Cluster c = r.getClusters().get(0);
        final PhylogenyNode picked = c.getRepresentative();
        // the DP-chosen medoid must achieve the brute-force minimum total distance
        final double picked_total = totalDistance(picked, c.getMembers());
        double min_total = Double.MAX_VALUE;
        for (final PhylogenyNode m : c.getMembers()) {
            min_total = Math.min(min_total, totalDistance(m, c.getMembers()));
        }
        if (Math.abs(picked_total - min_total) > 1e-9) {
            return fail("medoid total " + picked_total + " != brute-force min " + min_total);
        }
        // C is the unique minimum-total tip (not the first tip A) -> pins the DP pick, not just a tie
        if (!"C".equals(picked.getName())) {
            return fail("medoid should be C (the unique minimum-total tip), got " + picked.getName());
        }
        return true;
    }

    private static boolean longestBranch() {
        final Phylogeny phy = caterpillar();
        final SelectionResult r = RepresentativeTipSelector.selectByCutoff(phy, 10.0,
                RepresentativePick.LONGEST_BRANCH);
        final RepresentativeTipSelector.Cluster c = r.getClusters().get(0);
        // representative must have the maximal terminal branch length in the group
        double max_terminal = 0.0;
        for (final PhylogenyNode m : c.getMembers()) {
            max_terminal = Math.max(max_terminal, m.getDistanceToParent());
        }
        if (Math.abs(c.getRepresentative().getDistanceToParent() - max_terminal) > 1e-9) {
            return fail("longest-branch representative should have the max terminal branch");
        }
        return true;
    }

    private static boolean topologicalFallback() {
        // same shape as twoCherries but with NO branch lengths
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode x = new PhylogenyNode();
        final PhylogenyNode y = new PhylogenyNode();
        root.addAsChild(x);
        root.addAsChild(y);
        x.addAsChild(named("A"));
        x.addAsChild(named("B"));
        y.addAsChild(named("C"));
        y.addAsChild(named("D"));
        final Phylogeny phy = finish(root);
        if (RepresentativeTipSelector.hasUsableBranchLengths(phy)) {
            return fail("tree without branch lengths should report none usable");
        }
        // topological distance: within-cherry diameter = 2 edges, cross = 4 edges
        final SelectionResult at1 = RepresentativeTipSelector.selectByCutoff(phy, 1.0, RepresentativePick.MEDOID);
        if (!at1.usedTopologicalDistance()) {
            return fail("expected topological distance mode");
        }
        if (at1.getClusterCount() != 4) {
            return fail("topological cutoff 1 should give 4 singletons, got " + at1.getClusterCount());
        }
        if (RepresentativeTipSelector.selectByCutoff(phy, 2.0, RepresentativePick.MEDOID).getClusterCount() != 2) {
            return fail("topological cutoff 2 should give 2 groups");
        }
        if (!at1.summary().contains("topological")) {
            return fail("summary should mention topological distance");
        }
        return true;
    }

    private static boolean extractionDeterminism() {
        // the representatives selected on a tree and on its copy must be the same tips (by name), so the GUI
        // can safely recompute on a copy to build the extracted tree
        final Phylogeny phy = nested();
        final Phylogeny copy = phy.copy();
        final Set<String> a = repNames(
                RepresentativeTipSelector.selectByCutoff(phy, 0.52, RepresentativePick.MEDOID));
        final Set<String> b = repNames(
                RepresentativeTipSelector.selectByCutoff(copy, 0.52, RepresentativePick.MEDOID));
        if (!a.equals(b)) {
            return fail("representatives differ between a tree and its copy: " + a + " vs " + b);
        }
        return true;
    }

    private static boolean copyPreservesExternalOrder() {
        // the GUI extraction index-maps kept tips onto phy.copy(); that is correct only if copy() preserves
        // external-node order -- pin that property directly
        final Phylogeny phy = nested();
        final Phylogeny copy = phy.copy();
        final List<PhylogenyNode> a = phy.getExternalNodes();
        final List<PhylogenyNode> b = copy.getExternalNodes();
        if (a.size() != b.size()) {
            return fail("copy changed the number of external nodes");
        }
        for (int i = 0; i < a.size(); ++i) {
            final String an = a.get(i).getName();
            final String bn = b.get(i).getName();
            if ((an == null) ? (bn != null) : !an.equals(bn)) {
                return fail("copy did not preserve external-node order at index " + i + ": " + an + " vs " + bn);
            }
        }
        return true;
    }

    private static boolean accessors() {
        final Phylogeny phy = twoCherries();
        final SelectionResult r = RepresentativeTipSelector.selectByCutoff(phy, 0.02, RepresentativePick.MEDOID);
        if ((r.getKeptTips().size() != r.getKeptCount()) || (r.getKeptTips().size() != 2)) {
            return fail("accessors: getKeptTips should aggregate all kept members");
        }
        if (r.getPick() != RepresentativePick.MEDOID) {
            return fail("accessors: getPick");
        }
        if (r.getEffectiveCutoff() < 0.0) {
            return fail("accessors: effective cutoff must be non-negative");
        }
        for (final RepresentativeTipSelector.Cluster c : r.getClusters()) {
            if (c.getCladeRoot() == null) {
                return fail("accessors: getCladeRoot null");
            }
            if (c.getKeptMembers().isEmpty() || !c.getMembers().containsAll(c.getKeptMembers())) {
                return fail("accessors: kept members must be a subset of the group's members");
            }
        }
        // the all-singletons fallback in the target-count SEARCH branch must report cutoff 0, not the -1
        // sentinel (twoCherries target 3 is equidistant from counts 4 and 2, so it falls back to 4 singletons)
        final SelectionResult sentinel = RepresentativeTipSelector.selectByTargetCount(twoCherries(), 3,
                RepresentativePick.MEDOID);
        if ((sentinel.getEffectiveCutoff() != 0.0) || (sentinel.getClusterCount() != 4)) {
            return fail("accessors: all-singletons target path must report cutoff 0, got "
                    + sentinel.getEffectiveCutoff());
        }
        // formatDistance renders integers cleanly and strips trailing zeros
        if (!"0".equals(RepresentativeTipSelector.formatDistance(0.0))
                || !"2".equals(RepresentativeTipSelector.formatDistance(2.0))
                || !"0.05".equals(RepresentativeTipSelector.formatDistance(0.05))) {
            return fail("accessors: formatDistance");
        }
        return true;
    }

    private static boolean protection() {
        // (1) a protected tip is kept even when it is not the medoid; it replaces the group's representative
        final Phylogeny cat = caterpillar();
        final SelectionResult r1 = RepresentativeTipSelector.selectByCutoff(cat, 10.0, RepresentativePick.MEDOID,
                ids(tip(cat, "A")));
        if ((r1.getClusterCount() != 1) || (r1.getKeptCount() != 1) || (r1.getProtectedKeptCount() != 1)) {
            return fail("protection(1): expected one group, one kept, one protected");
        }
        if (!r1.representativeIds().equals(ids(tip(cat, "A")))) {
            return fail("protection(1): the protected outlier A must be the sole kept tip");
        }
        // grammar: a single kept/protected tip must read singular ("1 tip", not "1 tips")
        if (!r1.summary().contains("Keeping 1 tip, including 1 selected tip protected")) {
            return fail("protection(1): summary should be singular for one kept/protected tip");
        }
        // (2) several protected tips in one group are all kept
        final Phylogeny cat2 = caterpillar();
        final Set<Long> prot2 = ids(tip(cat2, "A"), tip(cat2, "C"));
        final SelectionResult r2 = RepresentativeTipSelector.selectByCutoff(cat2, 10.0, RepresentativePick.MEDOID,
                prot2);
        if ((r2.getClusterCount() != 1) || (r2.getKeptCount() != 2) || (r2.getProtectedKeptCount() != 2)) {
            return fail("protection(2): both protected tips should be kept in the single group");
        }
        if (!r2.representativeIds().equals(prot2)) {
            return fail("protection(2): kept set should be exactly the protected tips");
        }
        // (3) protection can push the kept count above the target
        final Phylogeny nes = nested();
        final Set<Long> prot3 = ids(tip(nes, "A"), tip(nes, "B"), tip(nes, "C"));
        final SelectionResult r3 = RepresentativeTipSelector.selectByTargetCount(nes, 1, RepresentativePick.MEDOID,
                prot3);
        if ((r3.getClusterCount() != 1) || (r3.getKeptCount() != 3) || (r3.getProtectedKeptCount() != 3)) {
            return fail("protection(3): kept count should exceed the target (protection wins)");
        }
        // (4) in a multi-group tree, protecting a non-representative swaps only that group's kept tip
        final Phylogeny tc = twoCherries();
        final SelectionResult r4 = RepresentativeTipSelector.selectByCutoff(tc, 0.02, RepresentativePick.MEDOID,
                ids(tip(tc, "B")));
        if ((r4.getClusterCount() != 2) || (r4.getKeptCount() != 2)) {
            return fail("protection(4): still two groups, two kept");
        }
        if (!r4.representativeIds().contains(tip(tc, "B").getId())
                || r4.representativeIds().contains(tip(tc, "A").getId())) {
            return fail("protection(4): protected B should replace default representative A");
        }
        // (5) an empty protection set behaves exactly like the un-protected call
        final Phylogeny nes2 = nested();
        final SelectionResult plain = RepresentativeTipSelector.selectByCutoff(nes2, 0.52, RepresentativePick.MEDOID);
        final SelectionResult empty = RepresentativeTipSelector.selectByCutoff(nes2, 0.52, RepresentativePick.MEDOID,
                new HashSet<>());
        if (!plain.representativeIds().equals(empty.representativeIds()) || (empty.getProtectedKeptCount() != 0)) {
            return fail("protection(5): empty protection must match no-protection");
        }
        return true;
    }

    private static boolean edgeCases() {
        // single-tip tree
        final PhylogenyNode only = new PhylogenyNode();
        only.setName("solo");
        final Phylogeny one = new Phylogeny();
        one.setRoot(only);
        one.externalNodesHaveChanged();
        if (RepresentativeTipSelector.selectByCutoff(one, 0.5, RepresentativePick.MEDOID).getClusterCount() != 1) {
            return fail("single-tip tree should give 1 group");
        }
        // two-tip cherry
        final PhylogenyNode root = new PhylogenyNode();
        child(root, "A", 0.1);
        child(root, "B", 0.1);
        final Phylogeny two = finish(root);
        if (RepresentativeTipSelector.selectByCutoff(two, 0.1, RepresentativePick.MEDOID).getClusterCount() != 2) {
            return fail("cherry below its diameter should give 2 singletons");
        }
        if (RepresentativeTipSelector.selectByCutoff(two, 0.2, RepresentativePick.MEDOID).getClusterCount() != 1) {
            return fail("cherry at its diameter should give 1 group");
        }
        // invalid inputs
        if (!throwsIae(() -> RepresentativeTipSelector.selectByCutoff(null, 0.1, RepresentativePick.MEDOID))) {
            return fail("null phylogeny should throw");
        }
        if (!throwsIae(() -> RepresentativeTipSelector.selectByCutoff(new Phylogeny(), 0.1, RepresentativePick.MEDOID))) {
            return fail("empty phylogeny should throw");
        }
        if (!throwsIae(() -> RepresentativeTipSelector.selectByCutoff(two, -1.0, RepresentativePick.MEDOID))) {
            return fail("negative cutoff should throw");
        }
        if (!throwsIae(() -> RepresentativeTipSelector.selectByTargetCount(two, 0, RepresentativePick.MEDOID))) {
            return fail("target 0 should throw");
        }
        return true;
    }

    // --- helpers -------------------------------------------------------------------------------------------

    private static int count(final Phylogeny phy, final int target) {
        return RepresentativeTipSelector.selectByTargetCount(phy, target, RepresentativePick.MEDOID)
                .getClusterCount();
    }

    private static PhylogenyNode tip(final Phylogeny phy, final String name) {
        for (final PhylogenyNode n : phy.getExternalNodes()) {
            if (name.equals(n.getName())) {
                return n;
            }
        }
        return null;
    }

    private static Set<Long> ids(final PhylogenyNode... nodes) {
        final Set<Long> s = new HashSet<>();
        for (final PhylogenyNode n : nodes) {
            s.add(n.getId());
        }
        return s;
    }

    private static PhylogenyNode child(final PhylogenyNode parent, final String name, final double dist) {
        final PhylogenyNode n = new PhylogenyNode();
        if (name != null) {
            n.setName(name);
        }
        parent.addAsChild(n);
        n.setDistanceToParent(dist);
        return n;
    }

    private static PhylogenyNode named(final String name) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName(name);
        return n;
    }

    private static Phylogeny finish(final PhylogenyNode root) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot(root);
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static double maxPairwise(final List<PhylogenyNode> members) {
        double m = 0.0;
        for (int i = 0; i < members.size(); ++i) {
            for (int j = i + 1; j < members.size(); ++j) {
                m = Math.max(m, PhylogenyMethods.calculateDistance(members.get(i), members.get(j)));
            }
        }
        return m;
    }

    private static double totalDistance(final PhylogenyNode a, final List<PhylogenyNode> members) {
        double t = 0.0;
        for (final PhylogenyNode b : members) {
            if (a != b) {
                t += PhylogenyMethods.calculateDistance(a, b);
            }
        }
        return t;
    }

    private static Set<String> repNames(final SelectionResult r) {
        final Set<String> names = new HashSet<>();
        for (final PhylogenyNode n : r.getRepresentatives()) {
            names.add(n.getName());
        }
        return names;
    }

    private static boolean throwsIae(final Runnable r) {
        try {
            r.run();
            return false;
        } catch (final IllegalArgumentException e) {
            return true;
        }
    }

    private static boolean fail(final String msg) {
        System.out.println("RepresentativeTipSelectorTest failed: " + msg);
        return false;
    }
}
