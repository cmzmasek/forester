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

package org.forester.pccx;

import java.util.SortedMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/*
 * Interface providing implementations of scoring methods used by
 * ExternalNodeBasedCoverageMethod.
 *
 * @author Christian M. Zmasek
 */
public interface ScoringMethodForExternalNode {

    /**
     * This calculates the coverage score for one external node.
     *
     *
     * @param external_node_scores
     *            SortedMap<PhylogenyNode, Double> in which the external node
     *            scores are stored (node->score)
     * @param phylogeny
     *            Phylogeny containing the external nodes to score
     * @param external_node
     *            PhylogenyNod for which to calculate the score
     * @param options
     *            CoverageCalculationOptions
     * @param annotate_phylogeny
     *
     */
    public void calculateScoreForExternalNode( final SortedMap<PhylogenyNode, Double> external_node_scores,
                                               final Phylogeny phylogeny,
                                               final PhylogenyNode external_node,
                                               final CoverageCalculationOptions options );

    /**
     * This returns a short description of this scoring method
     *
     * @return short description of this scoring method
     */
    public String getDesciption();

    /**
     * This calculates a normalization factor, so that a normalized score of 1.0
     * means complete coverage.
     *
     *
     * @param phylogeny
     *            Phylogeny containing the external nodes to score
     * @return normalization factor
     */
    public double getNormalizationFactor( final Phylogeny phylogeny );
}
