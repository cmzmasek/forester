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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/*
 *
 * @author Christian M. Zmasek
 */
public class BranchLengthBasedScoringMethod extends BranchCountingBasedScoringMethod {

    public static final double MIN_ALLOWED_BL_VALUE = 0.001;

    @Override
    double calculateScoreContributionPerExternalNode( final PhylogenyNode external_node,
                                                      final PhylogenyNode current_node ) {
        double score_contribution = 0.0;
        if ( current_node == external_node ) {
            score_contribution = external_node.getDistanceToParent();
            // This, of course, is completely /ad hoc/.
        }
        else {
            score_contribution = ModelingUtils.calculateBranchLengthSum( external_node, current_node );
        }
        return 1.0 / ( score_contribution > BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE ? score_contribution
                : BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE );
    }

    @Override
    public String getDesciption() {
        return "sum of 1/branch-length-sum [for self: 1/branch-length] [min branch length: "
                + BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE + "]";
    }

    @Override
    public double getNormalizationFactor( final Phylogeny phylogeny ) {
        double s = 0.0;
        double d = 0.0;
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorExternalForward(); iter.hasNext(); ) {
            d = iter.next().getDistanceToParent();
            s += ( 1.0 / ( d > BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE ? d
                    : BranchLengthBasedScoringMethod.MIN_ALLOWED_BL_VALUE ) );
        }
        return 1.0 / s;
    }
}
