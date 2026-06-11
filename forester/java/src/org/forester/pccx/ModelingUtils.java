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
import java.util.TreeMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/*
 * @author Christian M. Zmasek
 */
public final class ModelingUtils {

    static double calculateBranchLengthSum( final PhylogenyNode n1, final PhylogenyNode n2 ) {
        final PhylogenyNode lca = PhylogenyMethods.calculateLCA( n1, n2 );
        return ModelingUtils.calculateBranchLengthSumHelper( n1, lca )
                + ModelingUtils.calculateBranchLengthSumHelper( n2, lca );
    }

    private static double calculateBranchLengthSumHelper( final PhylogenyNode outer, final PhylogenyNode inner ) {
        PhylogenyNode my_outer = outer;
        double l = 0;
        while ( my_outer != inner ) {
            if ( my_outer.getDistanceToParent() > 0.0 ) {
                l += my_outer.getDistanceToParent();
            }
            my_outer = my_outer.getParent();
        }
        return l;
    }

    static int calculateBranchSum( final PhylogenyNode n1, final PhylogenyNode n2 ) {
        final PhylogenyNode lca = PhylogenyMethods.calculateLCA( n1, n2 );
        return ModelingUtils.calculateBranchSumHelper( n1, lca ) + ModelingUtils.calculateBranchSumHelper( n2, lca );
    }

    private static int calculateBranchSumHelper( final PhylogenyNode outer, final PhylogenyNode inner ) {
        PhylogenyNode my_outer = outer;
        int s = 0;
        while ( my_outer != inner ) {
            s++;
            my_outer = my_outer.getParent();
        }
        return s;
    }

    static SortedMap<PhylogenyNode, Double> setUpExternalCoverageHashMap( final Phylogeny phylogeny ) {
        final SortedMap<PhylogenyNode, Double> external_node_coverage = new TreeMap<PhylogenyNode, Double>();
        for( final PhylogenyNodeIterator iter = phylogeny.iteratorExternalForward(); iter.hasNext(); ) {
            external_node_coverage.put( iter.next(), 0.0 );
        }
        return external_node_coverage;
    }
}
