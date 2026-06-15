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

package org.forester.analysis;

import java.util.ArrayList;
import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public final class AncestralTaxonomyInference {

    public static void inferTaxonomyFromDescendents(final Phylogeny phy) throws AncestralTaxonomyInferenceException {
        for (final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if (!node.isExternal()) {
                inferTaxonomyFromDescendentsSimple(node);
            }
        }
    }

    private static void inferTaxonomyFromDescendentsSimple(final PhylogenyNode n) throws AncestralTaxonomyInferenceException {
        if (n.isExternal()) {
            throw new IllegalArgumentException("attempt to infer taxonomy from descendants of external node");
        }
        n.getNodeData().setTaxonomy(null);
        final PhylogenyNode n1 = n.getChildNode1();
        final PhylogenyNode n2 = n.getChildNode2();

        if (n1.getNodeData().isHasTaxonomy() && n2.getNodeData().isHasTaxonomy() && n1.getNodeData().getTaxonomy().getLineage() != null && n2.getNodeData().getTaxonomy().getLineage() != null) {
            final List<String> lin1 = n1.getNodeData().getTaxonomy().getLineage();
            final List<String> lin2 = n2.getNodeData().getTaxonomy().getLineage();

            final List<String> lin = new ArrayList<>();
            final int s = lin1.size() <= lin2.size() ? lin1.size() : lin2.size();

            for (int i = 0; i < s; ++i) {
                if (lin1.get(i).equalsIgnoreCase(lin2.get(i))) {
                    lin.add(lin1.get(i));
                } else {
                    break;
                }
            }
            final Taxonomy tax = new Taxonomy();
            tax.setLineage(lin);

            if ( lin.size() > 0) {
                final String last = lin.get(lin.size() - 1);
                tax.setScientificName(last);
                if (n1.getNodeData().getTaxonomy().getScientificName().equals(last) ) {
                    if (n1.isInternal()) {
                        n1.getNodeData().getTaxonomy().setScientificName(null);
                    }
                }
                if (n2.getNodeData().getTaxonomy().getScientificName().equals(last) ) {
                    if (n2.isInternal()) {
                        n2.getNodeData().getTaxonomy().setScientificName(null);
                    }
                }
            }
            n.getNodeData().setTaxonomy(tax);
        } else if (n1.getNodeData().isHasTaxonomy() && n1.getNodeData().getTaxonomy().getLineage() != null) {
            final Taxonomy tax = new Taxonomy();
            tax.setLineage(n1.getNodeData().getTaxonomy().getLineage());
            n.getNodeData().setTaxonomy(tax);
        } else if (n2.getNodeData().isHasTaxonomy() && n2.getNodeData().getTaxonomy().getLineage() != null) {
            final Taxonomy tax = new Taxonomy();
            tax.setLineage(n2.getNodeData().getTaxonomy().getLineage());
            n.getNodeData().setTaxonomy(tax);
        }
    }

}
