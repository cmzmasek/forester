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

import java.util.ArrayList;
import java.util.List;

import org.forester.util.ForesterConstants;

/**
 * The main literature citations for the algorithms implemented in Archaeopteryx, shown by the
 * Help &rarr; References dialog. Kept in one place (one {@link Reference} per algorithm) so adding a new
 * algorithm's citation is a one-line change and a later "make it consistent/pretty" pass edits only this
 * list and {@link #asText()} -- not any call site.
 *
 * <p>Most entries are literature citations; an algorithm with no separate publication (e.g. representative-tip
 * selection) carries a short prose description of the method instead.
 */
final class AlgorithmReferences {

    /** One algorithm and its main citation (or, for an unpublished method, a short description). */
    record Reference(String algorithm, String citation) {}

    static List<Reference> all() {
        final List<Reference> refs = new ArrayList<>();
        refs.add( new Reference( "MAD rooting (Tools → MAD-Root)",
                "Tria FDK, Landan G, Dagan T (2017): \"Phylogenetic rooting using minimal ancestor deviation\", "
                        + "Nature Ecology & Evolution 1:0193." ) );
        refs.add( new Reference( "Midpoint rooting (Tools → Midpoint-Root)",
                "Farris JS (1972): \"Estimating Phylogenetic Trees from Distance Matrices\", "
                        + "The American Naturalist 106(951):645–668." ) );
        refs.add( new Reference( "Gene tree / species tree reconciliation (Analysis → GSDI, GSDIR)",
                "Zmasek CM, Eddy SR (2001): \"A simple algorithm to infer gene duplication and speciation events "
                        + "on a gene tree\", Bioinformatics 17(9):821–828." ) );
        refs.add( new Reference( "Representative-tip selection / sequence clustering (Tools → Select Representative Tips)",
                "Reduces redundancy by grouping tips that are close in patristic (tree) distance and keeping one "
                        + "representative per group, either within a distance cutoff or down to a target count. The "
                        + "representative kept is the group's medoid (most central) or longest-branch (most divergent) "
                        + "tip, and selected or found tips can be protected from removal (implemented in Archaeopteryx)." ) );
        refs.add( new Reference( "phyloXML (file format)", ForesterConstants.PHYLO_XML_REFERENCE ) );
        return refs;
    }

    /** The references as grouped, human-readable text for the dialog: each algorithm, then its citation. */
    static String asText() {
        final StringBuilder sb = new StringBuilder();
        sb.append( "Main citations for the algorithms implemented in Archaeopteryx.\n\n" );
        for ( final Reference r : all() ) {
            sb.append( r.algorithm() ).append( ":\n    " ).append( r.citation() ).append( "\n\n" );
        }
        return sb.toString().stripTrailing();
    }

    private AlgorithmReferences() {
    }
}
