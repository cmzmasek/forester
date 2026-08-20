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
        refs.add( new Reference( "Tanglegram auto-untangle / crossing minimisation (Analysis → Create Tanglegram → Auto-untangle)",
                "Reorders the two trees of a tanglegram by rotating clades (reversing the child order at internal "
                        + "nodes -- a topology-preserving operation) to reduce the number of crossing connectors "
                        + "between their matched tips. Archaeopteryx uses a barycentre heuristic -- each clade's "
                        + "children are ordered by the mean vertical position of the tips they link to in the other "
                        + "tree -- applied alternately to both trees and iterated to convergence, combined with random "
                        + "restarts, keeping the arrangement with the fewest crossings (it never increases the "
                        + "crossing count). Minimising the crossings of a tanglegram is NP-hard, so this is a heuristic "
                        + "(implemented in Archaeopteryx). Barycentre crossing-reduction heuristic: Sugiyama K, Tagawa "
                        + "S, Toda M (1981): \"Methods for visual understanding of hierarchical system structures\", "
                        + "IEEE Transactions on Systems, Man, and Cybernetics 11(2):109–125. Tanglegram layout and "
                        + "crossing minimisation: Scornavacca C, Zickmann F, Huson DH (2011): \"Tanglegrams for rooted "
                        + "phylogenetic trees and networks\", Bioinformatics 27(13):i248–i256." ) );
        refs.add( new Reference( "Tanglegram entanglement / congruence score (Analysis → Create Tanglegram)",
                "A size-normalised measure, in [0,1], of how tangled a tanglegram's tip-to-tip connectors are: 0 means "
                        + "the two leaf orderings agree perfectly (no crossings), 1 means they are fully reversed. Two "
                        + "connectors cross exactly when their matched tips are in the opposite vertical order in the two "
                        + "trees, so the total number of crossings is the number of DISCORDANT tip pairs — an inversion "
                        + "count, equal to the Kendall rank-correlation (tau) distance between the two leaf orderings. "
                        + "Archaeopteryx divides this crossing count by the maximum possible number of pairs, n(n−1)/2, so "
                        + "the score is comparable across trees of different sizes, and computes it in O(n log n) via a "
                        + "counting merge sort (implemented in Archaeopteryx). Rank-correlation (concordant/discordant "
                        + "pairs) basis: Kendall MG (1938): \"A New Measure of Rank Correlation\", Biometrika 30(1–2):81–93. "
                        + "Tanglegram crossings as a measure of tree congruence: Scornavacca C, Zickmann F, Huson DH "
                        + "(2011): \"Tanglegrams for rooted phylogenetic trees and networks\", Bioinformatics 27(13):i248–"
                        + "i256. A related, leaf-position-based entanglement measure is implemented in the dendextend R "
                        + "package: Galili T (2015): \"dendextend: an R package for visualizing, adjusting and comparing "
                        + "trees of hierarchical clustering\", Bioinformatics 31(22):3718–3720." ) );
        refs.add( new Reference( "Geologic time axis (Settings → Overlays → Time Axis → Geologic)",
                "Draws the international geologic time scale beneath a dated (time-calibrated) tree as coloured, "
                        + "named interval bands (System/Period over Series/Epoch), so a node's position along the "
                        + "branch-length (= time) axis can be read directly against the named geologic intervals. The "
                        + "interval boundaries, names, and official colours are those of the International "
                        + "Chronostratigraphic Chart of the International Commission on Stratigraphy (ICS / IUGS), "
                        + "www.stratigraphy.org. " + GeologicTimeScale.REFERENCE ) );
        refs.add( new Reference( "Fossil range bars (Settings → Overlays → Data Overlays → Fossil Range Bars)",
                "On a dated (time-calibrated) tree, draws a stratigraphic-range bar at each fossil tip spanning its "
                        + "observed duration -- from its First Appearance Datum (FAD, oldest occurrence) to its Last "
                        + "Appearance Datum (LAD, youngest occurrence), read from the tip's phyloXML date min/max. The "
                        + "convention of plotting taxon stratigraphic ranges against a time-calibrated phylogeny follows "
                        + "the strap R package: Bell MA, Lloyd GT (2015): \"strap: an R package for plotting phylogenies "
                        + "against stratigraphy and assessing their stratigraphic congruence\", Palaeontology 58(2):379–389." ) );
        refs.add( new Reference( "Node age bars / spindles (HPD) (Settings → Overlays → Data Overlays)",
                "On a dated (time-calibrated) phylogram, draws each internal node's divergence-time uncertainty from its "
                        + "phyloXML date -- the point estimate (median/mean height) and the 95% Highest Posterior Density "
                        + "(HPD) interval (min/max), as produced by Bayesian MCMC dating (e.g. BEAST / TreeAnnotator "
                        + "height_95%_HPD). \"Node age shape\" selects the rendering: a flat BAR spanning the interval (the "
                        + "FigTree \"node bars\" convention), or a SPINDLE -- a tapered lens peaking at the point estimate "
                        + "and narrowing to zero at the HPD bounds, so the position of the estimate WITHIN its interval "
                        + "(the skew) is visible. NOTE: the spindle is a SCHEMATIC of the summarised uncertainty (the point "
                        + "estimate + the 95% HPD), NOT the raw posterior density -- a summary (MCC) tree does not carry "
                        + "the per-node posterior sample. Source of the dated ages and HPD intervals: Suchard MA, Lemey P, "
                        + "Baele G, Ayres DL, Drummond AJ, Rambaut A (2018): \"Bayesian phylogenetic and phylodynamic data "
                        + "integration using BEAST 1.10\", Virus Evolution 4(1):vey016. The node-bars visualisation "
                        + "convention is that of FigTree (A. Rambaut, http://tree.bio.ed.ac.uk/software/figtree/)." ) );
        refs.add( new Reference( "Auspice / Nextstrain JSON import (File → Read Tree from File → .json)",
                "Reads an Auspice / Nextstrain v2 dataset (dataset.json) -- the interchange format for dated, "
                        + "annotated pathogen phylogenies -- into Archaeopteryx's native model: node dates (num_date) and "
                        + "their confidence intervals become the calendar time axis + node-age spindles, the cumulative "
                        + "divergence (div) is kept for a time/divergence view, and each discrete trait (country, region, "
                        + "clade_membership, host, ...) becomes a colourable/searchable node property whose per-node "
                        + "confidence drives the ancestral-state pies. The map, entropy, and frequencies panels are not "
                        + "imported (Archaeopteryx is a tree viewer). Nextstrain / Auspice: Hadfield J, Megill C, Bell SM, "
                        + "Huddleston J, Potter B, Callender C, Sagulenko P, Bedford T, Neher RA (2018): \"Nextstrain: "
                        + "real-time tracking of pathogen evolution\", Bioinformatics 34(23):4121–4123." ) );
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
