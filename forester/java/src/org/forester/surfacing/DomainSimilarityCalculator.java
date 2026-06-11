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

package org.forester.surfacing;

import java.util.List;
import java.util.SortedSet;

interface DomainSimilarityCalculator {

    public SortedSet<DomainSimilarity> calculateSimilarities( final PairwiseDomainSimilarityCalculator pairwise_calculator,
                                                              final List<GenomeWideCombinableDomains> cdc_list,
                                                              final boolean ignore_domains_without_combinations_in_any_genome,
                                                              final boolean ignore_domains_specific_to_one_genome );;

                                                              public static enum Detailedness {
                                                                  BASIC, LIST_COMBINING_DOMAIN_FOR_EACH_SPECIES, PUNCTILIOUS
                                                              }

                                                              public static enum GoAnnotationOutput {
                                                                  ALL, NONE
                                                              }
}
