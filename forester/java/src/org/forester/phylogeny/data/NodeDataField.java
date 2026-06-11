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

package org.forester.phylogeny.data;

public enum NodeDataField {
    NODE_NAME,
    EVENT,
    SEQUENCE_NAME,
    GENE_NAME,
    SEQUENCE_SYMBOL,
    SEQUENCE_MOL_SEQ_FASTA,
    SEQUENCE_ACC,
    TAXONOMY_SCIENTIFIC_NAME,
    TAXONOMY_CODE,
    UNKNOWN,
    GO_TERM_IDS,
    SEQ_ANNOTATIONS,
    DOMAINS_ALL,
    DOMAINS_COLLAPSED_PER_PROTEIN;

    @Override
    public String toString() {
        switch ( this ) {
            case DOMAINS_ALL:
                return "Domains";
            case DOMAINS_COLLAPSED_PER_PROTEIN:
                return "Domains (collapsed per protein)";
            case EVENT:
                return "Events";
            case GENE_NAME:
                return "Gene Names";
            case GO_TERM_IDS:
                return "GO Term IDs";
            case NODE_NAME:
                return "Node Names";
            case SEQ_ANNOTATIONS:
                return "Sequence Annotations";
            case SEQUENCE_ACC:
                return "Sequence Accessors";
            case SEQUENCE_MOL_SEQ_FASTA:
                return "Molecular Sequences (Fasta)";
            case SEQUENCE_NAME:
                return "Sequence Names";
            case SEQUENCE_SYMBOL:
                return "Sequence Symbols";
            case TAXONOMY_CODE:
                return "Taxonomy Codes";
            case TAXONOMY_SCIENTIFIC_NAME:
                return "Scientific Names";
            case UNKNOWN:
                return "User Selected Data Fields";
            default:
                throw new IllegalArgumentException();
        }
    }
}