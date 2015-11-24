
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