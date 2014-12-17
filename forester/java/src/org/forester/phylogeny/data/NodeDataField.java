
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

    
    public String toString() {
        switch ( this ) {
            case DOMAINS_ALL:
                return "Domain";
            case DOMAINS_COLLAPSED_PER_PROTEIN:
                return "Domain (collapsed per protein)";
            case EVENT:
                return "Event";
            case GENE_NAME:
                return "Gene Name";
            case GO_TERM_IDS:
                return "GO Term ID";
            case NODE_NAME:
                return "Node Name";
            case SEQ_ANNOTATIONS:
                return "Sequence Annotation";
            case SEQUENCE_ACC:
                return "Sequence Accessor";
            case SEQUENCE_MOL_SEQ_FASTA:
                return "Molecular Sequence (Fasta)";
            case SEQUENCE_NAME:
                return "Sequence Name";
            case SEQUENCE_SYMBOL:
                return "Sequence Symbol";
            case TAXONOMY_CODE:
                return "Taxonomy Code";
            case TAXONOMY_SCIENTIFIC_NAME:
                return "Scientific Name";
            case UNKNOWN:
                return "User/UI Selected Data Field(s)";
            default:
                throw new IllegalArgumentException();
        }
    }
}