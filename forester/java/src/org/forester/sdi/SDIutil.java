
package org.forester.sdi;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class SDIutil {

    static String taxonomyToString( final PhylogenyNode n, final TaxonomyComparisonBase base ) {
        switch ( base ) {
            case ID:
                final Identifier id = n.getNodeData().getTaxonomy().getIdentifier();
                if ( id == null ) {
                    return null;
                }
                return id.getValuePlusProvider();
            case CODE:
                return n.getNodeData().getTaxonomy().getTaxonomyCode();
            case SCIENTIFIC_NAME:
                return n.getNodeData().getTaxonomy().getScientificName();
            default:
                throw new IllegalArgumentException( "unknown comparison base for taxonomies: " + base );
        }
    }

    public enum ALGORITHM {
        GSDIR, GSDI, SDI, SDIR
    }

    public enum TaxonomyComparisonBase {
        ID {

            @Override
            public String toString() {
                return "taxonomy id";
            }
        },
        CODE {

            @Override
            public String toString() {
                return "taxonomy code/mnemonic";
            }
        },
        SCIENTIFIC_NAME {

            @Override
            public String toString() {
                return "scientific name";
            }
        }
    }

    public final static TaxonomyComparisonBase determineTaxonomyComparisonBase( final Phylogeny gene_tree ) {
        int with_id_count = 0;
        int with_code_count = 0;
        int with_sn_count = 0;
        int max = 0;
        for( final PhylogenyNodeIterator iter = gene_tree.iteratorExternalForward(); iter.hasNext(); ) {
            final PhylogenyNode g = iter.next();
            if ( g.getNodeData().isHasTaxonomy() ) {
                final Taxonomy tax = g.getNodeData().getTaxonomy();
                if ( ( tax.getIdentifier() != null ) && !ForesterUtil.isEmpty( tax.getIdentifier().getValue() ) ) {
                    if ( ++with_id_count > max ) {
                        max = with_id_count;
                    }
                }
                if ( !ForesterUtil.isEmpty( tax.getTaxonomyCode() ) ) {
                    if ( ++with_code_count > max ) {
                        max = with_code_count;
                    }
                }
                if ( !ForesterUtil.isEmpty( tax.getScientificName() ) ) {
                    if ( ++with_sn_count > max ) {
                        max = with_sn_count;
                    }
                }
            }
        }
        if ( max == 0 ) {
            throw new IllegalArgumentException( "gene tree has no taxonomic data" );
        }
        else if ( max == 1 ) {
            throw new IllegalArgumentException( "gene tree has only one node with taxonomic data" );
        }
        else if ( max == with_id_count ) {
            return TaxonomyComparisonBase.ID;
        }
        else if ( max == with_sn_count ) {
            return TaxonomyComparisonBase.SCIENTIFIC_NAME;
        }
        else {
            return TaxonomyComparisonBase.CODE;
        }
    }
}
