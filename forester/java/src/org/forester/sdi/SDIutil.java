
package org.forester.sdi;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class SDIutil {

    public final static TaxonomyComparisonBase determineTaxonomyComparisonBase( final Phylogeny gene_tree )
            throws SDIException {
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
            throw new SDIException( "gene tree has no taxonomic data" );
        }
        else if ( max == 1 ) {
            throw new SDIException( "gene tree has only one node with taxonomic data" );
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

    public final static Phylogeny parseSpeciesTree( final Phylogeny gene_tree,
                                                    final File species_tree_file,
                                                    final boolean replace_undescores_in_nhx_trees,
                                                    final boolean ignore_quotes_in_nhx_trees,
                                                    final TAXONOMY_EXTRACTION taxonomy_extraction_in_nhx_trees )
                                                            throws FileNotFoundException, PhyloXmlDataFormatException, IOException, SDIException {
        Phylogeny species_tree;
        final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
        final PhylogenyParser p = ParserUtils.createParserDependingOnFileType( species_tree_file, true );
        if ( p instanceof PhyloXmlParser ) {
            species_tree = factory.create( species_tree_file, p )[ 0 ];
        }
        else {
            if ( p instanceof NHXParser ) {
                final NHXParser nhx = ( NHXParser ) p;
                nhx.setReplaceUnderscores( replace_undescores_in_nhx_trees );
                nhx.setIgnoreQuotes( ignore_quotes_in_nhx_trees );
                nhx.setTaxonomyExtraction( taxonomy_extraction_in_nhx_trees );
            }
            else if ( p instanceof NexusPhylogeniesParser ) {
                final NexusPhylogeniesParser nex = ( NexusPhylogeniesParser ) p;
                nex.setReplaceUnderscores( replace_undescores_in_nhx_trees );
                nex.setIgnoreQuotes( ignore_quotes_in_nhx_trees );
                nex.setTaxonomyExtraction( taxonomy_extraction_in_nhx_trees );
            }
            species_tree = factory.create( species_tree_file, p )[ 0 ];
            species_tree.setRooted( true );
            final TaxonomyComparisonBase comp_base = determineTaxonomyComparisonBase( gene_tree );
            switch ( comp_base ) {
                case SCIENTIFIC_NAME:
                    PhylogenyMethods
                    .transferNodeNameToField( species_tree,
                                              PhylogenyMethods.PhylogenyNodeField.TAXONOMY_SCIENTIFIC_NAME,
                                              true );
                    break;
                case CODE:
                    PhylogenyMethods.transferNodeNameToField( species_tree,
                                                              PhylogenyMethods.PhylogenyNodeField.TAXONOMY_CODE,
                                                              true );
                    break;
                case ID:
                    PhylogenyMethods.transferNodeNameToField( species_tree,
                                                              PhylogenyMethods.PhylogenyNodeField.TAXONOMY_ID,
                                                              true );
                    break;
                default:
                    throw new SDIException( "unable to determine comparison base" );
            }
        }
        return species_tree;
    }

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
}
