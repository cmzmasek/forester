// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: www.phylosoft.org/forester

package org.forester.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.UniProtTaxonomy;

public final class AncestralTaxonomyInference {

    public static void inferTaxonomyFromDescendents( final Phylogeny phy ) throws IOException,
            AncestralTaxonomyInferenceException {
        TaxonomyDataManager.clearCachesIfTooLarge();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( !node.isExternal() ) {
                inferTaxonomyFromDescendents( node );
            }
        }
    }

    private static void inferTaxonomyFromDescendents( final PhylogenyNode n ) throws IOException,
            AncestralTaxonomyInferenceException {
        if ( n.isExternal() ) {
            throw new IllegalArgumentException( "attempt to infer taxonomy from descendants of external node" );
        }
        n.getNodeData().setTaxonomy( null );
        final List<PhylogenyNode> descs = n.getDescendants();
        final List<String[]> lineages = new ArrayList<String[]>();
        int shortest_lin_length = Integer.MAX_VALUE;
        for( final PhylogenyNode desc : descs ) {
            if ( desc.getNodeData().isHasTaxonomy()
                    && ( TaxonomyDataManager.isHasAppropriateId( desc.getNodeData().getTaxonomy() )
                            || !ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getScientificName() )
                            || !ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getLineage() )
                            || !ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getTaxonomyCode() ) || !ForesterUtil
                            .isEmpty( desc.getNodeData().getTaxonomy().getCommonName() ) ) ) {
                final UniProtTaxonomy up_tax = TaxonomyDataManager.obtainUniProtTaxonomy( desc.getNodeData()
                        .getTaxonomy(), null, null );
                if ( ( up_tax == null ) && ForesterUtil.isEmpty( desc.getNodeData().getTaxonomy().getLineage() ) ) {
                    String desc_str = "";
                    if ( !ForesterUtil.isEmpty( desc.getName() ) ) {
                        desc_str = "\"" + desc.getName() + "\"";
                    }
                    else {
                        desc_str = "[" + desc.getId() + "]";
                    }
                    System.out.println( desc.getNodeData().getTaxonomy().toString() );
                    System.out.println( ForesterUtil.stringListToString( desc.getNodeData().getTaxonomy().getLineage(),
                                                                         "  >  " ) );
                    throw new AncestralTaxonomyInferenceException( "a taxonomy for node " + desc_str
                            + " could not be established from the database" );
                }
                String[] lineage = ForesterUtil.stringListToArray( desc.getNodeData().getTaxonomy().getLineage() );
                if ( ( lineage == null ) || ( lineage.length < 1 ) ) {
                    lineage = ForesterUtil.stringListToArray( up_tax.getLineage() );
                }
                if ( ( lineage == null ) || ( lineage.length < 1 ) ) {
                    throw new AncestralTaxonomyInferenceException( "a taxonomic lineage for node \""
                            + desc.getNodeData().getTaxonomy().toString() + "\" could not be established" );
                }
                if ( lineage.length < shortest_lin_length ) {
                    shortest_lin_length = lineage.length;
                }
                lineages.add( lineage );
            }
            else {
                String node = "";
                if ( !ForesterUtil.isEmpty( desc.getName() ) ) {
                    node = "\"" + desc.getName() + "\"";
                }
                else {
                    node = "[" + desc.getId() + "]";
                }
                throw new AncestralTaxonomyInferenceException( "node " + node
                        + " has no or inappropriate taxonomic information" );
            }
        }
        final List<String> last_common_lineage = new ArrayList<String>();
        String last_common = null;
        if ( shortest_lin_length > 0 ) {
            I: for( int i = 0; i < shortest_lin_length; ++i ) {
                final String lineage_0 = lineages.get( 0 )[ i ];
                for( int j = 1; j < lineages.size(); ++j ) {
                    if ( !lineage_0.equals( lineages.get( j )[ i ] ) ) {
                        break I;
                    }
                }
                last_common_lineage.add( lineage_0 );
                last_common = lineage_0;
            }
        }
        if ( last_common_lineage.isEmpty() ) {
            boolean saw_viruses = false;
            boolean saw_cellular_organism = false;
            for( final String[] lineage : lineages ) {
                if ( lineage.length > 0 ) {
                    if ( lineage[ 0 ].equalsIgnoreCase( UniProtTaxonomy.VIRUSES ) ) {
                        saw_viruses = true;
                    }
                    else if ( lineage[ 0 ].equalsIgnoreCase( UniProtTaxonomy.CELLULAR_ORGANISMS ) ) {
                        saw_cellular_organism = true;
                    }
                    if ( saw_cellular_organism && saw_viruses ) {
                        break;
                    }
                }
            }
            if ( saw_cellular_organism && saw_viruses ) {
                last_common_lineage.add( UniProtTaxonomy.CELLULAR_ORGANISMS );
                last_common = UniProtTaxonomy.CELLULAR_ORGANISMS;
            }
            else {
                String msg = "no common lineage for:\n";
                int counter = 0;
                for( final String[] strings : lineages ) {
                    msg += counter + ": ";
                    ++counter;
                    for( final String string : strings ) {
                        msg += string + " ";
                    }
                    msg += "\n";
                }
                throw new AncestralTaxonomyInferenceException( msg );
            }
        }
        final Taxonomy tax = new Taxonomy();
        n.getNodeData().setTaxonomy( tax );
        tax.setScientificName( last_common );
        final UniProtTaxonomy up_tax = TaxonomyDataManager.obtainUniProtTaxonomyFromLineage( last_common_lineage );
        if ( up_tax != null ) {
            if ( !ForesterUtil.isEmpty( up_tax.getRank() ) ) {
                try {
                    tax.setRank( up_tax.getRank().toLowerCase() );
                }
                catch ( final PhyloXmlDataFormatException ex ) {
                    tax.setRank( "" );
                }
            }
            if ( !ForesterUtil.isEmpty( up_tax.getId() ) ) {
                tax.setIdentifier( new Identifier( up_tax.getId(), "uniprot" ) );
            }
            if ( !ForesterUtil.isEmpty( up_tax.getCommonName() ) ) {
                tax.setCommonName( up_tax.getCommonName() );
            }
            if ( !ForesterUtil.isEmpty( up_tax.getSynonym() ) && !tax.getSynonyms().contains( up_tax.getSynonym() ) ) {
                tax.getSynonyms().add( up_tax.getSynonym() );
            }
            if ( up_tax.getLineage() != null ) {
                tax.setLineage( new ArrayList<String>() );
                for( final String lin : up_tax.getLineage() ) {
                    if ( !ForesterUtil.isEmpty( lin ) ) {
                        tax.getLineage().add( lin );
                    }
                }
            }
        }
        if ( ForesterUtil.isEmpty( tax.getLineage() ) ) {
            tax.setLineage( new ArrayList<String>() );
            for( final String lin : last_common_lineage ) {
                if ( !ForesterUtil.isEmpty( lin ) ) {
                    tax.getLineage().add( lin );
                }
            }
        }
        for( final PhylogenyNode desc : descs ) {
            if ( !desc.isExternal() && desc.getNodeData().isHasTaxonomy()
                    && desc.getNodeData().getTaxonomy().isEqual( tax ) ) {
                desc.getNodeData().setTaxonomy( null );
            }
        }
    }
}
