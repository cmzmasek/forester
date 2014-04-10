// $Id:
// forester -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2010 Christian M. Zmasek
// Copyright (C) 2008-2010 Burnham Institute for Medical Research
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx.webservices;

import java.util.ArrayList;
import java.util.List;

import org.forester.archaeopteryx.webservices.WebservicesManager.WsPhylogenyFormat;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.phylogeny.iterators.PreorderTreeIterator;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

public final class WebserviceUtil {

    public static final String PFAM_INST                       = "pfam";
    public static final String PFAM_NAME                       = "Pfam";
    public static final String PFAM_SERVER                     = "http://pfam.janelia.org";
    public static final String TOL_NAME                        = "Tree of Life (ToL)";
    public static final String TOL_URL_BASE                    = "http://tolweb.org/onlinecontributors/app?service=external&page=xml/TreeStructureService&node_id=";
    public static final String TOL_WEBSERVER                   = TOL_URL_BASE
                                                                       + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER;
    public static final String TREE_BASE_DESC                  = "This data set was downloaded from TreeBASE, a relational database of phylogenetic knowledge. TreeBASE has been supported by the NSF, Harvard University, Yale University, SDSC and UC Davis. Please do not remove this acknowledgment.";
    public static final String TREE_BASE_INST                  = "treebase";
    public static final String TREE_BASE_NAME                  = "TreeBASE";
    public static final String TREE_FAM_INST                   = "tree_fam";
    public static final String TREE_FAM_NAME                   = "TreeFam";
    public static final String TREE_FAM_URL_BASE               = "http://www.treefam.org/family/TF";
    public static final String TREEBASE_PHYLOWS_STUDY_URL_BASE = "http://purl.org/phylo/treebase/phylows/study/TB2:S";
    public static final String TREEBASE_PHYLOWS_TREE_URL_BASE  = "http://purl.org/phylo/treebase/phylows/tree/TB2:Tr";

    public static List<PhylogeniesWebserviceClient> createDefaultClients() {
        final List<PhylogeniesWebserviceClient> clients = new ArrayList<PhylogeniesWebserviceClient>();
        clients.add( new BasicPhylogeniesWebserviceClient( TREE_BASE_NAME,
                                                           "Read Tree(s) from TreeBASE Study...",
                                                           "Use TreeBASE to obtain evolutionary tree(s) from a study",
                                                           "Please enter a TreeBASE study (\"S\") identifier (without the \"S\")\n(Examples: 14909, 14525, 15613, 15632)",
                                                           WsPhylogenyFormat.TREEBASE_STUDY,
                                                           null,
                                                           TREEBASE_PHYLOWS_STUDY_URL_BASE
                                                                   + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER
                                                                   + "?format=nexus",
                                                           true,
                                                           "http://www.treebase.org",
                                                           TREE_BASE_INST ) );
        clients.add( new BasicPhylogeniesWebserviceClient( TREE_BASE_NAME,
                                                           "Read Tree from TreeBASE...",
                                                           "Use TreeBASE to obtain a evolutionary tree",
                                                           "Please enter a TreeBASE tree (\"Tr\") identifier (without the \"Tr\")\n(Examples: 2406, 422, 2654, 825, 4931, 2518, 4934)",
                                                           WsPhylogenyFormat.TREEBASE_TREE,
                                                           null,
                                                           TREEBASE_PHYLOWS_TREE_URL_BASE
                                                                   + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER
                                                                   + "?format=nexus",
                                                           true,
                                                           "http://www.treebase.org",
                                                           TREE_BASE_INST ) );
        clients.add( new BasicPhylogeniesWebserviceClient( PFAM_NAME,
                                                           "Read Gene Tree from Pfam...",
                                                           "Use  Pfam to obtain gene trees for seed alignments",
                                                           "Please enter a Pfam (PF) accession number\n(Examples: 01849 for NAC, 00452 for Bcl-2, 00046 for Homeobox)",
                                                           WsPhylogenyFormat.PFAM,
                                                           null,
                                                           PFAM_SERVER + "/family/PF"
                                                                   + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER
                                                                   + "/tree/download",
                                                           false,
                                                           PFAM_SERVER,
                                                           PFAM_INST ) );
        clients.add( new BasicPhylogeniesWebserviceClient( TREE_FAM_NAME,
                                                           "Read Gene Tree from TreeFam...",
                                                           "Use TreeFam to obtain a gene tree",
                                                           "Please enter a TreeFam (TF) accession number\n(Examples: 101004 for Cyclin D, 315938 for Hox, 105310 for Wnt)",
                                                           WsPhylogenyFormat.NHX,
                                                           null,
                                                           TREE_FAM_URL_BASE
                                                                   + PhylogeniesWebserviceClient.QUERY_PLACEHOLDER
                                                                   + "/tree/newick",
                                                           true,
                                                           "http://www.treefam.org",
                                                           TREE_FAM_INST ) );
        clients.add( new BasicPhylogeniesWebserviceClient( TOL_NAME,
                                                           "Read Tree from Tree of Life (ToL)...",
                                                           "Use ToL webservice to obtain a evolutionary tree",
                                                           "Please enter a Tree of Life node identifier\n(Examples: "
                                                                   + "14923 for ray-finned fishes, 19386 for Cephalopoda, 2461 for Cnidaria)",
                                                           WsPhylogenyFormat.TOL_XML_RESPONSE,
                                                           PhylogenyMethods.PhylogenyNodeField.TAXONOMY_SCIENTIFIC_NAME,
                                                           WebserviceUtil.TOL_WEBSERVER,
                                                           true,
                                                           "http://tolweb.org",
                                                           null ) );
        return clients;
    }

    public static void processInstructions( final PhylogeniesWebserviceClient client, final Phylogeny phylogeny )
            throws PhyloXmlDataFormatException {
        if ( client.getProcessingInstructions().equals( WebserviceUtil.TREE_FAM_INST ) ) {
            WebserviceUtil.processTreeFamTrees( phylogeny );
        }
        else if ( client.getProcessingInstructions().equals( WebserviceUtil.PFAM_INST ) ) {
            WebserviceUtil.extractSpTremblAccFromNodeName( phylogeny, "sptrembl" );
            PhylogenyMethods.transferInternalNodeNamesToConfidence( phylogeny, "bootstrap" );
        }
        else if ( client.getProcessingInstructions().equals( WebserviceUtil.TREE_BASE_INST ) ) {
            if ( PhylogenyMethods.isInternalNamesLookLikeConfidences( phylogeny ) ) {
                PhylogenyMethods.transferInternalNodeNamesToConfidence( phylogeny, "" );
            }
            WebserviceUtil.processTreeBaseTrees( phylogeny );
        }
    }

    static void extractSpTremblAccFromNodeName( final Phylogeny phy, final String source ) {
        final PreorderTreeIterator it = new PreorderTreeIterator( phy );
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                final String name = n.getName();
                final int i = name.lastIndexOf( "/" );
                if ( i > 0 ) {
                    final String acc_str = name.substring( 0, i );
                    if ( !ForesterUtil.isEmpty( acc_str ) ) {
                        final Sequence seq = new Sequence();
                        final Accession acc = new Accession( acc_str, source );
                        seq.setAccession( acc );
                        n.getNodeData().setSequence( seq );
                    }
                }
            }
        }
    }

    static void processTreeBaseTrees( final Phylogeny phy ) {
        phy.setDescription( TREE_BASE_DESC );
        final PhylogenyNodeIterator it = phy.iteratorExternalForward();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                final Accession acc = SequenceAccessionTools.parseAccessorFromString( n.getName() );
                if ( acc != null ) {
                    if ( !n.getNodeData().isHasSequence() ) {
                        n.getNodeData().addSequence( new Sequence() );
                    }
                    final Sequence s = n.getNodeData().getSequence();
                    if ( s.getAccession() == null ) {
                        s.setAccession( acc );
                    }
                }
            }
        }
    }

    static void processTreeFamTrees( final Phylogeny phy ) {
        final PhylogenyNodeIterator it = phy.iteratorPostorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.isExternal() ) {
                n.getNodeData().setEvent( null );
                if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                    final Accession acc = SequenceAccessionTools.parseAccessorFromString( n.getName() );
                    if ( acc != null ) {
                        if ( !n.getNodeData().isHasSequence() ) {
                            n.getNodeData().addSequence( new Sequence() );
                        }
                        final Sequence s = n.getNodeData().getSequence();
                        if ( s.getAccession() == null ) {
                            s.setAccession( acc );
                        }
                    }
                }
            }
            else {
                if ( ( n.getBranchData() != null ) && n.getBranchData().isHasConfidences()
                        && ( n.getBranchData().getConfidence( 0 ) != null ) ) {
                    n.getBranchData().getConfidence( 0 ).setType( "bootstrap" );
                }
                if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                    if ( !n.getNodeData().isHasTaxonomy() ) {
                        n.getNodeData().addTaxonomy( new Taxonomy() );
                    }
                    final Taxonomy t = n.getNodeData().getTaxonomy();
                    if ( ForesterUtil.isEmpty( t.getScientificName() ) ) {
                        t.setScientificName( n.getName() );
                        n.setName( "" );
                    }
                }
            }
            if ( n.getNodeData().isHasTaxonomy() && ( n.getNodeData().getTaxonomy().getIdentifier() != null ) ) {
                n.getNodeData()
                        .getTaxonomy()
                        .setIdentifier( new Identifier( n.getNodeData().getTaxonomy().getIdentifier().getValue(),
                                                        "ncbi" ) );
            }
        }
    }
}
