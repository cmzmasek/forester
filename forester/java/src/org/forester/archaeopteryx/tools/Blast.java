// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx.tools;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import org.forester.archaeopteryx.AptxUtil;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

public final class Blast {

    final public static void openNcbiBlastWeb( final String query,
                                               final boolean is_nucleic_acids,
                                             
                                               final TreePanel p ) {
        //http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&PAGE=Proteins&DATABASE=swissprot&QUERY=gi|163848401
        final StringBuilder uri_str = new StringBuilder();
        uri_str.append( "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&DATABASE=nr&PAGE=" );
        if ( is_nucleic_acids ) {
            uri_str.append( "Nucleotide" );
        }
        else {
            uri_str.append( "Proteins" );
        }
        uri_str.append( "&QUERY=" );
        uri_str.append( query );
        try {
            AptxUtil.launchWebBrowser( new URI( uri_str.toString() ), "_aptx_blast" );
        }
        catch ( final IOException e ) {
            AptxUtil.showErrorMessage( p, e.toString() );
            e.printStackTrace();
        }
        catch ( final URISyntaxException e ) {
            AptxUtil.showErrorMessage( p, e.toString() );
            e.printStackTrace();
        }
    }

    final public static String obtainQueryForBlast( final PhylogenyNode node ) {
        String query = "";
        if ( node.getNodeData().isHasSequence() ) {
            if ( !ForesterUtil.isEmpty( node.getNodeData().getSequence().getMolecularSequence() ) ) {
                query = node.getNodeData().getSequence().getMolecularSequence();
            }
            if ( ForesterUtil.isEmpty( query ) && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() ) ) {
                final Accession id = SequenceAccessionTools.parseAccessorFromString( node.getNodeData().getSequence()
                                                                                     .getAccession().getValue() );
                if ( id != null ) {
                    query = id.getValue();
                }
            }
            if ( ForesterUtil.isEmpty( query ) && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getName() ) ) {
                final Accession id = SequenceAccessionTools.parseAccessorFromString( node.getNodeData().getSequence()
                                                                                     .getName() );
                if ( id != null ) {
                    query = id.getValue();
                }
            }
            if ( ForesterUtil.isEmpty( query ) && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getSymbol() ) ) {
                final Accession id = SequenceAccessionTools.parseAccessorFromString( node.getNodeData().getSequence()
                                                                                     .getSymbol() );
                if ( id != null ) {
                    query = id.getValue();
                }
            }
            if ( ForesterUtil.isEmpty( query )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getGeneName() ) ) {
                final Accession id = SequenceAccessionTools.parseAccessorFromString( node.getNodeData().getSequence()
                                                                                     .getGeneName() );
                if ( id != null ) {
                    query = id.getValue();
                }
            }
        }
        if ( ForesterUtil.isEmpty( query ) && !ForesterUtil.isEmpty( node.getName() ) ) {
            final Accession id = SequenceAccessionTools.parseAccessorFromString( node.getName() );
            if ( id != null ) {
                query = id.getValue();
            }
        }
        return query;
    }

    final public static boolean isContainsQueryForBlast( final PhylogenyNode node ) {
        return !ForesterUtil.isEmpty( obtainQueryForBlast( node ) );
    }
}
