// Exp $
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

package org.forester.archaeopteryx.tools;

import java.io.IOException;
import java.net.UnknownHostException;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.swing.JOptionPane;

import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.ws.uniprot.UniProtEntry;
import org.forester.ws.uniprot.UniProtWsTools;

public class UniProtSequenceObtainer implements Runnable {

    private final Phylogeny            _phy;
    private final MainFrameApplication _mf;
    private final TreePanel            _treepanel;

    public UniProtSequenceObtainer( final MainFrameApplication mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    private String getBaseUrl() {
        return UniProtWsTools.BASE_URL;
    }

    private void execute() {
        _mf.getMainPanel().getCurrentTreePanel().setWaitCursor();
        SortedSet<String> not_found = null;
        try {
            not_found = obtainSeqInformation( _phy );
        }
        catch ( final UnknownHostException e ) {
            _mf.getMainPanel().getCurrentTreePanel().setArrowCursor();
            JOptionPane.showMessageDialog( _mf,
                                           "Could not connect to \"" + getBaseUrl() + "\"",
                                           "Network error during taxonomic information gathering",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final IOException e ) {
            _mf.getMainPanel().getCurrentTreePanel().setArrowCursor();
            e.printStackTrace();
            JOptionPane.showMessageDialog( _mf,
                                           e.toString(),
                                           "Failed to obtain taxonomic information",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        finally {
            _mf.getMainPanel().getCurrentTreePanel().setArrowCursor();
        }
        _treepanel.setTree( _phy );
        _mf.showWhole();
        _treepanel.setEdited( true );
        if ( ( not_found != null ) && ( not_found.size() > 0 ) ) {
            int max = not_found.size();
            boolean more = false;
            if ( max > 20 ) {
                more = true;
                max = 20;
            }
            final StringBuffer sb = new StringBuffer();
            sb.append( "Not all identifiers could be resolved.\n" );
            if ( not_found.size() == 1 ) {
                sb.append( "The following identifier was not found:\n" );
            }
            else {
                sb.append( "The following identifiers were not found (total: " + not_found.size() + "):\n" );
            }
            int i = 0;
            for( final String string : not_found ) {
                if ( i > 19 ) {
                    break;
                }
                sb.append( string );
                sb.append( "\n" );
                ++i;
            }
            if ( more ) {
                sb.append( "..." );
            }
            try {
                JOptionPane.showMessageDialog( _mf,
                                               sb.toString(),
                                               "UniProt Sequence Tool Completed",
                                               JOptionPane.WARNING_MESSAGE );
            }
            catch ( final Exception e ) {
                // Not important if this fails, do nothing. 
            }
        }
        else {
            try {
                JOptionPane.showMessageDialog( _mf,
                                               "UniProt sequence tool successfully completed",
                                               "UniProt Sequence Tool Completed",
                                               JOptionPane.INFORMATION_MESSAGE );
            }
            catch ( final Exception e ) {
                // Not important if this fails, do nothing.
            }
        }
    }

    synchronized public static SortedSet<String> obtainSeqInformation( final Phylogeny phy ) throws IOException {
        final SortedSet<String> not_found = new TreeSet<String>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.getNodeData().isHasSequence() ) {
                //TODO  Do something
            }
            //  else if ( !ForesterUtil.isEmpty( node.getName() ) ) {
            //     not_found.add( node.getName() );
            //  }
            else if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                String query = node.getName();
                if ( query.indexOf( '/' ) > 0 ) {
                    query = query.substring( 0, query.indexOf( '/' ) );
                }
                if ( query.indexOf( '.' ) > 0 ) {
                    query = query.substring( 0, query.indexOf( '.' ) );
                }
                if ( query.indexOf( '_' ) > 0 ) {
                    query = query.substring( 0, query.indexOf( '_' ) );
                }
                final UniProtEntry upe = obtainUniProtEntry( query );
                if ( upe != null ) {
                    final Sequence seq = new Sequence();
                    final Taxonomy tax = new Taxonomy();
                    if ( !ForesterUtil.isEmpty( upe.getAc() ) ) {
                        seq.setAccession( new Accession( upe.getAc(), "uniprot" ) );
                    }
                    if ( !ForesterUtil.isEmpty( upe.getRecName() ) ) {
                        seq.setName( upe.getRecName() );
                    }
                    if ( !ForesterUtil.isEmpty( upe.getSymbol() ) ) {
                        seq.setSymbol( upe.getSymbol() );
                    }
                    if ( !ForesterUtil.isEmpty( upe.getOsScientificName() ) ) {
                        tax.setScientificName( upe.getOsScientificName() );
                    }
                    if ( !ForesterUtil.isEmpty( upe.getTaxId() ) ) {
                        tax.setIdentifier( new Identifier( upe.getTaxId(), "uniprot" ) );
                    }
                    node.getNodeData().setTaxonomy( tax );
                    node.getNodeData().setSequence( seq );
                }
                else {
                    not_found.add( node.getName() );
                }
                //}
            }
        }
        return not_found;
    }

    static UniProtEntry obtainUniProtEntry( final String query ) throws IOException {
        final List<String> lines = UniProtWsTools.queryUniprot( "uniprot/" + query + ".txt", 200 );
        return UniProtEntry.createInstanceFromPlainText( lines );
    }

    @Override
    public void run() {
        execute();
    }
}
