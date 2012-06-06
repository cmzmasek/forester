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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.UnknownHostException;
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
import org.forester.util.SequenceIdParser;
import org.forester.ws.uniprot.DatabaseTools;
import org.forester.ws.uniprot.SequenceDatabaseEntry;
import org.forester.ws.uniprot.UniProtWsTools;

public final class SequenceDataRetriver extends RunnableProcess {

    private final Phylogeny            _phy;
    private final MainFrameApplication _mf;
    private final TreePanel            _treepanel;
    private final static boolean       DEBUG = false;

    private enum Db {
        UNKNOWN, UNIPROT, EMBL, NCBI;
    }

    public SequenceDataRetriver( final MainFrameApplication mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    private String getBaseUrl() {
        return UniProtWsTools.BASE_URL;
    }

    private void execute() {
        start( _mf, "sequence data" );
        SortedSet<String> not_found = null;
        try {
            not_found = obtainSeqInformation( _phy );
        }
        catch ( final UnknownHostException e ) {
            JOptionPane.showMessageDialog( _mf,
                                           "Could not connect to \"" + getBaseUrl() + "\"",
                                           "Network error during taxonomic information gathering",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            JOptionPane.showMessageDialog( _mf,
                                           e.toString(),
                                           "Failed to obtain taxonomic information",
                                           JOptionPane.ERROR_MESSAGE );
            return;
        }
        finally {
            end( _mf );
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
            if ( not_found.size() == 1 ) {
                sb.append( "Data for the following sequence identifier was not found:\n" );
            }
            else {
                sb.append( "Data for the following sequence identifiers was not found (total: " + not_found.size()
                        + "):\n" );
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
                                               "Sequence Tool Completed",
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

    public static SortedSet<String> obtainSeqInformation( final Phylogeny phy ) throws IOException {
        final SortedSet<String> not_found = new TreeSet<String>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            Sequence seq = null;
            Taxonomy tax = null;
            if ( node.getNodeData().isHasSequence() ) {
                seq = node.getNodeData().getSequence();
            }
            else {
                seq = new Sequence();
            }
            if ( node.getNodeData().isHasTaxonomy() ) {
                tax = node.getNodeData().getTaxonomy();
            }
            else {
                tax = new Taxonomy();
            }
            String query = null;
            Db db = Db.UNKNOWN;
            if ( node.getNodeData().isHasSequence() && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() )
                    && node.getNodeData().getSequence().getAccession().getValue().toLowerCase().startsWith( "uniprot" ) ) {
                query = node.getNodeData().getSequence().getAccession().getValue();
                db = Db.UNIPROT;
            }
            else if ( node.getNodeData().isHasSequence()
                    && ( node.getNodeData().getSequence().getAccession() != null )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getSource() )
                    && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() )
                    && ( node.getNodeData().getSequence().getAccession().getValue().toLowerCase().startsWith( "embl" ) || node
                            .getNodeData().getSequence().getAccession().getValue().toLowerCase().startsWith( "ebi" ) ) ) {
                query = node.getNodeData().getSequence().getAccession().getValue();
                db = Db.EMBL;
            }
            else if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                if ( ( query = UniProtWsTools.parseUniProtAccessor( node.getName() ) ) != null ) {
                    db = Db.UNIPROT;
                }
                else if ( ( query = SequenceIdParser.parseGenbankAccessor( node.getName() ) ) != null ) {
                    db = Db.NCBI;
                }
            }
            if ( !ForesterUtil.isEmpty( query ) ) {
                SequenceDatabaseEntry db_entry = null;
                if ( db == Db.UNIPROT ) {
                    if ( DEBUG ) {
                        System.out.println( "uniprot: " + query );
                    }
                    try {
                        db_entry = UniProtWsTools.obtainUniProtEntry( query, 200 );
                    }
                    catch ( final FileNotFoundException e ) {
                        // Ignore.
                    }
                }
                if ( ( db == Db.EMBL ) || ( ( db == Db.UNIPROT ) && ( db_entry == null ) ) ) {
                    if ( DEBUG ) {
                        System.out.println( "embl: " + query );
                    }
                    try {
                        db_entry = UniProtWsTools.obtainEmblEntry( query, 200 );
                    }
                    catch ( final FileNotFoundException e ) {
                        // Ignore.
                    }
                    if ( ( db == Db.UNIPROT ) && ( db_entry != null ) ) {
                        db = Db.EMBL;
                    }
                }
                if ( ( db_entry != null ) && !db_entry.isEmpty() ) {
                    if ( !ForesterUtil.isEmpty( db_entry.getAccession() ) ) {
                        String type = null;
                        if ( db == Db.EMBL ) {
                            type = "embl";
                        }
                        else if ( db == Db.UNIPROT ) {
                            type = "uniprot";
                        }
                        seq.setAccession( new Accession( db_entry.getAccession(), type ) );
                    }
                    if ( !ForesterUtil.isEmpty( db_entry.getSequenceName() ) ) {
                        seq.setName( db_entry.getSequenceName() );
                    }
                    if ( !ForesterUtil.isEmpty( db_entry.getSequenceSymbol() ) ) {
                        seq.setSymbol( db_entry.getSequenceSymbol() );
                    }
                    if ( !ForesterUtil.isEmpty( db_entry.getTaxonomyScientificName() ) ) {
                        tax.setScientificName( db_entry.getTaxonomyScientificName() );
                    }
                    if ( !ForesterUtil.isEmpty( db_entry.getTaxonomyIdentifier() ) ) {
                        tax.setIdentifier( new Identifier( db_entry.getTaxonomyIdentifier(), "uniprot" ) );
                    }
                    node.getNodeData().setTaxonomy( tax );
                    node.getNodeData().setSequence( seq );
                }
                else {
                    not_found.add( node.getName() );
                }
            }
        }
        return not_found;
    }

    @Override
    public void run() {
        execute();
    }
}
