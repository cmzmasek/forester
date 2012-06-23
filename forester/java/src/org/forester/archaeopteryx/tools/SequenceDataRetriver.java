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
import org.forester.ws.seqdb.SequenceDatabaseEntry;
import org.forester.ws.seqdb.SequenceDbWsTools;

public final class SequenceDataRetriver extends RunnableProcess {

    public final static int            DEFAULT_LINES_TO_RETURN = 50;
    private final Phylogeny            _phy;
    private final MainFrameApplication _mf;
    private final TreePanel            _treepanel;
    private final static boolean       DEBUG                   = false;

    private enum Db {
        UNIPROT, EMBL, NCBI, NONE, REFSEQ;
    }

    public SequenceDataRetriver( final MainFrameApplication mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    private void execute() {
        start( _mf, "sequence data" );
        SortedSet<String> not_found = null;
        try {
            not_found = obtainSeqInformation( _phy, false, true );
        }
        catch ( final UnknownHostException e ) {
            final String what = "_"; //TODO FIXME 
            JOptionPane.showMessageDialog( _mf,
                                           "Could not connect to \"" + what + "\"",
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
                                               "Sequence tool successfully completed",
                                               "Sequence Tool Completed",
                                               JOptionPane.INFORMATION_MESSAGE );
            }
            catch ( final Exception e ) {
                // Not important if this fails, do nothing.
            }
        }
    }

    public static SortedSet<String> obtainSeqInformation( final Phylogeny phy,
                                                          final boolean ext_nodes_only,
                                                          final boolean allow_to_set_taxonomic_data )
            throws IOException {
        final SortedSet<String> not_found = new TreeSet<String>();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( ext_nodes_only && node.isInternal() ) {
                continue;
            }
            final Sequence seq = node.getNodeData().isHasSequence() ? node.getNodeData().getSequence() : new Sequence();
            final Taxonomy tax = node.getNodeData().isHasTaxonomy() ? node.getNodeData().getTaxonomy() : new Taxonomy();
            String query = null;
            Identifier id = null;
            Db db = Db.NONE;
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
                if ( ( query = SequenceDbWsTools.parseUniProtAccessor( node.getName() ) ) != null ) {
                    db = Db.UNIPROT;
                }
                else if ( ( id = SequenceIdParser.parse( node.getName() ) ) != null ) {
                    if ( id.getProvider().equalsIgnoreCase( Identifier.NCBI ) ) {
                        db = Db.NCBI;
                    }
                    else if ( id.getProvider().equalsIgnoreCase( Identifier.REFSEQ ) ) {
                        db = Db.REFSEQ;
                    }
                }
            }
            if ( db == Db.NONE ) {
                not_found.add( node.getName() );
            }
            SequenceDatabaseEntry db_entry = null;
            if ( !ForesterUtil.isEmpty( query ) ) {
                if ( db == Db.UNIPROT ) {
                    if ( DEBUG ) {
                        System.out.println( "uniprot: " + query );
                    }
                    db_entry = SequenceDbWsTools.obtainUniProtEntry( query, DEFAULT_LINES_TO_RETURN );
                }
                if ( ( db == Db.EMBL ) || ( ( db == Db.UNIPROT ) && ( db_entry == null ) ) ) {
                    if ( DEBUG ) {
                        System.out.println( "embl: " + query );
                    }
                    db_entry = SequenceDbWsTools.obtainEmblEntry( new Identifier( query ), DEFAULT_LINES_TO_RETURN );
                    if ( ( db == Db.UNIPROT ) && ( db_entry != null ) ) {
                        db = Db.EMBL;
                    }
                }
            }
            else if ( ( db == Db.REFSEQ ) && ( id != null ) ) {
                db_entry = SequenceDbWsTools.obtainRefSeqEntryFromEmbl( id, DEFAULT_LINES_TO_RETURN );
            }
            else if ( ( db == Db.NCBI ) && ( id != null ) ) {
                db_entry = SequenceDbWsTools.obtainEmblEntry( id, DEFAULT_LINES_TO_RETURN );
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
                    else if ( db == Db.NCBI ) {
                        type = "ncbi";
                    }
                    else if ( db == Db.REFSEQ ) {
                        type = "refseq";
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
                if ( allow_to_set_taxonomic_data && !ForesterUtil.isEmpty( db_entry.getTaxonomyIdentifier() ) ) {
                    tax.setIdentifier( new Identifier( db_entry.getTaxonomyIdentifier(), "uniprot" ) );
                }
                node.getNodeData().setTaxonomy( tax );
                node.getNodeData().setSequence( seq );
            }
            else if ( db != Db.NONE ) {
                not_found.add( node.getName() );
            }
            try {
                Thread.sleep( 10 );// Sleep for 10 ms
            }
            catch ( final InterruptedException ie ) {
            }
        }
        return not_found;
    }

    @Override
    public void run() {
        execute();
    }
}
