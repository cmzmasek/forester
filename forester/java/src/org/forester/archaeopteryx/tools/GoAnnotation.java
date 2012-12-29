// $Id:
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.archaeopteryx.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Arrays;
import java.util.List;

import org.forester.archaeopteryx.MainFrameApplication;
import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class GoAnnotation extends RunnableProcess {

    private static final String        SYMBOL   = "Symbol";
    private static final String        ASPECT   = "Aspect";
    private static final String        DB       = "DB";
    private static final String        EVIDENCE = "Evidence";
    private static final String        GO_NAME  = "GO Name";
    private static final String        GO_ID    = "GO ID";
    private final Phylogeny            _phy;
    private final MainFrameApplication _mf;
    private final TreePanel            _treepanel;

    public GoAnnotation( final MainFrameApplication mf, final TreePanel treepanel, final Phylogeny phy ) {
        _phy = phy;
        _mf = mf;
        _treepanel = treepanel;
    }

    private void annotate() {
        start( _mf, "GO annotate" );
        for( final PhylogenyNodeIterator iter = _phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( ( node.getNodeData().getSequences() != null ) && !node.getNodeData().getSequences().isEmpty() ) {
                for( final Sequence seq : node.getNodeData().getSequences() ) {
                    if ( ( ( seq.getAccession() != null ) && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) && ( seq
                            .getAnnotations() == null ) ) || seq.getAnnotations().isEmpty() ) {
                        final Accession acc = seq.getAccession();
                        try {
                            final URL url = new URL( "http://www.ebi.ac.uk/QuickGO/GAnnotation?protein="
                                    + acc.getValue() + "&format=tsv" );
                            final HttpURLConnection url_connection = ( HttpURLConnection ) url.openConnection();
                            final BufferedReader br = new BufferedReader( new InputStreamReader( url_connection.getInputStream() ) );
                            final List<String> columns = Arrays.asList( br.readLine().split( "\t" ) );
                            System.out.println( columns );
                            final int db_index = columns.indexOf( DB );
                            final int goid_index = columns.indexOf( GO_ID );
                            final int name_index = columns.indexOf( GO_NAME );
                            final int evidence_index = columns.indexOf( EVIDENCE );
                            final int taxon_index = columns.indexOf( "Taxon" );
                            final int qualifier_index = columns.indexOf( "Qualifier" );
                            final int reference_index = columns.indexOf( "Reference" );
                            final int symbol_index = columns.indexOf( SYMBOL );
                            final int splice_index = columns.indexOf( "Splice" );
                            final int with_index = columns.indexOf( "With" );
                            final int aspect_index = columns.indexOf( ASPECT );
                            final int source_index = columns.indexOf( "Source" );
                            String line;
                            while ( ( line = br.readLine() ) != null ) {
                                final String[] fields = line.split( "\t" );
                                final Annotation a = new Annotation( fields[ goid_index ] );
                                a.setDesc( name_index >= 0 ? fields[ name_index ] : "" );
                                a.setSource( db_index >= 0 ? fields[ db_index ] : "" );
                                a.setEvidence( evidence_index >= 0 ? fields[ evidence_index ] : "" );
                                a.setType( aspect_index >= 0 ? fields[ aspect_index ] : "" );
                                seq.addAnnotation( a );
                                if ( ForesterUtil.isEmpty( seq.getSymbol() ) && ( symbol_index >= 0 )
                                        && !ForesterUtil.isEmpty( fields[ symbol_index ] ) ) {
                                    seq.setSymbol( fields[ symbol_index ] );
                                }
                                System.out.println( DB + ": " + fields[ db_index ] );
                                System.out.println( GO_ID + ": " + fields[ goid_index ] );
                                System.out.println( GO_NAME + ": " + fields[ name_index ] );
                                System.out.println( EVIDENCE + ": " + fields[ evidence_index ] );
                                System.out.println( " taxon" + ": " + fields[ taxon_index ] );
                                System.out.println( " qualifier" + ": " + fields[ qualifier_index ] );
                                System.out.println( " reference" + ": " + fields[ reference_index ] );
                                System.out.println( SYMBOL + ": " + fields[ symbol_index ] );
                                System.out.println( " splice" + ": " + fields[ splice_index ] );
                                System.out.println( " with" + ": " + fields[ with_index ] );
                                System.out.println( ASPECT + ": " + fields[ aspect_index ] );
                                System.out.println( " source" + ": " + fields[ source_index ] );
                            }
                            br.close();
                        }
                        catch ( final IOException e ) {
                            // TODO Auto-generated catch block
                            e.printStackTrace();
                        }
                    }
                }
            }
        }
        end( _mf );
        _treepanel.repaint();
        _treepanel.setEdited( true );
    }

    @Override
    public void run() {
        annotate();
    }
}
