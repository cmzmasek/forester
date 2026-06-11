// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.go.GoId;
import org.forester.go.GoTerm;
import org.forester.go.GoUtils;
import org.forester.go.OBOparser;
import org.forester.go.PfamToGoMapping;
import org.forester.go.PfamToGoParser;

public class pfam2go_extractor {

    final static private String PRG_NAME = "pfam2go_extractor";

    public static void main( final String args[] ) {
        if ( args.length < 3 ) {
            printHelp();
        }
        final PfamToGoParser p = new PfamToGoParser( new File( args[ 0 ] ) );
        List<PfamToGoMapping> pfam2go = null;
        try {
            pfam2go = p.parse();
        }
        catch ( final IOException e ) {
            printHelp();
            e.printStackTrace();
        }
        final OBOparser parser = new OBOparser( new File( args[ 1 ] ), OBOparser.ReturnType.BASIC_GO_TERM );
        List<GoTerm> all_go_terms = null;
        try {
            all_go_terms = parser.parse();
        }
        catch ( final IOException e ) {
            printHelp();
            e.printStackTrace();
        }
        final Map<GoId, GoTerm> goid_to_term_map = GoUtils.createGoIdToGoTermMap( all_go_terms );
        System.out.println( "# pfam2go : " + args[ 0 ] );
        System.out.println( "# OBO file: " + args[ 1 ] );
        final GoId[] queries = new GoId[ args.length - 2 ];
        for( int i = 2; i < args.length; ++i ) {
            queries[ i - 2 ] = new GoId( args[ i ] );
            System.out.println( "# " + ( i - 2 ) + ": " + queries[ i - 2 ].getId() + " = "
                    + goid_to_term_map.get( queries[ i - 2 ] ).getName() + " ("
                    + goid_to_term_map.get( queries[ i - 2 ] ).getDefinition() + ")" );
        }
        final SortedSet<String> pfams = new TreeSet<String>();
        for( final PfamToGoMapping pfam_to_go_mapping : pfam2go ) {
            final String domain_id = pfam_to_go_mapping.getKey();
            final GoId go_id = pfam_to_go_mapping.getValue();
            final Set<GoId> supers = GoUtils.getAllSuperGoIds( go_id, goid_to_term_map );
            supers.add( go_id );
            for( final GoId querie : queries ) {
                if ( supers.contains( querie ) ) {
                    pfams.add( domain_id.toString() );
                }
            }
        }
        for( final String pfam : pfams ) {
            System.out.println( pfam );
        }
    }

    private static void printHelp() {
        System.out.println();
        System.out.println( PRG_NAME
                            + " <pfam2go mapping file> <file with all GO terms, in 'obo' format> <GO id> [more GO ids]" );
        System.out.println();
    }
}
