// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2009 Christian M. Zmasek
// Copyright (C) 2009 Burnham Institute for Medical Research
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

package org.forester.application;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;

import org.forester.go.PfamToGoMapping;
import org.forester.go.PfamToGoParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class pfam_go {

    private static final String ALLOW_DUPLICATES_OPTION = "d";
    final static private String HELP_OPTION_1           = "help";
    final static private String HELP_OPTION_2           = "h";
    final static private String PRG_NAME                = "pfam_go";
    final static private String PRG_VERSION             = "1.10";
    final static private String PRG_DATE                = "2011.06.26";
    final static private String E_MAIL                  = "czmasek@burnham.org";
    final static private String WWW                     = "www.phylosoft.org";

    private static void process( final File pfams_file,
                                 final List<PfamToGoMapping> mappings,
                                 final boolean allow_duplicates ) throws IOException {
        final BufferedReader reader = ForesterUtil.obtainReader( pfams_file );
        String line = "";
        int found_count = 0;
        int not_found_count = 0;
        Collection<String> encountered_domains = null;
        if ( allow_duplicates ) {
            encountered_domains = new ArrayList<String>();
        }
        else {
            encountered_domains = new HashSet<String>();
        }
        while ( ( line = reader.readLine() ) != null ) {
            line = line.trim();
            if ( ForesterUtil.isEmpty( line ) || line.startsWith( "##" ) ) {
                continue;
            }
            else if ( line.startsWith( "#" ) ) {
                encountered_domains.clear();
                line = line.replace( '#', '>' );
                System.out.println( line );
            }
            else {
                if ( allow_duplicates || !encountered_domains.contains( line ) ) {
                    encountered_domains.add( line );
                    boolean found = false;
                    for( final PfamToGoMapping mapping : mappings ) {
                        if ( mapping.getKey().getId().equals( line ) ) {
                            System.out.println( mapping.getValue() );
                            found = true;
                        }
                    }
                    if ( found ) {
                        found_count++;
                    }
                    else {
                        not_found_count++;
                    }
                }
                else {
                    System.err.println( "# duplicate domain: " + line );
                }
            }
        }
        System.out.println( "# pfams with mapping to GO   : " + found_count );
        System.out.println( "# pfams without mapping to GO: " + not_found_count );
        reader.close();
    }

    public static void main( final String args[] ) {
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length == 0 ) ) {
            printHelp();
            System.exit( 0 );
        }
        final List<String> allowed_options = new ArrayList<String>();
        allowed_options.add( ALLOW_DUPLICATES_OPTION );
        if ( ( cla.getNumberOfNames() != 2 ) && ( cla.getNumberOfNames() != 3 ) ) {
            printHelp();
            System.exit( -1 );
        }
        final String dissallowed_options = cla.validateAllowedOptionsAsString( allowed_options );
        if ( dissallowed_options.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + dissallowed_options );
        }
        final File pfam2go_file = cla.getFile( 0 );
        final File pfams_file = cla.getFile( 1 );
        boolean allow_duplicates = false;
        if ( cla.isOptionSet( ALLOW_DUPLICATES_OPTION ) ) {
            allow_duplicates = true;
        }
        final PfamToGoParser pfam2go_parser = new PfamToGoParser( pfam2go_file );
        List<PfamToGoMapping> mappings = null;
        try {
            mappings = pfam2go_parser.parse();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        try {
            process( pfams_file, mappings, allow_duplicates );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        System.out.println();
    }

    private static void printHelp() {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE, E_MAIL, WWW );
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " [-" + ALLOW_DUPLICATES_OPTION
                + " to allow duplicates] <pfam2go file> <file with pfams>" );
        System.out.println();
        System.out.println();
    }
}