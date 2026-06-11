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

package org.forester.test.examples;

import java.io.File;
import java.io.IOException;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.util.ForesterUtil;

public class Example4 {

    public static void main( final String[] args ) {
        // Reading-in of (a) tree(s) from a file.
        final File treefile = new File( "/home/czmasek/tol_117_TEST.xml" );
        PhylogenyParser parser = null;
        try {
            parser = ParserUtils.createParserDependingOnFileType( treefile, true );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        Phylogeny[] phys = null;
        try {
            phys = PhylogenyMethods.readPhylogenies( parser, treefile );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        // Writing trees to a file.
        final File outfile = new File( "/home/czmasek/tol_117_TEST_out.xml" );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( phys, 0, outfile, ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
        }
    }
}
