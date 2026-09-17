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

package org.forester.io.parsers.util;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.json.AuspiceJsonParser;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.nhx.NHXParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;

/**
 * Headless tests for how {@link ParserUtils} picks a parser: by file name where the name is unambiguous, by the first
 * line where it is not -- and, deliberately, NOT by the two suffixes that tools disagree about.
 */
public final class ParserUtilsTest {

    public static boolean test() {
        return testBySuffix() && testByFirstLine() && testAmbiguousSuffixesAreSniffed();
    }

    private static boolean testBySuffix() {
        final Object[][] cases = { { "aln.NXS", NexusPhylogeniesParser.class }, { "t.nex", NexusPhylogeniesParser.class }, { "t.nexus", NexusPhylogeniesParser.class },
                { "t.nx", NexusPhylogeniesParser.class }, { "t.nwk", NHXParser.class }, { "t.nhx", NHXParser.class },
                { "t.newick", NHXParser.class }, { "t.xml", PhyloXmlParser.class }, { "auspice_tree.json", AuspiceJsonParser.class } };
        for( final Object[] c : cases ) {
            final PhylogenyParser p = ParserUtils.createParserDependingOnSuffixForTest( ( String ) c[ 0 ] );
            if ( ( p == null ) || ( p.getClass() != c[ 1 ] ) ) {
                return fail( c[ 0 ] + " must be read by " + ( ( Class<?> ) c[ 1 ] ).getSimpleName() + ", got "
                        + ( p == null ? "no parser" : p.getClass().getSimpleName() ) );
            }
        }
        // the name decides nothing here: a producer's habit is not a format (".trees" is BEAST's Nexus AND a common name
        // for a Newick gene-tree or bootstrap list; ".tre" and ".t" likewise)
        for( final String name : new String[] { "tree.tre", "run1.t", "tree.txt", "tree", "sequences.contree", "run.trees",
                "genes.trees", "run.con", "run.con.tre" } ) {
            if ( ParserUtils.createParserDependingOnSuffixForTest( name ) != null ) {
                return fail( name + " must be left to the first-line sniff" );
            }
        }
        return true;
    }

    /** ONE ladder for files and URLs. */
    private static boolean testByFirstLine() {
        final Object[][] cases = { { "#NEXUS", NexusPhylogeniesParser.class }, { "  #nexus  ", NexusPhylogeniesParser.class },
                { "# NEXUS", NexusPhylogeniesParser.class }, { "nexus", NexusPhylogeniesParser.class },
                { "Begin trees;", NexusPhylogeniesParser.class }, { "<?xml version=\"1.0\"?>", PhyloXmlParser.class },
                { "  <phyloxml>", PhyloXmlParser.class }, { "((A,B),C);", NHXParser.class }, { "", NHXParser.class },
                { null, NHXParser.class } };
        for( final Object[] c : cases ) {
            final PhylogenyParser p = ParserUtils.createParserFromFirstLine( ( String ) c[ 0 ], false );
            if ( ( p == null ) || ( p.getClass() != c[ 1 ] ) ) {
                return fail( "a first line of '" + c[ 0 ] + "' must be read by " + ( ( Class<?> ) c[ 1 ] ).getSimpleName()
                        + ", got " + ( p == null ? "no parser" : p.getClass().getSimpleName() ) );
            }
        }
        return true;
    }

    /** Why ".tre", ".t", ".trees" and ".con" are not in the suffix list: the SAME name, either format, both read right. */
    private static boolean testAmbiguousSuffixesAreSniffed() {
        try {
            for( final String suffix : new String[] { ".tre", ".t", ".trees", ".con", ".con.tre" } ) {
                final File newick = File.createTempFile( "parserutils", suffix );
                final File nexus = File.createTempFile( "parserutils", suffix );
                newick.deleteOnExit();
                nexus.deleteOnExit();
                Files.write( newick.toPath(), "((A:1,B:1):1,C:2);".getBytes( StandardCharsets.UTF_8 ) );
                Files.write( nexus.toPath(), "#NEXUS\nbegin trees;\n\ttree t1 = ((A:1,B:1):1,C:2);\nend;\n".getBytes( StandardCharsets.UTF_8 ) );
                if ( !( ParserUtils.createParserDependingOnFileType( newick, false ) instanceof NHXParser ) ) {
                    return fail( "a Newick file named *" + suffix + " must be read as Newick" );
                }
                if ( !( ParserUtils.createParserDependingOnFileType( nexus, false ) instanceof NexusPhylogeniesParser ) ) {
                    return fail( "a Nexus file named *" + suffix + " must be read as Nexus" );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "sniffing threw: " + e );
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "ParserUtils test failed: " + msg );
        return false;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }
}
