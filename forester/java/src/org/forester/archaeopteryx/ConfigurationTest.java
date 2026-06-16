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

package org.forester.archaeopteryx;

import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;

/**
 * Unit tests for {@link Configuration}. Archaeopteryx no longer reads configuration files, so the
 * (now sole, no-argument) constructor must simply produce the built-in defaults. This is headless
 * and runs as part of the suite via {@link #test()}.
 */
public final class ConfigurationTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "Configuration: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return testDefaults() && testClickToOptionAlignment();
    }

    /**
     * The click-to options use a fragile parallel-array scheme: each {@code final static int} index
     * must address the matching row of {@link Configuration#clickto_options}. Guards that coupling
     * (notably after "Change Node Font(s)" was removed and the trailing indices shifted down).
     */
    private static boolean testClickToOptionAlignment() {
        final Configuration c = new Configuration();
        if ( !"Display Node Data".equals( c.getClickToTitle( Configuration.display_node_data ) ) ) {
            return false;
        }
        if ( !"Colorize Node(s)".equals( c.getClickToTitle( Configuration.color_node_font ) ) ) {
            return false;
        }
        if ( !"Colorize Subtree(s)".equals( c.getClickToTitle( Configuration.color_subtree ) ) ) {
            return false;
        }
        if ( !"List Node Data".equals( c.getClickToTitle( Configuration.get_ext_desc_data ) ) ) {
            return false;
        }
        if ( !"Order Subtree".equals( c.getClickToTitle( Configuration.order_subtree ) ) ) {
            return false;
        }
        // the removed "Change Node Font(s)" option must no longer appear anywhere in the table
        for ( final String[] option : Configuration.clickto_options ) {
            if ( "Change Node Font(s)".equals( option[ 0 ] ) ) {
                return false;
            }
        }
        // the last index must still address the final row (no off-by-one after the shift)
        if ( Configuration.order_subtree != ( Configuration.clickto_options.length - 1 ) ) {
            return false;
        }
        return true;
    }

    /** The no-arg constructor yields the built-in defaults (no file is read). */
    private static boolean testDefaults() {
        final Configuration c = new Configuration();
        // phyloXML XSD validation defaults on and can no longer be switched off (the only switch
        // was the removed config-file parser), so it must be true.
        if ( !c.isValidatePhyloXmlAgainstSchema() ) {
            return false;
        }
        // the display-color map is still initialized (TreeColorSet reads it) but, with no parser to
        // populate it, must be empty.
        if ( ( c.getDisplayColors() == null ) || !c.getDisplayColors().isEmpty() ) {
            return false;
        }
        if ( c.getCladogramType() != AptxConstants.CLADOGRAM_TYPE_DEFAULT ) {
            return false;
        }
        if ( c.getDefaultNodeShapeSize() != AptxConstants.DEFAULT_NODE_SHAPE_SIZE_DEFAULT ) {
            return false;
        }
        if ( c.getTaxonomyExtraction() != TAXONOMY_EXTRACTION.NO ) {
            return false;
        }
        if ( c.isMidpointReroot() ) {
            return false;
        }
        return true;
    }

    private ConfigurationTest() {
    }
}
