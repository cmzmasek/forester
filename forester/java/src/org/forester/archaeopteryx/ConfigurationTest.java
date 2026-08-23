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
        return testDefaults() && testOptionEnums();
    }

    /**
     * The click-to actions and "Display Data" options are now the type-safe enums {@link ClickToOption} /
     * {@link DisplayOption} (they replaced the old fragile positional {@code String[][]} arrays + parallel
     * int index constants, so the previous alignment guard is obsolete). Spot-check that each enum carries a
     * title and that a couple of representative default states are correct, and that the removed actions/
     * options ("Launch BLAST", "Open PDB", "Binary Characters", "Vector Data", ...) are gone.
     */
    private static boolean testOptionEnums() {
        if ( !"Display Node Data".equals( ClickToOption.DISPLAY_NODE_DATA.title() )
                || !"Ladderize Subtree".equals( ClickToOption.ORDER_SUBTREE.title() )
                || !"Open Taxonomy DB".equals( ClickToOption.OPEN_TAX_WEB.title() ) ) {
            return false;
        }
        for ( final ClickToOption o : ClickToOption.values() ) {
            if ( ( o.title() == null ) || o.title().isEmpty() ) {
                return false;
            }
            // the removed BLAST / PDB actions must be gone
            if ( "Launch BLAST".equals( o.title() ) || "Open PDB".equals( o.title() ) ) {
                return false;
            }
        }
        // DefaultState replaces the old "yes"/"no"/"?" third column
        if ( !DisplayOption.SHOW_NODE_NAMES.isCheckedByDefault()
                || DisplayOption.SHOW_DOMAIN_ARCHITECTURES.isCheckedByDefault()
                || !DisplayOption.DISPLAY_AS_PHYLOGRAM.isGuessedByDefault() ) {
            return false;
        }
        for ( final DisplayOption o : DisplayOption.values() ) {
            if ( ( o.title() == null ) || o.title().isEmpty() ) {
                return false;
            }
            // the removed display features must be gone
            final String t = o.title();
            if ( "Binary Characters".equals( t ) || "Binary Char Counts".equals( t ) || "Vector Data".equals( t )
                    || "Relation Confidence".equals( t ) || "Multiple Seq Alignment".equals( t )
                    || "Taxonomy Images".equals( t ) ) {
                return false;
            }
        }
        return true;
    }

    /**
     * The no-arg constructor yields the built-in defaults (no file is read). Configuration no longer seeds
     * {@link Options} (those display defaults now live solely in {@code Options.init()} / the enums), so this
     * only exercises the surface Configuration still owns: XSD validation, the tip-color map, the parse-time
     * taxonomy-extraction knob, the default click-to action, and the editable / draw-as-phylogram convenience.
     */
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
        if ( c.getTaxonomyExtraction() != TAXONOMY_EXTRACTION.NO ) {
            return false;
        }
        // the click-to dropdown defaults to "Display Node Data"
        if ( c.getDefaultDisplayClicktoOption() != ClickToOption.DISPLAY_NODE_DATA ) {
            return false;
        }
        // trees are editable by default; the default display type is a cladogram (phylogram is only GUESSed)
        if ( !c.isEditable() || c.isDrawAsPhylogram() ) {
            return false;
        }
        return true;
    }

    private ConfigurationTest() {
    }
}
