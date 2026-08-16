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

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.util.ForesterUtil;

/**
 * Headless support for colouring a tanglegram's connectors by a categorical tip attribute: which fields are available
 * to colour by (taxonomy scientific name / code, and any categorical node property the tips carry), and the
 * value&rarr;colour palette for a chosen field. The palette reuses {@link AptxUtil#assignDistinctColors} (the same
 * deterministic, reproducible spread the rank colorizer uses).
 */
final class TanglegramColoring {

    /** A categorical field connectors can be coloured by (the value comes from a connector's tip). */
    static final class Field {

        enum Kind {
            SCIENTIFIC_NAME,
            TAXONOMY_CODE,
            PROPERTY
        }

        private final String _label;
        private final Kind   _kind;
        private final String _ref; // the property ref, for Kind.PROPERTY

        Field( final String label, final Kind kind, final String ref ) {
            _label = label;
            _kind = kind;
            _ref = ref;
        }

        String label() {
            return _label;
        }

        /** The tip's categorical value for this field, or "" when it carries none. */
        String valueFor( final PhylogenyNode node ) {
            switch ( _kind ) {
                case SCIENTIFIC_NAME:
                    return node.getNodeData().isHasTaxonomy()
                            ? trimOrEmpty( node.getNodeData().getTaxonomy().getScientificName() ) : "";
                case TAXONOMY_CODE:
                    return node.getNodeData().isHasTaxonomy()
                            ? trimOrEmpty( node.getNodeData().getTaxonomy().getTaxonomyCode() ) : "";
                case PROPERTY:
                    return trimOrEmpty( PropertyColorScheme.valueFor( node, _ref ) );
                default:
                    return "";
            }
        }

        private static String trimOrEmpty( final String s ) {
            return ForesterUtil.isEmpty( s ) ? "" : s.trim();
        }
    }

    /** The categorical fields the tree's tips carry, in a stable order (taxonomy first, then property refs). */
    static List<Field> availableFields( final Phylogeny tree ) {
        boolean any_sci = false;
        boolean any_code = false;
        for( final PhylogenyNode tip : TanglegramLinker.externalTipsInDisplayOrder( tree ) ) {
            if ( tip.getNodeData().isHasTaxonomy() ) {
                if ( !ForesterUtil.isEmpty( tip.getNodeData().getTaxonomy().getScientificName() ) ) {
                    any_sci = true;
                }
                if ( !ForesterUtil.isEmpty( tip.getNodeData().getTaxonomy().getTaxonomyCode() ) ) {
                    any_code = true;
                }
            }
        }
        final List<Field> fields = new ArrayList<>();
        if ( any_sci ) {
            fields.add( new Field( "Taxonomy: Scientific Name", Field.Kind.SCIENTIFIC_NAME, null ) );
        }
        if ( any_code ) {
            fields.add( new Field( "Taxonomy: Code", Field.Kind.TAXONOMY_CODE, null ) );
        }
        // categorical property refs = the app's colorable refs minus the numeric (gradient) ones -- reusing
        // PropertyColorScheme's filters so constant and per-tip-unique categorical columns are skipped, exactly as
        // the main "Color by" list does (colouring by a one-value-per-tip column is just noise)
        final Set<String> numeric = new HashSet<>( PropertyColorScheme.numericRefs( tree ) );
        for( final String ref : PropertyColorScheme.colorableRefs( tree ) ) {
            if ( !numeric.contains( ref ) ) {
                fields.add( new Field( ref, Field.Kind.PROPERTY, ref ) );
            }
        }
        return fields;
    }

    /** A distinct colour per non-empty value of {@code field} over the connectors' left tips (deterministic order). */
    static Map<String, Color> colorMap( final List<TanglegramLinker.Link> links, final Field field ) {
        final TreeSet<String> values = new TreeSet<>();
        for( final TanglegramLinker.Link link : links ) {
            final String value = field.valueFor( link.getA() );
            if ( !value.isEmpty() ) {
                values.add( value );
            }
        }
        return AptxUtil.assignDistinctColors( values );
    }

    private TanglegramColoring() {
    }
}
