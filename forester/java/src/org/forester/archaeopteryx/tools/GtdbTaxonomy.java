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

package org.forester.archaeopteryx.tools;

import java.util.LinkedHashMap;
import java.util.Map;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Parses a GTDB (Genome Taxonomy Database) taxonomy string -- the 7-rank, prefix-tagged lineage that GTDB and GTDB-Tk
 * emit, e.g.
 * <pre>d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli</pre>
 * and writes it onto a tree tip so the genome-based bacterial/archaeal taxonomy can drive Archaeopteryx's existing
 * visualizations OFFLINE (no network lookup): each rank becomes a categorical {@code gtdb:<rank>} node PROPERTY (so
 * "Color by" / "Annotation Columns" / search work at any rank), and the most specific rank present becomes a
 * {@code <taxonomy>} (so Display Data: Taxonomy shows the GTDB name).
 *
 * <p>GTDB uses the standardized seven ranks; the prefixes map as d__=domain, p__=phylum, c__=class, o__=order,
 * f__=family, g__=genus, s__=species. Empty ranks (e.g. a bare {@code s__} for an unnamed species) are skipped.
 * Pure and Swing-free.</p>
 *
 * @see <a href="https://gtdb.ecogenomic.org">https://gtdb.ecogenomic.org</a> (Parks et al.)
 */
public final class GtdbTaxonomy {

    /** The GTDB single-letter rank prefixes, most general first, mapped to their full rank names. */
    private static final char[]   PREFIXES = { 'd', 'p', 'c', 'o', 'f', 'g', 's' };
    private static final String[] RANKS    = { "domain", "phylum", "class", "order", "family", "genus", "species" };
    /** The property namespace for the imported GTDB ranks. User-visible (not the internal {@code aptx:} namespace). */
    public static final String    REF_PREFIX = "gtdb:";

    /** The full rank name for a GTDB single-letter prefix (d/p/c/o/f/g/s), or {@code null} if unrecognized. */
    public static String rankForPrefix( final char prefix ) {
        for( int i = 0; i < PREFIXES.length; ++i ) {
            if ( PREFIXES[ i ] == prefix ) {
                return RANKS[ i ];
            }
        }
        return null;
    }

    /** True if {@code s} looks like a GTDB lineage string -- it carries at least one {@code X__Name} rank tag at the
     *  START of a {@code ';'}-delimited field. Delegates to {@link #parse} so the detector can't drift from the parser:
     *  a stray value that merely CONTAINS a letter+{@code "__"} mid-token (e.g. {@code "photo__album"}, {@code
     *  "sample_g__bin"}) is NOT mistaken for a classification. Used to auto-detect the classification column. */
    public static boolean looksLikeGtdb( final String s ) {
        return !parse( s ).isEmpty();
    }

    /**
     * Parses a GTDB classification string into an ordered {@code rank -> name} map (domain first, species last),
     * skipping ranks with no name. Tolerant of surrounding whitespace, a trailing ';', and unrecognized fragments
     * (which are ignored). Never null; empty if {@code classification} carries no named rank.
     */
    public static Map<String, String> parse( final String classification ) {
        final Map<String, String> lineage = new LinkedHashMap<String, String>();
        if ( classification == null ) {
            return lineage;
        }
        for( final String raw : classification.split( ";" ) ) {
            final String part = raw.trim();
            // a field is "X__Name": a single prefix letter, then "__", then the (possibly empty) name
            if ( ( part.length() < 3 ) || ( part.charAt( 1 ) != '_' ) || ( part.charAt( 2 ) != '_' ) ) {
                continue;
            }
            final String rank = rankForPrefix( part.charAt( 0 ) );
            if ( rank == null ) {
                continue;
            }
            final String name = part.substring( 3 ).trim();
            if ( name.length() > 0 ) {
                lineage.put( rank, name );
            }
        }
        return lineage;
    }

    /**
     * Applies a GTDB classification to a tip: writes a {@code gtdb:<rank>} categorical property for each named rank and,
     * for the most specific rank present, a {@code <taxonomy>} (scientific name + rank) so taxonomy display shows the
     * GTDB name. ALL existing {@code gtdb:*} properties are cleared first and rewritten, so a re-import fully REPLACES
     * the GTDB annotation -- re-importing the same table is idempotent, and re-importing a shorter/different lineage
     * leaves no stale deeper-rank properties. Returns the number of ranks written (0 if nothing parsed, in which case
     * the node is left untouched).
     */
    public static int applyToNode( final PhylogenyNode node, final String classification )
            throws PhyloXmlDataFormatException {
        final Map<String, String> lineage = parse( classification );
        if ( lineage.isEmpty() ) {
            return 0;
        }
        PropertiesList props = node.getNodeData().getProperties();
        if ( props == null ) {
            props = new PropertiesList();
            node.getNodeData().setProperties( props );
        }
        removeAllGtdbProperties( props ); // full replace on re-import (no duplicates, no stale deeper ranks)
        String most_specific_rank = null;
        String most_specific_name = null;
        for( final Map.Entry<String, String> e : lineage.entrySet() ) {
            final String ref = REF_PREFIX + e.getKey();
            props.addProperty( new Property( ref, e.getValue(), "", "xsd:string", AppliesTo.NODE ) );
            most_specific_rank = e.getKey(); // lineage is ordered domain..species, so the last is the most specific
            most_specific_name = e.getValue();
        }
        // the most specific rank present becomes the tip's taxonomy (so Display Data: Taxonomy shows the GTDB name)
        final Taxonomy tax = new Taxonomy();
        tax.setScientificName( most_specific_name );
        tax.setRank( most_specific_rank );
        node.getNodeData().setTaxonomy( tax );
        return lineage.size();
    }

    /**
     * Applies GTDB classifications to a whole tree by TIP NAME: every external tip whose name is a key in
     * {@code classificationByTipName} gets {@link #applyToNode(PhylogenyNode, String)}. Returns the number of tips
     * actually annotated (a matched-but-unparseable classification counts for nothing). Pure.
     */
    public static int applyByTipName( final org.forester.phylogeny.Phylogeny phy,
                                      final Map<String, String> classificationByTipName )
            throws PhyloXmlDataFormatException {
        int annotated = 0;
        if ( ( phy == null ) || phy.isEmpty() || ( classificationByTipName == null ) ) {
            return 0;
        }
        for( final PhylogenyNode tip : phy.getExternalNodes() ) {
            final String cls = classificationByTipName.get( tip.getName() );
            if ( ( cls != null ) && ( applyToNode( tip, cls ) > 0 ) ) {
                ++annotated;
            }
        }
        return annotated;
    }

    /** Remove every {@code gtdb:*} property (the backing list is {@code getProperties()}), so a re-import fully
     *  replaces the GTDB annotation rather than leaving stale ranks from a previous, longer classification. */
    private static void removeAllGtdbProperties( final PropertiesList props ) {
        props.getProperties().removeIf( p -> ( p.getRef() != null ) && p.getRef().startsWith( REF_PREFIX ) );
    }

    private GtdbTaxonomy() {
    }
}
