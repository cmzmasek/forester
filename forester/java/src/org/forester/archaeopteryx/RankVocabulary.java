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

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * The taxonomic ranks "Annotate Clades by Rank" and "Colorize Subtrees via Taxonomic Rank" work with -- the joint
 * vocabulary of Archaeopteryx and Archaeopteryx.js (ported from forester.js, "Taxonomic ranks: clades by rank";
 * Christian, 2026-10-19: "i go with your recommendation, with one exception"). {@link #RANKS} is phyloXML's
 * vocabulary less "other" and "unknown" down to species, in {@code TaxonomyUtil.TAXONOMY_RANKS_LIST}'s order
 * (section and subsection under the genus), then the terms viruses and microbes are grouped by below species. These
 * are not a strict hierarchy across viruses; their order decides only a chooser's listing and the automatic default.
 * <p>
 * A rank is handled by its KEY ({@link #key}): lower case with spaces, hyphens, underscores and dots removed, so
 * "Sub-Species", "sub_species" and "ssp." are one rank; a few abbreviations fold onto their term. Refused, by his
 * word: strain, substrain and isolate (a strain is one tip), and phyloXML's "other" and "unknown". A rank outside
 * the list is still a rank when a {@code <taxonomy>} names it; a property is a rank only by the list or by naming
 * a clade or lineage scheme of the file's own ({@link #propertyRank}). {@code TaxonomyUtil.TAXONOMY_RANKS_LIST}
 * itself is untouched: it is the phyloXML editing vocabulary, and cladinator reads it. Ported from forester.js at
 * cc86fb6 (after Christian's review of the JS side: varietas folds onto variety, 67 ranks; a scheme is a key
 * ENDING in clade or lineage, a level suffix allowed, plus clade_membership by name).
 */
final class RankVocabulary {

    /** The ranks in Linnaean order, then the terms below species; 67 entries. */
    static final List<String> RANKS = Collections.unmodifiableList( Arrays.asList( "domain", "superkingdom", "kingdom",
            "subkingdom", "branch", "infrakingdom", "superphylum", "phylum", "subphylum", "infraphylum", "microphylum",
            "superdivision", "division", "subdivision", "infradivision", "superclass", "class", "subclass", "infraclass",
            "superlegion", "legion", "sublegion", "infralegion", "supercohort", "cohort", "subcohort", "infracohort",
            "superorder", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "supertribe", "tribe",
            "subtribe", "infratribe", "genus", "subgenus", "section", "subsection", "superspecies", "species",
            "subspecies", "variety", "subvariety", "form", "subform", "forma specialis", "cultivar",
            "biovar", "serovar", "pathovar", "morphovar", "serogroup", "serotype", "genogroup", "genotype",
            "subgenotype", "type", "subtype", "clade", "subclade", "lineage", "sublineage", "variant" ) );

    private static final Map<String, String>  ALIASES = new HashMap<String, String>();
    private static final Set<String>          REFUSED = new HashSet<String>( Arrays.asList( "", "other", "unknown",
            "strain", "substrain", "isolate" ) );
    private static final Map<String, Integer> INDEX   = new HashMap<String, Integer>(); // key -> place in RANKS
    private static final Map<String, String>  LABEL   = new HashMap<String, String>();  // key -> canonical spelling
    static {
        ALIASES.put( "ssp", "subspecies" );
        ALIASES.put( "subsp", "subspecies" );
        ALIASES.put( "var", "variety" );
        ALIASES.put( "varietas", "variety" );
        ALIASES.put( "subvar", "subvariety" );
        ALIASES.put( "forma", "form" );
        ALIASES.put( "fsp", "formaspecialis" );
        ALIASES.put( "formasp", "formaspecialis" );
        ALIASES.put( "cv", "cultivar" );
        for( int i = 0; i < RANKS.size(); ++i ) {
            INDEX.put( fold( RANKS.get( i ) ), i );
            LABEL.put( fold( RANKS.get( i ) ), RANKS.get( i ) );
        }
    }

    /** A rank's name folded: trimmed, lower-cased, spaces, hyphens, underscores and dots removed; "" for null. */
    static String fold( final String rank ) {
        return ( rank == null ) ? "" : rank.trim().toLowerCase( Locale.ROOT ).replaceAll( "[\\s_.\\-]+", "" );
    }

    /** A rank's key, or null when it is none, or refused (strain, substrain, isolate, other, unknown, blank). */
    static String key( final String rank ) {
        String k = fold( rank );
        final String alias = ALIASES.get( k );
        if ( alias != null ) {
            k = alias;
        }
        return REFUSED.contains( k ) ? null : k;
    }

    /** Whether {@code key} (from {@link #key}) is one of {@link #RANKS}. */
    static boolean isListed( final String key ) {
        return ( key != null ) && INDEX.containsKey( key );
    }

    /** The canonical spelling of a listed rank's key ("forma specialis" for formaspecialis), null for any other. */
    static String label( final String key ) {
        return ( key == null ) ? null : LABEL.get( key );
    }

    /** The local name of a property reference: what follows the last colon ("vipr:Genus" -> "Genus"). */
    static String localName( final String ref ) {
        return ( ref == null ) ? "" : ref.substring( ref.lastIndexOf( ':' ) + 1 ).trim();
    }

    /**
     * The rank a property reference names, as a key: its local name when that is a listed rank, in ANY namespace
     * ("vipr:Genus" -> "genus", "BVBRC:species", "repseq:family", "ncbi:genus", "tol:Domain"); or a clade or lineage
     * naming scheme of the file's own -- a key ENDING in "clade" or "lineage", a level suffix allowed
     * ("ird:H5Clade" -> h5clade, "vipr:PANGO_Lineage_L1" -> pangolineagel1), plus "clade_membership" by name -- as a
     * rank of its own, under the file's spelling. Null for anything else: a property named host, country or genome
     * type is a field for the Color-by, not a rank, and so are clade_size and lineage_count; and for a refused name
     * (a strain field is one tip).
     */
    static String propertyRank( final String ref ) {
        if ( ref == null ) {
            return null;
        }
        final String key = key( localName( ref ) );
        if ( key == null ) {
            return null;
        }
        if ( INDEX.containsKey( key ) ) {
            return key;
        }
        return ( SCHEME.matcher( key ).matches() || key.equals( "clademembership" ) ) ? key : null;
    }

    /** A clade or lineage naming scheme of the file's own: a key ending in clade or lineage, with at most a level
     *  suffix of a letter and digits after it (h5clade, pangolineagel1, h5clade2). */
    private static final java.util.regex.Pattern SCHEME = java.util.regex.Pattern.compile( ".*(clade|lineage)[a-z]?[0-9]{0,2}" );

    /** The listed ranks first, in their order, then the rest alphabetically (keys). */
    static int compare( final String a, final String b ) {
        final Integer ia = INDEX.get( a );
        final Integer ib = INDEX.get( b );
        if ( ( ia != null ) || ( ib != null ) ) {
            if ( ia == null ) {
                return 1;
            }
            if ( ib == null ) {
                return -1;
            }
            return ia - ib;
        }
        return a.compareTo( b );
    }

    private RankVocabulary() {
    }
}
