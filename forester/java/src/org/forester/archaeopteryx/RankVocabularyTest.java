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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The joint rank vocabulary ({@link RankVocabulary}), pinned to the same literals Archaeopteryx.js pins in its
 * test/rank_test.js: the 68 ranks and their order, the key folding and aliases, the refused names (strain,
 * substrain, isolate by Christian's exception; other, unknown), a property's rank in any namespace or by a clade /
 * lineage naming scheme of the file's own, and the ordering of keys.
 */
public final class RankVocabularyTest {

    public static void main( final String[] args ) {
        System.out.println( test() ? "RankVocabularyTest: OK." : "RankVocabularyTest: FAILED." );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        final List<String> R = RankVocabulary.RANKS;
        check( ok, R.size() == 67, "67 ranks, got " + R.size() );
        check( ok, !R.contains( "varietas" ) && "variety".equals( RankVocabulary.key( "varietas" ) ), "varietas folds onto variety" );
        check( ok, !R.contains( "other" ) && !R.contains( "unknown" ) && !R.contains( "strain" ) && !R.contains( "substrain" )
                && !R.contains( "isolate" ), "other, unknown, strain, substrain and isolate are not listed" );
        check( ok, R.indexOf( "species" ) < R.indexOf( "subspecies" ), "species before subspecies" );
        check( ok, R.indexOf( "subtype" ) < R.indexOf( "clade" ), "subtype before clade" );
        check( ok, R.indexOf( "subgenus" ) < R.indexOf( "section" ) && R.indexOf( "subsection" ) < R.indexOf( "species" ),
               "section and subsection sit under the genus" );
        check( ok, R.indexOf( "variant" ) == ( R.size() - 1 ), "variant last" );
        check( ok, R.indexOf( "domain" ) == 0 && R.indexOf( "genus" ) < R.indexOf( "species" ), "domain first, genus before species" );
        // keys
        check( ok, "subspecies".equals( RankVocabulary.key( "Sub-Species" ) ) && "subspecies".equals( RankVocabulary.key( "sub_species" ) )
                && "subspecies".equals( RankVocabulary.key( " ssp. " ) ) && "subspecies".equals( RankVocabulary.key( "subsp" ) ),
               "Sub-Species, sub_species, ssp. and subsp are one rank" );
        check( ok, "formaspecialis".equals( RankVocabulary.key( "Forma specialis" ) ) && "formaspecialis".equals( RankVocabulary.key( "f. sp." ) ),
               "forma specialis and f. sp. fold together" );
        check( ok, "variety".equals( RankVocabulary.key( "var." ) ) && "cultivar".equals( RankVocabulary.key( "cv" ) )
                && "subtype".equals( RankVocabulary.key( "Sub_Type" ) ) && "h5clade".equals( RankVocabulary.key( "H5Clade" ) ),
               "var., cv, Sub_Type, H5Clade" );
        check( ok, "forma specialis".equals( RankVocabulary.label( "formaspecialis" ) ) && "genus".equals( RankVocabulary.label( "genus" ) )
                && ( RankVocabulary.label( "h5clade" ) == null ), "labels: the canonical spelling of a listed key, null otherwise" );
        // refused
        for( final String r : new String[] { "strain", "Sub-Strain", "substrain", "isolate", "other", "unknown", "", null } ) {
            check( ok, RankVocabulary.key( r ) == null, "refused: " + r );
        }
        check( ok, !RankVocabulary.isListed( null ) && RankVocabulary.isListed( "genus" ) && !RankVocabulary.isListed( "h5clade" ), "isListed" );
        // a property's rank
        check( ok, "genus".equals( RankVocabulary.propertyRank( "vipr:Genus" ) ) && "species".equals( RankVocabulary.propertyRank( "BVBRC:species" ) )
                && "family".equals( RankVocabulary.propertyRank( "repseq:family" ) ) && "genus".equals( RankVocabulary.propertyRank( "ncbi:genus" ) )
                && "domain".equals( RankVocabulary.propertyRank( "tol:Domain" ) ) && "genus".equals( RankVocabulary.propertyRank( "Genus" ) ),
               "a listed rank in any namespace, or none" );
        check( ok, "h5clade".equals( RankVocabulary.propertyRank( "ird:H5Clade" ) ) && "pangolineagel1".equals( RankVocabulary.propertyRank( "vipr:PANGO_Lineage_L1" ) )
                && "clademembership".equals( RankVocabulary.propertyRank( "ns:clade_membership" ) ),
               "a clade or lineage naming scheme of the file's own is a rank" );
        check( ok, ( RankVocabulary.propertyRank( "x:clade_size" ) == null ) && ( RankVocabulary.propertyRank( "x:lineage_count" ) == null )
                && "h5clade2".equals( RankVocabulary.propertyRank( "ird:H5_Clade_2" ) ), "clade_size and lineage_count are fields; a level suffix is allowed" );
        check( ok, ( RankVocabulary.propertyRank( "BVBRC:host_group" ) == null ) && ( RankVocabulary.propertyRank( "x:genome_type" ) == null )
                && ( RankVocabulary.propertyRank( "vipr:Strain" ) == null ) && ( RankVocabulary.propertyRank( "y:Isolate" ) == null )
                && ( RankVocabulary.propertyRank( "x:sub-strain" ) == null ) && ( RankVocabulary.propertyRank( null ) == null )
                && ( RankVocabulary.propertyRank( "style:font_color" ) == null ),
               "host, genome type, strain, isolate, substrain, a style: a field is not a rank" );
        check( ok, "Genus".equals( RankVocabulary.localName( "vipr:Genus" ) ) && "Genus".equals( RankVocabulary.localName( "Genus" ) )
                && "".equals( RankVocabulary.localName( null ) ), "local names" );
        // ordering: listed first in their order, then the rest alphabetically
        final List<String> keys = new ArrayList<String>( Arrays.asList( "subtype", "zzz", "species", "aaa", "domain", "genus", "h5clade", "subclade" ) );
        keys.sort( RankVocabulary::compare );
        check( ok, String.join( " ", keys ).equals( "domain genus species subtype subclade aaa h5clade zzz" ), "compare: " + keys );
        return ok[ 0 ];
    }

    private static void check( final boolean[] ok, final boolean c, final String what ) {
        if ( !c ) {
            System.out.println( "  [RankVocabularyTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private RankVocabularyTest() {
    }
}
