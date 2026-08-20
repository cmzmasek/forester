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

package org.forester.demo;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Regenerates the synthetic Auspice / Nextstrain v2 demo dataset {@code forester/demo/nextstrain-ncov.json} (read back
 * by {@link org.forester.io.parsers.json.AuspiceJsonParser}). Dev-only, run once and commit the JSON (users need no
 * build step); excluded from the shipped jar. A hand-authored 20-tip SARS-CoV-2-like tree with rich per-node metadata:
 * <ul>
 * <li>{@code num_date} (increasing root&rarr;tips) &rarr; the Calendar axis + a "Color by" date gradient; internal-node
 *     {@code num_date.confidence} &rarr; the Node Age spindles.</li>
 * <li>{@code div} (cumulative divergence, increasing root&rarr;tips) &rarr; the Time&harr;Divergence branch-length toggle.</li>
 * <li>{@code region} (with internal-node confidence &rarr; the Ancestral-State Pies), {@code country}, {@code host}
 *     (mostly human, plus a mink and a cat, so "Color by host" is interesting), and {@code clade_membership}.</li>
 * <li>{@code branch_attrs.labels.clade} &rarr; the Nextstrain clade name on each internal branch.</li>
 * </ul>
 * Synthetic (no data-licensing issues); dates/divisions are deliberately monotonic (every child &ge; its parent).
 */
public final class AuspiceDemoGenerator {

    public static void main( final String[] args ) throws IOException {
        final File dir = new File( System.getProperty( "user.dir" ), "forester/demo" );
        final N root = buildTree();
        final StringBuilder sb = new StringBuilder();
        sb.append( "{\n" );
        sb.append( "  \"version\": \"v2\",\n" );
        sb.append( "  \"meta\": {\n" );
        sb.append( "    \"title\": \"Synthetic SARS-CoV-2 phylodynamics (Archaeopteryx demo)\",\n" );
        sb.append( "    \"updated\": \"2026-01-01\"\n" );
        sb.append( "  },\n" );
        sb.append( "  \"tree\":\n" );
        emit( sb, root, "  " );
        sb.append( "\n}\n" );
        final File out = new File( dir, "nextstrain-ncov.json" );
        Files.writeString( out.toPath(), sb.toString() );
        System.out.println( "wrote " + countTips( root ) + "-tip Auspice demo -> " + out.getAbsolutePath() );
    }

    /** A 20-tip SARS-CoV-2-like clade tree. Dates and divergences are monotonic root->tips. */
    private static N buildTree() {
        return internal( "19A", 2019.90, 0.0, 2019.80, 2019.98, "Asia", conf( "Asia", 0.95, "Europe", 0.05 ), "19A",
                internal( "19B", 2019.93, 0.0002, 2019.86, 2020.00, "Asia", conf( "Asia", 0.9, "Europe", 0.1 ), "19B",
                        leaf( "China/Wuhan-Hu-1/2019", 2019.96, 0.0003, "Asia", "China", "Homo sapiens", "19B" ),
                        leaf( "Japan/TY/2020", 2020.05, 0.0004, "Asia", "Japan", "Homo sapiens", "19B" ),
                        leaf( "SouthKorea/KCDC/2020", 2020.10, 0.0005, "Asia", "South Korea", "Homo sapiens", "19B" ) ),
                internal( "20A", 2020.05, 0.0006, 2019.98, 2020.12, "Europe", conf( "Europe", 0.8, "Asia", 0.2 ), "20A",
                        leaf( "Italy/LOM/2020", 2020.18, 0.0008, "Europe", "Italy", "Homo sapiens", "20A" ),
                        leaf( "France/IDF/2020", 2020.22, 0.0009, "Europe", "France", "Homo sapiens", "20A" ),
                        leaf( "Spain/CAT/2020", 2020.25, 0.0010, "Europe", "Spain", "Homo sapiens", "20A" ),
                        leaf( "Belgium/ANT/2020", 2020.28, 0.0010, "Europe", "Belgium", "Homo sapiens", "20A" ),
                        internal( "20B", 2020.30, 0.0011, 2020.22, 2020.38, "Europe",
                                conf( "Europe", 0.85, "North America", 0.15 ), "20B",
                                leaf( "England/SANG/2020", 2020.40, 0.0013, "Europe", "United Kingdom", "Homo sapiens",
                                        "20B" ),
                                leaf( "Germany/BAV/2020", 2020.45, 0.0014, "Europe", "Germany", "Homo sapiens", "20B" ),
                                leaf( "Denmark/mink-cluster5/2020", 2020.55, 0.0016, "Europe", "Denmark",
                                        "Neovison vison", "20B" ) ),
                        internal( "20C", 2020.32, 0.0012, 2020.24, 2020.40, "North America",
                                conf( "North America", 0.9, "Europe", 0.1 ), "20C",
                                leaf( "USA/CA-CDC/2020", 2020.50, 0.0015, "North America", "USA", "Homo sapiens",
                                        "20C" ),
                                leaf( "USA/cat-NY/2020", 2020.52, 0.0015, "North America", "USA", "Felis catus",
                                        "20C" ),
                                leaf( "Canada/ON-PHL/2020", 2020.58, 0.0017, "North America", "Canada", "Homo sapiens",
                                        "20C" ) ) ),
                internal( "20I", 2020.55, 0.0018, 2020.45, 2020.65, "Europe",
                        conf( "Europe", 0.7, "Africa", 0.2, "Asia", 0.1 ), "20I (Alpha)",
                        leaf( "England/alpha-B117/2020", 2020.90, 0.0020, "Europe", "United Kingdom", "Homo sapiens",
                                "20I (Alpha)" ),
                        internal( "21J", 2021.20, 0.0022, 2021.05, 2021.35, "Asia",
                                conf( "Asia", 0.8, "North America", 0.2 ), "21J (Delta)",
                                leaf( "India/MH-delta/2021", 2021.35, 0.0024, "Asia", "India", "Homo sapiens",
                                        "21J (Delta)" ),
                                leaf( "India/KA-delta/2021", 2021.42, 0.0025, "Asia", "India", "Homo sapiens",
                                        "21J (Delta)" ) ),
                        internal( "21K", 2021.85, 0.0026, 2021.70, 2021.95, "Africa",
                                conf( "Africa", 0.6, "Europe", 0.2, "North America", 0.2 ), "21K (Omicron)",
                                leaf( "SouthAfrica/omicron-BA1/2021", 2021.90, 0.0027, "Africa", "South Africa",
                                        "Homo sapiens", "21K (Omicron)" ),
                                leaf( "USA/omicron-CA/2022", 2022.05, 0.0029, "North America", "USA", "Homo sapiens",
                                        "21K (Omicron)" ),
                                leaf( "England/omicron/2022", 2022.10, 0.0030, "Europe", "United Kingdom",
                                        "Homo sapiens", "21K (Omicron)" ),
                                leaf( "Australia/VIC-omicron/2022", 2022.30, 0.0032, "Oceania", "Australia",
                                        "Homo sapiens", "21K (Omicron)" ) ) ) );
    }

    // ---- node model + builders --------------------------------------------------------------------------------

    private static final class N {
        String              name;                        // leaf name (also used as the internal-node name here)
        double              date, div;
        double              date_lo, date_hi;            // internal-node num_date confidence (0 => none)
        String              region, country, host, clade; // per-node traits
        Map<String, Double> region_conf;                 // internal-node region posterior (null on leaves)
        String              clade_label;                 // branch_attrs.labels.clade (internal)
        final List<N>       children = new ArrayList<>();
    }

    private static N leaf( final String name, final double date, final double div, final String region,
                           final String country, final String host, final String clade ) {
        final N n = new N();
        n.name = name;
        n.date = date;
        n.div = div;
        n.region = region;
        n.country = country;
        n.host = host;
        n.clade = clade;
        return n;
    }

    private static N internal( final String name, final double date, final double div, final double lo, final double hi,
                               final String region, final Map<String, Double> region_conf, final String clade_label,
                               final N... kids ) {
        final N n = new N();
        n.name = name;
        n.date = date;
        n.div = div;
        n.date_lo = lo;
        n.date_hi = hi;
        n.region = region;
        n.region_conf = region_conf;
        n.clade_label = clade_label;
        for ( final N k : kids ) {
            n.children.add( k );
        }
        return n;
    }

    private static Map<String, Double> conf( final Object... kv ) {
        final Map<String, Double> m = new LinkedHashMap<>();
        for ( int i = 0; i < kv.length; i += 2 ) {
            m.put( (String) kv[ i ], ( (Number) kv[ i + 1 ] ).doubleValue() );
        }
        return m;
    }

    private static int countTips( final N n ) {
        if ( n.children.isEmpty() ) {
            return 1;
        }
        int t = 0;
        for ( final N c : n.children ) {
            t += countTips( c );
        }
        return t;
    }

    // ---- JSON emitter -----------------------------------------------------------------------------------------

    private static void emit( final StringBuilder sb, final N n, final String ind ) {
        sb.append( ind ).append( "{\n" );
        if ( n.name != null ) {
            sb.append( ind ).append( "  \"name\": " ).append( q( n.name ) ).append( ",\n" );
        }
        sb.append( ind ).append( "  \"node_attrs\": {\n" );
        sb.append( ind ).append( "    \"div\": " ).append( num( n.div ) );
        sb.append( ",\n" ).append( ind ).append( "    \"num_date\": { \"value\": " ).append( num( n.date ) );
        if ( n.date_hi > 0 ) {
            sb.append( ", \"confidence\": [" ).append( num( n.date_lo ) ).append( ", " ).append( num( n.date_hi ) )
                    .append( "]" );
        }
        sb.append( " }" );
        if ( n.region != null ) {
            sb.append( ",\n" ).append( ind ).append( "    \"region\": { \"value\": " ).append( q( n.region ) );
            if ( n.region_conf != null ) {
                sb.append( ", \"confidence\": " ).append( confMap( n.region_conf ) );
            }
            sb.append( " }" );
        }
        if ( n.country != null ) {
            sb.append( ",\n" ).append( ind ).append( "    \"country\": { \"value\": " ).append( q( n.country ) )
                    .append( " }" );
        }
        if ( n.host != null ) {
            sb.append( ",\n" ).append( ind ).append( "    \"host\": { \"value\": " ).append( q( n.host ) )
                    .append( " }" );
        }
        if ( n.clade != null ) {
            sb.append( ",\n" ).append( ind ).append( "    \"clade_membership\": { \"value\": " ).append( q( n.clade ) )
                    .append( " }" );
        }
        sb.append( "\n" ).append( ind ).append( "  }" );
        if ( n.clade_label != null ) {
            sb.append( ",\n" ).append( ind ).append( "  \"branch_attrs\": { \"labels\": { \"clade\": " )
                    .append( q( n.clade_label ) ).append( " } }" );
        }
        if ( !n.children.isEmpty() ) {
            sb.append( ",\n" ).append( ind ).append( "  \"children\": [\n" );
            for ( int i = 0; i < n.children.size(); ++i ) {
                emit( sb, n.children.get( i ), ind + "    " );
                sb.append( i < ( n.children.size() - 1 ) ? ",\n" : "\n" );
            }
            sb.append( ind ).append( "  ]" );
        }
        sb.append( "\n" ).append( ind ).append( "}" );
    }

    private static String confMap( final Map<String, Double> m ) {
        final StringBuilder sb = new StringBuilder( "{ " );
        int i = 0;
        for ( final Map.Entry<String, Double> e : m.entrySet() ) {
            if ( i++ > 0 ) {
                sb.append( ", " );
            }
            sb.append( q( e.getKey() ) ).append( ": " ).append( num( e.getValue() ) );
        }
        return sb.append( " }" ).toString();
    }

    /** A plain-decimal number (no scientific notation), so the JSON source reads cleanly (e.g. 0.0002, not 2.0E-4). */
    private static String num( final double d ) {
        if ( d == 0.0 ) {
            return "0";
        }
        return BigDecimal.valueOf( d ).stripTrailingZeros().toPlainString();
    }

    private static String q( final String s ) {
        return "\"" + s.replace( "\\", "\\\\" ).replace( "\"", "\\\"" ) + "\"";
    }

    private AuspiceDemoGenerator() {
    }
}
