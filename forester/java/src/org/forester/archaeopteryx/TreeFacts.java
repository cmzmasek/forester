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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.time.Instant;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

/**
 * Everything the Tree Properties window shows read-only, computed in one pass over the tree and returned as plain
 * data: a list of {@link Group}s (File, Structure, Branch lengths, Support values, Annotation coverage, Time
 * axis), each a list of key/value {@link Fact}s plus an optional {@link Histogram}. No Swing in here, so the
 * numbers are testable and the window is only rendering.
 */
final class TreeFacts {

    static final String FILE          = "File";
    static final String STRUCTURE     = "Structure";
    static final String BRANCH_LENGTHS = "Branch lengths";
    static final String SUPPORT       = "Support values";
    static final String COVERAGE      = "Annotation coverage";
    static final String TIME_AXIS     = "Time axis";
    /** Bins of the branch-length / support histograms. */
    static final int    HISTOGRAM_BINS = 12;

    private static final NumberFormat  INT  = NumberFormat.getIntegerInstance( Locale.US );
    private static final DecimalFormat DEC  = new DecimalFormat( "0.######", DecimalFormatSymbols.getInstance( Locale.US ) );
    private static final DateTimeFormatter DATE = DateTimeFormatter.ofPattern( "yyyy-MM-dd HH:mm" );

    private TreeFacts() {
    }

    /** One key/value line. */
    static final class Fact {

        final String key;
        final String value;

        Fact( final String key, final String value ) {
            this.key = key;
            this.value = value;
        }

        @Override
        public String toString() {
            return key + ": " + value;
        }
    }

    /** Equal-width bins over [min, max]; {@code counts[i]} covers {@code [edge(i), edge(i+1))} (the last bin is closed). */
    static final class Histogram {

        final double min;
        final double max;
        final int[]  counts;

        Histogram( final double min, final double max, final int[] counts ) {
            this.min = min;
            this.max = max;
            this.counts = counts;
        }

        double edge( final int i ) {
            return min + ( ( max - min ) * i ) / counts.length;
        }

        int maxCount() {
            int m = 0;
            for( final int c : counts ) {
                m = Math.max( m, c );
            }
            return m;
        }

        int total() {
            int t = 0;
            for( final int c : counts ) {
                t += c;
            }
            return t;
        }
    }

    /** A titled block of facts with a muted detail (shown next to the title) and an optional histogram. */
    static final class Group {

        final String     title;
        final String     detail;
        final List<Fact> facts;
        final Histogram  histogram;

        Group( final String title, final String detail, final List<Fact> facts, final Histogram histogram ) {
            this.title = title;
            this.detail = ( detail == null ) ? "" : detail;
            this.facts = Collections.unmodifiableList( facts );
            this.histogram = histogram;
        }

        /** The value of the fact with this key, or null. */
        String value( final String key ) {
            for( final Fact f : facts ) {
                if ( f.key.equals( key ) ) {
                    return f.value;
                }
            }
            return null;
        }
    }

    /** What the tree panel knows about the time axis of the displayed tree (the tree alone does not). */
    static final class TimeAxis {

        final Options.TIME_AXIS_TYPE type;
        final boolean                dated;
        final String                 unit;
        final double                 root_age_ma;
        final double                 present_date;

        TimeAxis( final Options.TIME_AXIS_TYPE type, final boolean dated, final String unit,
                  final double root_age_ma, final double present_date ) {
            this.type = type;
            this.dated = dated;
            this.unit = unit;
            this.root_age_ma = root_age_ma;
            this.present_date = present_date;
        }
    }

    // ------------------------------------------------------------------ entry point
    /**
     * All groups for {@code phy}. {@code file} may be null (a tree not from / not yet saved to a file);
     * {@code time} may be null (no panel, or the Time-axis group is not wanted). Groups that have nothing to say
     * (no branch lengths, no support values, no time axis) are left out.
     */
    static List<Group> compute( final Phylogeny phy, final File file, final boolean edited, final TimeAxis time ) {
        final List<Group> out = new ArrayList<>();
        out.add( fileGroup( file, edited ) );
        if ( ( phy == null ) || phy.isEmpty() ) {
            return out;
        }
        final Scan scan = new Scan( phy );
        out.add( structureGroup( phy, scan ) );
        final Group bl = branchLengthGroup( phy, scan );
        if ( bl != null ) {
            out.add( bl );
        }
        out.addAll( supportGroups( scan ) );
        out.add( coverageGroup( phy, scan ) );
        final Group ta = timeAxisGroup( time, scan );
        if ( ta != null ) {
            out.add( ta );
        }
        return out;
    }

    // ------------------------------------------------------------------ the one pass
    /** Everything counted in a single preorder walk. */
    static final class Scan {

        int                                 nodes;
        int                                 tips;
        int                                 collapsed;
        int                                 named_internal;
        int                                 tips_taxonomy;
        int                                 tips_tax_id;
        int                                 tips_sequence;
        int                                 tips_mol_seq;
        int                                 tips_domains;
        int                                 tips_date;
        int                                 nodes_date;
        int                                 tips_distribution;
        int                                 tips_reference;
        int                                 nodes_events;
        int                                 duplications;
        int                                 speciations;
        int                                 gene_losses;
        int                                 branches_with_length;
        int                                 zero_length;
        int                                 negative_length;
        final DescriptiveStatistics         lengths      = new BasicDescriptiveStatistics();
        final Map<String, DescriptiveStatistics> support = new LinkedHashMap<>();
        final Map<String, Integer>          property_refs = new TreeMap<>();
        final Set<Taxonomy>                 taxonomies;

        Scan( final Phylogeny phy ) {
            for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                ++nodes;
                final NodeData nd = n.getNodeData();
                final boolean tip = n.isExternal();
                if ( tip ) {
                    ++tips;
                }
                else {
                    if ( n.isCollapse() ) {
                        ++collapsed;
                    }
                    if ( !ForesterUtil.isEmpty( n.getName() ) ) {
                        ++named_internal;
                    }
                }
                if ( !n.isRoot() ) {
                    final double d = n.getDistanceToParent();
                    if ( d == PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
                        // no branch length
                    }
                    else if ( d < 0 ) {
                        ++negative_length;
                    }
                    else {
                        ++branches_with_length;
                        lengths.addValue( d );
                        if ( d == 0 ) {
                            ++zero_length;
                        }
                    }
                    if ( n.getBranchData().isHasConfidences() ) {
                        for( final Confidence c : n.getBranchData().getConfidences() ) {
                            if ( ( c != null ) && ( c.getValue() >= 0 ) ) {
                                final String type = ForesterUtil.isEmpty( c.getType() ) ? "" : c.getType().trim();
                                DescriptiveStatistics s = support.get( type );
                                if ( s == null ) {
                                    s = new BasicDescriptiveStatistics();
                                    support.put( type, s );
                                }
                                s.addValue( c.getValue() );
                            }
                        }
                    }
                }
                if ( tip ) {
                    if ( nd.isHasTaxonomy() && !nd.getTaxonomy().isEmpty() ) {
                        ++tips_taxonomy;
                        if ( ( nd.getTaxonomy().getIdentifier() != null )
                                && !ForesterUtil.isEmpty( nd.getTaxonomy().getIdentifier().getValue() ) ) {
                            ++tips_tax_id;
                        }
                    }
                    if ( nd.isHasSequence() ) {
                        ++tips_sequence;
                        boolean mol = false;
                        boolean dom = false;
                        for( final Sequence s : nd.getSequences() ) {
                            if ( s == null ) {
                                continue;
                            }
                            mol |= !ForesterUtil.isEmpty( s.getMolecularSequence() );
                            dom |= ( s.getDomainArchitecture() != null )
                                    && ( s.getDomainArchitecture().getNumberOfDomains() > 0 );
                        }
                        if ( mol ) {
                            ++tips_mol_seq;
                        }
                        if ( dom ) {
                            ++tips_domains;
                        }
                    }
                    if ( nd.isHasDistribution() ) {
                        ++tips_distribution;
                    }
                    if ( nd.isHasReference() ) {
                        ++tips_reference;
                    }
                }
                if ( nd.isHasDate() ) {
                    ++nodes_date;
                    if ( tip ) {
                        ++tips_date;
                    }
                }
                if ( nd.isHasEvent() ) {
                    final Event e = nd.getEvent();
                    if ( !e.isUnassigned() ) {
                        ++nodes_events;
                        duplications += Math.max( 0, e.getNumberOfDuplications() ); // the raw counts: a mixed
                        speciations += Math.max( 0, e.getNumberOfSpeciations() ); // event carries several
                        gene_losses += Math.max( 0, e.getNumberOfGeneLosses() );
                    }
                }
                if ( nd.isHasProperties() ) {
                    final Set<String> seen = new java.util.HashSet<>();
                    for( final Property p : nd.getProperties().getProperties() ) {
                        if ( ( p != null ) && !ForesterUtil.isEmpty( p.getRef() ) && seen.add( p.getRef() ) ) {
                            property_refs.merge( p.getRef(), 1, Integer::sum );
                        }
                    }
                }
            }
            taxonomies = AptxUtil.obtainAllDistinctTaxonomies( phy.getRoot() );
        }
    }

    // ------------------------------------------------------------------ groups
    static Group fileGroup( final File file, final boolean edited ) {
        final List<Fact> f = new ArrayList<>();
        String detail = "";
        if ( file == null ) {
            f.add( new Fact( "Location", "not saved to a file yet" ) );
        }
        else {
            String path;
            try {
                path = file.getCanonicalPath();
            }
            catch ( final IOException e ) {
                path = file.getAbsolutePath();
            }
            f.add( new Fact( "Path", path ) );
            f.add( new Fact( "Format", fileFormatLabel( file ) ) );
            if ( file.isFile() ) {
                f.add( new Fact( "Size", humanSize( file.length() ) ) );
                f.add( new Fact( "Modified", DATE.format( Instant.ofEpochMilli( file.lastModified() )
                        .atZone( ZoneId.systemDefault() ) ) ) );
            }
            detail = file.getName();
        }
        f.add( new Fact( "Unsaved changes", edited ? "yes" : "no" ) );
        return new Group( FILE, detail, f, null );
    }

    static Group structureGroup( final Phylogeny phy, final Scan s ) {
        final List<Fact> f = new ArrayList<>();
        final int internal = s.nodes - s.tips;
        f.add( new Fact( "Tips", INT.format( s.tips ) ) );
        f.add( new Fact( "Internal nodes", INT.format( internal ) ) );
        f.add( new Fact( "Total nodes", INT.format( s.nodes ) ) );
        f.add( new Fact( "Branches", INT.format( Math.max( 0, s.nodes - 1 ) ) ) ); // one per non-root node
        f.add( new Fact( "Rooted", phy.isRooted() ? "yes" : "no" ) );
        f.add( new Fact( "Rerootable", phy.isRerootable() ? "yes" : "no" ) );
        final int poly = PhylogenyMethods.countNumberOfPolytomies( phy );
        f.add( new Fact( "Branching", ( poly == 0 ) ? "fully binary"
                : ( INT.format( poly ) + ( poly == 1 ? " polytomy" : " polytomies" ) ) ) );
        f.add( new Fact( "Depth", INT.format( PhylogenyMethods.calculateMaxDepth( phy ) ) + " (root to deepest tip)" ) );
        if ( s.branches_with_length > 0 ) {
            f.add( new Fact( "Height", DEC.format( PhylogenyMethods.calculateMaxDistanceToRoot( phy ) )
                    + " (longest root-to-tip path)" ) );
        }
        if ( s.collapsed > 0 ) {
            f.add( new Fact( "Collapsed clades", INT.format( s.collapsed ) ) );
        }
        return new Group( STRUCTURE, INT.format( s.tips ) + ( s.tips == 1 ? " tip" : " tips" ), f, null );
    }

    /** Null when the tree has no branch lengths at all. */
    static Group branchLengthGroup( final Phylogeny phy, final Scan s ) {
        final int branches = Math.max( 0, s.nodes - 1 );
        if ( ( s.branches_with_length == 0 ) && ( s.negative_length == 0 ) ) {
            return null;
        }
        final List<Fact> f = new ArrayList<>();
        final DescriptiveStatistics st = s.lengths;
        f.add( new Fact( "Branches with lengths", ofText( s.branches_with_length, branches ) ) );
        if ( st.getN() > 0 ) {
            f.add( new Fact( "Median", DEC.format( st.median() ) ) );
            f.add( new Fact( "Mean", DEC.format( st.arithmeticMean() )
                    + ( st.getN() > 1 ? " ± " + DEC.format( st.sampleStandardDeviation() ) + " (sd)" : "" ) ) );
            f.add( new Fact( "Minimum", DEC.format( st.getMin() ) ) );
            f.add( new Fact( "Maximum", DEC.format( st.getMax() ) ) );
            f.add( new Fact( "Sum", DEC.format( st.getSum() ) + " (total tree length)" ) );
        }
        if ( s.zero_length > 0 ) {
            f.add( new Fact( "Zero-length", INT.format( s.zero_length ) ) );
        }
        if ( s.negative_length > 0 ) {
            f.add( new Fact( "Negative", INT.format( s.negative_length ) + " (not drawn to scale)" ) );
        }
        f.add( new Fact( "Ultrametric", AptxUtil.isUltrametric( phy ) ? "yes (all tips equidistant from the root)"
                : "no" ) );
        final String detail = ( st.getN() > 0 ) ? ( "median " + DEC.format( st.median() ) ) : "";
        return new Group( BRANCH_LENGTHS, detail, f, histogram( st ) );
    }

    /** One group per support type ("Support values", or "Support values (bootstrap)" when typed / several). */
    static List<Group> supportGroups( final Scan s ) {
        final List<Group> out = new ArrayList<>();
        for( final Map.Entry<String, DescriptiveStatistics> e : s.support.entrySet() ) {
            final DescriptiveStatistics st = e.getValue();
            final List<Fact> f = new ArrayList<>();
            final String type = e.getKey();
            if ( !type.isEmpty() ) {
                f.add( new Fact( "Type", type ) );
            }
            f.add( new Fact( "Branches with support", INT.format( st.getN() ) ) );
            f.add( new Fact( "Median", DEC.format( st.median() ) ) );
            f.add( new Fact( "Mean", DEC.format( st.arithmeticMean() )
                    + ( st.getN() > 1 ? " ± " + DEC.format( st.sampleStandardDeviation() ) + " (sd)" : "" ) ) );
            f.add( new Fact( "Minimum", DEC.format( st.getMin() ) ) );
            f.add( new Fact( "Maximum", DEC.format( st.getMax() ) ) );
            final String title = ( type.isEmpty() && ( s.support.size() == 1 ) ) ? SUPPORT
                    : ( SUPPORT + " (" + ( type.isEmpty() ? "untyped" : type ) + ")" );
            out.add( new Group( title, "median " + DEC.format( st.median() ), f, histogram( st ) ) );
        }
        return out;
    }

    static Group coverageGroup( final Phylogeny phy, final Scan s ) {
        final List<Fact> f = new ArrayList<>();
        final int t = s.tips;
        f.add( new Fact( "Taxonomy", ofText( s.tips_taxonomy, t ) + " tips" ) );
        if ( s.tips_taxonomy > 0 ) {
            f.add( new Fact( "Taxonomy identifiers", ofText( s.tips_tax_id, t ) + " tips" ) );
            final int distinct = ( s.taxonomies == null ) ? 0 : s.taxonomies.size();
            f.add( new Fact( "Distinct taxonomies", INT.format( distinct ) ) );
        }
        f.add( new Fact( "Sequences", ofText( s.tips_sequence, t ) + " tips" ) );
        if ( s.tips_sequence > 0 ) {
            f.add( new Fact( "Molecular sequences", ofText( s.tips_mol_seq, t ) + " tips" ) );
            f.add( new Fact( "Domain architectures", ofText( s.tips_domains, t ) + " tips" ) );
        }
        f.add( new Fact( "Dates", ofText( s.tips_date, t ) + " tips"
                + ( ( s.nodes_date > s.tips_date ) ? ", " + INT.format( s.nodes_date - s.tips_date ) + " internal nodes"
                        : "" ) ) );
        if ( s.tips_distribution > 0 ) {
            f.add( new Fact( "Distributions", ofText( s.tips_distribution, t ) + " tips" ) );
        }
        if ( s.tips_reference > 0 ) {
            f.add( new Fact( "References", ofText( s.tips_reference, t ) + " tips" ) );
        }
        final int internal = s.nodes - s.tips;
        if ( internal > 0 ) {
            f.add( new Fact( "Named internal nodes", ofText( s.named_internal, internal ) ) );
        }
        if ( s.nodes_events > 0 ) {
            final List<String> parts = new ArrayList<>();
            if ( s.duplications > 0 ) {
                parts.add( INT.format( s.duplications ) + ( s.duplications == 1 ? " duplication" : " duplications" ) );
            }
            if ( s.speciations > 0 ) {
                parts.add( INT.format( s.speciations ) + ( s.speciations == 1 ? " speciation" : " speciations" ) );
            }
            if ( s.gene_losses > 0 ) {
                parts.add( INT.format( s.gene_losses ) + ( s.gene_losses == 1 ? " gene loss" : " gene losses" ) );
            }
            if ( parts.isEmpty() ) {
                parts.add( INT.format( s.nodes_events ) + ( s.nodes_events == 1 ? " node" : " nodes" ) + " with events" );
            }
            f.add( new Fact( "Events", String.join( ", ", parts ) ) );
        }
        if ( s.property_refs.isEmpty() ) {
            f.add( new Fact( "Properties", "none" ) );
        }
        else {
            for( final Map.Entry<String, Integer> e : s.property_refs.entrySet() ) {
                f.add( new Fact( "Property " + e.getKey(), INT.format( e.getValue() )
                        + ( e.getValue() == 1 ? " node" : " nodes" ) ) );
            }
        }
        final int kinds = ( s.tips_taxonomy > 0 ? 1 : 0 ) + ( s.tips_sequence > 0 ? 1 : 0 ) + ( s.nodes_date > 0 ? 1 : 0 )
                + ( s.tips_distribution > 0 ? 1 : 0 ) + ( s.nodes_events > 0 ? 1 : 0 )
                + ( s.property_refs.isEmpty() ? 0 : 1 );
        final String detail = ( kinds == 0 ) ? "no annotations" : ( kinds + ( kinds == 1 ? " kind" : " kinds" ) );
        return new Group( COVERAGE, detail, f, null );
    }

    /** Null when there is no time axis to speak of (no panel info, type Off and no dates). */
    static Group timeAxisGroup( final TimeAxis time, final Scan s ) {
        if ( time == null ) {
            return null;
        }
        final boolean off = ( time.type == null ) || ( time.type == Options.TIME_AXIS_TYPE.NONE );
        if ( off && !time.dated && ( s.nodes_date == 0 ) ) {
            return null;
        }
        final List<Fact> f = new ArrayList<>();
        f.add( new Fact( "Axis", off ? "off" : time.type.toString() ) );
        if ( s.nodes_date > 0 ) {
            f.add( new Fact( "Dated nodes", INT.format( s.nodes_date )
                    + ( time.dated ? " (a dated tree)" : "" ) ) );
        }
        if ( !ForesterUtil.isEmpty( time.unit ) ) {
            f.add( new Fact( "Unit", time.unit ) );
        }
        if ( time.type == Options.TIME_AXIS_TYPE.CALENDAR ) {
            if ( time.present_date > 0 ) {
                f.add( new Fact( "Most recent date", DEC.format( time.present_date ) ) );
            }
        }
        else if ( time.root_age_ma > 0 ) {
            f.add( new Fact( "Root age", DEC.format( time.root_age_ma ) + " Ma" ) );
        }
        return new Group( TIME_AXIS, off ? "off" : time.type.toString(), f, null );
    }

    // ------------------------------------------------------------------ helpers
    /** {@link #HISTOGRAM_BINS} equal-width bins, or null when there are too few / too similar values to bin. */
    static Histogram histogram( final DescriptiveStatistics st ) {
        if ( ( st == null ) || ( st.getN() < 4 ) ) {
            return null;
        }
        final double min = st.getMin();
        final double max = st.getMax();
        if ( ( max - min ) < 0.0001 ) {
            return null;
        }
        final int[] counts = new int[ HISTOGRAM_BINS ];
        final double width = ( max - min ) / HISTOGRAM_BINS;
        for( final double v : st.getDataAsDoubleArray() ) {
            int i = (int) ( ( v - min ) / width );
            if ( i >= HISTOGRAM_BINS ) {
                i = HISTOGRAM_BINS - 1; // the maximum lands in the last (closed) bin
            }
            if ( i < 0 ) {
                i = 0;
            }
            counts[ i ]++;
        }
        return new Histogram( min, max, counts );
    }

    /** "12 of 40 (30%)". */
    static String ofText( final int n, final int of ) {
        if ( of <= 0 ) {
            return INT.format( n );
        }
        final int pct = (int) Math.round( ( 100.0 * n ) / of );
        return INT.format( n ) + " of " + INT.format( of ) + " (" + pct + "%)";
    }

    static String humanSize( final long bytes ) {
        if ( bytes < 1024 ) {
            return bytes + " B";
        }
        final double kb = bytes / 1024.0;
        if ( kb < 1024 ) {
            return oneDecimal( kb ) + " KB";
        }
        final double mb = kb / 1024.0;
        if ( mb < 1024 ) {
            return oneDecimal( mb ) + " MB";
        }
        return oneDecimal( mb / 1024.0 ) + " GB";
    }

    private static String oneDecimal( final double d ) {
        final String s = String.format( Locale.US, "%.1f", d );
        return s.endsWith( ".0" ) ? s.substring( 0, s.length() - 2 ) : s;
    }

    /**
     * A best-effort label of the file's tree format from its first non-blank line (the same sniff the readers use):
     * phyloXML, Nexus, NHX, Newick, JSON, or "unknown". A zipped file is noted. Never throws.
     */
    static String fileFormatLabel( final File file ) {
        if ( ( file == null ) || !file.isFile() ) {
            return "unknown";
        }
        final String lc = file.getName().toLowerCase( Locale.US );
        if ( lc.endsWith( ".zip" ) ) {
            return "zipped (" + formatFromSuffix( lc.substring( 0, lc.length() - 4 ) ) + ")";
        }
        String first = null;
        try ( BufferedReader r = new BufferedReader( new FileReader( file ) ) ) {
            String line;
            int lines = 0;
            while( ( ( line = r.readLine() ) != null ) && ( lines++ < 50 ) ) {
                if ( !line.trim().isEmpty() ) {
                    first = line.trim();
                    break;
                }
            }
        }
        catch ( final IOException | RuntimeException e ) {
            return formatFromSuffix( lc );
        }
        if ( first == null ) {
            return "empty file";
        }
        final String f = first.toLowerCase( Locale.US );
        if ( f.startsWith( "<" ) ) {
            return "phyloXML";
        }
        if ( f.startsWith( "#nexus" ) || f.startsWith( "nexus" ) || f.startsWith( "# nexus" ) || f.startsWith( "begin" ) ) {
            return "Nexus";
        }
        if ( f.startsWith( "{" ) || f.startsWith( "[{" ) ) {
            return "JSON";
        }
        if ( first.contains( "[&&NHX" ) ) {
            return "NHX";
        }
        if ( f.startsWith( "(" ) || f.startsWith( "tree" ) || lc.endsWith( ".nwk" ) || lc.endsWith( ".nh" )
                || lc.endsWith( ".newick" ) || lc.endsWith( ".tre" ) || lc.endsWith( ".tree" ) ) {
            return "Newick";
        }
        return formatFromSuffix( lc );
    }

    private static String formatFromSuffix( final String lc ) {
        if ( lc.endsWith( ".xml" ) || lc.endsWith( ".phyloxml" ) ) {
            return "phyloXML";
        }
        if ( lc.endsWith( ".nex" ) || lc.endsWith( ".nexus" ) || lc.endsWith( ".nxs" ) ) {
            return "Nexus";
        }
        if ( lc.endsWith( ".nhx" ) ) {
            return "NHX";
        }
        if ( lc.endsWith( ".nwk" ) || lc.endsWith( ".nh" ) || lc.endsWith( ".newick" ) || lc.endsWith( ".tre" )
                || lc.endsWith( ".tree" ) ) {
            return "Newick";
        }
        if ( lc.endsWith( ".json" ) ) {
            return "JSON";
        }
        return "unknown";
    }

    /** Integer formatting with thousands separators (for the window's subtitle and tests). */
    static String count( final long n ) {
        return INT.format( n );
    }

    /** Up-to-six-decimals number formatting (no trailing zeros). */
    static String number( final double d ) {
        return DEC.format( d );
    }
}
