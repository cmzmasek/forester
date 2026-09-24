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
import java.util.List;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterUtil;

/**
 * What the hover card says about a node, in the same wording as the Archaeopteryx.js tooltip (so the two viewers
 * read alike): first the facts about the node itself -- name, distance to parent, date, distribution, depth and,
 * for an internal node, the tips below it -- then the branch confidences, then one HEADED block per taxonomy /
 * sequence / events / properties. A {@link Row} with a null value is a section heading.
 * <p>
 * Every row either stands in that leading block or under a heading, and that is the point of the order: a row
 * printed after a headed section is read as part of it, so a stray "Distribution" or "Tips below" trailing the
 * SEQUENCE block looked like a sequence field, and the properties -- headed by nothing -- looked like more of
 * whatever came before them. Pure.
 */
final class NodeHoverText {

    /** One line of the card: a key/value pair, or (value == null) a section heading. */
    static final class Row {

        final String key;
        final String value;

        private Row( final String key, final String value ) {
            this.key = key;
            this.value = value;
        }

        static Row heading( final String text ) {
            return new Row( text, null );
        }

        static Row line( final String key, final String value ) {
            return new Row( key, value );
        }

        boolean isHeading() {
            return value == null;
        }

        @Override
        public String toString() {
            return isHeading() ? "[" + key + "]" : key + ": " + value;
        }
    }

    /** Pfam accessions are "PF" plus at least five digits -- five today, and {5,} will not truncate a sixth. */
    private static final java.util.regex.Pattern PFAM_ACC = java.util.regex.Pattern
            .compile( "(PF\\d{5,})(\\.\\d+)?" );

    /**
     * The rollover for a protein DOMAIN: what it is, how good the hit is, where it sits, and which tip's
     * protein it belongs to -- the architectures line up in a column of their own, far from the labels, so
     * the tip has to be named here just as it is for an annotation-column cell.
     */
    static List<Row> domainRows( final PhylogenyNode node, final ProteinDomain d, final int protein_length ) {
        final List<Row> rows = new ArrayList<Row>();
        if ( d == null ) {
            return rows;
        }
        // The domain's name is DATA, so it is a value row: a heading is drawn small-caps uppercase and would
        // misrepresent a name that is case-sensitive ("Bcl-2", "wnt").
        rows.add( Row.line( "Domain", ForesterUtil.isEmpty( d.getName() ) ? "(unnamed)" : d.getName() ) );
        if ( d.getConfidence() >= 0 ) { // ProteinDomain.asText's own test for "has a confidence"
            rows.add( Row.line( "E-value", formatEValue( d.getConfidence() ) ) );
        }
        final int len = d.getLength(); // not a second copy of the 1-based inclusive convention
        rows.add( Row.line( "Residues", d.getFrom() + "\u2013" + d.getTo() + " (" + len + " aa)" ) );
        if ( protein_length > 0 ) {
            rows.add( Row.line( "Protein length", protein_length + " aa" ) );
        }
        rows.add( Row.line( "Tip", AnnotationColumns.tipLabel( node ) ) );
        // The row is "Accession", not "Pfam": phyloXML's domain id is whatever the annotation source used,
        // and a SMART or CDD accession called "Pfam" would be wrong in the same way the 404 below is. It is
        // shown whatever it is -- dropping it because it is not Pfam would lose a real identifier -- while
        // only a Pfam accession earns the entry LINK. (Archaeopteryx.js draws the same distinction.)
        if ( !ForesterUtil.isEmpty( d.getId() ) ) {
            rows.add( Row.line( "Accession", d.getId().trim() ) );
        }
        final String pfam = pfamAccession( d );
        // Say WHICH promise the click makes. A Pfam accession addresses an entry; a Pfam IDENTIFIER does not
        // -- measured against InterPro, /entry/pfam/NB-ARC/ is a 404 while /entry/pfam/PF00931/ is a 200 --
        // so a domain that carries only its name can be looked up but not linked to, and a user told
        // "the Pfam entry" who lands on a list of search results has been misled.
        if ( pfam != null ) {
            rows.add( Row.line( "Click", "the Pfam entry" ) );
        }
        else if ( !ForesterUtil.isEmpty( d.getName() ) ) {
            rows.add( Row.line( "Click", "to look this domain up at InterPro" ) );
        }
        return rows;
    }

    /** Whether a domain can be linked at all, and which of the two routes applies. */
    static boolean hasLink( final ProteinDomain d ) {
        return ( d != null ) && ( ( pfamAccession( d ) != null ) || !ForesterUtil.isEmpty( d.getName() ) );
    }

    /**
     * The Pfam accession of a domain, or null. phyloXML puts it in the {@code <domain>} element's {@code id}
     * attribute; a file that names its domains but carries no accession (which is the common case) gets none,
     * and the rollover then simply has no Pfam row.
     */
    static String pfamAccession( final ProteinDomain d ) {
        if ( ( d == null ) || ForesterUtil.isEmpty( d.getId() ) ) {
            return null;
        }
        // ANCHORED, not a substring search: an id that merely CONTAINS an accession -- "SM00109_PF00018",
        // a pipe-joined composite -- is not a Pfam accession, and linking it to that entry would send the
        // user to a different family. A trailing ".25" version suffix is allowed, nothing else is.
        final java.util.regex.Matcher m = PFAM_ACC.matcher( d.getId().trim().toUpperCase( java.util.Locale.US ) );
        return m.matches() ? m.group( 1 ) : null;
    }

    /**
     * An E-value reads as a power of ten: 1.2e-40 rather than 0.00000000000000000000000000000000000000012.
     * Values at or above 0.001 are written plainly, since that is where they stop being exponents.
     */
    // Per call, deliberately. DecimalFormat is NOT thread-safe, and this is a package-private static on a
    // stateless class -- the suite already calls its neighbours off the EDT, and two threads inside format()
    // share its internal DigitList and return interleaved digits rather than failing. A rollover formats one
    // number per pointer move, so there is nothing here worth a shared instance or a ThreadLocal.
    private static java.text.DecimalFormat format( final String pattern ) {
        return new java.text.DecimalFormat( pattern,
                                            java.text.DecimalFormatSymbols.getInstance( java.util.Locale.US ) );
    }

    static String formatEValue( final double v ) {
        if ( v == 0 ) {
            return "0";
        }
        return ( ( Math.abs( v ) >= 0.001 ) && ( Math.abs( v ) < 1000 ) ) ? format( "0.###" ).format( v )
                                                                          : format( "0.##E0" ).format( v );
    }

    private NodeHoverText() {
        // not instantiable
    }

    /** The card's rows for {@code n}; empty when there is nothing worth a card. */
    static List<Row> rows( final PhylogenyNode n ) {
        return rows( n, false );
    }

    /**
     * The card's rows. With {@code root_free} (a tree declared unrooted, shown unrooted --
     * {@link Rerooting#hidesRootDependentValues}) nothing that depends on where the tree is stored as rooted is shown:
     * an internal node lists the tips on each of its sides instead of its distance to parent, depth and tips below,
     * and a tip shows its (single) branch length but no depth.
     */
    static List<Row> rows( final PhylogenyNode n, final boolean root_free ) {
        final List<Row> rows = new ArrayList<>();
        final NodeData nd = n.getNodeData();
        if ( !ForesterUtil.isEmpty( n.getName() ) ) {
            rows.add( Row.line( "Name", n.getName() ) );
        }
        if ( n.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            if ( !root_free ) {
                rows.add( Row.line( "Distance to parent", NodeDataDraft.formatNumber( n.getDistanceToParent() ) ) );
            }
            else if ( n.isExternal() ) {
                rows.add( Row.line( "Branch length", NodeDataDraft.formatNumber( n.getDistanceToParent() ) ) );
            }
        }
        if ( nd.isHasDate() ) {
            final String date = dateText( nd.getDate() );
            if ( !date.isEmpty() ) {
                rows.add( Row.line( "Date", date ) );
            }
        }
        if ( nd.isHasDistribution() && !ForesterUtil.isEmpty( nd.getDistribution().getDesc() ) ) {
            rows.add( Row.line( "Distribution", nd.getDistribution().getDesc() ) );
        }
        if ( !root_free ) {
            rows.add( Row.line( "Depth", String.valueOf( n.calculateDepth() ) ) );
            if ( !n.isExternal() ) {
                rows.add( Row.line( "Tips below", String.valueOf( countTips( n ) ) ) );
            }
        }
        else if ( !n.isExternal() ) {
            rows.add( Row.line( "Tips around", Rerooting.tipsAroundText( n ) ) );
        }
        if ( n.getBranchData().isHasConfidences() ) {
            for( final Confidence c : n.getBranchData().getConfidences() ) {
                if ( c.getValue() == Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                    continue;
                }
                final String key = ForesterUtil.isEmpty( c.getType() ) || "unknown".equals( c.getType() )
                        ? "Confidence" : "Confidence [" + c.getType() + "]";
                rows.add( Row.line( key, NodeDataDraft.formatNumber( c.getValue() ) ) );
                if ( ( c.getStandardDeviation() != Confidence.CONFIDENCE_DEFAULT_VALUE )
                        && ( c.getStandardDeviation() != 0 ) ) {
                    rows.add( Row.line( "stdev", NodeDataDraft.formatNumber( c.getStandardDeviation() ) ) );
                }
            }
        }
        if ( nd.isHasTaxonomy() ) {
            for( final Taxonomy t : nd.getTaxonomies() ) {
                if ( t == null ) {
                    continue;
                }
                final int at = rows.size();
                if ( ( t.getIdentifier() != null ) && !ForesterUtil.isEmpty( t.getIdentifier().getValue() ) ) {
                    final String p = t.getIdentifier().getProvider();
                    rows.add( Row.line( ForesterUtil.isEmpty( p ) ? "Id" : "Id [" + p + "]",
                                        t.getIdentifier().getValue() ) );
                }
                add( rows, "Code", t.getTaxonomyCode() );
                add( rows, "Scientific name", t.getScientificName() );
                add( rows, "Common name", t.getCommonName() );
                add( rows, "Rank", t.getRank() );
                if ( rows.size() > at ) {
                    rows.add( at, Row.heading( "Taxonomy" ) );
                }
            }
        }
        if ( nd.isHasSequence() ) {
            for( final Sequence s : nd.getSequences() ) {
                if ( s == null ) {
                    continue;
                }
                final int at = rows.size();
                if ( ( s.getAccession() != null ) && !ForesterUtil.isEmpty( s.getAccession().getValue() ) ) {
                    final String src = s.getAccession().getSource();
                    rows.add( Row.line( ForesterUtil.isEmpty( src ) ? "Accession" : "Accession [" + src + "]",
                                        s.getAccession().getValue() ) );
                    add( rows, "comment", s.getAccession().getComment() );
                }
                add( rows, "Symbol", s.getSymbol() );
                add( rows, "Name", s.getName() );
                add( rows, "Gene name", s.getGeneName() );
                add( rows, "Location", s.getLocation() );
                add( rows, "Type", s.getType() );
                if ( rows.size() > at ) {
                    rows.add( at, Row.heading( "Sequence" ) );
                }
            }
        }
        if ( nd.isHasEvent() ) {
            final Event ev = nd.getEvent();
            final int at = rows.size();
            if ( !ev.isUnassigned() && ( ev.getEventType() != Event.EventType.mixed )
                    && ( ev.getEventType() != Event.EventType.unassigned ) ) {
                rows.add( Row.line( "Type", ev.getEventType().toString() ) );
            }
            if ( ev.getNumberOfDuplications() > 0 ) {
                rows.add( Row.line( "Duplications", String.valueOf( ev.getNumberOfDuplications() ) ) );
            }
            if ( ev.getNumberOfSpeciations() > 0 ) {
                rows.add( Row.line( "Speciations", String.valueOf( ev.getNumberOfSpeciations() ) ) );
            }
            if ( ev.getNumberOfGeneLosses() > 0 ) {
                rows.add( Row.line( "Losses", String.valueOf( ev.getNumberOfGeneLosses() ) ) );
            }
            if ( rows.size() > at ) {
                rows.add( at, Row.heading( "Events" ) );
            }
        }
        if ( nd.isHasProperties() ) {
            final int at = rows.size();
            for( final Property p : nd.getProperties().getProperties() ) {
                if ( TreePanelUtil.isInternalPropertyRef( p.getRef() ) || TreePanelUtil.isVisualStylePropertyRef( p.getRef() )
                        || ForesterUtil.isEmpty( p.getValue() ) ) {
                    continue;
                }
                String ref = p.getRef();
                if ( ref.indexOf( ':' ) > 0 ) {
                    ref = ref.substring( ref.indexOf( ':' ) + 1 );
                }
                rows.add( Row.line( ref, ForesterUtil.isEmpty( p.getUnit() ) ? p.getValue()
                        : p.getValue() + " " + p.getUnit() ) );
            }
            if ( rows.size() > at ) {
                rows.add( at, Row.heading( "Properties" ) );
            }
        }
        // a card that would only say "Depth: n" (a bare unnamed node) is noise: show nothing
        if ( ( rows.size() == 1 ) && "Depth".equals( rows.get( 0 ).key ) ) {
            rows.clear();
        }
        return rows;
    }

    /** "6.5 [5 - 8] mya", "6.5 mya", "[5 - 8]", or the description alone when there is no number. */
    static String dateText( final Date d ) {
        final StringBuilder sb = new StringBuilder();
        if ( d.getValue() != null ) {
            sb.append( d.getValue().toPlainString() );
        }
        // only a real WIDTH is a range: an exactly dated tip whose bounds differ by float noise (TreeAnnotator writes
        // {9.0,9.000000000000004}) would otherwise read "11.15 [11.149999999999999 - 11.150000000000002]" -- the same
        // claim the overlays used to draw. One predicate for what is drawn and what is said.
        if ( AptxUtil.hasDateIntervalWidth( d ) ) {
            if ( sb.length() > 0 ) {
                sb.append( ' ' );
            }
            sb.append( '[' ).append( d.getMin().toPlainString() ).append( " - " ).append( d.getMax().toPlainString() )
                    .append( ']' );
        }
        if ( ( sb.length() > 0 ) && !ForesterUtil.isEmpty( d.getUnit() ) ) {
            sb.append( ' ' ).append( d.getUnit() );
        }
        if ( !ForesterUtil.isEmpty( d.getDesc() ) ) {
            if ( sb.length() > 0 ) {
                sb.append( " (" ).append( d.getDesc() ).append( ')' );
            }
            else {
                sb.append( d.getDesc() );
            }
        }
        return sb.toString();
    }

    /** The tips under {@code n}, counted directly (the node's cached count can be stale after an edit). */
    static int countTips( final PhylogenyNode n ) {
        if ( n.isExternal() ) {
            return 1;
        }
        int count = 0;
        final java.util.ArrayDeque<PhylogenyNode> stack = new java.util.ArrayDeque<>();
        stack.push( n );
        while ( !stack.isEmpty() ) {
            final PhylogenyNode x = stack.pop();
            if ( x.isExternal() ) {
                ++count;
            }
            else {
                for( int i = 0; i < x.getNumberOfDescendants(); ++i ) {
                    stack.push( x.getChildNode( i ) );
                }
            }
        }
        return count;
    }

    private static void add( final List<Row> rows, final String key, final String value ) {
        if ( !ForesterUtil.isEmpty( value ) ) {
            rows.add( Row.line( key, value ) );
        }
    }
}
