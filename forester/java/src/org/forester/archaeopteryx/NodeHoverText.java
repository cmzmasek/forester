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

    private NodeHoverText() {
        // not instantiable
    }

    /** The card's rows for {@code n}; empty when there is nothing worth a card. */
    static List<Row> rows( final PhylogenyNode n ) {
        final List<Row> rows = new ArrayList<>();
        final NodeData nd = n.getNodeData();
        if ( !ForesterUtil.isEmpty( n.getName() ) ) {
            rows.add( Row.line( "Name", n.getName() ) );
        }
        if ( n.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            rows.add( Row.line( "Distance to parent", NodeDataDraft.formatNumber( n.getDistanceToParent() ) ) );
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
        rows.add( Row.line( "Depth", String.valueOf( n.calculateDepth() ) ) );
        if ( !n.isExternal() ) {
            rows.add( Row.line( "Tips below", String.valueOf( countTips( n ) ) ) );
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
        if ( ( d.getMin() != null ) && ( d.getMax() != null ) ) {
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
