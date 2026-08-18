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
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods.NDF;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

/**
 * One selectable field of the redesigned search tool: the thing a query is matched against (a node-data field
 * such as the scientific name, a specific custom {@code <property>} by its {@code ref}, or a numeric built-in
 * such as branch length / support). A field knows whether it is {@link #isNumeric() numeric} (which decides
 * which {@link SearchMode}s apply) and how to extract its value(s) from a {@link PhylogenyNode} -- fields are
 * frequently multi-valued (several synonyms, several properties with the same ref, several domains), and a
 * match against any one value is a match.
 *
 * <p>This is the extensible registry that replaces the fixed set of search checkboxes: adding a future field
 * is one factory here, not a checkbox plus a branch in a monolithic matcher. The {@code ANY_TEXT} field
 * reproduces today's default "search every text field at once" behaviour.
 */
final class SearchField {

    /** The text node-data fields OR-ed together by the synthetic {@link #anyText()} field. Mirrors the fields
     *  the legacy default (all-fields, no properties, no domains) search covers. */
    private static final NDF[]     TEXT_NDFS    = { NDF.NodeName, NDF.TaxonomyCode, NDF.TaxonomyCommonName,
            NDF.TaxonomyScientificName, NDF.TaxonomyIdentifier, NDF.TaxonomySynonym, NDF.TaxonomicLineage,
            NDF.SequenceName, NDF.GeneName, NDF.SequenceSymbol, NDF.SequenceAccession, NDF.Annotation, NDF.CrossRef,
            NDF.BinaryCharacter };
    private static final double[]  NO_NUMBERS   = new double[ 0 ];

    enum Kind {
        ANY_TEXT, NDF_FIELD, PROPERTY, BRANCH_LENGTH, CONFIDENCE, CLADE_SIZE, NUM_CHILDREN, DEPTH, DISTANCE_TO_ROOT,
        NODE_TYPE
    }

    private final String  _label;
    private final boolean _numeric;
    private final Kind    _kind;
    private final NDF     _ndf;          // NDF_FIELD only
    private final String  _property_ref; // PROPERTY only

    private SearchField( final String label, final boolean numeric, final Kind kind, final NDF ndf,
                         final String property_ref ) {
        _label = label;
        _numeric = numeric;
        _kind = kind;
        _ndf = ndf;
        _property_ref = property_ref;
    }

    // ---- factories ------------------------------------------------------------------------------------------

    /** The synthetic "search every text field" field -- the default, reproducing today's all-fields search. */
    static SearchField anyText() {
        return new SearchField( "Any Text", false, Kind.ANY_TEXT, null, null );
    }

    /** A single existing node-data text field (node name, scientific name, sequence name, ...). */
    static SearchField ofNdf( final NDF ndf ) {
        return new SearchField( labelForNdf( ndf ), false, Kind.NDF_FIELD, ndf, null );
    }

    /** A custom phyloXML {@code <property>} identified by its {@code ref} (e.g. {@code aptx:host}); {@code numeric}
     *  decides whether the numeric operators apply (derive it from the property datatype -- see
     *  {@link #datatypeIsNumeric(String)}). */
    static SearchField property( final String ref, final boolean numeric ) {
        return new SearchField( ref, numeric, Kind.PROPERTY, null, ref );
    }

    /** The branch length (distance to parent), a numeric built-in. */
    static SearchField branchLength() {
        return new SearchField( "Branch Length", true, Kind.BRANCH_LENGTH, null, null );
    }

    /** The branch support / confidence value(s), a numeric built-in. */
    static SearchField confidence() {
        return new SearchField( "Support / Confidence", true, Kind.CONFIDENCE, null, null );
    }

    /** Topological: the number of external tips in the node's subtree (a leaf = 1), a numeric built-in. */
    static SearchField cladeSize() {
        return new SearchField( "Structure: Clade Size (tips)", true, Kind.CLADE_SIZE, null, null );
    }

    /** Topological: the number of the node's DIRECT children (a leaf = 0), a numeric built-in -- e.g. > 2 finds
     *  polytomies / unresolved nodes. */
    static SearchField numChildren() {
        return new SearchField( "Structure: Number of Children", true, Kind.NUM_CHILDREN, null, null );
    }

    /** Topological: the number of EDGES from the root to the node (root = 0), a numeric built-in. */
    static SearchField depth() {
        return new SearchField( "Structure: Depth from Root (edges)", true, Kind.DEPTH, null, null );
    }

    /** Metric: the node's distance from the root = the sum of the branch lengths on the path to the root; a numeric
     *  built-in, offered only when the tree has branch lengths. */
    static SearchField distanceToRoot() {
        return new SearchField( "Structure: Distance from Root", true, Kind.DISTANCE_TO_ROOT, null, null );
    }

    /** Categorical: whether a node is a {@code leaf}, an {@code internal} node, or the {@code root} -- a small,
     *  enumerable value set, so it pairs naturally with the value-autocomplete popup (pick "leaf"/"internal"/"root"
     *  from the list rather than knowing to type it). */
    static SearchField nodeType() {
        return new SearchField( "Structure: Node Type", false, Kind.NODE_TYPE, null, null );
    }

    /** The static list of string fields for the search-tool field selector, in menu order ("Any text field"
     *  first = the default). Phase 3 will make this per-tree (only fields the tree actually has, plus one entry
     *  per custom-property ref and the numeric built-ins); until then this fixed set covers every text field. */
    static SearchField[] stringMenuFields() {
        return new SearchField[] { anyText(), ofNdf( NDF.NodeName ), ofNdf( NDF.TaxonomyScientificName ),
                ofNdf( NDF.TaxonomyCommonName ), ofNdf( NDF.TaxonomyCode ), ofNdf( NDF.TaxonomyIdentifier ),
                ofNdf( NDF.TaxonomySynonym ), ofNdf( NDF.TaxonomicLineage ), ofNdf( NDF.SequenceName ),
                ofNdf( NDF.GeneName ), ofNdf( NDF.SequenceSymbol ), ofNdf( NDF.SequenceAccession ), ofNdf( NDF.Domain ),
                ofNdf( NDF.Annotation ), ofNdf( NDF.CrossRef ), ofNdf( NDF.BinaryCharacter ),
                ofNdf( NDF.MolecularSequence ), ofNdf( NDF.Properties ) };
    }

    /** The canonical order of the (non-node-name) text fields in the field selector. */
    private static final NDF[] TEXT_MENU_ORDER = { NDF.TaxonomyScientificName, NDF.TaxonomyCommonName,
            NDF.TaxonomyCode, NDF.TaxonomyIdentifier, NDF.TaxonomySynonym, NDF.TaxonomicLineage, NDF.SequenceName,
            NDF.GeneName, NDF.SequenceSymbol, NDF.SequenceAccession, NDF.Domain, NDF.Annotation, NDF.CrossRef,
            NDF.BinaryCharacter, NDF.MolecularSequence };

    /**
     * The fields to offer for {@code phy}, tailored to what the tree actually carries (one preorder pass): always
     * "Any text field" + "Node name"; then each text field present; then the numeric built-ins Branch length /
     * Support (when present); then one entry per distinct custom-property {@code ref}, typed numeric-or-string by
     * its {@code datatype} (falling back to "do all its values parse as numbers?" when the datatype is absent). This
     * is what makes the selector self-documenting -- the user sees exactly what is searchable in this tree.
     */
    static List<SearchField> availableFields( final Phylogeny phy ) {
        final List<SearchField> out = new ArrayList<>();
        out.add( anyText() );
        out.add( ofNdf( NDF.NodeName ) );
        if ( ( phy == null ) || phy.isEmpty() ) {
            return out;
        }
        final EnumSet<NDF> pending = EnumSet.noneOf( NDF.class );
        for ( final NDF n : TEXT_MENU_ORDER ) {
            pending.add( n );
        }
        final EnumSet<NDF> present = EnumSet.noneOf( NDF.class );
        boolean has_branch_length = false;
        boolean has_confidence = false;
        // per property ref: [saw non-numeric datatype, saw numeric datatype, saw a value, all values numeric so far]
        final Map<String, boolean[]> prop_stats = new LinkedHashMap<>();
        final List<String> tmp = new ArrayList<>();
        for ( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( !pending.isEmpty() ) {
                for ( final Iterator<NDF> pit = pending.iterator(); pit.hasNext(); ) {
                    final NDF n = pit.next();
                    tmp.clear();
                    collectNdfStrings( node, n, tmp );
                    if ( anyNonBlank( tmp ) ) {
                        present.add( n );
                        pit.remove();
                    }
                }
            }
            if ( !has_branch_length && ( node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) ) {
                has_branch_length = true;
            }
            if ( !has_confidence && hasRealConfidence( node ) ) {
                has_confidence = true;
            }
            if ( node.getNodeData().getProperties() != null ) {
                for ( final Property p : node.getNodeData().getProperties().getProperties() ) {
                    // skip internal (aptx:) properties -- they are Archaeopteryx metadata (e.g. the re-import
                    // profile), not user data, so they are not offered as searchable fields
                    if ( ( p == null ) || ForesterUtil.isEmpty( p.getRef() )
                            || TreePanelUtil.isInternalPropertyRef( p.getRef() ) ) {
                        continue;
                    }
                    final boolean[] st = prop_stats.computeIfAbsent( p.getRef(),
                                                                     k -> new boolean[] { false, false, false, true } );
                    if ( !ForesterUtil.isEmpty( p.getDataType() ) ) {
                        if ( datatypeIsNumeric( p.getDataType() ) ) {
                            st[ 1 ] = true;
                        }
                        else {
                            st[ 0 ] = true;
                        }
                    }
                    if ( !ForesterUtil.isEmpty( p.getValue() ) ) {
                        st[ 2 ] = true;
                        if ( SearchMatcher.parseFiniteDouble( p.getValue() ) == null ) {
                            st[ 3 ] = false;
                        }
                    }
                }
            }
        }
        for ( final NDF n : TEXT_MENU_ORDER ) {
            if ( present.contains( n ) ) {
                out.add( ofNdf( n ) );
            }
        }
        if ( has_branch_length ) {
            out.add( branchLength() );
        }
        if ( has_confidence ) {
            out.add( confidence() );
        }
        for ( final Map.Entry<String, boolean[]> e : prop_stats.entrySet() ) {
            final boolean[] st = e.getValue();
            final boolean numeric;
            if ( st[ 0 ] ) {
                numeric = false;                 // a non-numeric datatype was declared -> string (wins over numeric)
            }
            else if ( st[ 1 ] ) {
                numeric = true;                  // a numeric datatype was declared
            }
            else {
                numeric = st[ 2 ] && st[ 3 ];    // no datatype hint -> numeric only if every value parses as a number
            }
            out.add( property( e.getKey(), numeric ) );
        }
        // topological fields -- always available (every non-empty tree has structure), offered last; the metric
        // distance-from-root is offered only when the tree carries branch lengths
        out.add( cladeSize() );
        out.add( numChildren() );
        out.add( depth() );
        if ( has_branch_length ) {
            out.add( distanceToRoot() );
        }
        out.add( nodeType() );
        return out;
    }

    private static boolean anyNonBlank( final List<String> values ) {
        for ( final String s : values ) {
            if ( !ForesterUtil.isEmpty( s ) && !s.trim().isEmpty() ) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasRealConfidence( final PhylogenyNode node ) {
        if ( ( node.getBranchData() == null ) || !node.getBranchData().isHasConfidences() ) {
            return false;
        }
        for ( final Confidence c : node.getBranchData().getConfidences() ) {
            if ( ( c != null ) && ( c.getValue() != Confidence.CONFIDENCE_DEFAULT_VALUE ) ) {
                return true;
            }
        }
        return false;
    }

    // ---- description ----------------------------------------------------------------------------------------

    /** A short, user-facing label for the field selector. */
    String label() {
        return _label;
    }

    /** Whether this field's values are numbers (branch length, support, a numeric property) -- which decides
     *  whether the string or numeric {@link SearchMode}s apply. */
    boolean isNumeric() {
        return _numeric;
    }

    Kind kind() {
        return _kind;
    }

    /** The backing {@link NDF} for an {@link Kind#NDF_FIELD} field, else {@code null}. */
    NDF ndf() {
        return _ndf;
    }

    /** The property {@code ref} for a {@link Kind#PROPERTY} field, else {@code null}. */
    String propertyRef() {
        return _property_ref;
    }

    // ---- extraction -----------------------------------------------------------------------------------------

    /** The string value(s) of this field on {@code node} (empty for a numeric field or when the node has no
     *  such data). Never {@code null}; may contain empty strings, which the matcher treats as no-match. */
    List<String> stringValues( final PhylogenyNode node ) {
        final List<String> out = new ArrayList<>();
        if ( node == null ) {
            return out;
        }
        switch ( _kind ) {
            case ANY_TEXT:
                for ( final NDF ndf : TEXT_NDFS ) {
                    collectNdfStrings( node, ndf, out );
                }
                collectVisiblePropertyStrings( node, out ); // so the default also finds custom annotation values
                break;
            case NDF_FIELD:
                collectNdfStrings( node, _ndf, out );
                break;
            case PROPERTY:
                collectPropertyStrings( node, _property_ref, out );
                break;
            case NODE_TYPE:
                add( out, nodeTypeLabel( node ) ); // "leaf" / "internal" / "root"
                break;
            default:
                break; // numeric fields have no string values
        }
        return out;
    }

    /** The categorical value of the {@link Kind#NODE_TYPE} field: {@code root} for the root, {@code leaf} for an
     *  external non-root node, else {@code internal}. */
    static String nodeTypeLabel( final PhylogenyNode node ) {
        if ( node == null ) {
            return "";
        }
        if ( node.isRoot() ) {
            return "root";
        }
        return node.isExternal() ? "leaf" : "internal";
    }

    /**
     * The distinct, non-blank string values this text field takes across every node of {@code phy}, sorted (the
     * candidate set for the value-autocomplete popup). Empty for a numeric field, the synthetic {@code ANY_TEXT}
     * field (its value set is everything mixed together -- not a meaningful pick list), or an empty/null tree. One
     * preorder pass -- the popup recomputes this on open rather than caching it, so it can never go stale against a
     * tree edit.
     */
    static List<String> distinctValues( final Phylogeny phy, final SearchField field ) {
        if ( ( phy == null ) || phy.isEmpty() || ( field == null ) || field._numeric
                || ( field._kind == Kind.ANY_TEXT ) ) {
            return new ArrayList<>();
        }
        // a molecular sequence is (near-)unique per tip and can be megabytes long -- a value pick-list of full
        // sequences is useless AND expensive to build, so it gets no autocomplete (you still search it with contains).
        if ( ( field._kind == Kind.NDF_FIELD ) && ( field._ndf == NDF.MolecularSequence ) ) {
            return new ArrayList<>();
        }
        final Set<String> set = new TreeSet<>();
        for ( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            for ( final String v : field.stringValues( it.next() ) ) {
                if ( ( v != null ) && !v.trim().isEmpty() ) {
                    set.add( v.trim() );
                }
            }
        }
        return new ArrayList<>( set );
    }

    /** The numeric value(s) of this field on {@code node} (empty for a string field or when the node has no
     *  such data / an unparseable value). Never {@code null}. */
    double[] numericValues( final PhylogenyNode node ) {
        if ( ( node == null ) || !_numeric ) {
            return NO_NUMBERS;
        }
        switch ( _kind ) {
            case BRANCH_LENGTH:
                final double bl = node.getDistanceToParent();
                return ( bl != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) ? new double[] { bl } : NO_NUMBERS;
            case CONFIDENCE:
                return confidenceValues( node );
            case PROPERTY:
                return numbersFromStrings( stringValuesForProperty( node, _property_ref ) );
            case CLADE_SIZE:
                return new double[] { node.getNumberOfExternalNodes() }; // cached count of the subtree's tips
            case NUM_CHILDREN:
                return new double[] { node.getNumberOfDescendants() };   // direct children (a leaf = 0)
            case DEPTH:
                return new double[] { node.calculateDepth() };           // edges from the root (root = 0)
            case DISTANCE_TO_ROOT:
                return new double[] { node.calculateDistanceToRoot() };  // summed branch lengths to the root
            default:
                return NO_NUMBERS;
        }
    }

    /** Whether {@code node} carries a value for this field: a numeric value (numeric field) or a non-blank string
     *  value (text field). Used for the field-aware "Inverse" complement -- the inverse of a field-scoped search
     *  is over the nodes the field APPLIES to, not every data-bearing node. */
    boolean carries( final PhylogenyNode node ) {
        if ( node == null ) {
            return false;
        }
        if ( _numeric ) {
            return numericValues( node ).length > 0;
        }
        for ( final String v : stringValues( node ) ) {
            if ( !ForesterUtil.isEmpty( v ) && !v.trim().isEmpty() ) {
                return true;
            }
        }
        return false;
    }

    // ---- helpers --------------------------------------------------------------------------------------------

    private static double[] confidenceValues( final PhylogenyNode node ) {
        if ( ( node.getBranchData() == null ) || !node.getBranchData().isHasConfidences() ) {
            return NO_NUMBERS;
        }
        final List<Double> vals = new ArrayList<>();
        for ( final Confidence c : node.getBranchData().getConfidences() ) {
            if ( ( c != null ) && ( c.getValue() != Confidence.CONFIDENCE_DEFAULT_VALUE ) ) {
                vals.add( c.getValue() );
            }
        }
        return toArray( vals );
    }

    private static List<String> stringValuesForProperty( final PhylogenyNode node, final String ref ) {
        final List<String> out = new ArrayList<>();
        collectPropertyStrings( node, ref, out );
        return out;
    }

    private static void collectPropertyStrings( final PhylogenyNode node, final String ref, final List<String> out ) {
        if ( ForesterUtil.isEmpty( ref ) || ( node.getNodeData().getProperties() == null ) ) {
            return;
        }
        for ( final Property p : node.getNodeData().getProperties().getProperties() ) {
            if ( ( p != null ) && ref.equals( p.getRef() ) ) {
                add( out, p.getValue() );
            }
        }
    }

    /** Adds the values of the node's USER-VISIBLE custom properties (skipping the internal aptx: namespace, e.g.
     *  the re-import profile), so the default "Any text field" also searches annotation columns like data:host. */
    private static void collectVisiblePropertyStrings( final PhylogenyNode node, final List<String> out ) {
        if ( node.getNodeData().getProperties() == null ) {
            return;
        }
        for ( final Property p : node.getNodeData().getProperties().getProperties() ) {
            if ( ( p != null ) && !TreePanelUtil.isInternalPropertyRef( p.getRef() ) ) {
                add( out, p.getValue() );
            }
        }
    }

    private static void collectNdfStrings( final PhylogenyNode node, final NDF ndf, final List<String> out ) {
        switch ( ndf ) {
            case NodeName:
                add( out, node.getName() );
                break;
            case TaxonomyCode:
                if ( node.getNodeData().isHasTaxonomy() ) {
                    add( out, node.getNodeData().getTaxonomy().getTaxonomyCode() );
                }
                break;
            case TaxonomyCommonName:
                if ( node.getNodeData().isHasTaxonomy() ) {
                    add( out, node.getNodeData().getTaxonomy().getCommonName() );
                }
                break;
            case TaxonomyScientificName:
                if ( node.getNodeData().isHasTaxonomy() ) {
                    add( out, node.getNodeData().getTaxonomy().getScientificName() );
                }
                break;
            case TaxonomyIdentifier:
                if ( node.getNodeData().isHasTaxonomy() && ( node.getNodeData().getTaxonomy().getIdentifier() != null ) ) {
                    add( out, node.getNodeData().getTaxonomy().getIdentifier().getValue() );
                }
                break;
            case TaxonomySynonym:
                if ( node.getNodeData().isHasTaxonomy() ) {
                    for ( final String syn : node.getNodeData().getTaxonomy().getSynonyms() ) {
                        add( out, syn );
                    }
                }
                break;
            case TaxonomicLineage:
                if ( node.getNodeData().isHasTaxonomy() && ( node.getNodeData().getTaxonomy().getLineage() != null ) ) {
                    for ( final String lin : node.getNodeData().getTaxonomy().getLineage() ) {
                        add( out, lin );
                    }
                }
                break;
            case SequenceName:
                if ( node.getNodeData().isHasSequence() ) {
                    add( out, node.getNodeData().getSequence().getName() );
                }
                break;
            case GeneName:
                if ( node.getNodeData().isHasSequence() ) {
                    add( out, node.getNodeData().getSequence().getGeneName() );
                }
                break;
            case SequenceSymbol:
                if ( node.getNodeData().isHasSequence() ) {
                    add( out, node.getNodeData().getSequence().getSymbol() );
                }
                break;
            case SequenceAccession:
                if ( node.getNodeData().isHasSequence() && ( node.getNodeData().getSequence().getAccession() != null ) ) {
                    add( out, node.getNodeData().getSequence().getAccession().getValue() );
                }
                break;
            case MolecularSequence:
                if ( node.getNodeData().isHasSequence() ) {
                    add( out, node.getNodeData().getSequence().getMolecularSequence() );
                }
                break;
            case Domain:
                if ( node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getDomainArchitecture() != null ) ) {
                    final DomainArchitecture da = node.getNodeData().getSequence().getDomainArchitecture();
                    for ( int i = 0; i < da.getNumberOfDomains(); ++i ) {
                        add( out, da.getDomain( i ).getName() );
                    }
                }
                break;
            case Annotation:
                if ( node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getAnnotations() != null ) ) {
                    for ( final Annotation ann : node.getNodeData().getSequence().getAnnotations() ) {
                        add( out, ann.getDesc() );
                        add( out, ann.getRef() );
                    }
                }
                break;
            case CrossRef:
                if ( node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getCrossReferences() != null ) ) {
                    for ( final Accession x : node.getNodeData().getSequence().getCrossReferences() ) {
                        add( out, x.getComment() );
                        add( out, x.getSource() );
                        add( out, x.getValue() );
                    }
                }
                break;
            case BinaryCharacter:
                if ( node.getNodeData().getBinaryCharacters() != null ) {
                    for ( final String c : node.getNodeData().getBinaryCharacters().getPresentCharacters() ) {
                        add( out, c );
                    }
                    for ( final String c : node.getNodeData().getBinaryCharacters().getGainedCharacters() ) {
                        add( out, c );
                    }
                }
                break;
            case Properties:
                if ( node.getNodeData().getProperties() != null ) {
                    for ( final Property p : node.getNodeData().getProperties().getProperties() ) {
                        if ( p != null ) {
                            add( out, p.getValue() );
                        }
                    }
                }
                break;
            default:
                break;
        }
    }

    private static void add( final List<String> out, final String s ) {
        if ( s != null ) {
            out.add( s );
        }
    }

    private static double[] numbersFromStrings( final List<String> strings ) {
        final List<Double> vals = new ArrayList<>();
        for ( final String s : strings ) {
            final Double d = SearchMatcher.parseFiniteDouble( s );
            if ( d != null ) {
                vals.add( d );
            }
        }
        return toArray( vals );
    }

    private static double[] toArray( final List<Double> vals ) {
        if ( vals.isEmpty() ) {
            return NO_NUMBERS;
        }
        final double[] a = new double[ vals.size() ];
        for ( int i = 0; i < a.length; ++i ) {
            a[ i ] = vals.get( i );
        }
        return a;
    }

    /** The local names of the numeric XSD datatypes (the part after the namespace colon), lower-cased. */
    private static final Set<String> NUMERIC_XSD_TYPES = new HashSet<>( Arrays.asList( "decimal", "double", "float",
            "integer", "int", "long", "short", "byte", "unsignedint", "unsignedlong", "unsignedshort", "unsignedbyte",
            "nonnegativeinteger", "nonpositiveinteger", "negativeinteger", "positiveinteger" ) );

    /** Whether a phyloXML property {@code datatype} (e.g. {@code xsd:decimal}, {@code xsd:integer}) denotes a
     *  number -- used to offer numeric operators for a numeric custom-property field. Matches the datatype's local
     *  name against the known numeric XSD types (so a custom datatype that merely CONTAINS "int" etc., such as
     *  {@code x:footprint}, is NOT misclassified). An empty/unknown datatype is non-numeric (the caller may still
     *  fall back to inspecting the actual values). */
    static boolean datatypeIsNumeric( final String datatype ) {
        if ( ForesterUtil.isEmpty( datatype ) ) {
            return false;
        }
        final String d = datatype.toLowerCase( Locale.ROOT );
        final int colon = d.lastIndexOf( ':' );
        return NUMERIC_XSD_TYPES.contains( ( colon >= 0 ) ? d.substring( colon + 1 ) : d );
    }

    private static String labelForNdf( final NDF ndf ) {
        // Title Case, no colon groups -- consistent with the "Display Data" checkbox labels (e.g. "Taxonomy
        // Scientific", "Seq Name", "Gene Name"). The search-only Structure fields (see the factories) keep a
        // "Structure:" prefix since they have no Display-Data counterpart.
        switch ( ndf ) {
            case NodeName:
                return "Node Name";
            case TaxonomyCode:
                return "Taxonomy Code";
            case TaxonomyCommonName:
                return "Taxonomy Common";
            case TaxonomyScientificName:
                return "Taxonomy Scientific";
            case TaxonomyIdentifier:
                return "Taxonomy Identifier";
            case TaxonomySynonym:
                return "Taxonomy Synonym";
            case TaxonomicLineage:
                return "Taxonomy Lineage";
            case SequenceName:
                return "Seq Name";
            case GeneName:
                return "Gene Name";
            case SequenceSymbol:
                return "Gene Symbol"; // the phyloXML <sequence><symbol>, i.e. the (species-specific) gene symbol
            case SequenceAccession:
                return "Seq Accession";
            case MolecularSequence:
                return "Molecular Sequence";
            case Domain:
                return "Domain";
            case Annotation:
                return "Annotation";
            case CrossRef:
                return "Cross-Reference";
            case BinaryCharacter:
                return "Binary Character";
            case Properties:
                return "Any Property";
            default:
                return ndf.toString();
        }
    }
}
