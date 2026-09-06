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

import java.math.BigDecimal;
import java.net.URI;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.Set;
import java.util.regex.Pattern;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Point;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Reference;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;
import org.forester.util.ForesterUtil;
import org.forester.util.TaxonomyUtil;

/**
 * A plain, mutable, text-only draft of everything the node editor can change on one {@link PhylogenyNode}: the
 * view-model behind {@link NodeDataForm}. Reading a node into a draft ({@link #from}) has NO side effect on the
 * node; validation ({@link #validate}) reports human-readable problems keyed by field so the form can outline the
 * offending widget; and writing ({@link #writeTo}) is a DIFF against the draft the editor opened with -- only fields
 * the user actually changed are applied, so an untouched branch length keeps its full precision and untouched
 * sub-elements (sequence annotations, taxonomy lineage, distribution polygons, ...) survive a round trip intact.
 * <p>
 * Everything is a {@code String} exactly as it would appear in a text field ("" = absent), which makes dirtiness a
 * plain {@link #equals} against the baseline. Multi-instance data (confidences, sequences, properties) are lists of
 * small drafts; a sequence draft remembers the {@link Sequence} it came from so a write can mutate that object in
 * place and keep the extras the form does not edit.
 */
final class NodeDataDraft {

    // ---- field keys (a Problem names the field it is about with one of these; the form maps keys to widgets) ----
    static final String NAME                  = "name";
    static final String BRANCH_LENGTH         = "branch_length";
    static final String BRANCH_WIDTH          = "branch_width";
    static final String TAX_ID                = "taxonomy.id";
    static final String TAX_PROVIDER          = "taxonomy.provider";
    static final String TAX_CODE              = "taxonomy.code";
    static final String TAX_SCI_NAME          = "taxonomy.scientific_name";
    static final String TAX_AUTHORITY         = "taxonomy.authority";
    static final String TAX_COMMON_NAME       = "taxonomy.common_name";
    static final String TAX_SYNONYMS          = "taxonomy.synonyms";
    static final String TAX_RANK              = "taxonomy.rank";
    static final String TAX_URIS              = "taxonomy.uris";
    static final String EV_DUPLICATIONS       = "events.duplications";
    static final String EV_SPECIATIONS        = "events.speciations";
    static final String EV_GENE_LOSSES        = "events.gene_losses";
    static final String DATE_DESC             = "date.desc";
    static final String DATE_VALUE            = "date.value";
    static final String DATE_MIN              = "date.min";
    static final String DATE_MAX              = "date.max";
    static final String DATE_UNIT             = "date.unit";
    static final String DIST_DESC             = "distribution.desc";
    static final String DIST_DATUM            = "distribution.datum";
    static final String DIST_LAT              = "distribution.latitude";
    static final String DIST_LONG             = "distribution.longitude";
    static final String DIST_ALT              = "distribution.altitude";
    static final String DIST_ALT_UNIT         = "distribution.altitude_unit";
    static final String REF_DESC              = "reference.desc";
    static final String REF_DOI               = "reference.doi";
    /** Indexed keys: {@code confidence.<i>.value}, {@code sequence.<i>.symbol}, {@code property.<i>.ref}, ... */
    static String confidenceKey( final int i, final String field ) {
        return "confidence." + i + "." + field;
    }

    static String sequenceKey( final int i, final String field ) {
        return "sequence." + i + "." + field;
    }

    static String propertyKey( final int i, final String field ) {
        return "property." + i + "." + field;
    }

    static final String CONF_VALUE   = "value";
    static final String CONF_TYPE    = "type";
    static final String CONF_SD      = "sd";
    static final String SEQ_NAME     = "name";
    static final String SEQ_SYMBOL   = "symbol";
    static final String SEQ_GENE     = "gene_name";
    static final String SEQ_TYPE     = "type";
    static final String SEQ_ACC      = "accession";
    static final String SEQ_SOURCE   = "source";
    static final String SEQ_LOCATION = "location";
    static final String SEQ_MOL_SEQ  = "mol_seq";
    static final String SEQ_URIS     = "uris";
    static final String SEQ_ALIGNED  = "aligned";
    static final String PROP_REF     = "ref";
    static final String PROP_VALUE   = "value";
    static final String PROP_UNIT    = "unit";
    static final String PROP_TYPE    = "datatype";

    // ---- section names (used for the provenance sentence and the changed-section report) ----
    static final String SEC_BASIC        = "Basic";
    static final String SEC_TAXONOMY     = "Taxonomy";
    static final String SEC_SEQUENCES    = "Sequences";
    static final String SEC_EVENTS       = "Events";
    static final String SEC_DATE         = "Date";
    static final String SEC_DISTRIBUTION = "Distribution";
    static final String SEC_REFERENCE    = "Reference";
    static final String SEC_PROPERTIES   = "Properties";

    /** Confidence types offered in the editor's type combo (phyloXML leaves the vocabulary open; these are the usual). */
    static final List<String> CONFIDENCE_TYPES = Collections
            .unmodifiableList( Arrays.asList( "bootstrap", "probability", "posterior", "sh_alrt", "unknown" ) );
    /** Common xsd datatypes for a property, offered in the editor's datatype combo (any {@code prefix:name} is legal). */
    static final List<String> PROPERTY_DATATYPES = Collections.unmodifiableList( Arrays
            .asList( "xsd:string", "xsd:decimal", "xsd:integer", "xsd:int", "xsd:boolean", "xsd:date", "xsd:anyURI" ) );
    /** Which of the datatypes above carry a NUMBER (their value is validated as one). */
    private static final Set<String> NUMERIC_DATATYPES = new LinkedHashSet<>( Arrays
            .asList( "xsd:decimal", "xsd:integer", "xsd:int", "xsd:long", "xsd:short", "xsd:float", "xsd:double",
                     "xsd:nonNegativeInteger", "xsd:positiveInteger" ) );
    private static final String       DEFAULT_CONFIDENCE_TYPE = "unknown";
    private static final Pattern      LINE_BREAK              = Pattern.compile( "\\R" );
    /** What a molecular sequence may contain: letters plus the gap / terminator / unknown characters. */
    private static final Pattern      NOT_RESIDUE             = Pattern.compile( "[^A-Za-z\\-\\.\\*\\?]" );
    private static final DecimalFormat NUMBER_FORMAT;
    static {
        // Up to 9 decimals (what the phyloXML writer keeps), no trailing zeros, ALWAYS a leading zero ("0.5", not ".5").
        NUMBER_FORMAT = new DecimalFormat( "0.#########", DecimalFormatSymbols.getInstance( Locale.US ) );
    }

    /** One validation problem: which field ({@code key}, see the constants above) and what is wrong with it. */
    static final class Problem {

        final String key;
        final String message;

        Problem( final String key, final String message ) {
            this.key = key;
            this.message = message;
        }

        @Override
        public String toString() {
            return key + ": " + message;
        }
    }

    /** One {@code <confidence>} row. */
    static final class ConfidenceDraft {

        String value = "";
        String type  = "";
        String sd    = "";

        ConfidenceDraft() {
        }

        ConfidenceDraft( final String value, final String type, final String sd ) {
            this.value = nn( value );
            this.type = nn( type );
            this.sd = nn( sd );
        }

        boolean isBlank() {
            return value.isEmpty() && type.isEmpty() && sd.isEmpty();
        }

        ConfidenceDraft copy() {
            return new ConfidenceDraft( value, type, sd );
        }

        @Override
        public boolean equals( final Object o ) {
            if ( !( o instanceof ConfidenceDraft ) ) {
                return false;
            }
            final ConfidenceDraft c = (ConfidenceDraft) o;
            return value.equals( c.value ) && type.equals( c.type ) && sd.equals( c.sd );
        }

        @Override
        public int hashCode() {
            return Objects.hash( value, type, sd );
        }
    }

    /** One {@code <sequence>} card. {@link #origin} is the live object it was read from (null for a new card) --
     *  a write mutates it in place, so annotations, domain architecture, cross references etc. are kept. */
    static final class SequenceDraft {

        String         name     = "";
        String         symbol   = "";
        String         geneName = "";
        String         type     = "";
        String         accession = "";
        String         source   = "";
        String         location = "";
        String         molSeq   = "";
        boolean        aligned  = false;
        /** One URI per line. */
        String         uris     = "";
        /** The live object this was read from, or -- after a write created one -- written to; null for a new card. */
        Sequence       origin;

        SequenceDraft() {
            origin = null;
        }

        private SequenceDraft( final Sequence origin ) {
            this.origin = origin;
        }

        /** An EMPTY draft bound to {@code origin} (may be null) -- the identity for an in-place write, no values. */
        static SequenceDraft withOrigin( final Sequence origin ) {
            return new SequenceDraft( origin );
        }

        static SequenceDraft from( final Sequence s ) {
            final SequenceDraft d = new SequenceDraft( s );
            d.name = nn( s.getName() );
            d.symbol = nn( s.getSymbol() );
            d.geneName = nn( s.getGeneName() );
            d.type = nn( s.getType() );
            if ( s.getAccession() != null ) {
                d.accession = nn( s.getAccession().getValue() );
                d.source = nn( s.getAccession().getSource() );
            }
            d.location = nn( s.getLocation() );
            d.molSeq = nn( s.getMolecularSequence() );
            d.aligned = s.isMolecularSequenceAligned();
            d.uris = urisToLines( s.getUris() );
            return d;
        }

        boolean isBlank() {
            return name.isEmpty() && symbol.isEmpty() && geneName.isEmpty() && type.isEmpty() && accession.isEmpty()
                    && source.isEmpty() && location.isEmpty() && molSeq.isEmpty() && uris.isEmpty();
        }

        SequenceDraft copy() {
            final SequenceDraft d = new SequenceDraft( origin );
            d.name = name;
            d.symbol = symbol;
            d.geneName = geneName;
            d.type = type;
            d.accession = accession;
            d.source = source;
            d.location = location;
            d.molSeq = molSeq;
            d.aligned = aligned;
            d.uris = uris;
            return d;
        }

        /** Field equality only -- the origin object is deliberately NOT compared (it is provenance, not data). */
        @Override
        public boolean equals( final Object o ) {
            if ( !( o instanceof SequenceDraft ) ) {
                return false;
            }
            final SequenceDraft s = (SequenceDraft) o;
            return name.equals( s.name ) && symbol.equals( s.symbol ) && geneName.equals( s.geneName )
                    && type.equals( s.type ) && accession.equals( s.accession ) && source.equals( s.source )
                    && location.equals( s.location ) && molSeq.equals( s.molSeq ) && ( aligned == s.aligned )
                    && uris.equals( s.uris );
        }

        @Override
        public int hashCode() {
            return Objects.hash( name, symbol, geneName, type, accession, source, location, molSeq, aligned, uris );
        }
    }

    /** One {@code <property>} row. */
    static final class PropertyDraft {

        String    ref       = "";
        String    value     = "";
        String    unit      = "";
        String    datatype  = "xsd:string";
        AppliesTo appliesTo = AppliesTo.NODE;
        String    idRef     = "";

        PropertyDraft() {
        }

        PropertyDraft( final String ref, final String value, final String unit, final String datatype,
                       final AppliesTo applies_to ) {
            this.ref = nn( ref );
            this.value = nn( value );
            this.unit = nn( unit );
            this.datatype = nn( datatype );
            this.appliesTo = ( applies_to == null ) ? AppliesTo.NODE : applies_to;
        }

        static PropertyDraft from( final Property p ) {
            final PropertyDraft d = new PropertyDraft( p.getRef(), p.getValue(), p.getUnit(), p.getDataType(),
                                                       p.getAppliesTo() );
            d.idRef = nn( p.getIdRef() );
            return d;
        }

        boolean isBlank() {
            return ref.isEmpty() && value.isEmpty() && unit.isEmpty();
        }

        PropertyDraft copy() {
            final PropertyDraft d = new PropertyDraft( ref, value, unit, datatype, appliesTo );
            d.idRef = idRef;
            return d;
        }

        @Override
        public boolean equals( final Object o ) {
            if ( !( o instanceof PropertyDraft ) ) {
                return false;
            }
            final PropertyDraft p = (PropertyDraft) o;
            return ref.equals( p.ref ) && value.equals( p.value ) && unit.equals( p.unit )
                    && datatype.equals( p.datatype ) && ( appliesTo == p.appliesTo ) && idRef.equals( p.idRef );
        }

        @Override
        public int hashCode() {
            return Objects.hash( ref, value, unit, datatype, appliesTo, idRef );
        }
    }

    // ---- Basic ----
    String                      name         = "";
    String                      branchLength = "";
    /** "" = default width (1). */
    String                      branchWidth  = "";
    final List<ConfidenceDraft> confidences  = new ArrayList<>();
    // ---- Taxonomy ----
    String                      taxId          = "";
    String                      taxProvider    = "";
    String                      taxCode        = "";
    String                      taxSciName     = "";
    String                      taxAuthority   = "";
    String                      taxCommonName  = "";
    /** One synonym per line. */
    String                      taxSynonyms    = "";
    String                      taxRank        = "";
    /** One URI per line. */
    String                      taxUris        = "";
    // ---- Sequences ----
    final List<SequenceDraft>   sequences      = new ArrayList<>();
    // ---- Events (internal nodes only; "" = not set) ----
    String                      duplications   = "";
    String                      speciations    = "";
    String                      geneLosses     = "";
    // ---- Date ----
    String                      dateDesc       = "";
    String                      dateValue      = "";
    String                      dateMin        = "";
    String                      dateMax        = "";
    String                      dateUnit       = "";
    // ---- Distribution (description + the FIRST point) ----
    String                      distDesc       = "";
    String                      distDatum      = "";
    String                      distLat        = "";
    String                      distLong       = "";
    String                      distAlt        = "";
    String                      distAltUnit    = "";
    // ---- Reference ----
    String                      refDesc        = "";
    String                      refDoi         = "";
    // ---- Properties ----
    final List<PropertyDraft>   properties     = new ArrayList<>();

    NodeDataDraft() {
    }

    /** Reads {@code node} into a fresh draft. Pure: the node is not touched (no empty sub-elements are created). */
    static NodeDataDraft from( final PhylogenyNode node ) {
        final NodeDataDraft d = new NodeDataDraft();
        d.name = nn( node.getName() );
        if ( node.getDistanceToParent() != PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) {
            d.branchLength = formatNumber( node.getDistanceToParent() );
        }
        final BranchWidth bw = node.getBranchData().getBranchWidth();
        if ( ( bw != null ) && ( bw.getValue() != BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE ) ) {
            d.branchWidth = formatNumber( bw.getValue() );
        }
        if ( node.getBranchData().isHasConfidences() ) {
            for( final Confidence c : node.getBranchData().getConfidences() ) {
                if ( c.getValue() == Confidence.CONFIDENCE_DEFAULT_VALUE ) {
                    continue;
                }
                final String sd = ( c.getStandardDeviation() != Confidence.CONFIDENCE_DEFAULT_VALUE )
                        && ( c.getStandardDeviation() != 0 ) ? formatNumber( c.getStandardDeviation() ) : "";
                d.confidences.add( new ConfidenceDraft( formatNumber( c.getValue() ), c.getType(), sd ) );
            }
        }
        final NodeData nd = node.getNodeData();
        if ( nd.isHasTaxonomy() ) {
            final Taxonomy t = nd.getTaxonomy();
            if ( t.getIdentifier() != null ) {
                d.taxId = nn( t.getIdentifier().getValue() );
                d.taxProvider = nn( t.getIdentifier().getProvider() );
            }
            d.taxCode = nn( t.getTaxonomyCode() );
            d.taxSciName = nn( t.getScientificName() );
            d.taxAuthority = nn( t.getAuthority() );
            d.taxCommonName = nn( t.getCommonName() );
            d.taxSynonyms = linesOf( t.getSynonyms() );
            d.taxRank = nn( t.getRank() );
            d.taxUris = urisToLines( t.getUris() );
        }
        if ( nd.isHasSequence() ) {
            for( final Sequence s : nd.getSequences() ) {
                if ( s != null ) {
                    d.sequences.add( SequenceDraft.from( s ) );
                }
            }
        }
        if ( !node.isExternal() && nd.isHasEvent() ) {
            final Event e = nd.getEvent();
            d.duplications = countText( e.getNumberOfDuplications() );
            d.speciations = countText( e.getNumberOfSpeciations() );
            d.geneLosses = countText( e.getNumberOfGeneLosses() );
        }
        if ( nd.isHasDate() ) {
            final Date date = nd.getDate();
            d.dateDesc = nn( date.getDesc() );
            d.dateValue = decimalText( date.getValue() );
            d.dateMin = decimalText( date.getMin() );
            d.dateMax = decimalText( date.getMax() );
            d.dateUnit = nn( date.getUnit() );
        }
        if ( nd.isHasDistribution() ) {
            final Distribution dist = nd.getDistribution();
            d.distDesc = nn( dist.getDesc() );
            final Point p = firstPoint( dist );
            if ( p != null ) {
                d.distDatum = Point.UNKNOWN_GEODETIC_DATUM.equals( p.getGeodeticDatum() ) ? ""
                        : nn( p.getGeodeticDatum() );
                d.distLat = decimalText( p.getLatitude() );
                d.distLong = decimalText( p.getLongitude() );
                d.distAlt = decimalText( p.getAltitude() );
                d.distAltUnit = nn( p.getAltiudeUnit() );
            }
        }
        if ( nd.isHasReference() ) {
            d.refDesc = nn( nd.getReference().getDescription() );
            d.refDoi = nn( nd.getReference().getDoi() );
        }
        if ( nd.isHasProperties() ) {
            for( final Property p : nd.getProperties().getProperties() ) {
                d.properties.add( PropertyDraft.from( p ) );
            }
        }
        return d;
    }

    /** A deep copy (the baseline the editor compares against). */
    NodeDataDraft copy() {
        final NodeDataDraft d = new NodeDataDraft();
        d.name = name;
        d.branchLength = branchLength;
        d.branchWidth = branchWidth;
        for( final ConfidenceDraft c : confidences ) {
            d.confidences.add( c.copy() );
        }
        d.taxId = taxId;
        d.taxProvider = taxProvider;
        d.taxCode = taxCode;
        d.taxSciName = taxSciName;
        d.taxAuthority = taxAuthority;
        d.taxCommonName = taxCommonName;
        d.taxSynonyms = taxSynonyms;
        d.taxRank = taxRank;
        d.taxUris = taxUris;
        for( final SequenceDraft s : sequences ) {
            d.sequences.add( s.copy() );
        }
        d.duplications = duplications;
        d.speciations = speciations;
        d.geneLosses = geneLosses;
        d.dateDesc = dateDesc;
        d.dateValue = dateValue;
        d.dateMin = dateMin;
        d.dateMax = dateMax;
        d.dateUnit = dateUnit;
        d.distDesc = distDesc;
        d.distDatum = distDatum;
        d.distLat = distLat;
        d.distLong = distLong;
        d.distAlt = distAlt;
        d.distAltUnit = distAltUnit;
        d.refDesc = refDesc;
        d.refDoi = refDoi;
        for( final PropertyDraft p : properties ) {
            d.properties.add( p.copy() );
        }
        return d;
    }

    // ---- per-section emptiness (drives "start collapsed" in the form and "remove the element" on write) ----
    boolean hasTaxonomy() {
        return !( taxId.isEmpty() && taxProvider.isEmpty() && taxCode.isEmpty() && taxSciName.isEmpty()
                && taxAuthority.isEmpty() && taxCommonName.isEmpty() && taxSynonyms.isEmpty() && taxRank.isEmpty()
                && taxUris.isEmpty() );
    }

    boolean hasEvents() {
        return !( duplications.isEmpty() && speciations.isEmpty() && geneLosses.isEmpty() );
    }

    boolean hasDate() {
        return !( dateDesc.isEmpty() && dateValue.isEmpty() && dateMin.isEmpty() && dateMax.isEmpty()
                && dateUnit.isEmpty() );
    }

    boolean hasDistribution() {
        return !( distDesc.isEmpty() && !hasDistributionPoint() );
    }

    boolean hasDistributionPoint() {
        return !( distDatum.isEmpty() && distLat.isEmpty() && distLong.isEmpty() && distAlt.isEmpty()
                && distAltUnit.isEmpty() );
    }

    boolean hasReference() {
        return !( refDesc.isEmpty() && refDoi.isEmpty() );
    }

    // ---- validation ----
    /**
     * Every problem with the current values, in field order; empty when the draft can be written. {@code internal}
     * says whether the node is internal (events are only validated -- and only written -- for internal nodes).
     */
    List<Problem> validate( final boolean internal ) {
        final List<Problem> ps = new ArrayList<>();
        checkDouble( ps, BRANCH_LENGTH, branchLength, "Branch length must be a number" );
        if ( checkDouble( ps, BRANCH_WIDTH, branchWidth, "Branch width must be a number" ) && !branchWidth.isEmpty()
                && ( Double.parseDouble( branchWidth.trim() ) < 0 ) ) {
            ps.add( new Problem( BRANCH_WIDTH, "Branch width cannot be negative" ) );
        }
        for( int i = 0; i < confidences.size(); ++i ) {
            final ConfidenceDraft c = confidences.get( i );
            if ( c.isBlank() ) {
                continue;
            }
            if ( c.value.trim().isEmpty() ) {
                ps.add( new Problem( confidenceKey( i, CONF_VALUE ), "Confidence " + ( i + 1 ) + " needs a value" ) );
            }
            else {
                checkDouble( ps, confidenceKey( i, CONF_VALUE ), c.value, "Confidence " + ( i + 1 )
                        + " must be a number" );
            }
            if ( checkDouble( ps, confidenceKey( i, CONF_SD ), c.sd, "Confidence " + ( i + 1 )
                    + ": the standard deviation must be a number" ) && !c.sd.trim().isEmpty()
                    && ( Double.parseDouble( c.sd.trim() ) < 0 ) ) {
                ps.add( new Problem( confidenceKey( i, CONF_SD ), "Confidence " + ( i + 1 )
                        + ": the standard deviation cannot be negative" ) );
            }
        }
        if ( !taxCode.trim().isEmpty() && !PhyloXmlUtil.TAXOMONY_CODE_PATTERN.matcher( taxCode.trim() ).matches() ) {
            ps.add( new Problem( TAX_CODE, "Taxonomy code must be 3 to 5 upper-case letters or digits (e.g. HUMAN)" ) );
        }
        if ( !taxRank.trim().isEmpty()
                && !TaxonomyUtil.TAXONOMY_RANKS_SET.contains( taxRank.trim().toLowerCase( Locale.ROOT ) ) ) {
            ps.add( new Problem( TAX_RANK, "\"" + taxRank.trim() + "\" is not a phyloXML taxonomic rank" ) );
        }
        checkUris( ps, TAX_URIS, taxUris, "Taxonomy URI" );
        for( int i = 0; i < sequences.size(); ++i ) {
            final SequenceDraft s = sequences.get( i );
            final String label = "Sequence " + ( i + 1 );
            if ( !s.symbol.trim().isEmpty()
                    && !PhyloXmlUtil.SEQUENCE_SYMBOL_PATTERN.matcher( s.symbol.trim() ).matches() ) {
                ps.add( new Problem( sequenceKey( i, SEQ_SYMBOL ), label
                        + ": a symbol is 1 to 20 characters without whitespace" ) );
            }
            if ( !s.type.trim().isEmpty()
                    && !PhyloXmlUtil.SEQUENCE_TYPES.contains( s.type.trim().toLowerCase( Locale.ROOT ) ) ) {
                ps.add( new Problem( sequenceKey( i, SEQ_TYPE ), label + ": type must be dna, rna or protein" ) );
            }
            checkUris( ps, sequenceKey( i, SEQ_URIS ), s.uris, label + " URI" );
        }
        if ( internal ) {
            checkCount( ps, EV_DUPLICATIONS, duplications, "Duplications" );
            checkCount( ps, EV_SPECIATIONS, speciations, "Speciations" );
            checkCount( ps, EV_GENE_LOSSES, geneLosses, "Gene losses" );
        }
        checkDecimal( ps, DATE_VALUE, dateValue, "Date value" );
        checkDecimal( ps, DATE_MIN, dateMin, "Date min" );
        checkDecimal( ps, DATE_MAX, dateMax, "Date max" );
        if ( checkDecimal( ps, DIST_LAT, distLat, "Latitude" ) && !distLat.trim().isEmpty() ) {
            final BigDecimal v = new BigDecimal( distLat.trim() );
            if ( ( v.compareTo( BigDecimal.valueOf( -90 ) ) < 0 ) || ( v.compareTo( BigDecimal.valueOf( 90 ) ) > 0 ) ) {
                ps.add( new Problem( DIST_LAT, "Latitude must be between -90 and 90" ) );
            }
        }
        if ( checkDecimal( ps, DIST_LONG, distLong, "Longitude" ) && !distLong.trim().isEmpty() ) {
            final BigDecimal v = new BigDecimal( distLong.trim() );
            if ( ( v.compareTo( BigDecimal.valueOf( -180 ) ) < 0 )
                    || ( v.compareTo( BigDecimal.valueOf( 180 ) ) > 0 ) ) {
                ps.add( new Problem( DIST_LONG, "Longitude must be between -180 and 180" ) );
            }
        }
        if ( checkDecimal( ps, DIST_ALT, distAlt, "Altitude" ) && !distAlt.trim().isEmpty()
                && distAltUnit.trim().isEmpty() ) {
            ps.add( new Problem( DIST_ALT_UNIT, "An altitude needs a unit (e.g. m)" ) );
        }
        if ( !refDoi.trim().isEmpty() && !PhyloXmlUtil.LIT_REF_DOI_PATTERN.matcher( refDoi.trim() ).matches() ) {
            ps.add( new Problem( REF_DOI, "That does not look like a DOI (e.g. 10.1093/bioinformatics/btq243)" ) );
        }
        for( int i = 0; i < properties.size(); ++i ) {
            final PropertyDraft p = properties.get( i );
            if ( p.isBlank() ) {
                continue;
            }
            final String label = "Property " + ( i + 1 );
            if ( !isPrefixed( p.ref ) ) {
                ps.add( new Problem( propertyKey( i, PROP_REF ), label
                        + ": the reference needs a namespace prefix, like \"data:depth\"" ) );
            }
            if ( !p.unit.trim().isEmpty() && !isPrefixed( p.unit ) ) {
                ps.add( new Problem( propertyKey( i, PROP_UNIT ), label
                        + ": a unit needs a namespace prefix, like \"METRIC:m\"" ) );
            }
            if ( !isPrefixed( p.datatype ) ) {
                ps.add( new Problem( propertyKey( i, PROP_TYPE ), label
                        + ": the datatype needs a namespace prefix, like \"xsd:string\"" ) );
            }
            else if ( NUMERIC_DATATYPES.contains( p.datatype.trim() ) && !p.value.trim().isEmpty()
                    && !isDecimal( p.value ) ) {
                ps.add( new Problem( propertyKey( i, PROP_VALUE ), label + ": \"" + p.value.trim()
                        + "\" is not a number, but the datatype is " + p.datatype.trim() ) );
            }
        }
        return ps;
    }

    // ---- writing ----
    /**
     * Applies this draft to {@code node} as a DIFF against {@code baseline} (the draft the editor opened with, or
     * last wrote): a field equal to its baseline value is left alone, so untouched values keep their precision and
     * untouched sub-elements are not rebuilt. Sections that end up entirely empty are REMOVED from the node rather
     * than left as empty shells. The draft must validate first ({@link #validate}); a draft that does not throws
     * {@link IllegalStateException} before touching the node.
     *
     * @return the names of the sections that were changed (empty when nothing differed)
     */
    Set<String> writeTo( final PhylogenyNode node, final NodeDataDraft baseline ) {
        final List<Problem> problems = validate( !node.isExternal() );
        if ( !problems.isEmpty() ) {
            throw new IllegalStateException( "draft does not validate: " + problems.get( 0 ) );
        }
        // ONE computation of "what changed" (the same one the provenance sentence is built from) drives the write
        final Set<String> changed = changedSections( baseline, !node.isExternal() );
        final NodeData nd = node.getNodeData();
        if ( changed.contains( SEC_BASIC ) ) {
            writeBasics( node, baseline );
        }
        if ( changed.contains( SEC_TAXONOMY ) ) {
            writeTaxonomy( nd );
        }
        if ( changed.contains( SEC_SEQUENCES ) ) {
            writeSequences( nd );
        }
        if ( changed.contains( SEC_EVENTS ) ) {
            writeEvents( nd );
        }
        if ( changed.contains( SEC_DATE ) ) {
            writeDate( nd );
        }
        if ( changed.contains( SEC_DISTRIBUTION ) ) {
            writeDistribution( nd );
        }
        if ( changed.contains( SEC_REFERENCE ) ) {
            writeReference( nd );
        }
        if ( changed.contains( SEC_PROPERTIES ) ) {
            writeProperties( nd );
        }
        return changed;
    }

    /** The Basic section is written FIELD by field (an untouched branch length must keep its full precision). */
    private void writeBasics( final PhylogenyNode node, final NodeDataDraft baseline ) {
        if ( !name.equals( baseline.name ) ) {
            node.setName( name.trim() );
        }
        if ( !branchLength.equals( baseline.branchLength ) ) {
            node.setDistanceToParent( branchLength.trim().isEmpty() ? PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT
                    : Double.parseDouble( branchLength.trim() ) );
        }
        if ( !branchWidth.equals( baseline.branchWidth ) ) {
            final String w = branchWidth.trim();
            if ( w.isEmpty() || ( Double.parseDouble( w ) == BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE ) ) {
                node.getBranchData().setBranchWidth( null );
            }
            else {
                node.getBranchData().setBranchWidth( new BranchWidth( Double.parseDouble( w ) ) );
            }
        }
        if ( !confidencesEqual( baseline ) ) {
            final List<Confidence> list = node.getBranchData().getConfidences();
            list.clear();
            for( final ConfidenceDraft c : confidences ) {
                if ( c.isBlank() ) {
                    continue;
                }
                final String type = c.type.trim().isEmpty() ? DEFAULT_CONFIDENCE_TYPE : c.type.trim();
                final double sd = c.sd.trim().isEmpty() ? Confidence.CONFIDENCE_DEFAULT_VALUE : Double
                        .parseDouble( c.sd.trim() );
                list.add( new Confidence( Double.parseDouble( c.value.trim() ), type, sd ) );
            }
        }
    }

    /** The sections in which this draft differs from {@code baseline} (what {@link #writeTo} would change). */
    Set<String> changedSections( final NodeDataDraft baseline, final boolean internal ) {
        final Set<String> changed = new LinkedHashSet<>();
        if ( !name.equals( baseline.name ) || !branchLength.equals( baseline.branchLength )
                || !branchWidth.equals( baseline.branchWidth ) || !confidencesEqual( baseline ) ) {
            changed.add( SEC_BASIC );
        }
        if ( !taxonomyEquals( baseline ) ) {
            changed.add( SEC_TAXONOMY );
        }
        if ( !sequencesEqual( baseline ) ) {
            changed.add( SEC_SEQUENCES );
        }
        if ( internal && !eventsEquals( baseline ) ) {
            changed.add( SEC_EVENTS );
        }
        if ( !dateEquals( baseline ) ) {
            changed.add( SEC_DATE );
        }
        if ( !distributionEquals( baseline ) ) {
            changed.add( SEC_DISTRIBUTION );
        }
        if ( !referenceEquals( baseline ) ) {
            changed.add( SEC_REFERENCE );
        }
        if ( !propertiesEqual( baseline ) ) {
            changed.add( SEC_PROPERTIES );
        }
        return changed;
    }

    /** The provenance sentence appended to the tree description after a write (per the repo rule). Pure. */
    static String provenance( final String node_label, final Collection<String> sections ) {
        final StringBuilder sb = new StringBuilder( "Manually edited the node data (" );
        sb.append( String.join( ", ", sections ) ).append( ") of node \"" ).append( node_label ).append( "\"." );
        return sb.toString();
    }

    private void writeTaxonomy( final NodeData nd ) {
        if ( !hasTaxonomy() ) {
            // every field cleared -> drop the element, unless it still carries a lineage the form does not edit
            if ( nd.isHasTaxonomy() && ForesterUtil.isEmpty( nd.getTaxonomy().getLineage() ) ) {
                nd.setTaxonomies( null );
            }
            else if ( nd.isHasTaxonomy() ) {
                applyTaxonomyFields( nd.getTaxonomy() );
            }
            return;
        }
        if ( !nd.isHasTaxonomy() ) {
            nd.setTaxonomy( new Taxonomy() );
        }
        applyTaxonomyFields( nd.getTaxonomy() );
    }

    private void applyTaxonomyFields( final Taxonomy t ) {
        if ( taxId.trim().isEmpty() && taxProvider.trim().isEmpty() ) {
            t.setIdentifier( null );
        }
        else {
            t.setIdentifier( new Identifier( taxId.trim(), taxProvider.trim() ) );
        }
        try {
            t.setTaxonomyCode( taxCode.trim() );
            t.setRank( taxRank.trim().toLowerCase( Locale.ROOT ) );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            throw new IllegalStateException( "validate first: " + e.getMessage(), e );
        }
        t.setScientificName( taxSciName.trim() );
        t.setAuthority( taxAuthority.trim() );
        t.setCommonName( taxCommonName.trim() );
        t.getSynonyms().clear();
        t.getSynonyms().addAll( lines( taxSynonyms ) );
        t.setUris( linesToUris( taxUris, t.getUris() ) );
    }

    private void writeSequences( final NodeData nd ) {
        final List<Sequence> out = new ArrayList<>();
        for( final SequenceDraft d : sequences ) {
            final Sequence s = ( d.origin != null ) ? d.origin : new Sequence();
            if ( d.isBlank() && ( d.origin == null ) ) {
                continue; // an untouched, empty "+ Add sequence" card
            }
            if ( ( d.accession.trim().isEmpty() && d.source.trim().isEmpty() ) ) {
                s.setAccession( null );
            }
            else {
                final String comment = ( s.getAccession() != null ) ? s.getAccession().getComment() : null;
                s.setAccession( new Accession( d.accession.trim(), d.source.trim(), comment ) );
            }
            s.setName( d.name.trim() );
            try {
                s.setSymbol( d.symbol.trim() );
                s.setType( d.type.trim().toLowerCase( Locale.ROOT ) );
            }
            catch ( final PhyloXmlDataFormatException e ) {
                throw new IllegalStateException( "validate first: " + e.getMessage(), e );
            }
            s.setGeneName( d.geneName.trim() );
            s.setLocation( d.location.trim() );
            s.setMolecularSequence( sanitizeMolSeq( d.molSeq ) );
            s.setMolecularSequenceAligned( d.aligned );
            s.setUris( linesToUris( d.uris, s.getUris() ) );
            if ( d.isBlank() && s.isEmpty() ) {
                d.origin = null;
                continue; // an existing sequence the user emptied out, with nothing else left on it -> gone
            }
            d.origin = s;
            out.add( s );
        }
        nd.setSequences( out.isEmpty() ? null : out );
    }

    private void writeEvents( final NodeData nd ) {
        if ( !hasEvents() ) {
            if ( nd.isHasEvent() ) {
                // Every count cleared. An event whose type only ever reflected its counts ("mixed"/"unassigned")
                // and that carries no confidence has nothing left -> remove it. Otherwise rebuild it WITHOUT counts
                // (the count setters would stamp the type to "mixed", so a typed event cannot just be zeroed).
                final Event e = nd.getEvent();
                final boolean count_only = ( e.getEventType() == Event.EventType.mixed )
                        || ( e.getEventType() == Event.EventType.unassigned );
                if ( count_only && ( e.getConfidence() == null ) ) {
                    nd.setEvent( null );
                }
                else {
                    final Event bare = count_only ? new Event() : new Event( e.getEventType() );
                    bare.setConfidence( e.getConfidence() );
                    nd.setEvent( bare );
                }
            }
            return;
        }
        if ( !nd.isHasEvent() ) {
            nd.setEvent( new Event() );
        }
        final Event e = nd.getEvent();
        e.setDuplications( parseCount( duplications ) );
        e.setSpeciations( parseCount( speciations ) );
        e.setGeneLosses( parseCount( geneLosses ) );
    }

    private void writeDate( final NodeData nd ) {
        if ( !hasDate() ) {
            nd.setDate( null );
            return;
        }
        if ( nd.getDate() == null ) {
            nd.setDate( new Date() );
        }
        final Date date = nd.getDate();
        date.setDesc( dateDesc.trim() );
        date.setValue( decimalOrNull( dateValue ) );
        date.setMin( decimalOrNull( dateMin ) );
        date.setMax( decimalOrNull( dateMax ) );
        date.setUnit( dateUnit.trim() );
    }

    private void writeDistribution( final NodeData nd ) {
        final Distribution old = nd.isHasDistribution() ? nd.getDistribution() : null;
        final List<Point> points = new ArrayList<>();
        if ( ( old != null ) && ( old.getPoints() != null ) ) {
            points.addAll( old.getPoints() );
        }
        if ( hasDistributionPoint() ) {
            final Point p = new Point( distDatum.trim().isEmpty() ? Point.UNKNOWN_GEODETIC_DATUM : distDatum.trim(),
                                       decimalOrNull( distLat ),
                                       decimalOrNull( distLong ),
                                       decimalOrNull( distAlt ),
                                       distAltUnit.trim() );
            if ( points.isEmpty() ) {
                points.add( p );
            }
            else {
                points.set( 0, p );
            }
        }
        else if ( !points.isEmpty() ) {
            points.remove( 0 );
        }
        final boolean has_polygons = ( old != null ) && !ForesterUtil.isEmpty( old.getPolygons() );
        if ( distDesc.trim().isEmpty() && points.isEmpty() && !has_polygons ) {
            nd.setDistributions( null );
            return;
        }
        nd.setDistribution( new Distribution( distDesc.trim(), points, ( old != null ) ? old.getPolygons() : null ) );
    }

    private void writeReference( final NodeData nd ) {
        if ( !hasReference() ) {
            nd.setReferences( null );
            return;
        }
        if ( !nd.isHasReference() ) {
            nd.setReference( new Reference( "" ) );
        }
        nd.getReference().setValue( refDesc.trim() );
        try {
            nd.getReference().setDoi( refDoi.trim() );
        }
        catch ( final PhyloXmlDataFormatException e ) {
            throw new IllegalStateException( "validate first: " + e.getMessage(), e );
        }
    }

    private void writeProperties( final NodeData nd ) {
        final PropertiesList list = new PropertiesList();
        for( final PropertyDraft p : properties ) {
            if ( p.isBlank() ) {
                continue;
            }
            list.addProperty( new Property( p.ref.trim(), p.value.trim(), p.unit.trim(), p.datatype.trim(),
                                            p.appliesTo, p.idRef ) );
        }
        nd.setProperties( list.size() == 0 ? null : list );
    }

    // ---- the rows that COUNT: an untouched "+ Add" row/card is not data and must not make the draft dirty ----
    private List<ConfidenceDraft> activeConfidences() {
        final List<ConfidenceDraft> out = new ArrayList<>();
        for( final ConfidenceDraft c : confidences ) {
            if ( !c.isBlank() ) {
                out.add( c );
            }
        }
        return out;
    }

    private List<SequenceDraft> activeSequences() {
        final List<SequenceDraft> out = new ArrayList<>();
        for( final SequenceDraft d : sequences ) {
            if ( !( d.isBlank() && ( d.origin == null ) ) ) {
                out.add( d );
            }
        }
        return out;
    }

    private List<PropertyDraft> activeProperties() {
        final List<PropertyDraft> out = new ArrayList<>();
        for( final PropertyDraft p : properties ) {
            if ( !p.isBlank() ) {
                out.add( p );
            }
        }
        return out;
    }

    private boolean confidencesEqual( final NodeDataDraft b ) {
        return activeConfidences().equals( b.activeConfidences() );
    }

    private boolean sequencesEqual( final NodeDataDraft b ) {
        return activeSequences().equals( b.activeSequences() );
    }

    private boolean propertiesEqual( final NodeDataDraft b ) {
        return activeProperties().equals( b.activeProperties() );
    }

    // ---- section equality against a baseline ----
    private boolean taxonomyEquals( final NodeDataDraft b ) {
        return taxId.equals( b.taxId ) && taxProvider.equals( b.taxProvider ) && taxCode.equals( b.taxCode )
                && taxSciName.equals( b.taxSciName ) && taxAuthority.equals( b.taxAuthority )
                && taxCommonName.equals( b.taxCommonName ) && taxSynonyms.equals( b.taxSynonyms )
                && taxRank.equals( b.taxRank ) && taxUris.equals( b.taxUris );
    }

    private boolean eventsEquals( final NodeDataDraft b ) {
        return duplications.equals( b.duplications ) && speciations.equals( b.speciations )
                && geneLosses.equals( b.geneLosses );
    }

    private boolean dateEquals( final NodeDataDraft b ) {
        return dateDesc.equals( b.dateDesc ) && dateValue.equals( b.dateValue ) && dateMin.equals( b.dateMin )
                && dateMax.equals( b.dateMax ) && dateUnit.equals( b.dateUnit );
    }

    private boolean distributionEquals( final NodeDataDraft b ) {
        return distDesc.equals( b.distDesc ) && distDatum.equals( b.distDatum ) && distLat.equals( b.distLat )
                && distLong.equals( b.distLong ) && distAlt.equals( b.distAlt ) && distAltUnit.equals( b.distAltUnit );
    }

    private boolean referenceEquals( final NodeDataDraft b ) {
        return refDesc.equals( b.refDesc ) && refDoi.equals( b.refDoi );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( !( o instanceof NodeDataDraft ) ) {
            return false;
        }
        final NodeDataDraft b = (NodeDataDraft) o;
        return name.equals( b.name ) && branchLength.equals( b.branchLength ) && branchWidth.equals( b.branchWidth )
                && confidencesEqual( b ) && taxonomyEquals( b ) && sequencesEqual( b ) && eventsEquals( b )
                && dateEquals( b ) && distributionEquals( b ) && referenceEquals( b ) && propertiesEqual( b );
    }

    @Override
    public int hashCode() {
        return Objects.hash( name, branchLength, branchWidth, activeConfidences(), taxSciName, activeSequences(),
                             dateValue, distDesc, refDoi, activeProperties() );
    }

    // ---- text helpers ----
    static String nn( final String s ) {
        return ( s == null ) ? "" : s;
    }

    static String formatNumber( final double d ) {
        return NUMBER_FORMAT.format( d );
    }

    private static String countText( final int count ) {
        return ( count >= 0 ) ? String.valueOf( count ) : "";
    }

    private static String decimalText( final BigDecimal d ) {
        return ( d == null ) ? "" : d.toPlainString();
    }

    private static BigDecimal decimalOrNull( final String s ) {
        return s.trim().isEmpty() ? null : new BigDecimal( s.trim() );
    }

    private static int parseCount( final String s ) {
        return s.trim().isEmpty() ? Event.DEFAULT_VALUE : Integer.parseInt( s.trim() );
    }

    /** The non-blank, trimmed lines of a one-item-per-line text. */
    static List<String> lines( final String text ) {
        final List<String> out = new ArrayList<>();
        if ( text != null ) {
            for( final String line : LINE_BREAK.split( text ) ) {
                if ( !line.trim().isEmpty() ) {
                    out.add( line.trim() );
                }
            }
        }
        return out;
    }

    private static String linesOf( final List<String> items ) {
        if ( ForesterUtil.isEmpty( items ) ) {
            return "";
        }
        final List<String> non_empty = new ArrayList<>();
        for( final String s : items ) {
            if ( !ForesterUtil.isEmpty( s ) ) {
                non_empty.add( s );
            }
        }
        return String.join( "\n", non_empty );
    }

    private static String urisToLines( final List<Uri> uris ) {
        if ( ForesterUtil.isEmpty( uris ) ) {
            return "";
        }
        final List<String> out = new ArrayList<>();
        for( final Uri u : uris ) {
            if ( ( u != null ) && ( u.getValue() != null ) ) {
                out.add( u.getValue().toString() );
            }
        }
        return String.join( "\n", out );
    }

    /** Builds the URI list for the lines, REUSING an existing {@link Uri} whose value matches (so its description
     *  and type attributes, which the form does not edit, survive). Null when there are no lines. */
    private static List<Uri> linesToUris( final String text, final List<Uri> existing ) {
        final List<String> ls = lines( text );
        if ( ls.isEmpty() ) {
            return null;
        }
        final List<Uri> out = new ArrayList<>();
        for( final String line : ls ) {
            Uri reuse = null;
            if ( existing != null ) {
                for( final Uri u : existing ) {
                    if ( ( u != null ) && ( u.getValue() != null ) && u.getValue().toString().equals( line ) ) {
                        reuse = u;
                        break;
                    }
                }
            }
            try {
                out.add( ( reuse != null ) ? reuse : new Uri( parseAbsoluteUri( line ) ) );
            }
            catch ( final Exception e ) {
                throw new IllegalStateException( "not a URL (validate first): " + line, e );
            }
        }
        return out;
    }

    /** Keeps letters, gap/terminator characters ({@code - . * ?}); drops whitespace, digits (FASTA line numbers)
     *  and anything else, so a sequence pasted from a file or a web page comes out clean. */
    static String sanitizeMolSeq( final String s ) {
        return ( s == null ) ? "" : NOT_RESIDUE.matcher( s ).replaceAll( "" );
    }

    /** How many characters {@link #sanitizeMolSeq} would keep -- without building the cleaned string. */
    static int countResidues( final String s ) {
        if ( s == null ) {
            return 0;
        }
        int n = 0;
        for( int i = 0; i < s.length(); ++i ) {
            final char c = s.charAt( i );
            if ( ( ( c >= 'A' ) && ( c <= 'Z' ) ) || ( ( c >= 'a' ) && ( c <= 'z' ) ) || ( c == '-' ) || ( c == '.' )
                    || ( c == '*' ) || ( c == '?' ) ) {
                ++n;
            }
        }
        return n;
    }

    private static boolean isPrefixed( final String s ) {
        return ( s != null ) && ( s.trim().indexOf( ':' ) >= 1 );
    }

    private static boolean isDecimal( final String s ) {
        try {
            new BigDecimal( s.trim() );
            return true;
        }
        catch ( final NumberFormatException e ) {
            return false;
        }
    }

    private static Point firstPoint( final Distribution d ) {
        if ( ( d.getPoints() == null ) || d.getPoints().isEmpty() ) {
            return null;
        }
        return d.getPoints().get( 0 );
    }

    /** True when {@code text} is empty or a finite double (adds a problem and returns false otherwise). */
    private static boolean checkDouble( final List<Problem> ps, final String key, final String text,
                                        final String message ) {
        if ( text.trim().isEmpty() ) {
            return true;
        }
        try {
            final double d = Double.parseDouble( text.trim() );
            if ( Double.isNaN( d ) || Double.isInfinite( d ) ) {
                ps.add( new Problem( key, message ) );
                return false;
            }
            return true;
        }
        catch ( final NumberFormatException e ) {
            ps.add( new Problem( key, message ) );
            return false;
        }
    }

    private static boolean checkDecimal( final List<Problem> ps, final String key, final String text,
                                         final String what ) {
        if ( text.trim().isEmpty() || isDecimal( text ) ) {
            return true;
        }
        ps.add( new Problem( key, what + " must be a number" ) );
        return false;
    }

    private static void checkCount( final List<Problem> ps, final String key, final String text, final String what ) {
        if ( text.trim().isEmpty() ) {
            return;
        }
        try {
            if ( Integer.parseInt( text.trim() ) < 0 ) {
                ps.add( new Problem( key, what + " cannot be negative" ) );
            }
        }
        catch ( final NumberFormatException e ) {
            ps.add( new Problem( key, what + " must be a whole number" ) );
        }
    }

    /** Parses an ABSOLUTE URI (one with a scheme, like https://...); throws for anything else. */
    private static URI parseAbsoluteUri( final String s ) throws java.net.URISyntaxException {
        final URI u = new URI( s );
        if ( !u.isAbsolute() ) {
            throw new java.net.URISyntaxException( s, "missing scheme (e.g. https://)" );
        }
        return u;
    }

    private static void checkUris( final List<Problem> ps, final String key, final String text, final String what ) {
        for( final String line : lines( text ) ) {
            try {
                parseAbsoluteUri( line );
            }
            catch ( final Exception e ) {
                ps.add( new Problem( key, what + " \"" + line + "\" is not a valid URL" ) );
                return;
            }
        }
    }
}
