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

package org.forester.ws.seqdb;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;

/**
 * Enriches a phylogeny's nodes from sequence + taxonomy databases via clean, timeout/cancel-aware
 * clients ({@link SequenceFetcher} + {@link TaxonomyResolver}).
 *
 * <p>Per node it (1) fetches the sequence entry for the node's accession and writes the sequence +
 * organism fields + source-feature properties, then (2) resolves detailed taxonomy (rank, common name,
 * full lineage, NCBI tax-id) and refines the node's taxonomy. It polls a {@link CancelFlag} between
 * nodes and stops cleanly on the first transport error. The node-writers are pure (no network) so the
 * whole thing is unit-testable with fake clients. Network I/O happens here -- call off the EDT.
 */
public final class SequenceTaxonomyResolver {

    /** Outcome of a {@link #resolve} run. */
    public static final class Result {

        private final SortedSet<String> _not_found;
        private final int               _sequences_written;
        private final int               _taxonomies_written;
        private final boolean           _cancelled;
        private final String            _error;

        public Result( final SortedSet<String> not_found,
                       final int sequences_written,
                       final int taxonomies_written,
                       final boolean cancelled,
                       final String error ) {
            _not_found = not_found;
            _sequences_written = sequences_written;
            _taxonomies_written = taxonomies_written;
            _cancelled = cancelled;
            _error = error;
        }

        /** External nodes for which neither sequence nor taxonomy data could be obtained. */
        public SortedSet<String> getNotFound() {
            return _not_found;
        }

        public int getSequencesWritten() {
            return _sequences_written;
        }

        public int getTaxonomiesWritten() {
            return _taxonomies_written;
        }

        /** True if the user cancelled before the walk finished. */
        public boolean wasCancelled() {
            return _cancelled;
        }

        /** A transport-error message that aborted the walk, or {@code null} if none. */
        public String getError() {
            return _error;
        }
    }

    private final SequenceFetcher   _seq_fetcher;
    private final TaxonomyResolver  _tax_resolver;

    public SequenceTaxonomyResolver( final SequenceFetcher seq_fetcher, final TaxonomyResolver tax_resolver ) {
        _seq_fetcher = seq_fetcher;
        _tax_resolver = tax_resolver;
    }

    /**
     * Walks {@code phy}'s external nodes, enriching each from the databases. Stops promptly when
     * {@code cancel} is set, and aborts (returning an error) on the first transport failure. Mutates
     * the passed phylogeny in place.
     */
    public Result resolve( final Phylogeny phy, final CancelFlag cancel ) {
        final SortedSet<String> not_found = new TreeSet<String>();
        int seq_written = 0;
        int tax_written = 0;
        if ( ( phy == null ) || phy.isEmpty() ) {
            return new Result( not_found, 0, 0, false, null );
        }
        // Walk ALL nodes (postorder) -- the legacy tool enriched internal nodes too (e.g. accessions or
        // named clades on internal nodes); only external nodes are reported as "not found".
        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            if ( cancel.isCancelled() ) {
                return new Result( not_found, seq_written, tax_written, true, null );
            }
            final PhylogenyNode node = it.next();
            boolean wrote = false;
            try {
                String entry_tax_id = null;
                final Accession acc = SequenceAccessionTools.obtainAccessorFromDataFields( node );
                if ( ( acc != null ) && !ForesterUtil.isEmpty( acc.getValue() ) ) {
                    final SequenceEntry entry = _seq_fetcher.fetch( acc );
                    if ( ( entry != null ) && !entry.isEmpty() ) {
                        writeSequence( node, entry, acc.getSource() );
                        writeTaxonomyFromEntry( node, entry );
                        entry_tax_id = entry.getOrganismId();
                        ++seq_written;
                        wrote = true;
                    }
                }
                // Detail query: prefer the entry's authoritative NCBI tax-id (an exact efetch, no
                // esearch name-ambiguity that could clobber the organism with a wrong taxon); else the
                // node's scientific name / node name.
                final String query = !ForesterUtil.isEmpty( entry_tax_id ) ? entry_tax_id : taxonomyQuery( node );
                if ( !ForesterUtil.isEmpty( query ) ) {
                    final ResolvedTaxonomy rt = _tax_resolver.resolveTaxonomy( query );
                    if ( ( rt != null ) && !rt.isEmpty() ) {
                        writeTaxonomyDetail( node, rt );
                        ++tax_written;
                        wrote = true;
                    }
                }
            }
            catch ( final IOException e ) {
                // a transport failure will hit every remaining node -- abort the whole walk cleanly
                return new Result( not_found, seq_written, tax_written, false, friendly( e ) );
            }
            if ( !wrote && node.isExternal() ) {
                not_found.add( label( node ) );
            }
        }
        return new Result( not_found, seq_written, tax_written, false, null );
    }

    /** The most specific taxon query for a node: its taxonomy scientific name, else its node name. */
    static String taxonomyQuery( final PhylogenyNode node ) {
        if ( node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
            return node.getNodeData().getTaxonomy().getScientificName();
        }
        return ForesterUtil.isEmpty( node.getName() ) ? null : node.getName();
    }

    private static String label( final PhylogenyNode node ) {
        if ( !ForesterUtil.isEmpty( node.getName() ) ) {
            return node.getName();
        }
        if ( node.getNodeData().isHasTaxonomy()
                && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
            return node.getNodeData().getTaxonomy().getScientificName();
        }
        return "[" + node.getId() + "]";
    }

    private static String friendly( final IOException e ) {
        final String m = e.getLocalizedMessage();
        return ForesterUtil.isEmpty( m ) ? "could not connect to the database" : m;
    }

    // --- pure node writers (no network) -----------------------------------------------------------

    static void writeSequence( final PhylogenyNode node, final SequenceEntry entry, final String source ) {
        final boolean had_seq = node.getNodeData().isHasSequence();
        final Sequence seq = had_seq ? node.getNodeData().getSequence() : new Sequence();
        if ( !ForesterUtil.isEmpty( entry.getAccession() ) ) {
            seq.setAccession( new Accession( entry.getAccession(), source ) );
        }
        if ( !ForesterUtil.isEmpty( entry.getSequenceName() ) ) {
            seq.setName( entry.getSequenceName() );
        }
        if ( !ForesterUtil.isEmpty( entry.getGeneName() ) ) {
            seq.setGeneName( entry.getGeneName() );
        }
        // Only write the molecular sequence when the node has none -- never clobber a pre-existing
        // (possibly aligned/curated) sequence, and never reset its alignment flag on top of it.
        if ( !ForesterUtil.isEmpty( entry.getMolecularSequence() )
                && ForesterUtil.isEmpty( seq.getMolecularSequence() ) ) {
            seq.setMolecularSequence( entry.getMolecularSequence() );
            seq.setMolecularSequenceAligned( false );
            final String type = typeString( entry.getMoleculeType() );
            if ( type != null ) {
                try {
                    seq.setType( type );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    // leave the type unset if the vocabulary rejects it
                }
            }
        }
        for( final String go : entry.getGoIds() ) {
            seq.addAnnotation( new Annotation( go ) );
        }
        for( final String x : entry.getCrossReferences() ) {
            final int colon = x.indexOf( ':' );
            if ( colon > 0 ) {
                seq.addCrossReference( new Accession( x.substring( colon + 1 ), x.substring( 0, colon ) ) );
            }
        }
        if ( !had_seq ) {
            node.getNodeData().addSequence( seq );
        }
        writeProperties( node, entry );
    }

    /** Writes the entry's source-feature metadata (strain/host/country/...) as phyloXML node properties. */
    private static void writeProperties( final PhylogenyNode node, final SequenceEntry entry ) {
        if ( entry.getProperties().isEmpty() ) {
            return;
        }
        final boolean had_props = node.getNodeData().isHasProperties();
        final PropertiesList props = had_props ? node.getNodeData().getProperties() : new PropertiesList();
        for( final Map.Entry<String, String> e : entry.getProperties().entrySet() ) {
            // skip refs already present so a re-run does not accumulate duplicate properties
            if ( props.getPropertiesWithGivenRef( e.getKey() ).isEmpty() ) {
                props.addProperty( new Property( e.getKey(), e.getValue(), "", "xsd:string", Property.AppliesTo.NODE ) );
            }
        }
        if ( !had_props ) {
            node.getNodeData().setProperties( props );
        }
    }

    static void writeTaxonomyFromEntry( final PhylogenyNode node, final SequenceEntry entry ) {
        if ( ForesterUtil.isEmpty( entry.getOrganismName() ) && entry.getLineage().isEmpty()
                && ForesterUtil.isEmpty( entry.getOrganismId() ) ) {
            return;
        }
        final boolean had_tax = node.getNodeData().isHasTaxonomy();
        final Taxonomy tax = had_tax ? node.getNodeData().getTaxonomy() : new Taxonomy();
        if ( !ForesterUtil.isEmpty( entry.getOrganismName() ) ) {
            tax.setScientificName( entry.getOrganismName() );
        }
        if ( !entry.getLineage().isEmpty() ) {
            tax.setLineage( new ArrayList<String>( entry.getLineage() ) );
        }
        if ( !ForesterUtil.isEmpty( entry.getOrganismId() ) ) {
            tax.setIdentifier( new Identifier( entry.getOrganismId(), "ncbi" ) );
        }
        if ( !had_tax ) {
            node.getNodeData().addTaxonomy( tax );
        }
    }

    static void writeTaxonomyDetail( final PhylogenyNode node, final ResolvedTaxonomy rt ) {
        final boolean had_tax = node.getNodeData().isHasTaxonomy();
        final Taxonomy tax = had_tax ? node.getNodeData().getTaxonomy() : new Taxonomy();
        if ( !ForesterUtil.isEmpty( rt.getScientificName() ) ) {
            tax.setScientificName( rt.getScientificName() );
        }
        if ( !ForesterUtil.isEmpty( rt.getCommonName() ) ) {
            tax.setCommonName( rt.getCommonName() );
        }
        if ( !ForesterUtil.isEmpty( rt.getRank() ) ) {
            try {
                tax.setRank( rt.getRank().toLowerCase( Locale.ROOT ) );
            }
            catch ( final PhyloXmlDataFormatException e ) {
                // a rank outside the phyloXML vocabulary -- leave it unset
            }
        }
        if ( !rt.getLineage().isEmpty() ) {
            tax.setLineage( new ArrayList<String>( rt.getLineage() ) );
        }
        if ( !ForesterUtil.isEmpty( rt.getTaxId() ) ) {
            tax.setIdentifier( new Identifier( rt.getTaxId(), "ncbi" ) );
        }
        if ( !had_tax ) {
            node.getNodeData().addTaxonomy( tax );
        }
    }

    private static String typeString( final SequenceEntry.MoleculeType type ) {
        if ( type == null ) {
            return null;
        }
        switch ( type ) {
            case PROTEIN:
                return "protein";
            case DNA:
                return "dna";
            case RNA:
                return "rna";
            default:
                return null;
        }
    }
}
