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

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * The fields of one sequence-database entry needed to enrich a tree node, in a source-agnostic shape so
 * the UniProt and NCBI clients return the same type: accession + entry name, sequence (protein/product)
 * name, gene name, the molecular sequence + its type, organism scientific name + NCBI tax-id + lineage,
 * GO term ids, database cross-references, and free-form node properties keyed by phyloXML property ref
 * (e.g. the GenBank source-feature qualifiers strain/host/country/isolate/...). Immutable. Consumed by
 * {@code SequenceTaxonomyResolver}.
 *
 * <p>Any field may be {@code null}/empty; {@link #isEmpty()} is true when nothing usable was returned.
 */
public final class SequenceEntry {

    /** Molecule type of {@link #getMolecularSequence()} (maps to phyloXML sequence types). */
    public enum MoleculeType {
        PROTEIN, DNA, RNA
    }

    /** An empty entry (accession not found / nothing usable). */
    public static final SequenceEntry EMPTY = new SequenceEntry( null, null, null, null, null, null, null, null, null,
                                                                 null, null, null );
    private final String              _accession;
    private final String              _entry_name;
    private final String              _sequence_name;
    private final String              _gene_name;
    private final String              _molecular_sequence;
    private final MoleculeType        _molecule_type;
    private final String              _organism_name;
    private final String              _organism_id;
    private final List<String>        _lineage;
    private final List<String>        _go_ids;
    private final List<String>        _cross_references; // each formatted "db:id" (e.g. "refseq:NP_000001")
    private final Map<String, String> _properties;       // phyloXML property ref -> value (e.g. "source:country" -> "USA")

    public SequenceEntry( final String accession,
                          final String entry_name,
                          final String sequence_name,
                          final String gene_name,
                          final String molecular_sequence,
                          final MoleculeType molecule_type,
                          final String organism_name,
                          final String organism_id,
                          final List<String> lineage,
                          final List<String> go_ids,
                          final List<String> cross_references,
                          final Map<String, String> properties ) {
        _accession = accession;
        _entry_name = entry_name;
        _sequence_name = sequence_name;
        _gene_name = gene_name;
        _molecular_sequence = molecular_sequence;
        _molecule_type = molecule_type;
        _organism_name = organism_name;
        _organism_id = organism_id;
        _lineage = copy( lineage );
        _go_ids = copy( go_ids );
        _cross_references = copy( cross_references );
        _properties = ( properties == null ) ? Collections.<String, String> emptyMap()
                : Collections.unmodifiableMap( new LinkedHashMap<String, String>( properties ) );
    }

    private static List<String> copy( final List<String> in ) {
        return ( in == null ) ? Collections.<String> emptyList()
                : Collections.unmodifiableList( new ArrayList<String>( in ) );
    }

    public String getAccession() {
        return _accession;
    }

    public String getEntryName() {
        return _entry_name;
    }

    /** The protein/product name (UniProt protein_name / GenBank definition), or {@code null}. */
    public String getSequenceName() {
        return _sequence_name;
    }

    public String getGeneName() {
        return _gene_name;
    }

    /** The molecular sequence as a plain string, or {@code null}. */
    public String getMolecularSequence() {
        return _molecular_sequence;
    }

    public MoleculeType getMoleculeType() {
        return _molecule_type;
    }

    public String getOrganismName() {
        return _organism_name;
    }

    /** The organism's NCBI tax-id (numeric string), or {@code null}. */
    public String getOrganismId() {
        return _organism_id;
    }

    /** Unmodifiable lineage of scientific names, root&rarr;organism (never null). */
    public List<String> getLineage() {
        return _lineage;
    }

    /** Unmodifiable list of GO ids (e.g. "GO:0005737"), never null. */
    public List<String> getGoIds() {
        return _go_ids;
    }

    /** Unmodifiable list of cross-references, each "db:id" (e.g. "refseq:NP_000001"), never null. */
    public List<String> getCrossReferences() {
        return _cross_references;
    }

    /** Unmodifiable node properties keyed by phyloXML ref (e.g. "source:country" -&gt; "USA"), never null. */
    public Map<String, String> getProperties() {
        return _properties;
    }

    public boolean isEmpty() {
        return isBlank( _accession ) && isBlank( _molecular_sequence ) && isBlank( _sequence_name );
    }

    private static boolean isBlank( final String s ) {
        return ( s == null ) || s.trim().isEmpty();
    }

    @Override
    public String toString() {
        return "SequenceEntry[" + _accession + " (" + _entry_name + "), name=" + _sequence_name + ", gene=" + _gene_name
                + ", organism=" + _organism_name + " (" + _organism_id + "), type=" + _molecule_type + ", seqLen="
                + ( ( _molecular_sequence == null ) ? 0 : _molecular_sequence.length() ) + "]";
    }
}
