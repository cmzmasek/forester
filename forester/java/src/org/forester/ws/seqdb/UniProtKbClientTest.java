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

/**
 * Unit tests for {@link UniProtKbClient#parseTsv}, driven by captured TSV fixtures (no network). The
 * parser reads columns positionally in the requested field order, so the tests pin that contract.
 */
public final class UniProtKbClientTest {

    private static final String HEADER = String.join( "\t", "Entry", "Entry Name", "Protein names",
                                                      "Gene Names (primary)", "Gene Names", "Sequence", "Organism",
                                                      "Organism (ID)", "Taxonomic lineage", "Gene Ontology IDs",
                                                      "RefSeq", "PDB" );
    private static final String SEQ    = "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKK";

    public static void main( final String[] args ) {
        System.out.println( "UniProtKbClient: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        // a full row: gene_primary present, two GO ids, one RefSeq xref, an EMPTY trailing PDB column, and a
        // lineage in UniProt's real "Name (rank)" format (the "(rank)" must be stripped per entry)
        final String lineage = "Eukaryota (domain), Metazoa (kingdom), Chordata (phylum), Mammalia (class), Primates (order), Hominidae (family), Homo (genus)";
        final String row = String.join( "\t", "P99999", "CYC_HUMAN", "Cytochrome c", "CYCS", "CYCS CYC", SEQ,
                                        "Homo sapiens (Human)", "9606", lineage, "GO:0005739; GO:0005829",
                                        "NP_061820.1;", "" );
        final SequenceEntry d = UniProtKbClient.parseTsv( HEADER + "\n" + row + "\n" );
        if ( d.isEmpty() ) {
            return fail( "a full row must parse to a non-empty entry" );
        }
        if ( !"P99999".equals( d.getAccession() ) || !"CYC_HUMAN".equals( d.getEntryName() )
                || !"Cytochrome c".equals( d.getSequenceName() ) ) {
            return fail( "accession/entry/protein wrong: " + d );
        }
        if ( !"CYCS".equals( d.getGeneName() ) ) { // gene_primary wins over gene_names
            return fail( "gene_primary must win; got " + d.getGeneName() );
        }
        if ( !SEQ.equals( d.getMolecularSequence() ) || ( d.getMoleculeType() != SequenceEntry.MoleculeType.PROTEIN )
                || !"Homo sapiens".equals( d.getOrganismName() ) || !"9606".equals( d.getOrganismId() ) ) {
            return fail( "sequence/type/organism wrong: " + d );
        }
        if ( !d.getLineage().contains( "Eukaryota" ) || !"Homo".equals( d.getLineage().get( d.getLineage().size() - 1 ) )
                || d.getLineage().toString().contains( "(" ) ) {
            return fail( "lineage must split on ',' and strip each \"(rank)\", ending in a clean Homo; got "
                    + d.getLineage() );
        }
        if ( ( d.getGoIds().size() != 2 ) || !d.getGoIds().contains( "GO:0005739" )
                || !d.getGoIds().contains( "GO:0005829" ) ) {
            return fail( "GO ids must split on ';'; got " + d.getGoIds() );
        }
        // refseq present (trailing ';' dropped), pdb empty -> only the refseq xref, "db:id" formatted
        if ( ( d.getCrossReferences().size() != 1 ) || !d.getCrossReferences().contains( "refseq:NP_061820.1" ) ) {
            return fail( "cross-references wrong (trailing-empty PDB column must be handled); got "
                    + d.getCrossReferences() );
        }

        // gene-name fallback: empty gene_primary -> first token of gene_names
        final String row2 = String.join( "\t", "Q12345", "X_HUMAN", "Some protein", "", "ABC DEF", "MKV",
                                          "Homo sapiens", "9606", "Eukaryota", "", "", "" );
        if ( !"ABC".equals( UniProtKbClient.parseTsv( HEADER + "\n" + row2 + "\n" ).getGeneName() ) ) {
            return fail( "gene-name must fall back to the first gene_names token when gene_primary is empty" );
        }

        // edge cases: empty input and header-only must yield EMPTY (not throw)
        if ( !UniProtKbClient.parseTsv( "" ).isEmpty() || !UniProtKbClient.parseTsv( HEADER + "\n" ).isEmpty() ) {
            return fail( "empty / header-only TSV must yield an empty entry" );
        }
        // the search query is field-qualified by shape: an entry name (has '_') -> id:, else accession:.
        // UniProt 400s on a non-accession value in the accession: field, so the OR-both approach fails.
        if ( !"id:RL7_HUMAN".equals( UniProtKbClient.searchQuery( "RL7_HUMAN" ) ) ) {
            return fail( "a SwissProt entry name must query the id: field; got " + UniProtKbClient.searchQuery( "RL7_HUMAN" ) );
        }
        if ( !"accession:P12345".equals( UniProtKbClient.searchQuery( "P12345" ) ) ) {
            return fail( "a UniProt accession must query the accession: field; got " + UniProtKbClient.searchQuery( "P12345" ) );
        }

        // the taxonomy-only projection (accession, id, organism_name, organism_id) -> Organism, no record kept
        final String org_header = String.join( "\t", "Entry", "Entry Name", "Organism", "Organism (ID)" );
        final Organism org = UniProtKbClient
                .parseOrganismTsv( org_header + "\n" + String.join( "\t", "P99999", "CYC_HUMAN", "Homo sapiens (Human)",
                                                                    "9606" ) + "\n" );
        if ( org.isEmpty() || !"9606".equals( org.getTaxId() ) || !"Homo sapiens".equals( org.getScientificName() ) ) {
            return fail( "parseOrganismTsv must read the tax-id and strip the organism's common name; got " + org );
        }
        if ( !UniProtKbClient.parseOrganismTsv( org_header + "\n" ).isEmpty()
                || !UniProtKbClient.parseOrganismTsv( "" ).isEmpty() ) {
            return fail( "header-only / empty organism TSV must yield an empty Organism" );
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  UniProtKbClient test failed: " + msg );
        return false;
    }

    private UniProtKbClientTest() {
    }
}
