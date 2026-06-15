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
 * Unit tests for {@link NcbiSequenceClient#parseGbSeqXml}, driven by a captured GBSeq XML fixture (no
 * network) -- pins the field extraction including the gene name / tax-id / cross-references that live in
 * the feature-table qualifiers.
 */
public final class NcbiSequenceClientTest {

    private static final String GBSEQ = """
            <?xml version="1.0"?>
            <!DOCTYPE GBSet PUBLIC "-//NCBI//NCBI GBSeq/EN" "https://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd">
            <GBSet>
              <GBSeq>
                <GBSeq_locus>NP_061820</GBSeq_locus>
                <GBSeq_length>105</GBSeq_length>
                <GBSeq_moltype>AA</GBSeq_moltype>
                <GBSeq_definition>cytochrome c [Homo sapiens]</GBSeq_definition>
                <GBSeq_primary-accession>NP_061820</GBSeq_primary-accession>
                <GBSeq_organism>Homo sapiens</GBSeq_organism>
                <GBSeq_taxonomy>Eukaryota; Metazoa; Chordata; Mammalia; Primates; Hominidae; Homo</GBSeq_taxonomy>
                <GBSeq_feature-table>
                  <GBFeature>
                    <GBFeature_key>source</GBFeature_key>
                    <GBFeature_quals>
                      <GBQualifier><GBQualifier_name>organism</GBQualifier_name><GBQualifier_value>Homo sapiens</GBQualifier_value></GBQualifier>
                      <GBQualifier><GBQualifier_name>strain</GBQualifier_name><GBQualifier_value>K-12</GBQualifier_value></GBQualifier>
                      <GBQualifier><GBQualifier_name>host</GBQualifier_name><GBQualifier_value>Aedes aegypti</GBQualifier_value></GBQualifier>
                      <GBQualifier><GBQualifier_name>geo_loc_name</GBQualifier_name><GBQualifier_value>USA: California</GBQualifier_value></GBQualifier>
                      <GBQualifier><GBQualifier_name>isolate</GBQualifier_name><GBQualifier_value>iso-42</GBQualifier_value></GBQualifier>
                      <GBQualifier><GBQualifier_name>collection_date</GBQualifier_name><GBQualifier_value>2021-03-15</GBQualifier_value></GBQualifier>
                      <GBQualifier><GBQualifier_name>db_xref</GBQualifier_name><GBQualifier_value>taxon:9606</GBQualifier_value></GBQualifier>
                    </GBFeature_quals>
                  </GBFeature>
                  <GBFeature>
                    <GBFeature_key>CDS</GBFeature_key>
                    <GBFeature_quals>
                      <GBQualifier><GBQualifier_name>gene</GBQualifier_name><GBQualifier_value>CYCS</GBQualifier_value></GBQualifier>
                      <GBQualifier><GBQualifier_name>db_xref</GBQualifier_name><GBQualifier_value>taxon:9606</GBQualifier_value></GBQualifier>
                      <GBQualifier><GBQualifier_name>db_xref</GBQualifier_name><GBQualifier_value>GeneID:54205</GBQualifier_value></GBQualifier>
                    </GBFeature_quals>
                  </GBFeature>
                </GBSeq_feature-table>
                <GBSeq_sequence>mgdvekgkkifimkcsqchtvekggkhktgpnlhglfgrktgqapgysytaank</GBSeq_sequence>
              </GBSeq>
            </GBSet>""";

    public static void main( final String[] args ) {
        System.out.println( "NcbiSequenceClient: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        final SequenceEntry e = NcbiSequenceClient.parseGbSeqXml( GBSEQ );
        if ( e.isEmpty() ) {
            return fail( "GBSeq entry did not parse" );
        }
        if ( !"NP_061820".equals( e.getAccession() ) || !"NP_061820".equals( e.getEntryName() )
                || !"cytochrome c [Homo sapiens]".equals( e.getSequenceName() ) ) {
            return fail( "accession/locus/definition wrong: " + e );
        }
        if ( e.getMoleculeType() != SequenceEntry.MoleculeType.PROTEIN ) {
            return fail( "moltype AA must map to PROTEIN; got " + e.getMoleculeType() );
        }
        if ( !"Homo sapiens".equals( e.getOrganismName() ) || !"9606".equals( e.getOrganismId() ) ) {
            return fail( "organism/tax-id (from db_xref taxon:) wrong: " + e );
        }
        if ( !"CYCS".equals( e.getGeneName() ) ) {
            return fail( "gene name (from qualifiers) wrong: " + e.getGeneName() );
        }
        if ( !( e.getMolecularSequence() != null ) || !e.getMolecularSequence().startsWith( "MGDVEK" ) ) {
            return fail( "molecular sequence must be upper-cased; got " + e.getMolecularSequence() );
        }
        if ( !e.getLineage().contains( "Eukaryota" )
                || !"Homo sapiens".equals( e.getLineage().get( e.getLineage().size() - 1 ) ) ) {
            return fail( "lineage must split on ';' and end with the organism; got " + e.getLineage() );
        }
        // the GeneID db_xref becomes a cross-reference; the taxon db_xref must NOT (it is the tax-id)
        if ( !e.getCrossReferences().contains( "geneid:54205" ) ) {
            return fail( "GeneID db_xref must become a cross-reference; got " + e.getCrossReferences() );
        }
        // source-feature qualifiers -> phyloXML property refs (the viral/strain metadata path);
        // note geo_loc_name maps to "source:country"
        final java.util.Map<String, String> props = e.getProperties();
        if ( !"K-12".equals( props.get( "taxonomy:strain" ) )
                || !"USA: California".equals( props.get( "source:country" ) )
                || !"Aedes aegypti".equals( props.get( "source:host" ) )
                || !"iso-42".equals( props.get( "source:isolate" ) )
                || !"2021-03-15".equals( props.get( "source:collection_date" ) ) ) {
            return fail( "source-feature properties not extracted: " + props );
        }
        for( final String x : e.getCrossReferences() ) {
            if ( x.startsWith( "taxon:" ) ) {
                return fail( "the taxon db_xref must not be a cross-reference" );
            }
        }
        // edge cases
        if ( !NcbiSequenceClient.parseGbSeqXml( "" ).isEmpty()
                || !NcbiSequenceClient.parseGbSeqXml( "<GBSet></GBSet>" ).isEmpty()
                || !NcbiSequenceClient.parseGbSeqXml( "<broken" ).isEmpty() ) {
            return fail( "empty/no-GBSeq/malformed must yield an empty entry, not throw" );
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  NcbiSequenceClient test failed: " + msg );
        return false;
    }

    private NcbiSequenceClientTest() {
    }
}
