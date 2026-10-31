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

import java.util.Objects;

import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Tests for {@link NodeDataLinks}: every rule pinned by the LITERAL address it must give (each address was opened by
 * hand once, 2026-10-30), and every refusal by a null -- an unknown source or provider, a value that is not an
 * accession, an address that is not http(s).
 */
public final class NodeDataLinksTest {

    private static final String NCBI_TAX_ID   = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=";
    private static final String NCBI_TAX_NAME = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name=";
    private static final String UNIPROT_TAX   = "https://www.uniprot.org/taxonomy/";
    private static final String UNIPROT_KB    = "https://www.uniprot.org/uniprot/";
    private static final String NCBI_PROTEIN  = "https://www.ncbi.nlm.nih.gov/protein/";
    private static final String NCBI_NUCCORE  = "https://www.ncbi.nlm.nih.gov/nuccore/";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeDataLinks: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return taxonomyIds() & names() & menuChain() & ncbiAccessions() & otherSources() & withoutSource()
                    & placeholderSources() & doi() & web() & malformedUnicode();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean taxonomyIds() {
        boolean ok = eq( "ncbi id", NCBI_TAX_ID + "9606", NodeDataLinks.taxonomyIdUrl( "9606", "ncbi" ) );
        ok &= eq( "provider case", NCBI_TAX_ID + "9606", NodeDataLinks.taxonomyIdUrl( "9606", "NCBI" ) );
        ok &= eq( "ncbi_taxonomy", NCBI_TAX_ID + "9606", NodeDataLinks.taxonomyIdUrl( "9606", "ncbi_taxonomy" ) );
        ok &= eq( "no provider reads as ncbi (Identifier.normalizedProvider)", NCBI_TAX_ID + "9606",
                  NodeDataLinks.taxonomyIdUrl( "9606", "" ) );
        ok &= eq( "null provider reads as ncbi", NCBI_TAX_ID + "9606", NodeDataLinks.taxonomyIdUrl( "9606", null ) );
        ok &= eq( "id trimmed", NCBI_TAX_ID + "9606", NodeDataLinks.taxonomyIdUrl( " 9606 ", "ncbi" ) );
        ok &= eq( "uniprot id: the provider's own database", UNIPROT_TAX + "9606",
                  NodeDataLinks.taxonomyIdUrl( "9606", "uniprot" ) );
        ok &= eq( "uniprot.taxonomy", UNIPROT_TAX + "9606",
                  NodeDataLinks.taxonomyIdUrl( "9606", "UniProt.Taxonomy" ) );
        ok &= eq( "an id that is an address opens it", "https://example.org/tax/1",
                  NodeDataLinks.taxonomyIdUrl( "https://example.org/tax/1", "anything" ) );
        ok &= eq( "unknown provider: no link", null, NodeDataLinks.taxonomyIdUrl( "9606", "gtdb" ) );
        ok &= eq( "a non-numeric ncbi id: no link", null, NodeDataLinks.taxonomyIdUrl( "txid9606", "ncbi" ) );
        ok &= eq( "empty id", null, NodeDataLinks.taxonomyIdUrl( "", "ncbi" ) );
        ok &= eq( "null id", null, NodeDataLinks.taxonomyIdUrl( null, "ncbi" ) );
        return ok;
    }

    private static boolean names() {
        boolean ok = eq( "name, no id: NCBI by name", NCBI_TAX_NAME + "Homo+sapiens",
                         NodeDataLinks.scientificNameUrl( "Homo sapiens", null ) );
        ok &= eq( "name with an id: the id's page", UNIPROT_TAX + "9606",
                  NodeDataLinks.scientificNameUrl( "Homo sapiens", UNIPROT_TAX + "9606" ) );
        ok &= eq( "name query-encoded", NCBI_TAX_NAME + "Cercopithecus+l%27hoesti",
                  NodeDataLinks.scientificNameUrl( " Cercopithecus l'hoesti ", null ) );
        ok &= eq( "no name", null, NodeDataLinks.scientificNameUrl( " ", NCBI_TAX_ID + "1" ) );
        ok &= eq( "code: UniProt mnemonic", "https://www.uniprot.org/taxonomy?query=mnemonic:HUMAN",
                  NodeDataLinks.taxonomyCodeUrl( "HUMAN" ) );
        ok &= eq( "viral mnemonic", "https://www.uniprot.org/taxonomy?query=mnemonic:9HIV1",
                  NodeDataLinks.taxonomyCodeUrl( "9HIV1" ) );
        ok &= eq( "lower case is not a mnemonic", null, NodeDataLinks.taxonomyCodeUrl( "human" ) );
        ok &= eq( "too short", null, NodeDataLinks.taxonomyCodeUrl( "HU" ) );
        ok &= eq( "too long", null, NodeDataLinks.taxonomyCodeUrl( "HUMANS" ) );
        return ok;
    }

    /** The node menu's "Open Taxonomy DB": id, then scientific name, then code -- a page the window links; never
     *  the common name, which the window does not link. */
    private static boolean menuChain() throws Exception {
        boolean ok = eq( "id first", NCBI_TAX_ID + "9606",
                         NodeDataLinks.taxonomyUrl( tax( "9606", "ncbi", "Homo sapiens", "HUMAN", "human" ) ) );
        ok &= eq( "uniprot id", UNIPROT_TAX + "9606",
                  NodeDataLinks.taxonomyUrl( tax( "9606", "uniprot", "Homo sapiens", "HUMAN", "human" ) ) );
        ok &= eq( "unknown provider: the name", NCBI_TAX_NAME + "Homo+sapiens",
                  NodeDataLinks.taxonomyUrl( tax( "g__Homo", "gtdb", "Homo sapiens", "HUMAN", "human" ) ) );
        ok &= eq( "no id: the name before the code", NCBI_TAX_NAME + "Homo+sapiens",
                  NodeDataLinks.taxonomyUrl( tax( null, null, "Homo sapiens", "HUMAN", "human" ) ) );
        ok &= eq( "code only", "https://www.uniprot.org/taxonomy?query=mnemonic:HUMAN",
                  NodeDataLinks.taxonomyUrl( tax( null, null, null, "HUMAN", "human" ) ) );
        ok &= eq( "common name only: no link (the window links none)", null,
                  NodeDataLinks.taxonomyUrl( tax( null, null, null, null, "human" ) ) );
        ok &= eq( "an id that names nothing, a common name: still no link", null,
                  NodeDataLinks.taxonomyUrl( tax( "180092", "itis", null, null, "human" ) ) );
        ok &= eq( "an ncbi id alone opens it", NCBI_TAX_ID + "9606",
                  NodeDataLinks.taxonomyUrl( tax( "9606", "ncbi", null, null, null ) ) );
        ok &= eq( "a uniprot id alone opens it", UNIPROT_TAX + "9606",
                  NodeDataLinks.taxonomyUrl( tax( "9606", "uniprot", null, null, null ) ) );
        ok &= eq( "an id that is an address, alone", "https://example.org/tax/1",
                  NodeDataLinks.taxonomyUrl( tax( "https://example.org/tax/1", "x", null, null, null ) ) );
        ok &= eq( "an id and a code, no name: the id", NCBI_TAX_ID + "9606",
                  NodeDataLinks.taxonomyUrl( tax( "9606", "ncbi", null, "HUMAN", null ) ) );
        ok &= eq( "an id and a common name: the id", NCBI_TAX_ID + "9606",
                  NodeDataLinks.taxonomyUrl( tax( "9606", "ncbi", null, null, "human" ) ) );
        ok &= eq( "an id nobody can open, nothing else: no link", null,
                  NodeDataLinks.taxonomyUrl( tax( "g__Homo", "gtdb", null, null, null ) ) );
        ok &= eq( "no taxonomy", null, NodeDataLinks.taxonomyUrl( null ) );
        return ok;
    }

    private static boolean ncbiAccessions() {
        boolean ok = eq( "RefSeq protein", NCBI_PROTEIN + "YP_013377151.1",
                         NodeDataLinks.accessionUrl( "YP_013377151.1", "ncbi", "protein" ) );
        // the FORMAT decides, not the declared type: a protein sequence carrying its genome's accession
        ok &= eq( "a genome accession on a 'protein' sequence opens Nucleotide", NCBI_NUCCORE + "KJ473821",
                  NodeDataLinks.accessionUrl( "KJ473821", "ncbi", "protein" ) );
        ok &= eq( "a RefSeq genome on a 'protein' sequence", NCBI_NUCCORE + "NC_009021",
                  NodeDataLinks.accessionUrl( "NC_009021", "ncbi", "protein" ) );
        ok &= eq( "a WGS contig on a 'protein' sequence", NCBI_NUCCORE + "BTIZ01000004.1",
                  NodeDataLinks.accessionUrl( "BTIZ01000004.1", "ncbi", "protein" ) );
        ok &= eq( "GenBank protein, no type", NCBI_PROTEIN + "ASU90448",
                  NodeDataLinks.accessionUrl( "ASU90448", "ncbi", null ) );
        ok &= eq( "GenBank protein, seven digits", NCBI_PROTEIN + "WAA1234567.1",
                  NodeDataLinks.accessionUrl( "WAA1234567.1", "GenBank", "" ) );
        ok &= eq( "GenBank nucleotide", NCBI_NUCCORE + "MN908947.3",
                  NodeDataLinks.accessionUrl( "MN908947.3", "GenBank", null ) );
        ok &= eq( "RefSeq mRNA", NCBI_NUCCORE + "NM_007294", NodeDataLinks.accessionUrl( "NM_007294", "RefSeq", "rna" ) );
        ok &= eq( "RefSeq WGS", NCBI_NUCCORE + "NZ_SYN09000001.1",
                  NodeDataLinks.accessionUrl( "NZ_SYN09000001.1", "RefSeq", "dna" ) );
        ok &= eq( "RefSeq WP_", NCBI_PROTEIN + "WP_100002001.1",
                  NodeDataLinks.accessionUrl( "WP_100002001.1", "refseq", null ) );
        ok &= eq( "an accession with a segment label after a separator", NCBI_PROTEIN + "QEA08541.1",
                  NodeDataLinks.accessionUrl( "QEA08541.1:NS5", "ncbi", "protein" ) );
        ok &= eq( "a Swiss-Prot accession under ncbi: NCBI Protein, version kept", NCBI_PROTEIN + "Q5UQR0.1",
                  NodeDataLinks.accessionUrl( "Q5UQR0.1", "ncbi", null ) );
        ok &= eq( "[OPQ]+5 digits, declared dna: Nucleotide", NCBI_NUCCORE + "P38398",
                  NodeDataLinks.accessionUrl( "P38398", "ncbi", "dna" ) );
        ok &= eq( "[OPQ]+5 digits, no type: Protein", NCBI_PROTEIN + "P38398",
                  NodeDataLinks.accessionUrl( "P38398", "ncbi", null ) );
        ok &= eq( "an unplaced shape: the declared type", NCBI_PROTEIN + "XYZ12345678",
                  NodeDataLinks.accessionUrl( "XYZ12345678", "ncbi", "protein" ) );
        ok &= eq( "an unplaced shape, no type: no link", null, NodeDataLinks.accessionUrl( "XYZ12345678", "ncbi", null ) );
        ok &= eq( "an unplaced shape with a suffix: no link", null,
                  NodeDataLinks.accessionUrl( "XYZ12345678:x", "ncbi", "protein" ) );
        ok &= eq( "not an accession under ncbi: no link", null,
                  NodeDataLinks.accessionUrl( "fig|1000.9.peg.1", "ncbi", "protein" ) );
        ok &= eq( "a value glued to more characters: no link", null,
                  NodeDataLinks.accessionUrl( "KJ473821-x", "ncbi", null ) );
        ok &= eq( "gi", "https://www.ncbi.nlm.nih.gov/protein/gi:123456",
                  NodeDataLinks.accessionUrl( "123456", "gi", null ) );
        ok &= eq( "empty value", null, NodeDataLinks.accessionUrl( "  ", "ncbi", "protein" ) );
        ok &= eq( "gi must be a number", null, NodeDataLinks.accessionUrl( "hello", "gi", null ) );
        ok &= eq( "a trailing dot is not part of it", NCBI_NUCCORE + "KJ473821",
                  NodeDataLinks.accessionUrl( "KJ473821.", "ncbi", null ) );
        // RefSeq is NCBI's closed prefix list: a Pathoplexus id under an ncbi source is no NCBI page (the 32 such
        // records in the Orthoebolavirus reference trees opened a protein page that does not exist)
        ok &= eq( "RefSeq genome", NCBI_NUCCORE + "NC_045512.2", NodeDataLinks.accessionUrl( "NC_045512.2", "ncbi", "protein" ) );
        ok &= eq( "RefSeq WGS", NCBI_NUCCORE + "NZ_CP009072.1", NodeDataLinks.accessionUrl( "NZ_CP009072.1", "ncbi", null ) );
        ok &= eq( "RefSeq WP_", NCBI_PROTEIN + "WP_000117906.1", NodeDataLinks.accessionUrl( "WP_000117906.1", "ncbi", "dna" ) );
        ok &= eq( "not a RefSeq prefix: no link", null, NodeDataLinks.accessionUrl( "PP_00711R7.2", "ncbi", "dna" ) );
        // [OPQ] and five digits is a UniProt AND a GenBank nucleotide shape: the declared type decides
        ok &= eq( "either shape, rna", NCBI_NUCCORE + "P12345", NodeDataLinks.accessionUrl( "P12345", "ncbi", "rna" ) );
        ok &= eq( "either shape, dna", NCBI_NUCCORE + "P12345", NodeDataLinks.accessionUrl( "P12345", "ncbi", "dna" ) );
        ok &= eq( "either shape, protein", NCBI_PROTEIN + "P12345", NodeDataLinks.accessionUrl( "P12345", "ncbi", "protein" ) );
        ok &= eq( "either shape, no type: protein", NCBI_PROTEIN + "P12345", NodeDataLinks.accessionUrl( "P12345", "ncbi", null ) );
        for( final String src : new String[] { "genbank", "refseq", "insdc", "ddbj", "NCBI" } ) {
            ok &= eq( "an NCBI source: " + src, NCBI_NUCCORE + "KJ473821", NodeDataLinks.accessionUrl( "KJ473821", src, null ) );
        }
        return ok;
    }

    private static boolean otherSources() {
        boolean ok = eq( "UniProt", UNIPROT_KB + "P38398", NodeDataLinks.accessionUrl( "P38398", "UniProt", "protein" ) );
        ok &= eq( "UniProt drops the version", UNIPROT_KB + "Q5UQR0", NodeDataLinks.accessionUrl( "Q5UQR0.1", "uniprot", null ) );
        ok &= eq( "UniProtKB", UNIPROT_KB + "A0A096XNM6", NodeDataLinks.accessionUrl( "A0A096XNM6", "UniProtKB", null ) );
        ok &= eq( "an entry name", UNIPROT_KB + "BRCA1_HUMAN",
                  NodeDataLinks.accessionUrl( "BRCA1_HUMAN", "Swiss-Prot", null ) );
        for( final String src : new String[] { "sp", "tr", "trembl", "swissprot", "UniProtKB/Swiss-Prot" } ) {
            ok &= eq( "a UniProt source: " + src, UNIPROT_KB + "P38398", NodeDataLinks.accessionUrl( "P38398", src, null ) );
        }
        // a FASTA-style composite links neither its database tag nor anything else
        ok &= eq( "uniprot: a composite is not an accession", null,
                  NodeDataLinks.accessionUrl( "sp|P38398|BRCA1_HUMAN", "uniprot", null ) );
        ok &= eq( "trembl: a composite is not an accession", null,
                  NodeDataLinks.accessionUrl( "tr|A0A096XNM6|A0A096XNM6_9VIRU", "trembl", null ) );
        ok &= eq( "uniprot: a word is not one", null, NodeDataLinks.accessionUrl( "hello", "uniprot", null ) );
        ok &= eq( "ENA: a composite is not one", null,
                  NodeDataLinks.accessionUrl( "ENA|MN908947|MN908947.3", "ENA", null ) );
        ok &= eq( "embl: a word is not one", null, NodeDataLinks.accessionUrl( "hello", "embl", null ) );
        ok &= eq( "Pathoplexus: a word is not one", null, NodeDataLinks.accessionUrl( "hello", "pathoplexus", null ) );
        ok &= eq( "NCBI Gene", "https://www.ncbi.nlm.nih.gov/gene/1489668",
                  NodeDataLinks.accessionUrl( "1489668", "geneid", null ) );
        ok &= eq( "a gene id must be a number", null, NodeDataLinks.accessionUrl( "BRCA1", "geneid", null ) );
        ok &= eq( "InterPro", "https://www.ebi.ac.uk/interpro/entry/InterPro/IPR002552/",
                  NodeDataLinks.accessionUrl( "IPR002552", "interpro", null ) );
        ok &= eq( "not an InterPro accession", null, NodeDataLinks.accessionUrl( "PF00931", "interpro", null ) );
        ok &= eq( "Pfam (the tree's own domain link)", "https://www.ebi.ac.uk/interpro/entry/pfam/PF00931/",
                  NodeDataLinks.accessionUrl( "PF00931.25", "Pfam", null ) );
        ok &= eq( "an id that merely CONTAINS a Pfam accession is not one", null,
                  NodeDataLinks.accessionUrl( "SM00109_PF00018", "pfam", null ) );
        ok &= eq( "Ensembl drops the version", "https://www.ensembl.org/id/ENSG00000139618",
                  NodeDataLinks.accessionUrl( "ENSG00000139618.17", "Ensembl", null ) );
        ok &= eq( "ENA", "https://www.ebi.ac.uk/ena/browser/view/MN908947",
                  NodeDataLinks.accessionUrl( "MN908947", "ENA", null ) );
        ok &= eq( "EMBL", "https://www.ebi.ac.uk/ena/browser/view/X04370",
                  NodeDataLinks.accessionUrl( "X04370", "embl", null ) );
        ok &= eq( "Pathoplexus", "https://pathoplexus.org/seq/PP_006XHL9.4",
                  NodeDataLinks.accessionUrl( "PP_006XHL9.4", "pathoplexus", null ) );
        ok &= eq( "an unknown source is not guessed at", null,
                  NodeDataLinks.accessionUrl( "fig|1000.9.peg.1", "BV-BRC", null ) );
        ok &= eq( "an unknown source with a GenBank-shaped value: still no guess", null,
                  NodeDataLinks.accessionUrl( "KJ473821", "BV-BRC", null ) );
        return ok;
    }

    private static boolean withoutSource() {
        boolean ok = eq( "no source, UniProt accession", UNIPROT_KB + "P38398",
                         NodeDataLinks.accessionUrl( "P38398", "", null ) );
        ok &= eq( "no source, RefSeq", NCBI_NUCCORE + "NM_007294", NodeDataLinks.accessionUrl( "NM_007294", null, null ) );
        ok &= eq( "no source, GenBank protein", NCBI_PROTEIN + "AAY88866",
                  NodeDataLinks.accessionUrl( "AAY88866", " ", null ) );
        ok &= eq( "no source, GenBank nucleotide", NCBI_NUCCORE + "JX519213.1",
                  NodeDataLinks.accessionUrl( "JX519213.1", "", null ) );
        ok &= eq( "no source: an entry name is not an accession", null,
                  NodeDataLinks.accessionUrl( "BRCA1_HUMAN", "", null ) );
        ok &= eq( "no source: a composite is not one", null, NodeDataLinks.accessionUrl( "sp|P38398|BRCA1", "", null ) );
        ok &= eq( "no source: a word", null, NodeDataLinks.accessionUrl( "hello", "", "protein" ) );
        ok &= eq( "no source: UniProt drops the version", UNIPROT_KB + "Q5UQR0", NodeDataLinks.accessionUrl( "Q5UQR0.1", "", null ) );
        ok &= eq( "no source: a Pathoplexus id is not RefSeq", null, NodeDataLinks.accessionUrl( "PP_006XHL9.4", "", null ) );
        ok &= eq( "no source: AB_ is not RefSeq", null, NodeDataLinks.accessionUrl( "AB_123456", "", null ) );
        return ok;
    }

    /** forester writes "unknown" (phyloXML) and "?" (NHX AC=) for an accession without a source: read as none. */
    private static boolean placeholderSources() {
        boolean ok = eq( "unknown: as no source", UNIPROT_KB + "P38398", NodeDataLinks.accessionUrl( "P38398", "unknown", null ) );
        ok &= eq( "Unknown, any case", NCBI_NUCCORE + "KJ473821", NodeDataLinks.accessionUrl( "KJ473821", " Unknown ", null ) );
        ok &= eq( "?: as no source", UNIPROT_KB + "P38398", NodeDataLinks.accessionUrl( "P38398", "?", null ) );
        ok &= eq( "unknown: still only an accession", null, NodeDataLinks.accessionUrl( "hello", "unknown", null ) );
        return ok;
    }

    private static boolean doi() {
        boolean ok = eq( "bare DOI", "https://doi.org/10.1126/science.7545954",
                         NodeDataLinks.doiUrl( "10.1126/science.7545954" ) );
        ok &= eq( "doi: prefix", "https://doi.org/10.1093/bioinformatics/btq243",
                  NodeDataLinks.doiUrl( "doi:10.1093/bioinformatics/btq243" ) );
        ok &= eq( "already an address", "https://doi.org/10.1093/bioinformatics/btq243",
                  NodeDataLinks.doiUrl( "http://dx.doi.org/10.1093/bioinformatics/btq243" ) );
        ok &= eq( "path-escaped, not form-encoded",
                  "https://doi.org/10.1002/(SICI)1097-4636(199706)35:4%3C493::AID-JBM9%3E3.0.CO;2-J",
                  NodeDataLinks.doiUrl( "10.1002/(SICI)1097-4636(199706)35:4<493::AID-JBM9>3.0.CO;2-J" ) );
        ok &= eq( "a doi.org address", "https://doi.org/10.1093/bioinformatics/btq243",
                  NodeDataLinks.doiUrl( "https://doi.org/10.1093/bioinformatics/btq243" ) );
        ok &= eq( "a doi.org address, any case", "https://doi.org/10.1093/bioinformatics/btq243",
                  NodeDataLinks.doiUrl( "HTTPS://DOI.ORG/10.1093/bioinformatics/btq243" ) );
        // as an address bar copies it: its escapes are kept (re-quoted, doi.org answers "DOI not found")
        ok &= eq( "an escaped address stays as written",
                  "https://doi.org/10.1002/(SICI)1097-0258(19980430)17:8%3C857::AID-SIM777%3E3.0.CO;2-E",
                  NodeDataLinks.doiUrl( "https://doi.org/10.1002/(SICI)1097-0258(19980430)17:8%3C857::AID-SIM777%3E3.0.CO;2-E" ) );
        ok &= eq( "a '#' or '?' in a DOI is part of its name", "https://doi.org/10.1000/a%23b%3Fc",
                  NodeDataLinks.doiUrl( "10.1000/a#b?c" ) );
        ok &= eq( "a citation is not a DOI", null, NodeDataLinks.doiUrl( "Miki et al. 1994" ) );
        ok &= eq( "too short a registrant", null, NodeDataLinks.doiUrl( "10.1/x" ) );
        ok &= eq( "empty", null, NodeDataLinks.doiUrl( "" ) );
        return ok;
    }

    private static boolean web() {
        boolean ok = eq( "an address", "https://example.org/a?b=c", NodeDataLinks.webUrl( "https://example.org/a?b=c" ) );
        ok &= eq( "already escaped stays as written", "https://example.org/a%20b",
                  NodeDataLinks.webUrl( "https://example.org/a%20b" ) );
        ok &= eq( "a '|' is escaped", "https://www.bv-brc.org/view/Feature/fig%7C1000.9.peg.1",
                  NodeDataLinks.webUrl( "https://www.bv-brc.org/view/Feature/fig|1000.9.peg.1" ) );
        ok &= eq( "non-ASCII is escaped", "https://example.org/%C3%BCber", NodeDataLinks.webUrl( "https://example.org/über" ) );
        ok &= eq( "trimmed", "http://example.org/x", NodeDataLinks.webUrl( " http://example.org/x " ) );
        ok &= eq( "ftp: no", null, NodeDataLinks.webUrl( "ftp://example.org/f" ) );
        ok &= eq( "file: no", null, NodeDataLinks.webUrl( "file:///etc/passwd" ) );
        ok &= eq( "javascript: no", null, NodeDataLinks.webUrl( "javascript:alert(1)" ) );
        ok &= eq( "an address inside a sentence: no", null, NodeDataLinks.webUrl( "see https://example.org" ) );
        ok &= eq( "no host: no", null, NodeDataLinks.webUrl( "https://" ) );
        ok &= eq( "an empty host (it passes the address pattern): no", null, NodeDataLinks.webUrl( "https:///etc/passwd" ) );
        ok &= eq( "a '|' beside an escape: the escape kept", "https://example.org/a%20b%7Cc",
                  NodeDataLinks.webUrl( "https://example.org/a%20b|c" ) );
        ok &= eq( "a lone '%' is escaped", "https://example.org/100%25", NodeDataLinks.webUrl( "https://example.org/100%" ) );
        ok &= eq( "brackets in a path are escaped", "https://example.org/a%5B1%5D",
                  NodeDataLinks.webUrl( "https://example.org/a[1]" ) );
        ok &= eq( "a second '#' is escaped", "https://example.org/a#b%23c", NodeDataLinks.webUrl( "https://example.org/a#b#c" ) );
        ok &= eq( "a character beyond the BMP", "https://example.org/%F0%9F%98%80",
                  NodeDataLinks.webUrl( "https://example.org/\ud83d\ude00" ) );
        ok &= eq( "a number: no", null, NodeDataLinks.webUrl( "120" ) );
        return ok;
    }

    /** An unpaired surrogate (an Auspice JSON "\\ud83d" is read as one) makes {@code URI} throw a
     *  NullPointerException: every rule answers "no link" instead. */
    private static boolean malformedUnicode() {
        boolean ok = true;
        try {
            ok &= eq( "web address", null, NodeDataLinks.webUrl( "https://example.org/\ud83d" ) );
            ok &= eq( "DOI", null, NodeDataLinks.doiUrl( "10.1000/\ud83d" ) );
            ok &= eq( "taxonomy id", null, NodeDataLinks.taxonomyIdUrl( "https://example.org/\udc00x", "ncbi" ) );
            ok &= eq( "accession", null, NodeDataLinks.accessionUrl( "\ud83dP38398", "uniprot", null ) );
            ok &= eq( "scientific name (form-encoded, no URI)", NCBI_TAX_NAME + "a%3Fb",
                      NodeDataLinks.scientificNameUrl( "a\ud83db", null ) );
        }
        catch ( final RuntimeException e ) {
            System.out.println( "  FAIL an unpaired surrogate threw " + e );
            return false;
        }
        return ok;
    }

    private static Taxonomy tax( final String id, final String provider, final String sci, final String code,
                                 final String common ) throws Exception {
        final Taxonomy t = new Taxonomy();
        if ( id != null ) {
            t.setIdentifier( new Identifier( id, provider ) );
        }
        if ( sci != null ) {
            t.setScientificName( sci );
        }
        if ( code != null ) {
            t.setTaxonomyCode( code );
        }
        if ( common != null ) {
            t.setCommonName( common );
        }
        return t;
    }

    private static boolean eq( final String what, final Object expected, final Object actual ) {
        if ( !Objects.equals( expected, actual ) ) {
            System.out.println( "  FAIL " + what + ": expected [" + expected + "], got [" + actual + "]" );
            return false;
        }
        return true;
    }
}
