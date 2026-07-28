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

package org.forester.archaeopteryx.tools;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.ws.seqdb.SequenceEntry;

/**
 * Unit tests for {@link LabelDataExtractor}: the UniProt-header parser, the fill-only-empty node writer,
 * and the whole-tree extraction. No network, no GUI.
 */
public final class LabelDataExtractorTest {

    public static void main( final String[] args ) {
        System.out.println( "LabelDataExtractor: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        // a full TrEMBL header
        final SequenceEntry e = LabelDataExtractor.parseUniProtHeader(
                "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Radical s-adenosyl methionine domain-containing protein "
                        + "OS=Fusarium phyllophilum OX=47803 GN=FPHYL_5758 PE=4 SV=1" );
        if ( !"A0A8H5JZG0".equals( e.getAccession() ) || !"A0A8H5JZG0_9HYPO".equals( e.getEntryName() )
                || !"Radical s-adenosyl methionine domain-containing protein".equals( e.getSequenceName() )
                || !"FPHYL_5758".equals( e.getGeneName() ) || !"Fusarium phyllophilum".equals( e.getOrganismName() )
                || !"47803".equals( e.getOrganismId() ) ) {
            return fail( "full header parsed wrong: " + e );
        }
        // OS value may contain spaces and parentheses (a strain); OX must still be read after it
        final SequenceEntry strain = LabelDataExtractor.parseUniProtHeader(
                "tr|A0A1X7RZU2|A0A1X7RZU2_ZYMT9 Radical SAM core domain-containing protein "
                        + "OS=Zymoseptoria tritici (strain ST99CH_3D7) OX=1276538 GN=ZT3D7_G7874 PE=4 SV=1" );
        if ( !"Zymoseptoria tritici (strain ST99CH_3D7)".equals( strain.getOrganismName() )
                || !"1276538".equals( strain.getOrganismId() ) ) {
            return fail( "organism with parentheses/strain parsed wrong: " + strain );
        }
        // leading ">" tolerated; missing GN tolerated
        final SequenceEntry nogn = LabelDataExtractor
                .parseUniProtHeader( ">sp|P0AEX9|MALE_ECOLI Maltose-binding periplasmic protein "
                        + "OS=Escherichia coli OX=83333 PE=1 SV=1" );
        if ( !"P0AEX9".equals( nogn.getAccession() ) || ( nogn.getGeneName() != null )
                || !"Escherichia coli".equals( nogn.getOrganismName() ) || !"83333".equals( nogn.getOrganismId() ) ) {
            return fail( "\">\"-prefixed / GN-less header parsed wrong: " + nogn );
        }
        // a bare sp|ACC|ENTRY with no description parses (accession present) but has no name/organism
        final SequenceEntry bare = LabelDataExtractor.parseUniProtHeader( "sp|Q12345|SOME_RAT" );
        if ( bare.isEmpty() || !"Q12345".equals( bare.getAccession() ) || ( bare.getSequenceName() != null )
                || ( bare.getOrganismName() != null ) ) {
            return fail( "bare sp|ACC|ENTRY parsed wrong: " + bare );
        }
        // non-UniProt labels are ignored (EMPTY)
        if ( !LabelDataExtractor.parseUniProtHeader( "Homo_sapiens" ).isEmpty()
                || !LabelDataExtractor.parseUniProtHeader( "Felis catus" ).isEmpty()
                || !LabelDataExtractor.parseUniProtHeader( "" ).isEmpty() ) {
            return fail( "a non-UniProt label must parse to EMPTY" );
        }

        // applyToNode: fills empty fields, renames the node to the accession
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Some protein OS=Fusarium phyllophilum OX=47803 GN=FPHYL_5758 PE=4 SV=1" );
        if ( !LabelDataExtractor.applyToNode( n, LabelDataExtractor.parseUniProtHeader( n.getName() ) ) ) {
            return fail( "applyToNode must report a change for a header label" );
        }
        if ( !"A0A8H5JZG0".equals( n.getName() ) ) {
            return fail( "node name must be shortened to the accession; got " + n.getName() );
        }
        final Sequence ws = n.getNodeData().getSequence();
        if ( ( ws == null ) || ( ws.getAccession() == null ) || !"A0A8H5JZG0".equals( ws.getAccession().getValue() )
                || !"Some protein".equals( ws.getName() ) || !"FPHYL_5758".equals( ws.getGeneName() ) ) {
            return fail( "sequence fields not written: " + ws );
        }
        final Taxonomy wt = n.getNodeData().getTaxonomy();
        if ( ( wt == null ) || !"Fusarium phyllophilum".equals( wt.getScientificName() ) || ( wt.getIdentifier() == null )
                || !"47803".equals( wt.getIdentifier().getValue() ) || !"ncbi".equals( wt.getIdentifier().getProvider() ) ) {
            return fail( "taxonomy fields not written: " + wt );
        }

        // applyToNode must NOT clobber a curated sequence name or scientific name
        final PhylogenyNode cur = new PhylogenyNode();
        cur.setName( "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Some protein OS=Fusarium phyllophilum OX=47803 GN=FPHYL_5758 PE=4 SV=1" );
        final Sequence cs = new Sequence();
        cs.setName( "my curated name" );
        cur.getNodeData().addSequence( cs );
        final Taxonomy ct = new Taxonomy();
        ct.setScientificName( "Curated organism" );
        cur.getNodeData().addTaxonomy( ct );
        LabelDataExtractor.applyToNode( cur, LabelDataExtractor.parseUniProtHeader( cur.getName() ) );
        if ( !"my curated name".equals( cur.getNodeData().getSequence().getName() )
                || !"Curated organism".equals( cur.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "applyToNode must not overwrite curated fields" );
        }
        // ...but it still fills the empties (gene, tax-id) and renames
        if ( !"FPHYL_5758".equals( cur.getNodeData().getSequence().getGeneName() )
                || ( cur.getNodeData().getTaxonomy().getIdentifier() == null ) || !"A0A8H5JZG0".equals( cur.getName() ) ) {
            return fail( "applyToNode must still fill the empty fields and rename" );
        }

        // whole-tree extract + the proactive-offer gate
        final Phylogeny tree = headerTree();
        // 2 of 3 tips are headers -> a majority -> the load-time offer fires
        if ( !LabelDataExtractor.mostLabelsParsable( tree ) ) {
            return fail( "mostLabelsParsable must be true when a majority of tips are UniProt headers" );
        }
        final int changed = LabelDataExtractor.extract( tree );
        if ( changed != 2 ) {
            return fail( "extract must change exactly the 2 header tips (the plain tip is left alone); got " + changed );
        }
        if ( LabelDataExtractor.extract( tree ) != 0 ) {
            return fail( "a second extract must be a no-op (everything already filled)" );
        }
        final Phylogeny plain = plainTree();
        if ( LabelDataExtractor.mostLabelsParsable( plain ) || ( LabelDataExtractor.extract( plain ) != 0 ) ) {
            return fail( "a tree of plain names must have nothing to extract and must not trigger the offer" );
        }
        // a single stray header among many plain names must NOT trip the proactive offer (minority)
        if ( LabelDataExtractor.mostLabelsParsable( minorityHeaderTree() ) ) {
            return fail( "a minority of header tips must not trigger the load-time offer" );
        }
        return testGenbank();
    }

    private static boolean testGenbank() {
        // protein defline: accession + description + organism from the trailing [..]
        final SequenceEntry p = LabelDataExtractor
                .parseGenbankDefline( ">NP_788278.1 death executioner Bcl-2 [Drosophila melanogaster]" );
        if ( !"NP_788278.1".equals( p.getAccession() ) || !"death executioner Bcl-2".equals( p.getSequenceName() )
                || !"Drosophila melanogaster".equals( p.getOrganismName() ) || ( p.getOrganismId() != null )
                || ( p.getGeneName() != null ) ) {
            return fail( "protein defline parsed wrong: " + p );
        }
        // 3-letter GenBank protein accession
        final SequenceEntry caa = LabelDataExtractor.parseGenbankDefline( ">CAA01637.1 spike [Canine coronavirus]" );
        if ( !"CAA01637.1".equals( caa.getAccession() ) || !"spike".equals( caa.getSequenceName() )
                || !"Canine coronavirus".equals( caa.getOrganismName() ) ) {
            return fail( "GenBank protein accession defline parsed wrong: " + caa );
        }
        // organism with strain (spaces + dots) and a description containing a comma
        final SequenceEntry strain = LabelDataExtractor.parseGenbankDefline(
                ">WP_482516138.1 carbohydrate kinase family protein, partial [Escherichia coli str. K-12 substr. MG1655]" );
        if ( !"carbohydrate kinase family protein, partial".equals( strain.getSequenceName() )
                || !"Escherichia coli str. K-12 substr. MG1655".equals( strain.getOrganismName() ) ) {
            return fail( "defline with strain/comma parsed wrong: " + strain );
        }
        // nucleotide defline: NO bracket -> accession + description (PREDICTED: stripped), NO organism guessed
        final SequenceEntry nuc = LabelDataExtractor.parseGenbankDefline(
                ">XM_006503748.2 PREDICTED: Mus musculus ketohexokinase (Khk), transcript variant X1, mRNA" );
        if ( !"XM_006503748.2".equals( nuc.getAccession() ) || ( nuc.getOrganismName() != null )
                || !"Mus musculus ketohexokinase (Khk), transcript variant X1, mRNA".equals( nuc.getSequenceName() ) ) {
            return fail( "nucleotide defline must give accession + description but no guessed organism: " + nuc );
        }
        // not a defline: a leading token that is not an NCBI accession
        if ( !LabelDataExtractor.parseGenbankDefline( "Homo sapiens cytochrome c" ).isEmpty()
                || !LabelDataExtractor.parseGenbankDefline( "just_a_label" ).isEmpty() ) {
            return fail( "a non-accession leading token must not parse as a GenBank defline" );
        }
        // parseHeader dispatches: UniProt header stays UniProt (has tax-id), GenBank defline routes to GenBank
        if ( !"47803".equals( LabelDataExtractor
                .parseHeader( "tr|A0A8H5JZG0|X_9HYPO P OS=Fusarium phyllophilum OX=47803 GN=G PE=4 SV=1" )
                .getOrganismId() ) ) {
            return fail( "parseHeader must keep a UniProt header on the UniProt path" );
        }
        if ( !"NP_788278.1".equals(
                LabelDataExtractor.parseHeader( "NP_788278.1 death executioner Bcl-2 [Drosophila melanogaster]" )
                        .getAccession() ) ) {
            return fail( "parseHeader must route a GenBank defline to the GenBank path" );
        }
        // applyToNode on a GenBank tip: renames to the accession, writes the RefSeq source + organism, no tax-id
        final PhylogenyNode n = leaf( ">NP_788278.1 death executioner Bcl-2 [Drosophila melanogaster]" );
        if ( !LabelDataExtractor.applyToNode( n, LabelDataExtractor.parseHeader( n.getName() ) )
                || !"NP_788278.1".equals( n.getName() ) ) {
            return fail( "GenBank applyToNode must rename the node to its accession" );
        }
        final Sequence gs = n.getNodeData().getSequence();
        if ( ( gs.getAccession() == null ) || !"refseq".equals( gs.getAccession().getSource() )
                || !"death executioner Bcl-2".equals( gs.getName() ) ) {
            return fail( "GenBank accession must be written with the refseq source + the description: " + gs );
        }
        final Taxonomy gt = n.getNodeData().getTaxonomy();
        if ( ( gt == null ) || !"Drosophila melanogaster".equals( gt.getScientificName() )
                || ( gt.getIdentifier() != null ) ) {
            return fail( "GenBank organism must be written as scientific name with NO tax-id: " + gt );
        }
        // a tree of GenBank deflines triggers the proactive offer
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( leaf( "NP_788278.1 death executioner Bcl-2 [Drosophila melanogaster]" ) );
        root.addAsChild( leaf( "CAA01637.1 spike [Canine coronavirus]" ) );
        final Phylogeny gb = new Phylogeny();
        gb.setRoot( root );
        gb.externalNodesHaveChanged();
        if ( !LabelDataExtractor.mostLabelsParsable( gb ) || ( LabelDataExtractor.extract( gb ) != 2 ) ) {
            return fail( "a GenBank-defline tree must be offered and fully extracted" );
        }
        return testNcbiCds();
    }

    private static boolean testNcbiCds() {
        // full CDS-from-genomic defline: accession = protein_id, description = protein, gene = gene, NO organism
        final SequenceEntry full = LabelDataExtractor.parseNcbiCdsDefline(
                "lcl|MT726045.1_prot_QTJ30158.1_7 [gene=ORF7a] [protein=ORF7a protein] [protein_id=QTJ30158.1] "
                        + "[location=27187..27543] [gbkey=CDS]" );
        if ( !"QTJ30158.1".equals( full.getAccession() ) || !"ORF7a protein".equals( full.getSequenceName() )
                || !"ORF7a".equals( full.getGeneName() ) || ( full.getOrganismName() != null )
                || ( full.getOrganismId() != null ) ) {
            return fail( "full lcl| CDS defline parsed wrong: " + full );
        }
        // extra qualifiers (locus_tag, db_xref) are ignored; a RefSeq protein_id is taken verbatim
        final SequenceEntry extra = LabelDataExtractor.parseNcbiCdsDefline(
                "lcl|NC_014470.1_prot_YP_003858589.1_8 [gene=ORF7a] [locus_tag=BtCoVBM48_gp7] "
                        + "[db_xref=GeneID:9714829] [protein=ORF7a protein] [protein_id=YP_003858589.1] "
                        + "[location=27188..27544] [gbkey=CDS]" );
        if ( !"YP_003858589.1".equals( extra.getAccession() ) || !"ORF7a protein".equals( extra.getSequenceName() )
                || !"ORF7a".equals( extra.getGeneName() ) ) {
            return fail( "lcl| CDS defline with extra qualifiers parsed wrong: " + extra );
        }
        // no [gene=] -> gene null, but still parses (accession + description)
        final SequenceEntry nogene = LabelDataExtractor.parseNcbiCdsDefline(
                "lcl|MZ190137.1_prot_QVN46564.1_8 [protein=ORF7a protein] [protein_id=QVN46564.1] "
                        + "[location=27116..27472] [gbkey=CDS]" );
        if ( !"QVN46564.1".equals( nogene.getAccession() ) || !"ORF7a protein".equals( nogene.getSequenceName() )
                || ( nogene.getGeneName() != null ) ) {
            return fail( "lcl| CDS defline without a gene parsed wrong: " + nogene );
        }
        // a leading ">" is tolerated
        if ( !"QVN46564.1".equals( LabelDataExtractor
                .parseNcbiCdsDefline( ">lcl|MZ190137.1_prot_QVN46564.1_8 [protein_id=QVN46564.1] [gbkey=CDS]" )
                .getAccession() ) ) {
            return fail( "a \">\"-prefixed lcl| CDS defline must parse" );
        }
        // an lcl| identifier with NO [key=value] qualifiers is not this dialect
        if ( !LabelDataExtractor.parseNcbiCdsDefline( "lcl|MT726045.1_prot_QTJ30158.1_7" ).isEmpty() ) {
            return fail( "a bare lcl| identifier (no qualifiers) must parse to EMPTY" );
        }
        // a plain GenBank protein defline is NOT an lcl| defline -> EMPTY here (so it does not shadow parseGenbankDefline)
        if ( !LabelDataExtractor.parseNcbiCdsDefline( "NP_788278.1 death executioner Bcl-2 [Drosophila melanogaster]" )
                .isEmpty() ) {
            return fail( "a non-lcl| defline must not parse as an NCBI CDS defline" );
        }
        // parseHeader routes an lcl| CDS defline to the CDS path (protein_id accession), NOT the mis-reading
        // GenBank fallback (which would have turned the [..] blob into a description and "[gbkey=CDS]" into an organism)
        final SequenceEntry viaHeader = LabelDataExtractor.parseHeader(
                "lcl|MT726045.1_prot_QTJ30158.1_7 [gene=ORF7a] [protein=ORF7a protein] [protein_id=QTJ30158.1] "
                        + "[location=27187..27543] [gbkey=CDS]" );
        if ( !"QTJ30158.1".equals( viaHeader.getAccession() ) || !"ORF7a protein".equals( viaHeader.getSequenceName() )
                || ( viaHeader.getOrganismName() != null ) ) {
            return fail( "parseHeader must route an lcl| CDS defline to the CDS path: " + viaHeader );
        }
        // applyToNode on a CDS tip: renames to the protein_id, writes the description + gene, writes NO taxonomy
        final PhylogenyNode n = leaf(
                "lcl|MT726045.1_prot_QTJ30158.1_7 [gene=ORF7a] [protein=ORF7a protein] [protein_id=QTJ30158.1] "
                        + "[location=27187..27543] [gbkey=CDS]" );
        if ( !LabelDataExtractor.applyToNode( n, LabelDataExtractor.parseHeader( n.getName() ) )
                || !"QTJ30158.1".equals( n.getName() ) ) {
            return fail( "CDS applyToNode must rename the node to its protein_id accession" );
        }
        final Sequence cs = n.getNodeData().getSequence();
        if ( ( cs == null ) || ( cs.getAccession() == null ) || !"QTJ30158.1".equals( cs.getAccession().getValue() )
                || !"ORF7a protein".equals( cs.getName() ) || !"ORF7a".equals( cs.getGeneName() ) ) {
            return fail( "CDS sequence fields (accession, description, gene) not written: " + cs );
        }
        if ( n.getNodeData().isHasTaxonomy() ) {
            return fail( "a CDS defline carries no organism, so no taxonomy must be written" );
        }
        // a tree of lcl| CDS deflines triggers the proactive offer and is fully extracted
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( leaf( "lcl|MT726045.1_prot_QTJ30158.1_7 [protein=ORF7a protein] [protein_id=QTJ30158.1] [gbkey=CDS]" ) );
        root.addAsChild( leaf( "lcl|PV584046.1_prot_XSG41275.1_7 [protein=ORF7a protein] [protein_id=XSG41275.1] [gbkey=CDS]" ) );
        final Phylogeny cds_tree = new Phylogeny();
        cds_tree.setRoot( root );
        cds_tree.externalNodesHaveChanged();
        if ( !LabelDataExtractor.mostLabelsParsable( cds_tree ) || ( LabelDataExtractor.extract( cds_tree ) != 2 ) ) {
            return fail( "an lcl| CDS-defline tree must be offered and fully extracted" );
        }
        return true;
    }

    // (header_a, header_b, plain): two UniProt-header tips and one ordinary name
    private static Phylogeny headerTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( leaf( "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO P1 OS=Fusarium phyllophilum OX=47803 GN=G1 PE=4 SV=1" ) );
        root.addAsChild( leaf( "sp|P0AEX9|MALE_ECOLI P2 OS=Escherichia coli OX=83333 GN=malE PE=1 SV=1" ) );
        root.addAsChild( leaf( "just_a_label" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    // one header tip among three plain ones -> a minority (must not trip the proactive offer)
    private static Phylogeny minorityHeaderTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( leaf( "sp|P0AEX9|MALE_ECOLI P OS=Escherichia coli OX=83333 PE=1 SV=1" ) );
        root.addAsChild( leaf( "alpha" ) );
        root.addAsChild( leaf( "beta" ) );
        root.addAsChild( leaf( "gamma" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static Phylogeny plainTree() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( leaf( "Homo sapiens" ) );
        root.addAsChild( leaf( "Mus musculus" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  LabelDataExtractor test failed: " + msg );
        return false;
    }

    private LabelDataExtractorTest() {
    }
}
