// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.ws.seqdb;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import org.forester.go.GoTerm;
import org.forester.phylogeny.data.Accession;
import org.forester.util.ForesterUtil;

public final class EbiDbEntry implements SequenceDatabaseEntry {

    public static SequenceDatabaseEntry createInstanceFromPlainText( final List<String> lines ) {
        final EbiDbEntry e = new EbiDbEntry();
        for( final String line : lines ) {
            if ( line.startsWith( "PA" ) ) {
                e.setPA( SequenceDbWsTools.extractFrom( line, "PA" ) );
            }
            else if ( line.startsWith( "DE" ) ) {
                e.setDe( SequenceDbWsTools.extractFrom( line, "DE" ) );
            }
            else if ( line.startsWith( "OS" ) ) {
                if ( line.indexOf( "(" ) > 0 ) {
                    e.setOs( SequenceDbWsTools.extractFromTo( line, "OS", "(" ) );
                }
                else {
                    e.setOs( SequenceDbWsTools.extractFrom( line, "OS" ) );
                }
            }
            else if ( line.startsWith( "OX" ) ) {
                if ( line.indexOf( "NCBI_TaxID=" ) > 0 ) {
                    e.setTaxId( SequenceDbWsTools.extractFromTo( line, "NCBI_TaxID=", ";" ) );
                }
            }
        }
        return e;
    }

    public static SequenceDatabaseEntry createInstanceFromPlainTextForRefSeq( final List<String> lines ) {
        final Pattern X_PATTERN = Pattern.compile( "^[A-Z]+" );
        final Pattern chromosome_PATTERN = Pattern.compile( "\\s+/chromosome=\"(\\w+)\"" );
        final Pattern map_PATTERN = Pattern.compile( "\\s+/map=\"([\\w+\\.])\"" );
        final Pattern gene_PATTERN = Pattern.compile( "\\s+/gene=\"(.+)\"" );
        final Pattern mim_xref_PATTERN = Pattern.compile( "\\s+/db_xref=\"MIM:(\\d+)\"" );
        final Pattern taxon_xref_PATTERN = Pattern.compile( "\\s+/db_xref=\"taxon:(\\d+)\"" );
        final Pattern interpro_PATTERN = Pattern.compile( "\\s+/db_xref=\"InterPro:(IP\\d+)\"" );
        final Pattern uniprot_PATTERN = Pattern.compile( "\\s+/db_xref=\"UniProtKB/TrEMBL:(\\w+)\"" );
        final EbiDbEntry e = new EbiDbEntry();
        final StringBuilder def = new StringBuilder();
        boolean in_def = false;
        boolean in_features = false;
        boolean in_source = false;
        boolean in_gene = false;
        boolean in_cds = false;
        boolean in_protein = false;
        for( final String line : lines ) {
            if ( line.startsWith( "ACCESSION " ) ) {
                e.setPA( SequenceDbWsTools.extractFrom( line, "ACCESSION" ) );
                in_def = false;
            }
            else if ( line.startsWith( "DEFINITION " ) ) {
                if ( line.indexOf( "[" ) > 0 ) {
                    def.append( SequenceDbWsTools.extractFromTo( line, "DEFINITION", "[" ) );
                }
                else if ( line.indexOf( "." ) > 0 ) {
                    def.append( SequenceDbWsTools.extractFromTo( line, "DEFINITION", "." ) );
                }
                else {
                    def.append( SequenceDbWsTools.extractFrom( line, "DEFINITION" ) );
                }
                in_def = true;
            }
            else if ( line.startsWith( "  ORGANISM " ) ) {
                if ( line.indexOf( "(" ) > 0 ) {
                    e.setOs( SequenceDbWsTools.extractFromTo( line, "  ORGANISM", "(" ) );
                }
                else {
                    e.setOs( SequenceDbWsTools.extractFrom( line, "  ORGANISM" ) );
                }
                //  in_def = false;
            }
            else if ( line.startsWith( " " ) && in_def ) {
                def.append( " " );
                if ( line.indexOf( "[" ) > 0 ) {
                    def.append( SequenceDbWsTools.extractTo( line, "[" ) );
                }
                else if ( line.indexOf( "." ) > 0 ) {
                    def.append( SequenceDbWsTools.extractTo( line, "." ) );
                }
                else {
                    def.append( line.trim() );
                }
            }
            else {
                in_def = false;
            }
            if ( X_PATTERN.matcher( line ).find() ) {
                in_features = false;
                in_source = false;
                in_gene = false;
                in_cds = false;
                in_protein = false;
                // in_def = false;
            }
            if ( line.startsWith( "FEATURES " ) ) {
                in_features = true;
            }
            if ( in_features && line.startsWith( "     source " ) ) {
                in_source = true;
                in_gene = false;
                in_cds = false;
                in_protein = false;
            }
            if ( in_features && line.startsWith( "     gene " ) ) {
                in_source = false;
                in_gene = true;
                in_cds = false;
                in_protein = false;
            }
            if ( in_features && line.startsWith( "     CDS " ) ) {
                in_source = false;
                in_gene = false;
                in_cds = true;
                in_protein = false;
            }
            if ( in_features && line.startsWith( "     Protein " ) ) {
                in_source = false;
                in_gene = false;
                in_cds = false;
                in_protein = true;
            }
        }
        if ( def.length() > 0 ) {
            e.setDe( def.toString().trim() );
        }
        return e;
    }
    // FIXME actually this is NCBI entry
    //http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/emb/AAR37336/
    private String               _pa;
    private String               _de;
    private String               _os;
    private String               _tax_id;
    private String               _symbol;
    private String               _provider;
    private ArrayList<Accession> _cross_references;
    private String               _gene_name;

    // TODO  PUBMED   15798186
    //TODO  (FEATURES) 
    // source /db_xref="taxon:9606"
    // gene            1..2881  
    // /gene="RBM39" 
    //
    // /db_xref="MIM:604739"  
    // CDS
    // /gene="RBM39"
    // /db_xref="MIM:604739"
    // /db_xref="InterPro:IPR002475"
    // /product="Bcl-2"
    // /db_xref="UniProtKB/TrEMBL:Q5J7V1" <- reparse?
    //
    // Protein
    /*
    LOCUS       NM_184234               2881 bp    mRNA    linear   PRI 16-JUN-2013
    DEFINITION  Homo sapiens RNA binding motif protein 39 (RBM39), transcript
            variant 1, mRNA.
    ACCESSION   NM_184234
    VERSION     NM_184234.2  GI:336176061
    KEYWORDS    RefSeq.
    SOURCE      Homo sapiens (human)
    ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
            Catarrhini; Hominidae; Homo.
    REFERENCE   1  (bases 1 to 2881)
    AUTHORS   Sillars-Hardebol,A.H., Carvalho,B., Belien,J.A., de Wit,M.,
            Delis-van Diemen,P.M., Tijssen,M., van de Wiel,M.A., Ponten,F.,
            Meijer,G.A. and Fijneman,R.J.
    TITLE     CSE1L, DIDO1 and RBM39 in colorectal adenoma to carcinoma
            progression
    JOURNAL   Cell Oncol (Dordr) 35 (4), 293-300 (2012)
    PUBMED   22711543
    REMARK    GeneRIF: Data show that CSE1L, DIDO1 and RBM39 mRNA expression
            levels correlated with chromosome 20q DNA copy number status.
    REFERENCE   2  (bases 1 to 2881)
    AUTHORS   Huang,G., Zhou,Z., Wang,H. and Kleinerman,E.S.
    TITLE     CAPER-alpha alternative splicing regulates the expression of
            vascular endothelial growth factor(1)(6)(5) in Ewing sarcoma cells
    JOURNAL   Cancer 118 (8), 2106-2116 (2012)
    PUBMED   22009261
    REMARK    GeneRIF: Increased VEGF(165) expression is secondary to the
            down-regulation of CAPER-alpha by EWS/FLI-1. CAPER-alpha mediates
            alternative splicing and controls the shift from VEGF(189) to
            VEGF(165) .
    REFERENCE   3  (bases 1 to 2881)
    AUTHORS   Han,B., Stockwin,L.H., Hancock,C., Yu,S.X., Hollingshead,M.G. and
            Newton,D.L.
    TITLE     Proteomic analysis of nuclei isolated from cancer cell lines
            treated with indenoisoquinoline NSC 724998, a novel topoisomerase I
            inhibitor
    JOURNAL   J. Proteome Res. 9 (8), 4016-4027 (2010)
    PUBMED   20515076
    REMARK    Erratum:[J Proteome Res. 2011 Apr 1;10(4):2128]
    REFERENCE   4  (bases 1 to 2881)
    AUTHORS   Zhang,J.Y., Looi,K.S. and Tan,E.M.
    TITLE     Identification of tumor-associated antigens as diagnostic and
            predictive biomarkers in cancer
    JOURNAL   Methods Mol. Biol. 520, 1-10 (2009)
    PUBMED   19381943
    REFERENCE   5  (bases 1 to 2881)
    AUTHORS   Dutta,J., Fan,G. and Gelinas,C.
    TITLE     CAPERalpha is a novel Rel-TAD-interacting factor that inhibits
            lymphocyte transformation by the potent Rel/NF-kappaB oncoprotein
            v-Rel
    JOURNAL   J. Virol. 82 (21), 10792-10802 (2008)
    PUBMED   18753212
    REMARK    GeneRIF: this study identifies CAPERalpha (RNA binding motif
            protein 39) as a new transcriptional coregulator for v-Rel and
            reveals an important role in modulating Rel's oncogenic activity.
    REFERENCE   6  (bases 1 to 2881)
    AUTHORS   Cazalla,D., Newton,K. and Caceres,J.F.
    TITLE     A novel SR-related protein is required for the second step of
            Pre-mRNA splicing
    JOURNAL   Mol. Cell. Biol. 25 (8), 2969-2980 (2005)
    PUBMED   15798186
    REFERENCE   7  (bases 1 to 2881)
    AUTHORS   Dowhan,D.H., Hong,E.P., Auboeuf,D., Dennis,A.P., Wilson,M.M.,
            Berget,S.M. and O'Malley,B.W.
    TITLE     Steroid hormone receptor coactivation and alternative RNA splicing
            by U2AF65-related proteins CAPERalpha and CAPERbeta
    JOURNAL   Mol. Cell 17 (3), 429-439 (2005)
    PUBMED   15694343
    REFERENCE   8  (bases 1 to 2881)
    AUTHORS   Sun,N.N., Fastje,C.D., Wong,S.S., Sheppard,P.R., Macdonald,S.J.,
            Ridenour,G., Hyde,J.D. and Witten,M.L.
    TITLE     Dose-dependent transcriptome changes by metal ores on a human acute
            lymphoblastic leukemia cell line
    JOURNAL   Toxicol Ind Health 19 (7-10), 157-163 (2003)
    PUBMED   15747776
    REMARK    GeneRIF: 10 genes were down-regulated following treatment of the
            T-ALL cells with 0.15 and 1.5 microg/mL of metal ores at 72 h
    REFERENCE   9  (bases 1 to 2881)
    AUTHORS   Jung,D.J., Na,S.Y., Na,D.S. and Lee,J.W.
    TITLE     Molecular cloning and characterization of CAPER, a novel
            coactivator of activating protein-1 and estrogen receptors
    JOURNAL   J. Biol. Chem. 277 (2), 1229-1234 (2002)
    PUBMED   11704680
    REMARK    GeneRIF: This paper describes the mouse gene.
    REFERENCE   10 (bases 1 to 2881)
    AUTHORS   Imai,H., Chan,E.K., Kiyosawa,K., Fu,X.D. and Tan,E.M.
    TITLE     Novel nuclear autoantigen with splicing factor motifs identified
            with antibody from hepatocellular carcinoma
    JOURNAL   J. Clin. Invest. 92 (5), 2419-2426 (1993)
    PUBMED   8227358
    COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The
            reference sequence was derived from DC346351.1, BC141835.1 and
            C75555.1.
            On Jun 16, 2011 this sequence version replaced gi:35493810.
            
            Summary: This gene encodes a member of the U2AF65 family of
            proteins. The encoded protein is found in the nucleus, where it
            co-localizes with core spliceosomal proteins. It has been shown to
            play a role in both steroid hormone receptor-mediated transcription
            and alternative splicing, and it is also a transcriptional
            coregulator of the viral oncoprotein v-Rel. Multiple transcript
            variants have been observed for this gene. A related pseudogene has
            been identified on chromosome X. [provided by RefSeq, Aug 2011].
            
            Transcript Variant: This variant (1) encodes the longest isoform
            (a, also called CC1.4).
            
            Publication Note:  This RefSeq record includes a subset of the
            publications that are available for this gene. Please see the Gene
            record to access additional publications.
            
            ##Evidence-Data-START##
            Transcript exon combination :: BC141835.1, L10911.1 [ECO:0000332]
            RNAseq introns              :: mixed/partial sample support
                                           ERS025081, ERS025082 [ECO:0000350]
            ##Evidence-Data-END##
            COMPLETENESS: complete on the 3' end.
    PRIMARY     REFSEQ_SPAN         PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP
            1-578               DC346351.1         3-580
            579-2872            BC141835.1         429-2722
            2873-2881           C75555.1           1-9                 c
    FEATURES             Location/Qualifiers
     source          1..2881
                     /organism="Homo sapiens"
                     /mol_type="mRNA"
                     /db_xref="taxon:9606"
                     /chromosome="20"
                     /map="20q11.22"
     gene            1..2881
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /note="RNA binding motif protein 39"
                     /db_xref="GeneID:9584"
                     /db_xref="HGNC:15923"
                     /db_xref="HPRD:09201"
                     /db_xref="MIM:604739"
     exon            1..396
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /inference="alignment:Splign:1.39.8"
     STS             35..262
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /standard_name="REN58946"
                     /db_xref="UniSTS:383746"
     misc_feature    221..223
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /note="upstream in-frame stop codon"
     STS             299..453
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /standard_name="G64285"
                     /db_xref="UniSTS:158667"
     exon            397..460
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /inference="alignment:Splign:1.39.8"
     CDS             410..2002
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /note="isoform a is encoded by transcript variant 1;
                     coactivator of activating protein-1 and estrogen
                     receptors; functional spliceosome-associated protein 59;
                     RNA-binding region (RNP1, RRM) containing 2;
                     hepatocellular carcinoma protein 1; splicing factor HCC1"
                     /codon_start=1
                     /product="RNA-binding protein 39 isoform a"
                     /protein_id="NP_909122.1"
                     /db_xref="GI:35493811"
                     /db_xref="CCDS:CCDS13266.1"
                     /db_xref="GeneID:9584"
                     /db_xref="HGNC:15923"
                     /db_xref="HPRD:09201"
                     /db_xref="MIM:604739"
                     /translation="MADDIDIEAMLEAPYKKDENKLSSANGHEERSKKRKKSKSRSRS
                     HERKRSKSKERKRSRDRERKKSKSRERKRSRSKERRRSRSRSRDRRFRGRYRSPYSGP
                     KFNSAIRGKIGLPHSIKLSRRRSRSKSPFRKDKSPVREPIDNLTPEERDARTVFCMQL
                     AARIRPRDLEEFFSTVGKVRDVRMISDRNSRRSKGIAYVEFVDVSSVPLAIGLTGQRV
                     LGVPIIVQASQAEKNRAAAMANNLQKGSAGPMRLYVGSLHFNITEDMLRGIFEPFGRI
                     ESIQLMMDSETGRSKGYGFITFSDSECAKKALEQLNGFELAGRPMKVGHVTERTDASS
                     ASSFLDSDELERTGIDLGTTGRLQLMARLAEGTGLQIPPAAQQALQMSGSLAFGAVAE
                     FSFVIDLQTRLSQQTEASALAAAASVQPLATQCFQLSNMFNPQTEEEVGWDTEIKDDV
                     IEECNKHGGVIHIYVDKNSAQGNVYVKCPSIAAAIAAVNALHGRWFAGKMITAAYVPL
                     PTYHNLFPDSMTATQLLVPSRR"
     misc_feature    413..415
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /experiment="experimental evidence, no additional details
                     recorded"
                     /note="N-acetylalanine; propagated from
                     UniProtKB/Swiss-Prot (Q14498.2); acetylation site"
     
     exon            461..510
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /inference="alignment:Splign:1.39.8"
    
     exon            1902..2874
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /inference="alignment:Splign:1.39.8"
     STS             1956..2182
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /standard_name="REN58786"
                     /db_xref="UniSTS:383586"
     STS             2104..2148
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /standard_name="D19S1033"
                     /db_xref="UniSTS:154759"
     STS             2145..2400
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
                     /standard_name="REN58785"
                     /db_xref="UniSTS:383585"
    
     polyA_signal    2851..2856
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
     polyA_site      2874
                     /gene="RBM39"
                     /gene_synonym="CAPER; CAPERalpha; FSAP59; HCC1; RNPC2"
    ORIGIN      
        1 atttggagct tggggcagct tctcgcgaga gcccgtgctg agggctctgt gaggccccgt
       61 gtgtttgtgt gtgtgtatgt gtgctggtga atgtgagtac agggaagcag cggccgccat
      121 ttcagggagc ttgtcgacgc tgtcgcaggg gtggatcctg agctgccgaa gccgccgtcc
      181 tgctctcccg cgtgggcttc tctaattcca ttgttttttt tagattctct cgggcctagc
      241 cgtccttgga acccgatatt cgggctgggc ggttccgcgg cctgggccta ggggcttaac
    
    
    
    */
    private EbiDbEntry() {
    }

    private void addCrossReference( final Accession accession ) {
        if ( _cross_references == null ) {
            _cross_references = new ArrayList<Accession>();
        }
        System.out.println( "XREF ADDED: " + accession );
        _cross_references.add( accession );
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    @Override
    public String getAccession() {
        return _pa;
    }

    @Override
    public List<Accession> getCrossReferences() {
        return _cross_references;
    }

    @Override
    public String getGeneName() {
        return _gene_name;
    }

    @Override
    public List<GoTerm> getGoTerms() {
        return null;
    }

    @Override
    public String getProvider() {
        return _provider;
    }

    @Override
    public String getSequenceName() {
        return _de;
    }

    @Override
    public String getSequenceSymbol() {
        return _symbol;
    }

    @Override
    public String getTaxonomyIdentifier() {
        return _tax_id;
    }

    @Override
    public String getTaxonomyScientificName() {
        return _os;
    }

    @Override
    public boolean isEmpty() {
        return ( ForesterUtil.isEmpty( getAccession() ) && ForesterUtil.isEmpty( getSequenceName() )
                && ForesterUtil.isEmpty( getTaxonomyScientificName() )
                && ForesterUtil.isEmpty( getTaxonomyIdentifier() ) && ForesterUtil.isEmpty( getSequenceSymbol() ) );
    }

    private void setDe( final String rec_name ) {
        if ( _de == null ) {
            _de = rec_name;
        }
    }

    private void setGeneName( final String gene_name ) {
        if ( _gene_name == null ) {
            _gene_name = gene_name;
        }
    }

    private void setOs( final String os ) {
        if ( _os == null ) {
            _os = os;
        }
    }

    private void setPA( final String pa ) {
        if ( _pa == null ) {
            _pa = pa;
        }
    }

    public void setProvider( final String provider ) {
        _provider = provider;
    }

    private void setTaxId( final String tax_id ) {
        if ( _tax_id == null ) {
            _tax_id = tax_id;
        }
    }
}
