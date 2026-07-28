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

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;
import org.forester.util.SequenceAccessionTools;
import org.forester.ws.seqdb.SequenceEntry;

/**
 * Pulls the structured data out of FASTA-header node labels and writes it into the proper node fields
 * (sequence accession + name + gene, taxonomy scientific name + NCBI tax-id), then shortens the node
 * name to the bare accession. Two header dialects are recognized:
 * <ul>
 * <li><b>UniProt</b>: {@code tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Radical s-adenosyl methionine ... protein
 * OS=Fusarium phyllophilum OX=47803 GN=FPHYL_5758 PE=4 SV=1} -- yields accession, protein name, gene,
 * organism scientific name <i>and</i> NCBI tax-id.</li>
 * <li><b>GenBank/RefSeq</b>: {@code NP_788278.1 death executioner Bcl-2 [Drosophila melanogaster]} --
 * yields accession + description, and the organism <i>only</i> from a trailing {@code [..]} bracket (no
 * tax-id, no gene). A nucleotide defline that embeds the organism in prose with no bracket (e.g.
 * {@code XM_006503748.2 PREDICTED: Mus musculus ketohexokinase (Khk), ... mRNA}) yields just accession +
 * description; the organism is left to the online resolver rather than guessed from the text.</li>
 * <li><b>NCBI CDS-from-genomic</b>: {@code lcl|MT726045.1_prot_QTJ30158.1_7 [gene=ORF7a] [protein=ORF7a
 * protein] [protein_id=QTJ30158.1] [location=27187..27543] [gbkey=CDS]} -- yields accession (the
 * {@code protein_id}), description (the {@code protein} qualifier) and gene. Carries no organism (it is not
 * on a CDS defline).</li>
 * </ul>
 *
 * <p>This is the offline counterpart of "Fetch Sequence &amp; Taxonomic Data": it parses what is already
 * in the label rather than hitting the network, turning an unwieldy, not-yet-publishable tree into one
 * whose data drives display, colorize-by-rank and clade annotation.
 *
 * <p>Existing curated fields are never overwritten -- only empty ones are filled -- so a second run is a
 * no-op and a partially-annotated tree is respected. The parsers are pure and the node-writer takes a
 * single node, so the whole thing is unit-testable with no GUI.
 */
public final class LabelDataExtractor {

    // A UniProtKB FASTA header: db|accession|entryName <rest>, where db is sp (Swiss-Prot) or tr (TrEMBL).
    private static final Pattern UNIPROT_HEADER = Pattern.compile( "^(?:sp|tr)\\|([^|\\s]+)\\|(\\S+)\\s*(.*)$",
                                                                   Pattern.CASE_INSENSITIVE );
    // A UniProt "KEY=" qualifier token (OS=, OX=, GN=, PE=, SV=), at the start or after whitespace.
    private static final Pattern KEY_TOKEN      = Pattern.compile( "(?:^|\\s)([A-Z]{2})=" );
    // A GenBank/RefSeq defline's trailing organism, e.g. "... [Drosophila melanogaster]".
    private static final Pattern TRAILING_TAXON = Pattern.compile( "\\[([^\\]]+)\\]\\s*$" );
    private static final String  PREDICTED      = "PREDICTED: ";
    // An NCBI "CDS from genomic" FASTA defline starts with a local identifier: lcl|<nuc>_prot_<prot>_<n>.
    private static final Pattern LCL_PREFIX     = Pattern.compile( "^lcl\\|", Pattern.CASE_INSENSITIVE );
    // ...followed by [key=value] qualifiers, e.g. [gene=ORF7a] [protein=ORF7a protein] [protein_id=QTJ30158.1].
    private static final Pattern BRACKET_KV     = Pattern.compile( "\\[([a-zA-Z_]+)=([^\\]]*)\\]" );

    private LabelDataExtractor() {
    }

    /**
     * Extracts data from every external node whose name is a UniProt header, writing the parsed fields and
     * shortening the name to its accession. Returns how many nodes were changed. Non-header names are left
     * untouched.
     */
    public static int extract( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return 0;
        }
        int changed = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( !ForesterUtil.isEmpty( node.getName() ) && applyToNode( node, parseHeader( node.getName() ) ) ) {
                ++changed;
            }
        }
        return changed;
    }

    /**
     * True if more than half the external nodes are UniProt headers -- the threshold for <i>proactively</i>
     * offering extraction when a tree is loaded (a few stray header-like names among ordinary labels
     * should not nag the user; the menu action has no such gate). Pure.
     */
    public static boolean mostLabelsParsable( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return false;
        }
        int total = 0;
        int hits = 0;
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            ++total;
            if ( !parseHeader( it.next().getName() ).isEmpty() ) {
                ++hits;
            }
        }
        return ( total > 0 ) && ( ( hits * 2 ) > total );
    }

    /**
     * Parses a FASTA-header label, trying the UniProt dialect first and falling back to the GenBank/RefSeq
     * defline. Pure. {@link SequenceEntry#EMPTY} if the label is neither.
     */
    public static SequenceEntry parseHeader( final String raw ) {
        final SequenceEntry uniprot = parseUniProtHeader( raw );
        if ( !uniprot.isEmpty() ) {
            return uniprot;
        }
        // the NCBI "CDS from genomic" lcl| defline must be tried before the plain GenBank defline: its leading
        // token would otherwise be mis-read as an accession and the [key=value] blob taken for a description.
        final SequenceEntry cds = parseNcbiCdsDefline( raw );
        return cds.isEmpty() ? parseGenbankDefline( raw ) : cds;
    }

    /**
     * Parses an NCBI "CDS from genomic" FASTA defline: a local identifier followed by {@code [key=value]}
     * qualifiers, e.g. {@code lcl|MT726045.1_prot_QTJ30158.1_7 [gene=ORF7a] [protein=ORF7a protein]
     * [protein_id=QTJ30158.1] [location=27187..27543] [gbkey=CDS]}. Yields the accession (the
     * {@code protein_id}, else an accession pulled from the {@code lcl|} token), the description (the
     * {@code protein} qualifier) and the gene. This dialect carries <b>no organism</b> -- the organism lives
     * on the genome record, not the CDS defline -- so none is guessed. Pure; {@link SequenceEntry#EMPTY} for
     * anything that is not an {@code lcl|} defline with {@code [key=value]} qualifiers.
     */
    public static SequenceEntry parseNcbiCdsDefline( final String raw ) {
        if ( ForesterUtil.isEmpty( raw ) ) {
            return SequenceEntry.EMPTY;
        }
        String s = raw.trim();
        if ( s.startsWith( ">" ) ) {
            s = s.substring( 1 ).trim();
        }
        if ( !LCL_PREFIX.matcher( s ).find() ) {
            return SequenceEntry.EMPTY; // not a local (lcl|) defline
        }
        String protein = null;    // -> sequence name (description)
        String gene = null;
        String protein_id = null; // -> accession
        boolean any_qualifier = false;
        final Matcher m = BRACKET_KV.matcher( s );
        while ( m.find() ) {
            any_qualifier = true;
            final String key = m.group( 1 ).toLowerCase( Locale.ROOT );
            final String value = m.group( 2 ).trim();
            if ( "protein".equals( key ) && ( protein == null ) ) {
                protein = value;
            }
            else if ( "gene".equals( key ) && ( gene == null ) ) {
                gene = value;
            }
            else if ( "protein_id".equals( key ) && ( protein_id == null ) ) {
                protein_id = value;
            }
        }
        if ( !any_qualifier ) {
            return SequenceEntry.EMPTY; // "lcl|..." with no [key=value] qualifiers -> not this dialect
        }
        // accession: prefer the protein_id qualifier; otherwise pull an accession out of the leading lcl| token
        String accession = emptyToNull( protein_id );
        if ( accession == null ) {
            final String token = s.split( "\\s+", 2 )[ 0 ];
            accession = SequenceAccessionTools.parseRefSeqAccessorFromString( token );
            if ( ForesterUtil.isEmpty( accession ) ) {
                accession = SequenceAccessionTools.parseGenbankAccessorFromString( token );
            }
            accession = emptyToNull( accession );
        }
        if ( ( accession == null ) && ( protein == null ) && ( gene == null ) ) {
            return SequenceEntry.EMPTY; // nothing usable
        }
        return new SequenceEntry( accession, null, emptyToNull( protein ), emptyToNull( gene ), null,
                                  SequenceEntry.MoleculeType.PROTEIN, null, null, null, null, null, null );
    }

    /**
     * Parses an NCBI GenBank/RefSeq defline ({@code ACCESSION.v description [Organism]}). The leading token
     * must be a RefSeq/GenBank accession (else {@link SequenceEntry#EMPTY}); the organism is taken only
     * from a trailing {@code [..]} bracket -- never guessed from a nucleotide record's prose -- and the
     * "PREDICTED: " RefSeq prefix is dropped from the description. No gene, no tax-id. Pure.
     */
    public static SequenceEntry parseGenbankDefline( final String raw ) {
        if ( ForesterUtil.isEmpty( raw ) ) {
            return SequenceEntry.EMPTY;
        }
        String s = raw.trim();
        if ( s.startsWith( ">" ) ) {
            s = s.substring( 1 ).trim();
        }
        final String[] parts = s.split( "\\s+", 2 );
        String accession = SequenceAccessionTools.parseRefSeqAccessorFromString( parts[ 0 ] );
        if ( ForesterUtil.isEmpty( accession ) ) {
            accession = SequenceAccessionTools.parseGenbankAccessorFromString( parts[ 0 ] );
        }
        if ( ForesterUtil.isEmpty( accession ) ) {
            return SequenceEntry.EMPTY; // the leading token is not an NCBI accession -> not a GenBank defline
        }
        String rest = ( parts.length > 1 ) ? parts[ 1 ].trim() : "";
        String organism = null;
        final Matcher b = TRAILING_TAXON.matcher( rest );
        if ( b.find() ) {
            organism = b.group( 1 ).trim();
            rest = rest.substring( 0, b.start() ).trim();
        }
        if ( rest.startsWith( PREDICTED ) ) {
            rest = rest.substring( PREDICTED.length() ).trim();
        }
        return new SequenceEntry( accession, null, emptyToNull( rest ), null, null, null, emptyToNull( organism ), null,
                                  null, null, null, null );
    }

    /**
     * Parses a UniProt FASTA header into a {@link SequenceEntry} (accession, entry name, protein name,
     * gene, organism scientific name + NCBI tax-id). Pure. Returns {@link SequenceEntry#EMPTY} for
     * anything that is not a {@code sp|}/{@code tr|} header, so non-UniProt labels are simply ignored.
     */
    public static SequenceEntry parseUniProtHeader( final String raw ) {
        if ( ForesterUtil.isEmpty( raw ) ) {
            return SequenceEntry.EMPTY;
        }
        String s = raw.trim();
        if ( s.startsWith( ">" ) ) {
            s = s.substring( 1 ).trim();
        }
        final Matcher m = UNIPROT_HEADER.matcher( s );
        if ( !m.matches() ) {
            return SequenceEntry.EMPTY;
        }
        final String accession = m.group( 1 );
        final String entry_name = m.group( 2 );
        final String rest = m.group( 3 ); // "ProteinName OS=... OX=... GN=... PE=... SV=..."
        final List<int[]> spans = new ArrayList<int[]>(); // {tokenStart, valueStart} for each KEY= token
        final List<String> keys = new ArrayList<String>();
        final Matcher k = KEY_TOKEN.matcher( rest );
        while ( k.find() ) {
            keys.add( k.group( 1 ).toUpperCase( Locale.ROOT ) );
            spans.add( new int[] { k.start(), k.end() } );
        }
        // the protein name is everything before the first qualifier token
        final String protein_name = spans.isEmpty() ? rest.trim() : rest.substring( 0, spans.get( 0 )[ 0 ] ).trim();
        String os = null;
        String ox = null;
        String gn = null;
        for( int i = 0; i < keys.size(); ++i ) {
            final int value_start = spans.get( i )[ 1 ];
            final int value_end = ( ( i + 1 ) < spans.size() ) ? spans.get( i + 1 )[ 0 ] : rest.length();
            final String value = rest.substring( value_start, value_end ).trim();
            if ( "OS".equals( keys.get( i ) ) && ( os == null ) ) {
                os = value;
            }
            else if ( "OX".equals( keys.get( i ) ) && ( ox == null ) ) {
                ox = value;
            }
            else if ( "GN".equals( keys.get( i ) ) && ( gn == null ) ) {
                gn = value;
            }
        }
        final String organism_id = ( ( ox != null ) && ox.matches( "\\d+" ) ) ? ox : null;
        return new SequenceEntry( accession, emptyToNull( entry_name ), emptyToNull( protein_name ), emptyToNull( gn ),
                                  null, SequenceEntry.MoleculeType.PROTEIN, emptyToNull( os ), organism_id, null, null,
                                  null, null );
    }

    /**
     * Writes {@code e} onto {@code node}, filling only empty fields (never clobbering curated data) and
     * shortening the node name to the accession. Returns true if anything changed.
     */
    static boolean applyToNode( final PhylogenyNode node, final SequenceEntry e ) {
        if ( ( node == null ) || ( e == null ) || e.isEmpty() ) {
            return false;
        }
        boolean changed = false;
        // sequence: fill only empty fields
        final boolean had_seq = node.getNodeData().isHasSequence();
        final Sequence seq = had_seq ? node.getNodeData().getSequence() : new Sequence();
        boolean seq_changed = false;
        if ( ( ( seq.getAccession() == null ) || ForesterUtil.isEmpty( seq.getAccession().getValue() ) )
                && !ForesterUtil.isEmpty( e.getAccession() ) ) {
            seq.setAccession( new Accession( e.getAccession(), accessionSource( e.getAccession() ) ) );
            seq_changed = true;
        }
        if ( ForesterUtil.isEmpty( seq.getName() ) && !ForesterUtil.isEmpty( e.getSequenceName() ) ) {
            seq.setName( e.getSequenceName() );
            seq_changed = true;
        }
        if ( ForesterUtil.isEmpty( seq.getGeneName() ) && !ForesterUtil.isEmpty( e.getGeneName() ) ) {
            seq.setGeneName( e.getGeneName() );
            seq_changed = true;
        }
        if ( seq_changed && !had_seq ) {
            node.getNodeData().addSequence( seq );
        }
        changed |= seq_changed;
        // taxonomy: fill only empty fields
        final boolean had_tax = node.getNodeData().isHasTaxonomy();
        final Taxonomy tax = had_tax ? node.getNodeData().getTaxonomy() : new Taxonomy();
        boolean tax_changed = false;
        if ( ForesterUtil.isEmpty( tax.getScientificName() ) && !ForesterUtil.isEmpty( e.getOrganismName() ) ) {
            tax.setScientificName( e.getOrganismName() );
            tax_changed = true;
        }
        if ( ( tax.getIdentifier() == null ) && !ForesterUtil.isEmpty( e.getOrganismId() ) ) {
            tax.setIdentifier( new Identifier( e.getOrganismId(), "ncbi" ) );
            tax_changed = true;
        }
        if ( tax_changed && !had_tax ) {
            node.getNodeData().addTaxonomy( tax );
        }
        changed |= tax_changed;
        // shorten the node name to the accession (the full header now lives in the structured fields)
        if ( !ForesterUtil.isEmpty( e.getAccession() ) && !e.getAccession().equals( node.getName() ) ) {
            node.setName( e.getAccession() );
            changed = true;
        }
        return changed;
    }

    private static String emptyToNull( final String s ) {
        return ForesterUtil.isEmpty( s ) ? null : s;
    }

    /** The accession's database source label ("refseq"/"ncbi"/"uniprot"/...), so a later Fetch queries the right DB. */
    private static String accessionSource( final String accession ) {
        final Accession a = SequenceAccessionTools.parseAccessorFromString_UniProtPriority( accession );
        return ( a != null ) ? a.getSource() : Accession.Source.NCBI.toString();
    }
}
