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
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Function;

import org.forester.io.writers.SequenceWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.PhylogenyDataUtil;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.ForesterUtil;

/**
 * Read-only export of per-tip data: the molecular sequences as FASTA, and the taxonomy / sequence /
 * branch / property fields as a tab-separated table (one row per external node). Pure (no GUI).
 *
 * <p>The TSV columns are dynamic -- a column appears only when at least one tip has a value for it (the
 * {@code name} key column is always present) -- and the {@code name} column is the key the annotation
 * IMPORT will match tips on, so this schema is the round-trip contract.
 */
public final class NodeDataExporter {

    private static final int    FASTA_LINE_WIDTH = 60;
    private static final String TAB              = "\t";

    /**
     * FASTA of every tip's molecular sequence(s). The header is a self-describing line -- the tip identifier
     * (the round-trip key, kept as the first whitespace-delimited token), then, when present, the sequence
     * accession, the organism (scientific or common name), and {@code " | "} + a sequence/protein
     * description. Empty if no tip carries a molecular sequence.
     */
    public static String toFasta( final Phylogeny phy ) {
        final StringBuilder sb = new StringBuilder();
        if ( phy == null ) {
            return sb.toString();
        }
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( !n.getNodeData().isHasSequence() ) {
                continue;
            }
            for( final Sequence seq : n.getNodeData().getSequences() ) {
                if ( ( seq == null ) || ForesterUtil.isEmpty( seq.getMolecularSequence() ) ) {
                    continue;
                }
                sb.append( SequenceWriter.toFasta( fastaHeader( n, seq ), seq.getMolecularSequence(),
                                                   FASTA_LINE_WIDTH ) );
                sb.append( ForesterUtil.LINE_SEPARATOR );
            }
        }
        return sb.toString();
    }

    /** The FASTA header id: the tip name, else the sequence accession, else the sequence name. */
    private static String fastaIdentifier( final PhylogenyNode n, final Sequence seq ) {
        if ( !ForesterUtil.isEmpty( n.getName() ) ) {
            return n.getName();
        }
        if ( ( seq.getAccession() != null ) && !ForesterUtil.isEmpty( seq.getAccession().getValue() ) ) {
            return seq.getAccession().getValue();
        }
        if ( !ForesterUtil.isEmpty( seq.getName() ) ) {
            return seq.getName();
        }
        return "";
    }

    /**
     * A self-describing FASTA header line: {@code id [accession] [organism] [| description]}. The id is the
     * first token (so downstream FASTA tools still key on it); accession, organism and description are added
     * only when present, and the accession is skipped when it is already the id.
     */
    private static String fastaHeader( final PhylogenyNode n, final Sequence seq ) {
        final List<String> prefix = new ArrayList<>();
        final String id = fastaIdentifier( n, seq );
        if ( !id.isEmpty() ) {
            prefix.add( id );
        }
        final String acc = ( seq.getAccession() != null ) ? seq.getAccession().getValue() : "";
        if ( !ForesterUtil.isEmpty( acc ) && !acc.equals( id ) ) {
            prefix.add( acc );
        }
        final String organism = organismName( n );
        if ( !organism.isEmpty() ) {
            prefix.add( organism );
        }
        String header = String.join( " ", prefix );
        final String description = sequenceDescription( seq );
        if ( !description.isEmpty() ) {
            header = header.isEmpty() ? description : header + " | " + description;
        }
        return sanitize( header );
    }

    /** Organism for the FASTA header: scientific name, else common name. */
    private static String organismName( final PhylogenyNode n ) {
        if ( n.getNodeData().isHasTaxonomy() ) {
            final Taxonomy t = n.getNodeData().getTaxonomy();
            if ( !ForesterUtil.isEmpty( t.getScientificName() ) ) {
                return t.getScientificName();
            }
            if ( !ForesterUtil.isEmpty( t.getCommonName() ) ) {
                return t.getCommonName();
            }
        }
        return "";
    }

    /** Human-readable sequence descriptor for the FASTA header: sequence name, else gene name, else symbol. */
    private static String sequenceDescription( final Sequence seq ) {
        if ( !ForesterUtil.isEmpty( seq.getName() ) ) {
            return seq.getName();
        }
        if ( !ForesterUtil.isEmpty( seq.getGeneName() ) ) {
            return seq.getGeneName();
        }
        if ( !ForesterUtil.isEmpty( seq.getSymbol() ) ) {
            return seq.getSymbol();
        }
        return "";
    }

    /**
     * Whether the external-node names can serve as the import key: every tip has a non-empty name and all the
     * names are distinct. When false, {@link #toNodeDataTsv} adds a {@code node_id} key column.
     */
    public static boolean tipNamesFormUniqueKey( final Phylogeny phy ) {
        return ( phy != null ) && namesFormUniqueKey( phy.getExternalNodes() );
    }

    private static boolean namesFormUniqueKey( final List<PhylogenyNode> tips ) {
        final Set<String> seen = new HashSet<>();
        for( final PhylogenyNode n : tips ) {
            if ( ForesterUtil.isEmpty( n.getName() ) || !seen.add( n.getName() ) ) {
                return false;
            }
        }
        return true;
    }

    /** Tab-separated tip-data table; columns with no value across every tip are omitted (except {@code name}). */
    public static String toNodeDataTsv( final Phylogeny phy ) {
        if ( phy == null ) {
            return "";
        }
        final List<PhylogenyNode> tips = phy.getExternalNodes();
        final LinkedHashMap<String, String[]> cols = new LinkedHashMap<>();
        // When tip names are blank or duplicated they cannot key the table back to tips on import, so prepend
        // a guaranteed-unique node_id column (the stable per-tip node id within this tree) in that case.
        if ( !namesFormUniqueKey( tips ) ) {
            addColumn( cols, "node_id", tips, n -> String.valueOf( n.getId() ), true );
        }
        addColumn( cols, "name", tips, n -> n.getName(), true );
        addColumn( cols, "taxonomy_scientific_name", tips, n -> tax( n, Taxonomy::getScientificName ), false );
        addColumn( cols, "taxonomy_common_name", tips, n -> tax( n, Taxonomy::getCommonName ), false );
        addColumn( cols, "taxonomy_code", tips, n -> tax( n, Taxonomy::getTaxonomyCode ), false );
        addColumn( cols, "taxonomy_id", tips, NodeDataExporter::taxId, false );
        addColumn( cols, "taxonomy_rank", tips, n -> tax( n, Taxonomy::getRank ), false );
        addColumn( cols, "sequence_name", tips, n -> seq( n, Sequence::getName ), false );
        addColumn( cols, "gene_name", tips, n -> seq( n, Sequence::getGeneName ), false );
        addColumn( cols, "sequence_symbol", tips, n -> seq( n, Sequence::getSymbol ), false );
        addColumn( cols, "sequence_accession", tips, NodeDataExporter::seqAccession, false );
        addColumn( cols, "sequence_type", tips, n -> seq( n, Sequence::getType ), false );
        addColumn( cols, "branch_length", tips, NodeDataExporter::branchLength, false );
        for( final String ref : propertyRefs( tips ) ) {
            addColumn( cols, ref, tips, n -> propertyValue( n, ref ), false );
        }
        final StringBuilder sb = new StringBuilder();
        sb.append( String.join( TAB, cols.keySet() ) );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        for( int i = 0; i < tips.size(); i++ ) {
            final List<String> row = new ArrayList<>( cols.size() );
            for( final String[] vals : cols.values() ) {
                row.add( vals[ i ] );
            }
            sb.append( String.join( TAB, row ) );
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        return sb.toString();
    }

    private static void addColumn( final Map<String, String[]> cols,
                                   final String name,
                                   final List<PhylogenyNode> tips,
                                   final Function<PhylogenyNode, String> extractor,
                                   final boolean force ) {
        final String[] vals = new String[ tips.size() ];
        boolean any = false;
        for( int i = 0; i < tips.size(); i++ ) {
            final String v = sanitize( extractor.apply( tips.get( i ) ) );
            vals[ i ] = v;
            if ( !v.isEmpty() ) {
                any = true;
            }
        }
        if ( any || force ) {
            cols.put( name, vals );
        }
    }

    private static String tax( final PhylogenyNode n, final Function<Taxonomy, String> f ) {
        return n.getNodeData().isHasTaxonomy() ? f.apply( n.getNodeData().getTaxonomy() ) : "";
    }

    private static String taxId( final PhylogenyNode n ) {
        if ( n.getNodeData().isHasTaxonomy() && ( n.getNodeData().getTaxonomy().getIdentifier() != null ) ) {
            return n.getNodeData().getTaxonomy().getIdentifier().getValue();
        }
        return "";
    }

    private static String seq( final PhylogenyNode n, final Function<Sequence, String> f ) {
        return n.getNodeData().isHasSequence() ? f.apply( n.getNodeData().getSequence() ) : "";
    }

    private static String seqAccession( final PhylogenyNode n ) {
        if ( n.getNodeData().isHasSequence() ) {
            final Accession acc = n.getNodeData().getSequence().getAccession();
            if ( acc != null ) {
                return acc.getValue();
            }
        }
        return "";
    }

    private static String branchLength( final PhylogenyNode n ) {
        final double d = n.getDistanceToParent();
        return ( d == PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT ) ? "" : String.valueOf( d );
    }

    private static SortedSet<String> propertyRefs( final List<PhylogenyNode> tips ) {
        final SortedSet<String> refs = new TreeSet<>();
        for( final PhylogenyNode n : tips ) {
            if ( n.getNodeData().getProperties() != null ) {
                for( final Property p : n.getNodeData().getProperties().getProperties() ) {
                    if ( !ForesterUtil.isEmpty( p.getRef() ) ) {
                        refs.add( p.getRef() );
                    }
                }
            }
        }
        return refs;
    }

    private static String propertyValue( final PhylogenyNode n, final String ref ) {
        if ( n.getNodeData().getProperties() == null ) {
            return "";
        }
        final List<Property> ps = n.getNodeData().getProperties().getPropertiesWithGivenRef( ref );
        return ( ( ps != null ) && !ps.isEmpty() ) ? ps.get( 0 ).getValue() : "";
    }

    /** Tabs / newlines in a value would break the TSV; collapse them to spaces. */
    private static String sanitize( final String s ) {
        if ( s == null ) {
            return "";
        }
        return s.replace( '\t', ' ' ).replace( '\n', ' ' ).replace( '\r', ' ' );
    }

    private NodeDataExporter() {
    }
}
