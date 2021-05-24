
package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class protein_seq_x {

    private final static String PRG_NAME    = "protein_seq_x";
    private static final String PRG_DATE    = "2021-05-24";
    private static final String PRG_VERSION = "1.0.0";
    private final static int    MAX_LENGTH  = 1273;
    public static void main( final String[] args ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( args.length != 3 ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <tree file> <aa fasta infile> <aa fasta outfile>\n" );
            System.out.println( "Example: " + PRG_NAME + " tree.xml all_s_proteins.fasta s_proteins.fasta\n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File seqs_infile = new File( args[ 1 ] );
        final File seqs_outfile = new File( args[ 2 ] );
        final String error0 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( error0 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error0 );
        }
        final String error1 = ForesterUtil.isReadableFile( seqs_infile );
        if ( !ForesterUtil.isEmpty( error1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error1 );
        }
        final String error2 = ForesterUtil.isWritableFile( seqs_outfile );
        if ( !ForesterUtil.isEmpty( error2 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error2 );
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, " Error reading from " + intree + ": " + e.getLocalizedMessage() );
        }
        List<MolecularSequence> seqs = null;
        try {
            seqs = FastaParser.parse( new FileInputStream( seqs_infile ) );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getLocalizedMessage() );
        }
        final List<MolecularSequence> out_seqs = new ArrayList<>();
        final BasicDescriptiveStatistics stats = new BasicDescriptiveStatistics();
        final SortedSet<String> accs = new TreeSet<>();
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                if ( node.isHasNodeData() && node.getNodeData().isHasSequence()
                        && ( node.getNodeData().getSequence().getAccession() != null )
                        && !ForesterUtil.isEmpty( node.getNodeData().getSequence().getAccession().getValue() ) ) {
                    final String acc = node.getNodeData().getSequence().getAccession().getValue();
                    if ( accs.contains( acc ) ) {
                        ForesterUtil.fatalError( PRG_NAME, "duplicate sequence accession found: " + acc );
                    }
                    else {
                        accs.add( acc );
                        boolean found = false;
                        for( final MolecularSequence seq : seqs ) {
                            final String name = seq.getIdentifier();
                            final String[] s = name.split( "\\|" );
                            if ( s.length < 3 ) {
                                ForesterUtil.fatalError( PRG_NAME, "illegal sequence name: " + name );
                            }
                            final String genome_acc = s[ 1 ];
                            if ( acc.equals( genome_acc ) ) {
                                out_seqs.add( BasicSequence.createAaSequence( genome_acc,
                                                                              seq.getMolecularSequenceAsString() ) );
                                if ( seq.getLength() > MAX_LENGTH ) {
                                    System.out.println( genome_acc + ": " + seq.getLength() );
                                }
                                stats.addValue( seq.getLength() );
                                found = true;
                                continue;
                            }
                        }
                        if ( !found ) {
                            ForesterUtil.fatalError( PRG_NAME,
                                                     "no sequence found for node: " + node.getName() + ": " + acc );
                        }
                    }
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "external node without sequence accession found" );
                }
            }
        }
        try {
            SequenceWriter.writeSeqs( out_seqs, seqs_outfile, SEQ_FORMAT.FASTA, 80 );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "Error writing to " + seqs_outfile + ": " + e.getLocalizedMessage() );
        }
        System.out.println( "Wrote sequences to     : " + seqs_outfile );
        System.out.println( "Number of sequences    : " + out_seqs.size() );
        System.out.println( "Minimum sequence length: " + stats.getMin() );
        System.out.println( "Maximum sequence length: " + stats.getMax() );
        System.out.println( "Mean sequence length   : " + stats.arithmeticMean() );
        System.out.println( "Median sequence length : " + stats.median() + "\n" );
    }
}
