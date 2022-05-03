
package org.forester.application;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;

public class asr_prepare {

    private final static String  PRG_NAME    = "asr_prepare";
    private static final String  PRG_DATE    = "2021-10-11";
    private static final String  PRG_VERSION = "1.0.0";
    private static final Pattern P1          = Pattern.compile( "\\d+" );
    public static void main( final String[] args ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        if ( args.length != 4 ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <in-msa> <out-tree> <out-msa>\n" );
            System.out.println( "Example: " + PRG_NAME + " tree.nh msa.fasta ct.nh cmsa.fasta\n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File inmsa = new File( args[ 1 ] );
        final File outtree = new File( args[ 2 ] );
        final File outmsa = new File( args[ 3 ] );
        final String error0 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( error0 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error0 );
        }
        final String error1 = ForesterUtil.isReadableFile( inmsa );
        if ( !ForesterUtil.isEmpty( error1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error1 );
        }
        final String error2 = ForesterUtil.isWritableFile( outtree );
        if ( !ForesterUtil.isEmpty( error2 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error2 );
        }
        final String error3 = ForesterUtil.isWritableFile( outmsa );
        if ( !ForesterUtil.isEmpty( error3 ) ) {
            ForesterUtil.fatalError( PRG_NAME, error3 );
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "failure to read input tree [" + e.getMessage() + "]" );
        }
        List<MolecularSequence> seqs = null;
        try {
            seqs = FastaParser.parse( new FileInputStream( inmsa ) );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failure to read input MSA [" + e.getMessage() + "]" );
        }
        final SortedSet<String> msa_ids = processSequences( seqs );
        final SortedSet<String> tree_ids = processTree( p );
        System.out.println( "MSA              : " + msa_ids.size() + " identifiers" );
        System.out.println( "Tree             : " + tree_ids.size() + " identifiers" );
        final SortedSet<String> msa_ids_orig = new TreeSet<>();
        msa_ids_orig.addAll( msa_ids );
        msa_ids.retainAll( tree_ids );
        System.out.println( "Intersection     : " + msa_ids.size() + " identifiers" );
        tree_ids.removeAll( msa_ids );
        msa_ids_orig.removeAll( msa_ids );
        System.out.println( "Removed from tree: " + tree_ids.size() + " identifiers" );
        System.out.println( "Removed from MSA : " + msa_ids_orig.size() + " identifiers" );
        System.out.println();
        System.out.println( "Removed from tree:" );
        for( final String id : tree_ids ) {
            System.out.println( id );
        }
        System.out.println();
        System.out.println( "Removed from MSA: " );
        for( final String id : msa_ids_orig ) {
            System.out.println( id );
        }
        final List<MolecularSequence> clean_seqs = new ArrayList<>();
        for( final MolecularSequence seq : seqs ) {
            if ( msa_ids.contains( seq.getIdentifier() ) ) {
                clean_seqs.add( seq );
            }
        }
        try {
            SequenceWriter.writeSeqs( clean_seqs, outmsa, SEQ_FORMAT.FASTA, 80 );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failure to write output MSA [" + e.getMessage() + "]" );
        }
        System.out.println();
        System.out.println( "Wrote " + clean_seqs.size() + " sequences to: " + outmsa );
        final List<PhylogenyNode> to_delete = new ArrayList<>();
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                if ( !msa_ids.contains( node.getName() ) ) {
                    to_delete.add( node );
                }
            }
        }
        for( final PhylogenyNode node : to_delete ) {
            p.deleteSubtree( node, true );
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toNewHampshire( p, true, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE, outtree );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failure to write output tree [" + e.getMessage() + "]" );
        }
        System.out.println( "Wrote tree to: " + outtree );
        System.out.println( "OK." );
    }

    private static SortedSet<String> processTree( final Phylogeny p ) {
        final SortedSet<String> tree_ids = new TreeSet<>();
        for( final PhylogenyNodeIterator iter = p.iteratorPreorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                final Matcher m1 = P1.matcher( node.getName() );
                if ( m1.matches() ) {
                    node.setName( "_" + node.getName() );
                }
                if ( ForesterUtil.isEmpty( node.getName() ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "tree node with empty identifier detected" );
                }
                else if ( tree_ids.contains( node.getName() ) ) {
                    ForesterUtil.fatalError( PRG_NAME, "non unique identifier in tree: " + node.getName() );
                }
                else {
                    tree_ids.add( node.getName() );
                }
            }
        }
        return tree_ids;
    }

    private static SortedSet<String> processSequences( final List<MolecularSequence> seqs ) {
        final SortedSet<String> msa_ids = new TreeSet<>();
        int l = -1;
        for( final MolecularSequence seq : seqs ) {
            if ( l == -1 ) {
                l = seq.getLength();
            }
            else if ( l != seq.getLength() ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "sequences of unequal length detected: " + l + " != " + seq.getLength() + ", "
                                                 + seq.getIdentifier() );
            }
            final Matcher m1 = P1.matcher( seq.getIdentifier() );
            if ( m1.matches() ) {
                ( ( BasicSequence ) seq ).setIdentifier( "_" + seq.getIdentifier() );
            }
            if ( ForesterUtil.isEmpty( seq.getIdentifier() ) ) {
                ForesterUtil.fatalError( PRG_NAME, "sequence with empty identifier in MSA" );
            }
            else if ( msa_ids.contains( seq.getIdentifier() ) ) {
                ForesterUtil.fatalError( PRG_NAME, "non unique sequence identifier in MSA: " + seq.getIdentifier() );
            }
            else {
                msa_ids.add( seq.getIdentifier() );
            }
        }
        return msa_ids;
    }
}
