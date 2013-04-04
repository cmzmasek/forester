
package org.forester.io.writers;

import java.io.IOException;
import java.io.Writer;
import java.util.List;

import org.forester.sequence.Sequence;
import org.forester.util.ForesterUtil;

public class SequenceWriter {

    public static enum SEQ_FORMAT {
        FASTA;
    }

    public static StringBuilder toFasta( final Sequence seq, final int width ) {
        return toFasta( seq.getIdentifier(), seq.getMolecularSequenceAsString(), width );
    }

    public static StringBuilder toFasta( final String name, final String mol_seq, final int width ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( ">" );
        sb.append( name );
        sb.append( ForesterUtil.LINE_SEPARATOR );
        if ( ( width < 1 ) || ( width >= mol_seq.length() ) ) {
            sb.append( mol_seq );
        }
        else {
            final int lines = mol_seq.length() / width;
            final int rest = mol_seq.length() - ( lines * width );
            for( int i = 0; i < lines; ++i ) {
                sb.append( mol_seq, i * width, ( i + 1 ) * width );
                if ( i < ( lines - 1 ) ) {
                    sb.append( ForesterUtil.LINE_SEPARATOR );
                }
            }
            if ( rest > 0 ) {
                sb.append( ForesterUtil.LINE_SEPARATOR );
                sb.append( mol_seq, lines * width, mol_seq.length() );
            }
        }
        return sb;
    }

    public static void toFasta( final Sequence seq, final Writer w, final int width ) throws IOException {
        w.write( ">" );
        w.write( seq.getIdentifier().toString() );
        w.write( ForesterUtil.LINE_SEPARATOR );
        if ( ( width < 1 ) || ( width >= seq.getLength() ) ) {
            w.write( seq.getMolecularSequence() );
        }
        else {
            final int lines = seq.getLength() / width;
            final int rest = seq.getLength() - ( lines * width );
            for( int i = 0; i < lines; ++i ) {
                w.write( seq.getMolecularSequence(), i * width, width );
                if ( i < ( lines - 1 ) ) {
                    w.write( ForesterUtil.LINE_SEPARATOR );
                }
            }
            if ( rest > 0 ) {
                w.write( ForesterUtil.LINE_SEPARATOR );
                w.write( seq.getMolecularSequence(), lines * width, rest );
            }
        }
    }

    public static void writeSeqs( final List<Sequence> seqs,
                                  final Writer writer,
                                  final SEQ_FORMAT format,
                                  final int width ) throws IOException {
        switch ( format ) {
            case FASTA:
                for( final Sequence s : seqs ) {
                    toFasta( s, writer, width );
                    writer.write( ForesterUtil.LINE_SEPARATOR );
                }
                break;
            default:
                throw new RuntimeException( "unknown format " + format );
        }
    }
}
