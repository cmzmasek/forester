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

package org.forester.io.writers;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;

public class SequenceWriter {

    public static enum SEQ_FORMAT {
        FASTA;
    }

    public static StringBuilder toFasta( final MolecularSequence seq, final int width ) {
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

    public static void toFasta( final MolecularSequence seq, final Writer w, final int width ) throws IOException {
        w.write( ">" );
        w.write( seq.getIdentifier() );
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

    public static void writeSeqs( final List<MolecularSequence> seqs,
                                  final File file,
                                  final SEQ_FORMAT format,
                                  final int width ) throws IOException {
        final Writer w = ForesterUtil.createBufferedWriter( file );
        SequenceWriter.writeSeqs( seqs, w, format, width );
        w.close();
    }

    public static void writeSeqs( final List<MolecularSequence> seqs,
                                  final Writer writer,
                                  final SEQ_FORMAT format,
                                  final int width ) throws IOException {
        switch ( format ) {
            case FASTA:
                for( final MolecularSequence s : seqs ) {
                    toFasta( s, writer, width );
                    writer.write( ForesterUtil.LINE_SEPARATOR );
                }
                break;
            default:
                throw new RuntimeException( "unknown format " + format );
        }
    }
}
