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

package org.forester.msa;

import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.forester.io.writers.SequenceWriter;
import org.forester.io.writers.SequenceWriter.SEQ_FORMAT;
import org.forester.sequence.BasicSequence;
import org.forester.sequence.MolecularSequence;
import org.forester.sequence.MolecularSequence.TYPE;
import org.forester.util.ForesterUtil;

public class BasicMsa implements Msa {

    private final char[][]    _data;
    private final String[]    _identifiers;
    private final Set<String> _identifiers_set;
    private final TYPE        _type;

    public BasicMsa( final int rows, final int columns, final TYPE type ) {
        if ( ( rows < 1 ) || ( columns < 1 ) ) {
            throw new IllegalArgumentException( "basic msa of size zero are illegal" );
        }
        _data = new char[ rows ][ columns ];
        _identifiers = new String[ rows ];
        _identifiers_set = new HashSet<String>();
        _type = type;
    }

    BasicMsa( final BasicMsa msa ) {
        _data = msa._data;
        _identifiers = msa._identifiers;
        _type = msa._type;
        _identifiers_set = msa._identifiers_set;
    }

    @Override
    public List<MolecularSequence> asSequenceList() {
        final List<MolecularSequence> seqs = new ArrayList<MolecularSequence>();
        for( int i = 0; i < getNumberOfSequences(); ++i ) {
            seqs.add( getSequence( i ) );
        }
        return seqs;
    }

    @Override
    public List<Character> getColumnAt( final int col ) {
        final List<Character> column = new ArrayList<Character>();
        for( int row = 0; row < getNumberOfSequences(); ++row ) {
            column.add( getResidueAt( row, col ) );
        }
        return column;
    }

    @Override
    public String getIdentifier( final int row ) {
        return _identifiers[ row ];
    }

    @Override
    public int getLength() {
        return _data[ 0 ].length;
    }

    @Override
    public int getNumberOfSequences() {
        return _identifiers.length;
    }

    @Override
    public char getResidueAt( final int row, final int col ) {
        return _data[ row ][ col ];
    }

    @Override
    public MolecularSequence getSequence( final int row ) {
        return new BasicSequence( getIdentifier( row ), _data[ row ], getType() );
    }

    @Override
    public MolecularSequence getSequence( final String id ) {
        for( int i = 0; i < getNumberOfSequences(); ++i ) {
            if ( getIdentifier( i ).equals( id ) ) {
                return getSequence( i );
            }
        }
        return null;
    }

    @Override
    public StringBuffer getSequenceAsString( final int row ) {
        final StringBuffer sb = new StringBuffer( getLength() );
        for( int col = 0; col < getLength(); ++col ) {
            sb.append( getResidueAt( row, col ) );
        }
        return sb;
    }

    @Override
    public TYPE getType() {
        return _type;
    }

    @Override
    public boolean isGapAt( final int row, final int col ) {
        return getResidueAt( row, col ) == MolecularSequence.GAP;
    }

    @Override
    public void setIdentifier( final int row, final String id ) {
        if ( ForesterUtil.isEmpty( id ) ) {
            throw new IllegalArgumentException( "illegal attempt to create msa with empty identifier" );
        }
        if ( _identifiers_set.contains( id ) ) {
            throw new IllegalArgumentException( "illegal attempt to create msa with non-unique identifiers [" + id
                    + "]" );
        }
        _identifiers_set.add( id );
        _identifiers[ row ] = id;
    }

    @Override
    public void setResidueAt( final int row, final int col, final char residue ) {
        _data[ row ][ col ] = residue;
    }

    @Override
    public String toString() {
        final Writer w = new StringWriter();
        try {
            write( w, MSA_FORMAT.PHYLIP );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
        return w.toString();
    }

    @Override
    public void write( final Writer w, final MSA_FORMAT format ) throws IOException {
        switch ( format ) {
            case PHYLIP:
                writeToPhylip( w );
                break;
            case FASTA:
                writeToFasta( w );
                break;
            case NEXUS:
                writeToNexus( w );
                break;
            default:
                throw new RuntimeException( "unknown format " + format );
        }
    }

    private short determineMaxIdLength() {
        short max = 0;
        for( int row = 0; row < getNumberOfSequences(); ++row ) {
            final short l = ( short ) getIdentifier( row ).length();
            if ( l > max ) {
                max = l;
            }
        }
        return max;
    }

    private void writeToFasta( final Writer w ) throws IOException {
        SequenceWriter.writeSeqs( asSequenceList(), w, SEQ_FORMAT.FASTA, 100 );
    }

    private void writeToNexus( final Writer w ) throws IOException {
        final int max = determineMaxIdLength() + 1;
        w.write( "Begin Data;" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "   Dimensions NTax=" + getNumberOfSequences() );
        w.write( " NChar=" + getLength() );
        w.write( ";" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "   Format DataType=Protein Interleave=No gap=-;" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "   Matrix" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        for( int row = 0; row < getNumberOfSequences(); ++row ) {
            final MolecularSequence seq = getSequence( row );
            final String s = seq.getMolecularSequenceAsString();
            w.write( "      " );
            w.write( ForesterUtil.pad( getIdentifier( row ).replace( ' ', '_' ), max, ' ', false ).toString() );
            w.write( " " );
            w.write( s );
            w.write( ForesterUtil.LINE_SEPARATOR );
        }
        w.write( "   ;" );
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( "End;" );
        w.write( ForesterUtil.LINE_SEPARATOR );
    }

    private void writeToPhylip( final Writer w ) throws IOException {
        final int max = determineMaxIdLength() + 1;
        w.write( getNumberOfSequences() + " " + getLength() );
        w.write( ForesterUtil.LINE_SEPARATOR );
        for( int row = 0; row < getNumberOfSequences(); ++row ) {
            w.write( ForesterUtil.pad( getIdentifier( row ).replace( ' ', '_' ), max, ' ', false ).toString() );
            for( int col = 0; col < getLength(); ++col ) {
                w.write( getResidueAt( row, col ) );
            }
            w.write( ForesterUtil.LINE_SEPARATOR );
        }
    }

    public static Msa createInstance( final List<MolecularSequence> seqs ) {
        if ( seqs.size() < 1 ) {
            throw new IllegalArgumentException( "cannot create msa from less than one sequence" );
        }
        final int length = seqs.get( 0 ).getLength();
        final BasicMsa msa = new BasicMsa( seqs.size(), length, seqs.get( 0 ).getType() );
        for( int row = 0; row < seqs.size(); ++row ) {
            final MolecularSequence seq = seqs.get( row );
            if ( seq.getLength() != length ) {
                throw new IllegalArgumentException( "illegal attempt to build msa from sequences of unequal length ["
                        + seq.getIdentifier() + "]" );
            }
            if ( seq.getType() != msa.getType() ) {
                throw new IllegalArgumentException( "illegal attempt to build msa from sequences of different type ["
                        + seq.getIdentifier() + "]" );
            }
            msa.setIdentifier( row, seq.getIdentifier() );
            for( int col = 0; col < length; ++col ) {
                msa._data[ row ][ col ] = seq.getResidueAt( col );
            }
        }
        return msa;
    }
}
