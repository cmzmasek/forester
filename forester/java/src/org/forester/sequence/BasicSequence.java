// $Id:
//
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

package org.forester.sequence;

import org.forester.util.ForesterUtil;

public class BasicSequence implements Sequence {

    private final char[] _mol_sequence;
    private String       _identifier;
    private final TYPE   _type;

    private BasicSequence( final String identifier, final String mol_sequence, final TYPE type ) {
        if ( ForesterUtil.isEmpty( identifier ) ) {
            throw new IllegalArgumentException( "identifier of sequence cannot be empty" );
        }
        if ( ForesterUtil.isEmpty( mol_sequence ) ) {
            throw new IllegalArgumentException( "molecular sequence cannot be empty" );
        }
        _mol_sequence = mol_sequence.toCharArray();
        _identifier = identifier;
        _type = type;
    }

    // Only use if you know what you are doing!
    public BasicSequence( final String identifier, final char[] mol_sequence, final TYPE type ) {
        if ( ForesterUtil.isEmpty( identifier ) ) {
            throw new IllegalArgumentException( "identifier of sequence cannot be empty" );
        }
        if ( ( mol_sequence == null ) || ( mol_sequence.length < 1 ) ) {
            throw new IllegalArgumentException( "molecular sequence cannot be empty" );
        }
        _mol_sequence = mol_sequence;
        _identifier = identifier;
        _type = type;
    }

    public void setIdentifier( final String id ) {
        _identifier = id;
    }

    @Override
    public String getIdentifier() {
        return _identifier;
    }

    @Override
    public int getLength() {
        return _mol_sequence.length;
    }

    @Override
    public char[] getMolecularSequence() {
        return _mol_sequence;
    }

    @Override
    public char getResidueAt( final int position ) {
        return _mol_sequence[ position ];
    }

    @Override
    public TYPE getType() {
        return _type;
    }

    @Override
    public int getNumberOfGapResidues() {
        int gaps = 0;
        for( final char element : _mol_sequence ) {
            if ( element == GAP ) {
                ++gaps;
            }
        }
        return gaps;
    }

    @Override
    public boolean equals( final Object obj ) {
        if ( obj == null ) {
            return false;
        }
        if ( obj.getClass() != getClass() ) {
            return false;
        }
        final Sequence other = ( Sequence ) obj;
        if ( getMolecularSequenceAsString().equals( other.getMolecularSequenceAsString() ) ) {
            return true;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return getMolecularSequenceAsString().hashCode();
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( _identifier.toString() );
        sb.append( ": " );
        sb.append( getMolecularSequenceAsString() );
        return sb.toString();
    }

    public static Sequence copySequence( final Sequence seq ) {
        final char[] s = new char[ seq.getMolecularSequence().length ];
        for( int i = 0; i < seq.getMolecularSequence().length; i++ ) {
            s[ i ] = seq.getMolecularSequence()[ i ];
        }
        return new BasicSequence( new String( seq.getIdentifier() ), s, seq.getType() );
    }

    public static Sequence createAaSequence( final String identifier, final String mol_sequence ) {
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                .replaceAll( AA_REGEXP, Character.toString( UNSPECIFIED_AA ) ), TYPE.AA );
    }

    public static Sequence createDnaSequence( final String identifier, final String mol_sequence ) {
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                .replaceAll( DNA_REGEXP, Character.toString( UNSPECIFIED_NUC ) ), TYPE.DNA );
    }

    public static Sequence createRnaSequence( final String identifier, final String mol_sequence ) {
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                .replaceAll( RNA_REGEXP, Character.toString( UNSPECIFIED_NUC ) ), TYPE.RNA );
    }

    @Override
    public String getMolecularSequenceAsString() {
        return new String( getMolecularSequence() );
    }

    @Override
    public boolean isGapAt( final int position ) {
        return getResidueAt( position ) == GAP;
    }
}
