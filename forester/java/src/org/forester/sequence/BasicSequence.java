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
// WWW: www.phylosoft.org/forester

package org.forester.sequence;

public class BasicSequence implements Sequence {

    private final char[] _mol_sequence;
    private final Object _identifier;
    private final TYPE   _type;

    private BasicSequence( final Object identifier, final String mol_sequence, final TYPE type ) {
        _mol_sequence = mol_sequence.toCharArray();
        _identifier = identifier;
        _type = type;
    }

    // Only use if you know what you are doing!
    public BasicSequence( final Object identifier, final char[] mol_sequence, final TYPE type ) {
        _mol_sequence = mol_sequence;
        _identifier = identifier;
        _type = type;
    }

    public Object getIdentifier() {
        return _identifier;
    }

    public int getLength() {
        return _mol_sequence.length;
    }

    public char[] getMolecularSequence() {
        return _mol_sequence;
    }

    public char getResidueAt( final int position ) {
        return _mol_sequence[ position ];
    }

    public TYPE getType() {
        return _type;
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( _identifier.toString() );
        sb.append( " " );
        sb.append( new String( _mol_sequence ) );
        return sb.toString();
    }

    public static Sequence createAaSequence( final Object identifier, final String mol_sequence ) {
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                .replaceAll( AA_REGEXP, Character.toString( UNSPECIFIED_AA ) ), TYPE.AA );
    }

    public static Sequence createDnaSequence( final Object identifier, final String mol_sequence ) {
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                .replaceAll( DNA_REGEXP, Character.toString( UNSPECIFIED_NUC ) ), TYPE.DNA );
    }

    public static Sequence createRnaSequence( final Object identifier, final String mol_sequence ) {
        return new BasicSequence( identifier, mol_sequence.toUpperCase().replaceAll( "\\.", GAP_STR )
                .replaceAll( RNA_REGEXP, Character.toString( UNSPECIFIED_NUC ) ), TYPE.RNA );
    }
}
