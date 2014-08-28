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

public interface MolecularSequence {

    public static final char   UNSPECIFIED_AA  = 'X';
    public static final char   UNSPECIFIED_NUC = 'N';
    public static final char   GAP             = '-';
    public static final String GAP_STR         = Character.toString( GAP );
    public static final char   TERMINATE       = '*';
    static final String        AA_REGEXP       = "[^ARNDBCQEZGHILKMFPSTWYVXUO\\-\\*]";
    static final String        DNA_REGEXP      = "[^ACGTRYMKWSN\\-\\*]";
    static final String        RNA_REGEXP      = "[^ACGURYMKWSN\\-\\*]";

    public abstract String getIdentifier();

    public abstract int getLength();

    public abstract int getNumberOfGapResidues();

    public abstract char[] getMolecularSequence();

    public abstract String getMolecularSequenceAsString();

    public abstract char getResidueAt( final int position );

    public abstract boolean isGapAt( final int position );

    public abstract TYPE getType();

    public enum TYPE {
        RNA, DNA, AA;
    }
}