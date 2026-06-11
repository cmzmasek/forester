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
        RNA, DNA, AA, GENERAL;
    }
}