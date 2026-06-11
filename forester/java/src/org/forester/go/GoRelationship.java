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

package org.forester.go;

public interface GoRelationship extends Comparable<GoRelationship> {

    public static final String PART_OF_STR              = "part_of";
    public static final String REGULATES_STR            = "regulates";
    public static final String NEGATIVELY_REGULATES_STR = "negatively_regulates";
    public static final String POSITIVELY_REGULATES_STR = "positively_regulates";
    public static final String HAS_PART_STR             = "has_part";
    public static final String OCCURS_IN_STR            = "occurs_in";
    public static final String HAPPENS_DURING_STR       = "happens_during";
    public static final String ENDS_DURING_STR       = "ends_during";

    public GoId getGoId();

    public Type getType();

    public static enum Type {
        PART_OF, REGULATES, NEGATIVELY_REGULATES, POSITIVELY_REGULATES, HAS_PART, OCCURS_IN, HAPPENS_DURING, ENDS_DURING;
    }
}
