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

public interface GoSubset extends Comparable<GoSubset> {

    public static final String GOSLIM_GENERIC_STR     = "goslim_generic";
    public static final String GOSLIM_GOA_STR         = "goslim_goa";
    public static final String GOSLIM_PIR_STR         = "goslim_pir";
    public static final String GOSUBSET_PROK_STR      = "gosubset_prok";
    public static final String GOSLIM_CANDIDA_STR     = "goslim_candida";
    public static final String GOSLIM_ASPERGILLUS_STR = "goslim_aspergillus";
    public static final String GOSLIM_PLANT_STR       = "goslim_plant";
    public static final String GOSLIM_YEAST_STR       = "goslim_yeast";
    public static final String GOSLIM_POMBE_STR       = "goslim_pombe";

    public Type getType();

    public static enum Type {
        GOSLIM_GENERIC,
        GOSLIM_GOA,
        GOSLIM_PIR,
        GOSUBSET_PROK,
        GOSLIM_CANDIDA,
        GOSLIM_ASPERGILLUS,
        GOSLIM_PLANT,
        GOSLIM_YEAST,
        GOSLIM_POMBE,
        OTHER;
    }
}
