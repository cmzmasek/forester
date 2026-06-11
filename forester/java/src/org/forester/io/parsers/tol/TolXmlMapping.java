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

package org.forester.io.parsers.tol;

public final class TolXmlMapping {

    public static final String PHYLOGENY               = "TREE";
    public static final String CLADE                   = "NODE";
    public static final String AUTHDATE                = "AUTHDATE";
    public static final String AUTHORITY               = "AUTHORITY";
    public static final String TAXONOMY_NAME           = "NAME";
    public static final String OTHERNAMES              = "OTHERNAMES";
    public static final String OTHERNAME               = "OTHERNAME";
    public static final String OTHERNAME_NAME          = "NAME";
    public static final String NODE_ID_ATTR            = "ID";
    public static final String NODE_ITALICIZENAME_ATTR = "ITALICIZENAME";
    public static final String TOL_TAXONOMY_ID_TYPE    = "tol";

    private TolXmlMapping() {
        // Hidden.
    }
}