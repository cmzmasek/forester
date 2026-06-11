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

public interface GoXRef extends Comparable<GoXRef> {

    public static final String EC_STR                = "EC";
    public static final String META_CYC_STR          = "MetaCyc";
    public static final String REACTOME_STR          = "Reactome";
    public static final String RESID_STR             = "RESID";
    public static final String UM_BBD_ENZYME_ID_STR  = "UM-BBD_enzymeID";
    public static final String UM_BBD_PATHWAY_ID_STR = "UM-BBD_pathwayID";
    public static final String UM_BBD_REACTIONID_STR = "UM-BBD_reactionID";
    public static final String TC_STR                = "TC";
    public static final String ARACYC_STR            = "AraCyc";
    public static final String XX_STR                = "XX";
    public static final String PMID_STR              = "PMID";
    public static final String IMG_STR               = "IMG";
    public static final String GOC_STR               = "GOC";
    public static final String WIKIPEDIA_STR         = "Wikipedia";
    public static final String KEGG_STR              = "KEGG";
    public static final String RHEA_STR              = "RHEA";
    public static final String NIF_SUBCELLULAR_STR   = "NIF_Subcellular";
    public static final String CORUM_STR             = "CORUM";
    public static final String UNIPATHWAY_STR        = "UniPathway";
    public static final String PO_STR                = "PO";
    public static final String SABIO_RK_STR          = "SABIO-RK";

    public Type getType();

    public String getXRef();

    public static enum Type {
        EC,
        META_CYC,
        REACTOME,
        RESID,
        UM_BBD_ENZYME_ID,
        UM_BBD_PATHWAY_ID,
        UM_BBD_REACTIONID,
        TC,
        ARACYC,
        XX,
        PMID,
        IMG,
        GOC,
        WIKIPEDIA,
        KEGG,
        RHEA,
        NIF_SUBCELLULAR,
        CORUM,
        UNIPATHWAY,
        PO,
        SABIO_RK,
        OTHER;
    }
}
