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

package org.forester.archaeopteryx;

import java.awt.Color;
import java.awt.Dimension;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.util.ForesterConstants;

public final class AptxConstants {

    public final static String PRG_NAME = "Archaeopteryx";
    final static String VERSION = "0.11.117";
    final static String PRG_DATE = "2026-09-01";
    // The first three are bundled and registered at startup (see FontResources), so they are always
    // present and give identical, reproducible figure type across platforms; the rest are fallbacks.
    final static String[] DEFAULT_FONT_CHOICES = {
            "Source Sans 3", "Liberation Sans", "Noto Sans", "Arial Unicode MS", "Dialog", "SansSerif", "Sans",
            "Arial", "Helvetica"};
    // Default tree (large/tip-label) font size; small font is this minus 2. Tuned up from 10 for a more
    // readable, publication-ready look out of the box (Source Sans 3 reads well here).
    final static int DEFAULT_TREE_FONT_SIZE = 12;
    final static int DOMAIN_STRUCTURE_DEFAULT_WIDTH = 100;
    final static String AUTHOR_EMAIL = "czmasek AT jcvi DOT org";
    final static int DOMAIN_STRUCTURE_E_VALUE_THR_DEFAULT_EXP = -3;
    final static float BUTTON_ZOOM_IN_FACTOR = 1.25f;
    final static float BUTTON_ZOOM_OUT_FACTOR = 1 / AptxConstants.BUTTON_ZOOM_IN_FACTOR;
    final static float BUTTON_ZOOM_IN_X_CORRECTION_FACTOR = 1.2f;
    final static float BUTTON_ZOOM_OUT_X_CORRECTION_FACTOR = 1 / AptxConstants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR;
    final static float WHEEL_ZOOM_IN_FACTOR = 1.08f;
    final static float WHEEL_ZOOM_OUT_FACTOR = 1 / AptxConstants.WHEEL_ZOOM_IN_FACTOR;
    final static float WHEEL_ZOOM_IN_X_CORRECTION_FACTOR = 1.085f;
    final static float WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR = 1 / AptxConstants.WHEEL_ZOOM_IN_X_CORRECTION_FACTOR;
    static final double EXT_NODE_INFO_LENGTH_MAX_RATIO = 0.95;
    static final Dimension NODE_PANEL_SPLIT_MINIMUM_SIZE = new Dimension(100, 50);
    static final Dimension NODE_PANEL_SIZE = new Dimension(500, 540);
    static final Dimension NODE_FRAME_SIZE = new Dimension(520, 640);
    static final int MAX_TREES_TO_LOAD = 100;
    final static float PDF_LINE_WIDTH_DEFAULT = 0.5f;
    // The Archaeopteryx website. "Archaeopteryx Home" and "Documentation" both point here (there is no separate
    // documentation site yet); APTX_GITHUB is the source repository. The JS version points at its own SITE, where
    // the viewer actually runs in the browser -- not at its source repository, which is what a reader following a
    // "run it online" menu item is least likely to want.
    final static String APTX_WEB_SITE = "https://cmzmasek.github.io/archaeopteryx/";
    final static String APTX_JS_WEB_SITE = "https://cmzmasek.github.io/archaeopteryx-js/";
    final static String APTX_DOC_SITE = "https://cmzmasek.github.io/archaeopteryx/";
    final static String APTX_GITHUB = "https://github.com/cmzmasek/forester";
    final static String PHYLOXML_REFERENCE_URL = "https://www.biomedcentral.com/1471-2105/10/356/";
    final static String PHYLOXML_REFERENCE = ForesterConstants.PHYLO_XML_REFERENCE;
    final static String PHYLOXML_REFERENCE_SHORT = "Han MV and Zmasek CM (2009), BMC Bioinformatics, 10:356";
    final static short NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT = 3;
    final static short NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT = 2;
    public static final boolean NH_PARSING_IGNORE_QUOTES_DEFAULT = false;
    static final CLADOGRAM_TYPE CLADOGRAM_TYPE_DEFAULT = CLADOGRAM_TYPE.LINED_UP;
    final static boolean VALIDATE_AGAINST_PHYLOXML_XSD_SCJEMA_DEFAULT = true;
    final static String BACKUP_FILE_SUFFIX = ".BAK";
    final static double MIN_NOT_COLLAPSE_DEFAULT = 50;
    public final static Color DOMAIN_BASE_COLOR_FOR_PDF = new Color(100,
            100,
            100);
    public final static Color DOMAIN_LABEL_COLOR_FOR_PDF = new Color(0,
            0,
            0);
    final static short DEFAULT_NODE_SHAPE_SIZE_DEFAULT = 5;
    final static int   TIP_IMAGE_SIZE_DEFAULT           = 40; // tip-image target height (px); user-adjustable
    final static int   MSA_COLUMN_WIDTH_DEFAULT         = 7;  // alignment cell width (px/residue); user-adjustable
    final static int   MSA_COLUMN_WIDTH_MIN             = 1;
    final static int   MSA_COLUMN_WIDTH_MAX             = 40;
    // "Export at a fixed size" (ExportSizeSpec): DPI + width/height bounds and the default journal-figure size.
    final static int    EXPORT_SIZE_DPI_DEFAULT         = 300; // publication default
    final static int    EXPORT_SIZE_DPI_MIN             = 72;
    final static int    EXPORT_SIZE_DPI_MAX             = 1200;
    final static double EXPORT_SIZE_DIM_MIN             = 0.1;    // min width/height, in the selected unit
    final static double EXPORT_SIZE_DIM_MAX             = 20000;  // max width/height, in the selected unit
    final static double EXPORT_SIZE_WIDTH_MM_DEFAULT    = 170;    // double-column journal figure width (mm)
    final static double EXPORT_SIZE_HEIGHT_MM_DEFAULT   = 120;
    static final int MAX_LENGTH_FOR_COLLAPSED_NAME = 8;
    // External node names longer than this (e.g. whole UniProt/NCBI FASTA-header descriptions pasted in
    // as labels) are shown head + "…" when "Shorten Labels" is on. Also the threshold at which that
    // (data-gated) checkbox appears at all -- trees with sane labels never see it. Display-only.
    static final int LONG_NODE_NAME_LIMIT = 60;
    // Diameters (in tree coordinate space, so they scale with zoom like node shapes) of the
    // internal-node support symbols: the smallest dot drawn for low support in SIZE_SCALED mode,
    // and the fixed dot for THRESHOLD_MARKS / full support.
    public final static float SUPPORT_SYMBOL_MIN_DIAMETER = 2.0f;
    public final static float SUPPORT_SYMBOL_MAX_DIAMETER = 8.0f;
}
