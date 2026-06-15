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
    final static String VERSION = "0.9.22";
    final static String PRG_DATE = "2026-06-15";
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
    static final boolean SPECIAL_CUSTOM = false;                                                                             //TODO remove me
    static final double EXT_NODE_INFO_LENGTH_MAX_RATIO = 0.95;
    static final Dimension NODE_PANEL_SPLIT_MINIMUM_SIZE = new Dimension(100, 50);
    static final Dimension NODE_PANEL_SIZE = new Dimension(500, 540);
    static final Dimension NODE_FRAME_SIZE = new Dimension(520, 640);
    static final int MAX_TREES_TO_LOAD = 100;
    final static float PDF_LINE_WIDTH_DEFAULT = 0.5f;
    final static String APTX_WEB_SITE = "https://sites.google.com/view/archaeopteryx/home";
    final static String APTX_JS_WEB_SITE = "https://sites.google.com/view/archaeopteryxjs";
    final static String APTX_DOC_SITE = "https://sites.google.com/view/cmzmasek/christian-zmasek/software/archaeopteryx/documentation";
    final static String PHYLOXML_REFERENCE_URL = "http://www.biomedcentral.com/1471-2105/10/356/";
    final static String APTX_REFERENCE_URL = "http://www.biomedcentral.com/bmcbioinformatics/";
    final static String APTX_REFERENCE = "Zmasek...";                                                                       //TODO
    final static String PHYLOXML_REFERENCE = ForesterConstants.PHYLO_XML_REFERENCE;
    final static String PHYLOXML_REFERENCE_SHORT = "Han MV and Zmasek CM (2009), BMC Bioinformatics, 10:356";
    final static short NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT = 3;
    final static short NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT = 2;
    public static final boolean NH_PARSING_IGNORE_QUOTES_DEFAULT = false;
    static final CLADOGRAM_TYPE CLADOGRAM_TYPE_DEFAULT = CLADOGRAM_TYPE.LINED_UP;
    final static boolean VALIDATE_AGAINST_PHYLOXML_XSD_SCJEMA_DEFAULT = true;
    final static String BACKUP_FILE_SUFFIX = ".BAK";
    final static double MIN_NOT_COLLAPSE_DEFAULT = 50;
    final static Color GUI_BACKGROUND_DEFAULT = new Color(32, 32, 32);
    final static Color CHECKBOX_TEXT_COLOR_DEFAULT = new Color(220,
            220,
            220);
    final static Color CHECKBOX_AND_BUTTON_ACTIVE_COLOR_DEFAULT = new Color(255, 0, 0);
    final static Color BUTTON_TEXT_COLOR_DEFAULT = new Color(255,
            255,
            255);
    final static Color BUTTON_BACKGROUND_COLOR_DEFAULT = new Color(64, 64, 64);
    final static Color MENU_BACKGROUND_COLOR_DEFAULT = new Color(0, 0, 0);
    final static Color MENU_TEXT_COLOR_DEFAULT = new Color(255,
            255,
            255);
    final static Color BUTTON_BORDER_COLOR_DEFAULT = new Color(0, 0, 0);
    final static Color TAB_LABEL_FOREGROUND_COLOR_SELECTED = new Color(0, 0, 0);
    public final static Color DOMAIN_BASE_COLOR_FOR_PDF = new Color(100,
            100,
            100);
    public final static Color DOMAIN_LABEL_COLOR_FOR_PDF = new Color(0,
            0,
            0);
    final static short DEFAULT_NODE_SHAPE_SIZE_DEFAULT = 7;
    static final int MAX_LENGTH_FOR_COLLAPSED_NAME = 8;
    // Diameters (in tree coordinate space, so they scale with zoom like node shapes) of the
    // internal-node support symbols: the smallest dot drawn for low support in SIZE_SCALED mode,
    // and the fixed dot for THRESHOLD_MARKS / full support.
    public final static float SUPPORT_SYMBOL_MIN_DIAMETER = 2.0f;
    public final static float SUPPORT_SYMBOL_MAX_DIAMETER = 8.0f;
}
