// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// Copyright (C) 2003-2007 Ethalinda K.S. Cannon
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

package org.forester.archaeopteryx;

import java.awt.Color;
import java.awt.Dimension;

import org.forester.archaeopteryx.Options.CLADOGRAM_TYPE;
import org.forester.util.ForesterConstants;

public final class Constants {

    final static boolean        __ALLOW_PHYLOGENETIC_INFERENCE                                = true;
    public final static boolean __RELEASE                                                     = false;                                                               // TODO remove me
    public final static boolean __SNAPSHOT_RELEASE                                            = false;                                                               // TODO remove me
    public final static boolean __SYNTH_LF                                                    = false;                                                               // TODO remove me
    public final static boolean ALLOW_DDBJ_BLAST                                              = false;
    public final static String  PRG_NAME                                                      = "Archaeopteryx";
    final static String         VERSION                                                       = "0.985 CC";
    final static String         PRG_DATE                                                      = "140218";
    final static String         DEFAULT_CONFIGURATION_FILE_NAME                               = "_aptx_configuration_file";
    final static String[]       DEFAULT_FONT_CHOICES                                          = { "Arial", "Helvetica",
            "Verdana", "Tahoma", "Dialog", "Lucida Sans", "SansSerif", "Sans-serif", "Sans"  };
    final static boolean        VERBOSE_DEFAULT                                               = false;
    final static int            DOMAIN_STRUCTURE_DEFAULT_WIDTH                                = 100;
    final static String         AUTHOR_EMAIL                                                  = "phyloxml@gmail.com";
    final static int            DOMAIN_STRUCTURE_E_VALUE_THR_DEFAULT_EXP                      = -3;
    final static float          BUTTON_ZOOM_IN_FACTOR                                         = 1.25f;
    final static float          BUTTON_ZOOM_OUT_FACTOR                                        = 1 / Constants.BUTTON_ZOOM_IN_FACTOR;
    final static float          BUTTON_ZOOM_IN_X_CORRECTION_FACTOR                            = 1.2f;
    final static float          BUTTON_ZOOM_OUT_X_CORRECTION_FACTOR                           = 1 / Constants.BUTTON_ZOOM_IN_X_CORRECTION_FACTOR;
    final static float          WHEEL_ZOOM_IN_FACTOR                                          = 1.08f;
    final static float          WHEEL_ZOOM_OUT_FACTOR                                         = 1 / Constants.WHEEL_ZOOM_IN_FACTOR;
    final static float          WHEEL_ZOOM_IN_X_CORRECTION_FACTOR                             = 1.085f;
    final static float          WHEEL_ZOOM_OUT_X_CORRECTION_FACTOR                            = 1 / Constants.WHEEL_ZOOM_IN_X_CORRECTION_FACTOR;
    static final boolean        SPECIAL_CUSTOM                                                = false;                                                               //TODO remove me
    static final double         EXT_NODE_INFO_LENGTH_MAX_RATIO                                = 0.95;
    static final Dimension      NODE_PANEL_SPLIT_MINIMUM_SIZE                                 = new Dimension( 100, 50 );
    static final Dimension      NODE_PANEL_SIZE                                               = new Dimension( 500, 600 );
    static final Dimension      NODE_FRAME_SIZE                                               = new Dimension( 520, 640 );
    static final String         APPLET_PARAM_NAME_FOR_URL_OF_TREE_TO_LOAD                     = "url_of_tree_to_load";
    static final String         APPLET_PARAM_NAME_FOR_URL_OF_SPECIES_TREE_TO_LOAD             = "url_of_species_tree_to_load";
    static final String         APPLET_PARAM_NAME_FOR_CONFIG_FILE_URL                         = "config_file";
    static final String         APPLET_PARAM_NAME_FOR_DEFAULT_QUERY_SEQUENCE                  = "homology_type_analysis_query_sequence";
    static final String         APPLET_PARAM_NAME_FOR_DEFAULT_SEQUENCE_RELATION_TYPE          = "homology_type_analysis_initial_relation_type";
    static final int            MAX_TREES_TO_LOAD                                             = 100;
    static final int            US_LETTER_SIZE_X                                              = 612;
    static final int            US_LETTER_SIZE_Y                                              = 792;
    static final int            A4_SIZE_X                                                     = 595;
    static final int            A4_SIZE_Y                                                     = 845;
    final static float          PDF_LINE_WIDTH_DEFAULT                                        = 0.5f;
    final static String         APTX_WEB_SITE                                                 = "https://sites.google.com/site/cmzmasek/home/software/archaeopteryx";
    final static String         PHYLOXML_WEB_SITE                                             = ForesterConstants.PHYLO_XML_LOCATION;
    final static String         PHYLOXML_REFERENCE_URL                                        = "http://www.biomedcentral.com/1471-2105/10/356/";
    final static String         APTX_REFERENCE_URL                                            = "http://www.biomedcentral.com/bmcbioinformatics/";
    final static String         APTX_REFERENCE                                                = "Zmasek...";                                                         //TODO
    final static String         PHYLOXML_REFERENCE                                            = ForesterConstants.PHYLO_XML_REFERENCE;
    final static String         PHYLOXML_REFERENCE_SHORT                                      = "Han MV and Zmasek CM (2009), BMC Bioinformatics, 10:356";
    final static short          NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT = 2;
    final static short          NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT    = 1;
    public static final boolean NH_PARSING_IGNORE_QUOTES_DEFAULT                              = false;
    static final CLADOGRAM_TYPE CLADOGRAM_TYPE_DEFAULT                                        = CLADOGRAM_TYPE.EXT_NODE_SUM_DEP;
    final static boolean        VALIDATE_AGAINST_PHYLOXML_XSD_SCJEMA_DEFAULT                  = true;
    final static String         BACKUP_FILE_SUFFIX                                            = ".BAK";
    final static double         MIN_NOT_COLLAPSE_DEFAULT                                      = 50;
    final static Color          GUI_BACKGROUND_DEFAULT                                        = new Color( 32, 32, 32 );
    final static Color          CHECKBOX_TEXT_COLOR_DEFAULT                                   = new Color( 220,
                                                                                                           220,
                                                                                                           220 );
    final static Color          CHECKBOX_AND_BUTTON_ACTIVE_COLOR_DEFAULT                      = new Color( 255, 0, 0 );
    final static Color          BUTTON_TEXT_COLOR_DEFAULT                                     = new Color( 255,
                                                                                                           255,
                                                                                                           255 );
    final static Color          BUTTON_BACKGROUND_COLOR_DEFAULT                               = new Color( 64, 64, 64 );
    final static Color          MENU_BACKGROUND_COLOR_DEFAULT                                 = new Color( 0, 0, 0 );
    final static Color          MENU_TEXT_COLOR_DEFAULT                                       = new Color( 255,
                                                                                                           255,
                                                                                                           255 );
    final static Color          BUTTON_BORDER_COLOR_DEFAULT                                   = new Color( 0, 0, 0 );
    final static Color          TAB_LABEL_FOREGROUND_COLOR_SELECTED                           = new Color( 0, 0, 0 );
    final static Color          DOMAIN_STRUCTURE_BASE_COLOR_DEFAULT                           = new Color( 32, 32, 32 );
    final static Color          DOMAIN_STRUCTURE_FONT_COLOR_DEFAULT                           = new Color( 144,
                                                                                                           144,
                                                                                                           144 );
    final static String         NCBI_ALL_DATABASE_SEARCH                                      = "http://www.ncbi.nlm.nih.gov/gquery/?term=";
    final static short          DEFAULT_NODE_SHAPE_SIZE_DEFAULT                               = 4;
}
