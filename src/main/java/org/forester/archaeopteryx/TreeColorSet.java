// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2010 Christian M. Zmasek
// Copyright (C) 2008-2010 Burnham Institute for Medical Research
// Copyright (C) 2003-2010 Ethalinda K.S. Cannon
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
import java.util.Map;

import org.forester.util.ForesterUtil;

public final class TreeColorSet {

    public static final String ANNOTATION                 = "Annotation";
    public static final String BACKGROUND                 = "Background";
    public static final String BACKGROUND_GRADIENT_BOTTOM = "Background Gradient Bottom";
    public static final String BINARY_DOMAIN_COMBINATIONS = "Binary Domain Combinations";
    public static final String BRANCH                     = "Branch";
    public static final String BRANCH_LENGTH              = "Branch Length";
    public static final String COLLAPSED                  = "Collapsed";
    public static final String CONFIDENCE                 = "Confidence";
    public static final String DOMAIN_LABEL               = "Domain Label";
    public static final String DOMAIN_BASE                = "Domain Base";
    public static final String DUPLICATION                = "Duplication";
    public static final String DUPLICATION_OR_SPECATION   = "Duplication or Specation";
    public static final String MATCHING_NODES_A           = "Matching A";
    public static final String MATCHING_NODES_A_AND_B     = "Matching A and B";
    public static final String MATCHING_NODES_B           = "Matching B";
    public static final String NODE_BOX                   = "Node Box";
    public static final String OVERVIEW                   = "Overview";
    public static final String SEQUENCE                   = "Sequence";
    public static final String SPECIATION                 = "Speciation";
    public static final String TAXONOMY                   = "Taxonomy";
    static final String[]      COLOR_FIELDS               = { BACKGROUND, BACKGROUND_GRADIENT_BOTTOM, SEQUENCE,
            TAXONOMY, CONFIDENCE, BRANCH_LENGTH, BRANCH, NODE_BOX, COLLAPSED, MATCHING_NODES_A, MATCHING_NODES_B,
            MATCHING_NODES_A_AND_B, DUPLICATION, SPECIATION, DUPLICATION_OR_SPECATION, DOMAIN_LABEL, DOMAIN_BASE,
            BINARY_DOMAIN_COMBINATIONS, ANNOTATION, OVERVIEW };
    static final String[]      SCHEME_NAMES               = { "Default", "Black", "Black & White", "Silver", "Green",
            "White & Blue", "Cyan", "Orange", "Blue", "Blue & White", "Neon" };
    private int                _color_scheme;
    private final Color[][]    _color_schemes             = { { new Color( 0, 0, 0 ), // background_color
            new Color( 0, 100, 100 ), // background_color_gradient_bottom
            new Color( 230, 230, 230 ), // sequence  __ Default (same as Black)
            new Color( 180, 180, 180 ), // taxonomy
            new Color( 180, 180, 180 ), // support
            new Color( 140, 140, 140 ), // branch_length_color
            new Color( 255, 255, 255 ), // branch_color
            new Color( 255, 255, 255 ), // box_color
            new Color( 255, 255, 255 ), // collapesed_fill_color
            new Color( 0, 255, 0 ), // found_color 0
            new Color( 255, 0, 0 ), // found_color 1
            new Color( 255, 255, 0 ), // found_color 1 + 2
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 230, 230, 230 ), // domain_label
            new Color( 100, 100, 100 ), // domains_base
            new Color( 65, 105, 255 ), // binary_domain_combinations_color
            new Color( 173, 255, 47 ) // annotation
            , new Color( 130, 130, 130 )                 // overview
            }, { new Color( 0, 0, 0 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 230, 230, 230 ), // sequence  __ Black
            new Color( 180, 180, 180 ), // taxonomy
            new Color( 180, 180, 180 ), // support
            new Color( 140, 140, 140 ), // branch_length_color
            new Color( 255, 255, 255 ), // branch_color
            new Color( 255, 255, 255 ), // box_color
            new Color( 255, 255, 255 ), // collapesed_fill_color
            new Color( 0, 255, 0 ), // found_color 0
            new Color( 255, 0, 0 ), // found_color 1
            new Color( 255, 255, 0 ), // found_color 1 + 2
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 230, 230, 230 ), // domain_label
            new Color( 100, 100, 100 ), // domains_base
            new Color( 65, 105, 255 ), // binary_domain_combinations_color
            new Color( 173, 255, 47 ) // annotation
            , new Color( 130, 130, 130 ) // ov
            }, { new Color( 255, 255, 255 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 0, 0, 0 ), // sequence  __ Black & White
            new Color( 0, 0, 0 ), // taxonomy
            new Color( 0, 0, 0 ), // support
            new Color( 0, 0, 0 ), // branch_length_color
            new Color( 0, 0, 0 ), // branch_color
            new Color( 0, 0, 0 ), // box_color
            new Color( 0, 0, 0 ), // collapesed_fill_color
            new Color( 255, 0, 0 ), // found_color 0
            new Color( 0, 255, 0 ), // found_color 1
            new Color( 0, 0, 255 ), // found_color 1 + 2
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 0, 0, 0 ), // domain_label
            new Color( 100, 100, 100 ), // domains_base
            new Color( 0, 0, 0 ), // binary_domain_combinations_color
            new Color( 0, 0, 0 ) // annotation
            , new Color( 220, 220, 220 ) // ov
            }, { new Color( 0, 0, 0 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 220, 220, 220 ), // sequence __ Silver
            new Color( 180, 180, 180 ), // taxonomy
            new Color( 140, 140, 140 ), // support
            new Color( 140, 140, 140 ), // branch_length_color
            new Color( 240, 240, 240 ), // branch_color
            new Color( 140, 140, 140 ), // box_color
            new Color( 240, 240, 240 ), // collapesed_fill_color
            new Color( 255, 0, 0 ), // found_color 0
            new Color( 0, 255, 0 ), // found_color 1
            new Color( 255, 255, 0 ), // found_color 1 + 2
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 230, 230, 230 ), // domain_label
            new Color( 100, 100, 100 ), // domains_base
            new Color( 180, 180, 180 ), // binary_domain_combinations_color
            new Color( 140, 140, 140 ) // annotation
            , new Color( 40, 40, 40 ) // ov
            }, { new Color( 0, 10, 0 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 0, 255, 0 ), // sequence __ the Matrix
            new Color( 30, 200, 30 ), // taxonomy
            new Color( 0, 155, 0 ), // support
            new Color( 0, 100, 0 ), // branch_length_color
            new Color( 0, 155, 0 ), // branch_color
            new Color( 0, 255, 0 ), // box_color
            new Color( 0, 155, 0 ), // collapesed_fill_color
            new Color( 255, 0, 0 ), // found_color 0
            new Color( 0, 255, 0 ), // found_color 1
            new Color( 255, 255, 0 ), // found_color 1 + 2
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 230, 230, 230 ), // domain_label
            new Color( 100, 100, 100 ), // domains_base
            new Color( 0, 235, 0 ), // binary_domain_combinations_color
            new Color( 0, 235, 0 ) // annotation
            , new Color( 40, 40, 40 ) // ov
            }, { new Color( 255, 255, 255 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 0, 0, 0 ), //sequence __ White & Blue
            new Color( 40, 40, 40 ), // taxonomy
            new Color( 0, 125, 0 ), // support
            new Color( 70, 70, 0 ), // branch_length_color
            new Color( 0, 20, 200 ), // branch_color
            new Color( 0, 20, 200 ), // box_color
            new Color( 0, 20, 200 ), // collapesed_fill_color
            new Color( 0, 255, 0 ), // found_color 0
            new Color( 255, 0, 0 ), // found_color 1
            new Color( 0, 0, 255 ), // found_color 0 + 1
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 0, 0, 0 ), // domain_label
            new Color( 50, 50, 50 ), // domains_base
            new Color( 65, 105, 225 ), // binary_domain_combinations_color
            new Color( 173, 255, 47 ) // annotation
            , new Color( 220, 220, 220 ) // ov
            }, { new Color( 0, 0, 0 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 255, 255, 255 ), // sequence __ Cyan
            new Color( 200, 200, 200 ), // taxonomy
            new Color( 255, 255, 255 ), // support
            new Color( 200, 200, 200 ), // branch_length_color
            new Color( 0, 255, 255 ), // branch_color
            new Color( 0, 255, 255 ), // box_color
            new Color( 0, 255, 255 ), // collapesed_fill_color
            new Color( 0, 255, 0 ), // found_color 0
            new Color( 0, 0, 255 ), // found_color 1
            new Color( 0, 255, 255 ), // found_color 0 + 1
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 230, 230, 230 ), // domain_label
            new Color( 100, 100, 100 ), // domains_base
            new Color( 65, 105, 225 ), // binary_domain_combinations_color
            new Color( 173, 255, 47 ) // annotation
            , new Color( 0, 120, 120 ) // ov
            }, { new Color( 0, 0, 0 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 255, 200, 0 ), // sequence __ Clockwork
            new Color( 255, 200, 0 ), // taxonomy
            new Color( 255, 200, 0 ), // support
            new Color( 255, 200, 0 ), // branch_length_color
            new Color( 255, 200, 0 ), // branch_color
            new Color( 255, 200, 0 ), // box_color
            new Color( 255, 200, 0 ), // collapesed_fill_color
            new Color( 255, 255, 0 ), // found_color 0
            new Color( 0, 255, 255 ), // found_color 1
            new Color( 255, 255, 255 ), // found_color 0 + 1
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 255, 200, 0 ), // domain_label
            new Color( 255, 200, 0 ), // domains_base
            new Color( 150, 150, 150 ), // binary_domain_combinations_color
            new Color( 150, 150, 150 ) // annotation
            , new Color( 150, 150, 150 ) // ov
            }, { new Color( 0, 0, 100 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 255, 255, 255 ), // sequence __ Blue
            new Color( 255, 255, 255 ), // taxonomy
            new Color( 255, 0, 0 ), // support
            new Color( 255, 0, 0 ), // branch_length_color
            new Color( 255, 0, 0 ), // branch_color
            new Color( 255, 0, 0 ), // box_color
            new Color( 255, 0, 0 ), // collapesed_fill_color
            new Color( 0, 255, 0 ), // found_color
            new Color( 255, 0, 0 ), // found_color 1
            new Color( 255, 255, 0 ), // found_color 1 + 2
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 255, 255, 255 ), // domain_label
            new Color( 100, 100, 100 ), // domains_base
            new Color( 255, 255, 255 ), // binary_domain_combinations_color
            new Color( 255, 255, 255 ) // annotation
            , new Color( 77, 77, 255 ) // ov
            }, { new Color( 0, 0, 0 ), // background_color
            new Color( 0, 255, 255 ), // background_color_gradient_bottom
            new Color( 255, 255, 255 ), // sequence __ blue &  white
            new Color( 255, 255, 255 ), // taxonomy
            new Color( 255, 255, 255 ), // support
            new Color( 0, 191, 255 ), // branch_length_color
            new Color( 0, 191, 255 ), // branch_color
            new Color( 0, 191, 255 ), // box_color
            new Color( 0, 191, 255 ), // collapesed_fill_color
            new Color( 255, 0, 0 ), // found_color 0
            new Color( 0, 255, 0 ), // found_color 1
            new Color( 255, 255, 0 ), // found_color 0 + 1
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 255, 255, 255 ), // domain_label
            new Color( 150, 150, 150 ), // domains_base
            new Color( 255, 255, 255 ), // binary_domain_combinations_color
            new Color( 255, 255, 255 ) // annotation
            , new Color( 170, 187, 204 ) // ov
            }, { new Color( 0, 0, 0 ), // background_color
            new Color( 255, 255, 0 ), // background_color_gradient_bottom
            new Color( 127, 255, 0 ), // sequence __ Neon
            new Color( 255, 110, 199 ), // taxonomy
            new Color( 234, 173, 234 ), // support
            new Color( 77, 77, 255 ), // branch_length_color
            new Color( 234, 173, 234 ), // branch_color
            new Color( 77, 77, 255 ), // box_color
            new Color( 234, 173, 234 ), // collapsed_fill_color
            new Color( 243, 243, 21 ), // found_color 0
            new Color( 255, 20, 147 ), // found_color 1
            new Color( 255, 255, 255 ), // found_color 1 + 2
            new Color( 255, 0, 0 ), // duplication_box_color
            new Color( 0, 255, 0 ), // speciation_box_color
            new Color( 255, 255, 0 ), // duplication_speciation_color
            new Color( 127, 255, 0 ), // domain_label
            new Color( 234, 173, 234 ), // domains_base
            new Color( 27, 255, 0 ), // binary_domain_combinations_color
            new Color( 27, 255, 0 ) // annotation
            , new Color( 77, 77, 255 ) // ov
            }                                            };
    private Color              annotation_color;
    private Color              background_color;
    private Color              background_color_gradient_bottom;
    private Color              binary_domain_combinations_color;
    private Color              bootstrap_color;
    private Color              box_color;
    private Color              branch_color;
    private Color              branch_length_color;
    private Color              collapse_fill_color;
    private Color              domain_label_color;
    private Color              domain_base_color;
    private Color              dup_box_color;
    private Color              duplication_or_specation_color;
    private Color              found_color_0;
    private Color              found_color_0_and_1;
    private Color              found_color_1;
    private Color              ov_color;
    // The drawing colors
    private Color              seq_color;
    private Color              spec_box_color;
    private Color              taxonomy_color;

    private TreeColorSet() {
        // Hidden constructor.
    }

    public Color getDomainBaseColor() {
        return domain_base_color;
    }

    public Color getDomainLabelColor() {
        return domain_label_color;
    }

    private void setColorForDefault( final int i, final Color color ) {
        _color_schemes[ 0 ][ i ] = color;
    }

    void cycleColorScheme() {
        if ( getCurrentColorScheme() >= ( _color_schemes.length - 1 ) ) {
            setColorSchema( 0 );
        }
        else {
            setColorSchema( getCurrentColorScheme() + 1 );
        }
    }

    Color getAnnotationColor() {
        return annotation_color;
    }

    Color getBackgroundColor() {
        return background_color;
    }

    Color getBackgroundColorGradientBottom() {
        return background_color_gradient_bottom;
    }

    Color getBinaryDomainCombinationsColor() {
        if ( Constants.SPECIAL_CUSTOM ) {
            return new Color( 50, 50, 50 );
        }
        return binary_domain_combinations_color;
    }

    Color getBoxColor() {
        return box_color;
    }

    Color getBranchColor() {
        return branch_color;
    }

    Color getBranchColorForPdf() {
        return Color.BLACK;
    }

    Color getBranchLengthColor() {
        return branch_length_color;
    }

    Color getCollapseFillColor() {
        return collapse_fill_color;
    }

    Color[][] getColorSchemes() {
        return _color_schemes;
    }

    Color getConfidenceColor() {
        return bootstrap_color;
    }

    int getCurrentColorScheme() {
        return _color_scheme;
    }

    String getCurrentColorSchemeName() {
        return SCHEME_NAMES[ getCurrentColorScheme() ];
    }

    Color getDuplicationBoxColor() {
        return dup_box_color;
    }

    Color getDuplicationOrSpeciationColor() {
        return duplication_or_specation_color;
    }

    Color getFoundColor0() {
        return found_color_0;
    }

    Color getFoundColor0and1() {
        return found_color_0_and_1;
    }

    Color getFoundColor1() {
        return found_color_1;
    }

    Color getGainedCharactersColor() {
        return Color.GREEN;
    }

    Color getLostCharactersColor() {
        return Color.RED;
    }

    Color getOvColor() {
        return ov_color;
    }

    Color getSequenceColor() {
        return seq_color;
    }

    Color getSpecBoxColor() {
        return spec_box_color;
    }

    Color getTaxonomyColor() {
        return taxonomy_color;
    }

    void setColorforDefault( final String color_field_name, final Color color ) {
        final String query = color_field_name.trim().replace( '_', ' ' );
        boolean found = false;
        int i = 0;
        for( final String cf : COLOR_FIELDS ) {
            if ( query.equalsIgnoreCase( cf ) ) {
                found = true;
                setColorForDefault( i, color );
                break;
            }
            ++i;
        }
        if ( !found ) {
            throw new IllegalArgumentException( "unknown color field name [" + color_field_name + "]" );
        }
    }

    /**
     * Switches colors between different schemes.
     */
    void setColorSchema( final int scheme ) {
        _color_scheme = scheme;
        background_color = _color_schemes[ scheme ][ 0 ];
        background_color_gradient_bottom = _color_schemes[ scheme ][ 1 ];
        seq_color = _color_schemes[ scheme ][ 2 ];
        taxonomy_color = _color_schemes[ scheme ][ 3 ];
        bootstrap_color = _color_schemes[ scheme ][ 4 ];
        branch_length_color = _color_schemes[ scheme ][ 5 ];
        branch_color = _color_schemes[ scheme ][ 6 ];
        box_color = _color_schemes[ scheme ][ 7 ];
        collapse_fill_color = _color_schemes[ scheme ][ 8 ];
        found_color_0 = _color_schemes[ scheme ][ 9 ];
        found_color_1 = _color_schemes[ scheme ][ 10 ];
        found_color_0_and_1 = _color_schemes[ scheme ][ 11 ];
        dup_box_color = _color_schemes[ scheme ][ 12 ];
        spec_box_color = _color_schemes[ scheme ][ 13 ];
        duplication_or_specation_color = _color_schemes[ scheme ][ 14 ];
        domain_label_color = _color_schemes[ scheme ][ 15 ];
        domain_base_color = _color_schemes[ scheme ][ 16 ];
        binary_domain_combinations_color = _color_schemes[ scheme ][ 17 ];
        annotation_color = _color_schemes[ scheme ][ 18 ];
        ov_color = _color_schemes[ scheme ][ 19 ];
    }

    void setCurrentColorScheme( final int color_scheme ) {
        _color_scheme = color_scheme;
    }

    static TreeColorSet createInstance() {
        final TreeColorSet tcs = new TreeColorSet();
        tcs.setColorSchema( 0 );
        return tcs;
    }

    static TreeColorSet createInstance( final Configuration configuration ) {
        final TreeColorSet tcs = new TreeColorSet();
        if ( ( configuration != null ) && ( configuration.getDisplayColors() != null )
                && ( configuration.getDisplayColors().size() > 0 ) ) {
            final Map<String, Color> colors = configuration.getDisplayColors();
            for( final String field : colors.keySet() ) {
                final Color color = colors.get( field );
                try {
                    tcs.setColorforDefault( field, color );
                }
                catch ( final IllegalArgumentException ex ) {
                    ForesterUtil.printWarningMessage( Constants.PRG_NAME, ex.getMessage() );
                }
            }
        }
        tcs.setColorSchema( 0 );
        return tcs;
    }
}
