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
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.prefs.Preferences;

import org.forester.io.parsers.nhx.NHXParser.TAXONOMY_EXTRACTION;
import org.forester.util.ForesterUtil;

/**
 * Holds the runtime GUI theme (look-and-feel preference + custom cross-platform colors), the tip-facing color
 * maps, and the handful of file-reading/overview knobs that {@link TreePanel} and the parsers read directly.
 *
 * <p>Archaeopteryx no longer reads configuration files. The <em>default</em> display settings that used to be
 * copied from here into {@link Options} now live solely in {@code Options.init()} (backed by
 * {@link AptxConstants}) and in the {@link DisplayOption} / {@link ClickToOption} enums, so this class is no
 * longer a defaults-carrier for {@code Options} -- it retains only what is genuinely consumed at runtime.
 */
public final class Configuration {

    public enum UI {
        CROSSPLATFORM, NATIVE, NIMBUS, FLAT_LIGHT, FLAT_DARK, UNKNOWN
    }

    private static final String PREFS_NODE = "org/forester/archaeopteryx";
    private static final String PREFS_KEY_UI = "ui";

    // The click-to actions and the "Display Data" options are the type-safe enums {@link ClickToOption}
    // and {@link DisplayOption} (they replaced the old positional String[][] clickto_options /
    // display_options arrays and their parallel int index constants).

    private static Hashtable<String, Color> _domain_colors;
    private static Hashtable<String, Color> _species_colors;
    private static String DEFAULT_FONT_FAMILY = "";

    // The click-to action selected by default in the dropdown.
    private final ClickToOption _default_clickto = ClickToOption.DISPLAY_NODE_DATA;

    private SortedMap<String, Color> _display_colors = null;

    private final Color _gui_background_color = AptxConstants.GUI_BACKGROUND_DEFAULT;
    private final Color _gui_button_background_color = AptxConstants.BUTTON_BACKGROUND_COLOR_DEFAULT;
    private final Color _gui_button_border_color = AptxConstants.BUTTON_BORDER_COLOR_DEFAULT;
    private final Color _gui_button_text_color = AptxConstants.BUTTON_TEXT_COLOR_DEFAULT;
    private final Color _gui_checkbox_and_button_active_color = AptxConstants.CHECKBOX_AND_BUTTON_ACTIVE_COLOR_DEFAULT;
    private final Color _gui_checkbox_text_color = AptxConstants.CHECKBOX_TEXT_COLOR_DEFAULT;
    private final Color _gui_menu_background_color = AptxConstants.MENU_BACKGROUND_COLOR_DEFAULT;
    private final Color _gui_menu_text_color = AptxConstants.MENU_TEXT_COLOR_DEFAULT;

    // Read directly by the NH/NHX/Nexus parse path (see Archaeopteryx.main).
    private final boolean _internal_number_are_confidence_for_nh_parsing = false;
    private final boolean _nh_parsing_replace_underscores = false;
    private final TAXONOMY_EXTRACTION _taxonomy_extraction = TAXONOMY_EXTRACTION.NO;

    // Read directly by TreePanel.
    private final int _min_base_font_size = 2;
    private final short _number_of_digits_after_comma_for_branch_length_values = AptxConstants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT;
    private final short _number_of_digits_after_comma_for_confidence_values = AptxConstants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT;
    private final short _ov_max_height = 80;
    private final short _ov_max_width = 80;

    private UI _ui = UI.UNKNOWN;


    static {
        for (final String font_name : AptxConstants.DEFAULT_FONT_CHOICES) {
            if (Arrays.binarySearch(AptxUtil.getAvailableFontFamiliesSorted(), font_name) >= 0) {
                DEFAULT_FONT_FAMILY = font_name;
                break;
            }
        }
        if (ForesterUtil.isEmpty(DEFAULT_FONT_FAMILY)) {
            DEFAULT_FONT_FAMILY = AptxConstants.DEFAULT_FONT_CHOICES[AptxConstants.DEFAULT_FONT_CHOICES.length - 1];
        }
    }

    public Configuration() {
        // Archaeopteryx no longer reads configuration files; all display settings come from the built-in
        // defaults (Options.init + the DisplayOption enum) and the Settings dialog at runtime.
        setDisplayColors(new TreeMap<String, Color>());
    }

    public void setDisplayColors(final SortedMap<String, Color> display_colors) {
        _display_colors = display_colors;
    }

    private void initSpeciesColors() {
        _species_colors = new Hashtable<String, Color>();
    }


    ClickToOption getDefaultDisplayClicktoOption() {
        return _default_clickto;
    }

    SortedMap<String, Color> getDisplayColors() {
        return _display_colors;
    }

    Map<String, Color> getDomainColors() {
        if (_domain_colors == null) {
            _domain_colors = new Hashtable<String, Color>();
        }
        return _domain_colors;
    }

    Color getGuiBackgroundColor() {
        return _gui_background_color;
    }

    Color getGuiButtonBackgroundColor() {
        return _gui_button_background_color;
    }

    Color getGuiButtonBorderColor() {
        return _gui_button_border_color;
    }

    Color getGuiButtonTextColor() {
        return _gui_button_text_color;
    }

    Color getGuiCheckboxAndButtonActiveColor() {
        return _gui_checkbox_and_button_active_color;
    }

    Color getGuiCheckboxTextColor() {
        return _gui_checkbox_text_color;
    }

    Color getGuiMenuBackgroundColor() {
        return _gui_menu_background_color;
    }

    Color getGuiMenuTextColor() {
        return _gui_menu_text_color;
    }

    static int getGuiFontSize() {
        return 11;
    }

    int getMinBaseFontSize() {
        return _min_base_font_size;
    }

    short getNumberOfDigitsAfterCommaForBranchLengthValues() {
        return _number_of_digits_after_comma_for_branch_length_values;
    }

    short getNumberOfDigitsAfterCommaForConfidenceValues() {
        return _number_of_digits_after_comma_for_confidence_values;
    }

    short getOvMaxHeight() {
        return _ov_max_height;
    }

    short getOvMaxWidth() {
        return _ov_max_width;
    }

    Hashtable<String, Color> getSpeciesColors() {
        if (_species_colors == null) {
            initSpeciesColors();
        }
        return _species_colors;
    }

    TAXONOMY_EXTRACTION getTaxonomyExtraction() {
        return _taxonomy_extraction;
    }

    /**
     * Convenience method.
     *
     * @return true if the tree should be drawn as a phylogram (vs. cladogram) by default
     */
    boolean isDrawAsPhylogram() {
        return DisplayOption.DISPLAY_AS_PHYLOGRAM.isCheckedByDefault();
    }

    boolean isEditable() {
        return true;
    }


    boolean isInternalNumberAreConfidenceForNhParsing() {
        return _internal_number_are_confidence_for_nh_parsing;
    }

    boolean isReplaceUnderscoresInNhParsing() {
        return _nh_parsing_replace_underscores;
    }

    /**
     * Returns the resolved look-and-feel selection. When no preference has been set yet
     * ({@code UNKNOWN}), the last theme the user chose at runtime (persisted via
     * {@link java.util.prefs.Preferences}) is used; if none was ever saved, the modern
     * FlatLaf light theme is the default.
     */
    UI getUi() {
        if (_ui == UI.UNKNOWN) {
            _ui = readUiPreference();
        }
        return _ui;
    }

    void setUi(final UI ui) {
        _ui = ui;
    }

    private static UI readUiPreference() {
        try {
            final String saved = Preferences.userRoot().node(PREFS_NODE).get(PREFS_KEY_UI, null);
            if (saved != null) {
                final UI ui = UI.valueOf(saved);
                if (ui != UI.UNKNOWN) {
                    return ui;
                }
            }
        } catch (final Exception e) {
            // an invalid or inaccessible preference simply falls through to the default
        }
        return UI.FLAT_LIGHT;
    }

    /**
     * Persists the user's runtime look-and-feel choice so it survives a restart.
     */
    static void saveUiPreference(final UI ui) {
        try {
            Preferences.userRoot().node(PREFS_NODE).put(PREFS_KEY_UI, ui.name());
        } catch (final Exception e) {
            // failing to persist the preference is non-fatal
        }
    }

    boolean isUseNativeUI() {
        return getUi() == UI.NATIVE;
    }

    /**
     * Whether the legacy, hand-themed cross-platform GUI colors and fonts should be
     * applied to the Swing components. The native and FlatLaf look-and-feels style the
     * components themselves, so custom colors are only applied for {@code CROSSPLATFORM}.
     */
    boolean isApplyCustomGuiColors() {
        return getUi() == UI.CROSSPLATFORM;
    }


    boolean isValidatePhyloXmlAgainstSchema() {
        return true;
    }


    static String getDefaultFontFamilyName() {
        return DEFAULT_FONT_FAMILY;
    }

}
