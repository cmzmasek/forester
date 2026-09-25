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
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Window;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.SwingConstants;

/**
 * "Control Panel Cheat Sheet" (Help menu): one row per control, each showing the control's OWN icon beside its own
 * description.
 * <p>
 * Both halves are read off the live {@link ControlPanel} rather than written out again here: the icon is the very
 * {@link Icon} the button is painting, and the text is the description {@code ControlPanel.describe} already
 * attaches to each button as its tooltip and accessible name. So the sheet cannot drift from the panel -- add a
 * button, or reword a tooltip, and this follows without being edited. A hand-maintained list of "what the buttons
 * mean" is exactly the kind of second copy that goes quietly stale.
 * <p>
 * Controls the panel hides for the current tree (the data-gated Display Data rows) are skipped, so the sheet
 * describes what is actually on screen.
 */
final class ControlPanelCheatSheet {

    /** One row: the control's icon (may be null -- a checkbox or a slider has none) and what it does. */
    record Entry(Icon icon, String label, String description) {
    }

    /**
     * Walks the panel top to bottom -- which is its visual order, since the control panel is a single GridBag
     * column -- and collects every VISIBLE control that carries a description.
     */
    static List<Entry> entries(final ControlPanel cp) {
        final List<Entry> out = new ArrayList<Entry>();
        collect(cp, out, new String[ 1 ]);
        return out;
    }

    /**
     * {@code recent_label} carries the last plain {@link JLabel} seen in the walk. A dropdown or a slider has no
     * text of its own, but the panel puts its name in a label immediately above it ("Color by:", "Font size:"),
     * so the row reads as the user sees it rather than as "(dropdown)".
     */
    private static void collect(final Container c, final List<Entry> out, final String[] recent_label) {
        for (final Component comp : c.getComponents()) {
            if (!comp.isVisible()) {
                continue; // a Display Data row this tree has no data for: not on screen, not on the sheet
            }
            if ((comp instanceof JLabel) && !(comp instanceof AbstractButton)) {
                final String lt = ((JLabel) comp).getText();
                if ((lt != null) && !lt.isBlank()) {
                    recent_label[ 0 ] = lt.trim();
                }
            }
            final String tip = describedText(comp);
            if (tip != null) {
                out.add(new Entry(iconOf(comp), labelOf(comp, recent_label[ 0 ]), tip));
                // A control is a LEAF for this walk. A JComboBox (and a JSlider under some look-and-feels) is
                // itself a Container whose inner parts inherit the tooltip, so descending into one listed every
                // dropdown twice -- once as the control, once as its editor/renderer.
                continue;
            }
            if (comp instanceof Container) {
                collect((Container) comp, out, recent_label);
            }
        }
    }

    /**
     * The control's own description, or null when {@code comp} is not a control at all (a spacer, a panel, a bare
     * label). The description may be EMPTY: several Display Data checkboxes carry no tooltip because their own
     * label already says it ("Node Name", "Taxonomy Code", "Domain Architectures"). Requiring a tooltip dropped
     * exactly those from the sheet -- a cheat sheet that omits controls is not a cheat sheet -- so a control with
     * a name but no tooltip is listed by name.
     */
    private static String describedText(final Component comp) {
        if (!((comp instanceof AbstractButton) || (comp instanceof JComboBox) || (comp instanceof JSlider))) {
            return null;
        }
        final javax.swing.JComponent jc = (javax.swing.JComponent) comp;
        final String tip = jc.getToolTipText();
        if ((tip != null) && !tip.isBlank()) {
            return tip;
        }
        // no tooltip: keep it only if it names itself (an icon-only button with neither is not describable)
        if (comp instanceof AbstractButton) {
            final String txt = ((AbstractButton) comp).getText();
            if ((txt != null) && !txt.isBlank()) {
                return "";
            }
        }
        return null;
    }

    private static Icon iconOf(final Component comp) {
        return (comp instanceof AbstractButton) ? ((AbstractButton) comp).getIcon() : null;
    }

    private static String labelOf(final Component comp, final String recent_label) {
        if (comp instanceof AbstractButton) {
            final String t = ((AbstractButton) comp).getText();
            return ((t == null) || t.isBlank()) ? "" : t;
        }
        // A dropdown or slider carries no text; the panel names it with the label just above it. That label
        // includes the control's CURRENT VALUE ("Font size: 6", "Tree share: 40%") because the panel doubles as
        // a readout -- but a cheat sheet says what a control IS, not where it happens to be set, and a printed
        // sheet showing "40%" would be wrong for every reader but the one who printed it. Keep the name only.
        if ((recent_label == null) || recent_label.isBlank()) {
            return "";
        }
        final int colon = recent_label.indexOf(':');
        return (colon < 0) ? recent_label.trim() : recent_label.substring(0, colon).trim();
    }

    /** The sheet as a component -- also what the PNG export renders, so the file IS what the window shows. */
    static JPanel buildSheet(final ControlPanel cp) {
        final JPanel sheet = new JPanel();
        sheet.setLayout(new BoxLayout(sheet, BoxLayout.Y_AXIS));
        sheet.setBackground(Color.WHITE);
        sheet.setBorder(BorderFactory.createEmptyBorder(14, 16, 14, 16));
        final JLabel title = new JLabel("Control Panel");
        title.setFont(title.getFont().deriveFont(Font.BOLD, title.getFont().getSize() + 4f));
        title.setForeground(Color.BLACK);
        title.setAlignmentX(Component.LEFT_ALIGNMENT);
        sheet.add(title);
        final JLabel sub = new JLabel("Every control, with the icon it draws and what it does.");
        sub.setForeground(Color.DARK_GRAY);
        sub.setAlignmentX(Component.LEFT_ALIGNMENT);
        sheet.add(sub);
        sheet.add(Box.createVerticalStrut(10));
        for (final Entry e : entries(cp)) {
            sheet.add(row(e));
            sheet.add(Box.createVerticalStrut(3));
        }
        return sheet;
    }

    private static JPanel row(final Entry e) {
        final JPanel r = new JPanel();
        r.setLayout(new BoxLayout(r, BoxLayout.X_AXIS));
        r.setBackground(Color.WHITE);
        r.setAlignmentX(Component.LEFT_ALIGNMENT);
        final JLabel glyph = new JLabel();
        if (e.icon() != null) {
            glyph.setIcon(e.icon()); // the button's OWN icon instance -- never a copy that could differ
        }
        glyph.setPreferredSize(new Dimension(30, 20));
        glyph.setMinimumSize(new Dimension(30, 20));
        glyph.setMaximumSize(new Dimension(30, 20));
        glyph.setHorizontalAlignment(SwingConstants.CENTER);
        r.add(glyph);
        r.add(Box.createHorizontalStrut(8));
        // A COLON between the name and what it does, never a dash (house style; see CLAUDE.md).
        final String sep = (e.label().isEmpty() || e.description().isEmpty()) ? "" : ": ";
        final JLabel text = new JLabel("<html><body style='width:520px'><b>" + escape(e.label()) + escape(sep)
                + "</b>" + escape(e.description()) + "</body></html>");
        text.setForeground(Color.BLACK);
        r.add(text);
        r.add(Box.createHorizontalGlue());
        r.setMaximumSize(new Dimension(Integer.MAX_VALUE, r.getPreferredSize().height));
        return r;
    }

    private static String escape(final String s) {
        return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;");
    }

    static JDialog buildDialog(final Window owner, final ControlPanel cp) {
        final JDialog dialog = new JDialog(owner, "Control Panel Cheat Sheet");
        final JPanel sheet = buildSheet(cp);
        final JScrollPane scroll = new JScrollPane(sheet);
        scroll.getVerticalScrollBar().setUnitIncrement(16);
        dialog.getContentPane().add(scroll);
        final JPanel buttons = new JPanel();
        final JButton png = new JButton("Save as PNG…");
        png.addActionListener(a -> savePng(dialog, sheet));
        buttons.add(png);
        final JButton close = new JButton("Close");
        close.addActionListener(a -> dialog.dispose());
        buttons.add(close);
        dialog.getContentPane().add(buttons, java.awt.BorderLayout.SOUTH);
        dialog.pack();
        final Dimension d = dialog.getSize();
        dialog.setSize(Math.min(d.width + 30, 760), Math.min(d.height + 20, 800));
        dialog.setLocationRelativeTo(owner);
        return dialog;
    }

    /** Renders the very panel the window shows, so the picture and the window can never disagree. */
    static BufferedImage render(final JPanel sheet) {
        final Dimension d = sheet.getPreferredSize();
        final int w = Math.max(1, d.width);
        final int h = Math.max(1, d.height);
        sheet.setSize(w, h);
        sheet.doLayout();
        layoutDeep(sheet);
        final BufferedImage img = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        final Graphics2D g = img.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, w, h);
        sheet.printAll(g);
        g.dispose();
        return img;
    }

    /** A panel that has never been shown lays out only when told to, all the way down. */
    private static void layoutDeep(final Container c) {
        c.doLayout();
        for (final Component comp : c.getComponents()) {
            if (comp instanceof Container) {
                layoutDeep((Container) comp);
            }
        }
    }

    private static void savePng(final Window owner, final JPanel sheet) {
        final JFileChooser fc = new JFileChooser();
        fc.setDialogTitle("Save Control Panel Cheat Sheet");
        fc.setSelectedFile(new File("archaeopteryx-control-panel.png"));
        if (fc.showSaveDialog(owner) != JFileChooser.APPROVE_OPTION) {
            return;
        }
        try {
            ImageIO.write(render(sheet), "png", fc.getSelectedFile());
        }
        catch (final Exception e) {
            JOptionPane.showMessageDialog(owner, "Could not write the image: " + e.getMessage(), "Save failed",
                    JOptionPane.ERROR_MESSAGE);
        }
    }

    static void show(final MainFrame mf) {
        final ControlPanel cp = (mf.getMainPanel() == null) ? null : mf.getMainPanel().getControlPanel();
        if (cp == null) {
            return;
        }
        buildDialog(mf, cp).setVisible(true);
    }

    private ControlPanelCheatSheet() {
        // not instantiable
    }
}
