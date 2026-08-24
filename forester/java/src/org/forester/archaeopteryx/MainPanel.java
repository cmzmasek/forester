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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JViewport;
import javax.swing.SwingConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

public class MainPanel extends JPanel implements ComponentListener {

    private static final long serialVersionUID = -2682765312661416435L;
    MainFrame _mainframe;
    List<TreePanel> _treepanels;
    ControlPanel _control_panel;
    private List<JScrollPane> _treegraphic_scroll_panes;
    private List<JPanel> _treegraphic_scroll_pane_panels;
    Configuration _configuration;
    private JTabbedPane _tabbed_pane;
    private TreeColorSet _colorset;
    private TreeFontSet _fontset;
    private Phylogeny _cut_or_copied_tree;
    private Set<Long> _copied_and_pasted_nodes;

    public MainPanel(final Configuration configuration, final MainFrame parent) {
        if (configuration == null) {
            throw new IllegalArgumentException("configuration is null");
        }
        addComponentListener(this);
        _configuration = configuration;
        _mainframe = parent;
        _treepanels = new ArrayList<TreePanel>();
        initialize();
        _control_panel = new ControlPanel(this, configuration);
        // Wrap the control panel so it scrolls vertically when the window is too short to show all
        // of its controls (it lays its rows out at their natural height rather than cramming them).
        final JScrollPane control_panel_scroller = new JScrollPane(_control_panel,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        // Keep the look-and-feel's own scroll-pane border (a thin line separating the control panel from
        // the tree canvas). Do NOT null it out: the L&F installs that border at construction, but a runtime
        // theme switch (updateComponentTreeUI) re-installs it -- so a null here made the border show up only
        // after the first theme toggle or tree load, leaving a fresh, tree-less startup looking "incomplete".
        control_panel_scroller.getViewport().setBackground(_control_panel.getBackground());
        control_panel_scroller.getVerticalScrollBar().setUnitIncrement(16);
        add(control_panel_scroller, BorderLayout.WEST);
        setupTreeGraphic(configuration, getControlPanel());
        getControlPanel().showWhole();
    }

    MainPanel() {
    }

    public void addPhylogenyInNewTab(final Phylogeny phy,
                                     final Configuration config,
                                     final String default_name,
                                     final String full_path) {
        final TreePanel treepanel = new TreePanel(phy, config, this);
        // if this tree was saved after an annotation import, restore the remembered profile so File -> Re-import works
        // across a save/reload (the profile rides along as a property on the root node)
        treepanel.setLastImportProfile(org.forester.archaeopteryx.tools.NodeDataImporter.readProfileFromTree(phy));
        // restore a saved per-tree Time-Axis config (aptx:time_axis on the root); a saved config wins over the
        // auto-derived default, and null (no/unparsable property) leaves this tab on auto-derive
        treepanel.applyTimeAxisConfig(TimeAxisConfig.readFromTree(phy));
        getControlPanel().phylogenyAdded(config);
        treepanel.setControlPanel(getControlPanel());
        _treepanels.add(treepanel);
        final String name = tabTitleFor(phy, default_name, getTabbedPane().getTabCount() + 1);
        final JScrollPane treegraphic_scroll_pane = new JScrollPane(treepanel);
        // The overview, scale and tree name are painted at getVisibleRect-relative (viewport-fixed) positions.
        // With the default BLIT_SCROLL_MODE, dragging a scrollbar blits the existing pixels to a shifted
        // position (carrying those fixed items to the wrong place) and only repaints the newly-exposed strip,
        // so the stale overlay flickers for a frame before a full repaint corrects it -- the jitter seen only
        // when scrolling via the scrollbars. SIMPLE_SCROLL_MODE repaints the whole (node-culled) visible area on
        // every scroll instead, so the fixed items are always drawn at the right spot. No blit, no jitter.
        treegraphic_scroll_pane.getViewport().setScrollMode(JViewport.SIMPLE_SCROLL_MODE);
        // (No scrollbar AdjustmentListener repaints are needed: SIMPLE_SCROLL_MODE already repaints the whole
        // visible area on every scroll, which keeps the viewport-fixed overview/scale/name in step.)
        treegraphic_scroll_pane.getHorizontalScrollBar().setUnitIncrement(10);
        treegraphic_scroll_pane.getHorizontalScrollBar().setBlockIncrement(200);
        treegraphic_scroll_pane.getVerticalScrollBar().setUnitIncrement(10);
        treegraphic_scroll_pane.getVerticalScrollBar().setBlockIncrement(200);
        final JPanel treegraphic_scroll_pane_panel = new JPanel();
        treegraphic_scroll_pane_panel.setLayout(new BorderLayout());
        treegraphic_scroll_pane_panel.add(treegraphic_scroll_pane, BorderLayout.CENTER);
        _treegraphic_scroll_pane_panels.add(treegraphic_scroll_pane_panel);
        _treegraphic_scroll_panes.add(treegraphic_scroll_pane);
        getTabbedPane().addTab(name, null, treegraphic_scroll_pane_panel, "");
        getTabbedPane().setSelectedIndex(getTabbedPane().getTabCount() - 1);
        getControlPanel().showWhole();
    }

    @Override
    public void componentHidden(final ComponentEvent e) {
        // Do nothing.
    }

    @Override
    public void componentMoved(final ComponentEvent e) {
        // Do nothing.
    }

    @Override
    public void componentResized(final ComponentEvent e) {
        if (getCurrentTreePanel() != null) {
            getCurrentTreePanel().updateOvSettings();
            getCurrentTreePanel().updateOvSizes();
        }
    }

    @Override
    public void componentShown(final ComponentEvent e) {
        // Do nothing.
    }

    public ControlPanel getControlPanel() {
        return _control_panel;
    }

    public Set<Long> getCopiedAndPastedNodes() {
        return _copied_and_pasted_nodes;
    }

    public TreePanel getCurrentTreePanel() {
        final int selected = getTabbedPane().getSelectedIndex();
        if (selected >= 0) {
            return _treepanels.get(selected);
        } else {
            if (_treepanels.size() == 1) {
                return _treepanels.get(0);
            } else {
                return null;
            }
        }
    }

    public Options getOptions() {
        return _mainframe.getOptions();
    }

    public JTabbedPane getTabbedPane() {
        return _tabbed_pane;
    }

    public TreeFontSet getTreeFontSet() {
        return _fontset;
    }

    public void setArrowCursor() {
        setCursor(TreePanel.ARROW_CURSOR);
        repaint();
    }

    public void setCopiedAndPastedNodes(final Set<Long> node_ids) {
        _copied_and_pasted_nodes = node_ids;
    }

    public void setWaitCursor() {
        setCursor(TreePanel.WAIT_CURSOR);
        repaint();
    }


    void adjustJScrollPane() {
        if (getTabbedPane() != null) {
            getCurrentScrollPanePanel().remove(getCurrentScrollPane());
            getCurrentScrollPanePanel().add(getCurrentScrollPane(), BorderLayout.CENTER);
        }
        getCurrentScrollPane().revalidate();
    }

    void closeCurrentPane() {
        final int index = getCurrentTabIndex();
        if ((index >= 0) && (getTabbedPane().getTabCount() > 0)) {
            getTabbedPane().remove(index);
            getTreePanels().remove(index);
            _treegraphic_scroll_panes.remove(index);
            _treegraphic_scroll_pane_panels.remove(index);
            getControlPanel().phylogenyRemoved(index);
        }
    }

    Configuration getConfiguration() {
        return _configuration;
    }

    Phylogeny getCurrentPhylogeny() {
        if (getCurrentTreePanel() == null) {
            return null;
        }
        return getCurrentTreePanel().getPhylogeny();
    }

    JScrollPane getCurrentScrollPane() {
        if (_treegraphic_scroll_panes.size() > 0) {
            final int selected = getTabbedPane().getSelectedIndex();
            if (selected >= 0) {
                return _treegraphic_scroll_panes.get(selected);
            } else {
                return _treegraphic_scroll_panes.get(0);
            }
        } else {
            return null;
        }
    }

    JPanel getCurrentScrollPanePanel() {
        final int selected = getTabbedPane().getSelectedIndex();
        if (selected >= 0) {
            return _treegraphic_scroll_pane_panels.get(selected);
        } else {
            return _treegraphic_scroll_pane_panels.get(0);
        }
    }

    int getCurrentTabIndex() {
        final int selected = getTabbedPane().getSelectedIndex();
        if (selected >= 0) {
            return selected;
        } else {
            return 0;
        }
    }

    Phylogeny getCutOrCopiedTree() {
        return _cut_or_copied_tree;
    }



    MainFrame getMainFrame() {
        return _mainframe;
    }

    Phylogeny getPhylogeny(final int index) {
        if (getCurrentTreePanel() == null) {
            return null;
        }
        return _treepanels.get(index).getPhylogeny();
    }

    Dimension getSizeOfViewport() {
        return getCurrentScrollPane().getViewport().getExtentSize();
    }

    TreeColorSet getTreeColorSet() {
        return _colorset;
    }

    List<TreePanel> getTreePanels() {
        return _treepanels;
    }

    void initialize() {
        setTreeFontSet(new TreeFontSet(this));
        getTreeFontSet().setBaseFont(getOptions().getBaseFont());
        setLayout(new BorderLayout());
        setTreeColorSet(TreeColorSet.createInstance(getConfiguration()));
        _treegraphic_scroll_panes = new ArrayList<JScrollPane>();
        _treegraphic_scroll_pane_panels = new ArrayList<JPanel>();
        _tabbed_pane = new JTabbedPane(SwingConstants.TOP);
        _tabbed_pane.addChangeListener(new ChangeListener() {

            // This method is called whenever the selected tab changes
            @Override
            public void stateChanged(final ChangeEvent evt) {
                getControlPanel().tabChanged();
            }
        });
        // A tab's title IS the tree name (save writes it back onto the tree). Double-clicking a tab opens the
        // editor for that name and the tree description; the tab right-click menu offers the same item (see
        // MainFrameApplication.createTabPopupMenu).
        _tabbed_pane.addMouseListener(new MouseAdapter() {

            @Override
            public void mouseClicked(final MouseEvent e) {
                if ((e.getClickCount() == 2) && (e.getButton() == MouseEvent.BUTTON1)) {
                    final int i = _tabbed_pane.indexAtLocation(e.getX(), e.getY());
                    if ((i >= 0) && (_mainframe != null)) {
                        _tabbed_pane.setSelectedIndex(i);
                        _mainframe.showTreeInfoDialog();
                    }
                }
            }
        });
        _tabbed_pane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
        add(_tabbed_pane, BorderLayout.CENTER);
    }

    void setCutOrCopiedTree(final Phylogeny cut_or_copied_tree) {
        _cut_or_copied_tree = cut_or_copied_tree;
    }


    void setTitleOfSelectedTab(final String title) {
        final int selected = getTabbedPane().getSelectedIndex();
        if (selected >= 0) {
            getTabbedPane().setTitleAt(selected, title);
        }
    }

    /**
     * The label a tab carries for {@code phy}: its name, else its identifier, else {@code fallback} (e.g. the
     * loaded file's base name), else a positional "[n]" placeholder. Shared by tab creation and by the
     * undo/redo tab-title re-sync, so an unnamed tree always gets a stable, never-empty label -- never a stale
     * one left over from a since-undone rename.
     */
    static String tabTitleFor(final Phylogeny phy, final String fallback, final int positional) {
        if (!ForesterUtil.isEmpty(phy.getName())) {
            return phy.getName();
        }
        if (phy.getIdentifier() != null) {
            return phy.getIdentifier().toString();
        }
        if (!ForesterUtil.isEmpty(fallback)) {
            return fallback;
        }
        return "[" + positional + "]";
    }

    void setTreeColorSet(final TreeColorSet colorset) {
        _colorset = colorset;
        for (final TreePanel p : getTreePanels()) {
            p.setBackground(colorset.getBackgroundColor());
        }
    }

    void setTreeFontSet(final TreeFontSet fontset) {
        _fontset = fontset;
    }

    void setupTreeGraphic(final Configuration config_settings, final ControlPanel control) {
        control.setSpeciesColors(config_settings.getSpeciesColors());
        RenderableDomainArchitecture.setColorMap(config_settings.getDomainColors());
    }

    void terminate() {
        for (final TreePanel atvtreepanel : _treepanels) {
            atvtreepanel.removeAllEditNodeJFrames();
        }
    }

}
