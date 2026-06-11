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

package org.forester.archaeopteryx.phylogeny.data;

import java.awt.Dimension;
import java.awt.Graphics2D;

import org.forester.archaeopteryx.TreePanel;
import org.forester.phylogeny.data.PhylogenyData;

public interface RenderablePhylogenyData extends PhylogenyData {

    public Dimension getOriginalSize();

    public Object getParameter();

    public Dimension getRenderingSize();

    /**
     * This can be used to render phylogeny data as graphics (for example,
     * display of the domain structure). In most Renderable implementations this
     * will do nothing (i.e. just return).
     *
     * @param g
     *            the Graphics to render to
     */
    public void render( final float x, final float y, final Graphics2D g, final TreePanel tree_panel, boolean to_pdf );

    public void setParameter( final double parameter );

    public void setRenderingHeight( final float rendering_height );
}
