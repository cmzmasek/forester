// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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
// WWW: www.phylosoft.org/forester

package org.forester.phylogeny;

public interface PhylogenyNodeI {

    public enum NH_CONVERSION_SUPPORT_VALUE_STYLE {
        NONE, IN_SQUARE_BRACKETS, AS_INTERNAL_NODE_NAMES;
    }

    public void addAsChild( PhylogenyNodeI node );

    public PhylogenyNode getChildNode( int i );

    public double getDistanceToParent();

    public int getId();

    public String getName();

    public void setDistanceToParent( double d );

    public void setName( String name );

    public void setParent( PhylogenyNode phylogenyNode );
}
