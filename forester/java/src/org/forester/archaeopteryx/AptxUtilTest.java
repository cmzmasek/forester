// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2026 Christian M. Zmasek
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
// Contact: phylosoft @ gmail . com

package org.forester.archaeopteryx;

import java.util.ArrayList;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

/**
 * Unit tests for {@link AptxUtil}. Lives in the {@code org.forester.archaeopteryx} package. Run
 * standalone via {@link #main(String[])}, or as part of the suite via {@link #test()}.
 */
public final class AptxUtilTest {

    public static void main( final String[] args ) {
        System.out.println( "AptxUtil: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return testHasAtLeastOneNodeWithDomainArchitecture();
    }

    private static boolean testHasAtLeastOneNodeWithDomainArchitecture() {
        // a tree whose leaves carry no sequence/domain data
        final Phylogeny no_domains = new Phylogeny();
        final PhylogenyNode root0 = new PhylogenyNode();
        root0.addAsChild( new PhylogenyNode() );
        root0.addAsChild( new PhylogenyNode() );
        no_domains.setRoot( root0 );
        no_domains.externalNodesHaveChanged();
        if ( AptxUtil.isHasAtLeastOneNodeWithDomainArchitecture( no_domains ) ) {
            return fail( "a tree without domain architectures must not be detected as having them" );
        }
        // give a single leaf a sequence with a domain architecture
        final Phylogeny with_domains = new Phylogeny();
        final PhylogenyNode root1 = new PhylogenyNode();
        final PhylogenyNode leaf = new PhylogenyNode();
        final ArrayList<PhylogenyData> domains = new ArrayList<PhylogenyData>();
        domains.add( new ProteinDomain( "d0", 10, 20 ) );
        domains.add( new ProteinDomain( "d1", 30, 40 ) );
        final Sequence seq = new Sequence();
        seq.setDomainArchitecture( new DomainArchitecture( domains, 100 ) );
        leaf.getNodeData().setSequence( seq );
        root1.addAsChild( leaf );
        root1.addAsChild( new PhylogenyNode() );
        with_domains.setRoot( root1 );
        with_domains.externalNodesHaveChanged();
        if ( !AptxUtil.isHasAtLeastOneNodeWithDomainArchitecture( with_domains ) ) {
            return fail( "a tree with a domain architecture must be detected" );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [AptxUtilTest] " + message );
        return false;
    }

    private AptxUtilTest() {
        // not instantiable
    }
}
