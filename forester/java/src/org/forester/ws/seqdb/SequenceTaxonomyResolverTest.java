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

package org.forester.ws.seqdb;

import java.io.IOException;
import java.util.Arrays;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Unit tests for {@link SequenceTaxonomyResolver}, driven by in-memory fake clients (no network). They
 * pin the node-writing (sequence + two-phase taxonomy), the not-found bookkeeping, cancellation, and
 * clean abort on a transport error.
 */
public final class SequenceTaxonomyResolverTest {

    public static void main( final String[] args ) {
        System.out.println( "SequenceTaxonomyResolver: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        // node "P99999" is a valid UniProt accession -> resolvable; "leaf_two" is not
        final SequenceTaxonomyResolver r = new SequenceTaxonomyResolver( new FakeSeqFetcher(), new FakeTaxResolver() );
        final Phylogeny tree = twoLeafTree();
        final SequenceTaxonomyResolver.Result res = r.resolve( tree, CancelFlag.NEVER );
        if ( res.wasCancelled() || ( res.getError() != null ) ) {
            return fail( "a normal run must not be cancelled or errored: " + res.getError() );
        }
        if ( ( res.getSequencesWritten() != 1 ) || ( res.getTaxonomiesWritten() < 1 ) ) {
            return fail( "expected 1 sequence + >=1 taxonomy written; got " + res.getSequencesWritten() + "/"
                    + res.getTaxonomiesWritten() );
        }
        if ( ( res.getNotFound().size() != 1 ) || !res.getNotFound().contains( "leaf_two" ) ) {
            return fail( "the unresolvable leaf must be reported not-found; got " + res.getNotFound() );
        }
        final PhylogenyNode n1 = tree.getNode( "P99999" );

        // sequence written
        final Sequence seq = n1.getNodeData().getSequence();
        if ( ( seq == null ) || !"Cytochrome c".equals( seq.getName() ) || !"CYCS".equals( seq.getGeneName() )
                || !"MGDVEK".equals( seq.getMolecularSequence() ) || !"protein".equals( seq.getType() ) ) {
            return fail( "sequence fields not written: " + seq );
        }
        if ( ( seq.getAccession() == null ) || !"P99999".equals( seq.getAccession().getValue() ) ) {
            return fail( "sequence accession not written" );
        }
        if ( seq.getAnnotations().isEmpty() || seq.getCrossReferences().isEmpty() ) {
            return fail( "GO annotation / cross-reference not written" );
        }

        // taxonomy: organism from the entry, then refined by the detail phase (rank + common name)
        final Taxonomy tax = n1.getNodeData().getTaxonomy();
        if ( ( tax == null ) || !"Homo sapiens".equals( tax.getScientificName() ) || !"species".equals( tax.getRank() )
                || !"human".equals( tax.getCommonName() ) ) {
            return fail( "taxonomy detail not written: " + tax );
        }
        if ( ( tax.getIdentifier() == null ) || !"9606".equals( tax.getIdentifier().getValue() )
                || !"ncbi".equals( tax.getIdentifier().getProvider() ) ) {
            return fail( "taxonomy identifier (NCBI tax-id) not written: " + tax.getIdentifier() );
        }

        // source-feature metadata written as a node property (the viral/strain path)
        if ( !n1.getNodeData().isHasProperties()
                || n1.getNodeData().getProperties().getPropertiesWithGivenRef( "source:country" ).isEmpty()
                || !"USA".equals( n1.getNodeData().getProperties().getPropertiesWithGivenRef( "source:country" ).get( 0 )
                        .getValue() ) ) {
            return fail( "source:country property not written" );
        }

        // cancellation: a flag that is already set stops the walk before any node is written
        final Phylogeny fresh = twoLeafTree();
        final SequenceTaxonomyResolver.Result cancelled = r.resolve( fresh, new CancelFlag() {

            @Override
            public boolean isCancelled() {
                return true;
            }
        } );
        if ( !cancelled.wasCancelled() || ( cancelled.getSequencesWritten() != 0 ) ) {
            return fail( "an already-cancelled run must write nothing" );
        }

        // transport error: the fetcher throws -> the walk aborts with an error, not a hang
        final SequenceTaxonomyResolver erroring = new SequenceTaxonomyResolver( new SequenceFetcher() {

            @Override
            public SequenceEntry fetch( final Accession acc ) throws IOException {
                throw new IOException( "connection reset" );
            }
        }, new FakeTaxResolver() );
        final SequenceTaxonomyResolver.Result err = erroring.resolve( twoLeafTree(), CancelFlag.NEVER );
        if ( err.getError() == null ) {
            return fail( "a transport failure must surface as an error, aborting the walk" );
        }
        return true;
    }

    private static Phylogeny twoLeafTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "P99999" ); // a valid UniProt accession -> obtainAccessorFromDataFields picks it up
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "leaf_two" ); // not an accession, unknown to the fake taxonomy DB
        root.addAsChild( a );
        root.addAsChild( b );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static final class FakeSeqFetcher implements SequenceFetcher {

        @Override
        public SequenceEntry fetch( final Accession acc ) {
            if ( ( acc != null ) && "P99999".equals( acc.getValue() ) ) {
                final java.util.Map<String, String> props = new java.util.LinkedHashMap<String, String>();
                props.put( "source:country", "USA" ); // a GenBank source-feature qualifier (viral metadata path)
                return new SequenceEntry( "P99999", "CYC_HUMAN", "Cytochrome c", "CYCS", "MGDVEK",
                                          SequenceEntry.MoleculeType.PROTEIN, "Homo sapiens", "9606",
                                          Arrays.asList( "Eukaryota", "Homo sapiens" ), Arrays.asList( "GO:0005739" ),
                                          Arrays.asList( "refseq:NP_061820" ), props );
            }
            return SequenceEntry.EMPTY;
        }
    }

    private static final class FakeTaxResolver implements TaxonomyResolver {

        @Override
        public ResolvedTaxonomy resolveTaxonomy( final String query ) {
            // the resolver now prefers the entry's authoritative tax-id ("9606") over the name
            if ( "9606".equals( query ) || "Homo sapiens".equals( query ) ) {
                return new ResolvedTaxonomy( "Homo sapiens", "species", "9606", "human",
                                             Arrays.asList( "Eukaryota", "Metazoa", "Homo sapiens" ) );
            }
            return ResolvedTaxonomy.EMPTY;
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  SequenceTaxonomyResolver test failed: " + msg );
        return false;
    }

    private SequenceTaxonomyResolverTest() {
    }
}
