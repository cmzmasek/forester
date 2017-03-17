
package org.forester.surfacing;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.protein.Domain;
import org.forester.protein.Protein;
import org.forester.species.BasicSpecies;
import org.forester.species.Species;

public final class MinimalDomainomeCalculator {

    public static void main( final String[] args ) {
        Set<String> a = new HashSet<String>();
        Set<String> b = new HashSet<String>();
        Set<String> c = new HashSet<String>();
        Set<String> d = new HashSet<String>();
        a.add( "x" );
        a.add( "b" );
        a.add( "c" );
        b.add( "a" );
        b.add( "b" );
        b.add( "c" );
        c.add( "a" );
        c.add( "b" );
        c.add( "c" );
        c.add( "c" );
        c.add( "f" );
        d.add( "a" );
        d.add( "c" );
        d.add( "d" );
        List<Set<String>> domains_per_genome_list = new ArrayList();
        domains_per_genome_list.add( a );
        domains_per_genome_list.add( b );
        domains_per_genome_list.add( c );
        domains_per_genome_list.add( d );
        Set<String> x = x( domains_per_genome_list );
        System.out.println( x );
    }

    static final public void calc( Phylogeny tre, SortedMap<Species, List<Protein>> protein_lists_per_species ) {
        for( final PhylogenyNodeIterator iter = tre.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isInternal() ) {
                System.out.println();
                if ( node.getNodeData().isHasTaxonomy() ) {
                    System.out.println( node.getNodeData().getTaxonomy().getScientificName() + ":" );
                }
                else {
                    System.out.println( node.getName() + ":" );
                }
                final List<PhylogenyNode> e = node.getAllExternalDescendants();
                final List<Set<String>> domains_per_genome_list = new ArrayList();
                for( PhylogenyNode en : e ) {
                    final String code = en.getNodeData().getTaxonomy().getTaxonomyCode();
                    System.out.print( code + " " );
                    //System.out.println( protein_lists_per_species );
                    final List<Protein> x = protein_lists_per_species.get( new BasicSpecies( code ) );
                    if ( x != null ) {
                        final Set<String> d = new HashSet<String>();
                        for( Protein protein : x ) {
                            List<Domain> domains = protein.getProteinDomains();
                            for( Domain domain : domains ) {
                                d.add( domain.getDomainId() );
                            }
                        }
                        domains_per_genome_list.add( d );
                    }
                }
                System.out.println();
                Set<String> x = x( domains_per_genome_list );
                System.out.println( x );
            }
        }
    }

    static final Set<String> x( List<Set<String>> domains_per_genome_list ) {
        Set<String> first = domains_per_genome_list.get( 0 );
        for( int i = 1; i < domains_per_genome_list.size(); ++i ) {
            first.retainAll( domains_per_genome_list.get( i ) );
        }
        return first;
    }
}
