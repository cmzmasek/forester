
package org.forester.surfacing;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.go.GoId;
import org.forester.protein.BinaryDomainCombination;
import org.forester.protein.BinaryDomainCombination.DomainCombinationType;
import org.forester.protein.Domain;
import org.forester.protein.Protein;
import org.forester.species.Species;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class BasicGenomeWideCombinableDomains implements GenomeWideCombinableDomains {

    private final static NumberFormat                  FORMATTER                                  = new DecimalFormat( "0.0E0" );
    private static final Comparator<CombinableDomains> DESCENDING_KEY_DOMAIN_COUNT_ORDER          = new Comparator<CombinableDomains>() {

                                                                                                      @Override
                                                                                                      public int compare( final CombinableDomains d1,
                                                                                                                          final CombinableDomains d2 ) {
                                                                                                          if ( d1.getKeyDomainCount() < d2
                                                                                                                  .getKeyDomainCount() ) {
                                                                                                              return 1;
                                                                                                          }
                                                                                                          else if ( d1
                                                                                                                  .getKeyDomainCount() > d2
                                                                                                                  .getKeyDomainCount() ) {
                                                                                                              return -1;
                                                                                                          }
                                                                                                          else {
                                                                                                              return d1
                                                                                                                      .getKeyDomain()
                                                                                                                      .compareTo( d2
                                                                                                                              .getKeyDomain() );
                                                                                                          }
                                                                                                      }
                                                                                                  };
    private static final Comparator<CombinableDomains> DESCENDING_KEY_DOMAIN_PROTEINS_COUNT_ORDER = new Comparator<CombinableDomains>() {

                                                                                                      @Override
                                                                                                      public int compare( final CombinableDomains d1,
                                                                                                                          final CombinableDomains d2 ) {
                                                                                                          if ( d1.getKeyDomainProteinsCount() < d2
                                                                                                                  .getKeyDomainProteinsCount() ) {
                                                                                                              return 1;
                                                                                                          }
                                                                                                          else if ( d1
                                                                                                                  .getKeyDomainProteinsCount() > d2
                                                                                                                  .getKeyDomainProteinsCount() ) {
                                                                                                              return -1;
                                                                                                          }
                                                                                                          else {
                                                                                                              return d1
                                                                                                                      .getKeyDomain()
                                                                                                                      .compareTo( d2
                                                                                                                              .getKeyDomain() );
                                                                                                          }
                                                                                                      }
                                                                                                  };
    private static final Comparator<CombinableDomains> DESCENDING_COMBINATIONS_COUNT_ORDER        = new Comparator<CombinableDomains>() {

                                                                                                      @Override
                                                                                                      public int compare( final CombinableDomains d1,
                                                                                                                          final CombinableDomains d2 ) {
                                                                                                          if ( d1.getNumberOfCombinableDomains() < d2
                                                                                                                  .getNumberOfCombinableDomains() ) {
                                                                                                              return 1;
                                                                                                          }
                                                                                                          else if ( d1
                                                                                                                  .getNumberOfCombinableDomains() > d2
                                                                                                                  .getNumberOfCombinableDomains() ) {
                                                                                                              return -1;
                                                                                                          }
                                                                                                          else {
                                                                                                              return d1
                                                                                                                      .getKeyDomain()
                                                                                                                      .compareTo( d2
                                                                                                                              .getKeyDomain() );
                                                                                                          }
                                                                                                      }
                                                                                                  };
    final private SortedMap<String, CombinableDomains> _combinable_domains_map;
    final private Species                              _species;
    final private DomainCombinationType                _dc_type;

    private BasicGenomeWideCombinableDomains( final Species species, final DomainCombinationType dc_type ) {
        _combinable_domains_map = new TreeMap<String, CombinableDomains>();
        _species = species;
        _dc_type = dc_type;
    }

    private void add( final String key, final CombinableDomains cdc ) {
        _combinable_domains_map.put( key, cdc );
    }

    @Override
    public boolean contains( final String key_id ) {
        return _combinable_domains_map.containsKey( key_id );
    }

    @Override
    public CombinableDomains get( final String key_id ) {
        return _combinable_domains_map.get( key_id );
    }

    @Override
    public SortedMap<String, CombinableDomains> getAllCombinableDomainsIds() {
        return _combinable_domains_map;
    }

    @Override
    public SortedSet<String> getAllDomainIds() {
        final SortedSet<String> domains = new TreeSet<String>();
        for( final String key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            final List<String> ds = cb.getAllDomains();
            for( final String d : ds ) {
                domains.add( d );
            }
        }
        return domains;
    }

    @Override
    public DomainCombinationType getDomainCombinationType() {
        return _dc_type;
    }

    @Override
    public SortedSet<String> getMostPromiscuosDomain() {
        final SortedSet<String> doms = new TreeSet<String>();
        final int max = ( int ) getPerGenomeDomainPromiscuityStatistics().getMax();
        for( final String key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            if ( cb.getNumberOfCombinableDomains() == max ) {
                doms.add( key );
            }
        }
        return doms;
    }

    @Override
    public DescriptiveStatistics getPerGenomeDomainPromiscuityStatistics() {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final String key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            stats.addValue( cb.getNumberOfCombinableDomains() );
        }
        return stats;
    }

    @Override
    public int getSize() {
        return _combinable_domains_map.size();
    }

    @Override
    public Species getSpecies() {
        return _species;
    }

    @Override
    public SortedSet<BinaryDomainCombination> toBinaryDomainCombinations() {
        final SortedSet<BinaryDomainCombination> binary_combinations = new TreeSet<BinaryDomainCombination>();
        for( final String key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            for( final BinaryDomainCombination b : cb.toBinaryDomainCombinations() ) {
                binary_combinations.add( b );
            }
        }
        return binary_combinations;
    }

    @Override
    public String toString() {
        return toStringBuilder( GenomeWideCombinableDomainsSortOrder.ALPHABETICAL_KEY_ID ).toString();
    }

    // Produces something like: 
    // 2-oxoacid_dh      5       5       2       4.8E-67   Biotin_lipoyl [4], E3_binding [3]
    @Override
    public StringBuilder toStringBuilder( final GenomeWideCombinableDomainsSortOrder sort_order ) {
        final StringBuilder sb = new StringBuilder();
        final List<CombinableDomains> combinable_domains = new ArrayList<CombinableDomains>();
        for( final String key : getAllCombinableDomainsIds().keySet() ) {
            final CombinableDomains cb = getAllCombinableDomainsIds().get( key );
            combinable_domains.add( cb );
        }
        if ( sort_order == GenomeWideCombinableDomainsSortOrder.KEY_DOMAIN_COUNT ) {
            Collections.sort( combinable_domains, BasicGenomeWideCombinableDomains.DESCENDING_KEY_DOMAIN_COUNT_ORDER );
        }
        else if ( sort_order == GenomeWideCombinableDomainsSortOrder.KEY_DOMAIN_PROTEINS_COUNT ) {
            Collections.sort( combinable_domains,
                              BasicGenomeWideCombinableDomains.DESCENDING_KEY_DOMAIN_PROTEINS_COUNT_ORDER );
        }
        else if ( sort_order == GenomeWideCombinableDomainsSortOrder.COMBINATIONS_COUNT ) {
            Collections.sort( combinable_domains, BasicGenomeWideCombinableDomains.DESCENDING_COMBINATIONS_COUNT_ORDER );
        }
        for( final CombinableDomains cb : combinable_domains ) {
            sb.append( ForesterUtil.pad( new StringBuffer( cb.getKeyDomain().toString() ), 18, ' ', false ) );
            sb.append( ForesterUtil.pad( new StringBuffer( "" + cb.getKeyDomainCount() ), 8, ' ', false ) );
            sb.append( ForesterUtil.pad( new StringBuffer( "" + cb.getKeyDomainProteinsCount() ), 8, ' ', false ) );
            sb.append( ForesterUtil.pad( new StringBuffer( "" + cb.getNumberOfCombinableDomains() ), 8, ' ', false ) );
            sb.append( ForesterUtil.pad( new StringBuffer( ""
                                                 + FORMATTER.format( cb.getKeyDomainConfidenceDescriptiveStatistics()
                                                         .median() ) ),
                                         10,
                                         ' ',
                                         false ) );
            sb.append( cb.getCombiningDomainIdsAsStringBuilder() );
            sb.append( ForesterUtil.getLineSeparator() );
        }
        return sb;
    }

    private static void countDomains( final Map<String, Integer> domain_counts,
                                      final Map<String, Integer> domain_protein_counts,
                                      final Map<String, DescriptiveStatistics> stats,
                                      final Set<String> saw_c,
                                      final String id_i,
                                      final double support ) {
        if ( domain_counts.containsKey( id_i ) ) {
            domain_counts.put( id_i, 1 + domain_counts.get( ( id_i ) ) );
            if ( !saw_c.contains( id_i ) ) {
                domain_protein_counts.put( id_i, 1 + domain_protein_counts.get( ( id_i ) ) );
            }
        }
        else {
            stats.put( id_i, new BasicDescriptiveStatistics() );
            domain_counts.put( id_i, 1 );
            domain_protein_counts.put( id_i, 1 );
        }
        stats.get( id_i ).addValue( support );
        saw_c.add( id_i );
    }

    public static BasicGenomeWideCombinableDomains createInstance( final List<Protein> protein_list,
                                                                   final boolean ignore_combination_with_same_domain,
                                                                   final Species species ) {
        return createInstance( protein_list,
                               ignore_combination_with_same_domain,
                               species,
                               null,
                               DomainCombinationType.BASIC,
                               null,
                               null );
    }

    public static BasicGenomeWideCombinableDomains createInstance( final List<Protein> protein_list,
                                                                   final boolean ignore_combination_with_same_domain,
                                                                   final Species species,
                                                                   final DomainCombinationType dc_type ) {
        return createInstance( protein_list, ignore_combination_with_same_domain, species, null, dc_type, null, null );
    }

    public static BasicGenomeWideCombinableDomains createInstance( final List<Protein> protein_list,
                                                                   final boolean ignore_combination_with_same_domain,
                                                                   final Species species,
                                                                   final Map<String, List<GoId>> domain_id_to_go_ids_map,
                                                                   final DomainCombinationType dc_type,
                                                                   final Map<String, DescriptiveStatistics> protein_length_stats_by_dc,
                                                                   final Map<String, DescriptiveStatistics> domain_number_stats_by_dc ) {
        final BasicGenomeWideCombinableDomains instance = new BasicGenomeWideCombinableDomains( species, dc_type );
        final Map<String, Integer> domain_counts = new HashMap<String, Integer>();
        final Map<String, Integer> domain_protein_counts = new HashMap<String, Integer>();
        final Map<String, DescriptiveStatistics> stats = new HashMap<String, DescriptiveStatistics>();
        for( final Protein protein : protein_list ) {
            if ( !protein.getSpecies().equals( species ) ) {
                throw new IllegalArgumentException( "species (" + protein.getSpecies()
                        + ") does not match species of combinable domains collection (" + species + ")" );
            }
            final Set<String> saw_i = new HashSet<String>();
            final Set<String> saw_c = new HashSet<String>();
            for( int i = 0; i < protein.getProteinDomains().size(); ++i ) {
                final Domain pd_i = protein.getProteinDomain( i );
                final String id_i = pd_i.getDomainId();
                final int current_start = pd_i.getFrom();
                BasicGenomeWideCombinableDomains.countDomains( domain_counts,
                                                               domain_protein_counts,
                                                               stats,
                                                               saw_c,
                                                               id_i,
                                                               pd_i.getPerSequenceEvalue() );
                if ( !saw_i.contains( id_i ) ) {
                    if ( dc_type == DomainCombinationType.BASIC ) {
                        saw_i.add( id_i );
                    }
                    CombinableDomains domain_combination = null;
                    if ( instance.contains( id_i ) ) {
                        domain_combination = instance.get( id_i );
                    }
                    else {
                        if ( dc_type == DomainCombinationType.DIRECTED_ADJACTANT ) {
                            domain_combination = new AdjactantDirectedCombinableDomains( pd_i.getDomainId(), species );
                        }
                        else if ( dc_type == DomainCombinationType.DIRECTED ) {
                            domain_combination = new DirectedCombinableDomains( pd_i.getDomainId(), species );
                        }
                        else {
                            domain_combination = new BasicCombinableDomains( pd_i.getDomainId(), species );
                        }
                        instance.add( id_i, domain_combination );
                    }
                    final Set<String> saw_j = new HashSet<String>();
                    if ( ignore_combination_with_same_domain ) {
                        saw_j.add( id_i );
                    }
                    Domain closest = null;
                    for( int j = 0; j < protein.getNumberOfProteinDomains(); ++j ) {
                        if ( ( dc_type != DomainCombinationType.BASIC )
                                && ( current_start >= protein.getProteinDomain( j ).getFrom() ) ) {
                            continue;
                        }
                        if ( i != j ) {
                            final String id = protein.getProteinDomain( j ).getDomainId();
                            if ( !saw_j.contains( id ) ) {
                                saw_j.add( id );
                                if ( dc_type != DomainCombinationType.DIRECTED_ADJACTANT ) {
                                    domain_combination
                                            .addCombinableDomain( protein.getProteinDomain( j ).getDomainId() );
                                }
                                else {
                                    if ( closest == null ) {
                                        closest = protein.getProteinDomain( j );
                                    }
                                    else {
                                        if ( protein.getProteinDomain( j ).getFrom() < closest.getFrom() ) {
                                            closest = protein.getProteinDomain( j );
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if ( ( dc_type == DomainCombinationType.DIRECTED_ADJACTANT ) && ( closest != null ) ) {
                        domain_combination.addCombinableDomain( closest.getDomainId() );
                    }
                    if ( protein_length_stats_by_dc != null ) {
                        final List<BinaryDomainCombination> dcs = domain_combination.toBinaryDomainCombinations();
                        for( final BinaryDomainCombination dc : dcs ) {
                            final String dc_str = dc.toString();
                            if ( !protein_length_stats_by_dc.containsKey( dc_str ) ) {
                                protein_length_stats_by_dc.put( dc_str, new BasicDescriptiveStatistics() );
                            }
                            protein_length_stats_by_dc.get( dc_str ).addValue( protein.getLength() );
                        }
                    }
                    if ( domain_number_stats_by_dc != null ) {
                        final List<BinaryDomainCombination> dcs = domain_combination.toBinaryDomainCombinations();
                        for( final BinaryDomainCombination dc : dcs ) {
                            final String dc_str = dc.toString();
                            if ( !domain_number_stats_by_dc.containsKey( dc_str ) ) {
                                domain_number_stats_by_dc.put( dc_str, new BasicDescriptiveStatistics() );
                            }
                            domain_number_stats_by_dc.get( dc_str ).addValue( protein.getNumberOfProteinDomains() );
                        }
                    }
                    //
                }
            }
        }
        for( final String key_id : domain_counts.keySet() ) {
            instance.get( key_id ).setKeyDomainCount( domain_counts.get( key_id ) );
            instance.get( key_id ).setKeyDomainProteinsCount( domain_protein_counts.get( key_id ) );
            instance.get( key_id ).setKeyDomainConfidenceDescriptiveStatistics( stats.get( key_id ) );
        }
        return instance;
    }
}
