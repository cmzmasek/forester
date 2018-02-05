
package org.forester.tools;

import org.forester.msa.Msa;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.sequence.MolecularSequence;
import org.forester.util.ForesterUtil;

public final class MutationInference {

    public final static void inferBranchEvents( final Phylogeny p, final String property_ref ) {
        final PhylogenyNodeIterator it = p.iteratorPreorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() ) {
                continue;
            }
            final PhylogenyNode parent = n.getParent();
            if ( n.getNodeData().isHasSequence() && parent.getNodeData().isHasSequence() ) {
                final String mol_seq_child = n.getNodeData().getSequence().getMolecularSequence();
                final String mol_seq_paren = parent.getNodeData().getSequence().getMolecularSequence();
                if ( !ForesterUtil.isEmpty( mol_seq_child ) && !ForesterUtil.isEmpty( mol_seq_paren ) ) {
                    if ( mol_seq_child.length() != mol_seq_paren.length() ) {
                        throw new IllegalArgumentException( "sequences must be of equal length" );
                    }
                    final PropertiesList custom_data = new PropertiesList();
                    for( int i = 0; i < mol_seq_child.length(); ++i ) {
                        if ( mol_seq_child.charAt( i ) != mol_seq_paren.charAt( i ) ) {
                            final Property prop = new Property( property_ref,
                                                                ( i + 1 ) + ": " + mol_seq_paren.charAt( i ) + " -> "
                                                                        + mol_seq_child.charAt( i ),
                                                                "",
                                                                "xsd:string",
                                                                AppliesTo.PARENT_BRANCH );
                            custom_data.addProperty( prop );
                        }
                    }
                    n.getNodeData().setProperties( custom_data );
                }
                else {
                    final String x;
                    if ( !ForesterUtil.isEmpty( mol_seq_child ) ) {
                        x = n.getName();
                    }
                    else {
                        x = parent.getName();
                    }
                    System.out.println( "warning: node \"" + x + "\" has no sequence associated" );
                }
            }
            else {
                final String x;
                if ( !n.getNodeData().isHasSequence() ) {
                    x = n.getName();
                }
                else {
                    x = parent.getName();
                }
                System.out.println( "warning: node \"" + x + "\" has no sequence associated" );
            }
        }
    }

    public final static void addMolSeqs( final Msa msa, final Phylogeny p ) {
        final PhylogenyNodeIterator it = p.iteratorPreorder();
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            final String name = n.getName();
            final MolecularSequence seq = msa.getSequence( name );
            if ( seq != null ) {
                final Sequence s = new Sequence();
                s.setMolecularSequenceAligned( true );
                s.setMolecularSequence( seq.getMolecularSequenceAsString() );
                n.getNodeData().setSequence( s );
            }
            else {
                throw new IllegalArgumentException( "no sequence for node \"" + name + "\" found in MSA" );
            }
        }
    }
}
