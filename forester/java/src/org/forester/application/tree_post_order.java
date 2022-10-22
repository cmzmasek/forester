
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class tree_post_order {

    private static final String  PRG_DATE     = "2022-06-02";
    private static final String  PRG_VERSION  = "0.0.1";
    private static final String  PRG_NAME     = "tree_post_order";
    private static final Pattern annotation_p = Pattern.compile( "_\\{(.+)\\}" );
    //private static final Pattern annotation_p = Pattern.compile( "(.+?)\\|.+" );
    private static final String  XSD_STRING   = "xsd:string";
    private static final String  REF          = "subspecies:clade";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        System.out.println();
        if ( ( args.length != 2 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <out-tree> \n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File outtree = new File( args[ 1 ] );
        final String e0 = ForesterUtil.isWritableFile( outtree );
        if ( !ForesterUtil.isEmpty( e0 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e0 );
        }
        final String e1 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( e1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e1 );
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "Could not read \"" + intree + "\" [" + e.getMessage() + "]" );
        }
        if ( !p.isRooted() ) {
            ForesterUtil.fatalError( PRG_NAME, "\"" + intree + "\" is not rooted" );
        }
        p.setRerootable( false );
        ForesterUtil
                .programMessage( PRG_NAME,
                                 "Successfully read in tree with " + p.getNumberOfExternalNodes() + " external nodes" );
        final int properties_added = 0;
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                if ( node.isHasNodeData() && ( node.getNodeData().getProperties() != null )
                        && ( node.getNodeData().getProperties().size() > 0 )
                        && ( PhylogenyMethods.getNodePropertyValues( node, REF ).size() == 1 ) ) {
                }
                else if ( !ForesterUtil.isEmpty( node.getName() ) ) {
                    final String name = node.getName();
                    String annotation = null;
                    final Matcher m = annotation_p.matcher( name );
                    if ( m.find() ) {
                        annotation = m.group( 1 );
                    }
                    else {
                        //
                        if ( name.length() < 14 ) {
                            annotation = name;
                        }
                        else {
                            ForesterUtil.fatalError( PRG_NAME, "No annotation found in " + name );
                        }
                        //
                        //ForesterUtil.fatalError( PRG_NAME, "No annotation found in " + name );
                    }
                    final Property prop = new Property( REF, annotation, "", XSD_STRING, Property.AppliesTo.NODE );
                    PropertiesList custom_data = node.getNodeData().getProperties();
                    if ( custom_data == null ) {
                        custom_data = new PropertiesList();
                    }
                    custom_data.addProperty( prop );
                    node.getNodeData().setProperties( custom_data );
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "No annotation found in node" + node.getId() );
                }
            }
            else {
                final List<PhylogenyNode> descs = node.getDescendants();
                final List<String> annotatons = new ArrayList<>();
                for( final PhylogenyNode desc : descs ) {
                    if ( desc.isHasNodeData() && ( desc.getNodeData().getProperties() != null )
                            && ( desc.getNodeData().getProperties().size() > 0 )
                            && ( PhylogenyMethods.getNodePropertyValues( desc, REF ).size() == 1 ) ) {
                        annotatons.add( PhylogenyMethods.getNodePropertyValues( desc, REF ).get( 0 ) );
                    }
                    else {
                        ForesterUtil.fatalError( PRG_NAME, "No annotation found in node " + node.getId() );
                    }
                }
                //   final String x = ForesterUtil.greatestCommonPrefix( annotatons, "." );
                final String x = ForesterUtil.greatestCommonPrefix( annotatons );
                final Property prop = new Property( REF, x, "", XSD_STRING, Property.AppliesTo.NODE );
                PropertiesList custom_data = node.getNodeData().getProperties();
                if ( custom_data == null ) {
                    custom_data = new PropertiesList();
                }
                custom_data.addProperty( prop );
                node.getNodeData().setProperties( custom_data );
                node.setName( x );
            }
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outtree );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outtree + "]: " + e.getLocalizedMessage() );
        }
        ///
        final SortedMap<String, Integer> counts = new TreeMap<>();
        for( final PhylogenyNodeIterator it = p.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getParent() != null ) {
                String parent_annotation = "";
                final PropertiesList custom_data_p = n.getParent().getNodeData().getProperties();
                if ( custom_data_p != null ) {
                    final List<Property> clades = custom_data_p.getProperties( REF );
                    if ( ( clades != null ) && ( clades.size() == 1 ) ) {
                        parent_annotation = clades.get( 0 ).getValue();
                    }
                }
                final PropertiesList custom_data = n.getNodeData().getProperties();
                if ( custom_data != null ) {
                    final List<Property> clades = custom_data.getProperties( REF );
                    if ( ( clades != null ) && ( clades.size() == 1 ) ) {
                        final String my_annotation = clades.get( 0 ).getValue();
                        if ( !parent_annotation.equals( my_annotation ) ) {
                            if ( counts.containsKey( my_annotation ) ) {
                                counts.put( my_annotation, counts.get( my_annotation ) + 1 );
                            }
                            else {
                                counts.put( my_annotation, 1 );
                            }
                        }
                    }
                }
            }
        }
        int count_monop = 0;
        int count_polyp = 0;
        int count_total = 0;
        System.out.println();
        System.out.println();
        for( final Map.Entry<String, Integer> entry : counts.entrySet() ) {
            ++count_total;
            final String key = entry.getKey();
            final Integer value = entry.getValue();
            if ( value > 1 ) {
                ++count_polyp;
                System.out.println( key + "\t" + value );
            }
            else {
                ++count_monop;
            }
        }
        System.out.println();
        System.out.println( "Mono :\t" + count_monop );
        System.out.println( "Poly :\t" + count_polyp );
        System.out.println( "Total:\t" + count_total );
        System.out.println();
        ////
        System.out.println();
        System.out.println( "Sum of properties added: " + properties_added );
        System.out.println( "Wrote outtree to       : " + outtree );
        System.out.println();
    }
}
