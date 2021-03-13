
package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.SortedMap;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

public class pangolin_x {

    private static final String XSD_STRING            = "xsd:string";
    private static final String VIPR_PANGOLIN_CLADE   = "vipr:PANGO_Lineage";
    private static final String VIPR_PANGOLIN_CLADE_0 = "vipr:PANGO_Lineage_L0";
    private static final String VIPR_PANGOLIN_CLADE_1 = "vipr:PANGO_Lineage_L1";
    private static final String VIPR_PANGOLIN_CLADE_2 = "vipr:PANGO_Lineage_L2";
    private final static String PRG_NAME              = "pangolin_x";
    public static void main( final String args[] ) {
        if ( args.length != 3 ) {
            System.out.println( "\nWrong number of arguments, expected: lineage_file intree outtree\n" );
            System.exit( -1 );
        }
        final File lineage_file = new File( args[ 0 ] );
        final File infile = new File( args[ 1 ] );
        final File outfile = new File( args[ 2 ] );
        if ( !lineage_file.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + lineage_file + "] does not exist" );
        }
        if ( !infile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + infile + "] does not exist" );
        }
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "[" + outfile + "] already exists" );
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( infile, true );
            p = factory.create( infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + infile + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        BasicTable<String> mapping_table = null;
        try {
            mapping_table = BasicTableParser.parse( lineage_file, ',', false, false );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to read [" + lineage_file + "] [" + e.getMessage() + "]" );
        }
        final SortedMap<String, String> id_to_clade = mapping_table.getColumnsAsMap( 0, 1 );
        final List<PhylogenyNode> ext_nodes = p.getExternalNodes();
        for( final PhylogenyNode ext_node : ext_nodes ) {
            final String name = ext_node.getName();
            if ( !ForesterUtil.isEmpty( name ) ) {
                if ( id_to_clade.containsKey( name ) ) {
                    final String clade = id_to_clade.get( name );
                    if ( !ForesterUtil.isEmpty( clade ) ) {
                        PropertiesList custom_data = ext_node.getNodeData().getProperties();
                        final String clade_split[] = clade.split( "\\." );
                        final String clade_0 = clade_split.length > 0 ? clade_split[ 0 ] : null;
                        String clade_1 = null;
                        String clade_2 = null;
                        if ( clade_split.length > 1 ) {
                            clade_1 = clade_0 + "." + clade_split[ 1 ];
                        }
                        if ( clade_split.length > 2 ) {
                            clade_2 = clade_1 + "." + clade_split[ 2 ];
                        }
                        if ( custom_data == null ) {
                            custom_data = new PropertiesList();
                        }
                        custom_data.addProperty( new Property( VIPR_PANGOLIN_CLADE,
                                                               clade,
                                                               "",
                                                               XSD_STRING,
                                                               AppliesTo.NODE ) );
                        if ( !ForesterUtil.isEmpty( clade_0 ) ) {
                            custom_data.addProperty( new Property( VIPR_PANGOLIN_CLADE_0,
                                                                   clade_0,
                                                                   "",
                                                                   XSD_STRING,
                                                                   AppliesTo.NODE ) );
                        }
                        if ( !ForesterUtil.isEmpty( clade_1 ) ) {
                            custom_data.addProperty( new Property( VIPR_PANGOLIN_CLADE_1,
                                                                   clade_1,
                                                                   "",
                                                                   XSD_STRING,
                                                                   AppliesTo.NODE ) );
                        }
                        if ( !ForesterUtil.isEmpty( clade_2 ) ) {
                            custom_data.addProperty( new Property( VIPR_PANGOLIN_CLADE_2,
                                                                   clade_2,
                                                                   "",
                                                                   XSD_STRING,
                                                                   AppliesTo.NODE ) );
                        }
                        ext_node.getNodeData().setProperties( custom_data );
                    }
                }
            }
        }
        addInternalCladeInformation( p, VIPR_PANGOLIN_CLADE );
        addInternalCladeInformation( p, VIPR_PANGOLIN_CLADE_0 );
        addInternalCladeInformation( p, VIPR_PANGOLIN_CLADE_1 );
        addInternalCladeInformation( p, VIPR_PANGOLIN_CLADE_2 );
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outfile );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outfile + "]: " + e.getMessage() );
        }
        System.out.println( "[" + PRG_NAME + "] wrote: [" + outfile + "]" );
        System.out.println( "[" + PRG_NAME + "] OK" );
        System.out.println();
    }

    private static void addInternalCladeInformation( final Phylogeny p, final String property ) {
        for( final PhylogenyNodeIterator it = p.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isInternal() ) {
                final List<PhylogenyNode> childs = n.getDescendants();
                boolean same = true;
                String clade = null;
                for( final PhylogenyNode child : childs ) {
                    if ( child.isFirstChildNode() ) {
                        final PropertiesList custom_data = child.getNodeData().getProperties();
                        if ( custom_data != null ) {
                            final List<Property> clades = custom_data.getProperties( property );
                            if ( ( clades != null ) && ( clades.size() == 1 ) ) {
                                clade = clades.get( 0 ).getValue();
                            }
                            else {
                                same = false;
                            }
                        }
                        else {
                            same = false;
                        }
                    }
                    else if ( !ForesterUtil.isEmpty( clade ) ) {
                        final PropertiesList custom_data = child.getNodeData().getProperties();
                        if ( custom_data != null ) {
                            final List<Property> clades = custom_data.getProperties( property );
                            if ( ( clades != null ) && ( clades.size() == 1 ) ) {
                                final String my_clade = clades.get( 0 ).getValue();
                                if ( ForesterUtil.isEmpty( my_clade ) || !clade.equals( my_clade ) ) {
                                    same = false;
                                }
                            }
                            else {
                                same = false;
                            }
                        }
                        else {
                            same = false;
                        }
                    }
                }
                if ( same ) {
                    if ( !ForesterUtil.isEmpty( clade ) ) {
                        PropertiesList custom_data = n.getNodeData().getProperties();
                        if ( custom_data == null ) {
                            custom_data = new PropertiesList();
                        }
                        custom_data.addProperty( new Property( property, clade, "", XSD_STRING, AppliesTo.NODE ) );
                        n.getNodeData().setProperties( custom_data );
                    }
                }
            }
        }
    }
}
