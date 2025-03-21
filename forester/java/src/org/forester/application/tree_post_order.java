
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

    private static final String PRG_DATE = "2022-06-02";
    private static final String PRG_VERSION = "0.0.1";
    private static final String PRG_NAME = "tree_post_order";
    private static final Pattern annotation_p = Pattern.compile("_\\{(.+)\\}");
    //private static final Pattern annotation_p = Pattern.compile( "(.+?)\\|.+" );
    private static final String XSD_STRING = "xsd:string";
    private static final String REF = "subspecies:clade";

    private static final String UNKNOWN = "unknown";
    private static final String HOST = "vipr:Host";
    private static final String CLADE = "vipr:Clade";
    private static final String COUNTRY = "vipr:Country";
    private static final String STRAIN = "vipr:Strain_Number";
    private static final String YEAR = "vipr:Year";
    private static final String SUBTYPE = "vipr:Subtype";
    private static final String REGION = "vipr:Region";
    private static final String STATE = "vipr:State";

    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        System.out.println();
        if ((args.length != 2)) {
            System.out.println(PRG_NAME + ": Wrong number of arguments.\n");
            System.out.println("Usage: " + PRG_NAME + " <in-tree> <out-tree> \n");
            System.exit(-1);
        }
        final File intree = new File(args[0]);
        final File outtree = new File(args[1]);
        final String e0 = ForesterUtil.isWritableFile(outtree);
        if (!ForesterUtil.isEmpty(e0)) {
            ForesterUtil.fatalError(PRG_NAME, e0);
        }
        final String e1 = ForesterUtil.isReadableFile(intree);
        if (!ForesterUtil.isEmpty(e1)) {
            ForesterUtil.fatalError(PRG_NAME, e1);
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(intree, true);
            p = factory.create(intree, pp)[0];
        } catch (final Exception e) {
            ForesterUtil.fatalError(PRG_NAME, "Could not read \"" + intree + "\" [" + e.getMessage() + "]");
        }
        if (!p.isRooted()) {
            ForesterUtil.fatalError(PRG_NAME, "\"" + intree + "\" is not rooted");
        }
        p.setRerootable(false);
        ForesterUtil
                .programMessage(PRG_NAME,
                        "Successfully read in tree with " + p.getNumberOfExternalNodes() + " external nodes");
        final int properties_added = 0;
        for (final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            procNode2(node);
        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML(p, 0, outtree);
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, "failed to write to [" + outtree + "]: " + e.getLocalizedMessage());
        }
        ///
       /* final SortedMap<String, Integer> counts = new TreeMap<>();
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

        */
        System.out.println();
        System.out.println("Sum of properties added: " + properties_added);
        System.out.println("Wrote outtree to       : " + outtree);
        System.out.println();
    }


    private static void procNode2(final PhylogenyNode node) {

        // 1. host
        // 2. country/state
        // 3. number
        // 4. year

        // Name: accn|OQ331006 Monkeypox virus isolate Monkeypox virus/Human/USA/CA-LACPHL-MA00397/2023, partial genome
        final Pattern PATTERN_1 = Pattern
                .compile("([A-za-z-\\s]*?)/([A-za-z-0-9\\s]*?)/([A-za-z-0-9\\s]*?)/(\\d{2,4})");

        final Pattern PATTERN_YEAR = Pattern
                .compile("[|/](\\d{4})");


        if (node.isExternal()) {
            final String name = node.getName();
            //final Matcher m0 = PATTERN_1.matcher(name);

            String host = "";
            String loc = "";
            String number = "";
            String year = "";
            String clade = "";

           /* if (m0.find()) {
                host = m0.group(1);
                loc = m0.group(2);
                number = m0.group(3);
                year = m0.group(4);

            }
            else {*/
            final Matcher m1 = PATTERN_YEAR.matcher(name);
            if (m1.find()) {
                year = m1.group(1);
            }


            final String name_lc = name.toLowerCase();
            if (name_lc.indexOf("nigeria") > 0) {
                loc = "Nigeria";
                //clade = "II";
            } else if (name_lc.indexOf("germany") > 0) {
                loc = "Germany";
            } else if (name_lc.indexOf("/chn") > 0) {
                loc = "China";
            } else if (name_lc.indexOf("ivoire") > 0) {
                loc = "Cote d'Ivoire";
                //clade = "II";
            } else if (name_lc.indexOf("democratic republic of the congo") > 0
                    || name.indexOf("DRC") > 0) {
                loc = "Democratic Republic of the Congo";
                //clade = "I";
            } else if (name_lc.indexOf("central african rep") > 0) {
                loc = "Central African Republic";
                //clade = "I";
            } else if (name_lc.indexOf("cameroon") > 0) {
                loc = "Cameroon";
            } else if (name_lc.indexOf("thailand") > 0) {
                loc = "Thailand";
            } else if (name_lc.indexOf("liberia") > 0) {
                loc = "Liberia";
                //clade = "II";
            } else if (name_lc.indexOf("sierra leone") > 0) {
                loc = "Sierra Leone";
            } else if (name_lc.indexOf("congo") > 0) {
                loc = "Congo";
                //clade = "I";
            } else if (name.indexOf("USA") > 0) {
                loc = "USA";
            } else if (name.indexOf("CAN") > 0) {
                loc = "Canada";
            } else if (name_lc.indexOf("australia") > 0) {
                loc = "Australia";
                //clade = "I";
            } else if (name_lc.indexOf("gabon") > 0) {
                loc = "Gabon";
                //clade = "I";
            }
            //else {
            ///    System.out.println(name + ":");
            // }
            if (name_lc.indexOf("sapiens") > 0 || name_lc.indexOf("human") > 0) {
                host = "Human";
            } else if (name_lc.indexOf("pan") > 0 || name_lc.indexOf("troglodytes") > 0) {
                host = "Chimpanzee";
            }


            //}

            if (loc.equals("CAN")) {
                loc = "Canada";
            }

            String new_name = name;

            new_name = new_name.replaceAll("\\|\\|", "|");

            new_name = new_name.replaceAll("\\s+\\|", "|");

            new_name = new_name.replaceAll("\\s+Monkeypox virus", "|Monkeypox virus");

            if (new_name.indexOf(", partial") > 1) {
                new_name = new_name.substring(0, new_name.indexOf(", partial"));
            }
            if (new_name.indexOf(", complete") > 1) {
               new_name = new_name.substring(0, new_name.indexOf(", complete"));
            }
            if (new_name.endsWith("|")) {
                new_name = new_name.substring(0, new_name.length() - 1);
            }


            if (new_name.startsWith("accn|")) {
                new_name = new_name.substring(5);
            }


            PhylogenyNode n = node;
            while (!n.isRoot()) {
                n = n.getParent();
                if (n.isInternal() && !ForesterUtil.isEmpty(n.getName())) {
                    clade = n.getName();
                    break;
                }
                /*

                if ( n.getName().equals("II.b.C.1") ) {
                    clade = "II.b.C.1";
                    break;
                }
                if ( n.getName().equals("II.b.B.1") ) {
                    clade = "II.b.B.1";
                    break;
                }
                if ( n.getName().equals("II.b") ) {
                    clade = "II.b";
                    break;
                }
                if ( n.getName().equals("II.a") ) {
                    clade = "II.a";
                    break;
                }
                if ( n.getName().equals("II") ) {
                    clade = "II";
                    break;
                }
                if ( n.getName().equals("I") ) {
                    clade = "I";
                    break;
                }*/

            }


            System.out.println(name + ":");
            System.out.println(new_name + ":");
            System.out.println(clade);
            System.out.println(host);
            System.out.println(loc);
            System.out.println(number);
            System.out.println(year);
            //new_name = new_name.substring(0, new_name.lastIndexOf('|'));

            //node.setName(clade + "|" + new_name);
            node.setName(new_name + "|" + clade);


            final PropertiesList custom_data = new PropertiesList();

            custom_data.addProperty(new Property(HOST, host, "", XSD_STRING, Property.AppliesTo.NODE));
            custom_data.addProperty(new Property(COUNTRY, loc, "", XSD_STRING, Property.AppliesTo.NODE));
            custom_data.addProperty(new Property(YEAR, year, "", XSD_STRING, Property.AppliesTo.NODE));
            custom_data.addProperty(new Property(STRAIN, number, "", XSD_STRING, Property.AppliesTo.NODE));
            custom_data.addProperty(new Property(CLADE, clade, "", XSD_STRING, Property.AppliesTo.NODE));


            node.getNodeData().setProperties(custom_data);

        }
    }

    private static void procNode1(final PhylogenyNode node) {
        if (node.isExternal()) {
            if (node.isHasNodeData() && (node.getNodeData().getProperties() != null)
                    && (node.getNodeData().getProperties().size() > 0)
                    && (PhylogenyMethods.getNodePropertyValues(node, REF).size() == 1)) {
            } else if (!ForesterUtil.isEmpty(node.getName())) {
                final String name = node.getName();
                String annotation = null;
                final Matcher m = annotation_p.matcher(name);
                if (m.find()) {
                    annotation = m.group(1);
                } else {
                    //
                    if (name.length() < 14) {
                        annotation = name;
                    } else {
                        ForesterUtil.fatalError(PRG_NAME, "No annotation found in " + name);
                    }
                    //
                    //ForesterUtil.fatalError( PRG_NAME, "No annotation found in " + name );
                }
                final Property prop = new Property(REF, annotation, "", XSD_STRING, Property.AppliesTo.NODE);
                PropertiesList custom_data = node.getNodeData().getProperties();
                if (custom_data == null) {
                    custom_data = new PropertiesList();
                }
                custom_data.addProperty(prop);
                node.getNodeData().setProperties(custom_data);
            } else {
                ForesterUtil.fatalError(PRG_NAME, "No annotation found in node" + node.getId());
            }
        } else {
            final List<PhylogenyNode> descs = node.getDescendants();
            final List<String> annotatons = new ArrayList<>();
            for (final PhylogenyNode desc : descs) {
                if (desc.isHasNodeData() && (desc.getNodeData().getProperties() != null)
                        && (desc.getNodeData().getProperties().size() > 0)
                        && (PhylogenyMethods.getNodePropertyValues(desc, REF).size() == 1)) {
                    annotatons.add(PhylogenyMethods.getNodePropertyValues(desc, REF).get(0));
                } else {
                    ForesterUtil.fatalError(PRG_NAME, "No annotation found in node " + node.getId());
                }
            }
            //   final String x = ForesterUtil.greatestCommonPrefix( annotatons, "." );
            final String x = ForesterUtil.greatestCommonPrefix(annotatons);
            final Property prop = new Property(REF, x, "", XSD_STRING, Property.AppliesTo.NODE);
            PropertiesList custom_data = node.getNodeData().getProperties();
            if (custom_data == null) {
                custom_data = new PropertiesList();
            }
            custom_data.addProperty(prop);
            node.getNodeData().setProperties(custom_data);
            node.setName(x);
        }
    }
}
