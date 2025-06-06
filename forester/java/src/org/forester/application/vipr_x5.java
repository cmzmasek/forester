
package org.forester.application;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sequence.MolecularSequence;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;
import org.forester.util.ViralUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class vipr_x5 {

    private static final boolean VERBOSE = false;
    private static final boolean DIE_ON_ERROR = false;
    private static final String PRG_DATE = "2025-04-06";
    private static final String PRG_VERSION = "1.0.2";

    private static final int COL_GENOME_ID = 0;
    private static final int COL_STRAIN = 15;
    private static final int COL_SEGMENT = 20;
    private static final int COL_SUBTYPE = 21;
    private static final int COL_GENBANK_ACC = 43;
    private static final int COL_ISOLATION_SOURCE = 65;
    private static final int COL_COLLECTION_DATE = 67;
    private static final int COL_HOST_NAME_SCI = 74;

    private static final String UNKNOWN = "unknown";
    private static final String XSD_STRING = "xsd:string";
    private static final String HOST = "vipr:Host";
    private static final String COUNTRY = "vipr:Country";
    private static final String STRAIN = "vipr:Strain_Number";
    private static final String YEAR = "vipr:Year";
    private static final String SUBTYPE = "vipr:Subtype";
    private static final String REGION = "vipr:Region";
    private static final String STATE = "vipr:State";
    private static final String GB_ACCESSTION = "vipr:GB_Accession";
    private static final String HOST_GROUP = "vipr:Host_Group";
    private static final String HOST_GROUP_DOMESTIC_WILD = "vipr:Host_Group_Domestic_vs_Wild";
    private static final String HOST_SCI = "vipr:Host_Scientific";
    private static final String YEAR_MONTH = "vipr:Year_Month";
    private static final String ISOLATION_SOURCE = "vipr:Isolation_Source";
    private static final String BV_BRC_ACC = "vipr:BVBRC_Accession";

    private final static String PRG_NAME = "vipr_x5";


    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 3 && args.length != 5) {
            System.out.println("\nWrong number of arguments, expected: intree mapfile outtree [name] [reroot node query]\n");
            System.out.println("    Examples: s4.xml map.txt seg_4.xml");
            System.out.println("              s4.xml map.txt seg_4.xml \"segment 4 2024-02-05\" \"/2003\"\n");

            System.exit(-1);
        }
        final File infile = new File(args[0]);
        final File mapfile = new File(args[1]);
        final File outfile = new File(args[2]);

        final String tree_name;
        final String reroot;
        if (args.length == 5) {
            tree_name = args[3].trim();
            reroot = args[4].trim();
        } else {
            tree_name = null;
            reroot = null;
        }

        if (!infile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile + "] does not exist");
        }
        if (!mapfile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + mapfile + "] does not exist");
        }
        if (outfile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + outfile + "] already exists");
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType(infile, true);
            p = factory.create(infile, pp)[0];
        } catch (final Exception e) {
            System.out.println("\nCould not read \"" + infile + "\" [" + e.getMessage() + "]\n");
            System.exit(-1);
        }

        BasicTable<String> map = null;
        try {
            map = BasicTableParser.parse(mapfile, '\t');
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, e.getMessage());
        }

        if (VERBOSE) {
            for (int r = 0; r < map.getNumberOfRows(); ++r) {
                String m_strain = map.getValue(COL_STRAIN, r).replaceAll("\"", "");
                System.out.println(m_strain);
                String m_segment = map.getValue(COL_SEGMENT, r).replaceAll("\"", "");
                System.out.println(m_segment);
                String m_subtype = map.getValue(COL_SUBTYPE, r).replaceAll("\"", "");
                System.out.println(m_subtype);
                String m_genbank = map.getValue(COL_GENBANK_ACC, r).replaceAll("\"", "");
                System.out.println(m_genbank);
                String m_isolationsource = map.getValue(COL_ISOLATION_SOURCE, r).replaceAll("\"", "");
                System.out.println(m_isolationsource);
                String m_coldate = map.getValue(COL_COLLECTION_DATE, r).replaceAll("\"", "");
                System.out.println(m_coldate);
                String m_hostname_sci = map.getValue(COL_HOST_NAME_SCI, r).replaceAll("\"", "");
                System.out.println(m_hostname_sci);
                System.out.println("---");
            }
            System.out.println();
        }

        final List<PhylogenyNode> ext_nodes = p.getExternalNodes();
        for (final PhylogenyNode ext_node : ext_nodes) {
            final String name = ext_node.getName();
            if (!ForesterUtil.isEmpty(name)) {
                final Matcher mg = ViralUtils.PATTERN_GB.matcher(name);
                final Matcher m0 = ViralUtils.PATTERN_0.matcher(name);
                final Matcher m00 = ViralUtils.PATTERN_00.matcher(name);
                final Matcher m1 = ViralUtils.PATTERN_1.matcher(name);
                final Matcher m2 = ViralUtils.PATTERN_2.matcher(name);
                final Matcher m3 = ViralUtils.PATTERN_3.matcher(name);
                final Matcher m1b = ViralUtils.PATTERN_1b.matcher(name);
                final Matcher mbv = ViralUtils.PATTERN_BVBRC_ACC.matcher(name);
                final Matcher mgi = ViralUtils.PATTERN_G.matcher(name);
                String type = UNKNOWN;
                String host = UNKNOWN;
                String location = UNKNOWN;
                String strain_number = UNKNOWN;
                String year = UNKNOWN;
                String subtype = UNKNOWN;
                String genbank_acc = "";
                String bv_acc = "";
                if (mg.find()) {
                    genbank_acc = mg.group(1).trim();
                }
                if (mbv.find()) {
                    bv_acc = mbv.group(1).trim();
                }

                if (m0.find()) {
                    // 1. type
                    // 2. host
                    // 3. country/state
                    // 4. strain number
                    // 5. year
                    // 6. subtype
                    type = m0.group(1).trim();
                    host = m0.group(2).trim();
                    location = m0.group(3).trim();
                    strain_number = m0.group(4).trim();
                    year = m0.group(5).trim();
                    subtype = m0.group(6).trim().toUpperCase();
                } else if (m1.find()) {
                    // 1. type
                    // 2. host
                    // 3. country/state
                    // 4. number
                    // 5. year
                    // 6. subtype
                    type = m1.group(1).trim();
                    host = m1.group(2).trim();
                    location = m1.group(3).trim();
                    strain_number = m1.group(4).trim();
                    year = m1.group(5).trim();
                    subtype = m1.group(6).trim().toUpperCase();
                } else if (m1b.find()) {
                    // 1. type
                    // 2. host
                    // 3. country/state
                    // 4. number
                    // 5. year
                    // 6. subtype
                    // 7. acc
                    type = m1b.group(1).trim();
                    host = m1b.group(2).trim();
                    location = m1b.group(3).trim();
                    strain_number = m1b.group(4).trim();
                    year = m1b.group(5).trim();
                    subtype = m1b.group(6).trim().toUpperCase();
                    genbank_acc = m1b.group(7).trim();

                } else if (m2.find()) {
                    // 1. type
                    // 2. country/state
                    // 3. number
                    // 4. year
                    // 5. subtype
                    host = "human";
                    type = m2.group(1).trim();
                    location = m2.group(2).trim();
                    strain_number = m2.group(3).trim();
                    year = m2.group(4).trim();
                    subtype = m2.group(5).trim().toUpperCase();
                } else if (m3.find()) {
                    // 1. type
                    // 2. host
                    // 3. country/state
                    // 4. number
                    // 5. year
                    type = "A";
                    host = m3.group(2).trim();
                    location = m3.group(3).trim();
                    strain_number = m3.group(4).trim();
                    year = m3.group(5).trim();
                    subtype = "H5N1";
                } else if (m00.find()) {
                    // 1. type
                    // 2. host
                    // 3. country/state
                    // 4. strain number
                    // 5. year
                    // 6. subtype
                    type = m00.group(1).trim();
                    host = m00.group(2).trim();
                    location = m00.group(3).trim();
                    strain_number = m00.group(4).trim();
                    year = m00.group(5).trim();
                    subtype = m00.group(6).trim().toUpperCase();
                } else if (mgi.find()) {
                    type = "A";
                    host = "human";
                    location = mgi.group(2).trim();
                    strain_number = mgi.group(3).trim();
                    year = mgi.group(4).trim();
                    subtype = "H5N1";
                    ext_node.getNodeData().setSequence(null);
                } else {
                    System.out.println("ERROR: name \"" + name + "\" could not be matched");
                    if (DIE_ON_ERROR) {
                        System.exit(-1);
                    } else {
                        if (name.startsWith("accn|")) {
                            ext_node.setName(name.substring(5));
                            System.out.println("     => " + ext_node.getName());
                        }
                        continue;
                    }
                }

                if (name.indexOf("/PHL-") > 0) {
                    ext_node.getNodeData().setSequence(null);
                }

                host = host.replace('_', ' ').trim();
                location = location.replace('_', ' ').trim();

                host = ViralUtils.cleanHost(host);

                host = ViralUtils.cleanHostString(host);
                location = ViralUtils.cleanLocationString(location);

                final String country = ViralUtils.determineCountry(location);
                final String state = ViralUtils.determineState(location);

                String m_hostname_sci = null;
                String m_coldate = null;
                String m_isolationsource = null;
                String m_bvbrc_acc = null;

                if (bv_acc.length() > 5) {
                    for (int r = 0; r < map.getNumberOfRows(); ++r) {
                        if (map.getValue(COL_GENOME_ID, r).replaceAll("\"", "").equals(bv_acc)) {
                            m_isolationsource = map.getValue(COL_ISOLATION_SOURCE, r).replaceAll("\"", "");
                            m_coldate = map.getValue(COL_COLLECTION_DATE, r).replaceAll("\"", "");
                            m_hostname_sci = map.getValue(COL_HOST_NAME_SCI, r).replaceAll("\"", "");
                            m_bvbrc_acc = map.getValue(COL_GENOME_ID, r).replaceAll("\"", "");
                            break;
                        }
                    }
                }

                if (VERBOSE) {
                    System.out.println();
                    System.out.println();
                    System.out.println("Name    : " + name);
                    System.out.println("Acc     : " + genbank_acc);
                    System.out.println("Type    : " + type);
                    System.out.println("Host    : " + host);
                    System.out.println("Location: " + location);
                    System.out.println("Country : " + country);
                    System.out.println("State   : " + state);
                    System.out.println("Number  : " + strain_number);
                    System.out.println("Year    : " + year);
                    System.out.println("Subtype : " + subtype);
                    if (!ForesterUtil.isEmpty(m_hostname_sci)) {
                        System.out.println("Host    : " + m_hostname_sci);
                    }
                    if (!ForesterUtil.isEmpty(m_coldate)) {
                        System.out.println("Coldate : " + m_coldate);
                    }
                    if (!ForesterUtil.isEmpty(m_isolationsource)) {
                        System.out.println("Source  : " + m_isolationsource);
                    }
                }

                year = ViralUtils.checkYear(year);

                final PropertiesList custom_data = new PropertiesList();
                if (genbank_acc.length() > 5) {
                    custom_data.addProperty(new Property(GB_ACCESSTION, genbank_acc, "", XSD_STRING, AppliesTo.NODE));
                }
                custom_data.addProperty(new Property(HOST, host, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(COUNTRY, country, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(YEAR, year, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(STRAIN, strain_number, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(SUBTYPE, subtype, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(STATE, state, "", XSD_STRING, AppliesTo.NODE));

                if (!ForesterUtil.isEmpty(m_hostname_sci)) {
                    m_hostname_sci = m_hostname_sci.replace('_', ' ');
                    custom_data.addProperty(new Property(HOST_SCI, m_hostname_sci, "", XSD_STRING, AppliesTo.NODE));
                }
                String year_month = null;
                if (!ForesterUtil.isEmpty(m_coldate)) {
                    year_month = makeYearMonth(year, m_coldate);
                    if (!ForesterUtil.isEmpty(year_month)) {
                        custom_data.addProperty(new Property(YEAR_MONTH, year_month, "", XSD_STRING, AppliesTo.NODE));
                    }
                }
                if (!ForesterUtil.isEmpty(m_isolationsource)) {
                    m_isolationsource = cleanIsolationSource(m_isolationsource);
                    custom_data.addProperty(new Property(ISOLATION_SOURCE, m_isolationsource, "", XSD_STRING, AppliesTo.NODE));
                }
                if (!ForesterUtil.isEmpty(m_bvbrc_acc)) {
                    custom_data.addProperty(new Property(BV_BRC_ACC, m_bvbrc_acc, "", XSD_STRING, AppliesTo.NODE));
                }

                ViralUtils.addRegion(country, custom_data, REGION);
                ViralUtils.addHostGroup(host, custom_data, HOST_GROUP, HOST_GROUP_DOMESTIC_WILD);

                ext_node.getNodeData().setProperties(custom_data);

                String new_name = null;

                if (!ForesterUtil.isEmpty(year_month)) {
                    year = year_month;
                }

                if (!ForesterUtil.isEmpty(m_isolationsource)) {
                    new_name = type + "/" + host + "/" + m_isolationsource + "/" + location + "/" + strain_number + "/" + year + "|" + subtype;
                } else {
                    new_name = type + "/" + host + "/" + location + "/" + strain_number + "/" + year + "|" + subtype;
                }

                if (genbank_acc.length() > 5) {
                    new_name += "|" + genbank_acc;
                }

                ext_node.setName(new_name);
            }
        }

        if (!ForesterUtil.isEmpty(tree_name)) {
            p.setName(tree_name);
        }
        if (!ForesterUtil.isEmpty(reroot)) {
            reRoot(reroot, p);
        }

        PhylogenyMethods.orderAppearance(p.getRoot(), true, true, PhylogenyMethods.DESCENDANT_SORT_PRIORITY.NODE_NAME);

        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML(p, 0, outfile);
        } catch (final IOException e) {
            ForesterUtil.fatalError(PRG_NAME, "failed to write to [" + outfile + "]: " + e.getMessage());
        }
        System.out.println();
        System.out.println();
        System.out.println("[" + PRG_NAME + "] wrote: [" + outfile + "]");
        System.out.println("[" + PRG_NAME + "] OK");
        System.out.println();
    }

    private static String cleanIsolationSource(String m_isolationsource) {
        m_isolationsource = m_isolationsource.replace('_', ' ').toLowerCase();
        if (m_isolationsource.equalsIgnoreCase("nasal swab")) {
            m_isolationsource = "Nasopharyngeal swab";
        }
        if (m_isolationsource.equalsIgnoreCase("brain")) {
            m_isolationsource = "Brain tissue";
        }
        m_isolationsource = m_isolationsource.substring(0, 1).toUpperCase() + m_isolationsource.substring(1);
        return m_isolationsource;
    }

    private static String makeYearMonth(String year, String m_coldate) {
        m_coldate = m_coldate.toLowerCase();
        String year_month = year + "_";
        if (m_coldate.indexOf("jan") > -1) {
            year_month += "01";
        } else if (m_coldate.indexOf("feb") > -1) {
            year_month += "02";
        } else if (m_coldate.indexOf("mar") > -1) {
            year_month += "03";
        } else if (m_coldate.indexOf("apr") > -1) {
            year_month += "04";
        } else if (m_coldate.indexOf("may") > -1) {
            year_month += "05";
        } else if (m_coldate.indexOf("jun") > -1) {
            year_month += "06";
        } else if (m_coldate.indexOf("jul") > -1) {
            year_month += "07";
        } else if (m_coldate.indexOf("aug") > -1) {
            year_month += "08";
        } else if (m_coldate.indexOf("sep") > -1) {
            year_month += "09";
        } else if (m_coldate.indexOf("oct") > -1) {
            year_month += "10";
        } else if (m_coldate.indexOf("nov") > -1) {
            year_month += "11";
        } else if (m_coldate.indexOf("dec") > -1) {
            year_month += "12";
        } else {
            year_month = null;
        }
        return year_month;
    }

    private static void reRoot(final String reroot, final Phylogeny p) {
        final Pattern pt = Pattern.compile(reroot);
        final List<PhylogenyNode> rnodes = p.getNodes(pt);
        if (rnodes.size() == 2) {
            final PhylogenyNode lca = PhylogenyMethods.calculateLCA(rnodes.get(0), rnodes.get(1));
            p.reRoot(lca);
        } else if (rnodes.size() == 1) {
            p.reRoot(rnodes.get(0));
        } else {
            ForesterUtil.fatalError(PRG_NAME, "cannot re-root with " + reroot);
        }
    }


}
