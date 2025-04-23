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
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;
import org.forester.util.ViralUtils;

import java.io.File;
import java.io.IOException;
import java.sql.Date;
import java.time.LocalDate;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class vipr_x6 {
    public final static Pattern PATTERN_YEAR = Pattern.compile("([12]\\d{3})");

    public final static Pattern PATTERN2_YEAR = Pattern.compile("[/\\-_]([12]\\d{3}),");

    public final static Pattern PATTERN3_YEAR = Pattern.compile("\\|([12]\\d{3})");

    public final static Pattern PATTERN4_YEAR = Pattern.compile("_([12]\\d{3})_");

    public final static Pattern GENOTYPE_PATTERN = Pattern.compile("/\\[([A-Z0-9]+)\\]");
    public final static Pattern GENOTYPE_PATTERN_2 = Pattern.compile("enotype\\s+([A-Z0-9]+)");
    public final static Pattern GENOTYPE_PATTERN_3 = Pattern.compile("ubgroup\\s+([A-Z0-9]+)");
    private static final boolean VERBOSE = true;
    private static final String PRG_DATE = "2025-04-23";
    private static final String PRG_VERSION = "1.0.4";

    private static final int COL_GENOME_ID = 0;
    private static final int COL_STRAIN = 15;
    private static final int COL_SEGMENT = 20;
    private static final int COL_SUBTYPE = 21;
    private static final int COL_LINEAGE = 28;
    private static final int COL_CLADE2 = 29;
    private static final int COL_SUBCLADE = 30;
    private static final int COL_GENBANK_ACC = 43;
    private static final int COL_ISOLATION_SOURCE = 65;
    private static final int COL_COLLECTION_DATE = 67;


    private static final int COL_ISOLATION_COUNTRY = 70;
    private static final int COL_HOST_NAME_SCI = 74;

    private static final String UNKNOWN = "unknown";
    private static final String XSD_STRING = "xsd:string";
    private static final String HOST = "vipr:Host";
    private static final String COUNTRY = "vipr:Country";
    private static final String STRAIN = "vipr:Strain_Number";
    private static final String YEAR = "vipr:Year";
    private static final String SUBTYPE = "vipr:Subtype";
    private static final String SUBCLADE = "vipr:Subclade";
    private static final String REGION = "vipr:Region";
    private static final String STATE = "vipr:State";
    private static final String GB_ACCESSTION = "vipr:GB_Accession";
    private static final String HOST_GROUP = "vipr:Host_Group";
    private static final String HOST_GROUP_DOMESTIC_WILD = "vipr:Host_Group_Domestic_vs_Wild";
    private static final String HOST_SCI = "vipr:Host_Scientific";
    private static final String YEAR_MONTH = "vipr:Year_Month";
    private static final String ISOLATION_SOURCE = "vipr:Isolation_Source";
    private static final String BV_BRC_ACC = "vipr:BVBRC_Accession";

    private final static String PRG_NAME = "vipr_x6";


    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 3 && args.length != 4 && args.length != 5) {
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
        } else if (args.length == 4) {
            tree_name = args[3].trim();
            reroot = null;
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
                System.out.println("Strain: " + m_strain);
                String m_segment = map.getValue(COL_SEGMENT, r).replaceAll("\"", "");
                System.out.println("Segment: " + m_segment);

                String m_subtype = map.getValue(COL_SUBTYPE, r).replaceAll("\"", "");
                System.out.println("Subtype: " + m_subtype);

                String m_lineage = map.getValue(COL_LINEAGE, r).replaceAll("\"", "");
                System.out.println("Lineage: " + m_lineage);

                String m_clade = map.getValue(COL_CLADE2, r).replaceAll("\"", "");
                System.out.println("Clade: " + m_clade);

                String m_subclade = map.getValue(COL_SUBCLADE, r).replaceAll("\"", "");
                System.out.println("Subclade: " + m_subclade);

                String m_genbank = map.getValue(COL_GENBANK_ACC, r).replaceAll("\"", "");
                System.out.println("Genbank: " + m_genbank);
                String m_isolationsource = map.getValue(COL_ISOLATION_SOURCE, r).replaceAll("\"", "");
                System.out.println("Source: " + m_isolationsource);
                String m_coldate = map.getValue(COL_COLLECTION_DATE, r).replaceAll("\"", "");
                System.out.println("Date: " + m_coldate);

                String m_isolation_country = map.getValue(COL_ISOLATION_COUNTRY, r).replaceAll("\"", "");
                System.out.println("Country: " + m_isolation_country);

                String m_hostname_sci = map.getValue(COL_HOST_NAME_SCI, r).replaceAll("\"", "");
                System.out.println("Host: " + m_hostname_sci);
                System.out.println("---");
            }
            System.out.println();
        }

        final List<PhylogenyNode> ext_nodes = p.getExternalNodes();
        for (final PhylogenyNode ext_node : ext_nodes) {
            final String name = ext_node.getName();
            if (!ForesterUtil.isEmpty(name)) {
                final Matcher mg = ViralUtils.PATTERN_GB.matcher(name);
                final Matcher mbv = ViralUtils.PATTERN_BVBRC_ACC.matcher(name);
                String type = null;
                String strain_number = null;

                String subtype = null;
                String genotype = null;
                String genbank_acc = "";
                String bv_acc = "";
                if (mg.find()) {
                    genbank_acc = mg.group(1).trim();
                }
                if (mbv.find()) {
                    bv_acc = mbv.group(1).trim();
                }

                final Matcher mgt = GENOTYPE_PATTERN.matcher(name);
                if (mgt.find()) {
                    genotype = mgt.group(1).trim();
                } else {
                    final Matcher mgt2 = GENOTYPE_PATTERN_2.matcher(name);
                    if (mgt2.find()) {
                        genotype = mgt2.group(1).trim();
                    } else {
                        final Matcher mgt3 = GENOTYPE_PATTERN_3.matcher(name);
                        if (mgt3.find()) {
                            genotype = mgt3.group(1).trim();
                        }
                    }
                }

                String m_hostname_sci = null;
                String m_coldate = null;
                String m_isolationsource = null;
                String m_bvbrc_acc = null;
                String m_strain = null;
                String m_segment = null;
                String m_subtype = null;
                String m_lineage = null;
                String m_clade = null;
                String m_subclade = null;
                String m_isolation_country = null;
                String m_genbank = null;

                if (genbank_acc.length() > 3) {
                    for (int r = 0; r < map.getNumberOfRows(); ++r) {
                        if (map.getValue(COL_GENBANK_ACC, r).replaceAll("\"", "").equals(genbank_acc)) {
                            m_isolationsource = map.getValue(COL_ISOLATION_SOURCE, r).replaceAll("\"", "");
                            m_coldate = map.getValue(COL_COLLECTION_DATE, r).replaceAll("\"", "");
                            m_hostname_sci = map.getValue(COL_HOST_NAME_SCI, r).replaceAll("\"", "");
                            m_bvbrc_acc = map.getValue(COL_GENOME_ID, r).replaceAll("\"", "");
                            m_strain = map.getValue(COL_STRAIN, r).replaceAll("\"", "");
                            m_segment = map.getValue(COL_SEGMENT, r).replaceAll("\"", "");
                            m_subtype = map.getValue(COL_SUBTYPE, r).replaceAll("\"", "");
                            m_lineage = map.getValue(COL_LINEAGE, r).replaceAll("\"", "");
                            m_clade = map.getValue(COL_CLADE2, r).replaceAll("\"", "");
                            m_subclade = map.getValue(COL_SUBCLADE, r).replaceAll("\"", "");
                            m_genbank = map.getValue(COL_GENBANK_ACC, r).replaceAll("\"", "");
                            m_isolation_country = map.getValue(COL_ISOLATION_COUNTRY, r).replaceAll("\"", "");
                            break;
                        }
                    }
                }

                if (ForesterUtil.isEmpty(m_isolation_country)) {
                    m_isolation_country = obtainCountryFromName(name);
                }

                if (!ForesterUtil.isEmpty(m_isolation_country)) {
                    m_isolation_country = ViralUtils.determineCountry(m_isolation_country);
                }

                if (ForesterUtil.isEmpty(m_hostname_sci)) {
                    m_hostname_sci = obtainHostFromName(name);
                }


                if (VERBOSE) {
                    System.out.println();
                    System.out.println();
                    System.out.println("Name    : " + name);
                    System.out.println("Strain  : " + m_strain);
                    System.out.println("Segment : " + m_segment);
                    System.out.println("Subtype : " + m_subtype);
                    System.out.println("Lineage : " + m_lineage);
                    System.out.println("Clade   : " + m_clade);
                    System.out.println("Genotype: " + genotype);
                    System.out.println("Subclade: " + m_subclade);
                    System.out.println("Country : " + m_isolation_country);
                    System.out.println("Genbank : " + m_genbank);
                    System.out.println("Acc     : " + genbank_acc);
                    System.out.println("Type    : " + type);
                    System.out.println("Host    : " + m_hostname_sci);
                    System.out.println("Number  : " + strain_number);
                    System.out.println("Date    : " + m_coldate);
                    System.out.println("Subtype : " + subtype);
                }

                String year = null;
                if (m_coldate != null) {
                    final Matcher my = PATTERN_YEAR.matcher(m_coldate);
                    if (my.find()) {
                        year = my.group(1);
                    }
                }

                if (ForesterUtil.isEmpty(year)) {
                    year = obtainYearFromName(name);
                }

                if (!ForesterUtil.isEmpty(year)) {
                    checkYear(year);
                }

                System.out.println("Year    : " + year);

                if (ForesterUtil.isEmpty(m_isolation_country)) {
                    System.out.println("WARNING: No country for " + name);
                }
                if (ForesterUtil.isEmpty(m_hostname_sci)) {
                    System.out.println("WARNING: No host for " + name);
                }
                if (ForesterUtil.isEmpty(year)) {
                    System.out.println("WARNING: No year for " + name);
                }

                final PropertiesList custom_data = new PropertiesList();
                if (genbank_acc.length() > 5) {
                    custom_data.addProperty(new Property(GB_ACCESSTION, genbank_acc, "", XSD_STRING, AppliesTo.NODE));
                }
                if (!ForesterUtil.isEmpty(m_isolation_country)) {
                    custom_data.addProperty(new Property(COUNTRY, m_isolation_country, "", XSD_STRING, AppliesTo.NODE));
                }
                if (!ForesterUtil.isEmpty(year)) {
                    custom_data.addProperty(new Property(YEAR, year, "", XSD_STRING, AppliesTo.NODE));
                }
                if (!ForesterUtil.isEmpty(strain_number)) {
                    custom_data.addProperty(new Property(STRAIN, strain_number, "", XSD_STRING, AppliesTo.NODE));
                }
                if (!ForesterUtil.isEmpty(subtype)) {
                    custom_data.addProperty(new Property(SUBTYPE, subtype, "", XSD_STRING, AppliesTo.NODE));
                }
                if (!ForesterUtil.isEmpty(genotype)) {
                    custom_data.addProperty(new Property(SUBCLADE, genotype, "", XSD_STRING, AppliesTo.NODE));
                }
                if (!ForesterUtil.isEmpty(m_hostname_sci)) {
                    m_hostname_sci = m_hostname_sci.replace('_', ' ');
                    if (m_hostname_sci.equalsIgnoreCase("human")) {
                        m_hostname_sci = "Homo sapiens";
                    }
                    custom_data.addProperty(new Property(HOST, m_hostname_sci, "", XSD_STRING, AppliesTo.NODE));
                }
                String year_month = null;
                if (!ForesterUtil.isEmpty(m_coldate) && !ForesterUtil.isEmpty(year)) {
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
                ext_node.getNodeData().setSequence(null); //TODO need to be option
                ext_node.getNodeData().setProperties(custom_data);

                if (!ForesterUtil.isEmpty(m_isolation_country)) {
                    ViralUtils.addRegion(m_isolation_country, custom_data, REGION);
                }

                String new_name = cleanName(name);

                boolean make_new_name = true;

                if (make_new_name && !ForesterUtil.isEmpty( genbank_acc )) {
                    new_name = genbank_acc;
                    if (!ForesterUtil.isEmpty( m_strain )) {
                        new_name = new_name + "|" + m_strain;
                    }
                }
                if (!ForesterUtil.isEmpty( new_name)) {
                    ext_node.setName(new_name);
                }
                System.out.println("New Name: " + new_name);
                System.out.println("---");
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


    private static final String obtainHostFromName(final String name) {
        final String name_lc = name.toLowerCase();
        String host = null;
        if (name_lc.indexOf("sapiens") > 0 || name_lc.indexOf("human") > 0) {
            host = "Human";
        } else if (name_lc.indexOf("pan") > 0 || name_lc.indexOf("troglodytes") > 0) {
            host = "Chimpanzee";
        }
        return host;
    }

    private static final String obtainYearFromName(final String name) {
        final Matcher m = PATTERN2_YEAR.matcher(name);
        if (m.find()) {
            return m.group(1);
        }
        final Matcher m2 = PATTERN3_YEAR.matcher(name);
        if (m2.find()) {
            return m2.group(1);
        }
        final Matcher m3 = PATTERN4_YEAR.matcher(name);
        if (m3.find()) {
            return m3.group(1);
        }
        return null;
    }


    private static final String obtainCountryFromName(final String name) {
        final String name_lc = name.toLowerCase();
        String loc = null;
        if (name_lc.indexOf("nigeria") > 0) {
            loc = "Nigeria";
        } else if (name_lc.indexOf("germany") > 0) {
            loc = "Germany";
        } else if (name_lc.indexOf("/chn") > 0) {
            loc = "China";
        } else if (name_lc.indexOf("ivoire") > 0) {
            loc = "Cote d'Ivoire";
        } else if (name_lc.indexOf("democratic republic of the congo") > 0 || name.indexOf("DRC") > 0 || name.indexOf("RDC") > 0) {
            loc = "Democratic Republic of the Congo";
        } else if (name_lc.indexOf("central african rep") > 0) {
            loc = "Central African Republic";
        } else if (name_lc.indexOf("cameroon") > 0) {
            loc = "Cameroon";
        } else if (name_lc.indexOf("thailand") > 0) {
            loc = "Thailand";
        } else if (name_lc.indexOf("liberia") > 0) {
            loc = "Liberia";
        } else if (name_lc.indexOf("sierra leone") > 0) {
            loc = "Sierra Leone";
        } else if (name_lc.indexOf("congo") > 0) {
            loc = "Congo";
        } else if (name.indexOf("USA") > 0) {
            loc = "USA";
        } else if (name.indexOf("CAN") > 0) {
            loc = "Canada";
        } else if (name_lc.indexOf("australia") > 0) {
            loc = "Australia";
        } else if (name_lc.indexOf("gabon") > 0) {
            loc = "Gabon";
        } else if (name_lc.indexOf("rwanda") > 0) {
            loc = "Rwanda";
        } else if (name_lc.indexOf("nycphl") > 0) {
            loc = "USA";
        } else if (name_lc.indexOf("wa-uw") > 0) {
            loc = "USA";
        } else if (name_lc.indexOf("ca-lacphl") > 0) {
            loc = "USA";
        } else if (name.indexOf("ROK") > 0) {
            loc = "South Korea";
        }
        return loc;
    }

    private final static String cleanName(final String name) {
        String new_name = name;

        new_name = new_name.replaceAll("\\|\\|", "|");

        new_name = new_name.replaceAll("\\s+\\|", "|");

        new_name = new_name.replaceAll("\\s+Monkeypox virus", "|Monkeypox virus");

        /*if (new_name.indexOf(", partial") > 1) {
            new_name = new_name.substring(0, new_name.indexOf(", partial"));
        }
        if (new_name.indexOf(", complete") > 1) {
            new_name = new_name.substring(0, new_name.indexOf(", complete"));
        }*/
        if (new_name.endsWith("|")) {
            new_name = new_name.substring(0, new_name.length() - 1);
        }


        if (new_name.startsWith("accn|")) {
            new_name = new_name.substring(5);
        }
        return new_name;
    }

    public static void checkYear(final String year) {
        int year_int = Integer.parseInt(year);
        if (year_int < 1900 || year_int > 2025) {
            System.out.println("Error year \"" + year + "\" is out of range");
            System.exit(-1);
        }
    }
}
