
package org.forester.application;

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
import org.forester.util.ForesterUtil;
import org.forester.util.ViralUtils;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;

public class vipr_x4 {

    private static final boolean VERBOSE = true;
    private static final boolean DIE_ON_ERROR = true;
    private static final String PRG_DATE = "2024-05-17";
    private static final String PRG_VERSION = "1.1.1";
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
    private final static String PRG_NAME = "vipr_x4";


    public static void main(final String args[]) {
        ForesterUtil.printProgramInformation(PRG_NAME, PRG_VERSION, PRG_DATE);
        if (args.length != 2) {
            System.out.println("\nWrong number of arguments, expected: intree outtree\n");
            System.exit(-1);
        }
        final File infile = new File(args[0]);
        final File outfile = new File(args[1]);
        if (!infile.exists()) {
            ForesterUtil.fatalError(PRG_NAME, "[" + infile + "] does not exist");
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
        final List<PhylogenyNode> ext_nodes = p.getExternalNodes();
        for (final PhylogenyNode ext_node : ext_nodes) {
            final String name = ext_node.getName();
            if (!ForesterUtil.isEmpty(name)) {
                final Matcher mg = ViralUtils.PATTERN_GB.matcher(name);
                final Matcher m0 = ViralUtils.PATTERN_0.matcher(name);
                final Matcher m1 = ViralUtils.PATTERN_1.matcher(name);
                final Matcher m2 = ViralUtils.PATTERN_2.matcher(name);
                final Matcher m1b = ViralUtils.PATTERN_1b.matcher(name);
                String type = UNKNOWN;
                String host = UNKNOWN;
                String location = UNKNOWN;
                String strain_number = UNKNOWN;
                String year = UNKNOWN;
                String subtype = UNKNOWN;
                String genbank_acc = "";
                if (mg.find()) {
                    genbank_acc = mg.group(1).trim();
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
                } else {
                    System.out.println("ERROR: name \"" + name + "\" could not be matched");
                    if (DIE_ON_ERROR) {
                        System.exit(-1);
                    } else {
                        continue;
                    }
                }

                host = ViralUtils.cleanHost(host);

                host = ViralUtils.cleanHostOrLocationString(host);
                location = ViralUtils.cleanHostOrLocationString(location);


                final String country = ViralUtils.determineCountry(location);
                final String state = ViralUtils.determineState(location);

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

                ext_node.getNodeData().setProperties(custom_data);

                ViralUtils.addRegion(country, custom_data, REGION);
                ViralUtils.addHostGroup(host, custom_data, HOST_GROUP, HOST_GROUP_DOMESTIC_WILD);

                String new_name = type + "/" + host + "/" + location + "/" + strain_number + "/" + year + "|" + subtype;

                if (genbank_acc.length() > 5) {
                    new_name += "|" + genbank_acc;
                }

                ext_node.setName(new_name);
            }
        }
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


}
