
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

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class vipr_x4 {

    private static final boolean VERBOSE = true;
    private static final boolean DIE_ON_ERROR = false;

    private static final String UNKNOWN = "unknown";

    private static final String XSD_STRING = "xsd:string";
    private static final String HOST = "vipr:Host";
    private static final String COUNTRY = "vipr:Country";
    private static final String STRAIN = "vipr:Strain_Number";
    private static final String YEAR = "vipr:Year";
    private static final String SUBTYPE = "vipr:Subtype";
    private static final String REGION = "vipr:Region";
    private static final String STATE = "vipr:State";
    private static final String GB_ACCESSTION = "vipr:gb_accession";
    private static final String H5_CLADE = "vipr:h5_clade";
    private static final String HOST_RANGE = "vipr:Host_Range";
    private static final String HOST_RANGE_DOMESTIC_WILD = "vipr:Host_Range_Domestic_vs_Wild";
    private final static String PRG_NAME = "vipr_x4";
    private static final String PRG_DATE = "2024-05-08";
    private static final String PRG_VERSION = "1.0.1";


    // 1. type
    // 2. host
    // 3. country/state
    // 4. number
    // 5. year
    // 6. subtype
    //  (A/Duck/Champasak/261/2022(H5N1))
    private final static Pattern PATTERN_0 = Pattern
            .compile("\\((.*?)/(.*?)/(.*?)/(.*?)/(.*?)\\((.*?)\\)\\)");

    // 1. type
    // 2. host
    // 3. country/state
    // 4. number
    // 5. year
    // 6. subtype

    //  |A/red-tailed_hawk/California/24-003714-001/2024|H5N1|
    private final static Pattern PATTERN_1 = Pattern
            .compile("\\|(.*?)/(.*?)/(.*?)/(.*?)/(.*?)\\|(.*?)\\|");


    // 1. type
    // 2. country/state
    // 3. number
    // 4. year
    // 5. subtype
    //AF084262.1_Influenza_A_virus_(A/HongKong/483/97(H5N1))_
    private final static Pattern PATTERN_2 = Pattern
            .compile("\\((.*?)/(.*?)/(.*?)/(.*?)\\((.*?)\\)\\)");


    private final static Pattern PATTERN_GB = Pattern
            .compile("^\\|?([A-Z][A-Z0-9.]+?)[_\\s|]");



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
                final Matcher mg = PATTERN_GB.matcher(name);
                final Matcher m0 = PATTERN_0.matcher(name);
                final Matcher m1 = PATTERN_1.matcher(name);
                final Matcher m2 = PATTERN_2.matcher(name);
                String type =UNKNOWN;
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
                    subtype = m0.group(6).trim();
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
                    subtype = m1.group(6).trim();
                } else if (m2.find()) {
                    // 1. type
                    // 2. country/state
                    // 3. number
                    // 4. year
                    // 5. subtype
                    type = m2.group(1).trim();
                    location = m2.group(2).trim();
                    strain_number = m2.group(3).trim();
                    year = m2.group(4).trim();
                    subtype = m2.group(5).trim();
                } else {
                    System.out.println("----> WARNING name \"" + name + "\" could not be matched");
                    if (DIE_ON_ERROR) {
                        System.exit(-1);
                    } else {
                        continue;
                    }
                }
                if (year.equals("97")) {
                    year = "1997";
                }
                if (year.equals("3500/2022")) {
                    year = "2022";
                }
                if (location.contentEquals("hongkong")) {
                    location = "Hong Kong";
                }

                host = procString(host);
                location = procString(location);

                final String country = determineCountry(location);
                final String state = determineState(location);

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
                final PropertiesList custom_data = new PropertiesList();

                if (genbank_acc.length() > 5 ) {
                    custom_data.addProperty(new Property(GB_ACCESSTION, genbank_acc, "", XSD_STRING, AppliesTo.NODE));
                }
                custom_data.addProperty(new Property(HOST, host, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(COUNTRY, country, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(YEAR, year, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(STRAIN, strain_number, "", XSD_STRING, AppliesTo.NODE));
                //custom_data.addProperty(new Property(H5_CLADE, "", "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(SUBTYPE, subtype, "", XSD_STRING, AppliesTo.NODE));
                custom_data.addProperty(new Property(STATE, state, "", XSD_STRING, AppliesTo.NODE));
                
                ext_node.getNodeData().setProperties(custom_data);

                addRegion(country, custom_data);
                addHostRange(host, custom_data);
                addHostRangeWildVsDomestic(host, custom_data);

                String new_name = type + "/" + host + "/" + location + "/" + strain_number + "/" + year + "|" + subtype;

                if (genbank_acc.length() > 5 ) {
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

    private static String procString(String s) {
        s = s.replaceAll("_", " ");
        s = s.substring(0, 1).toUpperCase() + s.substring(1);
        return s;
    }

    private static void addHostRange(final String host, final PropertiesList custom_data) {
        String hr = "";
        final String c = host.toLowerCase();
        if ( c.length() == 0 || c.equalsIgnoreCase(UNKNOWN) || c.equalsIgnoreCase( "na") || c.equals("?")) {
            hr = UNKNOWN;
        }
        else if (c.equals("cattle") || c.equals("cow") || c.equals("bovine") || c.equals("cat") || c.equals("feline") || c.equals("goat")
                || c.equals("swine") || c.equals("skunk") || c.equals("mountain lion") || c.equals("raccoon")
        ) {
            hr = "Mammalian";

        } else {
            hr = "Avian";
        }
        if (VERBOSE) {
            System.out.println(host + " -> " + hr);
        }
        if (!ForesterUtil.isEmpty(hr)) {
            custom_data.addProperty(new Property(HOST_RANGE, hr, "", XSD_STRING, AppliesTo.NODE));
        }
    }

    private static void addHostRangeWildVsDomestic(final String host, final PropertiesList custom_data) {
        String hr = "";
        final String c = host.toLowerCase();


        if ( c.length() == 0 || c.equalsIgnoreCase(UNKNOWN) || c.equalsIgnoreCase( "na") || c.equals("?")) {
            hr = UNKNOWN;
        }
        else if (c.equals("cattle") || c.equals("cow") || c.equals("bovine") || c.equals("cat") || c.equals("feline") || c.equals("goat")
                || c.equals("swine") || c.equals("skunk") || c.equals("mountain lion") || c.equals("raccoon")
        ) {
            if (c.equals("cattle") || c.equals("cow") || c.equals("bovine") || c.equals("cat") || c.equals("feline") || c.equals("goat")
                    || c.equals("swine")) {
                hr = "Mammalian (Domestic)";
            } else {
                hr = "Mammalian (Wild)";
            }
        } else {
            if (c.equals("chicken") || c.equals("turkey") ) {
                hr = "Avian (Domestic)";
            } else {
                hr = "Avian (Wild)";
            }
        }
        if (VERBOSE) {
            System.out.println(host + " -> " + hr);
        }
        if (!ForesterUtil.isEmpty(hr)) {
            custom_data.addProperty(new Property(HOST_RANGE_DOMESTIC_WILD, hr, "", XSD_STRING, AppliesTo.NODE));
        }
    }


    private static String determineCountry(final String country) {

        final String c = country.toLowerCase();
        if (c.equals("alabama")
                || c.equals("alaska")
                || c.equals("arizona")
                || c.equals("arkansas")
                || c.equals("california")
                || c.equals("colorado")
                || c.equals("connecticut")
                || c.equals("delaware")
                || c.equals("florida")
                || c.equals("georgia")
                || c.equals("hawaii")
                || c.equals("idaho")
                || c.equals("illinois")
                || c.equals("indiana")
                || c.equals("iowa")
                || c.equals("kansas")
                || c.equals("kentucky")
                || c.equals("louisiana")
                || c.equals("maine")
                || c.equals("maryland")
                || c.equals("massachusetts")
                || c.equals("michigan")
                || c.equals("minnesota")
                || c.equals("mississippi")
                || c.equals("missouri")
                || c.equals("montana")
                || c.equals("nebraska")
                || c.equals("nevada")
                || c.equals("new hampshire")
                || c.equals("new jersey")
                || c.equals("new mexico")
                || c.equals("new york")
                || c.equals("north carolina")
                || c.equals("north dakota")
                || c.equals("ohio")
                || c.equals("oklahoma")
                || c.equals("oregon")
                || c.equals("pennsylvania")
                || c.equals("rhode island")
                || c.equals("south carolina")
                || c.equals("south dakota")
                || c.equals("tennessee")
                || c.equals("texas")
                || c.equals("utah")
                || c.equals("vermont")
                || c.equals("virginia")
                || c.equals("washington")
                || c.equals("west virginia")
                || c.equals("wisconsin")
                || c.equals("wyoming")
        ) {
            return "USA";
        } else if (c.equals("ontario")) {
            return "Canada";
        } else if (c.equals("greenland")) {
            return "Denmark";
        }
        else if (c.equals("nanchang")) {
            return "China";
        }
        else {
            return country;
        }
    }

    private static String determineState(final String country) {

        final String c = country.toLowerCase();
        if (c.equals("alabama")
                || c.equals("alaska")
                || c.equals("arizona")
                || c.equals("arkansas")
                || c.equals("california")
                || c.equals("colorado")
                || c.equals("connecticut")
                || c.equals("delaware")
                || c.equals("florida")
                || c.equals("georgia")
                || c.equals("hawaii")
                || c.equals("idaho")
                || c.equals("illinois")
                || c.equals("indiana")
                || c.equals("iowa")
                || c.equals("kansas")
                || c.equals("kentucky")
                || c.equals("louisiana")
                || c.equals("maine")
                || c.equals("maryland")
                || c.equals("massachusetts")
                || c.equals("michigan")
                || c.equals("minnesota")
                || c.equals("mississippi")
                || c.equals("missouri")
                || c.equals("montana")
                || c.equals("nebraska")
                || c.equals("nevada")
                || c.equals("new hampshire")
                || c.equals("new jersey")
                || c.equals("new mexico")
                || c.equals("new york")
                || c.equals("north carolina")
                || c.equals("north dakota")
                || c.equals("ohio")
                || c.equals("oklahoma")
                || c.equals("oregon")
                || c.equals("pennsylvania")
                || c.equals("rhode island")
                || c.equals("south carolina")
                || c.equals("south dakota")
                || c.equals("tennessee")
                || c.equals("texas")
                || c.equals("utah")
                || c.equals("vermont")
                || c.equals("virginia")
                || c.equals("washington")
                || c.equals("west virginia")
                || c.equals("wisconsin")
                || c.equals("wyoming")
                || c.equals("ontario")
                || c.equals("greenland")
        ) {
            return country;
        } else {
            return "na";
        }
    }


    private static void addRegion(final String country, final PropertiesList custom_data) {
        String region = "";
        final String c = country.toLowerCase();
        if (c.equals("canada") || c.equals("usa") || c.equals("mexico")
        ) {
            region = "North America";
        } else if (c.equals("peru") || c.equals("ecuador") || c.equals("colombia") || c.equals("chile")
                || c.equals("brazil") || c.equals("argentina") || c.equals("guatemala") || c.equals("uruguay")
                || c.equals("venezuela")) {
            region = "South America";
        } else if (c.equals("denmark") || c.equals("finland") || c.equals("france") || c.equals("germany")
                || c.equals("netherlands") || c.equals("norway") || c.equals("united_kingdom")
                || c.equals("switzerland") || c.equals("austria") || c.equals("estonia") || c.equals("sweden")
                || c.equals("belgium") || c.equals("scotland")) {
            region = "Western Europe";
        } else if (c.equals("serbia") || c.equals("greece") || c.equals("malta") || c.equals("italy")
                || c.equals("spain") || c.equals("portugal")) {
            region = "Southern Europe";
        } else if (c.equals("poland")) {
            region = "Central Europe";
        } else if (c.equals("russia") || c.equals("belarus")) {
            region = "Eastern Europe";
        } else if (c.equals("japan") || c.equals("taiwan") || c.equals("hong_kong") || c.equals("hongkong") || c.equals("hong kong") || c.equals("south_korea")
                || c.equals("china") || c.equals("nanchang")) {
            region = "East Asia";
        } else if (c.equals("kazakhstan") || c.equals("uzbekistan") || c.equals("armenia")) {
            region = "Central Asia";
        } else if (c.equals("kuwait") || c.equals("jordan") || c.equals("bahrain") || c.equals("iraq")
                || c.equals("saudi_arabia") || c.equals("turkey") || c.equals("egypt") || c.equals("israel")
                || c.equals("west_bank") || c.equals("iran")) {
            region = "West Asia";
        } else if (c.equals("india")
                || ((c.equals("pakistan") | c.equals("bangladesh")) || c.equals("sri_lanka"))) {
            region = "South Asia";
        } else if (c.equals("laos") || c.equals("champasak") || c.equals("markeev") || c.equals("cambodia") || c.equals("thailand") || c.equals("malaysia")
                || c.equals("philippines") || c.equals("viet_nam") || c.equals("viet nam") || c.equals("vietnam")
                || c.equals("myanmar") || c.equals("timor_leste")) {
            region = "Southeast Asia";
        } else if (c.equals("mauritania") || c.equals("morocco") || c.equals("gambia") || c.equals("kenya") || c.equals("senegal")
                || c.equals("south_africa") || c.equals("tanzania") || c.equals("ghana") || c.equals("benin")
                || c.equals("tunisia") || c.equals("nigeria") || c.equals("libya") || c.equals("djibouti")
                || c.equals("sierra_leone") || c.equals("guinea") || c.equals("botswana")
                || c.equals("ethiopia") || c.equals("malawi") || c.equals("mali")) {
            region = "Africa";
        } else if (c.equals("australia") || c.equals("new_zealand")) {
            region = "Oceania";
        } else if (c.equals("dominican_republic") || c.equals("puerto_rico") || c.equals("jamaica")
                || c.equals("belize")) {
            region = "Caribbean";
        } else if (c.equals("na")) {
            region = "";
        } else {
            System.out.println("ERROR: unknown country \"" + c + "\"");
            System.exit(-1);
        }
        if (!ForesterUtil.isEmpty(region)) {
            custom_data.addProperty(new Property(REGION, region, "", XSD_STRING, AppliesTo.NODE));
        }
    }
}
