package org.forester.util;

import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;

import java.util.regex.Pattern;

public final class ViralUtils {


    // 1. type
    // 2. host
    // 3. country/state
    // 4. number
    // 5. year
    // 6. subtype
    //  (A/Duck/Champasak/261/2022(H5N1))
    public final static Pattern PATTERN_0 = Pattern
            .compile("\\((.*?)/([A-za-z-\\s'()]*?)/([A-za-z-\\s]*?)/(.*?)/(\\d{2,4}?)\\s*\\((.*?)\\)\\)");

    // 1. type
    // 2. host
    // 3. country/state
    // 4. number
    // 5. year
    // 6. subtype

    //  |A/red-tailed_hawk/California/24-003714-001/2024|H5N1|
    public final static Pattern PATTERN_1 = Pattern
            .compile("\\|(.*?)/([A-za-z-\\s'()]*?)/([A-za-z-\\s]*?)/(.*?)/(\\d{2,4}?)\\s*\\|(.*?)\\|");


    // 1. type
    // 2. host
    // 3. country/state
    // 4. number
    // 5. year
    // 6. subtype
    // 7. acc
    //  A/red-tailed_hawk/California/24-003714-001/2024|H5N1|acc
    public final static Pattern PATTERN_1b = Pattern
            .compile("(.*?)/([A-za-z-\\s'()]*?)/([A-za-z-\\s]*?)/(.*?)/(\\d{2,4}?)\\s*\\|(.*?)\\|(.+)");


    // 1. type
    // 2. country/state
    // 3. number
    // 4. year
    // 5. subtype
    //AF084262.1_Influenza_A_virus_(A/HongKong/483/97(H5N1))_
    public final static Pattern PATTERN_2 = Pattern
            .compile("\\((.*?)/(.*?)/(.*?)/(\\d{2,4}?)\\s*\\((.*?)\\)\\)");


    public final static Pattern PATTERN_GB = Pattern
            .compile("^(?:accn)?\\|?([A-Z][A-Z0-9.]+?)[_\\s|]");


    private static final String XSD_STRING = "xsd:string";

    private static final String UNKNOWN = "unknown";

    public static String cleanHost(final String host) {
        if (host.equalsIgnoreCase("duck")
                || host.equalsIgnoreCase("dk")
                || host.equalsIgnoreCase("mallard(anas platyrhynchos)")
                || host.equalsIgnoreCase("mallard duck")) {
            return "Duck";
        }
        if (host.equalsIgnoreCase("chicken")
                || host.equalsIgnoreCase("ck")) {
            return "Chicken";
        }
        if (host.equalsIgnoreCase("civet cat")) {
            return "Civet";
        }
        if (host.equalsIgnoreCase("domestic cat")) {
            return "Cat";
        }
        if (host.equalsIgnoreCase("canadian goose")) {
            return "Canada Goose";
        }
        if (host.equalsIgnoreCase("dairy cattle")
                ||host.equalsIgnoreCase("bovine")
                ||host.equalsIgnoreCase("cow")) {
            return "Cattle";
        }
        if (host.equalsIgnoreCase("cygnus cygnus")
                || host.equalsIgnoreCase("common swan")) {
            return "whooper swan";
        }
        return host;
    }

    public static String determineCountry(final String country) {

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
        } else if (c.equals("hokkaido")
                || c.equals("kyoto")
                || c.equals("miyagi")
                || c.equals("tottori")
                || c.equals("akita")
                || c.equals("chiba")) {
            return "Japan";
        } else if (c.equals("ontario") || c.equals("victoria")) {
            return "Canada";
        } else if (c.equals("greenland")) {
            return "Denmark";
        } else if (c.equals("nanchang")
                || c.equals("guiyang")
                || c.equals("shantou")
                || c.equals("jilin")
                || c.equals("wuhan")
                || c.equals("xinjiang")
                || c.equals("hunan")
                || c.equals("jiangsu")
                || c.equals("hubei")
                || c.equals("shandong")
                || c.equals("yunnan")
                || c.equals("sichuan")
                || c.equals("fujian")
                || c.equals("guangxi")
                || c.equals("tonghai")
                || c.equals("qinghai")
                || c.equals("jiangxi")
                || c.equals("zhejiang")
                || c.equals("hebei")
                || c.equals("huadong")
                || c.equals("ningxia")
                || c.equals("anhui")
                || c.equals("st")
                || c.equals("karakol lake")
                || c.equals("sheny")
                || c.equals("eastern china")
                || c.equals("beijing")) {
            return "China";
        } else if (c.equals("champasak") || c.equals("lao")) {
            return "Laos";
        } else if (c.equals("tyva") || c.equals("omsk")) {
            return "Russia";
        } else if (c.equals("arica y parinacota")) {
            return "Chile";
        } else if (c.equals("korea") ) {
            return "South Korea";
        } else if (c.equals("hongkong") || c.equals("hk")) {
            return "Hong Kong";
        } else if (c.equals("viet nam")
                || c.equals("vietnam hau giang")
        ) {
            return "Vietnam";
        } else if (c.equals("sidenreng rappang") || c.equals("sleman") || c.equals("klaten")
                || c.equals("majalengka")
                || c.equals("lamongan")
                || c.equals("westjava")
                || c.equals("east java")
                || c.equals("banten")
                || c.equals("pekalongan")
                || c.equals("banyuwangi")
                || c.equals("denpasar")) {
            return "Indonesia";
        } else if (c.equals("sagaing") || c.equals("yangon")) {
            return "Myanmar";
        } else if (c.equals("west bengal") || c.equals("sikkim")) {
            return "India";
        } else if (c.equals("gharbia") || c.equals("giza")) {
            return "Egypt";
        } else if (c.equals("england") || c.equals("scotland")) {
            return "United_Kingdom";
        } else if (c.equals("mangystau")) {
            return "Kazakhstan";
        }
        else if (c.equals("united states")) {
            return "USA";
        }
        else {
            return country;
        }
    }

    public static String determineState(final String country) {

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
                || c.equals("victoria")
                || c.equals("greenland")
        ) {
            return country;
        } else {
            return "";
        }
    }


    public static void addRegion(final String country,
                                 final PropertiesList custom_data,
                                 final String reg_ref) {
        String region = "";
        final String c = country.toLowerCase();
        if (c.equals("canada") || c.equals("usa") || c.equals("mexico")
        ) {
            region = "North America";
        } else if (c.equals("peru") || c.equals("ecuador") || c.equals("colombia") || c.equals("chile")
                || c.equals("brazil") || c.equals("argentina") || c.equals("guatemala") || c.equals("uruguay")
                || c.equals("venezuela")) {
            region = "South America";
        } else if (c.equals("france") || c.equals("germany")
                || c.equals("netherlands") || c.equals("united_kingdom")
                || c.equals("switzerland") || c.equals("austria") || c.equals("estonia")
                || c.equals("belgium")) {
            region = "Western Europe";
        } else if (c.equals("denmark") || c.equals("finland") | c.equals("norway") || c.equals("sweden")
                || c.equals("iceland")) {
            region = "Northern Europe";
        } else if (c.equals("serbia") || c.equals("greece") || c.equals("malta") || c.equals("italy")
                || c.equals("spain") || c.equals("portugal")) {
            region = "Southern Europe";
        } else if (c.equals("poland")) {
            region = "Central Europe";
        } else if (c.equals("russia") || c.equals("belarus")) {
            region = "Eastern Europe";
        } else if (c.equals("japan") || c.equals("taiwan") || c.equals("hong kong") || c.equals("south korea")
                || c.equals("tibet") || c.equals("china")) {
            region = "East Asia";
        } else if (c.equals("kazakhstan") || c.equals("uzbekistan") || c.equals("armenia")) {
            region = "Central Asia";
        } else if (c.equals("kuwait") || c.equals("jordan") || c.equals("bahrain") || c.equals("iraq")
                || c.equals("saudi arabia") || c.equals("turkey") || c.equals("egypt") || c.equals("israel")
                || c.equals("west bank") || c.equals("iran") || c.equals("lebanon")) {
            region = "West Asia";
        } else if (c.equals("india")
                || ((c.equals("pakistan") | c.equals("bangladesh")) || c.equals("sri lanka"))) {
            region = "South Asia";
        } else if (c.equals("laos") || c.equals("cambodia") || c.equals("thailand") || c.equals("malaysia")
                || c.equals("philippines") || c.equals("vietnam")
                || c.equals("myanmar") || c.equals("timor_leste") || c.equals("indonesia")) {
            region = "Southeast Asia";
        } else if (c.equals("mauritania") || c.equals("morocco") || c.equals("gambia") || c.equals("kenya") || c.equals("senegal")
                || c.equals("south africa") || c.equals("tanzania") || c.equals("ghana") || c.equals("benin")
                || c.equals("tunisia") || c.equals("nigeria") || c.equals("libya") || c.equals("djibouti")
                || c.equals("sierra leone") || c.equals("guinea") || c.equals("botswana") || c.equals("lesotho")
                || c.equals("ethiopia") || c.equals("namibia") || c.equals("malawi") || c.equals("mali") || c.equals("cameroon")) {
            region = "Africa";
        } else if (c.equals("australia") || c.equals("new_zealand")) {
            region = "Oceania";
        } else if (c.equals("dominican republic") || c.equals("puerto rico") || c.equals("jamaica")
                || c.equals("belize")) {
            region = "Caribbean";
        } else if (c.equals("na")) {
            region = "";
        } else {
            System.out.println("Error: unknown country \"" + c + "\"");
            System.exit(-1);
        }
        if (!ForesterUtil.isEmpty(region)) {
            custom_data.addProperty(new Property(reg_ref, region, "", XSD_STRING, Property.AppliesTo.NODE));
        }
    }


    public static String checkYear(final String year) {
        int year_int = Integer.parseInt(year);

        if (year_int >= 0) {
            if (year_int <= 24) {
                year_int += 2000;
            } else if (year_int <= 99) {
                year_int += 1900;
            }
        }

        if (year_int < 1800 || year_int > 2024) {
            System.out.println("Error year \"" + year + "\" is out of range");
            System.exit(-1);
        }
        return Integer.toString(year_int);
    }

    public static String cleanHostOrLocationString(final String s) {
        return s.replaceAll("_", " ").substring(0, 1).toUpperCase() + s.substring(1);
    }

    public static void addHostGroup(final String host, final PropertiesList custom_data,
                                    final String host_group_ref,
                                    final String host_group_dom_vs_wild_ref) {
        String hg1 = UNKNOWN;
        String hg2 = UNKNOWN;
        final String c = host.toLowerCase();
        if (c.length() == 0 || c.equalsIgnoreCase(UNKNOWN) || c.equalsIgnoreCase("na") || c.equals("?")) {
            hg1 = UNKNOWN;
            hg2 = UNKNOWN;
        } else if (c.equals("human")
                || c.equals("homo sapiens")) {
            hg1 = "Human";
            hg2 = "Human";
        } else if (c.equals("cattle")
                || c.equals("cow")
                || c.equals("bovine")
                || c.equals("equine")
                || c.equals("goat")
                || c.equals("sheep")
                || c.equals("swine")
                || c.equals("cat")
                || c.equals("feline")
                || c.equals("dog")
                || c.equals("canine")
                || c.equals("mink")
                || c.equals("ferret")
                || c.equals("stone marten")) {
            hg1 = "Non-Human Mammal";
            hg2 = "Non-Human Mammal (domestic)";
        } else if (c.equals("skunk")
                || c.equals("mountain lion")
                || c.equals("raccoon")
                || c.equals("pika")
                || c.equals("plateau pika")
                || c.equals("harbor seal")
                || c.equals("south american sea lion")
                || c.equals("tiger")
                || c.equals("owston's civet")
                || c.equals("civet")
        ) {
            hg1 = "Non-Human Mammal";
            hg2 = "Non-Human Mammal (wild)";

        } else if (c.equals("chicken")
                || c.equals("duck")
                || c.equals("mallard")
                || c.equals("turkey")
                || c.equals("goose")
                || c.equals("muscovy duck")
                || c.equals("poultry")) {
            hg1 = "Avian";
            hg2 = "Avian (domestic)";
        } else if (c.equals("openbill stork")
                || c.equals("pigeon")
                || c.equals("wild duck")
                || c.equals("peregrine falcon")
                || c.equals("common buzzard")
                || c.equals("bald eagle")
                || c.equals("bar-headed goose")
                || c.equals("whooper swan")
                || c.equals("quail")
                || c.equals("great crested grebe")
                || c.equals("crow")
                || c.equals("red-tailed hawk")
                || c.equals("blue-winged teal")
                || c.equals("american green-winged teal")
                || c.equals("cormorant")
                || c.equals("jungle crow")
                || c.equals("great black-backed gull")
                || c.equals("eurasian eagle owl")
                || c.equals("great horned owl")
                || c.equals("chukar")
                || c.equals("ruddy turnstone")
                || c.equals("teal")
                || c.equals("pacific black duck")
                || c.equals("gadwall duck")
                || c.equals("grackle")
                || c.equals("mute swan")
                || c.equals("canada goose")
                || c.equals("western gull")
                || c.equals("snow goose")
                || c.equals("common raven")
                || c.equals("harris-hawk")
                || c.equals("turkey vulture")
                || c.equals("american crow")
                || c.equals("redhead duck")
                || c.equals("american white pelican")
                || c.equals("wood duck")
                || c.equals("hooded merganser")
                || c.equals("blackbird")
        ) {
            hg1 = "Avian";
            hg2 = "Avian (wild)";
        } else if (c.equals("environment")) {
            hg1 = "Environment";
            hg2 = "Environment";
        } else {
            System.out.println("Error: Unknown host \"" + host + "\"");
            System.exit(-1);
        }

        if (!ForesterUtil.isEmpty(hg1)) {
            custom_data.addProperty(new Property(host_group_ref, hg1, "", XSD_STRING, Property.AppliesTo.NODE));
        }
        if (!ForesterUtil.isEmpty(hg2)) {
            custom_data.addProperty(new Property(host_group_dom_vs_wild_ref, hg2, "", XSD_STRING, Property.AppliesTo.NODE));
        }
    }
}
