package org.forester.util;

import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;

import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

    //A/Duck/Champasak/261/2022(H5N1)
    public final static Pattern PATTERN_00 = Pattern
            .compile("(.*?)/([A-za-z-\\s'()]*?)/([A-za-z-\\s]*?)/(.*?)/(\\d{2,4}?)\\s*\\((.*?)\\)");

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


    // 1. type
    // 2. host
    // 3. country/state
    // 4. number
    // 5. year
    //  [Influenza_A_virus_A/Domestic_Cat/USA/24-009038-001/2024|11320.588919]
    public final static Pattern PATTERN_3 = Pattern
            .compile("\\[(.*?)/([A-za-z-\\s'()]*?)/([A-za-z-\\s]*?)/(.*?)/(\\d{2,4}?)\\|");

    public final static Pattern PATTERN_GB = Pattern
            .compile("^(?:accn)?\\|?([A-Z][A-Z0-9.]+?)[_\\s|]");

    public final static Pattern PATTERN_BVBRC_ACC = Pattern
            .compile("\\|([0-9]+\\.[0-9]+)\\]$");


    // 1. type
    // 2. location
    // 3. number
    // 4. year
    //  |A/British_Columbia/PHL-2032/2024|2024-11-01|2.3.4.4b|PB1|A/British_Columbia/PHL-2032/2024_PB1|2|EPI3650020|A/British_Columbia/PHL-2032/2024|EPI_ISL_19548836
    public final static Pattern PATTERN_G = Pattern
            .compile("\\|([A-Z]+)/([A-Za-z_]+)/([A-Za-z_0-9-]+)/(\\d+)");


    private static final String XSD_STRING = "xsd:string";

    private static final String UNKNOWN = "unknown";

    public static String cleanHost(final String host) {
        final String h = host.toLowerCase();
        if (h.equals("duck")
                || h.equals("dk")
                || h.indexOf("anas platyrhynchos") >= 0
                || h.equals("mallard duck")) {
            return "Duck";
        } else if (h.equals("chicken")
                || h.equals("ck")) {
            return "Chicken";
        } else if (h.equals("civet cat")) {
            return "Civet";
        } else if (h.equals("domestic cat")) {
            return "Cat";
        } else if (h.equals("canadian goose") || h.equals("ganada goose")) {
            return "Canada goose";
        } else if (h.equals("dairy cattle")
                || h.equals("bovine")
                || h.equals("cattle milk product")
                || h.equals("dairy cow")
                || h.equals("cow")) {
            return "Cattle";
        } else if (h.indexOf("cygnus cygnus") >= 0
                || h.equals("common swan")) {
            return "Whooper swan";
        } else if (h.indexOf("thalasseus acuflavidus") >= 0) {
            return "Cabots tern";
        } else if (h.indexOf("sterna hirundo") >= 0) {
            return "Common tern";
        } else if (h.equals("pefa")) {
            return "Peregrine falcon";
        } else if (h.equals("cygnus olor")) {
            return "Swan";
        } else {
            return host;
        }
    }

    public static String determineCountry(final String location) {
        final String l = location.toLowerCase();
        if (l.equals("alabama")
                || l.equals("alaska")
                || l.equals("arizona")
                || l.equals("arkansas")
                || l.equals("california")
                || l.equals("ca")
                || l.equals("colorado")
                || l.equals("connecticut")
                || l.equals("delaware")
                || l.equals("florida")
                || l.equals("georgia")
                || l.equals("hawaii")
                || l.equals("idaho")
                || l.equals("illinois")
                || l.equals("indiana")
                || l.equals("iowa")
                || l.equals("kansas")
                || l.equals("kentucky")
                || l.equals("louisiana")
                || l.equals("maine")
                || l.equals("ma")
                || l.equals("me")
                || l.equals("va")
                || l.equals("wy")
                || l.equals("mn")
                || l.equals("co")
                || l.equals("or")
                || l.equals("ia")
                || l.equals("id")
                || l.equals("nm")
                || l.equals("tx")
                || l.equals("mi")
                || l.equals("sd")
                || l.equals("ok")
                || l.equals("wa")
                || l.equals("pa")
                || l.equals("fl")
                || l.equals("mt")
                || l.equals("hi")
                || l.equals("ut")
                || l.equals("in")
                || l.equals("il")
                || l.equals("oh")
                || l.equals("md")
                || l.equals("sc")
                || l.equals("mo")
                || l.equals("az")
                || l.equals("la")
                || l.equals("tn")
                || l.equals("ks")
                || l.equals("ak")
                || l.equals("nd")
                || l.equals("ny")
                || l.equals("wi")
                || l.equals("ne")
                || l.equals("ky")
                || l.equals("wv")
                || l.equals("nc")
                || l.equals("ar")
                || l.equals("ky")
                || l.equals("ne")
                || l.equals("ms")
                || l.equals("wv")
                || l.equals("nj")
                || l.equals("vt")
                || l.equals("ct")
                || l.equals("nv")
                || l.equals("de")
                || l.equals("maryland")
                || l.equals("massachusetts")
                || l.equals("michigan")
                || l.equals("minnesota")
                || l.equals("mississippi")
                || l.equals("missouri")
                || l.equals("montana")
                || l.equals("nebraska")
                || l.equals("nevada")
                || l.equals("new hampshire")
                || l.equals("nh")
                || l.equals("new jersey")
                || l.equals("new mexico")
                || l.equals("new york")
                || l.equals("north carolina")
                || l.equals("north dakota")
                || l.equals("ohio")
                || l.equals("oklahoma")
                || l.equals("oregon")
                || l.equals("pennsylvania")
                || l.equals("rhode island")
                || l.equals("ri")
                || l.equals("south carolina")
                || l.equals("south dakota")
                || l.equals("tennessee")
                || l.equals("texas")
                || l.equals("utah")
                || l.equals("vermont")
                || l.equals("virginia")
                || l.equals("washington")
                || l.equals("west virginia")
                || l.equals("wisconsin")
                || l.equals("wyoming")
                || l.equals("washington")
        ) {
            return "USA";
        } else if (l.equals("hokkaido")
                || l.equals("kyoto")
                || l.equals("miyagi")
                || l.equals("tottori")
                || l.equals("akita")
                || l.equals("fukuoka")
                || l.equals("chiba")) {
            return "Japan";
        } else if (l.equals("ontario") || l.equals("victoria")) {
            return "Canada";
        } else if (l.equals("greenland")) {
            return "Denmark";
        } else if (l.equals("nanchang")
                || l.equals("guiyang")
                || l.equals("shantou")
                || l.equals("jilin")
                || l.equals("wuhan")
                || l.equals("xinjiang")
                || l.equals("hunan")
                || l.equals("jiangsu")
                || l.equals("hubei")
                || l.equals("shandong")
                || l.equals("yunnan")
                || l.equals("sichuan")
                || l.equals("fujian")
                || l.equals("guangxi")
                || l.equals("tonghai")
                || l.equals("qinghai")
                || l.equals("jiangxi")
                || l.equals("zhejiang")
                || l.equals("hebei")
                || l.equals("huadong")
                || l.equals("ningxia")
                || l.equals("anhui")
                || l.equals("karakol lake")
                || l.equals("sheny")
                || l.equals("eastern china")
                || l.equals("chn")
                || l.equals("beijing")) {
            return "China";
        } else if (l.equals("champasak")
                || l.equals("lao")
                || l.equals("xiangkhouang")
                || l.equals("luang namtha")) {
            return "Laos";
        } else if (l.equals("tyva") || l.equals("omsk")) {
            return "Russia";
        } else if (l.equals("arica y parinacota")
                || l.equals("maule")) {
            return "Chile";
        } else if (l.equals("korea")) {
            return "South Korea";
        } else if (l.equals("hongkong") || l.equals("hk")) {
            return "Hong Kong";
        } else if (l.equals("viet nam")
                || l.equals("vietnam hau giang")
        ) {
            return "Vietnam";
        } else if (l.equals("sidenreng rappang") || l.equals("sleman") || l.equals("klaten")
                || l.equals("majalengka")
                || l.equals("lamongan")
                || l.equals("westjava")
                || l.equals("east java")
                || l.equals("banten")
                || l.equals("pekalongan")
                || l.equals("banyuwangi")
                || l.equals("denpasar")) {
            return "Indonesia";
        } else if (l.equals("sagaing") || l.equals("yangon")) {
            return "Myanmar";
        } else if (l.equals("west bengal") || l.equals("sikkim")) {
            return "India";
        } else if (l.equals("gharbia") || l.equals("giza")) {
            return "Egypt";
        } else if (l.equals("england") || l.equals("scotland")) {
            return "United Kingdom";
        } else if (l.equals("mangystau")) {
            return "Kazakhstan";
        } else if (l.equals("united states") || l.equals("us")) {
            return "USA";
        } else if (l.equals("espiritosanto") || l.equals("itapemirimbr") || l.equals("saofranciscodeitabapoanabr")
                || l.equals("macaebr") || l.equals("piumabr") || l.equals("riodasostrasbr") || l.equals("cabo friobr")
                || l.equals("bertiogabr")) {
            return "Brazil";
        } else if (l.equals("germany-nw")) {
            return "Germany";
        } else if (l.equals("veracruz")) {
            return "Mexico";
        } else if (l.equals("british columbia")) {
            return "Canada";
        } else if (l.startsWith("germany")) {
            return "Germany";
        } else if (l.equals("seedorf")) {
            return "Germany";
        } else {
            return location;
        }
    }

    public static String determineState(final String location) {
        final String l = location.toLowerCase();
        if (l.equals("nh")) {
            return "New Hampshire";
        }
        if (l.equals("ma")) {
            return "Massachusetts";
        }
        if (l.equals("me")) {
            return "Maine";
        }
        if (l.equals("vt")) {
            return "Vermont";
        }
        if (l.equals("wy")) {
            return "Wyoming";
        }
        if (l.equals("mn")) {
            return "Minnesota";
        }
        if (l.equals("ia")) {
            return "Iowa";
        }
        if (l.equals("co")) {
            return "Colorado";
        }
        if (l.equals("id")) {
            return "Idaho";
        }
        if (l.equals("mi")) {
            return "Michigan";
        }
        if (l.equals("tx")) {
            return "Texas";
        }
        if (l.equals("nm")) {
            return "New Mexico";
        }
        if (l.equals("sd")) {
            return "South Dakota";
        }
        if (l.equals("ny")) {
            return "New York";
        }

        if (l.equals("alabama")
                || l.equals("alaska")
                || l.equals("arizona")
                || l.equals("arkansas")
                || l.equals("california")
                || l.equals("colorado")
                || l.equals("connecticut")
                || l.equals("delaware")
                || l.equals("florida")
                || l.equals("georgia")
                || l.equals("hawaii")
                || l.equals("idaho")
                || l.equals("illinois")
                || l.equals("indiana")
                || l.equals("iowa")
                || l.equals("kansas")
                || l.equals("kentucky")
                || l.equals("louisiana")
                || l.equals("maine")
                || l.equals("maryland")
                || l.equals("massachusetts")
                || l.equals("michigan")
                || l.equals("minnesota")
                || l.equals("mississippi")
                || l.equals("missouri")
                || l.equals("montana")
                || l.equals("nebraska")
                || l.equals("nevada")
                || l.equals("new hampshire")
                || l.equals("new jersey")
                || l.equals("new mexico")
                || l.equals("new york")
                || l.equals("north carolina")
                || l.equals("north dakota")
                || l.equals("ohio")
                || l.equals("oklahoma")
                || l.equals("oregon")
                || l.equals("pennsylvania")
                || l.equals("rhode island")
                || l.equals("south carolina")
                || l.equals("south dakota")
                || l.equals("tennessee")
                || l.equals("texas")
                || l.equals("utah")
                || l.equals("vermont")
                || l.equals("virginia")
                || l.equals("washington")
                || l.equals("west virginia")
                || l.equals("wisconsin")
                || l.equals("wyoming")
                || l.equals("ontario")
                || l.equals("victoria")
                || l.equals("greenland")
        ) {
            return toTitleCase(location);
        } else {
            return "";
        }
    }

    public final static String toTitleCase(String words) {
        return Stream.of(words.trim().split("\\s"))
                .filter(word -> word.length() > 0)
                .map(word -> word.substring(0, 1).toUpperCase() + word.substring(1))
                .collect(Collectors.joining(" "));
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
        } else if (c.equals("australia") || c.equals("new_zealand")) {
            region = "Oceania";
        } else if (c.equals("dominican republic") || c.equals("puerto rico") || c.equals("jamaica")
                || c.equals("belize")) {
            region = "Caribbean";
        } else if (c.equals("antarctica")) {
            region = "Antarctica";
        } else if (c.equals("morocco") || c.equals("tunisia") || c.equals("libya")) {
            region = "North Africa";
        } else if (c.equals("mauritania") || c.equals("gambia") || c.equals("senegal") || c.equals("ghana") || c.equals("benin")
                || c.equals("liberia") || c.equals("nigeria") || c.equals("niger") || c.equals("togo") || c.equals("sierra leone")
                || c.equals("guinea") || c.equals("cote d'ivoire") || c.equals("mali")) {
            region = "West Africa";
        } else if (c.equals("kenya") || c.equals("tanzania") || c.equals("uganda") || c.equals("rwanda") || c.equals("ethiopia") || c.equals("djibouti")) {
            region = "East Africa";
        } else if (c.equals("angola") || c.equals("chad") || c.equals("gabon") || c.equals("central african republic") || c.equals("cameroon") || c.indexOf("congo") >= 0) {
            region = "Central Africa";
        } else if (c.equals("south africa") || c.equals("malawi") || c.equals("botswana") || c.equals("lesotho") || c.equals("namibia")) {
            region = "Southern Africa.";
        } else if (c.equals("na")) {
            region = "";
        } else {
            System.out.println("Warning: unknown country \"" + c + "\"");
            region = "";
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

        if (year_int < 1800 || year_int > 2025) {
            System.out.println("Error year \"" + year + "\" is out of range");
            System.exit(-1);
        }
        return Integer.toString(year_int);
    }

    public static String cleanHostString(final String s) {
        final String r = s.replaceAll("_", " ");
        return r.substring(0, 1).toUpperCase() + r.substring(1).toLowerCase();
    }

    public static String cleanLocationString(final String s) {
        final String r = s.replaceAll("_", " ");
        return r;
    }


    public static void addHostGroup(final String host, final PropertiesList custom_data,
                                    final String host_group_ref,
                                    final String host_group_dom_vs_wild_ref) {
        String hg1 = UNKNOWN;
        String hg2 = UNKNOWN;
        final String h = host.toLowerCase();
        if (h.length() == 0 || h.equalsIgnoreCase(UNKNOWN) || h.equalsIgnoreCase("na") || h.equals("?")) {
            hg1 = UNKNOWN;
            hg2 = UNKNOWN;
        } else if (h.equals("human")
                || h.equals("homo sapiens")) {
            hg1 = "Human";
            hg2 = "Human";
        } else if (h.equals("cattle")
                || h.equals("cow")
                || h.equals("bovine")
                || h.equals("equine")
                || h.equals("goat")
                || h.equals("sheep")
                || h.equals("swine")
                || h.equals("cat")
                || h.equals("feline")
                || h.equals("dog")
                || h.equals("canine")
                || h.equals("mink")
                || h.equals("pig")
                || h.equals("ferret")
                || h.equals("bovine milk")
                || h.equals("stone marten")) {
            hg1 = "Non-Human Mammal";
            hg2 = "Non-Human Mammal (domestic)";
        } else if (h.equals("skunk")
                || h.equals("mountain lion")
                || h.equals("raccoon")
                || h.equals("pika")
                || h.equals("plateau pika")
                || h.equals("harbor seal")
                || h.equals("south american sea lion")
                || h.equals("tiger")
                || h.equals("owston's civet")
                || h.equals("civet")
                || h.equals("alpaca")
                || h.equals("red fox")
                || h.equals("panthera pardus")
                || h.equals("panthera tigris")
                || h.equals("house mouse")
                || h.equals("bobcat")
                || h.equals("bottlenose dolphin")
                || h.startsWith("mustela")
                || h.equals("serval")
                || h.equals("fox")
                || h.equals("lion")
                || h.equals("ermine")
        ) {
            hg1 = "Non-Human Mammal";
            hg2 = "Non-Human Mammal (wild)";

        } else if (h.equals("chicken")
                || h.equals("domestic duck")
                || h.equals("mallard")
                || h.equals("turkey")
                || h.equals("domestic goose")
                || h.equals("muscovy duck")
                || h.equals("poultry")) {
            hg1 = "Avian";
            hg2 = "Avian (domestic)";
        } else if (h.equals("openbill stork")
                || h.equals("pigeon")
                || h.equals("wild duck")
                || h.equals("peregrine falcon")
                || h.equals("common buzzard")
                || h.equals("bald eagle")
                || h.equals("bar-headed goose")
                || h.equals("whooper swan")
                || h.equals("quail")
                || h.equals("great crested grebe")
                || h.equals("crow")
                || h.equals("red-tailed hawk")
                || h.equals("blue-winged teal")
                || h.equals("american green-winged teal")
                || h.equals("cormorant")
                || h.equals("jungle crow")
                || h.equals("great black-backed gull")
                || h.equals("eurasian eagle owl")
                || h.equals("great horned owl")
                || h.equals("chukar")
                || h.equals("ruddy turnstone")
                || h.equals("teal")
                || h.equals("pacific black duck")
                || h.equals("gadwall duck")
                || h.equals("grackle")
                || h.equals("mute swan")
                || h.equals("canada goose")
                || h.equals("western gull")
                || h.equals("snow goose")
                || h.equals("common raven")
                || h.equals("harris-hawk")
                || h.equals("turkey vulture")
                || h.equals("american crow")
                || h.equals("redhead duck")
                || h.equals("american white pelican")
                || h.equals("wood duck")
                || h.equals("hooded merganser")
                || h.equals("blackbird")
                || h.equals("cabots tern")
                || h.equals("aves")
                || h.equals("pelican")
                || h.equals("black vulture")
                || h.equals("sanderling")
                || h.equals("red tailed hawk")
                || h.equals("harris hawk")
                || h.equals("hawk")
                || h.equals("american wigeon")
                || h.equals("western sandpiper")
                || h.equals("common grackle")
                || h.equals("common tern")
                || h.equals("common loon")
                || h.equals("white-winged scoter")
                || h.equals("common eider")
                || h.equals("surf scoter")
                || h.equals("herring gull")
                || h.equals("scoter")
                || h.equals("brandt goose")
                || h.equals("black scoter")
                || h.equals("lesser scaup")
                || h.equals("house sparrow")
                || h.equals("wild-bird")
                || h.equals("wild goose")
                || h.equals("american robin")
                || h.equals("western kingbird")
                || h.equals("eurasian collared dove")
                || h.equals("peregrine falcon")
                || h.equals("south polar skua")
                || h.equals("horned grebe")
                || h.equals("emu")
                || h.equals("guineafowl")
                || h.equals("peafowl")
                || h.equals("seagull")
                || h.equals("dove")
                || h.equals("swan")
                || h.equals("eurasian goshawk")
                || h.equals("ostrich")
                || h.equals("pheasant")
                || h.equals("vulture")
                || h.equals("rock pigeon")
                || h.equals("northern shoveler")
                || h.equals("northern pintail")
                || h.equals("white-winged dove")
                || h.indexOf("stork") > -1
                || h.indexOf("duck") > -1
                || h.indexOf("pelican") > -1
                || h.indexOf("goose") > -1
                || h.indexOf("owl") > -1
                || h.indexOf("hawk") > -1
                || h.indexOf("teal") > -1
                || h.indexOf("heron") > -1
                || h.indexOf("falcon") > -1
                || h.indexOf("ibis") > -1
                || h.indexOf("gull") > -1
                || h.indexOf("swan") > -1
                || h.indexOf("egret") > -1
                || h.indexOf("marganser") > -1
                || h.indexOf("flamingo") > -1
        ) {
            hg1 = "Avian";
            hg2 = "Avian (wild)";
        } else if (h.equals("environment")) {
            hg1 = "Environment";
            hg2 = "Environment";
        } else {
            System.out.println("Warning: Unknown host \"" + host + "\"");
            hg1 = "";
            hg2 = "";
        }

        if (!ForesterUtil.isEmpty(hg1)) {
            custom_data.addProperty(new Property(host_group_ref, hg1, "", XSD_STRING, Property.AppliesTo.NODE));
        }
        if (!ForesterUtil.isEmpty(hg2)) {
            custom_data.addProperty(new Property(host_group_dom_vs_wild_ref, hg2, "", XSD_STRING, Property.AppliesTo.NODE));
        }
    }
}
