package edu.carnegiescience.dpb.evanslab;

import java.util.List;
import java.util.LinkedList;
import java.util.Set;
import java.util.HashSet;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Filters SNPs that appear in the first and second VCF files but NOT in the third VCF file.
 *
 * Parameters:
 *
 * -v1 --VCF1 VCF file 1, contains SNPs to consider
 * -v2 --VCF2 VCF file 2, contains SNPs that must match VCF1 location and ALT values
 * -v3 --VCF3 VCF file 3, contains SNPs that must NOT be present at VCF1 location or match VCF1 ALT values if present
 * -g  --GFF  GFF file containing regions (e.g. genes) that should not contain ANY matching SNPs from VCF3 (optional)
 *
 * @author Sam Hokin
 */
public class SNPFilter {

    /**
     * Main class does all the work.
     */
    public static void main(String[] args) throws Exception {

        Options options = new Options();

        Option vcf1Option = new Option("v1", "VCF1", true, "VCF file 1, contains SNPs to consider");
        vcf1Option.setRequired(true);
        options.addOption(vcf1Option);

        Option vcf2Option = new Option("v2", "VCF2", true, "contains SNPs that must match VCF1 location and ALT values");
        vcf2Option.setRequired(true);
        options.addOption(vcf2Option);

        Option vcf3Option = new Option("v3", "VCF3", true, "VCF file 3, contains SNPs that must NOT be present at VCF1 location or match VCF1 ALT values if present");
        vcf3Option.setRequired(true);
        options.addOption(vcf3Option);

        Option gffOption = new Option("g", "GFF", true, "GFF file containing regions (e.g. genes) that should not contain ANY matching SNPs from VCF3 (optional)");
        gffOption.setRequired(false);
        options.addOption(gffOption);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("SNPFilter", options);
            System.exit(1);
            return;
        }

        // filenames
        String vcf1Filename = cmd.getOptionValue("VCF1");
        String vcf2Filename = cmd.getOptionValue("VCF2");
        String vcf3Filename = cmd.getOptionValue("VCF3");
        String gffFilename = null; if (cmd.getOptionValue("GFF")!=null) gffFilename = cmd.getOptionValue("GFF");

        // output the parameters
        System.out.println("edu.carnegiescience.dpb.evanslab.SNPFilter");
        System.out.println("vcf1Filename:"+vcf1Filename);
        System.out.println("vcf2Filename:"+vcf2Filename);
        System.out.println("vcf3Filename:"+vcf3Filename);
        if (gffFilename!=null) System.out.println("gffFilename:"+gffFilename);
        System.out.println();

        // instantiate the VCF loaders
        VCFLoader vcf1Loader = new VCFLoader(vcf1Filename);
        VCFLoader vcf2Loader = new VCFLoader(vcf2Filename);
        VCFLoader vcf3Loader = new VCFLoader(vcf3Filename);

        // and the optional GFF loader
        GFFLoader gffLoader = null;
        if (gffFilename!=null) gffLoader = new GFFLoader(gffFilename);

        // list of "winning" locations for output
        List<String> winningLocations = new LinkedList<>();

        // list of "losing" locations for filtering out GFF regions
        List<String> losingLocations = new LinkedList<>();
        
        // spin through the source VCF file
        vcf1Loader.load();
        for (VariantContext vc1 : vcf1Loader.vcList) {
            // limit to SNPs
            if (vc1.isSNP()) {
                // get vc1 data
                String vc1ID = vc1.getID();
                String vc1Contig = vc1.getContig();
                int vc1Start = vc1.getStart();
                int vc1End = vc1.getEnd();
                Allele vc1Ref = vc1.getReference();
                String vc1RefString = vc1Ref.getBaseString();
                String vc1AltString = getAltString(vc1);
                // search for this location on VCF2
                List<VariantContext> vc2List = vcf2Loader.query(vc1Contig, vc1Start, vc1End).toList();
                for (VariantContext vc2 : vc2List) {
                    // get vc2 data
                    String vc2ID = vc2.getID();
                    String vc2Contig = vc2.getContig();
                    int vc2Start = vc2.getStart();
                    int vc2End = vc2.getEnd();
                    Allele vc2Ref = vc2.getReference();
                    String vc2RefString = vc2Ref.getBaseString();
                    String vc2AltString = getAltString(vc2);
                    // check for ALT match between VCF1 and VCF2
                    if (vc1AltString.equals(vc2AltString)) {
                        boolean winner = true;
                        // search for this location on VCF3
                        List<VariantContext> vc3List = vcf3Loader.query(vc1Contig, vc1Start, vc1End).toList();
                        for (VariantContext vc3 : vc3List) {
                            // get vc3 data
                            String vc3ID = vc3.getID();
                            String vc3Contig = vc3.getContig();
                            int vc3Start = vc3.getStart();
                            int vc3End = vc3.getEnd();
                            Allele vc3Ref = vc3.getReference();
                            String vc3RefString = vc3Ref.getBaseString();
                            String vc3AltString = getAltString(vc3);
                            // we allow VCF3 to have a DIFFERENT ALT allele at this position
                            winner = (!vc1AltString.equals(vc3AltString));
                        }
                        if (winner) {
                            winningLocations.add(vc1Contig+":"+vc1Start);
                        } else {
                            losingLocations.add(vc1Contig+":"+vc1Start);
                        }
                    }
                }
            }
        }

        // Optionally collect the genes that contain losing locations for further rejection
        Set<String> losingGenes = new HashSet<>();
        if (gffLoader!=null) {
            for (String loc : losingLocations) {
                String[] parts = loc.split(":");
                String contig = parts[0];
                Integer pos = Integer.parseInt(parts[1]);
                Location location = new Location(pos, pos);
                FeatureList overlapping = gffLoader.search(contig, location);
                for (FeatureI feature : overlapping) {
                    String geneID = feature.getAttribute("ID");
                    losingGenes.add(geneID);
                }
            }
        }

        // output header
        System.out.println("Contig\tPos\tREF\tALT\t"+
                           "Var1RF\tVar1RR\tVar1AF\tVar1AR\t"+
                           "Var2RF\tVar2RR\tVar2AF\tVar2AR");
        // run through the winners, outputting the relevant data
        for (String loc : winningLocations) {
            String[] parts = loc.split(":");
            String contig = parts[0];
            Integer pos = Integer.parseInt(parts[1]);
            // optional removal if this is on a losing gene
            String losingGene = null;
            if (gffLoader!=null) {
                Location location = new Location(pos, pos);
                FeatureList overlapping = gffLoader.search(contig, location);
                for (FeatureI feature : overlapping) {
                    String geneID = feature.getAttribute("ID");
                    if (losingGenes.contains(geneID)) losingGene = geneID;
                }
            }
            if (losingGene==null) {
                List<VariantContext> vc1List = vcf1Loader.query(contig, pos, pos).toList();
                for (VariantContext vc1 : vc1List) {
                    // get vc1 data
                    String vc1ID = vc1.getID();
                    String vc1Contig = vc1.getContig();
                    int vc1Start = vc1.getStart();
                    int vc1End = vc1.getEnd();
                    Allele vc1Ref = vc1.getReference();
                    String vc1RefString = vc1Ref.getBaseString();
                    String vc1AltString = getAltString(vc1);
                    List<Integer> vc1DP4 = vc1.getAttributeAsIntList("DP4", 0);
                    int vc1RefForward = vc1DP4.get(0);
                    int vc1RefReverse = vc1DP4.get(1);
                    int vc1AltForward = vc1DP4.get(2);
                    int vc1AltReverse = vc1DP4.get(3); 
                    List<VariantContext> vc2List = vcf2Loader.query(vc1Contig, vc1Start, vc1End).toList();
                    for (VariantContext vc2 : vc2List) {
                        // get vc2 data
                        String vc2ID = vc2.getID();
                        String vc2Contig = vc2.getContig();
                        int vc2Start = vc2.getStart();
                        int vc2End = vc2.getEnd();
                        Allele vc2Ref = vc2.getReference();
                        String vc2RefString = vc2Ref.getBaseString();
                        String vc2AltString = getAltString(vc2);
                        List<Integer> vc2DP4 = vc2.getAttributeAsIntList("DP4", 0);
                        int vc2RefForward = vc2DP4.get(0);
                        int vc2RefReverse = vc2DP4.get(1);
                        int vc2AltForward = vc2DP4.get(2);
                        int vc2AltReverse = vc2DP4.get(3); 
                        // output
                        System.out.println(vc1Contig+"\t"+vc1Start+"\t"+vc1RefString+"\t"+vc1AltString+"\t"+
                                           vc1RefForward+"\t"+vc1RefReverse+"\t"+vc1AltForward+"\t"+vc1AltReverse+"\t"+
                                           vc2RefForward+"\t"+vc2RefReverse+"\t"+vc2AltForward+"\t"+vc2AltReverse);
                    }
                }
            }
        }
    }

    /**
     * Return a comma-delimited ALT string from a VariantContext
     */
    static String getAltString(VariantContext vc) {
        List<Allele> vcAlts = vc.getAlternateAlleles();
        String vcAltString = "";
        for (Allele alt : vcAlts) {
            if (vcAltString.length()>0) vcAltString += ",";
            vcAltString += alt.getBaseString();
        }
        return vcAltString;
    }

}
