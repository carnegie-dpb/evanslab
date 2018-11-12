package edu.carnegiescience.dpb.evanslab;

import java.util.List;
import java.util.LinkedList;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

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
 * -v2 --VCF2 VCF file 2, contains SNPs that must match VCF1 location for output
 * -v3 --VCF3 VCF file 3, contains SNPs that must NOT match VCF1 for output
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

        Option vcf2Option = new Option("v2", "VCF2", true, "VCF file 2, contains SNPs that must match VCF1 location for output");
        vcf2Option.setRequired(true);
        options.addOption(vcf2Option);

        Option vcf3Option = new Option("v3", "VCF3", true, "VCF file 3, contains SNPs that must NOT match VCF1 for output");
        vcf3Option.setRequired(true);
        options.addOption(vcf3Option);


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

        // output the parameters
        System.out.println("edu.carnegiescience.dpb.evanslab.SNPFilter");
        System.out.println("vcf1Filename:"+vcf1Filename);
        System.out.println("vcf2Filename:"+vcf2Filename);
        System.out.println("vcf3Filename:"+vcf3Filename);
        System.out.println();

        // output header
        System.out.println("Contig\tPos\tREF\tALT\t"+
                           "Var1RF\tVar1RR\tVar1AF\tVar1AR\t"+
                           "Var2RF\tVar2RR\tVar2AF\tVar2AR");
        
        // instantiate the VCF loaders
        VCFLoader vcf1Loader = new VCFLoader(vcf1Filename);
        VCFLoader vcf2Loader = new VCFLoader(vcf2Filename);
        VCFLoader vcf3Loader = new VCFLoader(vcf3Filename);
        
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
                List<Allele> vc1Alts = vc1.getAlternateAlleles();
                String vc1AltString = "";
                for (Allele alt : vc1Alts) {
                    if (vc1AltString.length()>0) vc1AltString += ",";
                    vc1AltString += alt.getBaseString();
                }
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
                    List<Allele> vc2Alts = vc2.getAlternateAlleles();
                    String vc2AltString = "";
                    for (Allele alt : vc2Alts) {
                        if (vc2AltString.length()>0) vc2AltString += ",";
                        vc2AltString += alt.getBaseString();
                    }
                    // check for ALT match between VCF1 and VCF2
                    if (vc1AltString.equals(vc2AltString)) {
                        // search for this location on VCF3
                        List<VariantContext> vc3List = vcf3Loader.query(vc1Contig, vc1Start, vc1End).toList();
                        if (vc3List.size()==0) {
                            // We have a winner!
                            List<Integer> vc1DP4 = vc1.getAttributeAsIntList("DP4", 0);
                            int vc1RefForward = vc1DP4.get(0);
                            int vc1RefReverse = vc1DP4.get(1);
                            int vc1AltForward = vc1DP4.get(2);
                            int vc1AltReverse = vc1DP4.get(3); 
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

    //             Allele sourceuRef = sourceVC.getReference();
    //             List<Allele> sourceAlts = sourceVC.getAlternateAlleles();
    //             List<Integer> dp4List = sourceVC.getAttributeAsIntList("DP4", 0);
    //             int sourceRefForward = dp4List.get(0);
    //             int sourceRefReverse = dp4List.get(1);
    //             int sourceAltForward = dp4List.get(2);
    //             int sourceAltReverse = dp4List.get(3); 
    //             int sourceRefTotal = sourceRefForward + sourceRefReverse;
    //             int sourceAltTotal = sourceAltForward + sourceAltReverse;
    //             String sourceRefString = sourceRef.getBaseString();
    //             boolean sourceIsHet = sourceAlts.size()>1;
    //             String sourceAltString = "";
    //             for (Allele alt : sourceAlts) {
    //                 if (sourceAltString.length()>0) sourceAltString += ",";
    //                 sourceAltString += alt.getBaseString();
    //             }
    //             double sourceAltFraction = (double)(sourceAltTotal)/(double)(sourceRefTotal+sourceAltTotal);

    //             // source filtering NOTE: only homozygous calls allowed!
    //             boolean sourceOK =
    //                 (!sourceIsHet)
    //                 && (sourceAltTotal>sourceAltTotalMin)
    //                 && ((double)Math.min(sourceAltForward,sourceAltReverse)/(double)Math.max(sourceAltForward,sourceAltReverse)>sourceAltReadRatioMin)
    //                 && (sourceAltFraction>=sourceAltFractionMin);
                
    //             if (sourceOK) {

    //                 // search for the target gene(s) spanning this location on the source genome
    //                 Location sourceLocation = new Location(sourceStart,sourceStart);
    //                 FeatureList overlapping = remapGFFLoader.search(sourceContig, sourceLocation);
    //                 for (FeatureI feature : overlapping) {
    //                     String geneID = feature.getAttribute("ID");
    //                     // only process new genes
    //                     if (!targetGeneList.contains(geneID)) {

    //                         // find this gene on the target genome
    //                         FeatureList genes = targetGFFLoader.searchID(geneID);
    //                         for (FeatureI gene : genes) {
    //                             String chromosome = gene.seqname();
    //                             String type = gene.type();
    //                             Location loc = gene.location();
    //                             int start = loc.start();
    //                             int end = loc.end();
    //                             char strand = '+';
    //                             // if minus strand, indicate with "-" but make start<end
    //                             if (start<0) {
    //                                 strand = '-';
    //                                 int minusStart = -start;
    //                                 int minusEnd = -end;
    //                                 start = minusEnd;
    //                                 end = minusStart;
    //                             }
                                
    //                             // now search the target VCF for SNPs on the target genome
    //                             List<VariantContext> targetVCList = targetVCFLoader.query(chromosome, start, end).toList();
    //                             boolean targetHasSNPs = targetVCList.size()>0;
    //                             double minTargetRefFraction = 1.0;
    //                             for (VariantContext targetVC : targetVCList) {
    //                                 int targetStart = targetVC.getStart();
    //                                 Allele targetRef = targetVC.getReference();
    //                                 List<Allele> targetAlts = targetVC.getAlternateAlleles();
    //                                 List<Integer> targetDP4List = targetVC.getAttributeAsIntList("DP4", 0);
    //                                 int targetRefForward = targetDP4List.get(0);
    //                                 int targetRefReverse = targetDP4List.get(1);
    //                                 int targetAltForward = targetDP4List.get(2);
    //                                 int targetAltReverse = targetDP4List.get(3);
    //                                 int targetRefTotal = targetRefForward + targetRefReverse;
    //                                 int targetAltTotal = targetAltForward + targetAltReverse;
    //                                 String targetRefString = targetRef.getBaseString();
    //                                 String targetAltString = "";
    //                                 for (Allele alt : targetAlts) {
    //                                     if (targetAltString.length()>0) targetAltString += ",";
    //                                     targetAltString += alt.getBaseString();
    //                                 }
    //                                 double targetRefFraction = (double)targetRefTotal/(double)(targetRefTotal+targetAltTotal);
    //                                 minTargetRefFraction = Math.min(targetRefFraction,minTargetRefFraction);
    //                             }
                                    
    //                             // output record if passes target filter
    //                             boolean targetOK = (!targetHasSNPs) || (minTargetRefFraction>=targetRefFractionMin);
    //                             if (targetOK) {
    //                                 System.out.println(geneID+"\t"+chromosome+"\t"+start+"\t"+end+"\t"+strand+"\t"+minTargetRefFraction);
    //                                 targetGeneList.add(geneID);
    //                             }
                                    
    //                         }
                            
    //                     }
    //                 }
    //             }
    //         }
    //     }

    }

}
