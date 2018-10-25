package edu.carnegiescience.dpb.evanslab;

import java.util.List;
import java.util.LinkedList;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import htsjdk.samtools.util.CloseableIterator;

public class SNPComparer {

    // parameter defaults
    static int SOURCE_ALT_TOTAL_MIN = 100;
    static double SOURCE_ALT_READ_RATIO_MIN = 0.5;
    static double SOURCE_REF_FRACTION_MAX = 0.5;
    static double TARGET_REF_FRACTION_MIN =  0.5;

    /**
     * Main class does all the work.
     */
    public static void main(String[] args) throws Exception {
        if (args.length<4) {
            System.out.println("Usage: SNPComparer <source-vcf-file> <remap-gff-file> <target-gff-file> <target-vcf-file> "+
                               "[sourceAltTotalMin("+SOURCE_ALT_TOTAL_MIN+")] "+
                               "[sourceAltReadRatioMin("+SOURCE_ALT_READ_RATIO_MIN+")] "+
                               "[sourceRefFractionMax("+SOURCE_REF_FRACTION_MAX+")] " +
                               "[targetRefFractionMin("+TARGET_REF_FRACTION_MIN+")]");
            System.exit(0);
        }

        // defaults
        int sourceAltTotalMin = SOURCE_ALT_TOTAL_MIN;
        double sourceAltReadRatioMin = SOURCE_ALT_READ_RATIO_MIN;
        double sourceRefFractionMax = SOURCE_REF_FRACTION_MAX;
        double targetRefFractionMin = TARGET_REF_FRACTION_MIN;

        String sourceVCFFilename = args[0];
        String remapGFFFilename = args[1];
        String targetGFFFilename = args[2];
        String targetVCFFilename = args[3];

        if (args.length>=5) sourceAltTotalMin = Integer.parseInt(args[4]);
        if (args.length>=6) sourceAltReadRatioMin = Double.parseDouble(args[5]);
        if (args.length>=7) sourceRefFractionMax = Double.parseDouble(args[6]);
        if (args.length>=8) targetRefFractionMin = Double.parseDouble(args[7]);

        VCFLoader sourceVCFLoader = new VCFLoader(sourceVCFFilename);
        GFFLoader remapGFFLoader = new GFFLoader(remapGFFFilename);
        GFFLoader targetGFFLoader = new GFFLoader(targetGFFFilename);
        VCFLoader targetVCFLoader = new VCFLoader(targetVCFFilename);

        // output parameters
        System.out.println("edu.carnegiescience.dpb.evanslab.SNPComparer");
        System.out.println("sourceVCFFilename:"+sourceVCFFilename);
        System.out.println("remapGFFFilename:"+remapGFFFilename);
        System.out.println("targetGFFFilename:"+targetGFFFilename);
        System.out.println("targetVCFFilename:"+targetVCFFilename);
        System.out.println("sourceAltTotalMin="+sourceAltTotalMin);
        System.out.println("sourceAltReadRatioMin="+sourceAltReadRatioMin);
        System.out.println("sourceRefFractionMax="+sourceRefFractionMax);
        System.out.println("targetRefFractionMin="+targetRefFractionMin);
        System.out.println();
        
        // output header
        System.out.println("SrcContig\tSrcPos\tSrcRef\tSrcAlt\tSrcRF\tSrcRR\tSrcAF\tSrcAR\tSrcRefFrac\tGeneID\tMinTargetRefFrac");
        
        // load the source VCF file
        sourceVCFLoader.load();

        for (VariantContext sourceVC : sourceVCFLoader.vcList) {
            if (sourceVC.isSNP()) {

                // source VCF values
                String sourceID = sourceVC.getID();
                String sourceContig = sourceVC.getContig();
                int sourceStart = sourceVC.getStart();
                Allele sourceRef = sourceVC.getReference();
                List<Allele> sourceAlts = sourceVC.getAlternateAlleles();
                List<Integer> dp4List = sourceVC.getAttributeAsIntList("DP4", 0);
                int sourceRefForward = dp4List.get(0);
                int sourceRefReverse = dp4List.get(1);
                int sourceAltForward = dp4List.get(2);
                int sourceAltReverse = dp4List.get(3); 
                int sourceRefTotal = sourceRefForward + sourceRefReverse;
                int sourceAltTotal = sourceAltForward + sourceAltReverse;
                String sourceRefString = sourceRef.getBaseString();
                String sourceAltString = "";
                for (Allele alt : sourceAlts) {
                    if (sourceAltString.length()>0) sourceAltString += ",";
                    sourceAltString += alt.getBaseString();
                }
                double sourceRefFraction = (double)(sourceRefTotal)/(double)(sourceRefTotal+sourceAltTotal);

                // filtering!
                boolean sourceOK = (sourceAltTotal>sourceAltTotalMin)
                    && ((double)Math.min(sourceAltForward,sourceAltReverse)/(double)Math.max(sourceAltForward,sourceAltReverse)>sourceAltReadRatioMin)
                    && (sourceRefFraction<sourceRefFractionMax);
                
                if (sourceOK) {

                    // search for the target gene(s) spanning this location on the source genome
                    Location location = new Location(sourceStart,sourceStart);
                    FeatureList overlapping = remapGFFLoader.search(sourceContig, location);
                    for (FeatureI feature : overlapping) {
                        String geneID = feature.getAttribute("ID");

                        // find this gene on the target genome
                        FeatureList genes = targetGFFLoader.searchID(geneID);
                        for (FeatureI gene : genes) {
                            String chromosome = gene.seqname();
                            String type = gene.type();
                            Location loc = gene.location();
                            int start = loc.start();
                            int end = loc.end();
                            char strand = '+';
                            if (start<0) {
                                strand = '-';
                                int minusStart = -start;
                                int minusEnd = -end;
                                start = minusEnd;
                                end = minusStart;
                            }
                            // now search the target VCF for SNPs on the target genome
                            double minRefFraction = 1.0;
                            CloseableIterator<VariantContext> iterator = targetVCFLoader.query(chromosome, start, end);
                            while (iterator.hasNext()) {
                                VariantContext targetVC = iterator.next();
                                // target VCF values
                                int targetStart = targetVC.getStart();
                                Allele targetRef = targetVC.getReference();
                                List<Allele> targetAlts = targetVC.getAlternateAlleles();
                                List<Integer> targetDP4List = targetVC.getAttributeAsIntList("DP4", 0);
                                int targetRefForward = targetDP4List.get(0);
                                int targetRefReverse = targetDP4List.get(1);
                                int targetAltForward = targetDP4List.get(2);
                                int targetAltReverse = targetDP4List.get(3);
                                int targetRefTotal = targetRefForward + targetRefReverse;
                                int targetAltTotal = targetAltForward + targetAltReverse;
                                String targetRefString = targetRef.getBaseString();
                                String targetAltString = "";
                                for (Allele alt : targetAlts) {
                                    if (targetAltString.length()>0) targetAltString += ",";
                                    targetAltString += alt.getBaseString();
                                }
                                double targetRefFraction = (double)targetRefTotal/(double)(targetRefTotal+targetAltTotal);
                                if (targetRefFraction<minRefFraction) minRefFraction = targetRefFraction;
                            }
                            // output record if passes filter
                            if (minRefFraction>targetRefFractionMin) {
                                System.out.println(sourceContig+"\t"+sourceStart+
                                                   "\t"+sourceRefString+"\t"+sourceAltString+
                                                   "\t"+sourceRefForward+"\t"+sourceRefReverse+
                                                   "\t"+sourceAltForward+"\t"+sourceAltReverse+
                                                   "\t"+sourceRefFraction+
                                                   "\t"+geneID+"\t"+minRefFraction);
                            }
                        }
                    }
                }
            }
        }
    }

}
