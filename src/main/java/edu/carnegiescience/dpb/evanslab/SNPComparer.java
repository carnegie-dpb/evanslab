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

    /**
     * Main class does all the work.
     */
    public static void main(String[] args) throws Exception {
        if (args.length!=4) {
            System.out.println("Usage SNPComparer <source-vcf-file> <remap-gff-file> <target-gff-file> <target-vcf-file>");
            System.exit(0);
        }

        String sourceVCFFilename = args[0];
        String remapGFFFilename = args[1];
        String targetGFFFilename = args[2];
        String targetVCFFilename = args[3];
        
        VCFLoader sourceVCFLoader = new VCFLoader(sourceVCFFilename);
        GFFLoader remapGFFLoader = new GFFLoader(remapGFFFilename);
        GFFLoader targetGFFLoader = new GFFLoader(targetGFFFilename);
        VCFLoader targetVCFLoader = new VCFLoader(targetVCFFilename);
        
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
                int sourceRefTotal = sourceRefForward + sourceRefReverse;
                int sourceAltForward = dp4List.get(2);
                int sourceAltReverse = dp4List.get(3);
                int sourceAltTotal = sourceAltForward + sourceAltReverse;
                String sourceRefString = sourceRef.getBaseString();
                String sourceAltString = "";
                for (Allele alt : sourceAlts) {
                    if (sourceAltString.length()>0) sourceAltString += ",";
                    sourceAltString += alt.getBaseString();
                }

                // filtering!
                boolean sourceOK = (sourceAltTotal>100)
                    && ((double)sourceAltForward/(double)sourceAltTotal>0.2)
                    && (sourceRefTotal/sourceAltTotal<0.5);

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
                                start = -start;
                                end = -end;
                            }


                            // now search the target VCF for SNPs on the target genome
                            double maxAltRatio = 0.0;
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
                                int targetRefTotal = targetRefForward + targetRefReverse;
                                int targetAltForward = targetDP4List.get(2);
                                int targetAltReverse = targetDP4List.get(3);
                                int targetAltTotal = targetAltForward + targetAltReverse;
                                String targetRefString = targetRef.getBaseString();
                                String targetAltString = "";
                                for (Allele alt : targetAlts) {
                                    if (targetAltString.length()>0) targetAltString += ",";
                                    targetAltString += alt.getBaseString();
                                }
                                double altRatio = (double)targetAltTotal/(double)sourceAltTotal;
                                if (altRatio>maxAltRatio) maxAltRatio = altRatio;
                            }
                            // output
                            if (maxAltRatio<0.5) {
                                System.out.println(sourceContig+"\t"+sourceStart+
                                                   "\t"+sourceRefString+"\t"+sourceAltString+
                                                   "\t"+sourceRefForward+"\t"+sourceRefReverse+
                                                   "\t"+sourceAltForward+"\t"+sourceAltReverse+
                                                   "\t"+geneID+"\t"+maxAltRatio);
                            }
                        }
                    }
                }
            }
        }
    }

}
