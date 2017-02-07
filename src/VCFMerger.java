package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Map;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Merge consecutive SNPs from a VCF file for which there is no LD amongst the individuals.
 */
public class VCFMerger {

    static DecimalFormat lz = new DecimalFormat("000000000");
    static DecimalFormat df = new DecimalFormat("00.00");
    
    public static void main(String[] args) {

        if (args.length!=9) {
            System.err.println("Usage: VCFMerger file.vcf tab|csv alleles|zygotes minDifference maxMissing minQ minAvgQ minHom minHet");
            System.exit(1);
        }

        String filename = args[0];
        boolean tabFormat = args[1].equals("tab");
        boolean csvFormat = args[1].equals("csv");
        boolean alleleOutput = args[2].equals("alleles");
        boolean zygoteOutput = args[2].equals("zygotes");
        double minDifference = Double.parseDouble(args[3]);
        int maxMissing = Integer.parseInt(args[4]);
        int minQ = Integer.parseInt(args[5]);
        double minAvgQ = Double.parseDouble(args[6]);
        int minHom = Integer.parseInt(args[7]);
        int minHet = Integer.parseInt(args[8]);

        try {
            
            File file = new File(filename);
            if (!file.exists()) {
                System.err.println("Error: file "+file.getName()+" does not exist.");
                System.exit(1);
            }

            // keep track of last output so we only have the first
            String oldOutput = "";
            String oldContig = "";
            int oldPos = 0;

            // sum differences, which is essentially cM length of the chromosome
            double totalDifference = 0.0;

            VCFFileReader reader = new VCFFileReader(file, false);
            VCFHeader header = reader.getFileHeader();
            List<String> samples = header.getSampleNamesInOrder();

            if (csvFormat) {
                System.out.print("PHENO,");
                for (String sample : samples) System.out.print(","+sample);
                System.out.println("");
            }
            
            Iterator<VariantContext> vcIterator = reader.iterator();
            while (vcIterator.hasNext()) {

                // add up individual calls
                int countRef = 0;
                int countHom = 0;
                int countHet = 0;
                int countMissing = 0;

                double avgQ = 0.0;
                int nQ = 0;

                VariantContext vc = vcIterator.next();
                String contig = vc.getContig();
                int pos = vc.getStart();

                String output = "";
                        
                // always start with first locus at each chromosome
                if (!contig.equals(oldContig)) {

                    if (tabFormat && totalDifference>0.0) {
                        System.out.println("\t\ttotal\t"+df.format(totalDifference));
                        System.out.println("");
                    }
                    
                    oldOutput = "";
                    oldContig = contig;
                    oldPos = -10000000;
                    totalDifference = 0.0;
                }

                Allele ref = vc.getReference();
                
                for (String sample : samples) {
                    Genotype g = vc.getGenotype(sample);
                    if (g.hasGQ()) {
                        nQ++;
                        avgQ += (double) g.getGQ();
                    }
                    if (g.hasGQ() && g.getGQ()<minQ) {
                        countMissing++;
                        output += "-";
                    } else if (g.isHomRef()) {
                        countRef++;
                        if (zygoteOutput) output += "A";
                        if (alleleOutput) output += ref.getBaseString().toUpperCase();
                    } else if (g.isHet()) {
                        countHet++;
                        if (zygoteOutput) output += "H";
                        if (alleleOutput) output += g.getAllele(1).getBaseString().toLowerCase(); // 0=ref 1=alt
                    } else if (g.isHom()) {
                        countHom++;
                        if (zygoteOutput) output += "B";
                        if (alleleOutput) output += g.getAllele(0).getBaseString().toUpperCase(); // 0=alt
                    } else {
                        countMissing++;
                        output += "-";
                    }
                }

                // average Q value
                if (nQ>0) avgQ = avgQ/nQ;
                
                // must meet filter criteria

                if (countHet>=minHet && countHom>=minHom && countMissing<=maxMissing && avgQ>=minAvgQ) {

                    String marker = contig+"_"+lz.format(pos);
                    
                    if (oldOutput.length()==0) {
                        // new chromosome
                        if (tabFormat) System.out.println(marker+"\t"+df.format(avgQ)+"\t\t"+countHet+"\t"+countHom+"\t"+output);
                        if (csvFormat) {
                            System.out.print(marker+","+contig);
                            char[] chars = output.toCharArray();
                            for (int i=0; i<chars.length; i++) System.out.print(","+chars[i]);
                            System.out.println("");
                        }
                        oldOutput = output;
                    } else {
                        // calculate difference between old and new output, units are essentially cM
                        double difference = 0.00;
                        char[] oldChars = oldOutput.toCharArray();
                        char[] newChars = output.toCharArray();
                        for (int i=0; i<oldChars.length; i++) {
                            if (oldChars[i]!=newChars[i] && oldChars[i]!='-' && newChars[i]!='-') difference += 1.0;
                        }
                        difference = difference/samples.size()*100;

                        if (difference>=minDifference) {
                            
                            totalDifference += difference;
                            if (tabFormat) System.out.println(marker+"\t"+df.format(avgQ)+"\t"+df.format(difference)+"\t"+countHet+"\t"+countHom+"\t"+output);
                            if (csvFormat) {
                                System.out.print(marker+","+contig);
                                char[] chars = output.toCharArray();
                                for (int i=0; i<chars.length; i++) System.out.print(","+chars[i]);
                                System.out.println("");
                            }
                            oldOutput = output;
                        }
                    }

                }
                
            }
            
        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
