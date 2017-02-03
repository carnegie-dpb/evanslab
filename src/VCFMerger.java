package edu.carnegiescience.dpb.evanslab;

import java.io.File;
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

    public static void main(String[] args) {

        if (args.length!=4) {
            System.err.println("Usage: VCFMerger <VCF file> minSep minHom minHet");
            System.exit(1);
        }

        String filename = args[0];
        int minSep = Integer.parseInt(args[1]);
        int minHom = Integer.parseInt(args[2]);
        int minHet = Integer.parseInt(args[3]);

        try {
            
            File file = new File(filename);
            if (!file.exists()) {
                System.err.println("Error: file "+file.getName()+" does not exist.");
                System.exit(1);
            }

            // keep track of last output so we only have the first
            String lastOutput = "";
            String lastContig = "";
            int lastPos = 0;

            VCFFileReader reader = new VCFFileReader(file, false);
            VCFHeader header = reader.getFileHeader();
            List<String> samples = header.getSampleNamesInOrder();
            
            Iterator<VariantContext> vcIterator = reader.iterator();
            while (vcIterator.hasNext()) {

                // keep track of Het and Hom calls, don't output if below minHom and minHet
                int countHom = 0;
                int countHet = 0;

                String output = "";
                        
                VariantContext vc = vcIterator.next();
                String contig = vc.getContig();
                int pos = vc.getStart();

                // always start with first locus at each chromosome
                if (!contig.equals(lastContig)) {
                    lastPos = -10000000;
                    lastContig = contig;
                }
                
                for (String sample : samples) {
                    Genotype g = vc.getGenotype(sample);
                    if (g.isHomRef()) {
                        output += "A";
                    } else if (g.isHet()) {
                        countHet++;
                        output += "H";
                    } else if (g.isHom()) {
                        countHom++;
                        output += "B";
                    } else {
                        output += "-";
                    }
                }
                
                if ((pos-lastPos)>=minSep && countHet>=minHet && countHom>=minHom) {
                    // substitute A for missing and see if they match up
                    String outputPlusA = output.replaceAll("-","A");
                    String lastOutputPlusA = lastOutput.replaceAll("-","A");
                    boolean matchWithA = outputPlusA.equals(lastOutputPlusA);
                    // substitute H for missing and see if they match up
                    String outputPlusH = output.replaceAll("-","H");
                    String lastOutputPlusH = lastOutput.replaceAll("-","H");
                    boolean matchWithH = outputPlusH.equals(lastOutputPlusH);
                    // substitute B for missing and see if they match up
                    String outputPlusB = output.replaceAll("-","B");
                    String lastOutputPlusB = lastOutput.replaceAll("-","B");
                    boolean matchWithB = outputPlusB.equals(lastOutputPlusB);
                    if (output.equals(lastOutput) || matchWithA || matchWithH || matchWithB) {
                        // do nothing, stick with previous lastOutput
                    } else {
                        System.out.println(vc.getContig()+":"+vc.getStart()+"\t"+output);
                        lastOutput = output;
                        lastPos = pos;
                    }
                }
                
            }
            
        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
