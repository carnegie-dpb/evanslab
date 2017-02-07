package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.LinkedHashMap;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Slice out a specific range of locations and samples with the specified zygotes and display in a stacked list of sequences.
 * Useful for seeing what the exact bases are for, say, the homozygous calls of a sample subgroup.
 */
public class VCFSlicer {

    static DecimalFormat lz = new DecimalFormat("000000000");
    
    public static void main(String[] args) {

        if (args.length!=6) {
            System.err.println("Usage: VCFSlicer <REF fasta> <VCF file> <Samples file> contig start end");
            System.exit(1);
        }

        // input parameters
        // NOTE: vcf index will be assumed to be vcfFilename+".tbi" by the htsjdk routine
        String refFilename = args[0];
        String vcfFilename = args[1];
        String samplesFilename = args[2];
        String contig = args[3];
        int start = Integer.parseInt(args[4]);
        int end = Integer.parseInt(args[5]);

        try {

            File refFile = new File(refFilename);
            if (!refFile.exists()) {
                System.err.println("Error: file "+refFile.getName()+" does not exist.");
                System.exit(1);
            }

            File samplesFile = new File(samplesFilename);
            if (!samplesFile.exists()) {
                System.err.println("Error: file "+samplesFile.getName()+" does not exist.");
                System.exit(1);
            }
            
            File vcfFile = new File(vcfFilename);
            if (!vcfFile.exists()) {
                System.err.println("Error: file "+vcfFile.getName()+" does not exist.");
                System.exit(1);
            }

            // load desired samples as keys in a map of sequences 
            Map<String,String> sequenceMap = new LinkedHashMap<String,String>();
            BufferedReader samplesReader = new BufferedReader(new FileReader(samplesFile));
            String line = null;
            while ((line=samplesReader.readLine())!=null) {
                sequenceMap.put(line, "");
            }

            // get the reference contig sequence
            IndexedFastaSequenceFile refSequenceFile = new IndexedFastaSequenceFile(refFile);
            ReferenceSequence refSequence = refSequenceFile.getSubsequenceAt(contig, start, end);

            // load VCF file
            VCFFileReader vcfReader = new VCFFileReader(vcfFile, true);
            VCFHeader vcfHeader = vcfReader.getFileHeader();

            // spin through the position range
            for (int pos=start; pos<=end; pos++) {
            
                Iterator<VariantContext> vcIterator = vcfReader.query(contig, pos, pos);
                if (vcIterator.hasNext()) {

                    VariantContext vc = vcIterator.next();
                    
                    Allele ref = vc.getReference();
                    
                    for (String sample : sequenceMap.keySet()) {
                        
                        String sequence = sequenceMap.get(sample);
                        Genotype g = vc.getGenotype(sample);
                        if (g.isHet()) {
                            sequence += g.getAllele(1).getBaseString().toLowerCase(); // 0=ref, 1=alt
                        } else if (g.isHom()) {
                            sequence += g.getAllele(0).getBaseString().toUpperCase(); // 0=alt
                        } else {
                            sequence += "-";
                        }
                        sequenceMap.put(sample, sequence);
                    }
                    
                } else {

                    for (String sample : sequenceMap.keySet()) {
                        String sequence = sequenceMap.get(sample);
                        sequence += "-";
                        sequenceMap.put(sample, sequence);
                    }

                }

            }

            // now output it all
            System.out.println("Reference:\t"+refFilename);
            System.out.println("VCF file:\t"+vcfFilename);
            System.out.println("Samples:\t"+samplesFilename);
            System.out.println("Range:\t"+contig+":"+start+"-"+end);
            System.out.println("");
            System.out.println(contig+"\t"+start);
            System.out.println("\t"+"|");
            System.out.println("REF\t"+refSequence.getBaseString());
            for (String sample : sequenceMap.keySet()) {
                System.out.println(sample+"\t"+sequenceMap.get(sample));
            }
            
        } catch (Exception ex) {

            System.err.println(contig+":"+start+"-"+end);
            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
