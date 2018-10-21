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

import org.apache.commons.text.similarity.HammingDistance;

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
 * Useful for seeing what the exact bases are for the heterozygous and homozygous calls of a sample subgroup.
 *
 * @author Sam Hokin
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
            Map<String,String> sequenceMap = new LinkedHashMap<String,String>();         // store allele markers
            Map<String,String> fullSequenceMap = new LinkedHashMap<String,String>();     // actual base for calculation of distance
            BufferedReader samplesReader = new BufferedReader(new FileReader(samplesFile));
            String line = null;
            while ((line=samplesReader.readLine())!=null) {
                sequenceMap.put(line, "");
                fullSequenceMap.put(line, "");
            }

            // get the reference contig sequence
            IndexedFastaSequenceFile refSequenceFile = new IndexedFastaSequenceFile(refFile);
            ReferenceSequence refSequence = refSequenceFile.getSubsequenceAt(contig, start, end);
            String refSequenceString = refSequence.getBaseString();

            // load VCF file
            VCFFileReader vcfReader = new VCFFileReader(vcfFile, true);
            VCFHeader vcfHeader = vcfReader.getFileHeader();

            // spin through the position range
            for (int pos=start; pos<=end; pos++) {
            
                // query a single position at a time
                Iterator<VariantContext> vcIterator = vcfReader.query(contig, pos, pos);
                if (vcIterator.hasNext()) {

                    // this position has calls
                    VariantContext vc = vcIterator.next();
                    Allele ref = vc.getReference();
                    for (String sample : sequenceMap.keySet()) {
                        String sequence = sequenceMap.get(sample);
                        String fullSequence = fullSequenceMap.get(sample);

                        Genotype g = vc.getGenotype(sample);
                        if (g.isCalled()) {
                            if (g.isHom()) {
                                sequence += g.getAllele(0).getBaseString().toUpperCase();     // mark with upper case alt call
                                fullSequence += g.getAllele(0).getBaseString().toUpperCase();
                            } else if (g.isHet()) {
                                sequence += g.getAllele(1).getBaseString().toLowerCase();     // mark with lower case alt call
                                fullSequence += g.getAllele(1).getBaseString().toLowerCase();
                            } else {
                                sequence += ".";                                              // mark as not called
                                fullSequence += refSequenceString.charAt(pos-start);          // assume reference for distance calc
                            }
                        } else {
                            sequence += ".";                                              // mark as not called
                            fullSequence += refSequenceString.charAt(pos-start);          // assume reference for distance calc
                        }
                        // update the maps
                        sequenceMap.put(sample, sequence);
                        fullSequenceMap.put(sample, fullSequence);
                    }
                    
                } else {

                    // this position does NOT have calls
                    for (String sample : sequenceMap.keySet()) {
                        String sequence = sequenceMap.get(sample);
                        String fullSequence = fullSequenceMap.get(sample);
                        sequence += ".";                                                  // mark as not called
                        fullSequence += refSequenceString.charAt(pos-start);              // assume reference for distance calc
                        // update the maps
                        sequenceMap.put(sample, sequence);
                        fullSequenceMap.put(sample, fullSequence);
                    }

                }

            }

            // need to instantiate this for calculating distance from reference
            HammingDistance ld = new HammingDistance();
            
            // now output it all
            System.out.println("Reference:\t"+refFilename);
            System.out.println("VCF file:\t"+vcfFilename);
            System.out.println("Samples:\t"+samplesFilename);
            System.out.println("Range:\t"+contig+":"+start+"-"+end);
            System.out.println("");
            System.out.println(contig+"\t"+start);
            System.out.println("\t"+"|");
            System.out.println("REF\t"+refSequence.getBaseString()+" 0");
            for (String sample : sequenceMap.keySet()) {
                Integer distance = ld.apply(refSequenceString, fullSequenceMap.get(sample));
                System.out.println(sample+"\t"+sequenceMap.get(sample)+" "+distance);
            }
            
        } catch (Exception ex) {

            System.err.println(contig+":"+start+"-"+end);
            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
