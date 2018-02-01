package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;

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
 * Display the VCF file values for a given contig and position.
 *
 * Genotype methods:
 *
 * isAvailable()  true if the type of this genotype is set.
 * isCalled()     true if this genotype is comprised of any alleles that are not no-calls (even if some are).
 * isHet()        true if we're het (observed alleles differ); if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
 * isHetNonRef()  true if we're het (observed alleles differ) and neither allele is reference; if the ploidy is less than 2 or if any alleles are no-calls, this method will return false.
 * isHom()        true if all observed alleles are the same (regardless of whether they are ref or alt); if any alleles are no-calls, this method will return false
 * isHomRef()     true if all observed alleles are ref; if any alleles are no-calls, this method will return false.
 * isHomVar()     true if all observed alleles are alt; if any alleles are no-calls, this method will return false.
 * isMixed()      true if this genotype is comprised of both calls and no-calls.
 * isNoCall()     true if this genotype is not actually a genotype but a "no call" (e.g. './.' in VCF); if any alleles are not no-calls (even if some are), this method will return false.
 *
 */
public class VCFDisplayer {

    static DecimalFormat lz = new DecimalFormat("000000000");
    
    public static void main(String[] args) throws FileNotFoundException {

        if (args.length!=4) {
            System.err.println("Usage: VCFDisplayer <REF fasta> <VCF file> contig position");
            System.exit(1);
        }

        // input parameters
        // NOTE: vcf index will be assumed to be vcfFilename+".tbi" by the htsjdk routine
        String refFilename = args[0];
        String vcfFilename = args[1];
        String contig = args[2];
        int position = Integer.parseInt(args[3]);
        
        File refFile = new File(refFilename);
        if (!refFile.exists()) {
            System.err.println("Error: file "+refFile.getName()+" does not exist.");
            System.exit(1);
        }
        
        File vcfFile = new File(vcfFilename);
        if (!vcfFile.exists()) {
            System.err.println("Error: file "+vcfFile.getName()+" does not exist.");
            System.exit(1);
        }
        
        // get the reference base
        IndexedFastaSequenceFile refSequenceFile = new IndexedFastaSequenceFile(refFile);
        ReferenceSequence refSequence = refSequenceFile.getSubsequenceAt(contig, position, position);
        String refSequenceBase = refSequence.getBaseString();
        
        // load VCF file
        VCFFileReader vcfReader = new VCFFileReader(vcfFile, true);

        // meta data
        VCFHeader vcfHeader = vcfReader.getFileHeader();
        List<String> sampleNames = vcfHeader.getSampleNamesInOrder();

        // query the requested position
        Iterator<VariantContext> vcIterator = vcfReader.query(contig, position, position);
        if (vcIterator.hasNext()) {

            VariantContext vc = vcIterator.next();
            System.out.println(vc.toStringWithoutGenotypes());
            Map<String,Object> attributeMap = vc.getAttributes();
            for (String key : attributeMap.keySet()) {
                System.out.println(key+":"+attributeMap.get(key));
            }
            System.out.println("PhredScaledQual="+vc.getPhredScaledQual());

            System.out.print("REF="+vc.getReference().getBaseString()+"\t");
            System.out.print("ALT=");
            for (Allele alt : vc.getAlternateAlleles()) {
                System.out.print(alt.getBaseString()+" ");
            }
            System.out.println("");
            System.out.println("SAMP\tGQ\tCalled\tMixed\tHom\tHomVar\tHomRef\tHet\tHetNonRef");

            for (String sample : sampleNames) {
                Genotype g = vc.getGenotype(sample);
                System.out.print(g.getSampleName()+"\t");
                System.out.print(g.getGQ()+"\t");
                System.out.print(g.isCalled()+"\t");
                System.out.print(g.isMixed()+"\t");
                System.out.print(g.isHom()+"\t");
                System.out.print(g.isHomVar()+"\t");
                System.out.print(g.isHomRef()+"\t");
                System.out.print(g.isHet()+"\t");
                System.out.print(g.isHetNonRef()+"\t");
                System.out.println("");
            }
            
        }
        
    }

}
