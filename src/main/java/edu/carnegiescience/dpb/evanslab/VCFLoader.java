package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.util.List;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class VCFLoader {


    public static void main(String[] args) {

        if (args.length!=1) {
            System.out.println("Usage VCFLoader <vcf-file>");
            System.exit(0);
        }

        String vcfFilename = args[0];
        File vcfFile = new File(vcfFilename);

        VCFFileReader reader = new VCFFileReader(vcfFile);

        for (VariantContext vc : reader) {
            if (vc.isSNP()) {
                // values
                String id = vc.getID();
                String source = vc.getSource();
                String contig = vc.getContig();
                List<String> sampleNames = vc.getSampleNamesOrderedByName();
                int start = vc.getStart();
                Allele ref = vc.getReference();
                List<Allele> alts = vc.getAlternateAlleles();
                List<Integer> dp4List = vc.getAttributeAsIntList("DP4", 0);
                // output
                String altString = "";
                for (Allele alt : alts) {
                    if (altString.length()>0) altString += ",";
                    altString += alt.getBaseString();
                }                
                System.out.println(contig+"\t"+start+"\t"+ref.getBaseString()+"\t"+altString+"\t"+dp4List.get(0)+"\t"+dp4List.get(1)+"\t"+dp4List.get(2)+"\t"+dp4List.get(3));
            }
        }
    }

}
