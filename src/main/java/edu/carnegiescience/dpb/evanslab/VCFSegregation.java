package edu.carnegiescience.dpb.evanslab;

import java.io.File;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import org.mskcc.cbio.portal.stats.FisherExact;

/**
 * Loads a VCF file and computes segregation between the (first) two samples using Fisher's exact test.
 *
 * @author Sam Hokin
 */
public class VCFSegregation {

    /**
     * Main class outputs a tab-delimited summary of segregation per locus that has calls for both samples.
     */
    public static void main(String[] args) {
        if (args.length!=4) {
            System.out.println("Usage VCFSegregation <vcf-file> <sample1> <sample2> <maxSize>");
            System.exit(0);
        }

        String vcfFilename = args[0];
        String sample1 = args[1];
        String sample2 = args[2];
        int maxSize = Integer.parseInt(args[3]);

        FisherExact fisherExact = new FisherExact(maxSize);

        // output heading
        System.out.println("contig\tstart\tREF\tALT\ta\tb\tc\td\tsize\tp\tmlog10p\tsignif");

        VCFFileReader reader = new VCFFileReader(new File(vcfFilename));
        for (VariantContext vc : reader) {
            String id = vc.getID();
            String source = vc.getSource();
            String contig = vc.getContig();
            Set<String> sampleNames = vc.getSampleNames();
            if (sampleNames.contains(sample1) && sampleNames.contains(sample2)) {
                int start = vc.getStart();
                Allele ref = vc.getReference();
                List<Allele> alts = vc.getAlternateAlleles();
                // use only single-ALT calls
                if (alts.size()==1) {
                    String altString = alts.get(0).getBaseString();
                    // DP4 is deprecated in favor of ADF, ADR which are restricted to "high-quality" calls
                    // but we seem to need DP4 for htseq routines
                    //     RF1 RR1 AF1 AR1  RF2 RR2 AF2 AR2            ADF1:ADR1   ADF2:ADR2
                    // DP4=13, 17, 7,  0,   12, 29, 12, 0    ADF:ADR   13,6:17,0  12,11:29,0
                    // ADF .. Total allelic depths on the forward strand (Number=R,Type=Integer)
                    // ADR .. Total allelic depths on the reverse strand (Number=R,Type=Integer)
                    List<Integer> dp4List = vc.getAttributeAsIntList("DP4", 0); // FORWARD REF1,ALT1,REF2,ALT2
                    // restrict output to sites with full calls for both samples
                    if (dp4List.size()==8) {
                        int size = dp4List.get(0) + dp4List.get(1) + dp4List.get(2) + dp4List.get(3) +
                            dp4List.get(0) + dp4List.get(1) + dp4List.get(2) + dp4List.get(3);
                        int a = dp4List.get(0)+dp4List.get(1); // REF1
                        int b = dp4List.get(2)+dp4List.get(3); // ALT1
                        int c = dp4List.get(4)+dp4List.get(5); // REF2
                        int d = dp4List.get(6)+dp4List.get(7); // ALT2
                        double p = fisherExact.getP(a, b, c, d);
                        double minusLog10p = -Math.log10(p);
                        boolean significant = (p<0.01);
                        System.out.println(contig+"\t"+start+"\t"+ref.getBaseString()+"\t"+altString+"\t"+
                                           a+"\t"+b+"\t"+c+"\t"+d+"\t"+size+"\t"+p+"\t"+minusLog10p+"\t"+significant);
                    }
                }
            }
        }
    }
}
