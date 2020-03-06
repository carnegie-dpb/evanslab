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
 * This routine focuses on HET calls, which are defined by |log2(REF/ALT)| < 1.
 *
 * @author Sam Hokin
 */
public class VCFSegregation {

    static double LOG2 = Math.log(2.0);

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
        System.out.println("contig\tstart\tREF\tALT\tref1\talt1\tref2\talt2\tp");

        VCFFileReader reader = new VCFFileReader(new File(vcfFilename));
        for (VariantContext vc : reader) {
            String id = vc.getID();
            String source = vc.getSource();
            String contig = vc.getContig();
            Set<String> sampleNames = vc.getSampleNames();
            if (!sampleNames.contains(sample1) || !sampleNames.contains(sample2)) {
                System.err.println("ERROR: VCF does not contain both sample names "+sample1+" and "+sample2+".");
                System.exit(1);
            }
            int start = vc.getStart();
            Allele ref = vc.getReference();
            List<Allele> alts = vc.getAlternateAlleles();
            // use only single-ALT SNPs
            if (ref.getBaseString().length()==1 && alts.size()==1 && alts.get(0).getBaseString().length()==1) {
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
                    int ref1F = dp4List.get(0);
                    int ref1R = dp4List.get(1);
                    int alt1F = dp4List.get(2);
                    int alt1R = dp4List.get(3);
                    int ref2F = dp4List.get(4);
                    int ref2R = dp4List.get(5);
                    int alt2F = dp4List.get(6);
                    int alt2R = dp4List.get(7);
                    // if site has forward counts it must have reverse counts for same allele
                    boolean ref1ok = (ref1F==0 && ref1R==0) || (ref1F>0 && ref1R>0);
                    boolean alt1ok = (alt1F==0 && alt1R==0) || (alt1F>0 && alt1R>0);
                    boolean ref2ok = (ref2F==0 && ref2R==0) || (ref2F>0 && ref2R>0);
                    boolean alt2ok = (alt2F==0 && alt2R==0) || (alt2F>0 && alt2R>0);
                    if (ref1ok && alt1ok && ref2ok && alt2ok) {
                        int ref1 = ref1F + ref1R;
                        int alt1 = alt1F + alt1R;
                        int ref2 = ref2F + ref2R;
                        int alt2 = alt2F + alt2R;
                        // log2(ALT/REF) ratio is useful for downstream filtering
                        double ratio1 = (double)alt1 / (double)ref1;
                        double ratio2 = (double)alt2 / (double)ref2;
                        double p = fisherExact.getTwoTailedP(ref1, alt1, ref2, alt2);
                        System.out.println(contig+"\t"+start+"\t"+ref.getBaseString()+"\t"+altString+"\t"+
                                           ref1+"\t"+alt1+"\t"+ref2+"\t"+alt2+"\t"+p);
                    }
                }
            }
        }
    }
}
