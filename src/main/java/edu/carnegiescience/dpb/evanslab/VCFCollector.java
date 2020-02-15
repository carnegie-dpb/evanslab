package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
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

import org.mskcc.cbio.portal.stats.FisherExact;

/**
 * Collect REF vs HET counts across the samples in a VCF file according to their nontwin/twin condition.
 */
public class VCFCollector {

    static final String[] nts = {"VT1", "VT2", "VT4", "VT5", "VT6", "VT7", "VT8", "VT9",
                                      "VT10", "VT11", "VT12", "VT13", "VT14", "VT15", "VT16", "VT17", "VT18", "VT19",
                                      "VT20", "VT22", "VT23", "VT24", "VT26", "VT27", "VT28", "VT29",
                                      "VT30", "VT31", "VT32", "VT33", "VT34", "VT35", "VT36", "VT37", "VT38", "VT39",
                                      "VT40", "VT41", "VT42", "VT43", "VT44", "VT45", "VT46", "VT47", "VT48"};

    static final String[] tws = {"VT49",
                                   "VT50", "VT51", "VT52", "VT53", "VT54", "VT55", "VT56", "VT57", "VT58", "VT59",
                                   "VT60", "VT61", "VT62", "VT63", "VT64", "VT65", "VT66", "VT67", "VT68", "VT69",
                                   "VT70", "VT71", "VT72", "VT73", "VT74", "VT75", "VT76", "VT77", "VT78", "VT79",
                                   "VT80", "VT81", "VT82", "VT83", "VT84", "VT85", "VT86", "VT87", "VT88", "VT89",
                                   "VT90", "VT91", "VT92", "VT93", "VT94", "VT95", "VT96"};

    static Set<String> ntSet = new HashSet(Arrays.asList(nts));
    static Set<String> twSet = new HashSet(Arrays.asList(tws));

    public static void main(String[] args) {

        if (args.length!=1) {
            System.err.println("Usage: VCFCollector <VCF.gz file>");
            System.exit(1);
        }

        File vcfFile = new File(args[0]);

        int maxSize = nts.length + tws.length;
        FisherExact fisherExact = new FisherExact(maxSize);
        
        // header
        System.out.println("chr\tpos\tref\talt\t"+
                           "twRef\ttwHet\ttwHom\t" +
                           "ntRef\tntHet\tntHom\t" +
                           "p\tmlog10p\tOR");

        VCFFileReader reader = new VCFFileReader(vcfFile, false);
        VCFHeader header = reader.getFileHeader();
        List<String> samples = header.getSampleNamesInOrder();

        Iterator<VariantContext> vcIterator = reader.iterator();
        while (vcIterator.hasNext()) {
            VariantContext vc = vcIterator.next();
            String contig = vc.getContig();
            int pos = vc.getStart();
            Allele ref = vc.getReference();
            List<Allele> alts = vc.getAlternateAlleles();
            // limit ourselves to single variant SNPs
            if (alts.size()==1) {
                // TWINS
                VariantContext twContext = vc.subContextFromSamples(twSet);
                GenotypesContext twGC = vc.getGenotypes(twSet);
                int twRefCount = twContext.getHomRefCount();
                int twHetCount = twContext.getHetCount();
                int twHomCount = twContext.getHomVarCount();
                // NONTWINS
                VariantContext ntContext = vc.subContextFromSamples(ntSet);
                GenotypesContext ntGC = vc.getGenotypes(ntSet);
                int ntRefCount = ntContext.getHomRefCount();
                int ntHetCount = ntContext.getHetCount();
                int ntHomCount = ntContext.getHomVarCount();
                // STATS
                double p = fisherExact.getP(twRefCount, twHetCount, ntRefCount, ntHetCount);
                double mlog10p = -Math.log10(p);
                double or = (double)twHetCount*(double)ntRefCount / ( (double)ntHetCount*(double)twRefCount);
                // OUTPUT
                System.out.println(contig+"\t"+pos+"\t"+ref+"\t"+alts.get(0)+"\t"+
                                   twRefCount+"\t"+twHetCount+"\t"+twHomCount+"\t"+
                                   ntRefCount+"\t"+ntHetCount+"\t"+ntHomCount+"\t"+
                                   p+"\t"+mlog10p+"\t"+or);
            }
        }
    }
}
