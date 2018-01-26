package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Run through a combined VCF file, checking each non-called sample at each SNP to see if there are reference reads at that location.
 */
public class VCFRefCounter {

    public static void main(String[] args) {

        if (args.length!=3) {
            System.err.println("Usage: VCFRefCounter <VCF file> <BAM dir> <BAM suffix>");
            System.exit(1);
        }

        // input parameters
        String vcfFile = args[0];
        String bamDir = args[1];
        String bamSuffix = args[2];

        try {
            
            Map<String,SamReader> samReaders = new HashMap<String,SamReader>();
            List<String> samples = null;
            
            boolean first = true;
            VCFFileReader reader = new VCFFileReader(new File(vcfFile));
            Iterator<VariantContext> vcIterator = reader.iterator();
                        
            // loop over the variant lines
            while (vcIterator.hasNext()) {
                VariantContext vc = vcIterator.next();
                if (first) {
                    // load the samples and the SamReaders for all the BAM files
                    first = false;
                    samples = vc.getSampleNamesOrderedByName();
                    SamReaderFactory srf = SamReaderFactory.makeDefault();
                    for (String sample : samples) {
                        String bamFile = sample+"."+bamSuffix;
                        File samFile = new File(bamDir, bamFile);
                        if (samFile.exists()) {
                            SamReader sr = srf.open(samFile);
                            samReaders.put(sample, sr);
                            System.err.println("Loaded "+sr.getResourceDescription());
                        }
                    }
                }
                        
                // get the variants record
                String contig = vc.getContig();
                int start = vc.getStart();
                GenotypesContext gc = vc.getGenotypes();
                
                System.out.print(contig+"\t"+start);
                for (String sample : samples) {
                    Genotype g = gc.get(sample);
                    if (g.isNoCall()) {
                        SamReader sr = samReaders.get(sample);
                        if (sr!=null) {
                            SAMRecordIterator sri = sr.query(contig, start, start, false);
                            int count = 0;
                            double meanQual = 0.0;
                            while (sri.hasNext()) {
                                SAMRecord rec = sri.next();
                                count++;
                                meanQual += (double)rec.getMappingQuality();
                            }
                            if (count>0) {
                                meanQual = meanQual/count;
                                System.out.print("\tR:"+count+":"+Math.round(meanQual));
                            } else {
                                System.out.print("\tNC");
                            }
                            sri.close();
                        } else {
                            System.out.print("\tNB");
                        }
                    } else if (g.isHet()) {
                        System.out.print("\tHet");
                    } else if (g.isHom()) {
                        System.out.print("\tHom");
                    }
                }
                System.out.println("");
                
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }

}
