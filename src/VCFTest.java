package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Map;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Test various VCF classes and methods.
 */
public class VCFTest {

    public static void main(String[] args) {

        if (args.length>2) {
            System.err.println("Usage: VCFTest <VCF-file> [index-file]");
            System.exit(1);
        }

        try {
            
            File vcfFile = new File(args[0]);
            
            VCFFileReader reader;
            if (args.length==2) {
                File indexFile = new File(args[1]);
                reader = new VCFFileReader(vcfFile, indexFile, true);
            } else {
                reader = new VCFFileReader(vcfFile, false);
            }
            VCFHeader header = reader.getFileHeader();
            System.out.println(header.toString());

            Iterator<VariantContext> vcIterator = reader.iterator();
            while (vcIterator.hasNext()) {
                VariantContext vc = vcIterator.next();
                Allele ref = vc.getReference();
                if (vc.getStart()==vc.getEnd()) {
                    System.out.print(vc.getContig()+":"+vc.getStart()+"\t"+vc.getLog10PError()+"\t"+ref.getBaseString()+"->");
                } else {
                    System.out.print(vc.getContig()+":"+vc.getStart()+"-"+vc.getEnd()+":"+ref.getBaseString()+"->");
                }
                List<Allele> alts = vc.getAlternateAlleles();
                for (Allele alt : alts) System.out.print(alt.getBaseString()+" ");
                if (alts.size()>1) System.out.print("*");
                System.out.println("");
            }

            
            // SAMSequenceDictionary seqDictionary = VCFFileReader.getSequenceDictionary(vcfFile);
            // System.out.println("Sequence Dictionary: size="+seqDictionary.size()+"; ref length="+seqDictionary.getReferenceLength());

            // List<SAMSequenceRecord> seqRecordList = seqDictionary.getSequences();
            // for (SAMSequenceRecord rec : seqRecordList) {
            //     Set<Map.Entry<String,String>> attributes = rec.getAttributes();
            //     for (Map.Entry<String,String> map : attributes) {
            //         System.out.println(map.getKey()+":"+map.getValue());
            //     }
            // }

            // IntervalList intervalList = VCFFileReader.fromVcf(vcfFile);
            // System.out.println("Interval List: size="+intervalList.size()+"; base count="+intervalList.getBaseCount());
            // for (Interval interval : intervalList) {
            //     System.out.println(interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd()+":"+interval.toString());
            // }

            reader.close();

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
