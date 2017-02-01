package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;
import java.text.DecimalFormat;

/**
 * Combine SNPs which belong to the same exonic transcript(s).
 */
public class SNPCombine {

    static DecimalFormat df = new DecimalFormat("#.0");

    public static void main(String[] args) {

        if (args.length!=1) {
            System.out.println("Usage: SNPCombine <SNP file>");
            System.out.println("Example: SNPCombine VT1.snp");
            System.exit(1);
        }
        
        try {
            
            String fileName = args[0];

            File file = new File(fileName);
            if (!file.exists()) {
                System.err.println("Error: "+file.getName()+" does not exist.");
                System.exit(1);
            }

            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line = reader.readLine(); // header line
            String chr = "";
            String transcript = "";
            int pos1 = 0;
            int pos2 = 0;
            int homCount = 0;
            int hetCount = 0;
            int refDepth = 0;
            int altDepth = 0;
            double quality = 0.0;
            int mq = 0;
            while((line=reader.readLine())!=null) {
                SNPRecord rec = new SNPRecord(line);
                if (rec.isExonic()) {
                    if (rec.gene.equals(transcript)) {
                        // increment counter, etc.
                        pos2 = rec.pos2;
                        if (rec.isHomozygous()) homCount++;
                        if (rec.isHeterozygous()) hetCount++;
                        quality += rec.quality;
                        mq += rec.mq;
                        refDepth += rec.refDepth;
                        altDepth += rec.altDepth;
                    } else {
                        // output old stats if pure het or hom
                        boolean both = homCount>0 && hetCount>0;
                        if (mq>0 && !both && !chr.startsWith("B")) {
                            int count = homCount + hetCount;
                            String genotype = "hom";
                            if (homCount==0) genotype = "het";
                            SNPRecord out = new SNPRecord("exonic", transcript, chr, pos1, pos2, 'N', 'N', genotype, quality/count, mq/count, refDepth, altDepth);
                            System.out.println(out.toString());
                        }
                        // reset values for this SNP
                        chr = rec.chr;
                        transcript = rec.gene;
                        pos1 = rec.pos1;
                        pos2 = rec.pos2;
                        if (rec.isHomozygous()) { homCount = 1; hetCount=0; }
                        if (rec.isHeterozygous()) { hetCount = 1; homCount=0; }
                        refDepth = rec.refDepth;
                        altDepth = rec.altDepth;
                        quality = rec.quality;
                        mq = rec.mq;
                    }
                }
            }
            
        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
