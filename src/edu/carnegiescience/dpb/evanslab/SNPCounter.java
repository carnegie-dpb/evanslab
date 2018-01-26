package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;

/**
 * Run through all the SNP files and count the total number of Het vs Hom loci per sample.
 */
public class SNPCounter {

    public static void main(String[] args) {

        if (args.length!=6) {
            System.out.println("Usage: SNPCounter <SNP-dir> <prefix> <suffix> <minSample> <maxSample> <minDepth>");
            System.out.println("Example: SNPCounter ~/snps VT avinput.variant_function.snp 1 96 5");
            System.exit(1);
        }
        
        try {
            
            String snpDirName = args[0];
            String prefix = args[1];
            String suffix = args[2];
            int nMin = Integer.parseInt(args[3]);
            int nMax = Integer.parseInt(args[4]);
            int minDepth = Integer.parseInt(args[5]);

            File dir = new File(snpDirName);
            if (!dir.isDirectory()) {
                System.err.println("Error: "+dir.getName()+" is not a directory.");
                System.exit(1);
            }

            // header
            System.out.println("Samp"+"\t"+"Hom"+"\t"+"Het"+"\t"+"Ref"+"\t"+"Alt");

            // spin through the SNP files
            for (int i=nMin; i<=nMax; i++) {
                String fileName = prefix+String.valueOf(i)+"."+suffix;
                File file = new File(dir, fileName);
                if (file.exists()) {
                    int homCount = 0;
                    int hetCount = 0;
                    int refReads = 0;
                    int altReads = 0;
                    BufferedReader reader = new BufferedReader(new FileReader(file));
                    String line = reader.readLine(); // header line
                    while((line=reader.readLine())!=null) {
                        SNPRecord rec = new SNPRecord(line);
                        refReads += rec.refDepth;
                        altReads += rec.altDepth;
                        if (rec.isHomozygous() && rec.altDepth>=minDepth*2 && rec.refDepth==0) homCount++;
                        if (rec.isHeterozygous() && rec.altDepth>=minDepth && rec.refDepth>=minDepth) hetCount++;
                    }
                    System.out.println(i+"\t"+homCount+"\t"+hetCount+"\t"+refReads+"\t"+altReads);
                }
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
