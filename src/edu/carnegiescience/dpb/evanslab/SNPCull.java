package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.TreeMap;

/**
 * Run through all the SNP files and output a rotated CSV file (format="csvr") containing only loci for which there are significant calls throughout.
 */
public class SNPCull {

    static DecimalFormat cm = new DecimalFormat("#.000000");
    static DecimalFormat idx = new DecimalFormat("0000");

    static char NULL_CHAR = '-';

    public static void main(String[] args) {

        if (args.length!=7) {
            System.out.println("Usage: SNPCull <SNP-dir> <prefix> <suffix> <minSample> <maxSample> <minDepth> <minSamples>");
            System.out.println("Example: SNPCull ~/snps VT avinput.variant_function.snp 1 96 5 80");
            System.exit(1);
        }
        
        try {
            
            String snpDirName = args[0];
            String prefix = args[1];
            String suffix = args[2];
            int nMin = Integer.parseInt(args[3]);
            int nMax = Integer.parseInt(args[4]);
            int minDepth = Integer.parseInt(args[5]);
            int minSamples = Integer.parseInt(args[6]);

            File dir = new File(snpDirName);
            if (!dir.isDirectory()) {
                System.err.println("Error: "+dir.getName()+" is not a directory.");
                System.exit(1);
            }

            Map<String,Integer> saveMap = new TreeMap<String,Integer>();
            Map<String,Map<String,Character>> sampleMap = new TreeMap<String,Map<String,Character>>();

            // now spin through the SNP files and increment counter if locus contained in saveMap
            for (int i=nMin; i<=nMax; i++) {
                String fileName = prefix+String.valueOf(i)+"."+suffix;
                File file = new File(dir, fileName);
                if (file.exists()) {
                    System.err.println(file.getName());
                    Map<String,Character> genotypeMap = new TreeMap<String,Character>();
                    BufferedReader reader = new BufferedReader(new FileReader(file));
                    String line = reader.readLine(); // header line
                    while((line=reader.readLine())!=null) {
                        SNPRecord rec = new SNPRecord(line);
                        if (rec.isNumericChromosome()) {
                            String key = rec.chr+":"+String.format("%10d", rec.pos1); // format preserves alpha ordering
                            if (saveMap.containsKey(key)) {
                                // increment count
                                int count = saveMap.get(key) + 1;
                                saveMap.put(key, count);
                            } else {
                                // add new record
                                saveMap.put(key, 1);
                            }
                            if (rec.isReference() && rec.refDepth>=minDepth) genotypeMap.put(key, 'A');
                            if (rec.isHomozygous() && rec.altDepth>=minDepth) genotypeMap.put(key, 'B');
                            if (rec.isHeterozygous() && (rec.altDepth+rec.refDepth)>=minDepth) genotypeMap.put(key, 'H');
                        }
                    }
                    sampleMap.put(prefix+idx.format(i), genotypeMap);
                }
            }

            // output sample names
            System.out.print("PHENOTYPE,,");
            for (String sample : sampleMap.keySet()) {
                System.out.print(","+sample);
            }
            System.out.println("");

            // output a line for each locus; use 1 cM = 1 Mbp
            for (String key : saveMap.keySet()) {
                int count = saveMap.get(key);
                if (count>=minSamples) {
                    String[] parts = key.split(":");
                    String chr = parts[0];
                    int pos = Integer.parseInt(parts[1].trim());
                    String locus = chr+":"+pos;
                    System.out.print(locus+","+chr+","+cm.format((double)pos*1e-6));
                    for (String sample : sampleMap.keySet()) {
                        Character genotype = sampleMap.get(sample).get(key);
                        if (genotype==null) {
                            System.out.print(","+NULL_CHAR);
                        } else {
                            System.out.print(","+genotype);
                        }
                    }
                    System.out.println("");
                }
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
