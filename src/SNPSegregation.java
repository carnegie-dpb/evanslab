package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Map;
import java.util.TreeMap;

/**
 * Run through all the SNP files in the given directory and accumulate segregation statistics.
 */
public class SNPSegregation {

    public static void main(String[] args) {

        if (args.length!=10) {
            System.out.println("Usage: SNPSegregation hom|het|both <SNP-dir> <prefix> <suffix> <Amin> <Amax> <Bmin> <Bmax> <minDepth> <minLogSeg>");
            System.out.println("Example: SNPSegregation ~/snps VT avinput.variant_function.snp 1 48 49 96 1000 1.0");
            System.exit(1);
        }
        
        try {
            
            File dir = new File(args[1]);
            if (!dir.isDirectory()) {
                System.err.println("Error: "+args[0]+" is not a directory.");
                System.exit(1);
            }

            boolean hetOnly = args[0].equals("het");
            boolean homOnly = args[0].equals("hom");
            boolean both = !hetOnly && !homOnly;
            String prefix = args[2];
            String suffix = args[3];
            int aMin = Integer.parseInt(args[4]);
            int aMax = Integer.parseInt(args[5]);
            int bMin = Integer.parseInt(args[6]);
            int bMax = Integer.parseInt(args[7]);
            int minDepth = Integer.parseInt(args[8]);
            double minLogSeg = Double.parseDouble(args[9]);

            Map<String,Integer> aMap = new TreeMap<String,Integer>();
            Map<String,Integer> bMap = new TreeMap<String,Integer>();

            for (int i=aMin; i<=bMax; i++) {

                boolean isA = (aMin<=i && i<=aMax);
                boolean isB = (bMin<=i && i<=bMax);

                if (isA || isB) {

                    String fileName = prefix+String.valueOf(i)+"."+suffix;
                    File file = new File(dir, fileName);
                    if (file.exists()) {

                        BufferedReader reader = new BufferedReader(new FileReader(file));
                        String line = reader.readLine(); // header line
                        while((line=reader.readLine())!=null) {

                            SNPRecord rec = new SNPRecord(line);
                            String key = rec.chr+"\t"+rec.pos1+"\t"+rec.alt; // keep track of alleles separately

                            if ( (hetOnly && rec.genotype.equals("het")) || (homOnly && rec.genotype.equals("hom")) || both) {
                            
                                if (isA) {
                                    
                                    if (aMap.containsKey(key)) {
                                        // increment count
                                        int oldDepth = aMap.get(key);
                                        aMap.put(key, oldDepth+rec.altDepth);
                                    } else {
                                        // new SNP, don't care whether hom or het
                                        aMap.put(key, rec.altDepth);
                                    }
                                    
                                } else {
                                    
                                    if (bMap.containsKey(key)) {
                                        // increment count
                                        int oldDepth = bMap.get(key);
                                        bMap.put(key, oldDepth+rec.altDepth);
                                    } else {
                                        // new SNP, don't care whether hom or het
                                        bMap.put(key, rec.altDepth);
                                    }

                                }
                                
                            }

                        }
                        
                    }
                    
                }
                
            }

            // header line
            System.out.println("contig\tposition\talt\tlog10seg");

            // now determine segregated loci with SNPs in both groups and print out
            for (String key : aMap.keySet()) {
                if (bMap.containsKey(key)) {
                    int a = aMap.get(key);
                    int b = bMap.get(key);
                    if (a>minDepth || b>minDepth) {
                        double logSeg = Math.log10((double) b / (double) a);
                        if (Math.abs(logSeg)>minLogSeg) {
                            System.out.println(key+"\t"+logSeg);
                        }
                    }
                }
            }

            // now loci with SNPs in only A group, print out
            for (String key : aMap.keySet()) {
                if (!bMap.containsKey(key)) {
                    int a = aMap.get(key);
                    if (a>minDepth) {
                        System.out.println(key+"\tA ONLY:"+a);
                    }
                }
            }

            // now loci with SNPs in only B group, print out
            for (String key : bMap.keySet()) {
                if (!aMap.containsKey(key)) {
                    int b = bMap.get(key);
                    if (b>minDepth) {
                        System.out.println(key+"\tB ONLY:"+b);
                    }
                }
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
