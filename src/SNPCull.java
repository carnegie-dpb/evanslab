package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Map;
import java.util.TreeMap;

/**
 * Run through all the SNP files and output a CSV file containing only loci for which there are significant calls throughout.
 */
public class SNPCull {

    public static void main(String[] args) {

        if (args.length!=6) {
            System.out.println("Usage: SNPCull <SNP-dir> <prefix> <suffix> <minSample> <maxSample> <minDepth>");
            System.out.println("Example: SNPCull ~/snps VT avinput.variant_function.snp 1 96 5");
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
            // System.out.println("Samp"+"\t"+"Hom"+"\t"+"Het"+"\t"+"Ref"+"\t"+"Alt");

            Map<String,Integer> saveMap = new TreeMap<String,Integer>();
            Map<String,Map<String,Character>> sampleMap = new TreeMap<String,Map<String,Character>>();

            // the first file defines the maximal set of loci
            String baseName = prefix+String.valueOf(nMin)+"."+suffix;
            File baseFile = new File(dir, baseName);
            if (!baseFile.exists()) {
                System.err.println("Error: "+baseFile.getName()+" does not exist.");
                System.exit(1);
            }
            BufferedReader baseReader = new BufferedReader(new FileReader(baseFile));
            String baseLine = baseReader.readLine(); // header
            while ((baseLine=baseReader.readLine())!=null) {
                SNPRecord rec = new SNPRecord(baseLine);
                String key = rec.chr+":"+String.format("%10d", rec.pos1); // format preserves alpha ordering
                if ((rec.isReference() && rec.refDepth>=minDepth*2 && rec.altDepth==0) ||
                    (rec.isHomozygous() && rec.altDepth>=minDepth*2 && rec.refDepth==0) ||
                    (rec.isHeterozygous() && rec.altDepth>=minDepth && rec.refDepth>=minDepth)) {
                    saveMap.put(key, 0);
                }
            }

            // now spin through the SNP files and increment counter if locus contained in saveMap
            int samples = 0;
            for (int i=nMin; i<=nMax; i++) {
                String fileName = prefix+String.valueOf(i)+"."+suffix;
                File file = new File(dir, fileName);
                if (file.exists()) {
                    samples++;
                    System.err.println(file.getName());
                    Map<String,Character> genoMap = new TreeMap<String,Character>();
                    BufferedReader reader = new BufferedReader(new FileReader(file));
                    String line = reader.readLine(); // header line
                    while((line=reader.readLine())!=null) {
                        SNPRecord rec = new SNPRecord(line);
                        String key = rec.chr+":"+String.format("%10d", rec.pos1); // format preserves alpha ordering
                        if (saveMap.containsKey(key)) {
                            if (rec.isReference() && rec.refDepth>=minDepth*2 && rec.altDepth==0) {
                                int count = saveMap.get(key) + 1;
                                saveMap.put(key, count);
                                genoMap.put(key, 'A');
                            } else if (rec.isHomozygous() && rec.altDepth>=minDepth*2 && rec.refDepth==0) {
                                int count = saveMap.get(key) + 1;
                                saveMap.put(key, count);
                                genoMap.put(key, 'B');
                            } else if (rec.isHeterozygous() && rec.altDepth>=minDepth && rec.refDepth>=minDepth) {
                                int count = saveMap.get(key) + 1;
                                saveMap.put(key, count);
                                genoMap.put(key, 'H');
                            }
                        }
                    }
                    sampleMap.put(prefix+String.valueOf(i), genoMap);
                }
            }

            // sample header
            System.out.print("PHENOTYPE"+","+"");
            for (String sample : sampleMap.keySet()) {
                System.out.print(","+sample);
            }
            System.out.println("");
            
            // spit out the loci for which we have a full complement of satisfactory calls
            for (String key : saveMap.keySet()) {
                int count = saveMap.get(key);
                if (count==samples) {
                    String[] parts = key.split(":");
                    String chr = parts[0];
                    int pos = Integer.parseInt(parts[1].trim());
                    System.out.print(chr+":"+pos+","+chr);
                    for (String sample : sampleMap.keySet()) {
                        Map genoMap = sampleMap.get(sample);
                        System.out.print(","+genoMap.get(key));
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
