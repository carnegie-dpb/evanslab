package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;

/**
 * Combine a bunch of SNP files into a single CSV for input into R/qtl as a csvr file.
 */
public class SNPToCSV {

    public static void main(String[] args) {

        if (args.length!=6) {
            System.out.println("Usage: SNPToCSV <SNP-dir> <prefix> <suffix> <minSample> <maxSample> <maxMissing>");
            System.out.println("Example: SNPToCSV ~/snps VT avinput.variant_function.snp 1 96 10");
            System.exit(1);
        }
        
        try {
            
            String snpDirName = args[0];
            String prefix = args[1];
            String suffix = args[2];
            int nMin = Integer.parseInt(args[3]);
            int nMax = Integer.parseInt(args[4]);
            int maxMissing = Integer.parseInt(args[5]);

            File dir = new File(snpDirName);
            if (!dir.isDirectory()) {
                System.err.println("Error: "+dir.getName()+" is not a directory.");
                System.exit(1);
            }

            // store all the individual sample data in this map for CSV output
            Map<Integer,Map<String,Character>> sampleMap = new TreeMap<Integer,Map<String,Character>>();

            // store all the SNP keys in a master list
            Map<String,String> keyMap = new TreeMap<String,String>();

            // spin through the SNP files and load data into geno maps
            for (int i=nMin; i<=nMax; i++) {
                String fileName = prefix+String.valueOf(i)+"."+suffix;
                File file = new File(dir, fileName);
                if (file.exists()) {
                    System.err.println(file.getName());
                    Map<String,Character> genoMap = new TreeMap<String,Character>();
                    BufferedReader reader = new BufferedReader(new FileReader(file));
                    String line = reader.readLine(); // header line
                    while((line=reader.readLine())!=null) {
                        SNPRecord rec = new SNPRecord(line);
                        String key = rec.chr+"|"+String.format("%10d", rec.pos1);
                        String gene = rec.gene.replaceAll("transcript:","").replaceAll(",",";");
                        if (rec.isHomozygous()) {
                            genoMap.put(key, 'B');
                            keyMap.put(key,gene);
                        } else if (rec.isHeterozygous()) {
                            genoMap.put(key, 'H');
                            keyMap.put(key,gene);
                        }
                    }
                    // key sampleMap by sample number
                    sampleMap.put(i, genoMap);
                }
            }

            // sample header
            System.out.print("PHENOTYPE"+","+"");
            for (Integer id : sampleMap.keySet()) {
                System.out.print(","+prefix+id);
            }
            System.out.println("");
            
            // spit out the loci for which we have a full complement of satisfactory calls
            for (String key : keyMap.keySet()) {
                String[] parts = key.split("\\|");
                String chr = parts[0];
                String gene = keyMap.get(key);
                // count up the missing genotypes as well as loading into a list for output
                int missing = 0;
                List<Character> genotypeList = new ArrayList<Character>();
                for (Integer id : sampleMap.keySet()) {
                    Map<String,Character> genoMap = sampleMap.get(id);
                    if (genoMap.containsKey(key)) {
                        genotypeList.add(genoMap.get(key));
                    } else {
                        genotypeList.add('-');
                        missing++;
                    }
                }
                // output a line if we've got enough genotype calls
                if (missing<=maxMissing) {
                    System.out.print(gene+","+chr);
                    for (Character genotype : genotypeList) {
                        System.out.print(","+genotype);
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
