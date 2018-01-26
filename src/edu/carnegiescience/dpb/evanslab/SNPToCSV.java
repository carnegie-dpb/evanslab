package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.text.DecimalFormat;
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

    static DecimalFormat lz = new DecimalFormat("0000000000"); // leading zeros
    static DecimalFormat cmf = new DecimalFormat("#.00");    // cM
    
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

            // store all the genes as keys in a master list, with location as the value, updating pos1-pos2 as it expands
            Map<String,String> geneMap = new TreeMap<String,String>();

            // spin through the SNP files and load data into geno maps
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
                        String gene = rec.gene.replaceAll("transcript:","").replaceAll(",",";");
                        if (geneMap.containsKey(gene)) {
                            String loc = geneMap.get(gene);
                            String parts[] = loc.split(":");
                            String positions[] =  parts[1].split("-");
                            int pos1 = Integer.parseInt(positions[0]);
                            int pos2 = Integer.parseInt(positions[1]);
                            if (rec.pos1<pos1 || rec.pos2>pos2) {
                                // update region
                                pos1 = Math.min(pos1, rec.pos1);
                                pos2 = Math.max(pos2, rec.pos2);
                                loc = rec.chr+":"+lz.format(pos1)+"-"+lz.format(pos2);
                                geneMap.put(gene, loc);
                            }
                        } else {
                            String loc = rec.chr+":"+lz.format(rec.pos1)+"-"+lz.format(rec.pos2);
                            geneMap.put(gene, loc);
                        }
                        if (rec.isReference()) {
                            genotypeMap.put(gene, 'A');
                        } else if (rec.isHomozygous()) {
                            genotypeMap.put(gene, 'B');
                        } else if (rec.isHeterozygous()) {
                            genotypeMap.put(gene, 'H');
                        } else {
                            genotypeMap.put(gene, '-');
                        }
                    }
                    // key sampleMap by sample number
                    sampleMap.put(i, genotypeMap);
                }
            }

            // sample header
            System.out.print("PHENOTYPE"+","+"CHR"+","+"POS");
            for (Integer id : sampleMap.keySet()) {
                System.out.print(","+prefix+id);
            }
            System.out.println("");
            
            // spit out the loci for which we have a full complement of satisfactory calls
            for (String gene : geneMap.keySet()) {
                String loc = geneMap.get(gene);
                String[] parts = loc.split(":");
                String chr = parts[0];
                String[] positions = parts[1].split("-");
                int pos1 = Integer.parseInt(positions[0]);
                int pos2 = Integer.parseInt(positions[1]);
                // count up the missing genotypes as well as loading into a list for output
                int missing = 0;
                List<Character> genotypeList = new ArrayList<Character>();
                for (Integer id : sampleMap.keySet()) {
                    Map<String,Character> genotypeMap = sampleMap.get(id);
                    if (genotypeMap.containsKey(gene)) {
                        genotypeList.add(genotypeMap.get(gene));
                    } else {
                        genotypeList.add('-');
                        missing++;
                    }
                }
                // output a line assuming H if missing; use 1 Mb = 1 cM
                if (missing<=maxMissing) {
                    System.out.print(gene+","+chr+","+cmf.format((double)pos1*1e-6));
                    for (Character genotype : genotypeList) {
                        if ((char)genotype == '-') {
                            System.out.print(","+'H');
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
