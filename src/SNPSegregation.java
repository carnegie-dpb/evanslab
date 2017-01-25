package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;

/**
 * Run through all the SNP files in the given directory and accumulate segregation statistics.
 */
public class SNPSegregation {

    public static void main(String[] args) {

        if (args.length!=11) {
            System.out.println("Usage: SNPSegregation hom|het|both <SNP-dir> <prefix> <suffix> <Amin> <Amax> <Bmin> <Bmax> <minDepth> <minSamples> tsv|wig");
            System.out.println("Example: SNPSegregation hom ~/snps VT avinput.variant_function.snp 1 48 49 96 5 20");
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
            int minDepth = Integer.parseInt(args[8]);       // minimum alt depth to count a sample
            int minSamples = Integer.parseInt(args[9]);     // minimum number of samples with SNP to be included in output
            boolean tsvOutput = args[10].toLowerCase().equals("tsv");
            boolean wigOutput = args[10].toLowerCase().equals("wig");

            if (!tsvOutput && !wigOutput) {
                System.err.println("Error: "+args[10]+" is neither tsv nor wig.");
                System.exit(1);
            }

            // count of samples per groups
            int A = 0;
            int B = 0;

            // holds sample counts
            Map<String,Integer> aSamplesMap = new TreeMap<String,Integer>();
            Map<String,Integer> bSamplesMap = new TreeMap<String,Integer>();

            // holds total reads
            Map<String,Integer> aReadsMap = new TreeMap<String,Integer>();
            Map<String,Integer> bReadsMap = new TreeMap<String,Integer>();

            for (int i=aMin; i<=bMax; i++) {

                boolean isA = (aMin<=i && i<=aMax);
                boolean isB = (bMin<=i && i<=bMax);

                if (isA || isB) {

                    String fileName = prefix+String.valueOf(i)+"."+suffix;
                    File file = new File(dir, fileName);
                    if (file.exists()) {

                        if (isA) A++;
                        if (isB) B++;

                        BufferedReader reader = new BufferedReader(new FileReader(file));
                        String line = reader.readLine(); // header line
                        while((line=reader.readLine())!=null) {

                            SNPRecord rec = new SNPRecord(line);
                            String key = rec.chr+"\t"+String.format("%10d", rec.pos1)+"\t"+rec.ref+"\t"+rec.alt; // keep track of individual alleles

                            // filtering
                            boolean keep = true;
                            if (hetOnly) keep = keep && rec.genotype.equals("het");
                            if (homOnly) keep = keep && rec.genotype.equals("hom");
                            keep = keep && rec.altDepth>=minDepth;

                            if (keep && isA) {
                                    
                                if (aSamplesMap.containsKey(key)) {
                                    // increment sample count
                                    int count = aSamplesMap.get(key);
                                    aSamplesMap.put(key, count+1);
                                    // increment alt reads
                                    int reads = aReadsMap.get(key);
                                    aReadsMap.put(key, reads+rec.altDepth);
                                } else {
                                    // new SNP
                                    aSamplesMap.put(key, 1);
                                    aReadsMap.put(key, rec.altDepth);
                                }
                                
                            } else if (keep && isB) {
                                        
                                if (bSamplesMap.containsKey(key)) {
                                    // increment sample count
                                    int count = bSamplesMap.get(key);
                                    bSamplesMap.put(key, count+1);
                                    // increment alt reads
                                    int reads = bReadsMap.get(key);
                                    bReadsMap.put(key, reads+rec.altDepth);
                                } else {
                                    // new SNP
                                    bSamplesMap.put(key, 1);
                                    bReadsMap.put(key, rec.altDepth);
                                }
                                
                            }
                            
                        }
                        
                    }
                    
                }
                
            }

            // TSVheader line
            if (tsvOutput) System.out.println("contig\tposition\tref\talt\tAcount\tBcount\tAreads\tBreads\tlogOR");

            // merge the A and B keys into a single sorted set, since some are only in one map and some are only in the other
            Set<String> keySet = new TreeSet<String>();
            keySet.addAll(aSamplesMap.keySet());
            keySet.addAll(bSamplesMap.keySet());
            
            // spit out line in wig format when new chromosome
            // NOTE: there will be duplicate lines for SNPs at the same position; wig file must be pruned before converting to BigWig.
            String chrom = "";

            // now run through all the keys, outputting the stats
            for (String key : keySet) {
                String[] pieces = key.split("\t");
                int position = Integer.parseInt(pieces[1].trim());
                if (!chrom.equals(pieces[0])) {
                    chrom = pieces[0];
                    if (wigOutput) System.out.println("variableStep\tchrom="+chrom);
                }
                if (aSamplesMap.containsKey(key) && bSamplesMap.containsKey(key)) {
                    int a = aSamplesMap.get(key);
                    int b = bSamplesMap.get(key);
                    if (a>=minSamples || b>=minSamples) {
                        // fudge case where a==A or b==B by subtracting 1 from both counts
                        if (a==A || b==B) {
                            a = a-1;
                            b = b-1;
                        }
                        double oddsRatio = ((double)b / (double)(B-b)) / ((double)a / (double)(A-a));
                        if (tsvOutput) {
                            System.out.println(key+"\t"+aSamplesMap.get(key)+"\t"+bSamplesMap.get(key)+"\t"+aReadsMap.get(key)+"\t"+bReadsMap.get(key)+"\t"+Math.log(oddsRatio));
                        } else if (wigOutput) {
                            System.out.println(position+"\t"+Math.log(oddsRatio));
                        }
                    }
                } else if (aSamplesMap.containsKey(key)) {
                    int a = aSamplesMap.get(key);
                    int b = 0;
                    if (a>=minSamples) {
                        // fudge case where b=0 by adding 1 to both counts
                        a++;
                        b++;
                        // more fudge
                        if (a>=A) a = A - 1;
                        double oddsRatio = ((double)b / (double)(B-b)) / ((double)a / (double)(A-a));
                        if (tsvOutput) {
                            System.out.println(key+"\t"+aSamplesMap.get(key)+"\t"+0+"\t"+aReadsMap.get(key)+"\t"+0+"\t"+Math.log(oddsRatio));
                        } else if (wigOutput) {
                            System.out.println(position+"\t"+Math.log(oddsRatio));
                        }
                    }
                } else {
                    int a = 0;
                    int b = bSamplesMap.get(key);
                    if (b>=minSamples) {
                        // fudge case where a=0 by adding 1 to both counts
                        a++;
                        b++;
                        // more fudge
                        if (b>=B) b = B - 1;
                        double oddsRatio = ((double)b / (double)(B-b)) / ((double)a / (double)(A-a));
                        if (tsvOutput) {
                            System.out.println(key+"\t"+0+"\t"+bSamplesMap.get(key)+"\t"+0+"\t"+bReadsMap.get(key)+"\t"+Math.log(oddsRatio));
                        } else if (wigOutput) {
                            System.out.println(position+"\t"+Math.log(oddsRatio));
                        }
                    }
                }
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
