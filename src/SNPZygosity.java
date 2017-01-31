package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;

/**
 * Run through all the REF/SNP files in the given directory and accumulate the counts of reference, homozygous and heterozygous samples per locus per group.
 * The SNP files are in the directory given by SNP-dir.
 *
 * prefix = file name prefix (e.g. VT).
 * suffix = file name suffix (e.g. avinput.variant_function.snp).
 * Amin, Amax = min, max index of the A group of samples. (They are assumed to be contiguous.)
 * Bmin, Bmax = min, max index of the B group. (They are assumed to be contiguous.)
 * minDepth = minimum alt read depth for a hom SNP; minimum alt+ref depth for a het SNP; minimum ref depth for a Ref call to be counted.
 * minSamples = minimum number of samples for a SNP to be output. 
 *
 * The tab-separated REF/SNP file format is as follows:
 *
 * annoPos  gene                            chr     pos1    pos2    ref     alt     genotype quality MQ      ref-depth  alt-depth
 * upstream transcript:Zm00001d027578_T001  1       8362870 8362870 C       A       het      26.8    25      1          5
 * .        .                               1       8362912 8362912 .       .       ref      42.2    .       23         .
 *
 */
public class SNPZygosity {

    static double MIN_QUALITY = 20.0;
    static DecimalFormat df = new DecimalFormat("0.0000");

    public static void main(String[] args) {

        if (args.length!=9) {
            System.out.println("Usage: SNPZygosity <SNP-dir> <prefix> <suffix> <Amin> <Amax> <Bmin> <Bmax> <minDepth> <minSamples>");
            System.out.println("Example: SNPZygosity ~/snps VT avinput.variant_function.snp 1 48 49 96 5 20");
            System.exit(1);
        }
        
        try {
            
            File dir = new File(args[0]);
            if (!dir.isDirectory()) {
                System.err.println("Error: "+args[0]+" is not a directory.");
                System.exit(1);
            }

            String prefix = args[1];
            String suffix = args[2];
            int aMin = Integer.parseInt(args[3]);
            int aMax = Integer.parseInt(args[4]);
            int bMin = Integer.parseInt(args[5]);
            int bMax = Integer.parseInt(args[6]);
            int minDepth = Integer.parseInt(args[7]);       // minimum alt depth to count a sample
            int minSamples = Integer.parseInt(args[8]);     // minimum number of samples with SNP to be included in output

            // number of samples per group
            int A = 0;
            int B = 0;

            // holds count of homozygous samples per SNP
            Map<String,Integer> aHomMap = new TreeMap<String,Integer>();
            Map<String,Integer> bHomMap = new TreeMap<String,Integer>();

            // holds count of heterozygous samples per SNP
            Map<String,Integer> aHetMap = new TreeMap<String,Integer>();
            Map<String,Integer> bHetMap = new TreeMap<String,Integer>();

            // holds count of reference samples per SNP
            Map<String,Integer> aRefMap = new TreeMap<String,Integer>();
            Map<String,Integer> bRefMap = new TreeMap<String,Integer>();

            // holds total read counts
            Map<String,Integer> aReadsMap = new TreeMap<String,Integer>();
            Map<String,Integer> bReadsMap = new TreeMap<String,Integer>();

            // A group
            for (int i=aMin; i<=aMax; i++) {
                String fileName = prefix+String.valueOf(i)+"."+suffix;
                File file = new File(dir, fileName);
                if (file.exists()) {
                    // increment sample count
                    A++;
                    // spin through the records
                    BufferedReader reader = new BufferedReader(new FileReader(file));
                    // String line = reader.readLine(); // header line
                    String line = null; // forgot to put header lines on snp-rec files
                    while((line=reader.readLine())!=null) {
                        SNPRecord rec = new SNPRecord(line);
                        String key = rec.chr+"\t"+String.format("%10d", rec.pos1); // NOT allele-specific since we don't have them for Ref calls
                        if (rec.quality>=MIN_QUALITY) {
                            if (rec.isHomozygous() && rec.altDepth>=minDepth) {
                                if (aHomMap.containsKey(key)) {
                                    // increment homozygous count
                                    int count = aHomMap.get(key) + 1;
                                    aHomMap.put(key, count);
                                    // increment reads count
                                    int reads = aReadsMap.get(key) + rec.altDepth;
                                    aReadsMap.put(key, reads);
                                } else {
                                    // new homozygous SNP
                                    aHomMap.put(key, 1);
                                    aReadsMap.put(key, rec.altDepth);
                                }
                            } else if (rec.isHeterozygous() && (rec.refDepth+rec.altDepth)>=minDepth) {
                                if (aHetMap.containsKey(key)) {
                                    // increment heterozygous count
                                    int count = aHetMap.get(key) + 1;
                                    aHetMap.put(key, count);
                                    // increment reads count
                                    int reads = aReadsMap.get(key) + rec.altDepth;
                                    aReadsMap.put(key, reads);
                                } else {
                                    // new heterozygous SNP
                                    aHetMap.put(key, 1);
                                    aReadsMap.put(key, rec.altDepth);
                                }
                            } else if (rec.isReference() && rec.refDepth>=minDepth) {
                                if (aRefMap.containsKey(key)) {
                                    // increment reference count
                                    int count = aRefMap.get(key) + 1;
                                    aRefMap.put(key, count);
                                    // increment reads count
                                    int reads = aReadsMap.get(key) + rec.refDepth;
                                    aReadsMap.put(key, reads);
                                } else {
                                    // new reference call
                                    aRefMap.put(key, 1);
                                    aReadsMap.put(key, rec.refDepth);
                                }
                            }
                        }
                    }
                }
            }

            // B group
            for (int i=bMin; i<=bMax; i++) {
                String fileName = prefix+String.valueOf(i)+"."+suffix;
                File file = new File(dir, fileName);
                if (file.exists()) {
                    // increment sample count
                    B++;
                    // spin through the records
                    BufferedReader reader = new BufferedReader(new FileReader(file));
                    String line = reader.readLine(); // header line
                    while((line=reader.readLine())!=null) {
                        SNPRecord rec = new SNPRecord(line);
                        String key = rec.chr+"\t"+String.format("%10d", rec.pos1); // NOT allele-specific since we don't have them for Ref calls
                        if (rec.quality>=MIN_QUALITY) {
                            if (rec.isHomozygous() && rec.altDepth>=minDepth) {
                                if (bHomMap.containsKey(key)) {
                                    // increment homozygous count
                                    int count = bHomMap.get(key) + 1;
                                    bHomMap.put(key, count);
                                    // increment reads count
                                    int reads = bReadsMap.get(key) + rec.altDepth;
                                    bReadsMap.put(key, reads);
                                } else {
                                    // new homozygous SNP
                                    bHomMap.put(key, 1);
                                    bReadsMap.put(key, rec.altDepth);
                                }
                            } else if (rec.isHeterozygous() && (rec.refDepth+rec.altDepth)>minDepth) {
                                if (bHetMap.containsKey(key)) {
                                    // increment heterozygous count
                                    int count = bHetMap.get(key) + 1;
                                    bHetMap.put(key, count);
                                    // increment reads count
                                    int reads = bReadsMap.get(key) + rec.altDepth;
                                    bReadsMap.put(key, reads);
                                } else {
                                    // new heterozygous SNP
                                    bHetMap.put(key, 1);
                                    bReadsMap.put(key, rec.altDepth);
                                }
                            } else if (rec.isReference() && rec.refDepth>=minDepth) {
                                if (bRefMap.containsKey(key)) {
                                    // increment reference count
                                    int count = bRefMap.get(key) + 1;
                                    bRefMap.put(key, count);
                                    // increment reads count
                                    int reads = bReadsMap.get(key) + rec.refDepth;
                                    bReadsMap.put(key, reads);
                                } else {
                                    // new reference call
                                    bRefMap.put(key, 1);
                                    bReadsMap.put(key, rec.refDepth);
                                }

                            }
                        }
                    }
                }
            }

            // output header line
            System.out.println("contig\tposition\tAreads\tAref\tAhom\tAhet\tlnAOR\tBreads\tBref\tBhom\tBhet\tlnBOR");

            // merge the A and B keys into a single sorted set, since some are only in one map and some are only in the other
            Set<String> keySet = new TreeSet<String>();
            keySet.addAll(aReadsMap.keySet());
            keySet.addAll(bReadsMap.keySet());
            
            // now run through all the keys, outputting the stats
            for (String key : keySet) {
                String[] pieces = key.split("\t");
                int position = Integer.parseInt(pieces[1].trim());
                int aReads = 0;
                int bReads = 0;
                int aRef = 0;
                int bRef = 0;
                int aHom = 0;
                int aHet = 0;
                int bHom = 0;
                int bHet = 0;
                if (aReadsMap.get(key)!=null) aReads = aReadsMap.get(key);
                if (bReadsMap.get(key)!=null) bReads = bReadsMap.get(key);
                if (aHomMap.get(key)!=null) aHom = aHomMap.get(key);
                if (bHomMap.get(key)!=null) bHom = bHomMap.get(key);
                if (aHetMap.get(key)!=null) aHet = aHetMap.get(key);
                if (bHetMap.get(key)!=null) bHet = bHetMap.get(key);
                if (aRefMap.get(key)!=null) aRef = aRefMap.get(key);
                if (bRefMap.get(key)!=null) bRef = bRefMap.get(key);
                // have to have at least minSamples from one of the Hom or Het counts (we don't care about mostly Ref cases)
                if (aHom>=minSamples || aHet>=minSamples || bHom>=minSamples || bHet>=minSamples || aRef>=minSamples || bRef>=minSamples) {
                    // hom/het odds ratio for each group
                    double lnAOR = 0.0;
                    if (aHom>0 && aHet>0 && aHom<A && aHet<A) {
                        // unfudged
                        lnAOR = Math.log((double)aHom/(double)aHet * (double)(A-aHet)/(double)(A-aHom));
                    } else {
                        // fudge: add 1 to both sides
                        lnAOR = Math.log((double)(aHom+1)/(double)(aHet+1) * (double)(A+1-aHet)/(double)(A+1-aHom));
                    }
                    double lnBOR = 0.0;
                    if (bHom>0 && bHet>0 && bHom<B && bHet<B) {
                        // unfudged
                        lnBOR = Math.log((double)bHom/(double)bHet * (double)(B-bHet)/(double)(B-bHom));
                    } else {
                        // fudge: add 1 to both sides
                        lnBOR = Math.log((double)(bHom+1)/(double)(bHet+1) * (double)(B+1-bHet)/(double)(B+1-bHom));
                    }
                    // output
                    System.out.println(key+"\t"+aReads+"\t"+aRef+"\t"+aHom+"\t"+aHet+"\t"+df.format(lnAOR)+"\t"+bReads+"\t"+bRef+"\t"+bHom+"\t"+bHet+"\t"+df.format(lnBOR));
                }
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
