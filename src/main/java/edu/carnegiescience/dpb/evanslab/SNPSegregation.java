package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.text.DecimalFormat;

import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;

import org.mskcc.cbio.portal.stats.FisherExact;

/**
 * Run through all the SNP files in the given directory and accumulate segregation statistics subject to a command-line filter.
 */
public class SNPSegregation {

    static DecimalFormat posFormat = new DecimalFormat("0000000000");

    public static void main(String[] args) throws FileNotFoundException, IOException {

        if (args.length!=8) {
            System.out.println("Usage: SNPSegregation <SNP-dir> <prefix> <suffix> <Amin> <Amax> <Bmin> <Bmax> <minDepth>");
            System.out.println("Example: SNPSegregation ~/snps VT avinput.variant_function.snp 1 48 49 96 5");
            System.exit(1);
        }
        
        String snpDirName = args[0];
        String prefix = args[1];
        String suffix = args[2];
        int aMin = Integer.parseInt(args[3]);
        int aMax = Integer.parseInt(args[4]);
        int bMin = Integer.parseInt(args[5]);
        int bMax = Integer.parseInt(args[6]);
        int minDepth = Integer.parseInt(args[7]);

        // A group sample count
        int A = 0;
        // total reads, Het sample counts, Hom sample counts
        Map<String,Integer> aRefReadsMap = new TreeMap<String,Integer>();
        Map<String,Integer> aAltReadsMap = new TreeMap<String,Integer>();
        Map<String,Integer> aHetMap = new TreeMap<String,Integer>();
        Map<String,Integer> aHomMap = new TreeMap<String,Integer>();
        for (int i=aMin; i<=aMax; i++) {
            String fileName = snpDirName+"/"+prefix+String.valueOf(i)+"/"+prefix+String.valueOf(i)+"."+suffix;
            File file = new File(fileName);
            if (file.exists()) {
                // increment sample counter
                A++;
                BufferedReader reader = new BufferedReader(new FileReader(file));
                String line = reader.readLine(); // header line
                while((line=reader.readLine())!=null) {
                    SNPRecord rec = new SNPRecord(line);
                    if (!rec.hasNumericContig()) {
                        continue;
                    }
                    String key = formKey(rec);
                    // store alt and ref reads no matter what
                    if (aRefReadsMap.containsKey(key)) {
                        // increment reads
                        int refReads = aRefReadsMap.get(key) + rec.refDepth;
                        aRefReadsMap.put(key, refReads);
                        int altReads = aAltReadsMap.get(key) + rec.altDepth;
                        aAltReadsMap.put(key, altReads);
                    } else {
                        // new overall SNP
                        aRefReadsMap.put(key, rec.refDepth);
                        aAltReadsMap.put(key, rec.altDepth);
                    }
                    // filtering
                    if ((rec.refDepth+rec.altDepth)>=minDepth) {
                        if (rec.isHomozygous()) {
                            if (aHomMap.containsKey(key)) {
                                // increment hom count here
                                int count = aHomMap.get(key) + 1;
                                aHomMap.put(key, count);
                            } else {
                                // new Hom SNP
                                aHomMap.put(key, 1);
                            }
                        } else if (rec.isHeterozygous()) {
                            if (aHetMap.containsKey(key)) {
                                // increment het count here
                                int count = aHetMap.get(key) + 1;
                                aHetMap.put(key, count);
                            } else {
                                // new Het SNP
                                aHetMap.put(key, 1);
                            }
                        } else {
                            System.err.println("No het/hom call on this line:");
                            System.err.println(line);
                            System.exit(1);
                        }
                    }
                }
            }
        }

        // B group sample count
        int B = 0;
        // total reads, Het sample counts, Hom sample counts
        Map<String,Integer> bRefReadsMap = new TreeMap<String,Integer>();
        Map<String,Integer> bAltReadsMap = new TreeMap<String,Integer>();
        Map<String,Integer> bHetMap = new TreeMap<String,Integer>();
        Map<String,Integer> bHomMap = new TreeMap<String,Integer>();
        for (int i=bMin; i<=bMax; i++) {
            String fileName = snpDirName+"/"+prefix+String.valueOf(i)+"/"+prefix+String.valueOf(i)+"."+suffix;
            File file = new File(fileName);
            if (file.exists()) {
                // increment sample counter
                B++;
                BufferedReader reader = new BufferedReader(new FileReader(file));
                String line = reader.readLine(); // header line
                while((line=reader.readLine())!=null) {
                    SNPRecord rec = new SNPRecord(line);
                    if (!rec.hasNumericContig()) {
                        continue;
                    }
                    String key = formKey(rec);
                    // count alt and ref reads no matter what
                    if (bRefReadsMap.containsKey(key)) {
                        // increment reads
                        int refReads = bRefReadsMap.get(key) + rec.refDepth;
                        bRefReadsMap.put(key, refReads);
                        int altReads = bAltReadsMap.get(key) + rec.altDepth;
                        bAltReadsMap.put(key, altReads);
                    } else {
                        // new overall SNP
                        bRefReadsMap.put(key, rec.refDepth);
                        bAltReadsMap.put(key, rec.altDepth);
                    }
                    // filtering
                    if ((rec.refDepth+rec.altDepth)>=minDepth) {
                        if (rec.isHomozygous()) {
                            if (bHomMap.containsKey(key)) {
                                // increment hom count here
                                int count = bHomMap.get(key) + 1;
                                bHomMap.put(key, count);
                            } else {
                                // new Hom SNP
                                bHomMap.put(key, 1);
                            }
                        } else if (rec.isHeterozygous()) {
                            if (bHetMap.containsKey(key)) {
                                // increment het count here
                                int count = bHetMap.get(key) + 1;
                                bHetMap.put(key, count);
                            } else {
                                // new Het SNP
                                bHetMap.put(key, 1);
                            }
                        } else {
                            System.err.println("No het/hom call on this line:");
                            System.err.println(line);
                            System.exit(1);
                        }
                    }
                }
            }
        }

        // header line
        System.out.println("contig\tpos\tAref\tAalt\tAhet\tAhom\tBref\tBalt\tBhet\tBhom\tOR\tp");
            
        // merge the A and B keys into a single sorted set, since some are only in one map and some are only in the other
        Set<String> keySet = new TreeSet<String>();
        keySet.addAll(aRefReadsMap.keySet());
        keySet.addAll(bRefReadsMap.keySet());

        // Fisher's exact test
        FisherExact fisherExact = new FisherExact(A+B);
        
        // now run through all the keys, outputting the stats
        for (String key : keySet) {
            String[] pieces = key.split("_");
            String contig = getContig(key);
            int position = getPosition(key);
            int aRefReads = 0;
            int aAltReads = 0;
            int aHet = 0;
            int aHom = 0;
            int bRefReads = 0;
            int bAltReads = 0;
            int bHet = 0;
            int bHom = 0;
            if (aRefReadsMap.containsKey(key)) aRefReads = aRefReadsMap.get(key);
            if (aAltReadsMap.containsKey(key)) aAltReads = aAltReadsMap.get(key);
            if (aHetMap.containsKey(key)) aHet = aHetMap.get(key);
            if (aHomMap.containsKey(key)) aHom = aHomMap.get(key);
            if (bRefReadsMap.containsKey(key)) bRefReads = bRefReadsMap.get(key);
            if (bAltReadsMap.containsKey(key)) bAltReads = bAltReadsMap.get(key);
            if (bHetMap.containsKey(key)) bHet = bHetMap.get(key);
            if (bHomMap.containsKey(key)) bHom = bHomMap.get(key);
            // output
            if (aHet>0 & bHet>0) {
                int aRef = A - aHet;
                int bRef = B - bHet;
                double or = ( (double)bHet/(double)bRef ) / ( (double)aHet/(double)aRef );
                double p = fisherExact.getP(aHet, aRef, bHet, bRef);
                System.out.println(contig+"\t"+position+"\t"+aRefReads+"\t"+aAltReads+"\t"+aHet+"\t"+aHom+"\t"+bRefReads+"\t"+bAltReads+"\t"+bHet+"\t"+bHom+"\t"+or+"\t"+p);
            }
        }
    }

    static String formKey(SNPRecord rec) {
        return rec.contig+":"+posFormat.format(rec.pos1);
    }
    
    static String getContig(String key) {
        String[] pieces = key.split(":");
        return pieces[0];
    }
    
    static int getPosition(String key) {
        String[] pieces = key.split(":");
        return Integer.parseInt(pieces[1].trim());
    }
}
