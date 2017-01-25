package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Map;
import java.util.TreeMap;
import java.util.Set;
import java.util.TreeSet;

/**
 * Prune values at the same locus from an input WIG file. Saved value is the one that is largest in the positive direction.
 */
public class WigPruner {

    public static void main(String[] args) {

        if (args.length!=1) {
            System.out.println("Usage: WigPruner <WIG file>");
            System.exit(1);
        }
        
        try {
            
            File file = new File(args[0]);
            if (!file.exists()) {
                System.err.println("ERROR: file "+args[0]+" not found.");
                System.exit(1);
            }

            BufferedReader reader = new BufferedReader(new FileReader(file));
            int oldPos = 0;
            double oldVal = 0.0;
            String line;
            while ((line=reader.readLine())!=null) {
                if (line.startsWith("variableStep")) {
                    if (oldPos!=0) System.out.println(oldPos+"\t"+oldVal);
                    String[] parts = line.split("\t");
                    if (parts[1].startsWith("chrom=B")) break;
                    System.out.println(line);
                    oldPos = 0;
                } else {
                    String[] parts = line.split("\t");
                    int pos = Integer.parseInt(parts[0]);
                    double val = Double.parseDouble(parts[1]);
                    if (pos==oldPos) {
                        if (val>oldVal) oldVal = val;
                    } else {
                        if (oldPos!=0) System.out.println(oldPos+"\t"+oldVal);
                        oldPos = pos;
                        oldVal = val;
                    }
                }
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }

}
