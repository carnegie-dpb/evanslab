package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.TreeMap;

/**
 * Simply counts similar reads and spits them out with their counts.
 *
 * @author Sam Hokin
 */
public class ReadCounter {

    public static void main(String[] args) throws Exception {
        if (args.length!=2) {
            System.out.println("Usage: ReadCounter <reads file> min-count");
            System.exit(0);
        }
        
        String filename = args[0];
        int minCount = Integer.parseInt(args[1]);
        
        TreeMap<String,Integer> readCounts = new TreeMap<>();
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = null;
        while ((line=reader.readLine())!=null) {
            if (readCounts.containsKey(line)) {
                int count = readCounts.get(line) + 1;
                readCounts.put(line, count);
            } else {
                readCounts.put(line, 1);
            }
        }
        reader.close();

        // output reads above a minimum count as a FASTA
        int num = 0;
        for (String seq : readCounts.keySet()) {
            
            if (readCounts.get(seq)>=minCount) {
                num++;
                System.out.println(">read"+num+"-"+readCounts.get(seq));
                System.out.println(seq);
            }
        }
    }
}
