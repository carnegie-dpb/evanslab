package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Vector;
import java.util.TreeMap;

/**
 * Remove multiple values at the same locus from an input WIG file.
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

            // store positions in a Vector per chromosome and only spit out the ones that occur once.
            Vector<Integer> positions = null;
            TreeMap<Integer,Double> values = null;

            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line;
            while ((line=reader.readLine())!=null) {
                String[] parts = line.split("\t");
                if (parts[0].equals("variableStep")) {
                    if (positions!=null) {
                        // output previous chromosome data
                        for (Integer position : values.keySet()) {
                            int i1 = positions.indexOf(position);
                            int i2 = positions.lastIndexOf(position);
                            if (i1==i2) {
                                // unique position
                                System.out.println(position+"\t"+values.get(position));
                            }
                        }
                    }
                    // initialize new chromosome
                    System.out.println(line);
                    positions = new Vector<Integer>();
                    values = new TreeMap<Integer,Double>();
                } else {
                    int pos = Integer.parseInt(parts[0]);
                    double val = Double.parseDouble(parts[1]);
                    positions.add(pos);
                    values.put(pos, val);
                }
            }

            // output last chromosome data
            for (Integer position : values.keySet()) {
                int i1 = positions.indexOf(position);
                int i2 = positions.lastIndexOf(position);
                if (i1==i2) {
                    // unique position
                    System.out.println(position+"\t"+values.get(position));
                }
            }


        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }

}
