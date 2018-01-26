package edu.carnegiescience.dpb.evanslab;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

/**
 * Encapsulates a Reference record from a custom tab-separated reference call file.
 *
 * contig  position S1    S2    S3    ...    Sn
 * 1       8362870  NC    Het   Hom   ...    R:4:44.2
 *
 */
public class RefRecord {

    public String contig;
    public int position;
    public Map<String,String> calls;

    /**
     * Create a default record, identified by position=0.
     */
    public RefRecord(List<String> samples) {
        contig = null;
        position = 0;
        calls = new HashMap<String,String>();
        for (String sample : samples) {
            calls.put(sample, null);
        }
    }
    
    /**
     * Instantiate from a SNP line by parsing the elements and matching with samples.
     */
    public RefRecord(String line, List<String> samples) {
        String[] fields = line.split("\t");
        contig = fields[0];
        position = Integer.parseInt(fields[1]);
        calls = new HashMap<String,String>();
        for (int i=0; i<fields.length-2; i++) {
            calls.put(samples.get(i), fields[i+2]);
        }
    }

    /**
     * Load the samples from the header line.
     */
    public static List<String> getSamples(String header) {
        String[] fields = header.split("\t");
        List<String> samples = new ArrayList<String>();
        for (int i=0; i<fields.length-2; i++) {
            samples.add(fields[i+2]);
        }
        return samples;
    }

    /**
     * String summary output.
     */
    public String toString() {
        String out = contig+"\t"+position;
        for (String sample : calls.keySet()) {
            out += "\t"+sample+":"+calls.get(sample);
        }
        return out;
    }

    /**
     * Is the call for the given sample a Ref call?
     */
    public boolean isRef(String sample) {
        String[] fields = calls.get(sample).split(":");
        return fields[0].equals("R");
    }
        
    /**
     * Get the number of Ref reads for a sample; 0 if it's not a Ref.
     */
    public int getRefReads(String sample) {
        String[] fields = calls.get(sample).split(":");
        if (fields.length==3) {
            return Integer.parseInt(fields[1]);
        } else {
            return 0;
        }
    }

    /**
     * Get the Ref read quality for a sample; 0.0 if it's not a Ref.
     */
    public double getRefQuality(String sample) {
        String[] fields = calls.get(sample).split(":");
        if (fields.length==3) {
            return Double.parseDouble(fields[2]);
        } else {
            return 0.0;
        }
    }
    

    
}
