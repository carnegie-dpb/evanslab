package edu.carnegiescience.dpb.evanslab;

import java.util.Arrays;
import java.util.List;

/**
 * Encapsulates an exon record from the B73v4 GFF3 file.
 */
public class GFFExonRecord {

    // main fields
    public String chrom;
    public String source;
    public String type;
    public int start;
    public int end;
    public char strand;
    public String attributes;

    // exon attributes
    public String name;
    public String parent;
    public int constituitive;
    public int ensemblEndPhase;
    public String exonId;
    public int rank;
    public int version;

    /**
     * Instantiate from a GFF line by parsing the elements.
     */
    public GFFExonRecord(String line) {

        // get tab-separated fields
        String[] fields = line.split("\t");
            
        chrom = fields[0];
        source = fields[1];
        type = fields[2];
        start = Integer.parseInt(fields[3]);
        end = Integer.parseInt(fields[4]);
        strand = fields[6].charAt(0);
        attributes = fields[8];
    
        // get semicolon-separated attributes for exon
        if (fields[2].equals("exon")) {
            List<String> chunks = Arrays.asList(attributes.split(";"));
            for (String chunk : chunks) {
                String[] sides = chunk.split("=");
                if (sides[0].equals("Name")) {
                    name = sides[1];
                } else if (sides[0].equals("Parent")) {
                    parent = sides[1];
                } else if (sides[0].equals("constituitive")) {
                    constituitive = Integer.parseInt(sides[1]);
                } else if (sides[0].equals("ensembl_end_phase")) {
                    ensemblEndPhase = Integer.parseInt(sides[1]);
                } else if (sides[0].equals("exon_id")) {
                    exonId = sides[1];
                } else if (sides[0].equals("rank")) {
                    rank = Integer.parseInt(sides[1]);
                } else if (sides[0].equals("version")) {
                    version = Integer.parseInt(sides[1]);
                }
            }
        }

    }

    /**
     * Return true if this is an exon record.
     */
    public boolean isExon() {
        return exonId!=null;
    }

    /**
     * Standardized string representation
     */
    public String toString() {

        return chrom+":"+start+"-"+end+"\t"+strand+"\t"+name;

    }
    
}
