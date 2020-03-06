package edu.carnegiescience.dpb.evanslab;

import java.util.Arrays;
import java.util.List;

/**
 * Encapsulates a SNP record from a tab-separated Novogene .snp file:
 *
 * 0       1       2       3       4       5       6       7          8       9       10              11
 * annoPos info    contig  pos1    pos2    ref     alt     genotype   quality MQ      ref-depth       alt-depth
 */
public class SNPRecord {

    public String annoPos;
    public String info;
    public String contig;
    public int pos1;
    public int pos2;
    public char ref;
    public char alt;
    public String genotype;
    public double quality;
    public int mq;
    public int refDepth;
    public int altDepth;

    /**
     * Instantiate from a SNP line by parsing the elements.
     */
    public SNPRecord(String line) {
        String[] fields = line.split("\t");
        annoPos = fields[0];
        info = fields[1];
        contig = fields[2];
        pos1 = Integer.parseInt(fields[3]);
        pos2 = Integer.parseInt(fields[4]);
        ref = fields[5].charAt(0);
        alt = fields[6].charAt(0);
        genotype = fields[7];
        try { quality = Double.parseDouble(fields[8]); } catch (Exception e) { quality = 0.0; }
        try { mq = Integer.parseInt(fields[9]); } catch (Exception e) { mq = 0; }
        try { refDepth = Integer.parseInt(fields[10]); } catch (Exception e) { refDepth = 0; }
        try { altDepth = Integer.parseInt(fields[11]); } catch (Exception e) { altDepth = 0; }
    }

    /**
     * Instantiate from individual values
     */
    public SNPRecord(String annoPos, String info, String contig, int pos1, int pos2, char ref, char alt, String genotype, double quality, int mq, int refDepth, int altDepth) {
        this.annoPos = annoPos;
        this.info = info;
        this.contig = contig;
        this.pos1 = pos1;
        this.pos2 = pos2;
        this.ref = ref;
        this.alt = alt;
        this.genotype = genotype;
        this.quality = quality;
        this.mq = mq;
        this.refDepth = refDepth;
        this.altDepth = altDepth;
    }

    /**
     * Return true if this is a homozygous SNP.
     */
    public boolean isHomozygous() {
        return genotype.equals("hom");
    }

    /**
     * Return true if this is a heterozygous SNP.
     */
    public boolean isHeterozygous() {
        return genotype.equals("het");
    }

    /**
     * Return true if this is actually a Ref record.
     */
    public boolean isReference() {
        return genotype.equals("ref");
    }

    /**
     * Return true if this record is flagged as exonic
     */
    public boolean isExonic() {
        return annoPos.equals("exonic");
    }

    /**
     * Output the standardized string representation as used in a snp file.
     */
    public String toString() {
        return annoPos+"\t"+info+"\t"+contig+"\t"+pos1+"\t"+pos2+"\t"+ref+"\t"+alt+"\t"+genotype+"\t"+quality+"\t"+mq+"\t"+refDepth+"\t"+altDepth;
    }

    /** 
     * Return true if contig is a number. Useful for not including scaffolds.
     */
    public boolean hasNumericContig() {
        try {
            int i = Integer.parseInt(contig);
            return true;
        } catch (Exception ex) {
            return false;
        }
    }
}
