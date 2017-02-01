package edu.carnegiescience.dpb.evanslab;

import java.util.Arrays;
import java.util.List;

/**
 * Encapsulates a SNP record from a tab-separated Novogene .snp file:
 *
 * annoPos  gene                            chr     pos1    pos2    ref     alt     genotype quality MQ      ref-depth  alt-depth
 * upstream transcript:Zm00001d027578_T001  1       8362870 8362870 C       A       het      26      25      1          5
 *
 */
public class SNPRecord {

    public String annoPos;
    public String gene;
    public String chr;
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

        // get tab-separated fields
        String[] fields = line.split("\t");
            
        annoPos = fields[0];
        gene = fields[1];
        chr = fields[2];
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
    public SNPRecord(String annoPos, String gene, String chr, int pos1, int pos2, char ref, char alt, String genotype, double quality, int mq, int refDepth, int altDepth) {
        this.annoPos = annoPos;
        this.gene = gene;
        this.chr = chr;
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
        return annoPos+"\t"+gene+"\t"+chr+"\t"+pos1+"\t"+pos2+"\t"+ref+"\t"+alt+"\t"+genotype+"\t"+quality+"\t"+mq+"\t"+refDepth+"\t"+altDepth;
    }
    
}
