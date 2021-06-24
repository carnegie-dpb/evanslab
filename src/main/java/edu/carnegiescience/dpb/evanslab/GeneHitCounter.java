package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.HashMap;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.genome.parsers.gff.Location;

/**
 * Counts hits against genes from a BreakpointAnalyzer output file.
 * Feature positions are given by a GFF file.
 *
 * @author Sam Hokin
 */
public class GeneHitCounter {

    public static void main(String[] args) throws Exception {
        if (args.length!=2) {
            System.out.println("Usage: FeatureHitCounter <BreakpointAnalyzer output> <GFF file>");
            System.exit(0);
        }
        
        String bpaFilename = args[0];
        String gffFilename = args[1];

        // grab genes from GFF
        final FeatureList genes = GFF3Reader.read(gffFilename).selectByType("gene");
        System.err.println("Read "+genes.size()+" genes from "+gffFilename+".");

        // store per-gene counts in a HashMap
        HashMap<FeatureI,Integer> geneCounts = new HashMap<>();

        // 0=read                                  1=chr 2=hitstart 3=hitend   4=hitstrand  5=otherstart  6=otherend  7=otherstrand
        // A00197:89:HHHH2DMXX:1:1101:21296:28792  2     206608590  206608734  -            128247770     128247918   +
        BufferedReader bpaReader = new BufferedReader(new FileReader(bpaFilename));
        String line = null;
        while ((line=bpaReader.readLine())!=null) {
            String[] parts = line.split("\t");
            String read = parts[0];
            String chr = parts[1];
            // read that hit an exon
            int start = Integer.parseInt(parts[2]);
            int end = Integer.parseInt(parts[3]);
            Location location = null;
            if (parts[4].equals("+")) {
                location = new Location(start, end);
            } else {
                location = new Location(-end, -start);
            }
            // strand-specific hit-gene overlap 
            FeatureList overlaps = genes.selectOverlapping(chr, location, false);
            for (FeatureI gene : overlaps) {
                if (geneCounts.containsKey(gene)) {
                    int count = geneCounts.get(gene) + 1;
                    geneCounts.put(gene, count);
                } else {
                    geneCounts.put(gene, 1);
                }
            }
        }
        bpaReader.close();

        // output
        for (FeatureI gene : geneCounts.keySet()) {
            System.out.println("Exon read\t"+geneCounts.get(gene)+"\t"+gene.getAttribute("ID")+"\t"+gene.seqname()+":"+gene.location().bioStart()+"-"+gene.location().bioEnd());
        }
    }
    
}
