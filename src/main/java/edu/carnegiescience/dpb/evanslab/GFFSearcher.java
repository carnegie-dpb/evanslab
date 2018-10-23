package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.genome.parsers.gff.Location;

/**
 * Searches a GFF file for features that span given locations.
 */
public class GFFSearcher {

    FeatureList allFeatures;

    /**
     * Construct from a GFF file.
     */
    public GFFSearcher(String gffFilename) throws IOException {
        allFeatures = GFF3Reader.read(gffFilename);
    }

    public static void main(String[] args) throws Exception {

        if (args.length!=2) {
            System.out.println("Usage GFFSearcher <gff-file> <location-file>");
            System.exit(0);
        }
        String gffFilename = args[0];
        String locFilename = args[1];

        GFFSearcher searcher = new GFFSearcher(gffFilename);

        BufferedReader reader = new BufferedReader(new FileReader(locFilename));
        String line;
        while ((line=reader.readLine())!=null) {
            String[] parts = line.split("\t");
            String seqname = parts[0];
            int pos = Integer.parseInt(parts[1]);
            Location location = new Location(pos,pos);
            FeatureList overlapping = searcher.search(seqname, location);
            if (overlapping.size()>0) {
                // ONLY SPIT OUT THE FIRST MATCH!!!
                FeatureI feature = overlapping.get(0);
                System.out.println(seqname+"\t"+pos+"\t"+feature.getAttribute("ID")+"\t"+feature.getAttribute("Name"));
            } else {
                System.out.println(seqname+"\t"+pos+"\t"+"\t");
            }                
        }
        
    }

    /**
     * Search for a given location
     */
    public FeatureList search(String seqname, Location location) throws Exception {
        return allFeatures.selectOverlapping(seqname, location, true);
    }

}
