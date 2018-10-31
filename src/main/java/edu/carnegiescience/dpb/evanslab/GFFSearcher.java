package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.FileReader;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;

/**
 * Searches a GFF file for features that span given locations.
 *
 * @author Sam Hokin
 */
public class GFFSearcher {

    /**
     * Find features that overlap locations in a provided location file.
     */
    public static void main(String[] args) throws Exception {

        if (args.length!=2) {
            System.out.println("Usage GFFSearcher <gff-file> <location-file>");
            System.exit(0);
        }
        String gffFilename = args[0];
        String locFilename = args[1];

        GFFLoader loader = new GFFLoader(gffFilename);

        BufferedReader reader = new BufferedReader(new FileReader(locFilename));
        String line;
        while ((line=reader.readLine())!=null) {
            String[] parts = line.split("\t");
            String seqname = parts[0];
            int pos = Integer.parseInt(parts[1]);
            Location location = new Location(pos,pos);
            FeatureList overlapping = loader.search(seqname, location);
            if (overlapping.size()>0) {
                // ONLY SPIT OUT THE FIRST MATCH!!!
                FeatureI feature = overlapping.get(0);
                System.out.println(seqname+"\t"+pos+"\t"+feature.getAttribute("ID")+"\t"+feature.getAttribute("Name"));
            } else {
                System.out.println(seqname+"\t"+pos+"\t"+"\t");
            }                
        }
        
    }

}
