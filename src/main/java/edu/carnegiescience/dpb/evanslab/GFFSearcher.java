package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.BiConsumer;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;

/**
 * Searches a GFF file for features that span given locations.
 * Uses ConcurrentHashMap to run the gene search in parallel.
 *
 * Location file from a Fisher run has columns:
 * contig	start	REF	ALT	a	b	c	d	size	p	mlog10p	signif
 * 0            1       2       3       4       5       6       7       8       9       10      11
 *
 * Outputs wig file format.
 *
 * @author Sam Hokin
 */
public class GFFSearcher implements BiConsumer<String,String> {

    String gffFilename, locFilename;
    GFFLoader loader;
    ConcurrentHashMap<String,String> locValues;
    ConcurrentHashMap<String,String> geneLocValues;
    
    GFFSearcher(String gffFilename, String locFilename) throws IOException {
        this.gffFilename = gffFilename;
        this.locFilename = locFilename;
        loader = new GFFLoader(gffFilename);
    }

    /**
     * Load the positions and values in from the locations file.
     */
    void readLocFile() throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new FileReader(locFilename));
        locValues = new ConcurrentHashMap<>();
        int oldPos = 0;
        String line;
        while ((line=reader.readLine())!=null) {
            if (line.startsWith("#")) {
                continue; // comment
            } else if (line.toLowerCase().startsWith("contig")) {
                continue; // heading
            }
            String[] parts = line.split("\t");
            String contig = parts[0];
            int pos = Integer.parseInt(parts[1]);
            String value = parts[10];
            if (pos!=oldPos) {
                oldPos = pos;
                locValues.put(contig+":"+pos, value);
            }
        }
        reader.close();
    }

    /**
     * Search for overlap between the locValues and the genes.
     */
    void search() {
        locValues.forEach(1, this);
    }

    /**
     * BiConsumer method to do the searching.
     */
    public void accept(String loc, String val) {
        // perform the search and add result to the geneLocValues map
        String[] parts = loc.split(":");
        String contig = parts[0];
        int pos = Integer.parseInt(parts[1]);
        Location location = new Location(pos,pos);
        try {
            FeatureList overlapping = loader.search(contig, location);
            if (overlapping.size()>0) {
                geneLocValues.put(loc, val);
            }
        } catch (Exception e) {
            // do nothing
        }
    }

    /**
     * Print the results in wig format.
     */
    void printWig() {
        String oldContig = "";
        for (String loc : geneLocValues.keySet()) {
            String value = geneLocValues.get(loc);
            String[] parts = loc.split(":");
            String contig = parts[0];
            int pos = Integer.parseInt(parts[1]);
            if (!contig.equals(oldContig)) {
                oldContig = contig;
                System.out.println("variableStep chrom="+contig);
            }
            System.out.println(pos+"\t"+value);
        }
    }
    
    /**
     * Command-line operation.
     */
    public static void main(String[] args) throws Exception {
        if (args.length!=2) {
            System.out.println("Usage GFFSearcher <gff-file> <locations-file>");
            System.exit(0);
        }
        String gffFilename = args[0];
        String locFilename = args[1];
        
        GFFSearcher searcher = new GFFSearcher(gffFilename, locFilename);
        searcher.readLocFile();
        searcher.search();
        searcher.printWig();
    }
}
