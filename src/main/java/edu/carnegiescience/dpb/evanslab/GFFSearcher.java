package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.FileReader;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;

/**
 * Searches a GFF file for features that span given locations.
 *
 * Location file has columns:
 * contig	start	REF	ALT	a	b	c	d	size	p	mlog10p	signif
 * 0            1       2       3       4       5       6       7       8       9       10      11
 *
 * Writes wig file format.
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
        int oldPos = 0;
        String oldContig = "";
        String line;
        while ((line=reader.readLine())!=null) {
            if (line.startsWith("#")) {
                continue; // comment
            } else if (line.toLowerCase().startsWith("contig")) {
                continue; // heading
            }
            String[] parts = line.split("\t");
            String contig = parts[0];
            if (!contig.equals(oldContig)) {
                System.out.println("variableStep chrom="+contig);
                oldContig = contig;
            }
            int pos = Integer.parseInt(parts[1]);
            String mlog10p = parts[10];
            if (pos!=oldPos) {
                oldPos = pos;
                Location location = new Location(pos,pos);
                FeatureList overlapping = loader.search(contig, location);
                if (overlapping.size()>0) {
                    FeatureI feature = overlapping.get(0);
                    System.out.println(pos+"\t"+mlog10p);
                }
            }
        }
        reader.close();
    }

}
