package edu.carnegiescience.dpb.evanslab;

import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;

/**
 * Test the BioJava GFFReader functionality.
 */
public class GFFReaderTest {

    public static void main(String[] args) {

        if (args.length!=1) {
            System.out.println("Usage: GFFReaderTest <GFF-file>");
            System.exit(1);
        }

        try {
            
            FeatureList features = GFF3Reader.read(args[0]);
            for (FeatureI feature : features) {
                System.out.println(feature.toString());
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }

}
