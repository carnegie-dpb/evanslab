package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Iterator;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GFF3Reader;
import org.biojava.nbio.genome.parsers.gff.Location;

/**
 * Scans a sorted BAM file to find read pairs where one is on an exon and one is not.
 * Exon positions are given by a GFF file.
 *
 * @author Sam Hokin
 */
public class BreakpointAnalyzer {

    public static void main(String[] args) throws Exception {
        if (args.length!=2) {
            System.out.println("Usage BreakpointAnalyzer <BAM file> <GFF file>");
            System.exit(0);
        }
        
        String bamFilename = args[0];
        String gffFilename = args[1];

        final SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFilename));
        final FeatureList exons = GFF3Reader.read(gffFilename).selectByType("exon");
        System.err.println("Read "+exons.size()+" exons from "+gffFilename+".");

        SAMRecord record1 = null;
        SAMRecord record2 = null;
        boolean first = true;
        for (Iterator i=reader.iterator(); i.hasNext(); ) {
            if (first) {
                record1 = (SAMRecord) i.next();
            } else {
                record2 = (SAMRecord) i.next();
            }
            if (!first) {
                // check that the reads match
                if (!record1.getReadName().equals(record2.getReadName())) {
                    System.err.println("Whoops! Read pairs don't have same name:");
                    System.err.println(record1.getReadName());
                    System.err.println(record2.getReadName());
                    System.exit(1);
                }
                // locations
                String seqname1 = record1.getContig();
                String seqname2 = record2.getContig();
                Location location1 = createLocation(record1);
                Location location2 = createLocation(record2);
                // find strand-specific overlaps
                FeatureList overlaps1 = exons.selectOverlapping(seqname1, location1, false);
                FeatureList overlaps2 = exons.selectOverlapping(seqname2, location2, false);
                boolean onExon1 = overlaps1.size() > 0;
                boolean onExon2 = overlaps2.size() > 0;
                if (onExon1 && !onExon2) {
                    String bases2 = new String(record2.getReadBases(), StandardCharsets.US_ASCII);
                    System.out.println(record1.getReadName()+"\t"+seqname1
                                       +"\t"+location1.bioStart()+"\t"+location1.bioEnd()+"\t"+location1.bioStrand()
                                       +"\t"+bases2);
                } else if (!onExon1 && onExon2) {
                    String bases1 = new String(record1.getReadBases(), StandardCharsets.US_ASCII);
                    System.out.println(record2.getReadName()+"\t"+seqname2
                                       +"\t"+location2.bioStart()+"\t"+location2.bioEnd()+"\t"+location2.bioStrand()
                                       +"\t"+bases1);
                }
            }
            // toggle read pair flag
            first = !first;
        }
        reader.close();
    }

    /**
     * Create a biojava Location from an htsjdk SAMRecord
     */
    static Location createLocation(SAMRecord record) {
        if (record.getReadNegativeStrandFlag()) {
            return new Location(-record.getEnd(), -record.getStart());
        } else {
            return new Location(record.getStart(), record.getEnd());
        }
    }
}
