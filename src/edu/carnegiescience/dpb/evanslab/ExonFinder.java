package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.FileReader;

/**
 * Run through the GFF file, searching for SNP records that fall within the exon boundaries. If so, update with the exon data.
 */
public class ExonFinder {

    public static void main(String[] args) {

        if (args.length!=1) {
            System.out.println("Usage: ExonFinder <gff-file>");
            System.exit(1);
        }

        try {

            // get database connection from properties file
            DB db1 = new DB();
            DB db2 = new DB();
            
            // spin through GFF file
            BufferedReader in = new BufferedReader(new FileReader(args[0]));
            String line;
            while ((line=in.readLine())!=null) {
                if (!line.startsWith("#")) {
                    GFFExonRecord rec = new GFFExonRecord(line);
                    if (rec.isExon()) {
                        db1.executeQuery("SELECT * FROM gbscounts WHERE contig='"+rec.chrom+"' AND loc>="+rec.start+" AND loc<="+rec.end);
                        while (db1.rs.next()) {
                            String contig = db1.rs.getString("contig");
                            int loc = db1.rs.getInt("loc");
                            db2.executeUpdate("UPDATE gbscounts SET exonid='"+rec.exonId+"',start="+rec.start+",\"end\"="+rec.end+",strand='"+rec.strand+"' WHERE contig='"+contig+"' AND loc="+loc);
                            System.out.println(rec.toString()+"\t"+contig+":"+loc);
                        }
                    }
                }
            }

            // not necessary but looks good
            db1.close();
            db2.close();

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }

}





    
