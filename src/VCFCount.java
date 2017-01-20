package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Map;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Run through all the VCF files in the given directory and count the alt alleles, incrementing a counter in a database table.
 */
public class VCFCount {

    public static void main(String[] args) {

        if (args.length>2) {
            System.err.println("Usage: VCFCount <VCF-dir> [min(-log10(p))]");
            System.exit(1);
        }

        try {
            
            File dir = new File(args[0]);
            if (!dir.isDirectory()) {
                System.err.println("Error: "+args[0]+" is not a directory.");
                System.exit(1);
            }

            double mlog10pmin = 0.00;
            boolean pFilter = args.length>1;
            if (pFilter) {
                mlog10pmin = Double.parseDouble(args[1]);
                System.out.println("Filtering out values with -log10(p)<"+mlog10pmin);
            }

            List<File> fileList = Arrays.asList(dir.listFiles());

            DB db = new DB();
            
            for (File file : fileList) {
                if (file.getName().endsWith("vcf")) {

                    VCFFileReader reader = new VCFFileReader(file, false);
                    String contig = "";

                    Iterator<VariantContext> vcIterator = reader.iterator();
                    while (vcIterator.hasNext()) {
                        
                        VariantContext vc = vcIterator.next();
                        Allele ref = vc.getReference();
                        List<Allele> alts = vc.getAlternateAlleles();

                        // diagnostic output
                        if (!vc.getContig().equals(contig)) {
                            contig = vc.getContig();
                            System.out.println(file.getName()+":"+contig);
                        }

                        // increment database record
                        if (alts.size()==1) {
                            double mlog10p = -vc.getLog10PError();
                            if (!pFilter || mlog10p>mlog10pmin) {
                                Allele alt = alts.get(0);
                                db.executeQuery("SELECT count(*) FROM gbscounts WHERE contig='"+vc.getContig()+"' AND loc="+vc.getStart()); // assume SNPs
                                db.rs.next();
                                if (db.rs.getInt("count")==0) {
                                    // insert new contig record, will initialize default values
                                    db.executeUpdate("INSERT INTO gbscounts (contig,loc,ref) VALUES ('"+vc.getContig()+"',"+vc.getStart()+",'"+ref.getBaseString()+"')");
                                }
                                // update alt allele count and min -log10p value
                                db.executeUpdate("UPDATE gbscounts SET "+alt.getBaseString()+"="+alt.getBaseString()+"+1, minmlog10p=LEAST(minmlog10p,"+mlog10p+") WHERE contig='"+vc.getContig()+"' AND loc="+vc.getStart());
                            }
                        }
                        
                    }
                    
                }
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
