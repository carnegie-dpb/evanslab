package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.text.DecimalFormat;
import java.util.List;


/**
 * Run through the ref counts file and fill in Ref lines in the requested SNP file for the requested sample ID.
 *
 * The tab-separated SNP file format is as follows:
 *
 * annoPos  gene                            contig  pos1    pos2    ref     alt     genotype quality MQ      ref-depth  alt-depth
 * upstream transcript:Zm00001d027578_T001  1       8362870 8362870 C       A       het      26      25      1          5
 *
 * The tab-seperated Ref file format is as follows:
 *
 * contig   pos   S1       S2       S3       ...       Sn
 * 1        33879 NC       Het      R:4:44.2 ...       Hom
 *
 * So for lines in the Ref file that are missing from the SNP file that have Ref calls for the desired sample, we insert:
 *
 * .        .                               1       33879  33879    .       .       ref      44.2    .       4         .
 *
 * Only the R calls are used here since the SNP file already has the het/hom calls.
 *
 */
public class RefSNPMerger {

    static DecimalFormat df = new DecimalFormat("0.0");

    public static void main(String[] args) {

        if (args.length!=4) {
            System.out.println("Usage: RefSNPMerger <Ref file> <SNP dir> <SNP file suffix> <Sample ID>");
            System.out.println("Example: RefSNPMerger ~/refcounts.txt ~/snps avinput.variant_function.snp VT1");
            System.exit(1);
        }
        
        try {
            
            // input parameters
            String refFileName = args[0];
            String snpDirName = args[1];
            String snpFileSuffix = args[2];
            String sampleId = args[3];

            File refFile = new File(refFileName);
            if (!refFile.exists()) {
                System.err.println("Error: "+refFile.getName()+" does not exist.");
                System.exit(1);
            }
            
            File snpDir = new File(snpDirName);
            if (!snpDir.isDirectory()) {
                System.err.println("Error: "+snpDir.getName()+" is not a directory.");
                System.exit(1);
            }

            File snpFile = new File(snpDir, sampleId+"."+snpFileSuffix);
            if (!snpFile.exists()) {
                System.err.println("Error: "+snpFile.getName()+" does not exist.");
                System.exit(1);
            }

            // load samples from the Ref header
            BufferedReader refReader = new BufferedReader(new FileReader(refFile));
            String refLine = refReader.readLine();
            List<String> refSamples = RefRecord.getSamples(refLine);

            // output the SNP header line
            BufferedReader snpReader = new BufferedReader(new FileReader(snpFile));
            String snpLine = snpReader.readLine();
            System.out.println(snpLine);

            // spin through the SNP records, spitting out intervening Ref records
            while ((snpLine=snpReader.readLine())!=null) {
                SNPRecord snpRec = new SNPRecord(snpLine);
                boolean found = false;
                while (!found && (refLine=refReader.readLine())!=null) {
                    RefRecord refRec = new RefRecord(refLine, refSamples);
                    found = (refRec.position==snpRec.pos1);
                    int reads = refRec.getRefReads(sampleId);
                    double quality = refRec.getRefQuality(sampleId);
                    if (!found && reads>0 && quality>0.0) {
                        System.out.println(".\t.\t"+refRec.contig+"\t"+refRec.position+"\t"+refRec.position+"\t.\t.\tref\t"+df.format(quality)+"\t.\t"+reads+"\t.");
                    }
                }
                // spit out SNP line
                System.out.println(snpLine);
            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
