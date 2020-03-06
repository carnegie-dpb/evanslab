package edu.carnegiescience.dpb.evanslab;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.text.DecimalFormat;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import org.mskcc.cbio.portal.stats.FisherExact;

/**
 * Count the HET, HOM, No-calls per sample groups (A and B) from a VCF file.
 */
public class SNPCounter {

    static DecimalFormat phredFormat = new DecimalFormat("#.00");
    static DecimalFormat oddsFormat = new DecimalFormat("0.00");

    public static void main(String[] args) throws FileNotFoundException, IOException {

        if (args.length!=4) {
            System.err.println("Usage: SNPCounter <VCF file> <A samples list> <B samples list> <maxSize>");
            System.exit(1);
        }

        // input parameters
        // NOTE: vcf index will be assumed to be vcfFilename+".tbi" by the htsjdk routine
        String vcfFilename = args[0];
        String aSamplesFilename = args[1];
        String bSamplesFilename = args[2];
        int maxSize = Integer.parseInt(args[3]);

        FisherExact fisherExact = new FisherExact(maxSize);

        File aSamplesFile = new File(aSamplesFilename);
        if (!aSamplesFile.exists()) {
            System.err.println("Error: file "+aSamplesFile.getName()+" does not exist.");
            System.exit(1);
        }
        
        File bSamplesFile = new File(bSamplesFilename);
        if (!bSamplesFile.exists()) {
            System.err.println("Error: file "+bSamplesFile.getName()+" does not exist.");
            System.exit(1);
        }
        
        File vcfFile = new File(vcfFilename);
        if (!vcfFile.exists()) {
            System.err.println("Error: file "+vcfFile.getName()+" does not exist.");
            System.exit(1);
        }

        // load each group of samples into lists
        List<String> aSamplesList = new ArrayList<String>();
        List<String> bSamplesList = new ArrayList<String>();
        BufferedReader samplesReader = new BufferedReader(new FileReader(aSamplesFile));
        String line = null;
        while ((line=samplesReader.readLine())!=null) {
            aSamplesList.add(line);
        }
        samplesReader = new BufferedReader(new FileReader(bSamplesFile));
        line = null;
        while ((line=samplesReader.readLine())!=null) {
            bSamplesList.add(line);
        }

        // load the VCF file into a reader and header
        VCFFileReader vcfReader = new VCFFileReader(vcfFile, true);
        VCFHeader vcfHeader = vcfReader.getFileHeader();

        // output header
        System.out.println("contig"+"\t"+"position"+"\t"+"qual"+"\t"+"ref"+"\t"+"alt"+"\t"+
                           "aHomRef"+"\t"+"aHomVar"+"\t"+"aHet"+"\t"+
                           "BHomRef"+"\t"+"BHomVar"+"\t"+"bHet"+"\t"+
                           "OR");
        
        // spin through the VCF file
        Iterator<VariantContext> vcIterator = vcfReader.iterator();
        while (vcIterator.hasNext()) {

            VariantContext vc = vcIterator.next();

            // position data
            String contig = vc.getContig();
            int position = vc.getStart();
            String ref = vc.getReference().getBaseString();
            String alt = "";
            for (Allele a : vc.getAlternateAlleles()) {
                alt += a.getBaseString();
            }

            // A group
            int aNC = 0;
            int aHomRef = 0;
            int aHomVar = 0;
            int aHet = 0;
            for (String sample : aSamplesList) {
                Genotype g = vc.getGenotype(sample);
                if (g.isHomRef()) {
                    aHomRef++;
                } else if (g.isHomVar()) {
                    aHomVar++;
                } else if (g.isHet()) {
                    aHet++;
                } else if (g.isNoCall()) {
                    aNC++;
                } else {
                    System.out.println(sample+" has UNKNOWN genotype");
                }
            }

            // B group
            int bNC = 0;
            int bHomRef = 0;
            int bHomVar = 0;
            int bHet = 0;
            for (String sample : bSamplesList) {
                Genotype g = vc.getGenotype(sample);
                if (g.isHomRef()) {
                    bHomRef++;
                } else if (g.isHomVar()) {
                    bHomVar++;
                } else if (g.isHet()) {
                    bHet++;
                } else if (g.isNoCall()) {
                    bNC++;
                } else {
                    System.out.println(sample+" has UNKNOWN genotype");
                }
            }

            // stats
            double oddsRatio = ( (double)bHet/(double)bHomRef ) / ( (double)aHet/(double)aHomRef );
            double p = fisherExact.getP(aHomRef, aHet, bHomRef, bHet);

            // output line
            System.out.println(contig+"\t"+position+"\t"+phredFormat.format(vc.getPhredScaledQual())+"\t"+ref+"\t"+alt+"\t"+
                               aHomRef+"\t"+aHomVar+"\t"+aHet+"\t"+
                               bHomRef+"\t"+bHomVar+"\t"+bHet+"\t"+
                               oddsFormat.format(oddsRatio)+"\t"+p);
        }

    }


}
