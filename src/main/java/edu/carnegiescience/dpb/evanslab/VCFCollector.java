package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.Map;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import org.mskcc.cbio.portal.stats.FisherExact;

/**
 * Collect REF vs HET counts across the samples in a VCF file according to their nontwin/twin condition.
 *
 * NOTE: we ignore any calls with zero forward or reverse reads. These are biased by sequencing algorithm.
 */
public class VCFCollector {

    static final String[] nts = {"VT1", "VT2", "VT4", "VT5", "VT6", "VT7", "VT8", "VT9",
                                      "VT10", "VT11", "VT12", "VT13", "VT14", "VT15", "VT16", "VT17", "VT18", "VT19",
                                      "VT20", "VT22", "VT23", "VT24", "VT26", "VT27", "VT28", "VT29",
                                      "VT30", "VT31", "VT32", "VT33", "VT34", "VT35", "VT36", "VT37", "VT38", "VT39",
                                      "VT40", "VT41", "VT42", "VT43", "VT44", "VT45", "VT46", "VT47", "VT48"};

    static final String[] tws = {"VT49",
                                   "VT50", "VT51", "VT52", "VT53", "VT54", "VT55", "VT56", "VT57", "VT58", "VT59",
                                   "VT60", "VT61", "VT62", "VT63", "VT64", "VT65", "VT66", "VT67", "VT68", "VT69",
                                   "VT70", "VT71", "VT72", "VT73", "VT74", "VT75", "VT76", "VT77", "VT78", "VT79",
                                   "VT80", "VT81", "VT82", "VT83", "VT84", "VT85", "VT86", "VT87", "VT88", "VT89",
                                   "VT90", "VT91", "VT92", "VT93", "VT94", "VT95", "VT96"};

    static List<String> ntList = Arrays.asList(nts);
    static List<String> twList = Arrays.asList(tws);

    static int ntCount = ntList.size();
    static int twCount = twList.size();

    public static void main(String[] args) {

        if (args.length!=1) {
            System.err.println("Usage: VCFCollector <VCF.gz file>");
            System.exit(1);
        }

        File vcfFile = new File(args[0]);

        FisherExact fisherExact = new FisherExact(ntCount+twCount);
        
        // header
        System.out.println("contig\tpos\tref\talt\t"+
                           "TWRF\tTWRR\tTWAF\tTWAR\t" +
                           "NTRF\tNTRR\tNTAF\tNTAR\t" +
                           "TWRef\tTWHet\tTWHom\t" +
                           "NTRef\tNTHet\tNTHom\t" +
                           "p\tmlog10p\tOR");

        VCFFileReader reader = new VCFFileReader(vcfFile, false);
        VCFHeader header = reader.getFileHeader();
        List<String> samples = header.getSampleNamesInOrder();

        Iterator<VariantContext> vcIterator = reader.iterator();
        while (vcIterator.hasNext()) {
            VariantContext vc = vcIterator.next();
            String contig = vc.getContig();
            int pos = vc.getStart();
            Allele ref = vc.getReference();
            List<Allele> alts = vc.getAlternateAlleles();
            Allele alt = alts.get(0);
            // limit ourselves to single variant SNPs
            if (alts.size()==1 && ref.toString().length()==2 && alt.toString().length()==1) {
                // TWINS
                int twHet = 0;
                int twHom = 0;
                int twRF = 0;
                int twRR = 0;
                int twAF = 0;
                int twAR = 0;
                for (String sampleName : twList) {
                    Genotype g = vc.getGenotype(sampleName);
                    // g won't have attributes if no call
                    if (g.hasExtendedAttribute("ADF") && g.hasExtendedAttribute("ADR")) {
                        String[] adfParts = ((String) g.getExtendedAttribute("ADF")).split(",");
                        String[] adrParts = ((String) g.getExtendedAttribute("ADR")).split(",");
                        int refFor = Integer.parseInt(adfParts[0]);
                        int refRev = Integer.parseInt(adrParts[0]);
                        int altFor = Integer.parseInt(adfParts[1]);
                        int altRev = Integer.parseInt(adrParts[1]);
                        // require reads on both strands if not both zero
                        boolean discardRef = (refFor+refRev)>0 && (refFor*refRev)==0;
                        boolean discardAlt = (altFor+altRev)>0 && (altFor*altRev)==0;
                        if (!discardRef && !discardAlt) {
                            twRF += refFor;
                            twRR += refRev;
                            twAF += altFor;
                            twAR += altRev;
                            if (g.isHet()) {
                                twHet++;
                            } else if (g.isHomVar()) {
                                twHom++;
                            } else {
                                System.err.println("ERROR at "+contig+":"+pos+":"+ref+"/"+alt+":"+sampleName+":"+g.toString());
                            }
                        }
                    }
                }
                // NONTWINS
                int ntHet = 0;
                int ntHom = 0;
                int ntRF = 0;
                int ntRR = 0;
                int ntAF = 0;
                int ntAR = 0;
                for (String sampleName : ntList) {
                    Genotype g = vc.getGenotype(sampleName);
                    // g won't have attributes if no call
                    if (g.hasExtendedAttribute("ADF") && g.hasExtendedAttribute("ADR")) {
                        String[] adfParts = ((String) g.getExtendedAttribute("ADF")).split(",");
                        String[] adrParts = ((String) g.getExtendedAttribute("ADR")).split(",");
                        int refFor = Integer.parseInt(adfParts[0]);
                        int refRev = Integer.parseInt(adrParts[0]);
                        int altFor = Integer.parseInt(adfParts[1]);
                        int altRev = Integer.parseInt(adrParts[1]);
                        boolean discardRef = (refFor+refRev)>0 && (refFor*refRev)==0;
                        boolean discardAlt = (altFor+altRev)>0 && (altFor*altRev)==0;
                        if (!discardRef && !discardAlt) {
                            ntRF += refFor;
                            ntRR += refRev;
                            ntAF += altFor;
                            ntAR += altRev;
                            if (g.isHet()) {
                                ntHet++;
                            } else if (g.isHomVar()) {
                                ntHom++;
                            } else {
                                System.err.println("ERROR at "+contig+":"+pos+":"+ref+"/"+alt+":"+sampleName+":"+g.toString());
                            }
                        }
                    }
                }
                if (twHet>0 || twHom>0 || ntHet>0 || ntHom>0) {
                    // STATS
                    int twRef = twCount - twHet - twHom;
                    int ntRef = ntCount - ntHet - ntHom;
                    double p = fisherExact.getP(twRef, twHet, ntRef, ntHet);
                    double mlog10p = -Math.log10(p);
                    double or = ( (double)twHet*(double)ntRef ) / ( (double)ntHet*(double)twRef );
                    // OUTPUT
                    System.out.println(contig+"\t"+pos+"\t"+ref+"\t"+alt+"\t"+
                                       twRF+"\t"+twRR+"\t"+twAF+"\t"+twAR+"\t" +
                                       ntRF+"\t"+ntRR+"\t"+ntAF+"\t"+ntAR+"\t" +
                                       twRef+"\t"+twHet+"\t"+twHom+"\t"+
                                       ntRef+"\t"+ntHet+"\t"+ntHom+"\t"+
                                       p+"\t"+mlog10p+"\t"+or);
                }
            }
        }
    }
}
