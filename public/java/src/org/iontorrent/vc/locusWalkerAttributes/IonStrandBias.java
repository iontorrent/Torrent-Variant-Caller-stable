package org.iontorrent.vc.locusWalkerAttributes;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: ionadmin
 * Date: 4/13/12
 * Time: 11:56 AM
 * To change this template use File | Settings | File Templates.
 */
public class IonStrandBias {
    Map<String, AlignmentContext> stratifiedContexts;
    VariantContext vc;
    VariantCallContext vcc;
    ReferenceContext refContext;

    // Contains the per strand read count.  String in pair is allele, string in map is sample name.
    Map<String,List<Pair<String,Pair<Integer, Integer>>>> stratifiedAllelesStrandReadCounts;

    public IonStrandBias(Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, VariantCallContext vcc, ReferenceContext refContext)  {
        this.stratifiedContexts = stratifiedContexts;
        this.vc = vc;
        this.vcc = vcc;
        this.refContext = refContext;
        this.stratifiedAllelesStrandReadCounts = new HashMap<String,List<Pair<String,Pair<Integer, Integer>>>>();
    }

    public void determineStrandReadCounts()
    throws Exception {
        // Position of the current locus
        long startPos = vcc.getStart();
        long endPos = vcc.getEnd();

        // Add reference window around variant
        long refWinStart = refContext.getWindow().getStart();
        String refBasesInWin = new String(refContext.getBases());
        // and determine reference string
        int startIndex = (int)(startPos - refWinStart);
        int endIndex = (int)(startIndex + (endPos - startPos) + 1);
        String locusContextRefBases = refBasesInWin.substring(startIndex,endIndex);
        determineStrandReadCounts(locusContextRefBases);
    }

    private Pair<Integer, Integer> getPosNegStrandCounts (List<GATKSAMRecord> reads) {
        if(reads==null) return new Pair<Integer, Integer>(0,0);
        Pair<Integer, Integer> posNegStrandCounts;
        Integer posStrandReads = 0;
        Integer negStrandReads = 0;
        for(SAMRecord read : reads) {
            if(read.getReadNegativeStrandFlag()) {
                posStrandReads++;
            } else {
                negStrandReads++;
            }
        }
        posNegStrandCounts = new Pair<Integer, Integer>(posStrandReads,negStrandReads);
        return posNegStrandCounts;
    }

    public void determineStrandReadCounts(String locusContextRefBases)
    throws Exception {
        if (locusContextRefBases=="") {
            determineStrandReadCounts();
            return;
        }
        for (Map.Entry<String,AlignmentContext> stratifiedContext : stratifiedContexts.entrySet()) {
            AlignmentContext ac = stratifiedContext.getValue();
            List<Pair<String,Pair<Integer, Integer>>> allelesStrandReadCounts;
            allelesStrandReadCounts = new ArrayList<Pair<String,Pair<Integer, Integer>>>();
            if (ac.hasBasePileup()) {
                List<Pair<Byte,List<GATKSAMRecord>>> basesNReads;
                basesNReads = ac.getBasePileup().getEventBaseWithReadList(locusContextRefBases.getBytes()[0]);
                for(Pair<Byte,List<GATKSAMRecord>> baseNReads : basesNReads) {
                    String allele = ""+((char)((byte)baseNReads.getFirst()));

                    Pair<Integer, Integer> posNegStrandCounts;
                    posNegStrandCounts = getPosNegStrandCounts(baseNReads.getSecond());
                    allelesStrandReadCounts.add(new Pair<String, Pair<Integer, Integer>>(allele,posNegStrandCounts));
                }
            } else {
                List<Pair<String, List<GATKSAMRecord>>> locusSeqNReads;
                byte [] ref = vcc.getReference().getBases();
                if (vcc.isSimpleInsertion()) {
                    //System.err.printf("1D ref length: " + ref.length + "\n");
                    ref = new byte[] {'1','D'};
                } if (vcc.isSimpleDeletion()) {
                    //System.err.printf("+X ref length: " + ref.length + "\n");
                    ref = new byte[] {'+',ref[ref.length-1]};
                }
                locusSeqNReads = ac.getExtendedEventPileup().getEventStringWithReadList(ref,true);
                for(Pair<String, List<GATKSAMRecord>> alleleNReads : locusSeqNReads) {
                    Pair<Integer, Integer> posNegStrandCounts;
                    posNegStrandCounts = getPosNegStrandCounts(alleleNReads.getSecond());
                    allelesStrandReadCounts.add(new Pair<String, Pair<Integer, Integer>>
                            (alleleNReads.getFirst(),posNegStrandCounts));
                }
            }
            stratifiedAllelesStrandReadCounts.put(stratifiedContext.getKey(), allelesStrandReadCounts);
        }
    }

    public String getStrandReadCountsString() {
        String strandCountString = "";
        for (Map.Entry<String,List<Pair<String,Pair<Integer, Integer>>>> allelesStrandReadCounts :
                stratifiedAllelesStrandReadCounts.entrySet() ) {
            strandCountString += allelesStrandReadCounts.getKey() + ":";
            for(Pair<String,Pair<Integer, Integer>> alleleNStrandReadCounts : allelesStrandReadCounts.getValue()) {
                // Print out allele, read counts for + strand, and read counts for - strand.
                strandCountString += alleleNStrandReadCounts.getFirst() + ",";
                strandCountString += alleleNStrandReadCounts.getSecond().getFirst() + ",";
                strandCountString += alleleNStrandReadCounts.getSecond().getSecond() + ",";
            }
            strandCountString += ":";
        }
        if (strandCountString.length()<=1) return strandCountString;
        return strandCountString.substring(0,strandCountString.length()-2);
    }
}
