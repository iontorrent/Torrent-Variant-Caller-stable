package org.iontorrent.vc.scoring;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.IonSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.iontorrent.sam2flowgram.flowalign.FlowgramAlignment;
import org.iontorrent.sam2flowgram.util.SamToFlowgramAlignUtil;
import org.iontorrent.sam2flowgram.util.Tuple;

import java.io.PrintStream;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/20/12
 * Time: 2:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class BayesianScorer {

    /*
    Print stream for the reads' flowgram alignments in the context of the variant.
    */
    protected PrintStream flowAlignContextPileupStrm;

    /*
    The actual API call class
     */
    protected BayesianScorerAPI bayAPI;

    // Overall score
    // ToDo, it's the maximum right now, but perhaps what to do something with the
    // individual variant scores.
    Double bayesianLocusScore;
    
    // Internal thread id given by the API so that simul
    private Integer threadId;

    // For debugging can turn this off
    Boolean doAPICall = true;

    // Reference positions
    long refStartPosition;
    long refEndPosition;
    
    // Variant info
    long[] obsSeqPositions; // = new long[1000];
    int curVariantNum;  // 0 is for the reference

    // Alignment Context Needed
    int numFlowWindowSize = 20;

    /**
     * BayesianScorer extracts out the information from
     * @param flowAlignContextPileupStrm
     * @param doAPICall
     */
    public BayesianScorer(BayesianScorerAPI bayAPI, PrintStream flowAlignContextPileupStrm, Boolean doAPICall) {
        this.bayAPI = bayAPI;
        this.flowAlignContextPileupStrm = flowAlignContextPileupStrm;
        this.doAPICall = doAPICall;

    }

    public Double getBayesianLocusScore() throws Exception {
        if(bayesianLocusScore==null) {
            throw new Exception("No bayesianLocusScore, perhaps because there were no reads.");
        }
        return bayesianLocusScore;
    }

    public void reset() {
        curVariantNum = 1;
        obsSeqPositions = new long[1000];
    }
    
    /**
     * Adds reference context for the score
     * @param refStartPosition The genomic start position of that variant
     * @param refBases The actual reference sequence
     */
    void addReference(long refStartPosition, String refBases) {
        this.refStartPosition = refStartPosition;
        refEndPosition = refStartPosition + refBases.length() - 1;
        if(flowAlignContextPileupStrm!=null) {
            flowAlignContextPileupStrm.printf("addReference:" + refBases.length() + ";" + refBases + "\n");
        }
        if(doAPICall) {
            threadId = bayAPI.addRef(refBases.length(),refBases.getBytes());
            flowAlignContextPileupStrm.printf("addReference returned:" + threadId + "\n");
        } else {
            threadId = 0;
        }
    }

    private void checkThreadId() throws Exception {
        if (threadId==null) {
            throw new Exception("addReference must be called first before adding any additional information.\n");
        }
    }
    /**
     * A set of variants are considered together (b/c they are at or near the same position) for the score
     * @param obsSeqPosition The genomic position of the variant
     * @param observedSeq The observed sequence, can be reference or variant.
     * @throws Exception
     */
    void addObservedSeq(long obsSeqPosition, String observedSeq) throws Exception {
        checkThreadId();
        if(obsSeqPosition <refStartPosition) {
            throw new Exception(String.format("Given variant position %d is not in reference context given (%d-%d)\n", obsSeqPosition, refStartPosition, refEndPosition));
        }
        obsSeqPositions[curVariantNum++] = obsSeqPosition;

        int locPosRelRef = (int)(obsSeqPosition - refStartPosition);
        if(flowAlignContextPileupStrm!=null) {
            flowAlignContextPileupStrm.printf("addVariant:" + threadId + ";" + locPosRelRef + ";" + observedSeq + "\n");
        }
        //System.err.printf("addObservedSeq:" + locPosRelRef + ";" + observedSeq + "\n");
        observedSeq += '\000';  // API requires null termination
        if(doAPICall) bayAPI.addVariant(threadId,locPosRelRef,observedSeq.getBytes());
    }

    /**
     * For a given variant, add a read.  Call multiple times to get all the reads needed for the score.
     * @param variantNum The index to the variant, in order of the addVariant calls, reserving 0 for the reference
     * @param read The Ion specialized picard read object
     * @throws Exception
     */
    void addRead(int variantNum, IonSAMRecord read) throws Exception {
        checkThreadId();
        if(variantNum > curVariantNum ) {
            throw new Exception(String.format("Invalid curVariantNum\n"));
        }
        long variantPosition = obsSeqPositions[0];
        if(variantNum>0) variantPosition = obsSeqPositions[variantNum - 1];

        int strandInt = read.getReadNegativeStrandFlag() ? 0 : 1;  // 0 on neg. strand, 1 on positive
        FlowgramAlignment subAlign = read.flowAlign.extractPartialRegion(variantPosition, numFlowWindowSize,true);
        FlowgramAlignment something;
        
        if (flowAlignContextPileupStrm!=null) {
            flowAlignContextPileupStrm.printf("addRead:" + threadId + ";" + variantNum + ";" + subAlign.variantStart + ";" + subAlign.length + ";");

            int subFlowNum = 0;
            for(byte flowBase : subAlign.flowOrder) {
                flowAlignContextPileupStrm.print(SamToFlowgramAlignUtil.DNA[flowBase]);
                if(++subFlowNum==subAlign.length) break;
            }
            flowAlignContextPileupStrm.printf(";");
            subFlowNum = 0;
            for(int intensity : subAlign.qseq){
                flowAlignContextPileupStrm.print(intensity);
                flowAlignContextPileupStrm.print(',');
                if(++subFlowNum==subAlign.length) break;
            }
            subFlowNum = 0;
            for(char flowAlignSymbol : subAlign.aln) {
                flowAlignContextPileupStrm.print(flowAlignSymbol);
                flowAlignContextPileupStrm.print(',');
                if(++subFlowNum==subAlign.length) break;
            }

            flowAlignContextPileupStrm.printf(";" + strandInt);
            flowAlignContextPileupStrm.println();
        }

        // Make subalign arrays for the API
        byte[] flowBases = new byte[subAlign.length];
        int[] flowInten = new int[subAlign.length];
        byte[] flowAlignSymbols = new byte[subAlign.length];

        int subFlowNum = 0;
        for(byte flowBase : subAlign.flowOrder) {
            flowBases[subFlowNum] = (byte)SamToFlowgramAlignUtil.DNA[flowBase];
            flowAlignSymbols[subFlowNum] = (byte)subAlign.aln[subFlowNum];
            if(++subFlowNum==subAlign.length) break;
        }
        subFlowNum = 0;
        for(int intensity : subAlign.qseq){
            flowInten[subFlowNum] = intensity;
            if(++subFlowNum==subAlign.length) break;
        }
        subFlowNum = 0;
        for(char flowAlignSymbol : subAlign.aln) {
            flowAlignSymbols[subFlowNum] = (byte)flowAlignSymbol;
            if(++subFlowNum==subAlign.length) break;
        }

        //System.err.printf("flowBases/subAlign lengths: " + flowBases.length + "," + subAlign.length + "\n");
        if(doAPICall)
            bayAPI.addRead(threadId, variantNum,subAlign.variantStart,subAlign.length,
                flowBases,flowInten,flowAlignSymbols,strandInt);
    }

    public double calcScore() throws Exception {
        checkThreadId();

        double bayScore = 3.1415;  // dummy value if doAPICall is false for debugging
        if(doAPICall) {
            bayScore = bayAPI.calScore(threadId);
        }

        if (flowAlignContextPileupStrm!=null) {
            flowAlignContextPileupStrm.printf("calScore:" + threadId + ";" + bayScore + "\n");
        }
        return bayScore;
    }

    public void calcScoreNSetMax() throws Exception {
        double curVariantScore = calcScore();
        // Keep the maximum score
        if (bayesianLocusScore==null || curVariantScore>bayesianLocusScore) {
            bayesianLocusScore = curVariantScore;
        }
    }

    public void addVariantContext(Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, VariantCallContext vcc, ReferenceContext refContext)
    throws Exception {
        reset();
        // Position of the current locus
        long startPos = vcc.getStart();
        long endPos = vcc.getEnd();

        // Add reference window around variant
        long refWinStart = refContext.getWindow().getStart();
        String refBasesInWin = new String(refContext.getBases());
        addReference(refWinStart,refBasesInWin);
        // and determine reference string
        int startIndex = (int)(startPos - refWinStart);
        int endIndex = (int)(startIndex + (endPos - startPos) + 1);
        String locusContextRefBases = refBasesInWin.substring(startIndex,endIndex);

        // Loop through all the alleles
        bayesianLocusScore = null;
        Tuple<Integer> alignIndexAndOffset = new Tuple<Integer>(-1,-1);

        //TODO, The stratifedContext doesn't seem to get the reads right.  This is hacked upstream to use
        // rawContext from the walker, but of course this messes up separating things by sample.
        // This is ok, because BayesianScorer doesn't support multiple samples yet anyway.
        for (Map.Entry<String,AlignmentContext> stratifiedContext : stratifiedContexts.entrySet()) {
            //String sampleName = stratifiedContext.getKey();

            AlignmentContext ac = stratifiedContext.getValue();
            int variantNum = 0;  // 0 is actually the reference
            if (ac.hasBasePileup()) {
                List<Pair<Byte,List<GATKSAMRecord>>> basesNReads;
                basesNReads = ac.getBasePileup().getEventBaseWithReadList(locusContextRefBases.getBytes()[0]);
                Boolean isFirst = true;
                for(Pair<Byte,List<GATKSAMRecord>> baseNReads : basesNReads) {
                    String allele = ""+((char)((byte)baseNReads.getFirst()));
                    if (isFirst) {
                        isFirst = false;
                        obsSeqPositions[0] = startPos;
                    } else
                        this.addObservedSeq(startPos,allele);
                    addReadSet(baseNReads.getSecond(), variantNum);
                    variantNum++;
                    calcScoreNSetMax();
                }
            } else {            
                List<Pair<String, List<GATKSAMRecord>>> locusSeqNReads;
                //ToDo: locusContextRefBases.getBytes() is different than vcc.getReference().getBases();
                //ToDo: The translation below seems to work.  Couldn't find a util translation function for this
                //ToDo: so reversed engineered from the output.
                byte [] ref = vcc.getReference().getBases();
                if (vcc.isSimpleInsertion()) {
                    //System.err.printf("1D ref length: " + ref.length + "\n");
                    ref = new byte[] {'1','D',0};
                } else if (vcc.isSimpleDeletion()) {
                    //System.err.printf("+X ref length: " + ref.length + "\n");
                    ref = new byte[] {'+',ref[ref.length-1],0};
                }

                flowAlignContextPileupStrm.println("Total number of reads in pileup:" + ac.getExtendedEventPileup().getReads().size());
                locusSeqNReads = ac.getExtendedEventPileup().getEventStringWithReadList(ref,true);
                Boolean isFirst = true;
                for(Pair<String, List<GATKSAMRecord>> alleleNList : locusSeqNReads) {
                    if (isFirst) {
                        isFirst = false;
                        obsSeqPositions[0] = startPos;
                    } else
                        this.addObservedSeq(startPos,alleleNList.getFirst());
                    addReadSet(alleleNList.getSecond(), variantNum);
                    variantNum++;
                    calcScoreNSetMax();
                }
            }

        }
        flowAlignContextPileupStrm.printf("bayesianLocusScore = " + bayesianLocusScore + "\n");
        if(bayesianLocusScore==null) bayesianLocusScore = 0.0;
        flowAlignContextPileupStrm.printf("Finished:" + threadId + "\n");
        bayAPI.finished(threadId);
    }

    private void addReadSet(List<GATKSAMRecord> contextReads, int variantNum) throws Exception {
        if (contextReads==null) {
            return;
        }
        flowAlignContextPileupStrm.println("Size of contextReads " + contextReads.size() );
        for(GATKSAMRecord read : contextReads) {
            IonSAMRecord contextRead = (IonSAMRecord)read;
            addRead(variantNum,contextRead);
        }
    }


}
