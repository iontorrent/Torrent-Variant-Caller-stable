package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.StandardVCFWriter;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.IonSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.iontorrent.sam2flowgram.util.SamToFlowgramAlignUtil;
import org.iontorrent.sam2flowgram.util.Tuple;
import org.iontorrent.vc.locusWalkerAttributes.IonStrandBias;
import org.iontorrent.vc.scoring.BayesianScorer;
import org.iontorrent.vc.scoring.BayesianScorerFactory;

import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/14/12
 * Time: 9:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class IonUnifiedGenotyperEngine extends UnifiedGenotyperEngine {
    
    protected PrintStream flowIntensityWriter;
    protected PrintStream flowAlignContextPileupStrm;
    protected StandardVCFWriter vcfWriter;
    
    //API calls out to scorer.
    protected BayesianScorerFactory bayesianScorerFactory;
    protected Integer numThreads = 10; // ToDo, get this from GATK somehow, but now set it to some large possible value
    
    public Boolean displayPerAllele = true;
    
    public IonUnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC) {
        super(toolkit, UAC);
    }

    public IonUnifiedGenotyperEngine(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter,
                                     VariantAnnotatorEngine engine, Set<String> samples, PrintStream flowIntensityWriter,
                                     PrintStream flowAlignContextPileupStrm) {
        super(toolkit, UAC, logger, verboseWriter, engine, samples);
        this.flowIntensityWriter = flowIntensityWriter;
        this.flowAlignContextPileupStrm = flowAlignContextPileupStrm;
        String startupMsg = UAC.IGNORE_FLOW_INTENSITIES ? "Running regular genotyper ignoring flow intensities." :
                "Running flow intensity enhanced genotyper.";
        logger.info(startupMsg);

        //Add VCF line to debug output
        vcfWriter = new StandardVCFWriter(this.flowIntensityWriter, toolkit.getMasterSequenceDictionary(), false);
        bayesianScorerFactory = new BayesianScorerFactory(flowAlignContextPileupStrm,numThreads);
        bayesianScorerFactory.addDebugStream(flowIntensityWriter);
    }
    
    public void writeHeader(VCFHeader vcfHeader) {
        if (flowIntensityWriter!=null)
            vcfWriter.writeHeader(vcfHeader);
    }

    public VariantCallContext calculateLikelihoodsAndGenotypesIon(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        final GenotypeLikelihoodsCalculationModel.Model model = getCurrentGLModel(tracker, refContext, rawContext );
        if( model == null ) {
            return (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(tracker, refContext, null, rawContext) : null);
        }

        //Looks like the String here is coming from  AlignmentContextUtils.splitContextBySampleName(pileup);
        //From tests, it was the same string as the SM tag of the BAM file.
        Map<String, AlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext, model);
        if ( stratifiedContexts == null ) {
            return (UAC.OutputMode == OUTPUT_MODE.EMIT_ALL_SITES && UAC.GenotypingMode == GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES ? generateEmptyContext(tracker, refContext, stratifiedContexts, rawContext) : null);
        }

        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, null, true, model);
        if ( vc == null )
            return null;

        VariantCallContext vcc =
                calculateGenotypes(tracker,refContext,rawContext,stratifiedContexts,vc,model);
        if (vcc == null) return null;

        if(!vcc.shouldEmit) return vcc; //omit vcc that shouldn't be emitted.

        hackStratifiedContexts(stratifiedContexts,rawContext);

        callIonStrandBias(stratifiedContexts,vc,vcc,refContext);
        //debug output
        if(flowIntensityWriter!=null)
            outputLocusFlowInfo(stratifiedContexts, vc, vcc, refContext, rawContext);

        callBayesianScorer(stratifiedContexts,vc,vcc,refContext);

        return vcc;
    }

    //TODO, The stratifedContext doesn't seem to get the reads right.  This is a hack just to use the rawContext which
    //TODO, is bad for multiple samples, but we don't really support mutiple samples yet for TS 2.2.
    private void hackStratifiedContexts(Map<String, AlignmentContext> stratifiedContexts,AlignmentContext rawContext) {
        for (Map.Entry<String,AlignmentContext> stratifiedContext : stratifiedContexts.entrySet()) {
            if(flowIntensityWriter!=null) {
                int numReads = 0;
                numReads = stratifiedContext.getValue().hasBasePileup() ?
                        stratifiedContext.getValue().getBasePileup().getReads().size() :
                        stratifiedContext.getValue().getExtendedEventPileup().getReads().size();
                flowIntensityWriter.println("Original context for '" + stratifiedContext.getKey() + "' had " + numReads + " reads.");
            }
            stratifiedContexts.put(stratifiedContext.getKey(),rawContext);
        }
    }

    public void callBayesianScorer(Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, VariantCallContext vcc, ReferenceContext refContext) {
        //BayesianScorer scorer = new BayesianScorer(flowAlignContextPileupStrm, numThreads);
        //BayesianScorer debug = null;
        try {
            BayesianScorer bayesianScorer = bayesianScorerFactory.createBayesianScorer();
            bayesianScorer.addVariantContext(stratifiedContexts, vc, vcc, refContext);
            if(flowIntensityWriter!=null) {
                BayesianScorer debugScorer =  bayesianScorerFactory.createDebugScorer();
                debugScorer.addVariantContext(stratifiedContexts,vc,vcc,refContext);
            }
            vcc.putAttribute("Bayesian_Score",bayesianScorer.getBayesianLocusScore());
        } catch (Throwable ex) {
            System.out.println("Error in scorer.\n");
            ex.printStackTrace();
        }
    }

    public void callIonStrandBias(Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, VariantCallContext vcc,
                                  ReferenceContext refContext) {
        IonStrandBias strandBias = new IonStrandBias(stratifiedContexts,vc,vcc,refContext);
        try {
            strandBias.determineStrandReadCounts();
            vcc.putAttribute("Strand_Counts",strandBias.getStrandReadCountsString());
         } catch (Exception e) {
            System.out.println("Error in calling determining strand read counts.\n");
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
    public void outputLocusFlowInfo(Map<String, AlignmentContext> stratifiedContexts, VariantContext vc, VariantCallContext vcc,
                                    ReferenceContext refContext, AlignmentContext rawContext) {

        //Call is actually made
        // vcc has the alleles (and if they are called)
        // stratifiedContexts has pointers to the read.
        // So will map through all stratifiedContexts, look up allele in vcc

        // output vcf record
        vcfWriter.add(vcc);

        // Add custom output
        String contig = vcc.getChr();
        long startPos = vcc.getStart();
        long endPos = vcc.getEnd();
        flowIntensityWriter.printf("vcc position: " + contig + ":" + startPos + "-" + endPos + "\n");
        
        long refWinStart = refContext.getWindow().getStart();
        long refWinEnd =  refContext.getWindow().getStop();
        String refBasesInWin =   new String(refContext.getBases());
        flowIntensityWriter.printf("ref window [" + refWinStart + "-" + refWinEnd +
                "]: " + refBasesInWin + "\n");

        int refSlack = 5;
        int startIndex = (int)(startPos - refWinStart);
        int endIndex = (int)(startIndex + (endPos - startPos) + 1);
        String locusContextRefBases = refBasesInWin.substring(startIndex,endIndex);
        String locusContextRefBasesNSlack = refBasesInWin.substring(startIndex-refSlack,endIndex+refSlack);

        flowIntensityWriter.print("ref locus [" + refContext.getLocus().getStart() + "-" + refContext.getLocus().getStop() +
                "]: " + locusContextRefBases + "  Immediate context: " + locusContextRefBasesNSlack + "\n");
        String varType = "";
        if(!vcc.isIndel()) {
            if(vcc.isSNP()) {
                varType = "SNP"; //vcc.isSNP();
            } else if(vcc.isSymbolic()) {
                varType = "Symbolic";
            } else if(vcc.isMNP()) {
                varType = "MNP";
            }
        } else if(vcc.isSimpleInsertion()) {
            varType = "InDel.Insertion.Simple";
        } else if (vcc.isSimpleDeletion()) {
            varType = "InDel.Deletion.Simple";
        } else {
            varType = "InDel.Complex";
        }
        flowIntensityWriter.printf("varType: " + varType + "\n");
        flowIntensityWriter.printf("alleles: [");
        Boolean isFirst = true;
        for (Allele allele : vcc.getAlleles()) {
            if(!isFirst) flowIntensityWriter.printf(",");
            flowIntensityWriter.printf(allele.getDisplayString());
            isFirst = false;
        }
        flowIntensityWriter.printf("]\n");
        Tuple<Integer> alignIndexAndOffset = new Tuple<Integer>(-1,-1);
        for (Map.Entry<String,AlignmentContext> stratifiedContext : stratifiedContexts.entrySet()) {
            String sampleName = stratifiedContext.getKey();
            AlignmentContext ac = stratifiedContext.getValue();  //w/o hack to make them the same, rawContext actually had the complete set of reads for some reason.
            flowIntensityWriter.printf("  " + sampleName + "\t" + (ac.hasExtendedEventPileup() ? "Extended" : "") + "\n");
            
            if (!this.displayPerAllele) {
                List<GATKSAMRecord> contextReads =  ac.hasExtendedEventPileup() ?
                        ac.getExtendedEventPileup().getReads() :
                        ac.getBasePileup().getReads();
                outputReadSet(contextReads,startPos,endPos,alignIndexAndOffset);
                ac.getBasePileup();
                continue;
            }
            if (ac.hasBasePileup()) {
                List<Pair<Byte,List<GATKSAMRecord>>> basesNReads;
                basesNReads = ac.getBasePileup().getEventBaseWithReadList(locusContextRefBases.getBytes()[0]);
                for(Pair<Byte,List<GATKSAMRecord>> baseNReads : basesNReads) {
                    flowIntensityWriter.printf("Reads for base: " + (char)((byte)baseNReads.getFirst()) + "\n");
                    outputReadSet(baseNReads.getSecond(),startPos,endPos,alignIndexAndOffset);
                }
            } else {            
                List<Pair<String, List<GATKSAMRecord>>> locusSeqNReads;

                //ToDo: locusContextRefBases.getBytes() is different than vcc.getReference().getBases();
                //ToDo: The translation below seems to work.  Couldn't find a util translation function for this
                //ToDo: so reversed engineered from the output.
                byte [] ref = vcc.getReference().getBases();
                if (vcc.isSimpleInsertion()) {
                    ref = new byte[] {'1','D'};
                } else if (vcc.isSimpleDeletion()) {
                    ref = new byte[] {'+',ref[ref.length-1]};
                }

                flowIntensityWriter.println("Total number of reads in pileup:" + ac.getExtendedEventPileup().getReads().size());
                locusSeqNReads = ac.getExtendedEventPileup().getEventStringWithReadList(ref,true);
                for(Pair<String, List<GATKSAMRecord>> alleleNList : locusSeqNReads) {
                    flowIntensityWriter.printf("Reads for locus Seq: '" + alleleNList.getFirst() + "'\n");
                    outputReadSet(alleleNList.getSecond(),startPos,endPos,alignIndexAndOffset);
                }
            }
        }
        return;
    }    
    
    private void outputReadSet(List<GATKSAMRecord> contextReads, long startPos, long endPos,
                          Tuple<Integer> alignIndexAndOffset) {
        if (contextReads != null)
            flowIntensityWriter.println("Size of contextReads " + contextReads.size() );

        // Loop around positions of variant
        int posBuf = 1;
        for(long position = startPos-posBuf; position<endPos+1+posBuf; position++) {
            for(Boolean posStrandToOutput : new Boolean[] {true, false}) {
                flowIntensityWriter.printf("   " + position + (posStrandToOutput?" (+)":" (-)") + " [");
                Boolean isFirst = true;
                if (contextReads==null) {
                    flowIntensityWriter.printf("No reads.]\n");
                    return;
                }
                for(GATKSAMRecord read : contextReads) {
                    IonSAMRecord contextRead = (IonSAMRecord)read;
                    if(posStrandToOutput && contextRead.getReadNegativeStrandFlag()) continue;
                    if(!posStrandToOutput && !contextRead.getReadNegativeStrandFlag()) continue;
                    if(!isFirst) flowIntensityWriter.printf(",");

                    contextRead.flowAlign.getIndexAndOffset(position, alignIndexAndOffset);
                    int alignIndex = alignIndexAndOffset.one;

                    if(alignIndex>-1) {
                        flowIntensityWriter.printf(""+contextRead.flowAlign.alignment.aln[alignIndex]);
                        flowIntensityWriter.print("/"+SamToFlowgramAlignUtil.DNA[contextRead.flowAlign.alignment.flowOrder[alignIndex]]);
                        flowIntensityWriter.printf("/"+contextRead.flowAlign.getSignal(position));
                    }
                    isFirst = false;
                }
                flowIntensityWriter.printf("]\n");
            }
        }
    }
}
