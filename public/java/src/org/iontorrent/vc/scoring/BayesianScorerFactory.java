package org.iontorrent.vc.scoring;

import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/27/12
 * Time: 5:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class BayesianScorerFactory {

    /*
    The actual API call class
     */
    protected BayesianScorerAPI bayAPI;

    /*
    Print stream for the reads' flowgram alignments in the context of the variant.
    */
    protected PrintStream flowAlignContextPileupStrm;
    protected PrintStream debugStrm;
    
    Boolean doFlowAlignContextOutput = true;
    // For debugging can turn this off
    Boolean doAPICall = true;

    // Number of threads
    private Integer numThreads;

    public BayesianScorerFactory(PrintStream flowAlignContextPileupStrm, Integer numThreads) {
        this.flowAlignContextPileupStrm = flowAlignContextPileupStrm;
        this.numThreads = numThreads;
        if (numThreads==null) {
            doAPICall = false;
        }
        init();
    }

    public void addDebugStream(PrintStream debugStrm) {
        this.debugStrm = debugStrm;
    }

    // init just once to load up .so
    public void init() {
        bayAPI = new BayesianScorerAPI();
        if(doFlowAlignContextOutput)
            flowAlignContextPileupStrm.printf("rescorer_init:" + numThreads + "\n");
        if(doAPICall) bayAPI.rescorer_init(numThreads);
    }

    public BayesianScorer createBayesianScorer() {
        BayesianScorer bayesianScorer = new BayesianScorer(bayAPI, flowAlignContextPileupStrm, doAPICall);
        return bayesianScorer;
    }

    public BayesianScorer createDebugScorer() {
        BayesianScorer debugScorer = new BayesianScorer(bayAPI, debugStrm, false);
        return debugScorer;
    }

}
