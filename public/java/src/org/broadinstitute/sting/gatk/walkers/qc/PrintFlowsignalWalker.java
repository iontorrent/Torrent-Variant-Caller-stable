package org.broadinstitute.sting.gatk.walkers.qc;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/8/12
 * Time: 2:36 PM
 * To change this template use File | Settings | File Templates.
 */

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.IonSAMRecord;



// Flow align specific libraries
import org.iontorrent.sam2flowgram.util.*;
//import org.iontorrent.sam2fs.io.*;
import org.iontorrent.sam2flowgram.flowalign.*;

import net.sf.samtools.SAMRecord;
//import net.sf.samtools.*;
//import net.sf.samtools.util.*;
//import net.sf.picard.cmdline.*;
//import net.sf.picard.io.IoUtil;

import java.io.*;
//import java.util.*;
//import java.lang.Math;
// End of flow align specific libraries

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;

/**
 * At each locus in the input data set, prints the reference base, genomic location, and
 * all aligning reads in a compact but human-readable form.
 */
public class PrintFlowsignalWalker extends LocusWalker<AlignmentContext, Integer> implements TreeReducible<Integer> {
    @Output
    private PrintStream out;

    public void initialize() {
    }

    public AlignmentContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        long position = context.getPosition();
        List<Integer> signals = new ArrayList();
        int i;

        // The the reads for a particular AlignmentContext
        if(context.hasReads()) {
            List<GATKSAMRecord> reads;
            if (context.hasExtendedEventPileup()) {
                //reads = context.getExtendedEventPileup().getReads();
                return context; // ignore
            } else {
                reads = context.getBasePileup().getReads();
            }
            for(i=0;i<reads.size();i++) {
                IonSAMRecord read = (IonSAMRecord)reads.get(i);
                try {
                    if(null == read.flowAlign) {
                        continue;
                    }
                    Tuple<Integer> t = read.flowAlign.getSignalAndOffset(position);
                    if(null != t && 0 == t.two) {
                        signals.add(t.one);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("Error in read " + read.getSAMString());
                    System.exit(1);
                }
            }
        }

        if(0 < signals.size()) {
            out.printf( "In map: ref = %s, loc = %s, reads = %s%n", 
                    ref.getBaseAsChar(),
                    context.getLocation(),
                    Arrays.deepToString(signalsToArray(signals)) );
        }

        return context;
    }


    public Integer reduceInit() { return 0; }

    public Integer reduce(AlignmentContext context, Integer sum) {
        return sum + 1;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    @Override
    public boolean generateExtendedEvents() {
        return true;
    }

    public void onTraversalDone(Integer sum) {
    }

    private Integer[] signalsToArray(List<Integer> signals) {
        int i;
        Integer s[] = null;
        s = new Integer[signals.size()];
        for(i=0;i<signals.size();i++) {
            s[i] = signals.get(i);
        }
        return s;
    }
}
