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
import java.util.Random;

/**
 * At each locus in the input data set, prints the reference base, genomic location, and
 * all aligning reads in a compact but human-readable form.
 */
public class PrintFlowgramAlignWalker extends LocusWalker<AlignmentContext, Integer> implements TreeReducible<Integer> {
    @Output
    private PrintStream out;
    String flowOrderStr = "TACGTACGTCTGAGCATCGATCGATGTACAGC";
    String keySeq = "TCAG";
    FlowOrder flowOrder = new FlowOrder(flowOrderStr,keySeq);
    int phase_penalty = 1;   // The phase penalty when converting an alignment to flow space.
    int refBaseOffset = 10;  // Extra reference bases to consider during alignment.
    ReferenceSequence referenceSequence = null;
    Random generator = new Random();
    int curContigIdx = 0;

    int numReadsIterated = 0;
    int numUniqReadsFlowAligned = 0;

    public void initialize() {
        String refFN = "/ldata/ion_variant_calling/referenceLibrary/tmap-f2/hg19/hg19.fasta";
        File refFile = new File(refFN);
        out.printf( "Using reference file '" + refFN + "'\n");
        try {
            //referenceSequence = new ReferenceSequence(refFile);
            //referenceSequence.moveTo(curContigIdx);
        }  catch (Exception e) {
            e.printStackTrace();
            System.err.println("Please report bugs to eric.tsung@lifetech.com");
            System.exit(1);
        }
    }

    public AlignmentContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String eventString = "";

        //Check to see if moved to different contig/
        int lastContigIdx = curContigIdx;
        curContigIdx = ref.getLocus().getContigIndex();
        if(lastContigIdx!=curContigIdx) {
            try {
                //referenceSequence.moveTo(curContigIdx);
            } catch (Exception e) {
                e.printStackTrace();
                System.err.println("Error when moving to contig index " + curContigIdx + " that has a name of '" +
                        referenceSequence.getDictionary().toString() + "'.\n");
            }
        }

        boolean validateAlignments = false;

        // The the reads for a particular AlignmentContext
        if(context.hasReads()) {
            List<GATKSAMRecord> reads;
            if (context.hasExtendedEventPileup()) {
                reads = context.getExtendedEventPileup().getReads();
            } else {
                reads = context.getBasePileup().getReads();
            }
            for(int i = 0; i < reads.size(); i++) {
                IonSAMRecord read = (IonSAMRecord)reads.get(i);
                try {
                    numReadsIterated++;
                    if(read.flowAlign == null) {
                        numUniqReadsFlowAligned++;
                        //read.flowAlign = new FlowAlignRecord((SAMRecord)read,0,flowOrder);
                        //read.flowAlign.setAlignment(referenceSequence,phase_penalty,refBaseOffset,validateAlignments);

                    }
                    if(generator.nextInt(10000)<5) { // Randomly select some of alignments for printing
                        // P = 5/10000 * averageReadLength
                        FlowAlignRecord flowAlign = read.flowAlign;
                        System.out.print( "Total reads: " + numReadsIterated + ", uniq reads: " + numUniqReadsFlowAligned + "\n");
                        System.out.print( flowAlign.getSAMString() );
                        System.out.println(flowAlign.alignment.getAlignmentString());
                        System.out.println(flowAlign.positionStart+","+flowAlign.positionEnd+","+flowAlign.alignment.getScore());

                    }
                } catch (Exception e) {
                    //System.err.println("Error in read");
                    e.printStackTrace();
                    System.err.println("Error in read " + read.getSAMString());
                    //System.err.println("Please report bugs to eric.tsung@lifetech.com");
                    //System.exit(1);
                }
            }
        }

        if (context.hasExtendedEventPileup()) {
            context.getExtendedEventPileup().getReads();
            List<Pair<String,Integer>> events= context.getExtendedEventPileup().getEventStringsWithCounts();
            for (int i = 0; i < events.size(); i++) {
                eventString += "{" + events.get(i).getFirst() + ", " + events.get(i).getSecond() + "} ";
            }
        }
//        out.printf( "In map: ref = %s, loc = %s %s %s, reads = %s%n", ref.getBaseAsChar(),
//                                                                   context.getLocation(),
//                                                                   context.hasExtendedEventPileup() ? "[extended!]" : "",
//                                                                   eventString,
//                                                                   Arrays.deepToString( getReadNames(context.getReads()) ) );

//        if(context.hasExtendedEventPileup()) {
//            ReadBackedExtendedEventPileup pileup = context.getExtendedEventPileup();
//            pileup.getEventStringsWithCounts();
//        } else {
//            ReadBackedPileup pileup = context.getBasePileup();
//        }

        return context;
    }



    public Integer reduceInit() { return 0; }

    public Integer reduce(AlignmentContext context, Integer sum) {
        return sum + 1;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    private String[] getReadNames( List<GATKSAMRecord> reads ) {
        String[] readNames = new String[ reads.size() ];
        for( int i = 0; i < reads.size(); i++ ) {
            readNames[i] = String.format("%nname = %s, start = %d, end = %d", reads.get(i).getReadName(), reads.get(i).getAlignmentStart(), reads.get(i).getAlignmentEnd());
        }
        //Arrays.sort(readNames);
        return readNames;
    }

    @Override
    public boolean generateExtendedEvents() {
        return true;
    }

     public void onTraversalDone(Integer sum) {
         System.out.print( "Total reads: " + numReadsIterated + ", uniq reads: " + numUniqReadsFlowAligned + "\n");
     }
}
