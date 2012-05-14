package org.broadinstitute.sting.utils.sam;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordFactory;
import net.sf.samtools.BAMRecord;

import org.iontorrent.sam2flowgram.util.FlowAlignRecord;
import org.iontorrent.sam2flowgram.flowalign.FlowOrder;
import org.iontorrent.sam2flowgram.util.ReferenceSequence;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/8/12
 * Time: 4:34 PM
 * To change this template use File | Settings | File Templates.
 */

public class IonSamRecordFactory extends GATKSamRecordFactory {

    /** Ion Torrent specific attributes */
    // Key and flow order sequences
    String keySeq = "TCAG";
    String flowOrderStr = "TACGTACGTCTGAGCATCGATCGATGTACAGC";
    FlowOrder flowOrder = new FlowOrder(flowOrderStr,keySeq);
    
    // Flow alignment parameters
    int phase_penalty = 1;
    int refBaseOffset = 10;
    boolean validateAlignments = false;

    // Reference sequence
    // TODO, need to unify this with GATK's method for handling the reference
    ReferenceSequence referenceSequence = null;
    int lastContigIdx = -1;  //See if moveTo needs to be called.
    //TODO, threading ramifications of lastContigIdx?
    
    // Initialization of parameters
    public void setIonSequences(String keySeq, String flowOrderStr) {
        if(keySeq!=null) this.keySeq = keySeq;
        if(flowOrderStr!=null) this.flowOrderStr = flowOrderStr;
        this.flowOrder = new FlowOrder(flowOrderStr, keySeq);    
    }
    
    public void setFlowAlignSettings(Integer phase_penalty, Integer refBaseOffset) {
        if(phase_penalty!=null) this.phase_penalty=phase_penalty;
        if(refBaseOffset!=null) this.refBaseOffset=refBaseOffset;
    }
    
    public void setReferenceSequence(IndexedFastaSequenceFile referenceSequenceFile) {
        File refFile = new File(referenceSequenceFile.toString());
        System.out.printf("Using reference file '" + referenceSequenceFile.toString() + "'\n");
        try {
            referenceSequence = new ReferenceSequence(refFile,false);
            //referenceSequence.moveTo(lastContigIdx);
        }  catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error when trying to read in reference from " + referenceSequenceFile);
            System.err.println("Please report bugs to eric.tsung@lifetech.com");
            System.exit(1);
        }            
    }
    // Check for change in reference contig
    private void checkChangeInContig (SAMRecord newSAMRecord) {
        // Check if contig is different then that last one
        int curContigIdx = newSAMRecord.getReferenceIndex();
        if(lastContigIdx!=curContigIdx) {
            try {
                referenceSequence.moveTo(curContigIdx);
            } catch (Exception e) {
                e.printStackTrace();
                System.err.println("Error when moving to contig index " + curContigIdx + " that has a name of '" +
                        referenceSequence.getDictionary().toString() + "'.\n");
            }
            lastContigIdx = curContigIdx;
        }
    }
    
    //Perform flow alignment
    public void doFlowAlign(IonSAMRecord read) {
        try {
            read.flowAlign = new FlowAlignRecord((SAMRecord)read,0,flowOrder);
            read.flowAlign.setAlignment(referenceSequence,phase_penalty,refBaseOffset,validateAlignments);
        } catch (Exception e) {
            //System.err.println("Error in read");
            e.printStackTrace();
            System.err.println("Error in read " + read.getSAMString());
            //System.err.println("Please report bugs to eric.tsung@lifetech.com");
            //System.exit(1);
        }
    }
    
    /** Create a new BAM Record. */
    public BAMRecord createBAMRecord(final SAMFileHeader header,
                                     final int referenceSequenceIndex,
                                     final int alignmentStart,
                                     final short readNameLength,
                                     final short mappingQuality,
                                     final int indexingBin,
                                     final int cigarLen,
                                     final int flags,
                                     final int readLen,
                                     final int mateReferenceSequenceIndex,
                                     final int mateAlignmentStart,
                                     final int insertSize,
                                     final byte[] variableLengthBlock) {
        //System.out.print("Making new SAMRecord object of IonSamRecordFactory!\n");
        
        IonSAMRecord newSAMRecord =
                new IonSAMRecord(header,
                        referenceSequenceIndex,
                        alignmentStart,
                        readNameLength,
                        mappingQuality,
                        indexingBin,
                        cigarLen,
                        flags,
                        readLen,
                        mateReferenceSequenceIndex,
                        mateAlignmentStart,
                        insertSize,
                        variableLengthBlock);
        checkChangeInContig(newSAMRecord);
        doFlowAlign(newSAMRecord);
        return newSAMRecord;
    }
}
