package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.iontorrent.sam2flowgram.util.FlowAlignRecord;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/8/12
 * Time: 2:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class IonSAMRecord extends GATKSAMRecord {

    // Contains the flow alignment record
    public FlowAlignRecord flowAlign = null;

    // Defining constructors, picard SAMRecord doesn't have default, so have to define it
    // in all subclasses.
    public IonSAMRecord(final SAMFileHeader header) {
        super(header);
    }
    public IonSAMRecord(final SAMRecord read) {
        super(read);
    }

    public IonSAMRecord(final SAMFileHeader header,
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
        super(header, referenceSequenceIndex, alignmentStart, readNameLength, mappingQuality, indexingBin, cigarLen,
                flags, readLen, mateReferenceSequenceIndex, mateAlignmentStart, insertSize, variableLengthBlock);
    }

    // Yuck, basically a copy of the superclass, except there's a new IonSAMRecord instead.
    public static GATKSAMRecord emptyRead(IonSAMRecord read) {
        GATKSAMRecord emptyRead = new IonSAMRecord(read.getHeader(),
                read.getReferenceIndex(),
                0,
                (short) 0,
                (short) 0,
                0,
                0,
                read.getFlags(),
                0,
                read.getMateReferenceIndex(),
                read.getMateAlignmentStart(),
                read.getInferredInsertSize(),
                null);
        System.out.print("In GATKSAMRecord\n");
        emptyRead.setCigarString("");
        emptyRead.setReadBases(new byte[0]);
        emptyRead.setBaseQualities(new byte[0]);

        SAMReadGroupRecord samRG = read.getReadGroup();
        emptyRead.clearAttributes();
        if (samRG != null) {
            GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord(samRG);
            emptyRead.setReadGroup(rg);
        }

        return emptyRead;
    }

}
