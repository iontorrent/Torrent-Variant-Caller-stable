package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileReader;

import org.broadinstitute.sting.gatk.DownsamplingMethod;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.resourcemanagement.ThreadAllocation;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.sam.IonSamRecordFactory;
import org.broadinstitute.sting.utils.recalibration.BaseRecalibration;

import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/12/12
 * Time: 3:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class IonSAMDataSource extends SAMDataSource {
    final private static IonSamRecordFactory factory = new IonSamRecordFactory();
    public IonSamRecordFactory getFactory() { return factory; }

    public IonSAMDataSource(Collection<SAMReaderID> samFiles, ThreadAllocation threadAllocation, Integer numFileHandles, GenomeLocParser genomeLocParser) {
        super(samFiles,threadAllocation,numFileHandles,genomeLocParser);
    }

    public IonSAMDataSource(
            Collection<SAMReaderID> samFiles,
            ThreadAllocation threadAllocation,
            Integer numFileHandles,
            GenomeLocParser genomeLocParser,
            boolean useOriginalBaseQualities,
            SAMFileReader.ValidationStringency strictness,
            Integer readBufferSize,
            DownsamplingMethod downsamplingMethod,
            ValidationExclusion exclusionList,
            Collection<ReadFilter> supplementalFilters,
            boolean includeReadsWithDeletionAtLoci,
            boolean generateExtendedEvents) {
        super(samFiles,threadAllocation,numFileHandles,genomeLocParser,useOriginalBaseQualities,strictness,readBufferSize,
                downsamplingMethod,exclusionList,supplementalFilters,includeReadsWithDeletionAtLoci,generateExtendedEvents);
    }

    public IonSAMDataSource(
            Collection<SAMReaderID> samFiles,
            ThreadAllocation threadAllocation,
            Integer numFileHandles,
            GenomeLocParser genomeLocParser,
            boolean useOriginalBaseQualities,
            SAMFileReader.ValidationStringency strictness,
            Integer readBufferSize,
            DownsamplingMethod downsamplingMethod,
            ValidationExclusion exclusionList,
            Collection<ReadFilter> supplementalFilters,
            boolean includeReadsWithDeletionAtLoci,
            boolean generateExtendedEvents,
            BAQ.CalculationMode cmode,
            BAQ.QualityMode qmode,
            IndexedFastaSequenceFile refReader,
            BaseRecalibration bqsrApplier,
            byte defaultBaseQualities) {
        super(samFiles,threadAllocation,numFileHandles,genomeLocParser,useOriginalBaseQualities,strictness,readBufferSize,
                downsamplingMethod,exclusionList,supplementalFilters,includeReadsWithDeletionAtLoci,generateExtendedEvents,
                cmode, qmode, refReader, bqsrApplier, defaultBaseQualities);
        System.out.printf("In IonSAMDataSource constructor.\n");
        factory.setReferenceSequence(refReader);

    }
}
