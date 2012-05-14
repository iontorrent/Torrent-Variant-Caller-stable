
package org.broadinstitute.sting.gatk.walkers.torrent;

/**
 * Created by IntelliJ IDEA.
 * User: brinzadc
 * Date: 3/15/12
 * Time: 8:48 AM
 * To change this template use File | Settings | File Templates.
 */

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reference.ReferenceDataSource;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.Hashtable;
import java.util.Vector;


public class IndelAssembly extends ReadWalker<Integer, Integer> {

    @Output
    public PrintStream out;

    private int WINDOW_SIZE = 1000;
    private int WINDOW_PREFIX = 300;
    private int MIN_NUM_SOFTREADS = 10;
    private double MIN_RATIO_SOFT = 0.1;
    private int MIN_ASSEMBLY_WINDOW = 20;              // max indel length is MIN_ASSEMBLY_WINDOW - ANCHOR_LEN 
    private int KMER_LEN = 11;                         // fixed for now
    private int SHORT_SUFFIX_MATCH = 4;
    private int [] softwindow = new int[WINDOW_SIZE];  // implement as list later
    private int [] covwindow = new int[WINDOW_SIZE];   // implement as recursive list
    private int curLeft = 0;
    private String curChrom = "";
    private int softStart = 0;
    private int softLen = 0;
    private int softCov = 0;
    private int preSoftCov = 0;
    private Vector<GATKSAMRecord> ReadsBuffer = new Vector<GATKSAMRecord>();
    private IndexedFastaSequenceFile ancestralAlignments = null;


    @Override
    public Integer map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {

             if(curChrom.compareTo(read.getReferenceName())!=0)
        {
            for(int i=0;i<WINDOW_SIZE;i++)
            {
                DetectCandidateRegions(i);
            }

            FlushBuffer(Integer.MAX_VALUE);

            curChrom = read.getReferenceName();
            curLeft = read.getSoftStart()>WINDOW_PREFIX?(read.getSoftStart() - WINDOW_PREFIX):0;
            softStart = curLeft;
            softLen = 0;
            preSoftCov = 0;
            softCov = 0;
            for(int i=0;i<WINDOW_SIZE;i++) {softwindow[i]=0; covwindow[i]=0;}
        }

        if(read.getSoftEnd() >= curLeft + WINDOW_SIZE)
        {
            int delta = read.getSoftEnd() - curLeft - WINDOW_SIZE + 1;

            for(int i=0;i<(delta>WINDOW_SIZE?WINDOW_SIZE:delta);i++)
            {
                DetectCandidateRegions(i);
            }

            if(delta>=WINDOW_SIZE)   {
                curLeft = read.getSoftStart()>WINDOW_PREFIX?(read.getSoftStart() - WINDOW_PREFIX):0;
                for(int i=0;i<WINDOW_SIZE;i++) {softwindow[i]=0; covwindow[i]=0;}
            }
            else
            {
                for(int i=0;i<WINDOW_SIZE-delta;i++) {softwindow[i]=softwindow[i+delta]; covwindow[i]=covwindow[i+delta];}
                for(int i=WINDOW_SIZE-delta;i<WINDOW_SIZE;i++) {softwindow[i]=0; covwindow[i]=0;}
                curLeft += delta;
            }

            FlushBuffer(curLeft);
        }

        ReadsBuffer.add(read);
        //out.println(curChrom + " " + Integer.toString(curLeft) + " " + Integer.toString(read.getSoftStart()) + "-" + Integer.toString(read.getSoftEnd()) + "  " + read.getCigarString());
        AddSoftClipCounts(read.getSoftStart() - curLeft, read.getCigar());
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Integer reduceInit() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    private void DetectCandidateRegions( int i)
    {
            if(softwindow[i] > MIN_NUM_SOFTREADS && (covwindow[i]==0 || softwindow[i] > covwindow[i]*MIN_RATIO_SOFT))
            {
                softLen++;
                softCov+=softwindow[i];
            }
            else
            {
                if(softLen>MIN_ASSEMBLY_WINDOW)
                {
                    if(softCov/softLen > MIN_RATIO_SOFT*Math.max(preSoftCov, covwindow[i])) SegmentAssembly();
                }
                softStart = curLeft + i;
                softLen = 0;
                softCov = 0;
                preSoftCov = covwindow[i];
            }

    }

    private class TKmer{
        int Freq;
        int PosInReference;
     public TKmer()
     {
         Freq = 0;
         PosInReference = -1;
     }
    }

    private void BuildKMerSpectrum(Hashtable<String, TKmer> spectrum)
    {
        for(int i=0;i<ReadsBuffer.size();i++)
        {
            final GATKSAMRecord read = ReadsBuffer.get(i);
            int startOffset = read.getSoftStart() - softStart;
            int readOffset = 0;
            Cigar cigar = read.getCigar();
                                                              
            for(int j=0;j<cigar.numCigarElements();j++)
            {
                CigarOperator cgo = cigar.getCigarElement(j).getOperator();
                int cgl = cigar.getCigarElement(j).getLength();

                if(cgo == CigarOperator.SOFT_CLIP)
                {
                    if(startOffset + readOffset + cgl > 0 && startOffset + readOffset < softLen)
                    {
                      int startPos = (startOffset + readOffset >= 0)?(readOffset - KMER_LEN):(-startOffset -  KMER_LEN);
                      int  stopPos = readOffset + cgl + KMER_LEN + 1;
                      final String softAln = read.getReadString();  
                      final String novelSeq = softAln.substring(startPos<0?0:startPos, stopPos>softAln.length()?softAln.length():stopPos);
                        
                     for(int x = 0; x <= novelSeq.length() - KMER_LEN; x++)
                      {
                         final String kmer = novelSeq.substring(x,x+KMER_LEN);
                         if(spectrum.containsKey(kmer)) spectrum.get(kmer).Freq++;
                         else spectrum.put(novelSeq.substring(x,x+KMER_LEN),new TKmer());
                      }
                    }
                           readOffset+=cgl;
                }
                else if(cgo != CigarOperator.DELETION && cgo != CigarOperator.HARD_CLIP) readOffset+=cgl;
            }
        
        }
    }

    private void PrintVCF(String contig, int pos, String ref, String var, 
                          int varLen, int type, int qual, int totCounts, int refCounts, int varCounts)
    {
        if(qual==30 && type%2==1) qual=5;
        else
        if(varCounts<50) qual=7;
        else
        if(qual==30 && varLen<15) qual=9;

        if(totCounts<=0) totCounts = 1;
        if(refCounts<=0) refCounts = 0;
        float varfq = totCounts<=varCounts? (float) 1.0 :((float)varCounts)/((float)totCounts);
        float reffq = totCounts<=refCounts? (float) 1.0 :((float)refCounts)/((float)totCounts);
        
        String genotype = reffq > 0.2 ? "0/1" : "1/1";
        
        out.println(contig                + "\t" +
                    Integer.toString(pos) + "\t" +
                    "."                   + "\t" +
                    ref                   + "\t" +
                    var                   + "\t" +
                    Integer.toString(qual*100)+ "\t" +
                    "."                   + "\t" +
                    "Bayesian_Score=" + Integer.toString(qual)+";"+
                    "Variant-freq=" + Float.toString(varfq)+   ";"+
                    "Num-spanning-ref-reads=" + Integer.toString(refCounts)+ ";"+
                    "Num-variant-reads=" + Integer.toString(varCounts)+ ";"+
                    "Variant-length=" + Integer.toString(varLen) + "\t" +
                    "GT:AD:DP:FA:GQ:MQ0:PL\t"+
                    genotype + ":" + Integer.toString(refCounts) + "," + Integer.toString(varCounts) + ":" +
                    Integer.toString(totCounts) + ":" + Float.toString(varfq) + ":99:0:" + Integer.toString(qual*100));
    }

    private void getPath(String anchorKMer, int anchorFreq, Hashtable<String, TKmer> spectrum, int genStart, final String reference)
    {
      double cutFreq = 0.2*softCov/softLen;
      int startPos = spectrum.get(anchorKMer).PosInReference + KMER_LEN;
      String varSeq = "";
      String nextKmer = anchorKMer;
      int firstKMerFreq = 0;
      int lastKMerFreq = 0;

      while(varSeq.length()<WINDOW_PREFIX)
      {
      nextKmer = nextKmer.substring(1);

        String tmp = nextKmer + "A";
        int a = spectrum.containsKey(tmp)?spectrum.get(tmp).Freq:0;

        tmp = nextKmer + "C";
        int c = spectrum.containsKey(tmp)?spectrum.get(tmp).Freq:0;

        tmp = nextKmer + "G";
        int g = spectrum.containsKey(tmp)?spectrum.get(tmp).Freq:0;

        tmp = nextKmer + "T";
        int t = spectrum.containsKey(tmp)?spectrum.get(tmp).Freq:0;

        int maxFreq = Math.max(a,Math.max(c,Math.max(g,t)));

        if(maxFreq>cutFreq)
        {
            int nCandPath = 0;
            if(2*a >maxFreq) nCandPath++;
            if(2*c >maxFreq) nCandPath++;
            if(2*g >maxFreq) nCandPath++;
            if(2*t >maxFreq) nCandPath++;
            
            if(nCandPath == 1)
            {
                if(varSeq.length()==0) { firstKMerFreq = maxFreq; lastKMerFreq = maxFreq; }
                if(a==maxFreq) { varSeq+="A"; nextKmer+="A";}
                else
                if(c==maxFreq) { varSeq+="C"; nextKmer+="C";}
                else
                if(g==maxFreq) { varSeq+="G"; nextKmer+="G";}
                else
                if(t==maxFreq) { varSeq+="T"; nextKmer+="T";}

              int stopPos = spectrum.get(nextKmer).PosInReference;
              if( stopPos >= startPos )
              {
                  int anchor2Freq = spectrum.get(nextKmer).Freq;
                  int refCov = 0;
                  if((startPos<=KMER_LEN) && (stopPos>=softLen+KMER_LEN)) refCov = (anchorFreq + anchor2Freq)/2;
                  if(startPos>KMER_LEN) refCov = anchorFreq;
                  if(stopPos<softLen+KMER_LEN) refCov += anchor2Freq;

               if(varSeq.length() - KMER_LEN > stopPos - startPos)
                   // INSERTION  type 0
                    PrintVCF(curChrom, genStart + startPos - 1,  anchorKMer.substring(anchorKMer.length()-1), anchorKMer.substring(anchorKMer.length()-1)+varSeq.substring(0,varSeq.length() - KMER_LEN),varSeq.length() - KMER_LEN,0,100,refCov, refCov-Math.max(firstKMerFreq,lastKMerFreq)/2, Math.max(firstKMerFreq,lastKMerFreq)/2);
                  else
                  if(varSeq.length() - KMER_LEN <= stopPos - startPos)
                   // DELETION type 1
                   PrintVCF(curChrom, genStart + startPos - 1,  reference.substring(startPos - 1, stopPos), reference.substring(startPos - 1, startPos),stopPos - startPos,1,100,refCov, refCov-Math.max(firstKMerFreq,lastKMerFreq)/2, Math.max(firstKMerFreq,lastKMerFreq)/2);
              return;
              }
                lastKMerFreq = maxFreq;
            }
            else
             {
                 //out.println("MULTI " + Integer.toString(startPos) +  " " + varSeq);
                 return;
             }
        }  else
        {
            if(varSeq.length() >= SHORT_SUFFIX_MATCH)
            {
            String shortSuffix = varSeq.substring(varSeq.length()-SHORT_SUFFIX_MATCH);
            String shortRef = reference.substring(startPos,startPos + KMER_LEN>reference.length()?reference.length():startPos + KMER_LEN);
            if(shortRef.contains(shortSuffix))
            {
              int delta = 0; 
              for(int i=0; i<shortRef.length() - SHORT_SUFFIX_MATCH;i++)
              {
                if(shortRef.regionMatches(i,shortSuffix,0,SHORT_SUFFIX_MATCH)) {delta = i; break;}
              }
              if(varSeq.length() < delta + SHORT_SUFFIX_MATCH)
              {
                  // DELETION type 3
                  PrintVCF(curChrom, genStart + startPos - 1,  reference.substring(startPos - 1, startPos + delta + SHORT_SUFFIX_MATCH - varSeq.length()),reference.substring(startPos - 1, startPos),delta + SHORT_SUFFIX_MATCH - varSeq.length(),3,30,anchorFreq,anchorFreq-Math.max(firstKMerFreq,lastKMerFreq)/2, Math.max(firstKMerFreq,lastKMerFreq)/2);
              }
                else
              {
                  // INSERTION  type 2
                  PrintVCF(curChrom, genStart + startPos - 1,  anchorKMer.substring(anchorKMer.length()-1), anchorKMer.substring(anchorKMer.length()-1)+varSeq.substring(0,varSeq.length() - delta - SHORT_SUFFIX_MATCH),varSeq.length() - delta - SHORT_SUFFIX_MATCH,2,30,anchorFreq, anchorFreq-Math.max(firstKMerFreq,lastKMerFreq)/2, Math.max(firstKMerFreq,lastKMerFreq)/2);
              }
                return;
             }
            }
           // out.println("INCOM " + Integer.toString(startPos) +  " " + varSeq);
            return;
        }
      }
    }
    
    private void DetectIndel(int genStart, int genStop, final String reference, Hashtable<String, TKmer> spectrum)
    {
       String anchorKMer = null;
       int numHits = 0; 
       int firstAnchor = -1;
        int anchorFreq = 0;
       for(int i=0;i<reference.length()-KMER_LEN+1;i++)
       {
         final String refKMer= reference.substring(i,i + KMER_LEN);
         if(spectrum.containsKey(refKMer)) 
         {
             final TKmer kmer = spectrum.get(refKMer);
             //if(kmer.PosInReference!=-1) out.println("WARNING: Assembly is circularized!!!");  else
             {
                  kmer.PosInReference = i;
                  numHits++;
             }   
             if(kmer.Freq > 0.2*softCov/softLen && (firstAnchor == -1 || i-firstAnchor == 1))
             {
                 firstAnchor = i;
                 anchorKMer = refKMer;
                 anchorFreq = kmer.Freq;
             }    
         }
       }

       // out.println(" most left anchor: " + Integer.toString(firstAnchor) + " " + anchorKMer + "  # hits " + Integer.toString(numHits));
        if(anchorKMer!=null) getPath(anchorKMer, anchorFreq, spectrum, genStart, reference);

    }

    private void SegmentAssembly()
    {
        if(ancestralAlignments == null) ancestralAlignments = new ReferenceDataSource(getToolkit().getArguments().referenceFile).getReference();
        int genStart = softStart>KMER_LEN?softStart-KMER_LEN:0;
        int endOfChr = ancestralAlignments.getSequenceDictionary().getSequence(curChrom).getSequenceLength();
        int genStop = softStart + softLen + KMER_LEN + 1;
            if(genStop>endOfChr) genStop = endOfChr;

        String reference = new String(ancestralAlignments.getSubsequenceAt(curChrom, genStart, genStop).getBases());
       // out.println(curChrom + ":" + Integer.toString(softStart) + " " + Integer.toString(softLen) + " " + ReadsBuffer.size() + " " + Double.toString(softCov/softLen) + " " + reference);
        Hashtable spectrum = new Hashtable<String, TKmer>();
        
        BuildKMerSpectrum(spectrum);
        DetectIndel(genStart, genStop, reference, spectrum);
        
       // Vector<GATKSAMRecord> AssemblyBuffer = new Vector<GATKSAMRecord>();
    }

    private void FlushBuffer(int LeftMostPos)
    {
      while(ReadsBuffer.size()>0 && ReadsBuffer.get(0).getSoftEnd()<LeftMostPos){
          // out.println(ReadsBuffer.get(0).getSAMString());
          ReadsBuffer.remove(0);
      }
    }

    private void AddSoftClipCounts(int offset, Cigar cigar)
    {
      if(offset<0)
      {
        //  out.println("A read started upfront window! Increase WINDOW_PREFIX. Read skipped.");
          return;
      }
      for(int i=0;i<cigar.numCigarElements();i++)
      {
          CigarOperator cgo = cigar.getCigarElement(i).getOperator();
          int cgl = cigar.getCigarElement(i).getLength();

          if(cgo == CigarOperator.SOFT_CLIP)
          {
           for(int j=cgl+offset-1;j>=offset;j--) softwindow[j]++;
           offset+=cgl;
          }
          else if(cgo != CigarOperator.INSERTION && cgo != CigarOperator.HARD_CLIP)
          {
          for(int j=cgl+offset-1;j>=offset;j--) covwindow[j]++;
          offset+=cgl;         
          }
              
     }

    }
}



