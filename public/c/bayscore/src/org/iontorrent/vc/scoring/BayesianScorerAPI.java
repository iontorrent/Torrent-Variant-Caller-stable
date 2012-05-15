package org.iontorrent.vc.scoring;

/**
 * Created by IntelliJ IDEA.
 * User: ionadmin
 * Date: 3/26/12
 * Time: 11:09 AM
 * To change this template use File | Settings | File Templates.
 */
public class BayesianScorerAPI {

    // Interface to the C API
    static {
        //System.setProperty("java.library.path","/home/ionadmin/src/TorrentVariantCaller/source/bayscore/fromZhengDir");
       System.loadLibrary("BayesianScorerAPI");
    }

    // To make it thread safe, pid is returned from addRef, and needs to be used
    // for all calls, until finished is called, which releases it.
    
    // init specifying number of thread spots to reserve 
    // (basically max. number of simultenous calls)
    public native void rescorer_init(int numThread);

    // Calculates the score.  Done on a per variant basis.
    public native double calScore(int pid);

    public native void rescorer_end();
    //public native void addRef(int numRefBases, char* refBases);

    // Add reference window around variant, returns pid for future calls
    // Called first for each locus, a is the total length of the window
    public native int addRef(int a, byte [] refbase);
    //public native void addVariant(int varLocInRef, char *predicted_varaint, char *ref_allele);

    // Add variants found at the working locus
    // v is relative position of variant in the ref window given
    // pv, is variant sequence or representation, i.e. G / +T / 2D for SNP, insertion, deletion
    // ref, if used, need pv to coorespond, ie. ref=A, pv=AT, instead of +T (not used for now)
    public native void addVariant(int pid, int v, byte [] pv, byte [] ref);
    public native void addVariant(int pid, int v, byte [] pv);

    //Add reads for the variant or reference.  pid=0 is for the reference
    //Meaning of arguments:
    // int pid, int varNum, int varFlowIndex, int numFlow, char *flowbases, int *flowIntensities
    public native void addRead(int pid, int v, int v2, int f, byte [] fl, int [] fs, byte [] al, int ispositive);

    // Called when done with this locus
    public native void finished(int pid);

    public static void main(String[] args) {
    BayesianScorerAPI api = new BayesianScorerAPI();
        api.rescorer_init(1);
        byte [] a = new byte[10];
        a[0] = 'A';
        a[1] = 'C';
        a[2] = 'G';
        a[3] = 'T';
        a[4] = 0;
        int pid = api.addRef(5, a);
        int [] b = new int[5];
        b[0]=  2;
        b[1] = 1;
        b[2] = 2;
        b[3] = 3;
	byte [] align = new byte[10];
	align[0] = '|';
	align[1] = '|';
	align[2] = '+';
	align[3] = '-';
        api.addRead(pid, 1,2,4,a,b, align, 1);
	String var = "+G";
	var += '\000';
	api.addVariant(pid, 2,var.getBytes());
    }
}
