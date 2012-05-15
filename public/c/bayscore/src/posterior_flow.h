/*
 * Copyright (c) 2011 Life Technologies Corporation. All rights reserved.
 */

#ifndef POSTERIOR_PROB_H_
#define POSTERIOR_PROB_H_

#define MAX_LEN 80
#define NUM_HHH 500000
#define MINI_PEN 0.6
#define DEL_PENALTY 30.0
#define CONTEXT 3
#define LF 1

//penalty  is -10*log10(pro)

class flow_list 
{
    public:
	flow_list() {init(NUM_HHH);reset();}
	flow_list(int n) {
		init(n);
		reset();
	}
	~flow_list() {
	    delete [] seq[0];
	    delete [] seq;
	    delete [] flow[0];
	    delete [] flow;
	    delete [] prior;
	    delete [] length;
   	}
	void reset() {num_h = 0; v_type = 0;}
	void add_list(unsigned char *s, int *f, int len) {
	    if (num_h >= NUM_H) {return; } //
	    length[num_h] = len;
	    int i;
	    for (i = 0; i < len; i++) {seq[num_h][i] = s[i]; flow[num_h][i] = f[i];}
	    num_h++;
	}
	void add_list(unsigned char *s, int *f, int len, double p) {
	    if (num_h >= NUM_H) {return; }
	    prior[num_h] = p;
	    add_list(s, f, len);
	}
	int num_list() { return num_h;}
	int get_list(int n, unsigned char *&s, int *& f, int &le, double &p) {
	    if (n >= num_h) return 0;
	    s = seq[n]; le = length[n]; p = prior[n]; f = flow[n];
	    return 1;
	}
        int get_list(int n, unsigned char *&s, int *& f, int &le) {
            if (n >= num_h) return 0;
            s = seq[n]; le = length[n]; f = flow[n];
            return 1;
        }
	void output_hyp(int nh) {
	    int i;
	    char code[] = "ACGT";
	    for (i = 0; i < length[nh]; i++) {
		int j;
		for (j = 0; j < flow[nh][i]; j++) {
		    printf("%c", code[seq[nh][i]]);
		}
	    }
	    printf("\n");
	}
	void set_is_indel() { v_type = 1; }
	int is_indel() { return (v_type==1);}
    protected:
	void init(int n) {
	    NUM_H = n;
	    seq = new unsigned char*[NUM_H];
	    seq[0] = new unsigned char[NUM_H*MAX_LEN];
	    flow = new int*[NUM_H];
	    flow[0] = new int[NUM_H*MAX_LEN];
	    length = new int[NUM_H];
	    prior = new double[NUM_H];
	    int i;
	    for (i = 1; i < NUM_H; i++) {
		seq[i] = seq[i-1]+MAX_LEN;
		flow[i] = flow[i-1]+MAX_LEN;
	    }
	}
	int num_h;
	int NUM_H;
	int v_type;
	unsigned char **seq;
	int **flow;
	int *length;
	double *prior; // or strand
};

class flow_poster_prob_calc
{
    public:
	flow_poster_prob_calc() {
	    unsetRef();
	    do_split = 0;
	}
	~flow_poster_prob_calc(){ ;}
	void set_dosplit(int i) {do_split = i;}
	double calc_post_prob(int num, unsigned char *flow_order, int *flow_sig) {
	    return calc_post_prob(num, flow_order, flow_sig, 1);
	}
	double calc_post_prob(int num, unsigned char *flow_order, int *flow_sig, int dir);
	double calc_post_prob(int num_ref, unsigned char *ref_seq, int *ref_mer, int num, unsigned char *flow_order, int *flow_sig) {
	    setRef(num_ref, ref_seq, ref_mer);
	    double x= calc_post_prob(num, flow_order, flow_sig);
	    unsetRef();
	    return x;
	}
	double best_hyp(flow_list *alter_hyp, flow_list *read_list, int cov, float vh, char *line);
	void setRef(int n, unsigned char *rs, int *rm) { 
	    ref = rs; rlen = n; ref_r = rm;
	    int i;
	    for (i = 0; i < n; i++) {
		ref_rev[i] = 3-ref[n-i-1];
		ref_r_rev[i] = ref_r[n-i-1];
	    }
	}
	void unsetRef() { setRef(0, NULL, NULL);}
    protected:
	unsigned char *ref, ref_rev[1000];
	int *ref_r, ref_r_rev[1000];
	int rlen;
	double prob_match(int ref_len, int flow_len, unsigned char b);
	double prob_match(unsigned char refb, int ref_len, unsigned char flow_base, int flow_sign);
	double prob_ins(unsigned char refbase, int ref_len, unsigned char insbase, int flow_len);
	double prob_del(int len);
	double temp[MAX_LEN];
	double split[4];
	int c_sign[4];
	double length_factor(int l) {
	    if (l <= 5)
		//return ((double) (l))*0.005+0.055;
		//return 0.06;
		return 0.07-0.006*(l-1); 
	    if (l > 10) return 0.005;
	    return 0.02;
	}
	/*double length_factor(int l) {
	    double x = 0.075-((double) l)*0.005;
	    if (l >=5) x-=0.2;
	    if (x < 0.01) x = 0.01;
	    return x;
	}
	*/
	int do_split;
};

int rescorer_main(int argc, char **argv);
void rescorer_init(int nt);
double calScore(int pid);
void rescorer_end();
int addRef(int numRefBases, char* refBases);
void addVariant(int pid, int varLocInRef, char *predicted_varaint, char *ref_allele);
void addVariant(int pid, int vc, char *pr);
void addRead(int pid, int varNum, int varFlowIndex, int numFlow, char *flowbase, int *fs, char *al, int ispos);
void finished(int pid);


#endif
