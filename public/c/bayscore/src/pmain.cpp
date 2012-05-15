#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include "fasta-io.h"
#include "zutil.h"

#include "posterior_flow.h"

static void getlist(char *line, int *fs, int nf)
{
	int i;
	for (i = 0; i < nf; i++) {
		fs[i] = atoi(line);
		line = strchr(line, ',');
		if (!line) break;
		line++;
	}
	if (i < nf-1) fatal("to few flow");
}

int main(int argc, char **argv)
{
	char *file = argv[1];
	FILE *fp = fopen(file, "r");
	FILE *fvcf = fopen(argv[2], "r");
	FILE *outp = fopen(argv[3], "w");
	char line[1000000];
	char ref[100000];
	rescorer_init(1);
	char fb[100000];
    	int fs[100000];
	while (fgets(line, sizeof line, fp)) {
		fprintf(outp, "%s", line);
		char c = line[5];
		char *s = strchr(line, ':');
		if (!s) continue;
		if (strstr(line, "return")) continue;
		if (c == 'f') {
			finished(0);
			char line2[10000];
			while (fgets(line2, sizeof line2, fvcf)) {
			    if (line2[0] == '#') continue;
			    printf("%s", line2);
			    fprintf(outp, "%s", line2);
			    break;
			}
			s++;
			int len;
			sscanf(s, "%d;%s", &len, ref);
			printf("%d %s\n", len, ref);
			if (len != strlen(ref)) fatal("length of the reference is not correct %d %d %s\n", len, strlen(ref), ref);
			int x = addRef(strlen(ref), ref);
			if (x != 0) fatal("pid not right\n");
		} else if (c == 'r') {
			int pos;
			char pv[10000];
			sscanf(s+1, "%*d;%d;%s", &pos, pv);
			printf("%d %s\n", pos, pv);
			addVariant(0, pos, pv);
		} else if (c == 'o') {
			double bs = calScore(0);
                        printf("Bayesian score=%f\n", (float) bs);
                        fprintf(outp, "Bayesian score=%f\n", (float) bs);
		} else if (c == 'a') {
			int vn, vf, nf;
			char *seq;
			sscanf(s+1, "%*d;%d;%d;%d;", &vn, &vf, &nf);
			s = strchr(s, ';');
			s = strchr(s+1, ';');
			s = strchr(s+1, ';');
                        s = strchr(s+1, ';');
			seq = s+1;
			char *t = strchr(seq, ';');
			*t= 0;
			t++;
			getlist(t, fs, nf);
			t = strchr(t, ';');
			t++;
			char *align = t;
			t = strchr(t, ';');
			*t = 0;
			int ispos = atoi(t+1);
			addRead(0, vn, vf, nf, seq, fs, align, ispos);
		} 
	}
	rescorer_end();	
}
