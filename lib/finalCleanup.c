//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
#include "Python.h"

#define FINALCLEANUP_VERSION "0.0"



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct read {
        char *r_seq, *r_qual, *r_header;
        unsigned int r_len;

};

struct reads {
        struct read r1;
        struct read r2;
};



struct stats {

        unsigned long long int A;
        unsigned long long int C;
        unsigned long long int T;
        unsigned long long int G;
        unsigned long long int N;

        //Used double to get decimal points     
        double basePairs;
        double qualTotal;

        //discarded reads
        unsigned long long int se_discarded;
        unsigned long long int r1_discarded;
        unsigned long long int r2_discarded;

        unsigned long long int se_kept;
        unsigned long long int pe_kept;
        unsigned long long int numForcedPairs;

        unsigned long long int polyATrimmed;
        unsigned long long int polyTTrimmed;


        unsigned long int R1_length[700];
        unsigned long int R2_length[700];
        unsigned long int SE_length[700];

};

void statsConstruct(struct stats *s) {

        int i = 0;
        s->A = 0;
        s->T = 0;
        s->G = 0;
        s->C = 0;
        s->N = 0;

        s->basePairs = 0;
        s->qualTotal = 0;

        s->se_discarded = 0;
        s->r1_discarded = 0;
        s->r2_discarded = 0;
        s->se_kept = 0;
        s->pe_kept = 0;
        s->numForcedPairs = 0;

        s->polyATrimmed = 0;
        s->polyTTrimmed = 0;

        for (i = 0; i < 700; i++) {
                s->R1_length[i] = 0;
                s->R2_length[i] = 0;
                s->SE_length[i] = 0;
        }

}

int grabFour(FILE *f, struct read *r) {
        char data[567];
        int i = 0;
        r->r_header = NULL;

        while (fgets(data, 567, f) != NULL) {
                if (i == 0) {
                        r->r_header = strdup(data);
                } else if (i == 1) {
                        r->r_seq = strdup(data);
                } else if (i == 3) {
                        r->r_qual = strdup(data);
                        r->r_len = strlen(data);
                        return 0;
                }
                i++;
        }

        return 1;

}


int trim_poly_a(char *seq, int min_cutoff, int trim_errors, size_t start, size_t end) {

        int cut_off = 0, errors = 0, tmp_add = 0, index = 0;

        for (index = start; index >= end; index--) {
                if (seq[index] == 'A') {
                        cut_off++;
                        cut_off += tmp_add;
                        tmp_add = 0;
                } else {
                        errors++;
                        tmp_add++;
                }

                if (errors >= trim_errors) {
                        break;
                }
        }

        if (min_cutoff > cut_off) {
                return 0;
        } else {
                return cut_off;
        }
}


int trim_poly_t(char *seq, int min_cutoff, int trim_errors, size_t start, size_t end) {

        int cut_off = 0, errors = 0, tmp_add = 0, index;

        for (index = start; index < end; index++) {
                if (seq[index] == 'T') {
                        cut_off++;
                        cut_off += tmp_add;
                        tmp_add = 0;
                } else {
                        errors++;
                        tmp_add++;
                }

                if (errors >= trim_errors) {
                        break;
                }
        }

        if (min_cutoff > cut_off) {
                return 0;
        } else {
                return cut_off;
        }
}



void split(char *strTest, char *delimiter, struct reads *r) {

        //header 1 , read 1, read 2, qual 1, qual 2
        char **results = (char **)malloc(sizeof(char *) * 5);
        char *token;
        int count = 0;
        int i = 0;

        while ((token = strsep(&strTest, delimiter)) != NULL) {
                results[count] = strdup(token);
                count++;
        }


        if (count == 5) {
                        (r->r1).r_header = strdup(results[0]);
                        (r->r1).r_seq = strdup(results[1]);
                        (r->r1).r_qual = strdup(results[2]);
                        (r->r1).r_len = strlen((r->r1).r_seq);

                        (r->r2).r_header = strdup(results[0]);
                        (r->r2).r_seq = strdup(results[3]);
                        (r->r2).r_qual = strdup(results[4]);
                        (r->r2).r_len = strlen((r->r2).r_seq);

        } else {
                (r->r1).r_header = strdup(results[0]);
                (r->r1).r_seq = strdup(results[1]);
                (r->r1).r_qual = strdup(results[2]);
                (r->r1).r_len = strlen((r->r1).r_seq);
                (r->r2).r_header = NULL;

        }


        //free(token);


        for (i = 0; i < count; i++) {
               free(results[i]);
        }

}

int minLen = 30;




void readStats(struct read r, struct stats *s) {
        int i = 0;
        for (i = 0; i < r.r_len + 1; i++) {
                s->qualTotal += r.r_qual[i] - 33;
                s->basePairs++;

                char c = (r.r_seq)[i];

                if (c == 'A') { (s->A)++; }
                else if (c == 'T') { (s->T)++; }
                else if (c == 'C') { (s->C)++; }
                else if (c == 'G') { (s->G)++; }
                else if (c == 'N') { (s->N)++;}

        }

}


//Checks length first and then take length stats
void getStats(struct reads *r, struct stats *s) {

        if ((r->r1).r_len < minLen && (r->r2).r_header == NULL) {
                (s->se_discarded)++;
                (r->r1).r_header = NULL;
        } else if ((r->r1).r_len < minLen) {
                (s->r1_discarded)++;
                (r->r1).r_header = NULL;
		(r->r2).r_header = NULL;
        }

        if ((r->r2).r_header != NULL && (r->r2).r_len < minLen) {
                (s->r2_discarded)++;
                (r->r1).r_header = NULL;
                (r->r2).r_header = NULL;
        }


        if ((r->r1).r_header != NULL) {
                readStats(r->r1, s);
	}

        if ((r->r2).r_header != NULL) {
                readStats(r->r1, s);
	}




}



void PolyATCuts(struct reads *r, int minLenOfTail, int errorsAllowed, struct stats *s) {
        int cut = trim_poly_t((r->r1).r_seq, minLenOfTail, errorsAllowed, 0, (r->r1).r_len);

        if (cut != 0) {
                s->polyTTrimmed += 1;
                (r->r1).r_qual += cut;
                (r->r1).r_seq += cut;
                (r->r1).r_len -= cut;
                //Magic number for sure but if it gets close to the end - poly T can only be R2 if it spans R1
                if ((r->r1).r_len < 5 && (r->r2).r_header != NULL) {
			(r->r1) = (r->r2);
			(r->r2).r_header = NULL;
                        cut = trim_poly_t((r->r1).r_seq, minLenOfTail, errorsAllowed, 0, (r->r1).r_len);
                        if (cut != 0) {
                                s->polyTTrimmed += cut;
                                (r->r1).r_qual += cut;
                                (r->r1).r_seq += cut;
                                (r->r1).r_len -= cut;
                        }
                }

        }

	if ((r->r2).r_header != NULL) {
	        cut = trim_poly_a((r->r2).r_seq, minLenOfTail, errorsAllowed, (r->r2).r_len, 0);
        	if (cut != 0) {

                	s->polyATrimmed += 1;
	                (r->r2).r_seq[(r->r2).r_len - cut + 1] = '\0';
			(r->r2).r_qual[(r->r2).r_len - cut + 1] = '\n';

	                (r->r2).r_seq[(r->r2).r_len - cut + 2] = '\0';
			(r->r2).r_qual[(r->r2).r_len - cut + 2] = '\0';
			(r->r2).r_len -= cut;

			if ((r->r2).r_len < 5) {
				(r->r2).r_header = NULL;
				cut = trim_poly_a((r->r1).r_seq, minLenOfTail, errorsAllowed, (r->r1).r_len, 0);

				if (cut != 0) {
					s->polyATrimmed += cut;
					(r->r1).r_seq[(r->r1).r_len - cut + 1] = '\0';
					(r->r1).r_seq[(r->r1).r_len - cut + 1] = '\0';
					(r->r1).r_len -= cut;
                       		}
                	}
        	}
	}

}

int PolyATTrim = 0;


int grabTab(FILE *f, struct reads *r, struct stats *s) {
        char data[4096];
        int i = 0;

        if (fgets(data, sizeof(data), f) != NULL) {
                split(data, "\t", r);
		if (PolyATTrim) {
                	PolyATCuts(r, 10, 3, s);
                }
		if ((r->r1).r_header != NULL) {
			getStats(r, s);
		}

	        return 1;
        } else {
                return 0;
        }



}

char compliment(char c) {
        switch (c) {
                case 'A':
                        return 'T';
                        break;
                case 'T':
                        return 'A';
                        break;
                case 'G':
                        return 'C';
                        break;
                case 'C':
                        return 'G';
                        break;
                default:
                        return c;
                        break;
        }
}

char * reverseComp(char *read) {
        char *strReverse = strdup(read);
        int i = 0, end = strlen(read)-1;

        while (strReverse[i] != '\n' &&  strReverse[i] != '\0' && strReverse[i] != '\t') {
		strReverse[i] = compliment(read[end]);
                i++;
                end--;
        }

        return strReverse;


}

char *reverse(char *read) {

        char *strReverse = strdup(read);
        int i = 0, end = strlen(read)-2;

        while (strReverse[i] != '\n' &&  strReverse[i] != '\0' && strReverse[i] != '\t') {
                strReverse[i] = read[end];
                i++;
                end--;
        }

	if (strReverse[i] !=  '\n') {
		strReverse[i] = '\n';
	}
	strReverse[++i] = '\0';

        return strReverse;


}



int clean(char *devFile, char *logFile, char *strR1, char *strR2, char *strSE, int tmpforcePairs, int tmpPolyATTrim) {
        int forcePairs = tmpforcePairs;
	//PolyATTrim = tmpPolyATTrim;
	//PolyATTrim = tmpPolyATTrim;

        FILE *f = fopen(devFile, "r");
        //FILE *f = stderr;
        FILE *R1 = NULL;
        FILE *R2 = NULL;
        FILE *SE = NULL;
	FILE *log = NULL;

        struct reads r;
        struct stats s;
        statsConstruct(&s);

        while (grabTab(f, &r, &s)) {
                if ((r.r2).r_header != NULL && (r.r1).r_header != NULL) {
                        s.pe_kept++;

			if (strlen((r.r1).r_qual) == 0) {
				printf("2. im going to destory you\n");
			}
			
			if (R1 == NULL) {
        			R1 = fopen(strR1, "w");
			}		
			if (R2 == NULL) {
        			R2 = fopen(strR2, "w");
			}
			fprintf(R1, "@N%s\n%s\n+\n%s\n", (r.r1).r_header, (r.r1).r_seq, (r.r1).r_qual);
			fprintf(R2, "@N%s\n%s\n+\n%s", (r.r2).r_header, (r.r2).r_seq, (r.r2).r_qual);
                } else if (forcePairs && (r.r1).r_header != NULL) {
                        s.numForcedPairs++;
                        int loc = (strlen((r.r1).r_seq))/2;
                        char cSeq = (r.r1).r_seq[loc];
                        char cQual = (r.r1).r_qual[loc];

			if (R1 == NULL) {
        			R1 = fopen(strR1, "w");
			}		
			if (R2 == NULL) {
        			R2 = fopen(strR2, "w");
			}

                        (r.r1).r_seq[loc] = '\0';
                        (r.r1).r_qual[loc] = '\0';
                        fprintf(R1, "@%s\n%s\n+\n%s\n", (r.r1).r_header, (r.r1).r_seq, (r.r1).r_qual);
                        (r.r1).r_seq[loc] = cSeq;
                        (r.r1).r_qual[loc] = cQual;

                        fprintf(R2, "@%s\n%s\n+\n%s", (r.r1).r_header, reverseComp(&((r.r1).r_seq)[loc]), reverse(&((r.r1).r_qual)[loc]));
			
                } else if ((r.r2).r_header != NULL) {
			if (SE == NULL) {
                		SE = fopen(strSE, "w");
			}
                        s.se_kept++;
                        fprintf(SE, "@%s\n%s\n+\n%s", (r.r2).r_header, (r.r2).r_seq, (r.r2).r_qual);
                } else if ((r.r1).r_header != NULL) {
			if (SE == NULL) {
                		SE = fopen(strSE, "w");
			}
                        s.se_kept++;
                        fprintf(SE, "@%s\n%s\n+\n%s", (r.r1).r_header, (r.r1).r_seq, (r.r1).r_qual);
                }

        }

	if (log == NULL) {
		log = fopen(logFile, "w");
	}

        fprintf(log, "A\tT\tG\tC\tN\tPolyA_Removed_Reads\tPolyT_Removed_Reads\tShort_discarded\tPE_Kept\tSE_Kept\tForced_Pairs\tAverageQual\n");
        fprintf(log, "%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%.2f\n", s.A, s.T, s.G, s.C, s.N, s.polyATrimmed, s.polyTTrimmed, s.r1_discarded + s.r2_discarded + s.se_discarded, s.pe_kept, s.se_kept, s.numForcedPairs,
		(float)((float)(s.qualTotal)/(float)(s.A + s.T + s.C + s.G + s.N)));

	if (f != NULL) {
        	fclose(f);
        }

	if (R1 != NULL) {
		fclose(R1);
        	fclose(R2);
	}


        if (SE != NULL) {
                fclose(SE);
        }

}
   

static PyObject *
	finalCleanup(PyObject *self, PyObject *args) {
	const char *logFile, *devFile, *R1, *R2, *SE;
	int logFileSize, devFileSize, r1_size, r2_size, se_size;
	int polyAT = 0, split = 0;

	if (!PyArg_ParseTuple(args, "iiis#s#s#s#s#", &polyAT, &split, &minLen, &logFile, &logFileSize, &devFile, &devFileSize, &R1, &r1_size, &R2, &r2_size, &SE, &se_size)) {
		return NULL;		
	}
	clean(devFile, logFile, R1, R2, SE, split, polyAT);
	return args;
}



PyDoc_STRVAR(finalCleanup_doc,
"bounded_edit_distance(a, b, k, m) -> int, int\n\
    Calculates the bounded Levenshtein's edit distance between strings \"a\" and \"b\" with bound \"k\" and \"m\" matching bases at end\n");

static PyMethodDef finalCleanup_methods[] = {
    {   "finalCleanup", (PyCFunction)finalCleanup,
        METH_VARARGS,   finalCleanup_doc       },
    {   NULL, NULL, 0, NULL }  /* sentinel */
};

PyDoc_STRVAR(module_doc, "Calculate Hamming distance, Levenshtein's edit distance, and a edge bounded Levenshtein's edit distance.\n");

PyMODINIT_FUNC
initfinalCleanup(void)
{
    PyObject *m;

    m = Py_InitModule3("finalCleanup", finalCleanup_methods, module_doc);

    if (m == NULL) {
	printf("Null value in init\n");
	return NULL;
    }

    PyModule_AddStringConstant(m, "__version__", FINALCLEANUP_VERSION);
}

