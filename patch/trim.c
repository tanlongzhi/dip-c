#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include <stdio.h>
#include "kvec.h"

/********************
 * Global variables *
 ********************/

const char *lt_adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"; // Illumina 3'-end adapter
const char *lt_adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
const char *lt_oligo_for= "AGATGTGTATAAGAGACAG"; // 19bp transposon
const char *lt_oligo_rev= "CTGTCTCTTATACACATCT";


enum lt_type_e {
	LT_UNKNOWN = 0,
	LT_AMBI_BASE = 1,
	LT_SHORT_SEQ = 2,
	LT_SHORT_PE = 3,
	LT_SHORT_PE_SWAP = 4,
	LT_MERGED = 21,
	LT_NO_MERGE = 22,
};

typedef struct {
	int n_threads;
	int chunk_size;
	int min_seq_len;
	int max_qual;
	int max_ovlp_pen, min_ovlp_len;
	int max_adap_pen, min_adap_len;
	int max_olig_pen, min_olig_len;
	int tab_out;
} lt_opt_t;

static void lt_opt_init(lt_opt_t *opt)
{
	memset(opt, 0, sizeof(lt_opt_t));
	opt->n_threads = 2;
	opt->chunk_size = 10000000;
	opt->max_qual = 50;
	opt->min_seq_len = 40;
	opt->max_ovlp_pen = 2;
	opt->min_ovlp_len = 8;
	opt->max_adap_pen = 1;
	opt->min_adap_len = 3;
	opt->max_olig_pen = 2;
	opt->min_olig_len = 10; // total length = 19
}

/**********************
 * Reverse complement *
 **********************/

char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

void lt_seq_rev(int l, const char *f, char *r)
{
	int i;
	for (i = 0; i < l; ++i)
		r[l - i - 1] = f[i];
	r[l] = 0;
}

void lt_seq_revcomp(int l, const char *f, char *r)
{
	int i;
	for (i = 0; i < l; ++i)
		r[l - i - 1] = (uint8_t)f[i] >= 128? 'N' : comp_tab[(uint8_t)f[i]];
	r[l] = 0;
}

/**********************
 * Ungapped extension *
 **********************/

#define LT_QUAL_THRES 53 // =33+20
#define LT_HIGH_PEN 3
#define LT_LOW_PEN  1

int lt_ue_for1(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int min_len, int max_pen)
{
	int i, pen = 0;
	for (i = 0; i < l1 && i < l2; ++i) {
		if (s1[i] != s2[i]) {
			pen += q1[i] >= LT_QUAL_THRES && (q2 == 0 || q2[i] >= LT_QUAL_THRES)? LT_HIGH_PEN : LT_LOW_PEN;
			if (i <= min_len && pen > max_pen) break;
			if (i > min_len && pen * min_len > i * max_pen) break; // in effect: pen > max_pen * ((double)i / min_len)
		}
	}
	return i;
}

int lt_ue_rev1(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int min_len, int max_pen)
{
	int i, pen = 0;
	for (i = 0; i < l1 && i < l2; ++i) {
		if (s1[l1-1-i] != s2[l2-1-i]) {
			pen += q1[l1-1-i] >= LT_QUAL_THRES && (q2 == 0 || q2[l2-1-i] >= LT_QUAL_THRES)? LT_HIGH_PEN : LT_LOW_PEN;
			if (i <= min_len && pen > max_pen) break;
			if (i > min_len && pen * min_len > i * max_pen) break;
		}
	}
	return i;
}

int lt_ue_for(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int max_pen, int min_len, int max_pos, uint64_t *pos)
{
	int i, n = 0;
	for (i = min_len; i <= l1; ++i) {
		int l;
		l = lt_ue_for1(i, s1 + l1 - i, q1 + l1 - i, l2, s2, q2, min_len, max_pen);
		if (l >= min_len && (l == i || l == l2)) {
			pos[n++] = (uint64_t)(l1 - i) << 32 | l;
			if (n == max_pos) return n;
		}
	}
	return n;
}

int lt_ue_rev(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int max_pen, int min_len, int max_pos, uint64_t *pos)
{
	int i, n = 0;
	for (i = min_len; i <= l1; ++i) {
		int l;
		l = lt_ue_rev1(i, s1, q1, l2, s2, q2, min_len, max_pen);
		if (l >= min_len && (l == i || l == l2)) {
			pos[n++] = (uint64_t)(l1 - i) << 32 | l;
			if (n == max_pos) return n;
		}
	}
	return n;
}

int lt_ue_contained(int l1, const char *s1, const char *q1, int l2, const char *s2, const char *q2, int max_pen, int max_pos, uint64_t *pos)
{
	int i, n = 0;
	for (i = 1; i < l2 - l1; ++i) {
		int l;
		l = lt_ue_for1(l1, s1, q1, l2 - i, s2 + i, q2 + i, l1, max_pen);
		if (l == l1) {
			pos[n++] = (uint64_t)i << 32 | l;
			if (n == max_pos) return n;
		}
	}
	return n;
}

/**********************
 * Batch FASTQ reader *
 **********************/

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	uint32_t l_seq:31, dbl_bind:1;
	enum lt_type_e type;
	int olig_pos_f, olig_pos_r;
	char *name, *seq, *qual, *bc_f, *bc_r;
} bseq1_t;

bseq1_t *bseq_read(kseq_t *ks, int chunk_size, int *n_)
{
	int size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->seq = strdup(ks->seq.s);
		s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		s->bc_f = 0;
		s->bc_r = 0;
		s->l_seq = ks->seq.l;
		s->dbl_bind = 0;
		s->type = LT_UNKNOWN;
		s->olig_pos_f = 0;
		s->olig_pos_r = 0;
		size += seqs[n++].l_seq;
		if (size >= chunk_size && (n&1) == 0) break;
	}
	*n_ = n;
	return seqs;
}

/*********************************
 * Core trimming/merging routine *
 *********************************/

typedef struct {
	lt_opt_t opt;
	kseq_t *ks;
} lt_global_t;

void lt_global_init(lt_global_t *g)
{
	memset(g, 0, sizeof(lt_global_t));
	lt_opt_init(&g->opt);
}

// trim a read by l basepairs from the 5' end
static inline void trim_bseq_5(bseq1_t *s, int l)
{
	memmove(s->seq, s->seq + l, s->l_seq - l);
	memmove(s->qual, s->qual + l, s->l_seq - l);
	s->l_seq -= l;
	s->seq[s->l_seq] = s->qual[s->l_seq] = 0;
}

static inline int merge_base(int max_qual, char fc, char fq, char rc, char rq)
{
	int y;
	if (fc == rc) {
		int q = fq > rq? (fq - 33) + (rq - 33) / 2 : (rq - 33) + (fq - 33) / 2;
		y = toupper(fc) | (33 + (q < max_qual? q : max_qual)) << 8;
	} else {
		if (fq > rq) y = toupper(fc) | (33 + (fq - rq)) << 8;
		else y = toupper(rc) | (33 + (rq - fq)) << 8;
	}
	return y;
}

static inline void trim_adap(bseq1_t *s, const char *adap, int is_5, int min_len, int max_pen, int allow_contained)
{
	int n_hits, l_adap;
	uint64_t hits[4];
	l_adap = strlen(adap);
	if (is_5) n_hits = lt_ue_rev(s->l_seq, s->seq, s->qual, l_adap, adap, 0, max_pen, min_len, 4, hits);
	else n_hits = lt_ue_for(s->l_seq, s->seq, s->qual, l_adap, adap, 0, max_pen, min_len, 4, hits);
	if (n_hits > 0 && (allow_contained || (hits[0]>>32) + (uint32_t)hits[0] == s->l_seq || (hits[n_hits-1]>>32) + (uint32_t)hits[n_hits-1] == s->l_seq)) {
		int len = s->l_seq - (hits[n_hits-1]>>32); // trim the longest hit
		if (is_5) {
			if (len > min_len) trim_bseq_5(s, len); // trim
			else memset(s->qual, 33+1, len); // reduce baseQ
		} else {
			if (len > min_len) s->l_seq -= len, s->seq[s->l_seq] = s->qual[s->l_seq] = 0; // trim
			else memset(s->qual + (s->l_seq - len), 33+1, len); // reduce baseQ
		}
	}
}

// modified from trim_adap, but for finding oligos anywhere in the sequence, store the trimming position and trimmed bases (if needed)
static inline void trim_olig(bseq1_t *s, const char *adap, int is_5, int min_len, int max_pen, int *trim_pos, char *trim_seq, int if_output)
{
	int n_hits, l_adap;
	uint64_t hits[16];
	l_adap = strlen(adap);
	if (is_5) n_hits = lt_ue_rev(s->l_seq, s->seq, s->qual, l_adap, adap, 0, max_pen, min_len, 16, hits);
	else n_hits = lt_ue_for(s->l_seq, s->seq, s->qual, l_adap, adap, 0, max_pen, min_len, 16, hits);
	if (n_hits > 0) {
		int len = s->l_seq - (hits[n_hits-1]>>32); // trim the longest hit
		if (is_5) {
			if (len > min_len) { // trim
				if (if_output) {
					strncpy(trim_seq, s->seq, len);
					trim_seq[len] = 0;
					*trim_pos = len;
				}
				trim_bseq_5(s, len); 
			} else {
				memset(s->qual, 33+1, len); // reduce baseQ
			}
		} else {
			if (len > min_len) {  // trim
				s->l_seq -= len, s->seq[s->l_seq] = s->qual[s->l_seq] = 0;
			} else {
				memset(s->qual + (s->l_seq - len), 33+1, len); // reduce baseQ
			}
		}
	}
}

void lt_process(const lt_global_t *g, bseq1_t s[2])
{
	int i, k, mlen;
	mlen = s[0].l_seq > s[1].l_seq? s[0].l_seq : s[1].l_seq;

	// trim heading and trailing N
	for (k = 0; k < 2; ++k) {
		bseq1_t *sk = &s[k];
		for (i = sk->l_seq - 1; i >= 0; --i) // trim trailing "N"
			if (sk->seq[i] != 'N') break;
		sk->l_seq = i + 1;
		sk->seq[sk->l_seq] = sk->qual[sk->l_seq] = 0;
		for (i = 0; i < sk->l_seq; ++i) // trim heading "N"
			if (sk->seq[i] != 'N') break;
		if (i) trim_bseq_5(sk, i);
	}
	
	// trim Illumina PE adapters
	trim_adap(&s[0], lt_adapter1, 0, g->opt.min_adap_len, g->opt.max_adap_pen, 1);
	trim_adap(&s[1], lt_adapter2, 0, g->opt.min_adap_len, g->opt.max_adap_pen, 1);
	
	// trim transposon sequences and store barcodes
	int olig_pos[2];
	char *bc[2];
	
	for (k = 0; k < 2; ++k) {
		olig_pos[k] = 0;
		bc[k] = (char*)alloca(mlen + 1);
		bc[k][0] = 0;
		trim_olig(&s[k], lt_oligo_for, 1, g->opt.min_olig_len, g->opt.max_olig_pen, &olig_pos[k], bc[k], 1);
		trim_olig(&s[k], lt_oligo_rev, 0, g->opt.min_olig_len, g->opt.max_olig_pen, 0, 0, 0);
        trim_adap(&s[k], lt_oligo_rev, 0, g->opt.min_adap_len, g->opt.max_adap_pen, 1);
	}

	
	// merge the two ends
	int n_fh;
	uint64_t fh[2];
	char *rseq, *rqual, *xseq, *xqual;
	rseq = (char*)alloca(mlen + 1);
	rqual = (char*)alloca(mlen + 1);
	xseq = (char*)alloca(s[0].l_seq + s[1].l_seq + 1);
	xqual = (char*)alloca(s[0].l_seq + s[1].l_seq + 1);
	// reverse the other read
	lt_seq_revcomp(s[1].l_seq, s[1].seq, rseq);
	lt_seq_rev(s[1].l_seq, s[1].qual, rqual);
	// find overlaps
	n_fh = lt_ue_for(s[0].l_seq, s[0].seq, s[0].qual, s[1].l_seq, rseq, rqual, g->opt.max_ovlp_pen, g->opt.min_ovlp_len, 2, fh);
	if (n_fh == 1) {
    	int l = (uint32_t)fh[0], st = fh[0]>>32;
        if (st + l == s[0].l_seq) { // good stitch
            //printf("%s\n%s\n%s\n%s\n%d,%d,%d,%d\n",s[0].seq,s[0].qual,rseq,rqual,st,l,s[0].l_seq,s[1].l_seq);
    		int x = 0;
    		s[0].type = s[1].type = LT_MERGED;
    		for (i = 0; i < st; ++i)
    			xseq[x] = s[0].seq[i], xqual[x++] = s[0].qual[i];
    		for (i = 0; i < l; ++i) {
    			int j = st + i, y;
    			y = merge_base(g->opt.max_qual, s[0].seq[j], s[0].qual[j], rseq[i], rqual[i]);
    			xseq[x] = (uint8_t)y, xqual[x++] = y>>8;
    		}
    		for (i = l; i < s[1].l_seq; ++i)
    			xseq[x] = rseq[i], xqual[x++] = rqual[i];
            //printf("%s\n%s\n%d\n",xseq,xqual,x);
    		xseq[x] = xqual[x] = 0;
    		free(s[0].seq); free(s[0].qual);
    		s[0].seq = strdup(xseq);
    		s[0].qual = strdup(xqual);
    		s[0].l_seq = x;
            s[1].seq[0] = 0;
    		s[1].l_seq = 0;
        } else {
            s[0].type = s[1].type = LT_NO_MERGE;
        }
	} else {
		s[0].type = s[1].type = LT_NO_MERGE;
	}
	
	// check whether read length is too short
    if (s[0].l_seq < g->opt.min_seq_len) {
        if (s[0].type == LT_MERGED) { // single-end
            s[0].type = s[1].type = LT_SHORT_SEQ;
        } else if (s[1].l_seq >= g->opt.min_seq_len) { // read 1 is short, but read 2 is not => swap
			char *tmp;
			tmp = s[0].seq, s[0].seq = s[1].seq, s[1].seq = tmp;
			tmp = s[0].qual, s[0].qual = s[1].qual, s[1].qual = tmp;
			tmp = bc[0], bc[0] = bc[1], bc[1] = tmp;
			i = olig_pos[0], olig_pos[0] = olig_pos[1], olig_pos[1] = i;
			i = s[0].l_seq, s[0].l_seq = s[1].l_seq, s[1].l_seq = i;
            s[0].type = s[1].type = LT_SHORT_PE_SWAP;
        } else {
            s[0].type = s[1].type = LT_SHORT_SEQ;
        }
    } else if (s[0].type != LT_MERGED && s[1].l_seq < g->opt.min_seq_len) { // read 2 is short, but read 1 is not
        s[0].type = s[1].type = LT_SHORT_PE;
    }
    
    // record the barcodes
	for (k = 0; k < 2; ++k) {
		s[k].olig_pos_f = olig_pos[0];
		s[k].olig_pos_r = olig_pos[1];
		s[k].bc_f = strdup(bc[0]);
		s[k].bc_r = strdup(bc[1]);
	}
}

/**********************
 * Callback functions *
 **********************/

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	lt_global_t *g;
} data_for_t;

static void worker_for(void *_data, long i, int tid)
{
	data_for_t *data = (data_for_t*)_data;
	lt_process(data->g, &data->seqs[i<<1]);
}

static void *worker_pipeline(void *shared, int step, void *_data)
{
	int i;
	lt_global_t *g = (lt_global_t*)shared;
	if (step == 0) {
		data_for_t *ret;
		ret = calloc(1, sizeof(data_for_t));
		ret->seqs = bseq_read(g->ks, g->opt.chunk_size, &ret->n_seqs);
		assert((ret->n_seqs&1) == 0);
		ret->g = g;
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		data_for_t *data = (data_for_t*)_data;
		kt_for(g->opt.n_threads, worker_for, data, data->n_seqs>>1);
		return data;
	} else if (step == 2) {
		data_for_t *data = (data_for_t*)_data;
		if (g->opt.tab_out) { // tabular output
			for (i = 0; i < data->n_seqs; i += 2) {
				bseq1_t *s = &data->seqs[i];
				bseq1_t *s_r = &data->seqs[i+1];        
				printf("%s\t%d\t%d\t%d\t%zd\t%zd\n", s->name, s->type, s->olig_pos_f, s->olig_pos_r, strlen(s->seq),strlen(s_r->seq));
			}
		} else { // FASTQ output (FASTA not supported yet)
			for (i = 0; i < data->n_seqs; ++i) {
				bseq1_t *s = &data->seqs[i];
				if (s->l_seq > 0 && (s->type == LT_NO_MERGE || s->type == LT_MERGED || (s->type == LT_SHORT_PE && ~i&1)) ) {
					putchar(s->qual? '@' : '>'); fputs(s->name, stdout);
					if (s->type == LT_NO_MERGE) {
						putchar('/'); putchar("12"[i&1]);
					}
					printf(" YT:i:%d", s->type);
					printf("\tPF:i:%d", s->olig_pos_f);
					printf("\tPR:i:%d", s->olig_pos_r);
					if (s->bc_f) { fputs("\tBF:Z:", stdout); fputs(s->bc_f[0] == 0? "*" : s->bc_f, stdout); }
					if (s->bc_r) { fputs("\tBR:Z:", stdout); fputs(s->bc_r[0] == 0? "*" :s->bc_r, stdout); }
					putchar('\n');
					puts(s->seq);
					if (s->qual) { puts("+"); puts(s->qual); }
				}
			}
		}
		for (i = 0; i < data->n_seqs; ++i) { // deallocate
			bseq1_t *s = &data->seqs[i];
			free(s->bc_f); free(s->bc_r); free(s->seq); free(s->qual); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

#include <unistd.h>

int main_trim(int argc, char *argv[])
{
	int c;
	lt_global_t g;
	gzFile fp;

	lt_global_init(&g);
	while ((c = getopt(argc, argv, "Tt:b:l:")) >= 0) {
		if (c == 't') g.opt.n_threads = atoi(optarg);
		else if (c == 'T') g.opt.tab_out = 1;
		else if (c == 'l') g.opt.min_seq_len = atoi(optarg);
	}
	if (argc - optind < 1) {
		fprintf(stderr, "Usage: seqtk mergepe <read1.fq> <read2.fq> | lianti trim [options] -\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t INT     number of threads [%d]\n", g.opt.n_threads);
		fprintf(stderr, "  -l INT     min read/fragment length to output [%d]\n", g.opt.min_seq_len);
		fprintf(stderr, "  -T         tabular output for debugging\n");
		return 1;
	}

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	g.ks = kseq_init(fp);

	kt_pipeline(2, worker_pipeline, &g, 3);

	kseq_destroy(g.ks);
	gzclose(fp);
	return 0;
}
