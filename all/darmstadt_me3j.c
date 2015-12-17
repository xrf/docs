#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int min(int x, int y)
{
    return x < y ? x : y;
}

void swap(int *x, int *y)
{
    const int z = *x;
    *x = *y;
    *y = z;
}

/* Given two ranges x1 to x2 and y1 to y2, find their intersection and union.
   The intersection is stored in x1 and x2 and the union in y1 and y2 */
void intersect_range(int *x1, int *x2, int *y1, int *y2)
{
    if (*x1 < *y1) {
        swap(x1, y1);
    }
    if (*x2 > *y2) {
        swap(x2, y2);
    }
}

/** Find the range of J formed by coupling j1 with j2. */
void range_twoJ(int twoj1, int twoj2, int *twoJMin, int *twoJMax)
{
    /* Use the triangle condition. */
    *twoJMin = abs(twoj1 - twoj2);
    *twoJMax = (twoj1 + twoj2);
}

/** Find the common range of J formed by coupling j1 with j2 and coupling jj1
    with jj2. */
void range_twoJ_both(int twoj1, int twoj2, int twojj1, int twojj2,
                     int *twoJMin, int *twoJMax)
{
    int twoJJMin, twoJJMax;
    range_twoJ(twoj1, twoj2, twoJMin, twoJMax);
    range_twoJ(twojj1, twojj2, &twoJJMin, &twoJJMax);
    intersect_range(twoJMin, twoJMax, &twoJJMin, &twoJJMax);
}

/** Find the range of J formed by coupling j1 with j2, where j1 and j2 must be
    either both half-integers or both integers. */
void range_J(int twoj1, int twoj2, int *JMin, int *JMax)
{
    range_twoJ(twoj1, twoj2, JMin, JMax);
    *JMin /= 2;
    *JMax /= 2;
}

/** Iterate over `n`, `l`, `j` of all single-particle states.  The function
    `f` is called each iteration with an arbitary argument `env`. */
void iterate_nlj(int eMax, int nMax, int lMax,
                 void (*f)(void *env, int e, int n, int l, int twoj),
                 void *env)
{
    int e, n, l, lMin, twoj, twojMin, twojMax;
    for (e = 0; e <= eMax; e++) {
        lMin = e % 2;
        for (l = lMin; l <= min(e, lMax); l += 2) {
            n = (e - l) / 2;
            if (n > nMax)
                continue;
            range_twoJ(2 * l, 1, &twojMin, &twojMax);
            for (twoj = twojMin; twoj <= twojMax; twoj += 2) {
                (*f)(env, e, n, l, twoj);
            }
        }
    }
}

static void get_nljDim_iteratee(void *env, int e, int n, int l, int twoj)
{
    int *nlj = (int *)env;
    ++*nlj;
}

int get_nljDim(int eMax, int nMax, int lMax)
{
    int nlj = 0;
    iterate_nlj(eMax, nMax, lMax, &get_nljDim_iteratee, &nlj);
    return nlj;
}

struct nlj_table_iteratee_env {
    int nlj, *e_nlj, *n_nlj, *l_nlj, *twoj_nlj;
};

static void nlj_table_iteratee(void *env, int e, int n, int l, int twoj)
{
    struct nlj_table_iteratee_env *self = (struct nlj_table_iteratee_env *)env;
    self->e_nlj[self->nlj] = e;
    self->n_nlj[self->nlj] = n;
    self->l_nlj[self->nlj] = l;
    self->twoj_nlj[self->nlj] = twoj;
    ++self->nlj;
}

/** Initialize arrays that map `nlj` to `e`, `n`, `l`, `twoj` respectively. */
void nlj_table(int eMax, int nMax, int lMax,
               int *e_nlj, int *n_nlj, int *l_nlj, int *twoj_nlj)
{
    struct nlj_table_iteratee_env self;
    self.nlj = 0;
    self.e_nlj = e_nlj;
    self.n_nlj = n_nlj;
    self.l_nlj = l_nlj;
    self.twoj_nlj = twoj_nlj;
    iterate_nlj(eMax, nMax, lMax, &nlj_table_iteratee, &self);
}

struct nlj_table {
    int nljDim;
    int *e_nlj, *n_nlj, *l_nlj, *twoj_nlj;
};

int init_nlj_table(struct nlj_table *self, int eMax, int nMax, int lMax)
{
    self->nljDim = get_nljDim(eMax, nMax, lMax);
    self->e_nlj = (int *)malloc(sizeof(*self->e_nlj) *
                                (size_t)self->nljDim * 4);
    if (!self->e_nlj) {
        return 1;
    }
    self->n_nlj = self->e_nlj + self->nljDim;
    self->l_nlj = self->n_nlj + self->nljDim;
    self->twoj_nlj = self->l_nlj + self->nljDim;
    nlj_table(eMax, nMax, lMax, self->e_nlj,
              self->n_nlj, self->l_nlj, self->twoj_nlj);
    return 0;
}

void free_nlj_table(struct nlj_table *self)
{
    free(self->e_nlj);
}

int violates_parity_3(int nlj1, int nlj2, int nlj3,
                      int nnlj1, int nnlj2, int nnlj3,
                      const int *l_nlj)
{
    return (l_nlj[nlj1] + l_nlj[nlj2] + l_nlj[nlj3] -
            l_nlj[nnlj1] - l_nlj[nnlj2] - l_nlj[nnlj3]) % 2;
}

/** Iterate over all 3-body matrix elements in the Darmstadt ME3J order.  The
    function `f` is called each iteration with an arbitary argument `env`. */
void iterate_me3j(int nljDim, int EMax,
                  const int *e_nlj, const int *l_nlj, const int *twoj_nlj,
                  void (*f)(void *env,
                            int nlj1, int nlj2, int nlj3,
                            int nnlj1, int nnlj2, int nnlj3,
                            int J12, int JJ12, int J,
                            int T12, int TT12, int T),
                  void *env)
{
    int nlj1, nlj2, nlj3, nnlj1, nnlj2, nnlj3,
        J12, JJ12, J12Min, J12Max, JJ12Min, JJ12Max,
        twoJ, twoJMin, twoJMax, T12, TT12, twoT;
    for (nlj1 = 0; nlj1 < nljDim; nlj1++)
    for (nlj2 = 0; nlj2 <= nlj1; nlj2++)
    for (nlj3 = 0; nlj3 <= nlj2; nlj3++) {
        if (e_nlj[nlj1] + e_nlj[nlj2] + e_nlj[nlj3] > EMax)
            break;
        for (nnlj1 = 0; nnlj1 <= nlj1; nnlj1++)
        for (nnlj2 = 0; nnlj2 <= (nlj1 == nnlj1 ? nlj2 : nnlj1); nnlj2++)
        for (nnlj3 = 0; nnlj3 <= (nlj1 == nnlj1 && nlj2 == nnlj2 ?
                                  nlj3 : nnlj2); nnlj3++) {
            if (e_nlj[nnlj1] + e_nlj[nnlj2] + e_nlj[nnlj3] > EMax)
                break;
            if (violates_parity_3(nlj1, nlj2, nlj3, nnlj1, nnlj2, nnlj3, l_nlj))
                continue;
            range_J(twoj_nlj[nlj1], twoj_nlj[nlj2], &J12Min, &J12Max);
            range_J(twoj_nlj[nnlj1], twoj_nlj[nnlj2], &JJ12Min, &JJ12Max);
            for (J12 = J12Min; J12 <= J12Max; J12++)
            for (JJ12 = JJ12Min; JJ12 <= JJ12Max; JJ12++) {
                range_twoJ_both(2 * J12, twoj_nlj[nlj3],
                                2 * JJ12, twoj_nlj[nnlj3],
                                &twoJMin, &twoJMax);
                for (twoJ = twoJMin; twoJ <= twoJMax; twoJ += 2)
                for (T12 = 0; T12 <= 1; T12++)
                for (TT12 = 0; TT12 <= 1; TT12++)
                for (twoT = 1; twoT <= 2 * min(T12, TT12) + 1; twoT += 2) {
                    (*f)(env, nlj1, nlj2, nlj3, nnlj1, nnlj2, nnlj3,
                         J12, JJ12, twoJ, T12, TT12, twoT);
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Example program */
/* ------------------------------------------------------------------------ */

static void main_me3j_iteratee(void *env,
                               int nlj1, int nlj2, int nlj3,
                               int nnlj1, int nnlj2, int nnlj3,
                               int J12, int JJ12, int J,
                               int T12, int TT12, int T)
{
    FILE *stream = (FILE *)env;
    fprintf(stream, "%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
            nlj1, nlj2, nlj3, nnlj1, nnlj2, nnlj3,
            J12, JJ12, J, T12, TT12, T);
}

static void usage(const char *name, FILE *stream)
{
    fprintf(stream,
            "usage: %s [-o file] [-E EMax] [-l lMax] [-n nMax] eMax\n"
            "EMax/lMax/nMax default to the largest possible value\n",
            name);
}

static int get_options(int *argc, char ***argv, char *opt, const char **val,
                       const char *optstring)
{
    if (!*argc) {
        return 0;
    }
    if (!(**argv)[0] || (**argv)[0] != '-') {
        *opt = '\0';
        *val = **argv;
        ++*argv;
        --*argc;
        return -1;
    }
    *opt = (**argv)[1];
    if (!*opt) {
        fprintf(stderr, "invalid option: %s\n", **argv);
        return 1;
    }
    if (!strchr(optstring, *opt)) {
        fprintf(stderr, "unknown option: %s\n", **argv);
        return 1;
    }
    if ((**argv)[2]) {
        fprintf(stderr, "invalid option: %s\n", **argv);
        return 1;
    }
    ++*argv;
    --*argc;
    if (!*argc) {
        fprintf(stderr, "missing value for -%c\n", *opt);
        return 1;
    }
    *val = **argv;
    ++*argv;
    --*argc;
    return -1;
}

static int parse_args(int argc, char **argv,
                      const char **name, const char **filename,
                      int *eMax, int *nMax, int *lMax, int *EMax)
{
    const char *val;
    char opt;
    int e, eMax_set = 0;
    *name = "_";
    if (argc) {
        *name = argv[0];
        ++argv;
        --argc;
    }
    *filename = "darmstadt_me3j_quantum_numbers.txt";
    *nMax = *lMax = *EMax = INT_MAX;
    while ((e = get_options(&argc, &argv, &opt, &val, "Elno")) == -1) {
        switch (opt) {
        case 'E':
            *EMax = atoi(val);
            break;
        case 'l':
            *lMax = atoi(val);
            break;
        case 'n':
            *nMax = atoi(val);
            break;
        case 'o':
            *filename = val;
            break;
        case '\0':
            if (eMax_set) {
                fprintf(stderr, "unexpected argument: %s\n", val);
                return 1;
            }
            *eMax = atoi(val);
            eMax_set = 1;
            break;
        }
    }
    if (!eMax_set) {
        fprintf(stderr, "missing argument: eMax\n");
        return 1;
    }
    return e;
}

int main(int argc, char **argv)
{
    struct nlj_table table;
    const char *name, *filename;
    int eMax, nMax, lMax, EMax;
    FILE *stream;

    if (parse_args(argc, argv, &name, &filename, &eMax, &nMax, &lMax, &EMax)) {
        usage(name, stderr);
        return EXIT_FAILURE;
    }

    printf("Input limits:\n"
           "  eMax = %i\n"
           "  nMax = %i\n"
           "  lMax = %i\n"
           "  EMax = %i\n\n",
           eMax, nMax, lMax, EMax);

    if (init_nlj_table(&table, eMax, nMax, lMax)) {
        fprintf(stderr, "init_nlj_table failed\n");
        return EXIT_FAILURE;
    }

    printf("Quantum numbers of 3-body matrix elements will be written to:\n"
           "  %s\nWriting... ", filename);
    fflush(stdout);

    stream = fopen(filename, "w");
    if (!stream) {
        fprintf(stderr, "Can't open: %s\n", filename);
        free_nlj_table(&table);
        return EXIT_FAILURE;
    }

    /* write the quantum numbers out to file */
    fprintf(stream,
            "nlj1\tnlj2\tnlj3\tnnlj1\tnnlj2\tnnlj3\t"
            "J12\tJJ12\ttwoJ\tT12\tTT12\ttwoT\n");
    iterate_me3j(table.nljDim, EMax, table.e_nlj, table.l_nlj, table.twoj_nlj,
                 &main_me3j_iteratee, stream);

    printf("done.\n");
    fclose(stream);
    free_nlj_table(&table);
    return 0;
}
