#include <assert.h>
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

/** Find the common range of J formed by coupling j1 with j2 and coupling jj1
    with jj2, where j1 and j2 must be either both half-integers or both
    integers, and similarly for jj1 and jj2. */
void range_J_both(int twoj1, int twoj2, int twojj1, int twojj2,
                  int *JMin, int *JMax)
{
    range_twoJ_both(twoj1, twoj2, twojj1, twojj2, JMin, JMax);
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

int shell_number(int n, int l)
{
    return 2 * n + l;
}

int pack_nlj(int e, int twoj)
{
    return (e * (e + 1) + twoj) / 2;
}

int unpack_nlj(int nlj, int *n, int *l, int *twoj)
{
    const int e = ((int)sqrt(8 * nlj + 1) - 1) / 2;
    const int erem = nlj - e * (e + 1) / 2;
    *twoj = 1 + erem * 2;
    *n = (e - erem) / 2;
    *l = e - 2 * *n;
    return e;
}

int violates_parity_2(int nlj1, int nlj2, int nnlj1, int nnlj2,
                      const int *l_nlj)
{
    return (l_nlj[nlj1] + l_nlj[nlj2] - l_nlj[nnlj1] - l_nlj[nnlj2]) % 2;
}

/** Iterate over all 2-body matrix elements in the Darmstadt ME2J order.  The
    function `f` is called each iteration with an arbitary argument `env`. */
void iterate_me2j(int nljDim, int EMax,
                  const int *e_nlj, const int *l_nlj, const int *twoj_nlj,
                  void (*f)(void *env, int nlj1, int nlj2, int nnlj1, int nnlj2,
                            int J, int T, int MT),
                  void *env)
{
    int nlj1, nlj2, nnlj1, nnlj2, J, JMin, JMax, T, MT;
    for (nlj1 = 0; nlj1 < nljDim; nlj1++)
    for (nlj2 = 0; nlj2 <= nlj1; nlj2++) {
        if (e_nlj[nlj1] + e_nlj[nlj2] > EMax)
            break;
        for (nnlj1 = 0; nnlj1 <= nlj1; nnlj1++)
        for (nnlj2 = 0; nnlj2 <= (nnlj1 == nlj1 ? nlj2 : nnlj1); nnlj2++) {
            if (e_nlj[nnlj1] + e_nlj[nnlj2] > EMax)
                break;
            if (violates_parity_2(nlj1, nlj2, nnlj1, nnlj2, l_nlj))
                continue;
            range_J_both(twoj_nlj[nlj1], twoj_nlj[nlj2],
                         twoj_nlj[nnlj1], twoj_nlj[nnlj2],
                         &JMin, &JMax);
            for (J = JMin; J <= JMax; J++)
            for (T = 0; T <= 1; T++)
            for (MT = -T; MT <= T; MT++) {
                (*f)(env, nlj1, nlj2, nnlj1, nnlj2, J, T, MT);
            }
        }
    }
}

static void get_me2jDim_iteratee(void *env, int nlj1, int nlj2,
                                 int nnlj1, int nnlj2, int J, int T, int MT)
{
    int *counter = (int *)env;
    ++*counter;
}

int get_me2jDim(int nljDim, int EMax, const int *e_nlj,
                const int *l_nlj, const int *twoj_nlj)
{
    int counter = 0;
    iterate_me2j(nljDim, EMax, e_nlj, l_nlj, twoj_nlj,
                 &get_me2jDim_iteratee, &counter);
    return counter;
}

/** Read an .me2j matrix element text file into the given array. */
int read_matrixElems(const char *filename, int me2jDim, double *matrixElems)
{
    char line[1024];
    int i;

    FILE *stream = fopen(filename, "r");
    if (!stream) {
        return 1;
    }

    /* skip one line */
    if (!fgets(line, sizeof(line), stream) ||
        strlen(line) == sizeof(line) - 1) {
        goto invalid_format;
    }

    for (i = 0; i < me2jDim; i++) {
        double num;
        const int status = fscanf(stream, "%lf", &num);
        if (status != 1) {
            if (ferror(stream)) {
                goto invalid_format;
            }
            break;
        }
        matrixElems[i] = num;
    }

    fclose(stream);

    if (i < me2jDim) {
        return 3;
    }

    return 0;

invalid_format:
    fclose(stream);
    return 2;
}

struct me2j {
    int me2jDim;
    double *matrixElems;
    struct nlj_table table;
};

void free_me2j(struct me2j *self)
{
    free_nlj_table(&self->table);
    free(self->matrixElems);
}

/** Read an .me2j matrix element text file. */
int read_me2j(struct me2j *self, const char *filename,
              int eMax, int nMax, int lMax, int EMax)
{
    int e;

    if (init_nlj_table(&self->table, eMax, nMax, lMax)) {
        return -1;
    }

    self->me2jDim = get_me2jDim(self->table.nljDim, EMax, self->table.e_nlj,
                                self->table.l_nlj, self->table.twoj_nlj);
    self->matrixElems =
        (double*)malloc(sizeof(*self->matrixElems) * (size_t)self->me2jDim);
    if (!self->matrixElems) {
        free_nlj_table(&self->table);
        return -1;
    }

    e = read_matrixElems(filename, self->me2jDim, self->matrixElems);
    if (e) {
        free_nlj_table(&self->table);
        free(self->matrixElems);
        return e;
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Example program */
/* ------------------------------------------------------------------------ */

static void main_nlj_iteratee(void *env, int e, int n, int l, int twoj)
{
    int e2, n2, l2, twoj2, *nlj = (int *)env;
    printf("%i\t%i\t%i\t%i\t%.1f\n", *nlj, e, n, l, twoj / 2.);

    /* make sure pack_nlj and unpack_nlj are inverses */
    e2 = unpack_nlj(*nlj, &n2, &l2, &twoj2);
    assert(pack_nlj(e2, twoj2) == *nlj);

    ++*nlj;
}

static void main_me2j_iteratee(void *env,
                               int nlj1, int nlj2,
                               int nnlj1, int nnlj2,
                               int J, int T, int MT)
{
    FILE *stream = (FILE *)env;
    fprintf(stream, "%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
            nlj1, nlj2, nnlj1, nnlj2, J, T, MT);
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
    *filename = "darmstadt_me2j_quantum_numbers.txt";
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
    int nlj, eMax, nMax, lMax, EMax;
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

    /* print the single-particle states */
    printf("nlj\te\tn\tl\tj\n");
    nlj = 0;
    iterate_nlj(eMax, nMax, lMax, &main_nlj_iteratee, &nlj);
    printf("Total of %i single-particle state(s).\n\n", table.nljDim);

    printf("Quantum numbers of 2-body matrix elements will be written to:\n"
           "  %s\nWriting... ", filename);
    fflush(stdout);

    stream = fopen(filename, "w");
    if (!stream) {
        fprintf(stderr, "Can't open: %s\n", filename);
        free_nlj_table(&table);
        return EXIT_FAILURE;
    }

    /* write the quantum numbers out to file */
    fprintf(stream, "nlj1\tnlj2\tnnlj1\tnnlj2\tJ\tT\tMT\n");
    iterate_me2j(table.nljDim, EMax, table.e_nlj, table.l_nlj, table.twoj_nlj,
                 &main_me2j_iteratee, stream);

    printf("done.\n");
    fclose(stream);
    free_nlj_table(&table);
    return 0;
}
