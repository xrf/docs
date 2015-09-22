#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int min(int x, int y)
{
    return x < y ? x : y;
}

static int max(int x, int y)
{
    return x > y ? x : y;
}

/** Iterate over `n`, `l`, `j` of all single-particle states.  The function
    `f` is called each iteration with an arbitary argument `env`.  Returns the
    total number of `nlj` values.  */
int iterate_nlj(int eMax, int nMax, int lMax,
                void f(void *env, int nlj, int e, int n, int l, int twoj),
                void *env)
{
    int e, n, l, lMin, twoj, twojMin, twojMax, nlj;

    nlj = 0;

    for (e = 0; e <= eMax; e++) {
        lMin = e % 2;

        for (l = lMin; l <= min(e, lMax); l += 2) {
            n = (e - l) / 2;
            if (n > nMax)
                continue;

            twojMin = abs(2 * l - 1);
            twojMax = 2 * l + 1;

            for (twoj = twojMin; twoj <= twojMax; twoj += 2) {
                (*f)(env, nlj, e, n, l, twoj);
                nlj++;
            }
        }
    }

    return nlj;
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

    /* while looping, we must ensure that both nnlj1 >= nnlj2 and
       (nlj1, nlj2) >= (nnlj1, nnlj2) lexicographically */

    /* *** bra loops *** */
    for (nlj1 = 0; nlj1 < nljDim; nlj1++)
    for (nlj2 = 0; nlj2 <= nlj1; nlj2++) {

        /* ensure that e1 + e2 <= EMax */
        if (e_nlj[nlj1] + e_nlj[nlj2] > EMax)
            break;

        /* *** ket loops *** */
        for (nnlj1 = 0; nnlj1 <= nlj1; nnlj1++)
        for (nnlj2 = 0; nnlj2 <= (nnlj1 == nlj1 ? nlj2 : nnlj1); nnlj2++) {

         /* ensure that ee1 + ee2 <= EMax */
            if (e_nlj[nnlj1] + e_nlj[nnlj2] > EMax)
                break;

            /* ensure parity conservation */
            if ((l_nlj[nlj1] + l_nlj[nlj2] - l_nlj[nnlj1] - l_nlj[nnlj2]) % 2)
                continue;

            /* compute the range of J allowed by triangular condition */
            JMin = max(abs(twoj_nlj[nlj1] - twoj_nlj[nlj2]),
                       abs(twoj_nlj[nnlj1] - twoj_nlj[nnlj2])) / 2;
            JMax = min((twoj_nlj[nlj1] + twoj_nlj[nlj2]),
                       (twoj_nlj[nnlj1] + twoj_nlj[nnlj2])) / 2;

            /* *** diagonal loops *** */

            /* loop is automatically skipped if JMin > JMax */
            for (J = JMin; J <= JMax; J++) {
                for (T = 0; T <= 1; T++) {
                    for (MT = -T; MT <= T; MT++) {

                        /* antisymmetry forbidden states (swap phase
                           for nlj1 == nlj2 must be +1, same for
                           nnlj1 == nnlj2) are not removed in order to
                           have constant-size T-MT-blocks */

                        (*f)(env, nlj1, nlj2, nnlj1, nnlj2, J, T, MT);
                    }
                }
            }
        }
    }
}

struct nlj_table_iteratee_env {
    int *e_nlj, *n_nlj, *l_nlj, *twoj_nlj;
};

static void nlj_table_iteratee(void *env, int nlj,
                               int e, int n, int l, int twoj)
{
    struct nlj_table_iteratee_env *self = (struct nlj_table_iteratee_env *)env;
    self->e_nlj[nlj] = e;
    self->n_nlj[nlj] = n;
    self->l_nlj[nlj] = l;
    self->twoj_nlj[nlj] = twoj;
}

/** Initialize arrays that map `nlj` to `e`, `n`, `l`, `twoj` respectively. */
void nlj_table(int eMax, int nMax, int lMax,
               int *e_nlj, int *n_nlj, int *l_nlj, int *twoj_nlj)
{
    struct nlj_table_iteratee_env self;
    self.e_nlj = e_nlj;
    self.n_nlj = n_nlj;
    self.l_nlj = l_nlj;
    self.twoj_nlj = twoj_nlj;
    iterate_nlj(eMax, nMax, lMax, &nlj_table_iteratee, &self);
}

static void get_nljDim_iteratee(void *env, int nlj,
                                int e, int n, int l, int twoj)
{
}

int get_nljDim(int eMax, int nMax, int lMax)
{
    return iterate_nlj(eMax, nMax, lMax, &get_nljDim_iteratee, NULL);
}

static void get_me2jDim_iteratee(void *env, int nlj1, int nlj2,
                                 int nnlj1, int nnlj2, int J, int T, int MT)
{
    int *counter = (int *)env;
    ++*counter;
}

int get_me2jDim(int nljDim, int EMax, const int *e_nlj, const int *l_nlj, const *twoj_nlj)
{
    int counter = 0;
    iterate_me2j(nljDim, EMax, e_nlj, l_nlj, twoj_nlj,
                 &get_me2jDim_iteratee, &counter);
    return counter;
}

struct nlj_table {
    int nljDim;
    int *e_nlj, *n_nlj, *l_nlj, *twoj_nlj;
};

int init_nlj_table(struct nlj_table *self, int eMax, int nMax, int lMax)
{
    self->nljDim = get_nljDim(eMax, nMax, lMax);
    self->e_nlj = (int *)malloc(sizeof(*self->e_nlj) * self->nljDim * 4);
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
        (double*)malloc(sizeof(*self->matrixElems) * self->me2jDim);
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

static void basic_example_nlj_iteratee(void *env, int nlj,
                                int e, int n, int l, int twoj)
{
    int e2, n2, l2, twoj2;

    printf("%i\t%i\t%i\t%i\t%.1f\n", nlj, e, n, l, twoj / 2.);

    e2 = unpack_nlj(nlj, &n2, &l2, &twoj2);
    if (e != e2 || n != n2 || l != l2 || twoj != twoj2)
        abort();
    if (pack_nlj(e, twoj) != nlj)
        abort();
}

static void basic_example_me2j_iteratee(void *env, int nlj1, int nlj2,
                                        int nnlj1, int nnlj2,
                                        int J, int T, int MT)
{
    FILE *stream = (FILE *)env;
    fprintf(stream, "%i\t%i\t%i\t%i\t%i\t%i\t%i\n",
            nlj1, nlj2, nnlj1, nnlj2, J, T, MT);
}

void basic_example(void)
{
    static const char *filename = "darmstadt-me2j-quantum-numbers.txt";
    struct nlj_table table;
    FILE *stream;

    int eMax = 5;
    int nMax = 99999;
    int lMax = 99999;
    int EMax = 99999;

    if (init_nlj_table(&table, eMax, nMax, lMax)) {
        abort();
    }

    /* print the single-particle states */
    printf("nlj\te\tn\tl\tj\n");
    iterate_nlj(eMax, nMax, lMax, &basic_example_nlj_iteratee, NULL);
    printf("Total of %i single-particle state(s).\n", table.nljDim);

    /* write the quantum numbers out to file */
    printf((
        "\n"
        "Quantum numbers of 2-body matrix elements will be written to:\n"
        "  %s\n"
        "Writing... "
    ), filename);
    fflush(stdout);
    stream = fopen(filename, "w");
    if (!stream) {
        abort();
    }
    fprintf(stream, "nlj1\tnlj2\tnnlj1\tnnlj2\tJ\tT\tMT\n");
    iterate_me2j(table.nljDim, EMax, table.e_nlj, table.l_nlj,
                 table.twoj_nlj, &basic_example_me2j_iteratee, stream);
    printf("done.\n");

    free_nlj_table(&table);
}

int main(void)
{
    basic_example();
    return 0;
}
