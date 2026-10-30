/* Determinant of a square integer polynomial matrix modulo each prime of a list (FLINT 3.5.0; certificate
   cyl_stratification_w1014, which compiles this file).
   Usage: detmodp INFILE OUTFILE
   INFILE (whitespace separated): r; then the r*r entries row-major, each "L c_0 ... c_{L-1}" (L = 0 for a zero entry, decimal
          integers of any size, coefficients low -> high); then K and K primes (each checked with n_is_prime).
   OUTFILE: one line per prime, "p D c_0 ... c_D" with det mod p = sum c_k t^k (D = -1 for the zero polynomial).  Primes already in
          OUTFILE are skipped, so a killed run resumes; each line is written in one call and flushed when its prime is done.
   Method: entries reduced mod p (fmpz_poly_get_nmod_poly), then nmod_poly_mat_det_interpolate (evaluation and interpolation over
          Z/pZ; exact mod p because p exceeds the number of evaluation points, which is about the degree bound + 1). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

int main(int argc, char **argv)
{
    if (argc != 3) { fprintf(stderr, "usage: detmodp INFILE OUTFILE\n"); return 2; }
    FILE *in = fopen(argv[1], "r");
    if (!in) { perror("INFILE"); return 2; }
    slong r;
    if (fscanf(in, "%ld", &r) != 1 || r <= 0) { fprintf(stderr, "bad r\n"); return 3; }
    fmpz_poly_mat_t A;
    fmpz_poly_mat_init(A, r, r);
    fmpz_t c;
    fmpz_init(c);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < r; j++)
        {
            slong L;
            if (fscanf(in, "%ld", &L) != 1 || L < 0) { fprintf(stderr, "bad entry length\n"); return 3; }
            for (slong k = 0; k < L; k++)
            {
                if (fmpz_fread(in, c) <= 0) { fprintf(stderr, "bad coefficient\n"); return 3; }
                fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(A, i, j), k, c);
            }
        }
    slong K;
    if (fscanf(in, "%ld", &K) != 1 || K <= 0) { fprintf(stderr, "bad K\n"); return 3; }
    ulong *ps = flint_malloc(K * sizeof(ulong));
    for (slong k = 0; k < K; k++)
    {
        if (fscanf(in, "%lu", &ps[k]) != 1) { fprintf(stderr, "bad prime\n"); return 3; }
        if (!n_is_prime(ps[k])) { fprintf(stderr, "not prime: %lu\n", ps[k]); return 4; }
    }
    fclose(in);

    /* primes already done (complete lines only are expected; the Python driver checks every line) */
    int *done = flint_calloc(K, sizeof(int));
    FILE *old = fopen(argv[2], "r");
    if (old)
    {
        ulong p, x;
        slong D;
        while (fscanf(old, "%lu %ld", &p, &D) == 2)
        {
            for (slong k = 0; k < K; k++) if (ps[k] == p) done[k] = 1;
            for (slong k = 0; k <= D; k++) if (fscanf(old, "%lu", &x) != 1) break;
        }
        fclose(old);
    }
    FILE *out = fopen(argv[2], "a");
    if (!out) { perror("OUTFILE"); return 2; }
    for (slong k = 0; k < K; k++)
    {
        if (done[k]) continue;
        ulong p = ps[k];
        nmod_poly_mat_t Ap;
        nmod_poly_mat_init(Ap, r, r, p);
        for (slong i = 0; i < r; i++)
            for (slong j = 0; j < r; j++)
                fmpz_poly_get_nmod_poly(nmod_poly_mat_entry(Ap, i, j), fmpz_poly_mat_entry(A, i, j));
        nmod_poly_t D;
        nmod_poly_init(D, p);
        nmod_poly_mat_det_interpolate(D, Ap);
        slong deg = nmod_poly_degree(D);
        size_t cap = 64 + 24 * (size_t) (deg + 2);
        char *buf = flint_malloc(cap);
        size_t pos = (size_t) snprintf(buf, cap, "%lu %ld", p, deg);
        for (slong t = 0; t <= deg; t++)
            pos += (size_t) snprintf(buf + pos, cap - pos, " %lu", nmod_poly_get_coeff_ui(D, t));
        pos += (size_t) snprintf(buf + pos, cap - pos, "\n");
        fwrite(buf, 1, pos, out);
        fflush(out);
        flint_free(buf);
        nmod_poly_clear(D);
        nmod_poly_mat_clear(Ap);
    }
    fclose(out);
    flint_free(done);
    flint_free(ps);
    fmpz_clear(c);
    fmpz_poly_mat_clear(A);
    return 0;
}
