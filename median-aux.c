#include "ecm.h"
#include <time.h>
#include <stdlib.h>

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

enum { WANT_TNAIF = 1,
       WANT_TKARA = 2,
       WANT_TT3 = 4
};


void check_result (listz_t a, unsigned int m, listz_t c, unsigned int l,
                  listz_t b, unsigned int n);

#define POLYFORM

void
print_list (listz_t p, unsigned int n)
{
  unsigned int i;

  for (i = 0; i < n; i++)
    {
#ifdef POLYFORM
      if (i > 0 && mpz_cmp_ui (p[i], 0) >= 0)
	printf ("+");
#endif
      mpz_out_str (stdout, 10, p[i]);
#ifdef POLYFORM
      printf ("*x^%u", i);
#else
      printf (" ");
#endif
    }
  putchar ('\n');
}

void TMulNaif (listz_t b, unsigned int n,
	 listz_t a, unsigned int m, listz_t c, unsigned int l)
{
    unsigned int i, j;
    
    for (i = 0; i <= n; i++) {
        mpz_set_ui (b[i], 0);
        for (j = 0; (j <= m) && (i + j <= l); j++) {
            mpz_addmul (b[i], a[j], c[i + j]);
        }
    }
}

int
do_test (int m, int l, int n, int rand_style)
{
  listz_t a, b, c, t;
  int elapsed_time;
  int parano;
  int magic = 42;
  int h;
  unsigned int used_space, forseen_space;

  mpz_t max_rand;
  int i;
  gmp_randstate_t bla;
  a = init_list (m + 42);
  b = init_list (n + 42);
  parano = 20 * MAX(n,m);
  t = init_list (parano);
  c = init_list (l + 42);

  gmp_randinit_default (bla);

  mpz_init (max_rand);
  mpz_set_ui (max_rand, 10);


  if (rand_style == 0) {
  for (i = 0; i <= m; i++)
    {
      mpz_urandomm (a[i], bla, max_rand);
    }

  for (i = 0; i <= l; i++)
    {
      mpz_urandomm (c[i], bla, max_rand);
    }

  for (i = 0; i <= n; i++)
    {
      mpz_urandomm (b[i], bla, max_rand);
    }
  }
  else {
  for (i = 0; i <= m; i++)
    {
      mpz_random (a[i], rand_style);
    }

  for (i = 0; i <= l; i++)
    {
      mpz_random (c[i], rand_style);
    }

  for (i = 0; i <= n; i++)
    {
      mpz_random (b[i], rand_style);
    }
  }

  for (i = 0; i < parano; i++)
      mpz_set_ui (t[i], magic);

  elapsed_time = cputime();
  TToomCookMul (b, n, a, m, c, l, t);
  i = 0;
  while ((mpz_cmp_ui (t[i], magic) != 0) && (i < parano))
      i++;

  elapsed_time = cputime() - elapsed_time;
  printf ("Test effectué avec m = %d, l = %d, n = %d.\n",
          m, l, n);
  printf ("Temps passé dans TToomCookMul: %d.\n", elapsed_time);
  printf ("On a visiblement utilisé jusque %d d'espace temporaire.\n", i);
  used_space = i;

  h = MAX(m,n);

  if (i >= 4 * h)
      printf ("ATTENTION, ON DEBORDE.\n");

  forseen_space = TToomCookMul_space (n, m, l);
  if (i == parano)
      printf ("Oula, on a même attend tout l'espace.\n");
  printf ("A comparer avec %d d'espace prévu.\n", forseen_space);

  if (used_space > forseen_space)
      abort ();

  elapsed_time = cputime();
  check_result (b, n, a, m, c, l);
  elapsed_time = cputime() - elapsed_time;
  printf ("Temps passé dans la vérification naïve: %d.\n", elapsed_time);
  clear_list (a, m + 42);
  clear_list (b, n + 42);
  clear_list (t, parano);
  clear_list (c, l + 42);
  return 0;
}


void check_result (listz_t b, unsigned int n, listz_t a, unsigned int m,
                  listz_t c, unsigned int l)
{
    int i, j;
    mpz_t temoin;
    mpz_init (temoin);
    
    for (i = 0; i <= n; i++) {
        mpz_set_ui (temoin, 0);
        for (j = 0; (j <= m) && (i + j <= l); j++) {
            mpz_addmul (temoin, a[j], c[i + j]);
        }

        if (mpz_cmp (b[i], temoin) != 0) {
            printf ("On s'est vautré à la position %d.\n", i);
            printf ("On attendait: ");
            mpz_out_str (NULL, 10, temoin);
            printf (" et on a eu ");
            mpz_out_str (NULL, 10, b[i]);
            printf ("\n");
            exit (42);
        }
    }
}
        
void show_result (listz_t b, unsigned int n, listz_t a, unsigned int m,
                  listz_t c, unsigned int l)
{
    printf ("RES: TKarmul (");
    print_list (a, m + 1);
    printf (", ");
    print_list (c, l + 1);
    printf (") = ");
    print_list (b, n + 1);
    printf ("ENDRES\n");
    check_result (b, n, a, m, c, l);
}

unsigned int tt_theo (unsigned int k)
{
    unsigned int h;
    if (k <= 5)
        return k;

    h = k / 3 + 1;
    return 4 * tt_theo (h - 1) + tt_theo (k - 2 * h) + k;
}


int do_bench (int k, int size, int bench_mask)
{
  listz_t a, b, c, t;
  int elapsed_time;
  int tmptime;
  int rand_style;
  int nb_runs = 0;
  int msec_threshold = 2000;

  mpz_t max_rand;
  int i;
  gmp_randstate_t bla;
  c = init_list (k * 2);
  b = init_list (k + 1 + 42);
  t = init_list (5 * k + 42);
  a = init_list (k);


  for (i = 0; i < 2 * k - 1; i++)
  {
      mpz_random (c[i], size);
  }

  for (i = 0; i < k; i++)
  {
      mpz_random (a[i], size);
  }

  printf ("k | ");
  if (bench_mask & WANT_TT3)
      printf ("TT3-runs | TT3-time ");

  if (bench_mask & WANT_TKARA)
      printf ("TKara-runs | TKara-time ");

  if (bench_mask & WANT_TNAIF)
      printf ("TNaif-runs | TNaif-time ");

  printf ("\n%d ", k);

  if (bench_mask & WANT_TT3)
  {
    tmptime = 0;
    elapsed_time = cputime ();
    while (tmptime < msec_threshold)
    {
        TToomCookMul (b, k - 1, a, k - 1, c, 2 * k - 2, t);
        nb_runs++;
        tmptime = cputime() - elapsed_time;
    }
    printf ("%d %d %lf ", nb_runs, tmptime, (tmptime + 0.) / (nb_runs+ 0.));
  }

  nb_runs = 0;

  if (bench_mask & WANT_TKARA)
  {
      tmptime = 0;
      elapsed_time = cputime ();
      while (tmptime < msec_threshold)
      {
        TKarMul (b, k - 1, a, k - 1, c, 2 * k - 2, t);
        nb_runs++;
        tmptime = cputime() - elapsed_time;
      }
    printf ("%d %d %lf ", nb_runs, tmptime, (tmptime + 0.) / (nb_runs+ 0.));
  }

  nb_runs = 0;

  if (bench_mask & WANT_TNAIF)
  {
      tmptime = 0;
      elapsed_time = cputime ();
      while (tmptime < msec_threshold)
      {
        check_result (b, k - 1, a, k - 1, c, 2 * k - 2);
        nb_runs++;
        tmptime = cputime() - elapsed_time;
      }
    printf ("%d %d %lf ", nb_runs, tmptime, (tmptime + 0.) / (nb_runs+ 0.));
  }
  printf ("\n");
}


int random_int (int max)
{
    double tmp = (rand() + 1.) / RAND_MAX;
    return tmp * max;
}

int main (int argc, char *argv[])
{
    if (atoi(argv[1]) == 0) {

        do_test (atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
        return 0;
    }
    
    if (atoi(argv[1]) == 1) {
        int nb;
        int nbmax = atoi(argv[2]);
        int perc = 0;
        int lastperc = 0;
        int n, m, l;
        srand (time(0));
        for (nb = 0; nb < nbmax; nb++) {
            m = random_int (10000) + 100;
            n = random_int (10000) + 100;
            l = random_int (10000) + 100;
            do_test (m, l, n,
                     random_int(10) + 1);
            perc = 100 * nb / nbmax;
            if (perc != lastperc) {
                printf ("%d %% effectue.\n", perc);
                lastperc = perc;
            }
        }
        return 0;
    }
    if (atoi(argv[1]) == 2) {
        do_bench (atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
        return 0;
    }
    printf ("usage:\nmedian 0 m l n size\nmedian 1 nbtests\nmedian 2 k nblimbs \n");
}
