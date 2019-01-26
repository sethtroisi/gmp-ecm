/* Simple expression parser for GMP-ECM.

Copyright 2003, 2004, 2005, 2006, 2007, 2008, 2012 Jim Fougeron, Paul Zimmermann and Alexander Kruppa.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ecm-ecm.h"
#include "getprime_r.h"

#ifdef HAVE_STRINGS_H
# include <strings.h> /* for strncasecmp */
#endif

#ifdef HAVE_CTYPE_H
# include <ctype.h>
#endif


/*****************************************************************
 *   Syntax for this VERY simple recursive expression parser:	 *
 *                                                               *
 *   ( or [ or { along with ) or ] or } are valid for grouping   *
 *   Normal "simple" operators:  + - * /  (. can be used for *)  *
 *   Module:                     n%m    345%11                   *
 *   Unary minus is supported:   -n     -500                     *
 *   Exponentation:              n^m    2^500                    *
 *   Simple factorial:           n!     53!  == 1*2*3*4...*52*53 *
 *   Multi-factorial:            n!m    15!3 == 15.12.9.6.3      *
 *   Simple Primorial:           n#     11# == 2*3*5*7*11        *
 *   Reduced Primorial:          n#m    17#5 == 5.7.11.13.17     *
 *                                                               *
 * Adding (working on these at least:				 *
 *   Phi(x,n)							 *
 *   U(p,q,n)							 *
 *   primU(p,q,n)						 *
 *                                                               *
 * NOTE Lines ending in a \ character are "joined"               *
 * NOTE C++ // single line comments (rest of line is a comment)  *
 *                                                               *
 ****************************************************************/

/* value only used by the expression parser */
static mpz_t t, mpOne;
static char *expr_str;

static void eval_power (mpz_t prior_n, mpz_t n,char op);
static void eval_product (mpz_t prior_n, mpz_t n,char op);
static void eval_sum (mpz_t prior_n, mpz_t n,char op);
static int  eval_Phi (mpz_t *params, mpz_t n);
static int  eval_U (mpz_t *params, mpz_t n);
static int  eval_primU (mpz_t *params, mpz_t n);
static int  eval_2 (int bInFuncParams);

#if 0 /* strncasecmp is a required function in configure.in */
#if defined (_MSC_VER) || defined (__MINGW32__)
#define strncasecmp strnicmp
#endif
#endif

/***************************************/
/* Main expression evaluation function */
/* This is the function that the app   */
/* calls to read the expression line   */
/***************************************/
int eval (mpcandi_t *n, FILE *fd, int primetest)
{
  int ret;
  int nMaxSize = 2000, nCurSize = 0;
  int c;
  char *expr = (char *) malloc (nMaxSize + 1);

  ASSERT_ALWAYS (expr != NULL);
JoinLinesLoop:
  c = fgetc (fd);
  if (0)
    {
ChompLine:
      do
        c = fgetc (fd);
      while (c != EOF && !IS_NEWLINE(c));
      if (IS_NEWLINE(c))
        goto JoinLinesLoop;
    }
  
  while (c != EOF && !IS_NEWLINE(c) && c != ';')
    {
      if (c == '/')
	{
	  /* This might be a C++ // comment or it might be a / division operator.  
	     Check it out, and if it is a comment, then "eat it" */
	  int peek_c = fgetc (fd);
	  if (peek_c == '/')
	    /* Got a C++ single line comment, so Chomp the line */
	    goto ChompLine;

	  /* Put the char back on the file, then allow the code to add the '/' char to the buffer */
	  ungetc (peek_c, fd);
        }

      /* strip space and tabs out here, and then we DON'T have to mess with them in the rest of the parser */
      if (!isspace (c) && c != '"' && c != '\'')
	expr[nCurSize++] = (char) c;

      if (nCurSize == nMaxSize)
      {
	char *cp;
	nMaxSize += nMaxSize / 2;
	cp = (char *) realloc (expr, nMaxSize + 1);
        ASSERT_ALWAYS (cp != NULL);
	expr = cp;
      }
      c = fgetc (fd);
    }
  expr[nCurSize] = 0;
  if (!nCurSize)
    ret = 0;
  else
    {
      if (c == ';')
	ungetc (c, fd);
      mpz_init (t);
      expr_str = expr;
      ret = eval_2 (0);
      if (ret)
	{
	  char *s;
	  char *cpTmpExpr = expr;
	  s = mpz_get_str (NULL, 10, t);
	  if (!strcmp(s, cpTmpExpr))
	    cpTmpExpr = NULL;
	  ret = mpcandi_t_add_candidate (n, t, cpTmpExpr, primetest);
	  free (s); /* size strlen (s) + 1 */
	}
      mpz_clear(t);
    }
  free(expr);
  return ret;
}

int eval_str (mpcandi_t *n, char *cp, int primetest, char **EndChar)
{
  int ret;
  int nMaxSize=2000, nCurSize=0;
  char *c;
  char *expr = (char *) malloc(nMaxSize+1);

  ASSERT_ALWAYS (expr != NULL);
  c = cp;
JoinLinesLoop:
  if (*c == '#')
    {
      do
        ++c;
      while (*c && !IS_NEWLINE(*c));
      if (IS_NEWLINE(*c))
        goto JoinLinesLoop;
    }
  
  while (*c && !IS_NEWLINE(*c) && *c != ';')
    {
      /* strip space and tabs out here, and then we DON'T have to mess with them in the rest of the parser */
      if (!isspace((int) *c) && *c != '"' && *c != '\'')
	expr[nCurSize++] = *c;
      if (nCurSize == nMaxSize)
      {
	char *cp;
	nMaxSize += nMaxSize / 2;
	cp = (char *) realloc (expr, nMaxSize + 1);
        ASSERT_ALWAYS (cp != NULL);
	expr = cp;
      }
      ++c;
    }
  expr[nCurSize] = 0;
  if (!nCurSize)
    ret = 0;
  else
    {
      if (*c != ';')
	++c;
      mpz_init(t);
      expr_str = expr;
      ret = eval_2(0);
      if (ret)
	{
	  char *s;
	  char *cpTmpExpr = expr;
	  s = mpz_get_str (NULL, 10, t);
	  if (!strcmp(s, cpTmpExpr))
	    cpTmpExpr = NULL;
	  ret = mpcandi_t_add_candidate(n, t, cpTmpExpr, primetest);
	  free (s); /* size strlen (s) + 1 */
	}
      mpz_clear(t);
    }
  free(expr);
  if (EndChar && *EndChar)
    *EndChar = c;
  return ret;
}

void eval_power (mpz_t prior_n, mpz_t n,char op)
{
#if defined (DEBUG_EVALUATOR)
  if ('#'==op || '^'==op || '!'==op || '@'==op || '$'==op)
    {
      fprintf (stderr, "eval_power ");
      mpz_out_str(stderr, 10, prior_n);
      fprintf (stderr, "%c", op);
      mpz_out_str(stderr, 10, n);
      fprintf (stderr, "\n");
    }
#endif

  if ('^'==op)
    mpz_pow_ui(n,prior_n,mpz_get_ui(n));
  else if ('!'==op)	/* simple factorial  (syntax n!    example: 7! == 1*2*3*4*5*6*7) */
    mpz_fac_ui(n,mpz_get_ui(n));
  else if ('@'==op)	/* Multi factorial   (syntax n!prior_n
                           Example: 15!3 == 15*12*9*6*3
                           Note: 15!3 is substituted into 15@3 by the parser */
    {
      long nCur;
      unsigned long nDecr;
      nCur = mpz_get_si(prior_n);
      nDecr = mpz_get_ui(n);
      mpz_set_ui(n,1);
      while (nCur > 1)
	{
	  /* This could be done much more efficiently (bunching mults using smaller "built-ins"), but I am not going to bother for now */
	  mpz_mul_ui(n,n,nCur);
	  nCur -= nDecr;
	}
    }
  else if ('#'==op)  /* simple primorial (syntax  n#   example:  11# == 2*3*5*7*11 */
    {
      unsigned long nMax;
      unsigned long p;
      prime_info_t prime_info;

      prime_info_init (prime_info);
      nMax = mpz_get_ui (n);
      mpz_set_ui (n, 1);
      for (p = 2; p <= nMax; p = getprime_mt (prime_info))
	/* This could be done much more efficiently (bunching mults using smaller "built-ins"), but I am not going to bother for now */
	mpz_mul_ui (n, n, p);
      prime_info_clear (prime_info); /* free the prime table */
    }
  else if ('$'==op)  /* reduced primorial (syntax  n#prior_n   example:  13#5 == (5*7*11*13) */
    {
      unsigned long p;
      unsigned long nMax;
      unsigned long nStart;
      prime_info_t prime_info;

      prime_info_init (prime_info);
      nMax = mpz_get_ui (prior_n);
      nStart = mpz_get_ui (n);
      mpz_set_ui (n, 1);
      /*printf ("Reduced-primorial  %ld#%ld\n", nMax, nStart);*/
      for (p = 2; p <= nMax; p = getprime_mt (prime_info))
	{
	  if (p >= nStart)
	    /* This could be done much more efficiently (bunching mults using smaller "built-ins"), but I am not going to bother for now */
	    mpz_mul_ui (n, n, p);
	}
      prime_info_clear (prime_info); /* free the prime table */
    }
}

void
eval_product (mpz_t prior_n, mpz_t n, char op)
{
#if defined (DEBUG_EVALUATOR)
  if ('*'==op || '.'==op || '/'==op || '%'==op)
    {
      fprintf (stderr, "eval_product ");
      mpz_out_str(stderr, 10, prior_n);
      fprintf (stderr, "%c", op);
      mpz_out_str(stderr, 10, n);
      fprintf (stderr, "\n");
    }
#endif
  if ('*' == op || '.' == op)
    mpz_mul (n, prior_n, n);
  else if ('/' == op)
    {
      mpz_t r;
      mpz_init (r);
      mpz_tdiv_qr (n, r, prior_n, n);
      if (mpz_cmp_ui (r, 0) != 0)
        {
          fprintf (stderr, "Parsing Error: inexact division\n");
          exit (EXIT_FAILURE);
        }
      mpz_clear (r);
    }
  else if ('%' == op)
    mpz_tdiv_r (n, prior_n, n);
}

void eval_sum (mpz_t prior_n, mpz_t n,char op)
{
#if defined (DEBUG_EVALUATOR)
  if ('+'==op || '-'==op)
    {
      fprintf (stderr, "eval_sum ");
      mpz_out_str(stderr, 10, prior_n);
      fprintf (stderr, "%c", op);
      mpz_out_str(stderr, 10, n);
      fprintf (stderr, "\n");
    }
#endif

  if ('+' == op)
    mpz_add(n,prior_n,n);
  else if ('-' == op)
    mpz_sub(n,prior_n,n);
}

int eval_Phi (mpz_t* params, mpz_t n)
{
  int factors[200];
  unsigned dwFactors=0, dw;
  unsigned long B;
  unsigned long p;
  mpz_t D, T, org_n;
  prime_info_t prime_info;
  
  if (mpz_cmp_ui (n, 1) == 0)
    {
      /* return value is 1 if b is composite, or b if b is prime */
      int isPrime = mpz_probab_prime_p (params[0], PROBAB_PRIME_TESTS);
      if (isPrime)
	mpz_set (n, params[0]);
      else
	mpz_set (n, mpOne);
      return 1;
    }
  if (mpz_cmp_si (n, -1) == 0)
    /* this is actually INVALID, but it is easier to simply */
    return 0;

  /* OK parse the Phi out now */
  if (mpz_cmp_ui (params[0], 0) <= 0)
    return 0;

  if (mpz_cmp_ui (params[0], 1) == 0)
    {
      if (mpz_cmp_ui (n, 1) != 0)
	mpz_sub_ui (n, n, 1);
      return 1;
    }

  /* Ok, do the real h_primative work, since we are not one of the trivial case */

  if (mpz_fits_ulong_p (params[0]) == 0)
    return 0;

  B = mpz_get_ui (params[0]);

  /* Obtain the factors of B */
  prime_info_init (prime_info);
  for (p = 2; p <= B; p = getprime_mt (prime_info))
    {
      if (B % p == 0)
	{
	  /* Add the factor one time */
	  factors[dwFactors++] = p;
	  /* but be sure to totally remove it */
	  do { B /= p; } while (B % p == 0);
        }
     }
  prime_info_clear (prime_info); /* free the prime tables */
  B = mpz_get_si (params[0]);

  mpz_init_set (org_n, n);
  mpz_set_ui (n, 1);
  mpz_init_set_ui (D, 1);
  mpz_init (T);
	      
  for(dw=0;(dw<(1U<<dwFactors)); dw++)
    {
      /* for all Mobius terms */
      int iPower=B;
      int iMobius=0;
      unsigned dwIndex=0;
      unsigned dwMask=1;
		  
      while(dwIndex < dwFactors)
	{
	  if(dw&dwMask)
	    {
	      /* printf ("iMobius = %d iPower = %d, dwIndex = %d ", iMobius, iPower, dwIndex); */
	      iMobius++;
	      iPower/=factors[dwIndex];
	      /* printf ("Then iPower = %d\n", iPower);  */
	    }
	  dwMask<<=1;
	  ++dwIndex;
	}
      /*
      fprintf (stderr, "Taking ");
      mpz_out_str(stderr, 10, org_n);
      fprintf (stderr, "^%d-1\n", iPower);
      */
      mpz_pow_ui(T, org_n, iPower);
      mpz_sub_ui(T, T, 1);
    
      if(iMobius&1)
      {
	/*
	fprintf (stderr, "Muling D=D*T  ");
	mpz_out_str(stderr, 10, D);
	fprintf (stderr, "*");
	mpz_out_str(stderr, 10, T);
	fprintf (stderr, "\n");
	*/
	mpz_mul(D, D, T);
      }
      else
      {
	/*
	fprintf (stderr, "Muling n=n*T  ");
	mpz_out_str(stderr, 10, n);
	fprintf (stderr, "*");
	mpz_out_str(stderr, 10, T);
	fprintf (stderr, "\n");
	*/
	mpz_mul(n, n, T); 
      }
  }
  mpz_divexact(n, n, D);
  mpz_clear(T);
  mpz_clear(org_n);
  mpz_clear(D);

  return 1;
}

int eval_U (mpz_t *params, mpz_t n)
/* params[0]=P, params[1]=Q */
{
  unsigned long N;
  mpz_t U1,U0,org_n,D,T; /* At each step U1 holds U(k), and U0 holds U(k-1) */
  long k,l;
  
  if (mpz_cmp_si (n, 0) < 0)
    return 0;
  if (mpz_cmp_ui (n, 1) == 0)
    {
      mpz_set_ui (n, 1);
      return 1;
    }
  if (mpz_cmp_ui (n, 0) == 0)
    {
      mpz_set_ui (n, 0);
      return 1;
    }
  if (mpz_fits_ulong_p (n) == 0)
    return 0;

  N = mpz_get_ui (n);
  if (mpz_cmp_ui (params[0], 0) == 0)
    {
      if( N%2==0 )
        {
          mpz_set_ui (n, 0);
	}
      else
        {
	  mpz_neg (params[1], params[1]);
	  mpz_pow_ui (n, params[1], (N-1)/2);
	  mpz_neg (params[1], params[1]);
	}
      return 1;
    }


  mpz_init_set (org_n, n);
  mpz_init_set_ui (U1, 1);
  mpz_init_set_ui (U0, 0);
  mpz_init (D);
  mpz_init (T);
  mpz_mul (D, params[0], params[0]);
  mpz_submul_ui (D, params[1], 4);
  k=1;

  for(l=mpz_sizeinbase(org_n,2)-2;l>=0;l--)
    {
       mpz_mul (U0, U0, U0);
       mpz_mul (U1, U1, U1);
       mpz_mul (U0, U0, params[1]);
       mpz_sub (U0, U1, U0); // U(2k-1)=U(k)^2-QU(k-1)^2
       mpz_pow_ui (T, params[1], k);
       mpz_mul (U1, U1, D);
       mpz_addmul_ui (U1, T, 2);
       mpz_addmul (U1, params[1], U0); // U(2k+1)=DU(k)^2+2Q^k+QU(2k-1)
       if (mpz_tstbit (org_n, l) )
         {
	    k=2*k+1;
	    mpz_mul (U0,U0,params[1]); // U0 is 2k, U1 is 2k+1
	    mpz_add (U0,U1,U0);
	    mpz_divexact (U0,U0,params[0]);
	 }
       else
         {
 	    k=2*k;
	    mpz_addmul (U1,U0,params[1]); // U0 is 2k-1, U1 is 2k
	    mpz_divexact (U1,U1,params[0]);
         }
	 /* gmp_printf("%d %Zd %Zd\n",k,U0,U1); */
    }
  mpz_set(n, U1);

  mpz_clear(U0);
  mpz_clear(U1);
  mpz_clear(org_n);
  mpz_clear(D);
  mpz_clear(T);

  return 1;
}

int eval_primU (mpz_t* params, mpz_t n)
{
  int factors[200];
  unsigned dwFactors=0, dw;
  unsigned long N;
  unsigned long p;
  mpz_t D, T;
  
  if (mpz_cmp_ui (n, 0) <= 0)
    return 0;
  if (mpz_cmp_ui (n, 1) == 0)
    {
	mpz_set_ui (n, 1);
      return 1;
    }

  /* Ignore the special cases where P^2=0,Q or 4Q*/
  if (mpz_cmp_ui (params[0], 0) == 0)
    {
      return 0;
    }
  mpz_init(D);
  mpz_mul(D, params[0], params[0]);
  if (mpz_cmp (D, params[1]) == 0)
    {
      return 0;
    }
  mpz_submul_ui(D, params[1], 4);
  if (mpz_cmp_ui (D, 0) == 0)
    {
      return 0;
    }  


  if (mpz_fits_ulong_p (n) == 0)
    return 0;

  N = mpz_get_ui (n);

   /* Obtain the factors of N */
  for (p = 2; p <= N; p++)
    {
      if (N % p == 0)
	{
	  /* Add the factor one time */
	  factors[dwFactors++] = p;
	  /* but be sure to totally remove it */
	  do { N /= p; } while (N % p == 0);
        }
     }
  
  
  N = mpz_get_ui (n);

  mpz_set_ui (n, 1);
  mpz_set_ui (D, 1);
  mpz_init (T);
	      
  for(dw=0;(dw<(1U<<dwFactors)); dw++)
    {
      /* for all Mobius terms */
      int iPower=N;
      int iMobius=0;
      unsigned dwIndex=0;
      unsigned dwMask=1;
		  
      while(dwIndex < dwFactors)
	{
	  if(dw&dwMask)
	    {
	      /* printf ("iMobius = %d iPower = %d, dwIndex = %d ", iMobius, iPower, dwIndex); */
	      iMobius++;
	      iPower/=factors[dwIndex];
	      /* printf ("Then iPower = %d\n", iPower);  */
	    }
	  dwMask<<=1;
	  ++dwIndex;
	}
 /*     
      gmp_fprintf (stderr, "Taking U(%Zd,%Zd,%d)\n", P,Q,iPower);
 */     
      mpz_set_ui(T,iPower);
      eval_U(params, T);
    
      if(iMobius&1)
      {
/*	
	fprintf (stderr, "Muling D=D*T  ");
	mpz_out_str(stderr, 10, D);
	fprintf (stderr, "*");
	mpz_out_str(stderr, 10, T);
	fprintf (stderr, "\n");
*/	
	mpz_mul(D, D, T);
      }
      else
      {
/*	
	fprintf (stderr, "Muling n=n*T  ");
	mpz_out_str(stderr, 10, n);
	fprintf (stderr, "*");
	mpz_out_str(stderr, 10, T);
	fprintf (stderr, "\n");
*/	
	mpz_mul(n, n, T); 
      }
  }
  mpz_divexact(n, n, D);
  mpz_clear(T);
  mpz_clear(D);

  return 1;
}

/* A simple partial-recursive decent parser */
int eval_2 (int bInFuncParams)
{
  mpz_t n_stack[5];
  mpz_t n;
  mpz_t param_stack[5];
  int i,j;
  int num_base;
  char op_stack[5];
  char op;
  char negate;
  const int num_of_funcs=3;
  const char *func_names[]={"Phi","U","primU"};
  const int func_num_params[]={2,3,3};
  typedef int (*fptr)(mpz_t *,mpz_t);
  const fptr func_ptrs[]={eval_Phi,eval_U,eval_primU};
  char *paren_position;
  char tentative_func_name[20];
  int func_id;
  for (i=0;i<5;i++)
    {
      op_stack[i]=0;
      mpz_init(n_stack[i]); 
      mpz_init(param_stack[i]); 
    }
  mpz_init(n);
  op = 0;
  negate = 0;

  for (;;)
    {
      if ('-'==(*expr_str))
	{
	  expr_str++;
	  negate=1;
	}
      else
	negate=0;
      if ('('==(*expr_str) || '['==(*expr_str) || '{'==(*expr_str))    /* Number is subexpression */
	{
	  expr_str++;
	  eval_2 (bInFuncParams);
	  mpz_set(n, t);
	}
      else            /* Number is decimal value */
	{
	  for (i=0;isdigit((int) expr_str[i]);i++)
	    ;
	  if (!i)         /* No digits found */
	    {
	      /* check for a valid "function" */
	      paren_position=strchr(&expr_str[i], '(');
	      if (NULL==paren_position)
	        {
		  /* No parentheses found */
		  fprintf (stderr, "\nError - invalid number [%c]\n", expr_str[i]);
		  exit (EXIT_FAILURE);
		}
	      strncpy(tentative_func_name,&expr_str[i],paren_position-&expr_str[i]);
	      tentative_func_name[paren_position-&expr_str[i]]='\0';
	      for (func_id=0;func_id<num_of_funcs;func_id++)
	        {
		  if (!strcasecmp (tentative_func_name, func_names[func_id]))
		    break;
		}
	      if(func_id==num_of_funcs)	/* No matching function name found */
	        {
		  fprintf (stderr, "Error, Unknown function %s()\n", tentative_func_name);
		  exit (EXIT_FAILURE);
		}
	      /* Now we can actually process existing functions */
	      expr_str=paren_position+1;
	      for(j=0;j<func_num_params[func_id]-1;j++)
	        {
		  /* eval the first parameters.  NOTE we pass a 1 since we ARE in parameter mode, 
		     and this causes the ',' character to act as the end of expression */		
		  if(eval_2 (1) != 2)
		    {
		      fprintf (stderr, "Error, Function %s() requires %d parameters\n", func_names[func_id], func_num_params[func_id]);
                      exit (EXIT_FAILURE);
		    }
		  mpz_set(param_stack[j], t);
		}
	      /* Now eval the last parameter NOTE we pass a 0 since we are NOT expecting a ','
		 character to end the expression, but are expecting a ) character to end the function */	
	      if (eval_2 (0))
		{
		  mpz_set(n, t);
		  if( (func_ptrs[func_id])(param_stack, n) == 0 )
		    {
		      fprintf (stderr, "\nParsing Error -  Invalid "
                               "parameter passed to the %s function\n", func_names[func_id]);
		      exit (EXIT_FAILURE);
		    }
		}
	      goto MONADIC_SUFFIX_LOOP;
	    }
	  /* Now check for a hex number.  If so, handle it as such */
	  num_base=10;  /* assume base 10 */
	  if (i == 1 && !strncasecmp(expr_str, "0x", 2))
	    {
		num_base = 16;	/* Kick up to hex */
		expr_str += 2;	/* skip the 0x string of the number */
		for (i=0;isxdigit((int) expr_str[i]);i++)
	          ;
	    }
	  op=expr_str[i];
	  expr_str[i]=0;
	  mpz_set_str(n,expr_str,num_base);
	  expr_str+=i;
	  *expr_str=op;
      }
      if (negate) 
	mpz_neg(n,n);

/* This label is needed for "normal" primorials and factorials, since they are evaluated Monadic suffix 
   expressions.  Most of this parser assumes Dyadic operators  with the only exceptino being the
   Monadic prefix operator of "unary minus" which is handled by simply ignoring it (but remembering),
   and then fixing the expression up when completed. */
/* This is ALSO where functions should be sent.  A function should "act" like a stand alone number.
   We should NOT start processing, and expecting a number, but we should expect an operator first */
MONADIC_SUFFIX_LOOP:
        op=*expr_str++;
	    
      if (0==op || ')'==op || ']'==op || '}'==op || (','==op&&bInFuncParams))
	{
	  eval_power (n_stack[2],n,op_stack[2]);
	  eval_product (n_stack[1],n,op_stack[1]);
	  eval_sum (n_stack[0],n,op_stack[0]);
	  mpz_set(t, n);
	  mpz_clear(n);
	  for (i=0;i<5;i++) 
	    mpz_clear(n_stack[i]);
	  /* Hurray! a valid expression (or sub-expression) was parsed! */
	  return ','==op?2:1;
	}
      else
	{
	  if ('^' == op)
	    {
	      eval_power (n_stack[2],n,op_stack[2]);
	      mpz_set(n_stack[2], n);
	      op_stack[2]='^';
	    }
	  else if ('!' == op)
	    {
	      if (!isdigit((int) *expr_str))
		{
		  /* If the next char is not a digit, then this is a simple factorial, and not a "multi" factorial */
		  mpz_set(n_stack[2], n);
		  op_stack[2]='!';
		  goto MONADIC_SUFFIX_LOOP;
		}
	      eval_power (n_stack[2],n,op_stack[2]);
	      mpz_set(n_stack[2], n);
	      op_stack[2]='@';
	    }
	  else if ('#' == op)
	    {
	      if (!isdigit((int) *expr_str))
		{
  		  /* If the next char is not a digit, then this is a simple primorial, and not a "reduced" primorial */
		  mpz_set(n_stack[2], n);
		  op_stack[2]='#';
		  goto MONADIC_SUFFIX_LOOP;
		}
	      eval_power (n_stack[2],n,op_stack[2]);
	      mpz_set(n_stack[2], n);
	      op_stack[2]='$';
	    }
	  else
	    {
	      if ('.'==op || '*'==op || '/'==op || '%'==op)
		{
		  eval_power (n_stack[2],n,op_stack[2]);
		  op_stack[2]=0;
		  eval_product (n_stack[1],n,op_stack[1]);
		  mpz_set(n_stack[1], n);
		  op_stack[1]=op;
		}
	      else
		{
		  if ('+'==op || '-'==op)
		    {
		      eval_power (n_stack[2],n,op_stack[2]);
		      op_stack[2]=0;
		      eval_product (n_stack[1],n,op_stack[1]);
		      op_stack[1]=0;
		      eval_sum (n_stack[0],n,op_stack[0]);
		      mpz_set(n_stack[0], n);
		      op_stack[0]=op;
		    }
		  else    /* Error - invalid operator */
		    {
		      fprintf (stderr, "\nError - unknown operator: '%c'\n", op);
                      exit (EXIT_FAILURE);
		     }
		}
	    }
	}
    }
}

void init_expr(void)
{
  mpz_init_set_ui(mpOne, 1);
}

void free_expr(void)
{
  mpz_clear(mpOne);
}
