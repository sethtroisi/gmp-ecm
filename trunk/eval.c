/* Simple expression parser for GMP-ECM.

  Copyright 2003, 2004, 2005 Jim Fougeron, Paul Zimmermann and Alexander Kruppa.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* for strncasecmp */
#if !defined (_MSC_VER)
#include <unistd.h>
#endif
#include <ctype.h>
#include <time.h>
#include "gmp.h"
#include "ecm.h"
#include "ecm-ecm.h"

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
 *                                                               *
 * NOTE Lines ending in a \ character are "joined"               *
 * NOTE Lines starting with #are comments                        *
 * NOTE C++ // single line comments (rest of line is a comment)  *
 *                                                               *
 ****************************************************************/

/* value only used by the expression parser */
static mpz_t t, mpOne;
static char *expr_str;

static void eval_power (mpz_t prior_n, mpz_t n,char op);
static void eval_product (mpz_t prior_n, mpz_t n,char op);
static void eval_sum (mpz_t prior_n, mpz_t n,char op);
static int  eval_Phi (mpz_t prior_n, mpz_t n, int ParamCnt);
static int  eval_2 (int bInFuncParams);

#if defined (_MSC_VER) || defined (__MINGW32__)
#define strncasecmp strnicmp
#endif

/**************************************/
/* Main expression evalation function */
/* This is the function that the app  */
/* calls to read the expression line  */
/**************************************/
int eval (mpcandi_t *n, FILE *fd, int primetest)
{
  int ret;
  int nMaxSize = 2000, nCurSize = 0;
  int c;
  char *expr = (char *) malloc (nMaxSize + 1);

  if (expr == NULL)
    {
      fprintf (stderr, "Error: not enough memory\n");
      exit (EXIT_FAILURE);
    }
/* Lines ending in '\\' are "joined" as a single longer line */
JoinLinesLoop:;
  c = fgetc (fd);
  if (c == '#')
    {
ChompLine:;
      do
        c = fgetc (fd);
      while (c != EOF && c != '\n');
      if (c == '\n')
        goto JoinLinesLoop;
    }
  
  while (c != EOF && c != '\n' && c != ';')
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
      if (!isspace (c))
	expr[nCurSize++] = c;

      if (nCurSize == nMaxSize)
      {
	char *cp;
	nMaxSize += 5000;
	cp = (char *) realloc (expr, nMaxSize + 1);
	if (!cp)
	{
	  free (expr);
	  fprintf (stderr, "Severe warning!, out of core memory reading number!\n");
	  exit (EXIT_FAILURE);
	}
	expr = cp;
      }
      c = fgetc (fd);
    }
  expr[nCurSize] = 0;
  if (!nCurSize)
    ret = 0;
  else
    {
      if (expr[nCurSize-1] == '\\')
	{
	  /* remove the '\\' char, and then process the next line */
	  expr[--nCurSize] = 0;
	  goto JoinLinesLoop;
	}
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
	  FREE (s, strlen (s) + 1);
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

  if (expr == NULL)
    {
      fprintf (stderr, "Error: not enough memory\n");
      exit (EXIT_FAILURE);
    }
/* Lines ending in '\\' are "joined" as a single longer line */
  c = cp;
JoinLinesLoop:;
  if (*c == '#')
    {
      do
        ++c;
      while (*c && *c != '\n');
      if (*c == '\n')
        goto JoinLinesLoop;
    }
  
  while (*c && *c != '\n' && *c != ';')
    {
      /* strip space and tabs out here, and then we DON'T have to mess with them in the rest of the parser */
      if (!isspace(*c))
	expr[nCurSize++] = *c;
      if (nCurSize == nMaxSize)
      {
	char *cp;
	nMaxSize += 5000;
	cp = (char *) realloc (expr, nMaxSize + 1);
	if (!cp)
	{
	  free(expr);
	  fprintf(stderr, "Severe warning!, out of core memory reading number!\n");
	  exit (EXIT_FAILURE);
	}
	expr = cp;
      }
      ++c;
    }
  expr[nCurSize] = 0;
  if (!nCurSize)
    ret = 0;
  else
    {
      if (expr[nCurSize-1] == '\\')
	{
	  /* remove the '\\' char, and then process the next line */
	  expr[--nCurSize] = 0;
	  goto JoinLinesLoop;
	}
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
	  FREE (s, strlen (s) + 1);
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
  else if ('@'==op)	/* Multi factorial   (syntax n!prior_n.  example: 15!3 == 15*12*9*6*3) */
    {
      long nCur;
      unsigned long nDecr;
      nCur = mpz_get_si(prior_n);
      nDecr = mpz_get_ui(n);
      mpz_set_ui(n,1);
      /*printf ("Multi-factorial  %ld!%ld\n", nCur, nDecr);*/
      while (nCur > 1)
	{
	  /* This could be done much more efficiently (bunching mults using smaller "built-ins"), but I am not going to bother for now */
	  mpz_mul_ui(n,n,nCur);
	  nCur -= nDecr;
	}
    }
  else if ('#'==op)  /* simple primorial (syntax  n#   example:  11# == 2*3*5*7*11 */
    {
      long nMax;
      double p;
      nMax = mpz_get_si(n);
      mpz_set_ui(n,1);
      getprime (FREE_PRIME_TABLE);  /* free the prime tables, and reinitialize */
      for (p = 2.0; p <= nMax; p = getprime (p))
	/* This could be done much more efficiently (bunching mults using smaller "built-ins"), but I am not going to bother for now */
	mpz_mul_ui(n,n,(unsigned)p);
    }
  else if ('$'==op)  /* reduced primorial (syntax  n#prior_n   example:  13#5 == (5*7*11*13) */
    {
      double p;
      long nMax;
      unsigned long nStart;
      nMax = mpz_get_si(prior_n);
      nStart = mpz_get_ui(n);
      mpz_set_ui(n,1);
      getprime (FREE_PRIME_TABLE);  /* free the prime tables, and reinitialize */
      p = getprime (nStart);
      /*printf ("Reduced-primorial  %ld#%ld\n", nMax, nStart);*/
      for (; p <= nMax; p = getprime (p))
	{
	  /* Unfortunately, the SoE within GMP-ECM does not always start
	     correctly, so we have to skip the low end stuff by hand */
	  if (p >= nStart)
	    /* This could be done much more efficiently (bunching mults using smaller "built-ins"), but I am not going to bother for now */
	    mpz_mul_ui(n,n,(unsigned)p);
	}
    }
}

void eval_product (mpz_t prior_n, mpz_t n,char op)
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
  if ('*'==op || '.'==op)
    mpz_mul(n,prior_n,n);
  else if ('/'==op)
    mpz_tdiv_q(n,prior_n,n);
  else if ('%'==op)
    mpz_tdiv_r(n,prior_n,n);
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

int eval_Phi (mpz_t b, mpz_t n, int ParamCnt)
{
  int factors[200];
  unsigned dwFactors=0, dw;
  int N, B;
  double p;
  mpz_t D, T, org_n;

  if (ParamCnt == 0)
    {
      fprintf (stderr, "\nParsing Error - the Phi function (in ECM) requires 2 parameters\n");
      return 0;
    }

  if (mpz_cmp_ui(n, 1) == 0)
    {
      /* return value is 1 if b is composite, or b if b is prime */
      int isPrime = mpz_probab_prime_p (b, PROBAB_PRIME_TESTS);
      if (isPrime)
	mpz_set(n, b);
      else
	mpz_set(n, mpOne);
      return 1;
    }
  if (mpz_cmp_si(n, -1) == 0)
    {
      /* this is actually INVALID, but it is easier to simply */
      fprintf (stderr, "\nParsing Error - Invalid parameter passed to the Phi function\n");
      return 0;
    }
  /* OK parse the Phi out now */
  if (mpz_cmp_ui(b, 0) == 0)
    {
      /* this is valid, but return that it is NOT */
      mpz_set(n, mpOne);
      return 0;
    }
  if (mpz_cmp_ui(b, 1) == 0)
    {
      if (mpz_cmp_ui(n, 1) != 0)
	mpz_sub_ui(n, n, 1);
      return 1;
    }

  /* Ok, do the real h_primative work, since we are not one of the trivial case */

  N = mpz_get_si(n);
  B = mpz_get_si(b);

  if (mpz_cmp_ui(b, B))
    {
      fprintf (stderr, "\nParsing Error - Invalid parameter passed to the Phi function (first param B too high)\n");
      return 0;
    }

  /* Obtain the factors of B */
  getprime (FREE_PRIME_TABLE);  /* free the prime tables, and reinitialize */
  for (p = 2.0; p <= B; p = getprime (p))
    {
      if (B % (int) p == 0)
	{
	  /* Add the factor one time */
	  factors[dwFactors++] = (int) p;
	  /* but be sure to totally remove it */
	  do { B /= (int) p; } while (B % (int) p == 0);
        }
     }
  B = mpz_get_si(b);

  mpz_init_set(org_n, n);
  mpz_set_ui(n, 1);
  mpz_init_set_ui(D, 1);
  mpz_init(T);

	      
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

/* A simple partial-recursive decent parser */
int eval_2 (int bInFuncParams)
{
  mpz_t n_stack[4];
  mpz_t n;
  int i;
  int num_base;
  char op_stack[4];
  char op;
  char negate;
  for (i=0;i<4;i++)
    {
      op_stack[i]=0;
      mpz_init(n_stack[i]); 
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
	  for (i=0;isdigit(expr_str[i]);i++)
	    ;
	  if (!i)         /* No digits found */
	    {
	      /* check for a valid "function" */
	      if (!strncasecmp (&expr_str[i], "phi(", 4))
		{
		  /* Process the phi(B,N) function */
		  expr_str+=4;
		  /* eval the first parameter.  NOTE we pass a 1 since we ARE in parameter mode, 
		     and this causes the ',' character to act as the end of expression */
		  if (eval_2 (1) != 2)
		    {
		      fprintf (stderr, "Error, Function Phi() requires 2 parameters in GMP-ECM syntax\n");
		      mpz_clear(n);
		      for (i=0;i<4;i++) 
			mpz_clear(n_stack[i]);
		      return 0;
		    }
		  /* Save off the parameter */
		  mpz_set(n_stack[3], t);
		  /* Now eval the second parameter NOTE we pass a 0 since we are NOT expecting a ','
		     character to end the expression, but are expecting a ) character to end the function */
		  if (eval_2 (0))
		    {
		      mpz_set(n, t);
		      if (eval_Phi (n_stack[3], n, 1) == 0)
			{
			  mpz_clear(n);
			  for (i=0;i<4;i++) 
			    mpz_clear(n_stack[i]);
			  return 0;
			}
		    }
		  goto MONADIC_SUFFIX_LOOP;
		}
	      else
		{
		  fprintf (stderr, "\nError - invalid number [%c]\n", expr_str[i]);
		  mpz_clear(n);
		  for (i=0;i<4;i++) 
		    mpz_clear(n_stack[i]);
		  return 0;
		}
	    }
	  /* Now check for a hex number.  If so, handle it as such */
	  num_base=10;  /* assume base 10 */
	  if (i == 1 && !strncasecmp(expr_str, "0x", 2))
	    {
		num_base = 16;	/* Kick up to hex */
		expr_str += 2;	/* skip the 0x string of the number */
		for (i=0;isxdigit(expr_str[i]);i++)
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
MONADIC_SUFFIX_LOOP:;
        op=*expr_str++;
	    
      if (0==op || ')'==op || ']'==op || '}'==op || (','==op&&bInFuncParams))
	{
	  eval_power (n_stack[2],n,op_stack[2]);
	  eval_product (n_stack[1],n,op_stack[1]);
	  eval_sum (n_stack[0],n,op_stack[0]);
	  mpz_set(t, n);
	  mpz_clear(n);
	  for (i=0;i<4;i++) 
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
	      if (!isdigit(*expr_str))
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
	      if (!isdigit(*expr_str))
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
		      mpz_clear(n);
		      for (i=0;i<4;i++) 
			mpz_clear(n_stack[i]);
		      return 0;
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
