
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#if !defined (_MSC_VER)
#include <unistd.h>
#else
/* MinGW built GMP causes grief when gmp_*printf() is used.  This "works" around the problem */
#undef isascii
int isascii(int a) {return __isascii(a); }
#endif

#include "gmp.h"
#include "ecm.h"

/* Options for using an external PRPer rather than internal GMP */
char *externalprp = NULL; /* call out to external program  */  
char *externallog = NULL; /* where to find output */
char *externalinputprpfile = NULL; /* Used to place the n value (a temp file). Is deleted after system */
char *externalisprp = "-PRP!"; /* what to match output against */
char *externaliscomposite = "is composite"; /* ditto */
int externalprplen=800; /* length where we go external */
int externalprpval=-1; /* exit value meaning it's prp, -1 is impossible */

static int warned=0;

int smart_probab_prime_p(mpz_t const n, int c)
{
  int numlen, bufflen;
  char const *read;
  char *command, *build;
  int doneN, doneL, doneT, returnval;

  numlen=mpz_sizeinbase(n, 10);
  if(numlen<=externalprplen)
    {
      return mpz_probab_prime_p(n, c<3?3:c); /* for tiny numbers perform more MR's.  Perform at least 3, no matter what PROBAB_PRIME_TESTS is defined at in ecm.h */
    }
  
  if(!externalprp || warned)
    return mpz_probab_prime_p(n, c);

  if (externalprpval == -1 && externallog == NULL)
    {
      fprintf(stderr, "Usage of external PRP requires either return value listed, or external log file.\nNeither were present, so using GMP\n");
      warned=1;
      return mpz_probab_prime_p(n, c);
    }
  if(!strstr(externalprp,"%n") && !externalinputprpfile)
    {
      fprintf(stderr, "PRP program contains no %%n parameter and no external input file was listed, using GMP\n");
      warned=1;
      return mpz_probab_prime_p(n, c);
    }
  
  bufflen=50+strlen(externalprp)+numlen+(externallog&&strstr(externalprp,"%l")?strlen(externallog):0)+(externalinputprpfile&&strstr(externalprp,"%t")?strlen(externalinputprpfile):0);
  command=malloc(bufflen);
  if(!command)
    {
      fprintf (stderr, "Error: not enough memory\n");
      exit (EXIT_FAILURE);
    }
      
  read=externalprp;
  build=command;
  doneL=doneN=doneT=0;
  
  while(*read)
    {
      if(*read=='%')
        {
          if(read[1]=='n' && !doneN)
            {
              build+=gmp_sprintf(build, "%Zi", n);
              doneN=1;
              ++read;
            }
          else if(read[1]=='l' && !doneL)
            {
              build+=sprintf(build, "%s", externallog);
              doneL=1;
              ++read;
            }
          else if(read[1]=='t' && !doneT)
            {
              build+=sprintf(build, "%s", externalinputprpfile);
              doneT=1;
              ++read;
            }
          else if(read[1]=='x')
            {
              int chVal=0,chCnt=0;
              ++read;
              while (isxdigit(read[1]) && chCnt < 2)
                {

                   chVal *= 16;
                   chVal += (read[1]>='0' && read[1]<='9') ? read[1]-'0' : tolower(read[1])-'a'+10;
                   ++read;
                   ++chCnt;
                }
              build+=sprintf(build, "%c", chVal);
            }
          else 
            *build++='%';
        }
      else
        *build++=*read;
      if (*read) 
        ++read;
    }
  *build='\0';
  /*fprintf(stderr, "Running %s\n", command);*/
  
  /* WARNING - could delete important data! */
  if(externallog)
    unlink(externallog);
  if (externalinputprpfile)
    {
      FILE *out = fopen(externalinputprpfile, "w");
      if (!out)
        {
          fprintf(stderr, "Error creating temp PRP file %s so using GMP\n", externalinputprpfile);
          free (command);
          return mpz_probab_prime_p(n, c);
        } 
      gmp_fprintf(out, "%Zi", n);
      fclose(out);
    }
  returnval=system(command);
  /* WARNING - could delete important data! */
  if (externalinputprpfile)
    unlink(externalinputprpfile);

  if (returnval>=0)
    {
      /* First choice is to look at the return value, if that's useful */
      if (externalprpval >= 0)
        {
          returnval = (returnval==externalprpval);
        }
      /* Second choice is to plough through any log file */
      else if (externallog)
        {
	  FILE*lf;
          returnval=-1;
          lf=fopen(externallog,"r");
          if(lf)
            {
              while(fgets(command, bufflen, lf))
                {
                  if(externalisprp && strstr(command, externalisprp))
                    {
                      returnval=1;
                      break;
                    }
                  if(externaliscomposite && strstr(command, externaliscomposite))
                    {
                      returnval=0;
                      break;
                    }
                }
              fclose(lf);
            }
        }
        /* Else, the default's unknown state - fall back to internal */
       else
          returnval=-1;
    }

  if(externallog)
    unlink(externallog);
 
  if (returnval<0)
    {
      fprintf(stderr, "External PRP failed, using internal GMP\n");
      returnval = mpz_probab_prime_p(n, c);
    }              
          
  free (command);
  return returnval;
}
