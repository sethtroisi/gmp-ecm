/* Functions for reading a writing resume file lines.

  Copyright 2001, 2002, 2003 Alexander Kruppa and Paul Zimmermann.

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
#include <time.h>
#include <string.h>
#include "gmp.h"
#include "ecm.h"

#define DEBUG

/* Reads a string of characters from fd while they match the string s.
   Returns the number of matching characters that were read. 
*/

int 
facceptstr (FILE *fd, char *s)
{
  int c;
  unsigned i = 0;
  
  while (s[i] != 0 && (c = fgetc (fd)) != EOF)
    {
      if (c != s[i++])
        {
          ungetc (c, fd);
          return i-1;
        }
    }
  
  return i;
}

/* Reads a string from fd until the character "delim" is seen, or "len" 
   characters have been written to s (including terminating null), or 
   EOF is reached. If s is NULL, characters are read from fd but not written
   anywhere.
   Returns the number of characters read. */
int 
freadstrn (FILE *fd, char *s, char delim, unsigned int len)
{
  unsigned int i = 0;
  char c;
  
  while (i + 1 < len && (c = fgetc (fd)) != EOF)
    if (c == delim)
      {
        ungetc (c, fd);
        break;
      }
    else
      if (s != NULL)
        s[i++] = c;
  
  if (i < len && s != NULL)
    s[i++] = 0;
  
  return i;
}

int 
read_resumefile_line (int *method, mpz_t x, mpz_t n, mpz_t sigma, mpz_t A, 
        mpz_t x0, double *b1, char *program, char *who, char *rtime, 
        char *comment, FILE *fd)
{
  int a, c;
  int have_method, have_x, have_n, have_sigma, have_a, have_b1, have_checksum;
  unsigned int saved_checksum;
  
  while (!feof (fd))
    {
      /* Ignore empty lines */
      if (facceptstr (fd, "\n"))
        continue;
      
      /* Ignore lines beginning with '#'*/
      if (facceptstr (fd, "#"))
        {
          while ((c = fgetc (fd)) != EOF && c != '\n');
          continue;
        }
      
      if (feof (fd))
        break;
      
      have_method = have_x = have_n = have_sigma = have_a = have_b1 = 0;
      have_checksum = 0;

      /* Set optional fields to zero */
      mpz_set_ui (sigma, 0);
      mpz_set_ui (A, 0);
      if (program != NULL)
        program[0] = 0;
      if (who != NULL)
        who[0] = 0;
      if (rtime != NULL)
        rtime[0] = 0;
      if (comment != NULL)
        comment[0] = 0;

      while (!facceptstr (fd, "\n") && !feof (fd))
        {
          if ((a = facceptstr (fd, "METHOD=")))
            {
              if (a != 7)
                goto error;

              if (facceptstr (fd, "ECM") == 3)
                {
                  *method = EC_METHOD;
                }
              else if (facceptstr (fd, "P"))
                {
                  if ((a = facceptstr (fd, "M1") == 2))
                    {
                      *method = PM1_METHOD;
                    }
                  else if (a == 0 && facceptstr (fd, "P1") == 2)
                    {
                      *method = PP1_METHOD;
                    }
                  else 
                    goto error;
                }
              else
                goto error;

              have_method = 1;
            }
          
          else if (facceptstr (fd, "X"))
            {
              if (facceptstr (fd, "="))
                {
                  mpz_inp_str (x, fd, 0);
                  have_x = 1;
                }
              else if (facceptstr (fd, "0=") == 2)
                {
                  mpz_inp_str (x0, fd, 0);
                }
              else
                goto error;
            }
          
          else if (facceptstr (fd, "C"))
            {
              if ((a = facceptstr (fd, "HECKSUM=")))
                {
                  if (a != 8)
                    goto error;

                  fscanf (fd, "%u", &saved_checksum);
                  have_checksum = 1;
                }
              else if ((a = facceptstr (fd, "OMMENT=")))
                {
                  if (a != 7)
                    goto error;
                  
                  freadstrn (fd, comment, ';', 255);
                }
              else
                goto error;
            }

          else if ((a = facceptstr (fd, "N=")))
            {
              if (a != 2)
                goto error;

              mpz_inp_str (n, fd, 0);
              have_n = 1;
            }
          
          else if ((a = facceptstr (fd, "SIGMA=")))
            {
              if (a != 6)
                goto error;

              mpz_inp_str (sigma, fd, 0);
              have_sigma = 1;
            }
          
          else if ((a = facceptstr (fd, "A=")))
            {
              if (a != 2)
                goto error;

              mpz_inp_str (A, fd, 0);
              have_a = 1;
            }
          
          else if ((a = facceptstr (fd, "B1=")))
            {
              if (a != 3)
                goto error;

              fscanf (fd, "%lf", b1);
              have_b1 = 1;
            }

          else if ((a = facceptstr (fd, "PROGRAM=")))
            {
              if (a != 8)
                goto error;

              freadstrn (fd, program, ';', 255);
            }

          else if ((a = facceptstr (fd, "WHO=")))
            {
              if (a != 4)
                goto error;

              freadstrn (fd, who, ';', 255);
            }

          else if ((a = facceptstr (fd, "TIME=")))
            {
              if (a != 5)
                goto error;

              freadstrn (fd, rtime, ';', 255);
            }
          else /* Not a tag we know about */
            goto error;
         
          if (!facceptstr (fd, ";"))
            goto error;
          
          while (facceptstr (fd, " "));
        }
      
      /* Finished reading tags */
      
#ifdef DEBUG
      if (*method != EC_METHOD && (have_sigma || have_a))
        {
          fprintf (stderr, "Save file line has ");
          if (have_sigma)
            {
              fprintf (stderr, "SIGMA");
              mpz_set_ui (sigma, 0);
            }
          if (have_sigma && have_a)
            fprintf (stderr, " and ");
          if (have_a)
            {
              fprintf (stderr, "A");
              mpz_set_ui (A, 0);
            }
          fprintf (stderr, " value for method other than ECM.\n");
        }
#endif
      
      if (have_method && have_x && have_n && have_b1 &&
           (method != EC_METHOD || have_sigma || have_a))
        {
          if (have_checksum)
            {
              mpz_t checksum;
              
              mpz_init (checksum);
              mpz_set_d (checksum, *b1);
              if (have_sigma)
                mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
              if (have_a)
                mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (A, CHKSUMMOD));
              mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n, CHKSUMMOD));
              mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (x, CHKSUMMOD));
              if (mpz_fdiv_ui (checksum, CHKSUMMOD) != saved_checksum)
                {
                  fprintf (stderr, "Resume file line has bad checksum %u, expected %lu\n", 
                           saved_checksum, mpz_fdiv_ui (checksum, CHKSUMMOD));
                  mpz_clear (checksum);
                  continue;
                }
              mpz_clear (checksum);
            }
          return 1;
        }
      
      fprintf (stderr, "Save file line lacks fields, have_method = %d, "
                       "have_x = %d, have_n = %d\n",
               have_method, have_x, have_n);
      continue;
      
error:
      /* In case of error, read rest of line and try next line */
      c = fgetc (fd);
      while (c != EOF && c != '\n')
        c = fgetc (fd);
    }
    
    /* We hit EOF without reading a proper save line */
    return 0;
}

void 
write_resumefile_line (FILE *fd, int method, double B1, mpz_t sigma, mpz_t A, 
	mpz_t x, mpz_t n, mpz_t x0, char *comment)
{
  mpz_t checksum;
  time_t t;
  char timestring[256];

#ifdef DEBUG
  if (fd == NULL)
    {
      fprintf (stderr, "write_resumefile_line: fd == NULL\n");
      exit (EXIT_FAILURE);
    }
#endif
  
  mpz_init (checksum);
  mpz_set_d (checksum, B1);
  fprintf (fd, "METHOD=");
  if (method == PM1_METHOD)
    fprintf (fd, "PM1");
  else if (method == PP1_METHOD)
    fprintf (fd, "PP1");
  else 
    {
      fprintf (fd, "ECM; ");
      if (mpz_sgn (sigma) != 0)
        {
          fprintf (fd, "SIGMA=");
          mpz_out_str (fd, 10, sigma);
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
        }
      if (mpz_sgn (A) != 0)
        {
          fprintf (fd, "A=");
          mpz_out_str (fd, 10, A);
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (A, CHKSUMMOD));
        }
    }
  
  fprintf (fd, "; B1=%.0f; N=", B1);
  mpz_out_str (fd, 10, n);
  fprintf (fd, "; X=0x");
  mpz_out_str (fd, 16, x);
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n, CHKSUMMOD));
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (x, CHKSUMMOD));
  fprintf (fd, "; CHECKSUM=%lu; PROGRAM=GMP-ECM %s;",
           mpz_fdiv_ui (checksum, CHKSUMMOD), ECM_VERSION);
  mpz_clear (checksum);
  
  if (mpz_sgn (x0) != 0)
    {
      fprintf (fd, " X0=0x");
      mpz_out_str (fd, 16, x0);
      fprintf (fd, ";");
    }

  if (comment[0] != 0)
    fprintf (fd, " COMMENT=%.255s;", comment);
  
  t = time (NULL);
  strncpy (timestring, ctime (&t), 255);
  timestring[255] = 0;
  timestring[strlen (timestring) - 1] = 0; /* Remove newline */
  fprintf (fd, " TIME=%s;", timestring);
  fprintf (fd, "\n");
}
