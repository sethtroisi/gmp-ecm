/* Functions for reading a writing resume file lines.

  Copyright 2001, 2002, 2003 Paul Zimmermann and Alexander Kruppa.

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
#if !defined (_MSC_VER)
#include <unistd.h>
#endif
#include <gmp.h>
#include "ecm.h"
#include "ecm-ecm.h"

#if defined (_MSC_VER) || defined (__MINGW32__)
/* needed to declare GetComputerName() for write_resumefile_line() */
#include <windows.h>
#endif

/* Reads a string of characters from fd while they match the string s.
   Returns the number of matching characters that were read. 
*/

static int 
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

/* Reads a string from fd until the character "delim" or newline is seen, or 
   "len" characters have been written to "s" (including terminating null), 
   or EOF is reached. The "delim" and newline characters are left on the 
   stream.
   If s is NULL, characters are read from fd but not written anywhere.
   Returns the number of characters read.
*/

static int 
freadstrn (FILE *fd, char *s, char delim, unsigned int len)
{
  unsigned int i = 0;
  int c;
  
  while (i + 1 < len && (c = fgetc (fd)) != EOF)
    if (c == delim || IS_NEWLINE(c))
      {
        ungetc (c, fd);
        break;
      }
    else
      if (s != NULL)
        s[i++] = (char) c;
  
  if (i < len && s != NULL)
    s[i++] = 0;
  
  return i;
}

/* Reads an assignment from a save file. Return 1 if an assignment was
   successfully read, 0 if there are no more lines to read (at EOF) 
*/

int 
read_resumefile_line (int *method, mpz_t x, mpcandi_t *n, mpz_t sigma, mpz_t A, 
        mpz_t x0, double *b1, char *program, char *who, char *rtime, 
        char *comment, FILE *fd)
{
  int a, c;
  int have_method, have_x, have_z, have_n, have_sigma, have_a, have_b1, 
      have_checksum, have_qx;
  unsigned int saved_checksum;
  char tag[16];
  mpz_t z;
  
  while (!feof (fd))
    {
      /* Ignore empty lines */
      if (facceptstr (fd, "\n"))
        {
          continue;
        }
      
      /* Ignore lines beginning with '#'*/
      if (facceptstr (fd, "#"))
        {
          while ((c = fgetc (fd)) != EOF && !IS_NEWLINE(c));
          continue;
        }
      
      if (feof (fd))
        break;
      
      have_method = have_x = have_z = have_n = have_sigma = have_a = 
                    have_b1 = have_qx = 0;
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
          freadstrn (fd, tag, '=', 16);
          
          if (!facceptstr (fd, "="))
            {
              printf ("Save file line has no equal sign after: %s\n", tag);
              goto error;
            }
          
          if (strcmp (tag, "METHOD") == 0)
            {
              if (facceptstr (fd, "ECM") == 3)
                {
                  *method = ECM_ECM;
                }
              else if (facceptstr (fd, "P"))
                {
                  a = facceptstr (fd, "-1");
                  if (a == 2)
                    {
                      *method = ECM_PM1;
                    }
                  else if (a == 0 && facceptstr (fd, "+1") == 2)
                    {
                      *method = ECM_PP1;
                    }
                  else 
                    goto error;
                }
              else
                goto error;

              have_method = 1;
            }
          else if (strcmp (tag, "X") == 0)
            {
              mpz_inp_str (x, fd, 0);
              have_x = 1;
            }
          else if (strcmp (tag, "Z") == 0)
            {
              mpz_init (z);
              mpz_inp_str (z, fd, 0);
              have_z = 1;
            }
          else if (strcmp (tag, "QX") == 0)
            {
              mpz_inp_str (x, fd, 0);
              have_qx = 1;
            }
          else if (strcmp (tag, "X0") == 0)
            {
              mpz_inp_str (x0, fd, 0);
            }
          else if (strcmp (tag, "CHECKSUM") == 0)
            {
              fscanf (fd, "%u", &saved_checksum);
              have_checksum = 1;
            }
          else if (strcmp (tag, "COMMENT") == 0)
            {
              freadstrn (fd, comment, ';', 255);
            }
          else if (strcmp (tag, "N") == 0)
            {
              /*mpz_inp_str (n, fd, 0);*/
	      /* we want to "maintain" any expressions, which were possibly stored in the file for N */
              have_n = read_number (n, fd, 0);
            }
          else if (strcmp (tag, "SIGMA") == 0)
            {
              mpz_inp_str (sigma, fd, 0);
              have_sigma = 1;
            }
          else if (strcmp (tag, "A") == 0)
            {
              mpz_inp_str (A, fd, 0);
              have_a = 1;
            }
          else if (strcmp (tag, "B1") == 0)
            {
              fscanf (fd, "%lf", b1);
              have_b1 = 1;
            }
          else if (strcmp (tag, "PROGRAM") == 0)
            {
              freadstrn (fd, program, ';', 255);
            }
          else if (strcmp (tag, "WHO") == 0)
            {
              freadstrn (fd, who, ';', 255);
            }
          else if (strcmp (tag, "TIME") == 0)
            {
              freadstrn (fd, rtime, ';', 255);
            }
          else /* Not a tag we know about */
            {
              printf ("Save file line has unknown tag: %s\n", tag);
              goto error;
            }
         
          /* Prime95 lines have no semicolon after SIGMA */
          if (!facceptstr (fd, ";") && ! (have_qx && have_n && have_sigma))
            {
              printf ("%s field not followed by semicolon\n", tag);
              goto error;
            }
          
          while (facceptstr (fd, " "));
        }
      
      /* Finished reading tags */
      
      /* Handle Prime95 v22 lines. These have no METHOD=ECM field and
         QX= instead of X= */
      
      if (have_qx)
        {
          if (have_n && have_sigma)
            {
              *method = ECM_ECM;
              *b1 = 1.0;
              strcpy (program, "Prime95");
              mpz_mod (x, x, n->n);
              return 1;
            }
          goto error;
        }

#ifdef DEBUG
      if (*method != ECM_ECM && (have_sigma || have_a || have_z))
        {
          int count = have_sigma + have_a + have_z;
          printf ("Warning: Save file line has");
          if (have_sigma)
            {
              printf (" SIGMA");
              mpz_set_ui (sigma, 0);
              if (--count > 1)
                printf (",");
              else if (count > 0)
                printf (" and");
            }
          if (have_a)
            {
              printf (" A");
              mpz_set_ui (A, 0);
              if (--count > 0)
                printf (" and");
            }
          if (have_z)
            {
              printf (" Z");
              mpz_clear (Z);
              have_z = 0;
            }
          printf (" value for method other than ECM.\n");
        }
#endif
      
      if (!have_method || !have_x || !have_n || !have_b1 ||
          (method == ECM_ECM && !have_sigma && !have_a))
        {
          fprintf (stderr, "Save file line lacks fields\n");
          continue;
        }

      if (have_checksum)
        {
          mpz_t checksum;
          
          mpz_init (checksum);
          mpz_set_d (checksum, *b1);
          if (have_sigma)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
          if (have_a)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (A, CHKSUMMOD));
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n->n, CHKSUMMOD));
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (x, CHKSUMMOD));
          if (have_z)
            mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (z, CHKSUMMOD));
          if (mpz_fdiv_ui (checksum, CHKSUMMOD) != saved_checksum)
            {
              fprintf (stderr, "Resume file line has bad checksum %u, expected %lu\n", 
                       saved_checksum, mpz_fdiv_ui (checksum, CHKSUMMOD));
              mpz_clear (checksum);
              continue;
            }
          mpz_clear (checksum);
        }

      mpz_mod (x, x, n->n);
      if (have_z)	/* Must normalize */
        {
          if (!mpz_invert (z, z, n->n)) /* Factor found? */
            {
              /* Oh great. What do we do with it now? */
              /* mpres_gcd (f, z, n); */
              printf ("Oops, factor found while reading from save file.\n");
            }
          mpz_mul (z, z, x);
          mpz_mod (x, z, n->n);
        }
        
      return 1;
      
error:
      /* In case of error, read rest of line and try next line */
      c = fgetc (fd);
      while (c != EOF && !IS_NEWLINE(c))
        c = fgetc (fd);
    }
    
    /* We hit EOF without reading a proper save line */
    return 0;
}

void 
write_temp_resumefile (int method, double B1, mpz_t sigma, mpz_t A, mpz_t x, 
                       mpz_t n, mpz_t orig_x0, int verbose)
{
  FILE *fd;
  mpcandi_t data;
  char *comment = "ECM_WIP_AutoSave";

  /* Use a 2 file save method.  It saves to a second file name, and when that is "correctly"
     done, it deletes the final file, then renames from the very temp name, to the real
     incremental name */
  fd = fopen ("gmpecm1.wip", "w");
  if (!fd)
    {
      /* failure to create a temp file is not critical.  GMP-ECM ALWAYS failed to 
         create this before this code was added ;) and it still ran just fine */
      if (verbose >= 2)
	fprintf (stderr, "error, can't open incremental save file gmpecm1.wip\n");
      return;
    }
  
  mpcandi_t_init (&data);
  mpcandi_t_add_candidate (&data, n, NULL, 0);
  write_resumefile_line (fd, method, B1, sigma, A, x, &data, orig_x0, comment);
  mpcandi_t_free (&data);
  fclose (fd);
  remove ("gmp_ecm.wip");
  rename ("gmpecm1.wip", "gmp_ecm.wip");
}

void kill_temp_resume_file (void)
{
  remove ("gmp_ecm.wip");
}

void 
write_resumefile_line (FILE *fd, int method, double B1, mpz_t sigma, mpz_t A, 
	mpz_t x, mpcandi_t *n, mpz_t x0, const char *comment)
{
  mpz_t checksum;
  time_t t;
  char text[256];
  char *uname, mname[32];

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
  if (method == ECM_PM1)
    fprintf (fd, "P-1");
  else if (method == ECM_PP1)
    fprintf (fd, "P+1");
  else 
    {
      fprintf (fd, "ECM");
      if (mpz_sgn (sigma) != 0)
        {
          fprintf (fd, "; SIGMA=");
          mpz_out_str (fd, 10, sigma);
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (sigma, CHKSUMMOD));
        }
      else if (mpz_sgn (A) != 0)
        {
          fprintf (fd, "; A=");
          mpz_out_str (fd, 10, A);
          mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (A, CHKSUMMOD));
        }
    }
  
  fprintf (fd, "; B1=%.0f; N=", B1);
  if (n->cpExpr)
    fprintf(fd, "%s", n->cpExpr);
  else
    mpz_out_str (fd, 10, n->n);
  fprintf (fd, "; X=0x");
  mpz_out_str (fd, 16, x);
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (n->n, CHKSUMMOD));
  mpz_mul_ui (checksum, checksum, mpz_fdiv_ui (x, CHKSUMMOD));
  fprintf (fd, "; CHECKSUM=%lu; PROGRAM=GMP-ECM %s;",
           mpz_fdiv_ui (checksum, CHKSUMMOD), VERSION);
  mpz_clear (checksum);
  
  if (mpz_sgn (x0) != 0)
    {
      fprintf (fd, " X0=0x");
      mpz_out_str (fd, 16, x0);
      fprintf (fd, ";");
    }
  
  /* Try to get the users and his machines name */
  /* TODO: how to make portable? */
  uname = getenv ("LOGNAME");
  if (uname == NULL)
    uname = getenv ("USERNAME");
  if (uname == NULL)
    uname = "";
  
#if defined (_MSC_VER) || defined (__MINGW32__)
  /* dummy block, so that the vars needed here don't need to
    "spill" over to the rest of the function. */
  {
    DWORD size, i;
    TCHAR T[MAX_COMPUTERNAME_LENGTH+2];
    size=MAX_COMPUTERNAME_LENGTH+1;
    if (!GetComputerName(T, &size))
      strcpy(mname, "localPC");
    else
      {
        for (i = 0; i < sizeof(mname)-1; ++i)
          mname[i] = T[i];
        mname[sizeof(mname)-1] = 0;
      }
  }
#else
  if (gethostname (mname, 32) != 0)
    mname[0] = 0;
  mname[31] = 0; /* gethostname() may omit trailing 0 if hostname >31 chars */
#endif
  
  if (uname[0] != 0 || mname[0] != 0)
    {
      fprintf (fd, " WHO=%.233s@%.32s;", uname, mname);
    }

  if (comment[0] != 0)
    fprintf (fd, " COMMENT=%.255s;", comment);
  
  t = time (NULL);
  strncpy (text, ctime (&t), 255);
  text[255] = 0;
  text[strlen (text) - 1] = 0; /* Remove newline */
  fprintf (fd, " TIME=%s;", text);
  fprintf (fd, "\n");
  fflush (fd);
}
