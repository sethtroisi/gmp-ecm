/* Program running ecm factoring processes in parallel, minimizing simultaneous
   step2.

  Copyright 2005 Torbjörn Granlund.

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

  FIXME: We should clean up its reporting and logging functions, as well
         as its error handling, if we release it.
*/

#include <sys/types.h>		/* for fork, wait4 */
#include <unistd.h>		/* for fork */
#include <fcntl.h>		/* for open */
#include <sys/wait.h>		/* for wait4 */
#include <sys/time.h>		/* for gettimeofday, wait4 */
#include <sys/resource.h>	/* for wait4 */
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

enum foo {NONE = 0, STEP1_RUNNING, STEP2_RUNNING, STEP1_DONE, STEP2_DONE};

struct jobinfo
{
  pid_t pid;
  enum foo state;
  char filename[50];
};

#define RUNNING(x) (x.state == STEP1_RUNNING || x.state == STEP2_RUNNING)

char *progname;

unsigned long
bush ()
{
  unsigned long ran;
  int fd;

  fd = open ("/dev/urandom", O_RDONLY);
  if (fd >= 0)
    {
      unsigned char buf[sizeof (long)];
      size_t nread;
      nread = read (fd, buf, sizeof (long));
      if (nread != sizeof (long))
	goto stupid;
      ran = (buf[0] << 24) + (buf[1] << 16) + (buf[2] << 8) + buf[3];
      if (sizeof (long) > 4)
	{
	  unsigned long ran2;
	  ran2 = (buf[4] << 24) + (buf[5] << 16) + (buf[6] << 8) + buf[7];
	  ran += (ran2 << 31) << 1;
	}
      close (fd);
      return ran;
    }
  else
    {
      static int flag = 0;
    stupid:
      if (flag == 0)
	{
	  struct timeval tp;
	  gettimeofday (&tp, NULL);
	  srand48 ((tp.tv_sec << 16) + tp.tv_usec + getpid ());
	  flag = 1;
	}
      ran = mrand48 ();
    }
  return ran;
}

char *
pathfind (const char *command)
{
  char *path, *p, *buf;
  int len, clen;

  clen = strlen (command);

  path = getenv ("PATH");
  if (path == NULL)
    abort ();

  buf = malloc (strlen (path) + 1);

  for (;;)
    {
      p = strchr (path, ':');
      if (p == NULL)
	len = strlen (path);
      else
	len = p - path;

      memcpy (buf, path, len);
      if (buf[len - 1] != '/')
	{
	  buf[len] = '/';
	  memcpy (buf + len + 1, command, clen + 1);
	}
      else
	{
	  memcpy (buf + len, command, clen + 1);
	}
      if (access (buf, X_OK) == 0)
	return buf;

      path += len + 1;
    }

  free (buf);
  return NULL;
}

#define BUFSIZE 65536

int
main (int argc, char *argv[], char *envp[])
{
  int nprocs, i;
  char *B1, *B2;
  struct jobinfo *jiv;
  pid_t pid;
  int wstat;
  char *tmpdir;
  char filename[50];
  int (*result_chan)[2];		/* for reading output from passed */
  int n_running_procs;
  int fd;
  char *ecmfactor;
  char buf[BUFSIZE];
  char sigma[20];
  int next_cofac_i;
  size_t nread;
  struct rusage rus;
  unsigned used_ms;

  nprocs = 1;			/* default */
  B1 = NULL;
  B2 = NULL;
  next_cofac_i = 0;

  progname = argv[0];
  argv++;
  argc--;

  while (argc >= 2 && argv[0][0] == '-')
    {
      if (strcmp ("-B1", argv[0]) == 0)
	{
	  B1 = argv[1];
	  argv += 2;
	  argc -= 2;
	}
      else if (strcmp ("-B2", argv[0]) == 0)
	{
	  B2 = argv[1];
	  argv += 2;
	  argc -= 2;
	}
      else if (strcmp ("-n", argv[0]) == 0)
	{
	  nprocs = strtoul (argv[1], 0, 0);
	  argv += 2;
	  argc -= 2;
	}
      else
	{
	  fprintf (stderr, "%s: unknown option: %s\n", progname, argv[0]);
	  exit (1);
	}
    }
  printf ("There seem to be %d cofactor files\n", argc);

  if (B1 == NULL)
    {
      fprintf (stderr, "%s: missing B1 value\n", progname);
      exit (1);
    }

  ecmfactor = pathfind ("ecmfactor");

  result_chan = malloc (nprocs * sizeof (int [2]));

  jiv = malloc (nprocs * sizeof (struct jobinfo));
  for (i = 0; i < nprocs; i++)
    {
      jiv[i].pid = 0;
      jiv[i].state = NONE;
      pipe (result_chan[i]);
    }

  tmpdir = getenv ("TMPDIR");
  if (tmpdir == NULL)
    tmpdir = "/tmp";

  n_running_procs = 0;

  for (;;)
    {
      int n_step2_procs = 0;
      for (i = 0; i < nprocs; i++)
	{
	  n_step2_procs += (jiv[i].state == STEP2_RUNNING);
	}
      for (i = 0; i < nprocs; i++)
	{
	  if (! RUNNING (jiv[i]))
	    {
	      if (jiv[i].state == STEP1_DONE)
		{
		  if (n_step2_procs * 2 >= nprocs)
		    continue;
		  fprintf (stderr, "STARTING NEW STEP 2 JOB\n");
		  pid = fork ();
		  if (pid == 0)
		    {
		      /* Child */
		      char *outv[6], **op = outv;
		      dup2 (result_chan[i][1], 1);
		      close (result_chan[i][0]);
		      close (result_chan[i][1]);
		      *op++ = "ecmfactor";
		      *op++ = "-resume";
		      *op++ = jiv[i].filename;
		      *op++ = B1;
		      if (B2 != NULL)
			*op++ = B2;
		      *op = NULL;
		      execve (ecmfactor, outv, envp);
		      fprintf (stderr, "cannot execute %s\n", ecmfactor);
		      abort ();
		    }
		  n_running_procs++;
		  jiv[i].pid = pid;
		  jiv[i].state = STEP2_RUNNING;
		  continue;
		}

	      for (;;)
		{
		  if (next_cofac_i == 0)
		    sprintf (sigma, "%lu", bush ());
		  strcpy (filename, argv[next_cofac_i]);
		  next_cofac_i = (next_cofac_i + 1) % argc;
		  fd = open (filename, O_RDONLY);
		  if (fd != -1)
		    break;
		  usleep (50000);
		}

	      sprintf (jiv[i].filename, "%s/ecm-save-%u", tmpdir, i);
	      unlink (jiv[i].filename);
	      fprintf (stderr, "STARTING NEW STEP 1 JOB\n");
	      pid = fork ();
	      if (pid == 0)
		{
		  /* Child */
		  char *outv[8], **op = outv;
		  dup2 (result_chan[i][1], 1);
		  close (result_chan[i][0]);
		  close (result_chan[i][1]);
		  dup2 (fd, 0);
		  close (fd);
		  *op++ = "ecmfactor";
		  *op++ = "-save";
		  *op++ = jiv[i].filename;
		  *op++ = "-sigma";
		  *op++ = sigma;
		  *op++ = B1;
		  *op++ = "1";
		  *op = NULL;
		  execve (ecmfactor, outv, envp);
		  fprintf (stderr, "cannot execute %s\n", ecmfactor);
		  abort ();
		}
	      close (fd);
	      n_running_procs++;
	      jiv[i].pid = pid;
	      jiv[i].state = STEP1_RUNNING;
	      continue;
	    }
	}

      fprintf (stderr, "ABOUT TO WAIT (%d jobs running)\n", n_running_procs);
      pid = wait4 (0, &wstat, 0,  &rus);
      if (pid == -1)
	{
	  fprintf (stderr, "wait returned error %d\n", errno);
	  abort ();
	}
      n_running_procs--;

      used_ms = rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;

      if (WIFSIGNALED (wstat))
	{
	  fprintf (stderr, "*** child got signal %d\n", WTERMSIG (wstat));
	  exit (1);
	}
      for (i = 0; i < nprocs; i++)
	{
	  if (jiv[i].pid == pid)
	    goto yee;
	}
      abort ();
    yee:
      if (jiv[i].state == STEP1_RUNNING)
	{
	  nread = read (result_chan[i][0], buf, BUFSIZE);
	  jiv[i].state = STEP1_DONE;
	  fprintf (stderr, "STEP 1 JOB %d FINISHED (used %u ms)\n", i, used_ms);
	}
      else if (jiv[i].state == STEP2_RUNNING)
	{
	  unlink (jiv[i].filename);
	  nread = read (result_chan[i][0], buf, BUFSIZE);
	  jiv[i].state = STEP2_DONE;
	  fprintf (stderr, "STEP 2 JOB %d FINISHED (used %u ms)\n", i, used_ms);
	}
      else
	abort ();

      if (WEXITSTATUS (wstat) == 0)
	{
	  FILE *fs;
	  fs = popen ("mail -s \"NEW FACTOR\" tg@swox.se", "w");
	  fwrite (buf, nread, 1, fs);
	  pclose (fs);
	}
    }
}

/*
  To perform only step 1:  $ ./ecm -save toto B1 1 < composite
  Then to perform step 2:  $ ./ecm -resume toto B1 [B2]
*/
