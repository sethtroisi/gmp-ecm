/* This programs writes the string "1" to the file /proc/sys/vm/drop_caches.
   Under Linux, this causes the kernel to drop the page buffer, including
   cached file data from the disk devices. This is useful for consecutive 
   timing runs with a data set that is small enough that it would stay in 
   disk cache and thus not give a realistic estimate for the data transfer 
   speed. The drop_caches file is writeably only for root, so this program 
   must be executed with sudo. */

#include <stdlib.h>
#include <stdio.h>

int main()
{
  FILE *f;
  const char *fn = "/proc/sys/vm/drop_caches";
  
  f = fopen (fn, "w");
  if (f == NULL)
    {
      fprintf (stderr, "Could not open %s for writing\n", fn);
      exit (EXIT_FAILURE);
    }
  if (fprintf (f, "1\n") != 2)
    {
      fprintf (stderr, "Could not print \"1\" to %s\n", fn);
      exit (EXIT_FAILURE);
    }
  fclose (f);
  exit (EXIT_SUCCESS);    
}
