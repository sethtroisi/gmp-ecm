/* Memory usage utilities.

Copied from CADO-NFS.

This file is part of the ECM Library.

The ECM Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The ECM Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the ECM Library; see the file COPYING.LIB.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#ifndef _WIN32

#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include "ecm-impl.h"

/* Returns peak memory usage, in KB
 * This is the VmPeak field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure.
 */
long
PeakMemusage (void)
{
  pid_t pid = getpid ();

  char str[1024];
  char *truc;
  snprintf (str, 1024, "/proc/%d/status", (int) pid);

  FILE *file;
  file = fopen (str, "r");
  if (file == NULL)
    return -1; /* for example on Mac OS X */

  long mem;
  for(;;)
    {
      truc = fgets (str, 1023, file);
      if (truc == NULL)
	return -1; /* for example on FreeBSD */
      int ret = sscanf (str, "VmPeak: %ld", &mem);
      if (ret == 1)
        {
          fclose (file);
          return mem;
        }
    }
}

#else

#include <windows.h>
#include <psapi.h>

long
PeakMemusage (void)
{
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (long)(info.PeakWorkingSetSize >> 10);
}

#endif
