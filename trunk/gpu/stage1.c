/* Plain C stage 1 (using only the basic GMP operations for the critical loop).

  Copyright 2010 Julie Feltin and Paul Zimmermann.

  This file is part of the ECM Library.

  The ECM Library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or (at your
  option) any later version.

  The ECM Library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the ECM Library; see the file COPYING.LIB.  If not, write to
  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
  MA 02110-1301, USA.
*/

/* Example:
   $ ./stage1 29799904256775982671863388319999573561548825027149399972531599612392671227006866151136667908641695103422986028076864929902803267437351318167549013218980573566942647077444419419003164546362008247462049 17 1000000
a=9122032453422309303975228503303119186594096764565570324571255124462257150200764885675475282720699595928510737885986837078892475054042869967855173419156369123531441860593128494129801744621300240933171
Starting point: x=15103147079049907765711344782576620836201216022308570309709685128202269764959718444195054567986529779670921068575492301131280304481783998012863228763537884832421573625653695166561750401066843178542817
After stage 1, x=17628830287311089492501926991334045672323900629597277022267290291143065118509513624278744509458318241685541031462343435507840406720158512915044990777024060216890526386080098985253934838627077872164414
*/

#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include "prototype.h"
#include "getprime.h"

int
main (int argc, char*argv[])
{
  if (argc != 4)
    printf ("Error in call function\n./stage1 N sigma B1\n");
  else
    {
      mpz_t N;
      mpz_t sigma;
      mpz_t B1;
      mpz_t x;

      mpz_init (N);
      mpz_init (sigma);
      mpz_init (B1);
      mpz_init (x);

      mpz_set_str (N, argv[1], 10); /* in base 10 */
      mpz_set_str (sigma, argv[2], 10);
      mpz_set_str (B1, argv[3], 10);

      /* check N is odd */
      if (mpz_divisible_ui_p (N, 2))
        {
          fprintf (stderr, "Error, N should be odd\n");
          exit (1);
        }
		
      stageOne (B1, sigma, N, x);
      gmp_printf ("After stage 1, x=%Zd\n",x);
		
      mpz_clear (sigma);
      mpz_clear (N);
      mpz_clear (B1);
      mpz_clear (x);
    }
}
