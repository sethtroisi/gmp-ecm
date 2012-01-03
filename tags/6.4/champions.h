/* champions.h: defines the keepers of Top-10 lists for P-1, P+1, and ECM
   factors, and the size that is currently needed to enter the Top-10 */

/* people keeping track of champions and corresponding url's: ECM, P-1, P+1 */
static char *champion_keeper[3] =
{ "Richard Brent <champs@rpbrent.com>",
  "Paul Zimmermann <zimmerma@loria.fr>",
  "Paul Zimmermann <zimmerma@loria.fr>"};
static char *champion_url[3] =
{"http://wwwmaths.anu.edu.au/~brent/ftp/champs.txt",
 "http://www.loria.fr/~zimmerma/records/Pminus1.html",
 "http://www.loria.fr/~zimmerma/records/Pplus1.html"};
/* minimal number of digits to enter the champions table for ECM, P-1, P+1 */
static unsigned int champion_digits[3] = { 67, 53, 48 };
