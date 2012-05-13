/* 
  Implementation of fast stage 2 for P-1 and P+1 as described in
  "Improved Stage 2 to $P\pm{}1$ Factoring Algorithms" by
  Peter L. Montgomery and Alexander Kruppa, ANTS 2008 (8th Algorithmic 
  Number Theory Symposium).
   
  Copyright 2007, 2008 Alexander Kruppa.
  NTT functions are based on code Copyright 2005 Dave Newman.

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

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "ecm-impl.h"
#include "sp.h"
#include <math.h>

/* Some useful PARI functions:
   sumset(a,b) = {local(i, j, l); l = listcreate (length(a) * length(b)); for (i = 1, length(a), for (j = 1, length(b), listput(l, a[i] + b[j]))); listsort (l, 1); l}
*/

const static uint64_t Pvalues[] = {
6ull, 10ull, 18ull, 30ull, 42ull, 34ull, 54ull, 66ull, 90ull, 102ull, 126ull, 
150ull, 210ull, 198ull, 270ull, 330ull, 390ull, 378ull, 462ull, 510ull, 
630ull, 690ull, 714ull, 750ull, 810ull, 870ull, 1050ull, 1170ull, 1230ull, 
1470ull, 1530ull, 1650ull, 1890ull, 2310ull, 2130ull, 2730ull, 2610ull, 
2670ull, 3150ull, 3570ull, 3990ull, 4290ull, 4410ull, 4830ull, 5250ull, 
5610ull, 6090ull, 6930ull, 7350ull, 8190ull, 8610ull, 9030ull, 9450ull, 
9570ull, 10710ull, 11550ull, 11970ull, 11730ull, 13650ull, 14490ull, 16170ull, 
17850ull, 19110ull, 19950ull, 20790ull, 24570ull, 24990ull, 25410ull, 
27930ull, 30030ull, 32130ull, 34650ull, 35490ull, 39270ull, 40950ull, 
43890ull, 46410ull, 48510ull, 51870ull, 53130ull, 53550ull, 57330ull, 
57750ull, 60690ull, 62370ull, 62790ull, 66990ull, 67830ull, 71610ull, 
79170ull, 82110ull, 84630ull, 85470ull, 90090ull, 94710ull, 99330ull, 
101010ull, 108570ull, 111930ull, 117810ull, 122430ull, 131670ull, 136290ull, 
139230ull, 140910ull, 150150ull, 154770ull, 155610ull, 159390ull, 188370ull, 
196350ull, 200970ull, 203490ull, 210210ull, 214830ull, 219450ull, 232050ull, 
237510ull, 253890ull, 256410ull, 259350ull, 265650ull, 270270ull, 274890ull, 
307230ull, 324870ull, 330330ull, 334950ull, 353430ull, 363090ull, 371910ull, 
390390ull, 395010ull, 417690ull, 431970ull, 450450ull, 466830ull, 468930ull, 
478170ull, 482790ull, 510510ull, 570570ull, 589050ull, 630630ull, 658350ull, 
667590ull, 690690ull, 696150ull, 746130ull, 750750ull, 810810ull, 824670ull, 
833910ull, 870870ull, 881790ull, 903210ull, 921690ull, 930930ull, 981750ull, 
990990ull, 1009470ull, 1067430ull, 1111110ull, 1138830ull, 1193010ull, 
1217370ull, 1231230ull, 1272810ull, 1291290ull, 1345890ull, 1351350ull, 
1360590ull, 1411410ull, 1452990ull, 1531530ull, 1591590ull, 1610070ull, 
1623930ull, 1688610ull, 1711710ull, 1771770ull, 1831830ull, 1845690ull, 
1891890ull, 2012010ull, 2072070ull, 2132130ull, 2192190ull, 2238390ull, 
2372370ull, 2492490ull, 2552550ull, 2612610ull, 2645370ull, 2672670ull, 
2709630ull, 2792790ull, 2852850ull, 2912910ull, 3028410ull, 3033030ull, 
3202290ull, 3333330ull, 3416490ull, 3453450ull, 3573570ull, 3652110ull, 
3693690ull, 3730650ull, 3873870ull, 3993990ull, 4354350ull, 4516050ull, 
4594590ull, 4654650ull, 4834830ull, 5135130ull, 5222910ull, 5555550ull, 
5615610ull, 6096090ull, 6156150ull, 6216210ull, 6276270ull, 6322470ull, 
6516510ull, 6636630ull, 6715170ull, 7417410ull, 7597590ull, 7657650ull, 
7777770ull, 7837830ull, 8128890ull, 8207430ull, 8378370ull, 8558550ull, 
8678670ull, 8978970ull, 9579570ull, 9699690ull, 10360350ull, 10720710ull, 
10840830ull, 11191950ull, 11321310ull, 11741730ull, 11981970ull, 12684210ull, 
12762750ull, 13063050ull, 13123110ull, 13783770ull, 14176470ull, 14264250ull, 
14504490ull, 14804790ull, 15405390ull, 15668730ull, 15825810ull, 16546530ull, 
17160990ull, 17687670ull, 18888870ull, 20030010ull, 20281170ull, 20930910ull, 
21111090ull, 21411390ull, 21637770ull, 21951930ull, 22972950ull, 23130030ull, 
23393370ull, 23993970ull, 24534510ull, 25555530ull, 26036010ull, 26816790ull, 
26996970ull, 27057030ull, 27606810ull, 28318290ull, 29099070ull, 29699670ull, 
30120090ull, 30240210ull, 31141110ull, 32222190ull, 32462430ull, 33183150ull, 
33663630ull, 34204170ull, 34804770ull, 35225190ull, 36246210ull, 37267230ull, 
39369330ull, 40330290ull, 42372330ull, 44414370ull, 45435390ull, 47477430ull, 
48498450ull, 49519470ull, 49639590ull, 51482970ull, 51561510ull, 52582530ull, 
53063010ull, 54624570ull, 55645590ull, 56666610ull, 58708650ull, 60090030ull, 
62792730ull, 63333270ull, 64234170ull, 64913310ull, 65615550ull, 65855790ull, 
67897830ull, 71981910ull, 74023950ull, 79129050ull, 81171090ull, 82192110ull, 
82732650ull, 85804950ull, 87297210ull, 88438350ull, 91861770ull, 94444350ull, 
100150050ull, 103633530ull, 104654550ull, 105555450ull, 105675570ull, 
106696590ull, 109759650ull, 110780670ull, 115825710ull, 116966850ull, 
118107990ull, 119969850ull, 120126930ull, 123813690ull, 126095970ull, 
129159030ull, 132222090ull, 133243110ull, 142432290ull, 144354210ull, 
145495350ull, 148918770ull, 152642490ull, 162852690ull, 164894730ull, 
169999830ull, 170600430ull, 174083910ull, 176125950ull, 182011830ull, 
184294110ull, 192462270ull, 196846650ull, 199609410ull, 203693490ull, 
205735530ull, 207777570ull, 215104890ull, 222071850ull, 223092870ull, 
237387150ull, 242492250ull, 245555310ull, 246576330ull, 248197950ull, 
249339090ull, 251681430ull, 261891630ull, 281291010ull, 293543250ull, 
300690390ull, 317026710ull, 320089770ull, 328077750ull, 339489150ull, 
340510170ull, 358888530ull, 363993630ull, 378287910ull, 380570190ull, 
387477090ull, 397687290ull, 406816410ull, 410960550ull, 417086670ull, 
433062630ull, 434444010ull, 436486050ull, 455885430ull, 457927470ull, 
458948490ull, 475284810ull, 481410930ull, 485555070ull, 494684190ull, 
497668710ull, 504894390ull, 511801290ull, 512942430ull, 514083570ull, 
528377850ull, 531990690ull, 533482950ull, 538047510ull, 547777230ull, 
551861310ull, 564293730ull, 572281710ull, 585554970ull, 591681090ull, 
606996390ull, 616786170ull, 622311690ull, 636605970ull, 648858210ull, 
649879230ull, 669278610ull, 688677990ull, 692762070ull, 695524830ull, 
695825130ull, 708077370ull, 716245530ull, 766275510ull, 786695910ull, 
805074270ull, 833662830ull, 843873030ull, 857146290ull, 863272410ull, 
902071170ull, 927596670ull, 940869930ull, 979668690ull, 999068070ull, 
1021530510ull, 1037866830ull, 1057266210ull, 1076665590ull, 1091980890ull, 
1096064970ull, 1115464350ull, 1141710570ull, 1193061870ull, 1220449230ull, 
1251260010ull, 1303332030ull, 1367656290ull, 1376845470ull, 1406455050ull, 
1444232790ull, 1456665210ull, 1503451950ull, 1514683170ull, 1542250710ull, 
1561650090ull, 1655583930ull, 1702550850ull, 1716845130ull, 1775043270ull, 
1794442650ull, 1819968150ull, 1866935070ull, 1902850950ull, 1949637690ull, 
1969037070ull, 1988436450ull, 2007835830ull, 2034082050ull, 2066033970ull, 
2085433350ull, 2104832730ull, 2124232110ull, 2172220050ull, 2279427150ull, 
2383571190ull, 2454021570ull, 2512219710ull, 2531619090ull, 2547955410ull, 
2570417850ull, 2663991330ull, 2706213510ull, 2783811030ull, 2847714870ull, 
2861408550ull, 2900207310ull, 2919606690ull, 3041108070ull, 3064591530ull, 
3094201110ull, 3191198010ull, 3229996770ull, 3275942670ull, 3307594290ull, 
3346393050ull, 3425131710ull, 3579185610ull, 3598584990ull, 3656783130ull, 
3661347690ull, 3745611870ull, 3753780030ull, 3792578790ull, 3908975070ull, 
3909996090ull, 3947773830ull, 4003929930ull, 4005971970ull, 4102968870ull, 
4186272090ull, 4219365150ull, 4238764530ull, 4374560190ull, 4426632210ull, 
4510355850ull, 4587953370ull, 4665550890ull, 4684950270ull, 4731917190ull, 
4781947170ull, 5107652550ull, 5111736630ull, 5169934770ull, 5344529190ull, 
5383327950ull, 5422126710ull, 5459904450ull, 5577321750ull, 5713117410ull, 
5788672890ull, 5907111210ull, 5965309350ull, 6023507490ull, 6101105010ull, 
6187891710ull, 6256300050ull, 6314498190ull, 6469693230ull, 6760683930ull, 
6818882070ull, 6838281450ull, 6915878970ull, 7032275250ull, 7090473390ull, 
7150713570ull, 7362064710ull, 7517259750ull, 7536659130ull, 7556058510ull, 
7594857270ull, 7643866230ull, 7808250450ull, 7924646730ull, 8118640530ull, 
8254436190ull, 8351433090ull, 8512754250ull, 8720021310ull, 9146807670ull, 
9592993410ull, 10039179150ull, 10407767370ull, 10485364890ull, 10555815270ull, 
10931550630ull, 11125544430ull, 11377736370ull, 11532931410ull, 
11797675890ull, 11823922110ull, 12095513430ull, 12328305990ull, 
12598876290ull, 12716293590ull, 12929686770ull, 13162479330ull, 
13220677470ull, 13467764310ull, 13608665070ull, 13960916970ull, 
14054850810ull, 14081097030ull, 14132448330ull, 14501036550ull, 
14641937310ull, 14714429730ull, 14908423530ull, 14923738830ull, 
14947222290ull, 15432206790ull, 15651726090ull, 15839593770ull, 
15936590670ull, 16003977990ull, 16285779510ull, 16596169590ull, 
16731965250ull, 16867760910ull, 17100553470ull, 17107700610ull, 
17158751610ull, 17178150990ull, 17624336730ull, 17740733010ull, 
18047039010ull, 18342113790ull, 18516708210ull, 18691302630ull, 
18846497670ull, 18962893950ull, 19021092090ull, 19409079690ull, 
19603073490ull, 19855265430ull, 19971661710ull, 20146256130ull, 
20534243730ull, 20747636910ull, 21174423270ull, 21349017690ull, 
21640008390ull, 21950398470ull, 22221989790ull, 22532379870ull, 
22978565610ull, 23347153830ull, 23754540810ull, 23870937090ull, 
24317122830ull, 24763308570ull, 24957302370ull, 25034899890ull, 
25209494310ull, 26160063930ull, 26548051530ull, 27285227970ull, 
27440423010ull, 28332794490ull, 28410392010ull, 28778980230ull, 
28972974030ull, 29225165970ull, 30098138070ull, 30563723190ull, 
31009908930ull, 31223302110ull, 31456094670ull, 31667445810ull, 
32348466150ull, 33240837630ull, 33376633290ull, 33687023370ull, 
34579394850ull, 34598794230ull, 35393027670ull, 35471766330ull, 
36286540290ull, 36984917970ull, 37796628870ull, 38789060310ull, 
39487437990ull, 39662032410ull, 40825995210ull, 41272180950ull, 
42397344990ull, 43600106550ull, 44143289190ull, 44725270590ull, 
44841666870ull, 45287852610ull, 45734038350ull, 47518781310ull, 
47809772010ull, 47964967050ull, 48411152790ull, 48857338530ull, 
49788508770ull, 52038836850ull, 52426824450ull, 52779076350ull, 
52873010190ull, 55550124630ull, 55627722150ull, 57664657050ull, 
57781053330ull, 58227239070ull, 58988379450ull, 59119610550ull, 
59565796290ull, 60477567150ull, 61040149170ull, 61641529950ull, 
62242910730ull, 62994381450ull, 64027653690ull, 64648433850ull, 
64920025170ull, 65812396650ull, 66103387350ull, 67150953870ull, 
67597139610ull, 68043325350ull, 71166625530ull, 72854371590ull, 
73397554230ull, 74289925710ull, 74736111450ull, 76074668670ull, 
77878811010ull, 78480191790ull, 79197968850ull, 80730519870ull, 
82321269030ull, 82767454770ull, 84106011990ull, 84668594010ull, 
86336940690ull, 89906426610ull, 90798798090ull, 92137355310ull, 
93669906330ull, 94368284010ull, 95260655490ull, 95920234410ull, 
97045398450ull, 100129899870ull, 100614884370ull, 103738184550ull, 
103796382690ull, 104630556030ull, 105522927510ull, 106415298990ull, 
107307670470ull, 108859620870ull, 109984784910ull, 110877156390ull, 
113360277030ull, 114485441070ull, 115339013790ull, 117569942490ull, 
118462313970ull, 118908499710ull, 122380988730ull, 122477985630ull, 
122924171370ull, 123816542850ull, 124708914330ull, 126862245510ull, 
130063143210ull, 130800319650ull, 131401700430ull, 133050647730ull, 
134525000610ull, 135300975810ull, 135863557830ull, 136309743570ull, 
137202115050ull, 140325415230ull, 143894901150ull, 144787272630ull, 
145233458370ull, 148240362270ull, 153710987430ull, 155495730390ull, 
156116510550ull, 156834287610ull, 157280473350ull, 161742330750ull, 
163080887970ull, 165680404890ull, 171112231290ull, 172896974250ull, 
173343159990ull, 173789345730ull, 174681717210ull, 176912645910ull, 
176932045290ull, 177358831650ull, 178251203130ull, 182266874790ull, 
183120447510ull, 186728732190ull, 192082961070ull, 197437189950ull, 
197747580030ull, 199221932910ull, 200560490130ull, 201006675870ull, 
201452861610ull, 206360904750ull, 213499876590ull, 218000532750ull, 
218563114770ull, 220192662690ull, 222869777130ull, 223762148610ull, 
224654520090ull, 226439263050ull, 228224006010ull, 228670191750ull, 
231347306190ull, 235440575370ull, 239378649510ull, 239824835250ull, 
242055763950ull, 246963807090ull, 248302364310ull, 250087107270ull, 
252318035970ull, 255887521890ull, 258564636330ull, 259010822070ull, 
260194184250ull, 262134122250ull, 265257422430ull, 269719279830ull, 
272396394270ull, 278196808890ull, 283551037770ull, 291136195350ull, 
297382795710ull, 304075581810ull, 317014968270ull, 322640788470ull, 
325046311590ull, 329954354730ull, 338431883790ull, 338878069530ull, 
342893741190ull, 354940756170ull, 355833127650ull, 357520873710ull, 
366541585410ull, 368772514110ull, 374960916330ull, 380373343350ull, 
381711900570ull, 387958500930ull, 390565164990ull, 393312729810ull, 
394651287030ull, 407590673490ull, 408036859230ull, 409841001570ull, 
420530059950ull, 421868617170ull, 426718462170ull, 429899960490ull, 
432788426070ull, 433469446410ull, 436514007930ull, 437485118070ull, 
447533996910ull, 450870690270ull, 453900056610ull, 456147321630ull, 
459348219330ull, 462161129430ull, 463363890990ull, 472287605790ull, 
478398410490ull, 484780806510ull, 485226992250ull, 487011735210ull, 
489165066390ull, 491027406870ull, 495916050630ull, 496123317690ull, 
498166378710ull, 503520607590ull, 504859164810ull, 508428650730ull, 
511105765170ull, 514481257290ull, 516553927890ull, 522900588210ull, 
530117157570ull, 531921299910ull, 536984538090ull, 539661652530ull, 
542047776270ull, 546354438630ull, 553047224730ull, 555724339170ull, 
557955267870ull, 565986611190ull, 574017954510ull, 575802697470ull, 
584241427770ull, 585172598010ull, 586064969490ull, 601681470390ull, 
602573841870ull, 612836113890ull, 615513228330ull, 618636528510ull, 
619121513010ull, 627560243310ull, 636561555630ull, 639607258290ull, 
642730558470ull, 649423344570ull, 652100459010ull, 653439016230ull, 
666378402690ull, 667716959910ull, 670840260090ull, 681102532110ull, 
685118203770ull, 692257175610ull, 698503775970ull, 700288518930ull, 
705196562070ull, 712335533910ull, 718135948530ull, 722597805930ull, 
731075334990ull, 734644820910ull, 739999049790ull, 753830807730ull, 
767662565670ull, 781494323610ull, 795772267290ull, 800680310430ull, 
821651040210ull, 833698055190ull, 834590426670ull, 847529813130ull, 
850653113310ull, 878316629190ull, 886347972510ull, 892148387130ull, 
899287358970ull, 905980145070ull, 912226745430ull, 947475418890ull, 
961307176830ull, 963984291270ull, 967922365410ull, 975138934770ull, 
976923677730ull, 1002802450650ull, 1015295651370ull, 1015741837110ull, 
1028681223570ull, 1030465966530ull, 1044297724470ull, 1054559996490ull, 
1064822268510ull, 1072562621130ull, 1080438769410ull, 1085792998290ull, 
1099624756230ull, 1124882748990ull, 1145135701710ull, 1163875502790ull, 
1179938189430ull, 1183953861090ull, 1196893247550ull, 1224110577690ull, 
1229523004710ull, 1265605851510ull, 1279437609450ull, 1280155386510ull, 
1289699881470ull, 1300408339230ull, 1312455354210ull, 1326287112150ull, 
1352612070810ull, 1378044657990ull, 1386483388290ull, 1390091672970ull, 
1390984044450ull, 1403923430910ull, 1416862817370ull, 1417755188850ull, 
1454342419530ull, 1461035205630ull, 1473082220610ull, 1486913978550ull, 
1510561822770ull, 1514577494430ull, 1520377909050ull, 1525285952190ull, 
1533317295510ull, 1610953614270ull, 1613203942350ull, 1625231557950ull, 
1639063315890ull, 1675650546570ull, 1692159418950ull, 1714468705950ull, 
1722053863530ull, 1727408092410ull, 1774703780850ull, 1787604368550ull, 
1791212653230ull, 1805044411170ull, 1832707927050ull, 1846539684990ull, 
1856801957010ull, 1874804581650ull, 1882680729930ull, 1908559502850ull, 
1939792504650ull, 1947377662230ull, 1960317048690ull, 1966563649050ull, 
1973256435150ull, 1984857264390ull, 1999135208070ull, 2012520780270ull, 
2040184296150ull, 2049205007850ull, 2076771526830ull, 2081679569970ull, 
2095511327910ull, 2109343085850ull, 2115589686210ull, 2128529072670ull, 
2133592310850ull, 2137006601730ull, 2149499802450ull, 2154407845590ull, 
2167347232050ull, 2187425590350ull, 2193226004970ull, 2206165391430ull, 
2258485519290ull, 2275324181130ull, 2296741096650ull, 2302987697010ull, 
2316819454950ull, 2361438028950ull, 2369023186530ull, 2387316801870ull, 
2400256188330ull, 2484585293190ull, 2502646115970ull, 2503771280010ull, 
2551959339930ull, 2555528825850ull, 2565791097870ull, 2607286371690ull, 
2633165144610ull, 2671983303990ull, 2676445161390ull, 2715709506510ull, 
2736680236290ull, 2753189108670ull, 2762559009210ull, 2814762740790ull, 
2856258014610ull, 2903767096230ull, 2917831646730ull, 2925416804310ull, 
2953080320190ull, 3008407351950ull, 3034286124870ull, 3045886954110ull, 
3060164897790ull, 3086043670710ull, 3111922443630ull, 3119061415470ull, 
3194466805530ull, 3215437535310ull, 3217687863390ull, 3243547236930ull, 
3271210752810ull, 3298874268690ull, 3306013240530ull, 3326537784570ull, 
3344831399910ull, 3409528332210ull, 3435407105130ull, 3448346491590ull, 
3491626508370ull, 3539814568290ull, 3549048673170ull, 3551861583270ull, 
3575509427490ull, 3577740356190ull, 3590679742650ull, 3616558515570ull, 
3672331733070ull, 3686163491010ull, 3722750721690ull, 3771831153090ull, 
3796817554530ull, 3810649312470ull, 3838312828350ull, 3865976344230ull, 
3901225017690ull, 3904348317870ull, 3932729610810ull, 3952982563530ull, 
3978861336450ull, 4031957439510ull, 4069437041670ull, 4134133973970ull, 
4170275018910ull, 4172952133350ull, 4194330250110ull, 4198830906270ull, 
4211770292730ull, 4225602050670ull, 4250588452110ull, 4253265566550ull, 
4267543510230ull, 4326440027910ull, 4341164157330ull, 4350087872130ull
};


/* All the prime factors that can appear in eulerphi(P) */
const unsigned long phiPfactors[] = {2UL, 3UL, 5UL, 7UL, 11UL, 13UL, 
				     17UL, 19UL, 23UL, 29UL, 31UL, 37UL,
				     41UL, 43UL, 47UL, 53UL};

/* returns Euler's totient phi function */
uint64_t
eulerphi64 (uint64_t n)
{
  uint64_t phi = 1, p;

  for (p = 2UL; p * p <= n; p += 2)
    {
      if (n % p == 0)
	{
	  phi *= p - 1;
	  n /= p;
	  while (n % p == 0)
	    {
	      phi *= p;
	      n /= p;
	    }
	}

      if (p == 2UL)
	p--;
    }

  /* now n is prime or 1 */
  return (n == 1) ? phi : phi * (n - 1);
}

/* Approximate amount of memory in bytes each coefficient in an NTT takes 
   so that NTT can do transforms up to length lmax with modulus, or
   with 2*modulus if twice != 0 */
static size_t
ntt_coeff_mem (const uint64_t lmax, const mpz_t modulus, const int twice)
{
  mpz_t t, ltmp;
  size_t n;
  
  mpz_init (t);
  mpz_init (ltmp);
  mpz_mul (t, modulus, modulus);
  mpz_set_uint64 (ltmp, lmax);
  mpz_mul (t, t, ltmp);
  if (twice)
    mpz_mul_2exp (t, t, 1UL);
  /* +4: +1 for rounding up, +3 for extra words due to ECRT */
  n = (mpz_sizeinbase (t, 2) - 1) / SP_NUMB_BITS + 4;
  mpz_clear (t);
  mpz_clear (ltmp);
  return n * sizeof (sp_t);
}

size_t
pm1fs2_memory_use (const uint64_t lmax, const mpz_t modulus, 
		   const int use_ntt)
{
  if (use_ntt)
    {
      /* We store lmax / 2 + 1 coefficients for the DCT-I of F and lmax 
	 coefficients for G in NTT ready format. Each coefficient in 
	 NTT-ready format occupies approx. 
	 ceil(log(lmax*modulus^2)/log(bits per sp_t)) + 3 words. */
      
      size_t n;
      
      n = ntt_coeff_mem (lmax, modulus, 0) * (size_t) (3 * lmax / 2 + 1);
      outputf (OUTPUT_DEVVERBOSE, "pm1fs2_memory_use: Estimated memory use "
	       "with lmax = %lu NTT is %lu bytes\n", lmax, n);
      return n;
    }
  else
    {
      /* F stores s_1/2 residues,
	 h stores s_1 mpz_t structs (residues get cloned from F)
	 g stores lmax residues, 
	 R stores lmax-s_1 residues, 
	 and tmp stores 3*lmax+list_mul_mem (lmax / 2) residues.
	 Assume s_1 is close to lmax/2.
	 Then we have 
	 lmax/4 + lmax/2 + lmax + lmax/2 + 3*lmax + list_mul_mem (lmax / 2)
	 = (5+1/4)*lmax + list_mul_mem (lmax / 2) residues, plus s_1 mpz_t.
      */
      
      size_t n;
      
      n = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);
      n *= 5 * lmax + lmax / 4 + list_mul_mem (lmax / 2);
      n += lmax / 2 * sizeof (mpz_t);
      /* Memory use due to temp space allocation in TMulKS appears to 
	 approximately triple the estimated memory use. This is hard to
	 estimate precisely, so let's go with the fudge factor of 3 here */
      n *= 3;
      outputf (OUTPUT_DEVVERBOSE, "pm1fs2_memory_use: Estimated memory use "
	       "with lmax = %lu is %lu bytes\n", lmax, n);
      return n;
    }
}

/* return the possible lmax for given memory use and modulus */

unsigned long
pm1fs2_maxlen (const size_t memory, const mpz_t modulus, const int use_ntt)
{
  if (use_ntt)
    {
      size_t n, lmax = 1;
  
      n = ntt_coeff_mem (lmax, modulus, 0);
      lmax = 1UL << ceil_log2 (memory / n / 3);
      return lmax;
    }
  else
    {
      size_t lmax, n;
      
      n = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);

      /* Guess an initial value of lmax for list_mul_mem (lmax / 2) */
      /* memory = n * 25/4 * lmax + lmax / 2 * sizeof (mpz_t); */
      /* Fudge factor of 3 for TMulKS as above */
      lmax = memory / (3 * 25 * n / 4 + 3 * sizeof (mpz_t) / 2);
      return lmax;
    }
}

size_t 
pp1fs2_memory_use (const uint64_t lmax, const mpz_t modulus, 
		   const int use_ntt, const int twopass)
{
  size_t n, m;
  
  m = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);
  if (use_ntt)
    {
      /* In one pass mode, we store h_x_ntt and h_y_ntt, each of length 
	 lmax/2(+1), and g_x_ntt and g_y_ntt, each of length lmax, all in 
	 NTT ready format. In two pass mode, we store h_x_ntt, h_y_ntt and 
	 g_x_ntt as before, plus R which is lmax - s_1 mpz_t. 
	 We assume s_1 ~= lmax/2.
      */

      n = ntt_coeff_mem (lmax, modulus, !twopass);
      if (twopass)
	return lmax * (2 * n + m / 2);
      else
	return lmax * 3 * n;
    }
  else
    {
      /* We allocate:
	 F: s_1/2 coefficients
	 fh_x, fh_y: s_1/2 coefficients
	 h_x, h_y: s_1 mpz_t's (cloned from fh_x and fh_y)
	 g_x, g_y: lmax coefficients
	 R_x, R_y: lmax - s_1 coefficients
	 tmp: 3UL * lmax + list_mul_mem (lmax / 2)
	 Assuming s_1 ~ lmax/2, that's
	 lmax/2 + 2*lmax/4 + 2*lmax + 2*lmax/2 * 3*lmax + 
           list_mul_mem (lmax / 2) =
	 7 + list_mul_mem (lmax / 2) coefficients and lmax mpz_t.
       */
      
      n = m * (7 * lmax + list_mul_mem (lmax / 2));
      n += lmax * sizeof (mpz_t);
      n = 5 * n / 2; /* A fudge factor again */
      return n;
    }
}

unsigned long 
pp1fs2_maxlen (const size_t memory, const mpz_t modulus, const int use_ntt, 
	       const int twopass)
{
  size_t n, m;
  
  m = mpz_size (modulus) * sizeof (mp_limb_t) + sizeof (mpz_t);
  if (use_ntt)
    {
      n = ntt_coeff_mem (1, modulus, !twopass);
      if (twopass)
	n = memory / (2 * n + m / 2);
      else
	n = memory / (3 * n);
      return 1UL << (ceil_log2 (n / 2)); /* Rounded down to power of 2 */
    }
  else
    {
      return memory / 5 / (m * 8 + sizeof (mpz_t)) * 2;
    }
}


/* Test if for given P, nr, B2min and B2 we can choose an m_1 so that the 
   stage 2 interval [B2min, B2] is covered. The effective B2min and B2
   are stored in effB2min and effB2 */

static int
test_P (const mpz_t B2min, const mpz_t B2, mpz_t m_1, const uint64_t P, 
	const uint64_t nr, mpz_t effB2min, mpz_t effB2)
{
  mpz_t m, Ptmp, nrtmp;
  /* We need B2min >= max(S_1 + S_2) + (m_1 - 1) * P + 1, or
     B2min - max(S_1 + S_2) - 1 >= m_1*P - P, or
     (B2min - max(S_1 + S_2) - 1)/P + 1 >= m_1
     Choose m_1 accordingly */
  
  mpz_init (m);
  mpz_init (Ptmp);
  mpz_init (nrtmp);
  sets_max (m, P);
  mpz_mul_2exp (m, m, 1UL); /* m = 2*max(S_1 + S_2) */

  mpz_sub (m_1, B2min, m);
  mpz_sub_ui (m_1, m_1, 1UL); /* m_1 = B2min - max(S_1 + S_2) - 1 */
  mpz_set_uint64 (Ptmp, P);
  mpz_fdiv_q (m_1, m_1, Ptmp);
  mpz_add_ui (m_1, m_1, 1);
  
  /* Compute effB2min = max(S_1 + S_2) + (m_1 - 1)*P + 1 */
  
  mpz_sub_ui (effB2min, m_1, 1UL);
  mpz_mul (effB2min, effB2min, Ptmp);
  mpz_add (effB2min, effB2min, m);
  mpz_add_ui (effB2min, effB2min, 1UL);
  ASSERT_ALWAYS (mpz_cmp (effB2min, B2min) <= 0);

  /* Compute the smallest value coprime to P at the high end of the stage 2
     interval that will not be covered: 
     min(S_1 + S_2) + (m_1 + nr)*P. 
     We assume min(S_1 + S_2) = -max(S_1 + S_2) */
  mpz_set_uint64 (nrtmp, nr);
  mpz_add (effB2, m_1, nrtmp);
  mpz_mul (effB2, effB2, Ptmp);
  mpz_sub (effB2, effB2, m);

  /* The effective B2 values is that value, minus 1 */
  mpz_sub_ui (effB2, effB2, 1UL);

  mpz_clear (m);
  mpz_clear (Ptmp);
  mpz_clear (nrtmp);
  return (mpz_cmp (B2, effB2) <= 0);
}


static void
factor_phiP (int *exponents, const uint64_t phiP)
{
    const int nrprimes = sizeof (phiPfactors) / sizeof (unsigned long);
    uint64_t cofactor = phiP;
    int i;
    
    ASSERT_ALWAYS (phiP > 0);

    for (i = 0; i < nrprimes; i++)
	for (exponents[i] = 0; cofactor % phiPfactors[i] == 0UL; exponents[i]++)
	    cofactor /= phiPfactors[i];

    ASSERT_ALWAYS (cofactor == 1UL);
}


static unsigned long 
pow_ul (const unsigned long b, const unsigned int e)
{
    unsigned long r = 1UL;
    unsigned int i;

    for (i = 0; i < e; i++)
	r *= b;

    return r;
}

static uint64_t
absdiff_ul (uint64_t a, uint64_t b)
{
    return (a > b) ? a - b : b - a;
}

/* Choose s_1 so that s_1 * s_2 = phiP, s_1 is positive and even, 
   s_2 >= min_s2 and s_2 is minimal and abs(s_1 - l) is minimal 
   under those conditions. If use_ntt == 1, we require s_1 < l.
   Returns 0 if no such choice is possible */

static uint64_t
choose_s_1 (const uint64_t phiP, const uint64_t min_s2,
	    const uint64_t l, const int use_ntt)
{
  const int nrprimes = sizeof (phiPfactors) / sizeof (unsigned long);
  /* Using [nrprimes] here makes the compiler complain about variable-sized
     arrays */
  int phiPexponents[sizeof (phiPfactors) / sizeof (unsigned long)], 
    exponents[sizeof (phiPfactors) / sizeof (unsigned long)];
  uint64_t s_1 = 0UL, s_2 = 0UL, trys_1;
  int i;

  ASSERT_ALWAYS (phiP > 0 && phiP % 2 == 0);

  /* We want only even s_1. We divide one 2 out of phiP here... */
  factor_phiP (phiPexponents, phiP / 2);
  for (i = 0; i < nrprimes; i++)
      exponents[i] = 0;

  do {
      trys_1 = 2; /* ... and add a 2 here */
      for (i = 0; i < nrprimes; i++)
	  trys_1 *= pow_ul (phiPfactors[i], exponents[i]);
#if 0
      printf ("choose_s_1: Trying trys_1 = %" PRId64 "\n", trys_1);
#endif
      /* See if it satisfies all the required conditions and is an 
	 improvement over the previous choice */
      if (phiP / trys_1 >= min_s2 && 
	  (s_2 == 0 || phiP / trys_1 < s_2) && 
	  absdiff_ul (trys_1, l) < absdiff_ul (s_1, l) &&
	  (use_ntt == 0 || trys_1 < l))
      {
#if 0
	  printf ("choose_s_1: New best s_1 for "
	          "phiP = %" PRId64 ", min_s2 = %" PRId64 
		  ", l = %" PRId64 " : %" PRId64 "\n", 
                  phiP, min_s2, l, trys_1);
#endif
	  s_1 = trys_1;
      }
      for (i = 0; i < nrprimes; i++)
      {
	  if (++(exponents[i]) <= phiPexponents[i])
	      break;
	  exponents[i] = 0;
      }
  } while (i < nrprimes);

  return s_1;
}


/* Approximate cost of stage 2. Cost with and without ntt are not 
   comparable. We have l > s_1 and s_1 * s_2 = eulerphi(P), hence
   s_2*l > eulerphi(P) and so cost (s_2, l) > eulerphi(P) for all P */
static uint64_t
est_cost (const uint64_t s_2, const uint64_t l, const int use_ntt,
          const int method)
{
  if (method == ECM_PM1)
    {
      /* The time for building f, h and DCT-I of h seems to be about 
         7/6 of the time of computing g, h*g and gcd with NTT, and 
         3/2 of the time of computing g, h*g and gcd without NTT */

      if (use_ntt)
        return (7 * l) / 6 + s_2 * l;
      else
        return (3 * l) / 2 + s_2 * l;
    }
  else if (method == ECM_PP1)
    {
      /* Building f is the same, building h and its forward transform is
         twice about as expensive as for P-1. Each multi-point evaluation
         is twice as expensive as for P-1.
         FIXME: The estimate for NTT assumes the "one-pass" variant, in 
         "two-pass" the multipoint evaluations are slower, so the optimum 
         shifts towards smaller s_2 some more */
      if (use_ntt)
        return (4 * l) / 5 + s_2 * l;
      else
        return (3 * l) / 4 + s_2 * l;
    }
  else
    abort (); /* Invalid value for method */
}

/* Choose P so that a stage 2 range from B2min to B2 can be covered with
   multipoint evaluations, each using a convolution of length at most lmax. 
   The parameters for stage 2 are stored in finalparams, the final effective
   B2min and B2 values in final_B2min and final_B2, respecively. Each of these
   may be NULL, in which case the value is not stored. It is permissible
   to let B2min and final_B2min, or B2 and final_B2 point at the same mpz_t. */

long
choose_P (const mpz_t B2min, const mpz_t B2, const uint64_t lmax,
	  const uint64_t min_s2, faststage2_param_t *finalparams, 
	  mpz_t final_B2min, mpz_t final_B2, const int use_ntt, 
	  const int method)
{
  /* Let S_1 + S_2 == (Z/PZ)* (mod P).

     Let F(x) = \prod_{k_1 \in S_1} (x - b_1^{k_1}).

     If we evaluate F(b_1^{k_2 + m P}) for all k_2 \in S_2 with 
     m_1 <= m < m_1+nr, we test all exponents k_2 + m P - k_1.
     The largest value coprime to P at the low end of the stage 2 interval 
     *not* covered will be 
       max(S_2) + (m_1 - 1)*P - min(S_1).
     The smallest value at the high end not covered will be
       min(S_2) + (m_1 + nr)*P - max(S_1).
     Assume S_1 and S_2 are symmetric around 0, so that max(S_1) = -min(S_1).
     Then the largest ... is:
       max(S_1) + max(S_2) + (m_1 - 1)*P
     The smallest ... is:
       -(max(S_1) + max(S_2)) + (m_1 + nr)*P
     The effective B2min = (max(S_1) + max(S_2)) + (m_1 - 1)*P + 1
     The effective B2max = -(max(S_1) + max(S_2)) + (m_1 + nr)*P - 1

     Then the difference effB2max - effB2min =
       -2*(max(S_1) + max(S_2)) + P*(nr + 1) - 2

     We obviously require B2max - B2min <= nr*P
     Since nr < lmax, B2max - B2min <= lmax*P or
     P >= ceil((B2max - B2min)/lmax)

     Hence we are looking for an odd P with s_1 * s_2 = eulerphi(P) so that
     s_1 ~= lmax / 2 and the whole stage 2 interval is covered. s_2 should 
     be small, as long as s_1 is small enough.
  */

  mpz_t B2l, m_1, effB2min, tryeffB2, effB2, lmin, Ptmp;
  /* The best parameters found so far, P == 0 means that no suitable P
     has been found yet: */
  uint64_t P = 0, s_1 = 0, s_2 = 0, cost = 0, l = 0;
  unsigned int i;
  const unsigned int Pvalues_len = sizeof (Pvalues) / sizeof (uint64_t);
  int r;

  outputf (OUTPUT_TRACE, "choose_P(B2min = %Zd, B2 = %Zd, ", B2min, B2);
  outputf (OUTPUT_TRACE, "lmax = %" PRId64 ", min_s2 = %" PRId64 
          ", use_ntt = %d, method = %d\n", lmax, min_s2, use_ntt, method);

  if (mpz_cmp (B2, B2min) < 0)
    return 0L;

  /* If we use the NTT, we allow only power-of-two transform lengths.
     In that case, the code below assumes that lmax is a power of two.
     If that is not the case, print error and return. */
  if (use_ntt && (lmax & (lmax - 1UL)) != 0)
    {
      outputf (OUTPUT_ERROR, 
               "choose_P: Error, lmax = %lu is not a power of two\n", lmax);
      return ECM_ERROR;
    }
  
  mpz_init (effB2);
  mpz_init (tryeffB2);
  mpz_init (effB2min);
  mpz_init (B2l);
  mpz_init (m_1);
  mpz_init (lmin);
  mpz_init (Ptmp);
  
  mpz_sub (B2l, B2, B2min);
  mpz_add_ui (B2l, B2l, 1UL); /* +1 due to closed interval */
  
  /* For each candidate P, check if [B2min, B2] can be covered at all,
     and if so, what the best parameters (minimizing the cost, maximizing 
     effB2) are. If they are better than the best parameters for the best P 
     so far, remember them. */

  for (i = 0 ; i < Pvalues_len; i++)
    {
      uint64_t tryP, tryphiP, trys_1, trys_2, trycost, tryl;
      
      tryP = Pvalues[i];
      tryphiP = eulerphi64 (tryP);
      
      outputf (OUTPUT_TRACE, 
	       "choose_P: trying P = %" PRId64 ", eulerphi(P) = "
	       "%" PRId64 "\n", tryP, tryphiP);
      
      /* If we have a good P already and this tryphiP >= cost, then 
	 there's no hope for this tryP, since cost(s_2, l) > eulerphi(P) */
      if (P != 0 && tryphiP >= cost)
	{
	  outputf (OUTPUT_TRACE, 
		   "choose_P: tryphiP > cost = %" PRId64 
		   ", this P is too large\n", cost);
	  continue;
	}
      
      /* We have nr < l and effB2-effB2min <= nr*P. Hence we need 
	 l >= B2l/P */
      mpz_set_uint64(Ptmp, tryP);
      mpz_cdiv_q (lmin, B2l, Ptmp);
      outputf (OUTPUT_TRACE, "choose_P: lmin = %Zd for ", lmin);
      outputf (OUTPUT_TRACE, "P = %" PRId64 "\n", tryP);
      if (mpz_cmp_ui (lmin, lmax) > 0)
	{
	  outputf (OUTPUT_TRACE, 
		   "choose_P: lmin > lmax, this P is too small\n");
	  continue;
	}
      
      /* Try all possible transform lengths and store parameters in 
	 P, s_1, s_2, l if they are better than the previously best ones */
       
      /* Keep reducing tryl to find best parameters. For NTT, we only have 
	 power of 2 lengths so far, so we can simply divide by 2. 
	 For non-NTT, we have arbitrary transform lengths so we can decrease 
	 in smaller steps... let's say by, umm, 25% each time? */
      for (tryl = lmax; mpz_cmp_ui (lmin, tryl) <= 0;
	   tryl = (use_ntt) ? tryl / 2 : 3 * tryl / 4)
	{
	  trys_1 = choose_s_1 (tryphiP, min_s2, tryl / 2, use_ntt);
	  if (trys_1 == 0)
	    {
	      outputf (OUTPUT_TRACE, 
		       "choose_P: could not choose s_1 for "
		       "P = %" PRId64 ", l = %lu\n", tryP, tryl);
	      continue;
	    }
	  ASSERT (tryphiP % trys_1 == 0UL);
	  trys_2 = tryphiP / trys_1;
	  outputf (OUTPUT_TRACE, 
	      	"choose_P: chose s_1 = %" PRId64 ", k = s_2 = %" PRId64 
		   " for P = %" PRId64 ", l = %" PRId64 "\n", 
		   trys_1, trys_2, tryP, tryl);
	  
	  if (test_P (B2min, B2, m_1, tryP, tryl - trys_1, effB2min, tryeffB2))
	    {
	      /* can't mix 64-bit types and mpz_t on win32 for some reason */
	      outputf (OUTPUT_TRACE, 
		       "choose_P: P = %" PRId64 ", l = %" PRId64 ", "
		       "s_1 = %" PRId64 ", k = s_2 = %" PRId64 " works, ",
		       tryP, tryl, trys_1, trys_2);
	      outputf (OUTPUT_TRACE, 
		       "m_1 = %Zd, effB2min = %Zd, effB2 = %Zd\n",
		       m_1, effB2min, tryeffB2);
	      /* We use these parameters if we 
		 1. didn't have any suitable ones yet, or 
		 2. these cover [B2min, B2] and are cheaper than the best 
                    ones so far, or 
		 3. they are as expensive but reach greater effB2. */
	      trycost = est_cost (trys_2, tryl, use_ntt, method);
	      ASSERT (tryphiP < trycost);
	      if (P == 0 || trycost < cost ||
		  (trycost == cost && mpz_cmp (tryeffB2, effB2) > 0))
		{
		  outputf (OUTPUT_TRACE, 
			   "choose_P: and is the new "
			   "optimum (cost = %" PRId64 ")\n", trycost);
		  P = tryP;
		  s_1 = trys_1;
		  s_2 = trys_2;
		  l = tryl;
		  cost = trycost;
		  mpz_set (effB2, tryeffB2);
		}
	    }
	}
  }
  
  if (P != 0) /* If we found a suitable P */
    {
      /* Compute m_1, effB2min, effB2 again */
      r = test_P (B2min, B2, m_1, P, l - s_1, effB2min, effB2);
      ASSERT_ALWAYS(r != 0);
      if (finalparams != NULL)
	{
	  finalparams->P = P;
	  finalparams->s_1 = s_1;
	  finalparams->s_2 = s_2;
	  finalparams->l = l;
	  mpz_set (finalparams->m_1, m_1);
	}
      if (final_B2min != NULL)
	mpz_set (final_B2min, effB2min);
      if (final_B2 != NULL)
	mpz_set (final_B2, effB2);
    }
  
  mpz_clear (effB2);
  mpz_clear (tryeffB2);
  mpz_clear (effB2min);
  mpz_clear (B2l);
  mpz_clear (m_1);
  mpz_clear (lmin);
  mpz_clear (Ptmp);

  return (P != 0) ? (long) P : ECM_ERROR;
}
