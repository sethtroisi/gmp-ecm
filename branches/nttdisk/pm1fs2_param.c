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

const uint64_t Pvalues[] = {
3ull, 5ull, 9ull, 15ull, 21ull, 17ull, 27ull, 33ull, 45ull, 
51ull, 63ull, 75ull, 105ull, 99ull, 135ull, 165ull, 195ull, 189ull, 
231ull, 255ull, 315ull, 345ull, 357ull, 375ull, 405ull, 435ull, 
525ull, 585ull, 615ull, 735ull, 765ull, 825ull, 945ull, 1155ull, 
1065ull, 1365ull, 1305ull, 1335ull, 1575ull, 1785ull, 1995ull, 
2145ull, 2205ull, 2415ull, 2625ull, 2805ull, 3045ull, 3465ull, 
3675ull, 4095ull, 4305ull, 4515ull, 4725ull, 4785ull, 5355ull, 
5775ull, 5985ull, 5865ull, 6825ull, 7245ull, 8085ull,

8925ull, 9555ull, 9975ull, 10395ull, 12285ull, 12495ull, 12705ull,
13965ull, 15015ull, 16065ull,

17325ull, 17745ull, 19635ull, 20475ull, 21945ull, 23205ull, 24255ull,
25935ull, 26565ull, 26775ull, 28665ull, 28875ull, 30345ull, 31185ull,
31395ull,

33495ull, 33915ull, 35805ull, 39585ull, 41055ull, 42315ull, 42735ull,
45045ull, 47355ull, 49665ull, 50505ull, 54285ull, 55965ull, 58905ull,
61215ull,

65835ull, 68145ull, 69615ull, 70455ull, 75075ull, 77385ull, 77805ull,
79695ull, 94185ull, 98175ull, 100485ull, 101745ull, 105105ull, 107415ull,
109725ull, 116025ull, 118755ull, 126945ull, 128205ull, 129675ull,

132825ull, 135135ull, 137445ull, 153615ull, 162435ull, 165165ull, 167475ull,
176715ull, 181545ull, 185955ull, 195195ull, 197505ull, 208845ull, 215985ull,
225225ull, 233415ull, 234465ull, 239085ull, 241395ull, 255255ull,

285285ull, 294525ull, 315315ull, 329175ull, 333795ull, 345345ull, 348075ull,
373065ull, 375375ull, 405405ull, 412335ull, 416955ull, 435435ull, 440895ull,
451605ull, 460845ull, 465465ull, 490875ull, 495495ull, 504735ull,

533715ull, 555555ull, 569415ull, 596505ull, 608685ull, 615615ull, 636405ull,
645645ull, 672945ull, 675675ull, 680295ull, 705705ull, 726495ull, 765765ull,
795795ull, 805035ull, 811965ull, 844305ull, 855855ull, 885885ull, 915915ull,
922845ull, 945945ull, 1006005ull, 1036035ull,

1066065ull, 1096095ull, 1119195ull, 1186185ull, 1246245ull, 1276275ull,
1306305ull, 1322685ull, 1336335ull, 1354815ull, 1396395ull, 1426425ull,
1456455ull, 1514205ull, 1516515ull, 1601145ull, 1666665ull, 1708245ull,
1726725ull, 1786785ull, 1826055ull, 1846845ull, 1865325ull, 1936935ull,
1996995ull,

2177175ull, 2258025ull, 2297295ull, 2327325ull, 2417415ull, 2567565ull,
2611455ull, 2777775ull, 2807805ull, 3048045ull, 3078075ull, 3108105ull,
3138135ull, 3161235ull, 3258255ull, 3318315ull, 3357585ull, 3708705ull,
3798795ull, 3828825ull, 3888885ull, 3918915ull, 4064445ull, 4103715ull,
4189185ull,

4279275ull, 4339335ull, 4489485ull, 4789785ull, 4849845ull, 5180175ull,
5360355ull, 5420415ull, 5595975ull, 5660655ull, 5870865ull, 5990985ull,
6342105ull, 6381375ull, 6531525ull, 6561555ull, 6891885ull, 7088235ull,
7132125ull, 7252245ull, 7402395ull, 7702695ull, 7834365ull, 7912905ull,
8273265ull,

8580495ull, 8843835ull, 9444435ull, 10015005ull, 10140585ull, 10465455ull,
10555545ull, 10705695ull, 10818885ull, 10975965ull, 11486475ull, 11565015ull,
11696685ull, 11996985ull, 12267255ull, 12777765ull, 13018005ull, 13408395ull,
13498485ull, 13528515ull, 13803405ull, 14159145ull, 14549535ull, 14849835ull,
15060045ull, 15120105ull, 15570555ull, 16111095ull, 16231215ull, 16591575ull,

16831815ull, 17102085ull, 17402385ull, 17612595ull, 18123105ull, 18633615ull,
19684665ull, 20165145ull, 21186165ull, 22207185ull, 22717695ull, 23738715ull,
24249225ull, 24759735ull, 24819795ull, 25741485ull, 25780755ull, 26291265ull,
26531505ull, 27312285ull, 27822795ull, 28333305ull, 29354325ull, 30045015ull,
31396365ull, 31666635ull, 32117085ull, 32456655ull, 32807775ull, 32927895ull,

33948915ull, 35990955ull, 37011975ull, 39564525ull, 40585545ull, 41096055ull,
41366325ull, 42902475ull, 43648605ull, 44219175ull, 45930885ull, 47222175ull,
50075025ull, 51816765ull, 52327275ull, 52777725ull, 52837785ull, 53348295ull,
54879825ull, 55390335ull, 57912855ull, 58483425ull, 59053995ull, 59984925ull,
60063465ull, 61906845ull, 63047985ull, 64579515ull, 66111045ull, 66621555ull,

71216145ull, 72177105ull, 72747675ull, 74459385ull, 76321245ull, 81426345ull,
82447365ull, 84999915ull, 85300215ull, 87041955ull, 88062975ull, 91005915ull,
92147055ull, 96231135ull, 98423325ull, 99804705ull, 101846745ull, 102867765ull,
103888785ull, 107552445ull, 111035925ull, 111546435ull, 118693575ull,
121246125ull, 122777655ull, 123288165ull, 124098975ull, 124669545ull,
125840715ull, 130945815ull,

140645505ull, 146771625ull, 150345195ull, 158513355ull, 160044885ull,
164038875ull, 169744575ull, 170255085ull, 179444265ull, 181996815ull,
189143955ull, 190285095ull, 193738545ull, 198843645ull, 203408205ull,
205480275ull, 208543335ull, 216531315ull, 217222005ull, 218243025ull,
227942715ull, 228963735ull, 229474245ull, 237642405ull, 240705465ull,
242777535ull, 247342095ull, 248834355ull, 252447195ull, 255900645ull,
256471215ull, 257041785ull, 264188925ull, 265995345ull, 266741475ull,

269023755ull, 273888615ull, 275930655ull, 282146865ull, 286140855ull,
292777485ull, 295840545ull, 303498195ull, 308393085ull, 311155845ull,
318302985ull, 324429105ull, 324939615ull, 334639305ull, 344338995ull,
346381035ull, 347762415ull, 347912565ull, 354038685ull, 358122765ull,
383137755ull, 393347955ull, 402537135ull, 416831415ull, 421936515ull,
428573145ull, 431636205ull, 451035585ull, 463798335ull, 470434965ull,
489834345ull, 499534035ull, 510765255ull, 518933415ull, 528633105ull,

538332795ull, 545990445ull, 548032485ull, 557732175ull, 570855285ull,
596530935ull, 610224615ull, 625630005ull, 651666015ull, 683828145ull,
688422735ull, 703227525ull, 722116395ull, 728332605ull, 751725975ull,
757341585ull, 771125355ull, 780825045ull, 827791965ull, 851275425ull,
858422565ull, 887521635ull, 897221325ull, 909984075ull, 933467535ull,
951425475ull, 974818845ull, 984518535ull, 994218225ull, 1003917915ull,
1017041025ull, 1033016985ull, 1042716675ull, 1052416365ull, 1062116055ull,

1086110025ull, 1139713575ull, 1191785595ull, 1227010785ull, 1256109855ull,
1265809545ull, 1273977705ull, 1285208925ull, 1331995665ull, 1353106755ull,
1391905515ull, 1423857435ull, 1430704275ull, 1450103655ull, 1459803345ull,
1520554035ull, 1532295765ull, 1547100555ull, 1595599005ull, 1614998385ull,
1637971335ull, 1653797145ull, 1673196525ull, 1712565855ull, 1789592805ull,
1799292495ull, 1828391565ull, 1830673845ull, 1872805935ull, 1876890015ull,
1896289395ull, 1954487535ull, 1954998045ull, 1973886915ull, 2001964965ull,
2002985985ull, 2051484435ull, 2093136045ull, 2109682575ull, 2119382265ull,

2187280095ull, 2213316105ull, 2255177925ull, 2293976685ull, 2332775445ull,
2342475135ull, 2365958595ull, 2390973585ull, 2553826275ull, 2555868315ull,
2584967385ull, 2672264595ull, 2691663975ull, 2711063355ull, 2729952225ull,
2788660875ull, 2856558705ull, 2894336445ull, 2953555605ull, 2982654675ull,
3011753745ull, 3050552505ull, 3093945855ull, 3128150025ull, 3157249095ull,
3234846615ull, 3380341965ull, 3409441035ull, 3419140725ull, 3457939485ull,
3516137625ull, 3545236695ull, 3575356785ull, 3681032355ull, 3758629875ull,
3768329565ull, 3778029255ull, 3797428635ull, 3821933115ull, 3904125225ull,
3962323365ull, 4059320265ull, 4127218095ull, 4175716545ull, 4256377125ull,

4360010655ull, 4573403835ull, 4796496705ull, 5019589575ull, 5203883685ull,
5242682445ull, 5277907635ull, 5465775315ull, 5562772215ull, 5688868185ull,
5766465705ull, 5898837945ull, 5911961055ull, 6047756715ull, 6164152995ull,
6299438145ull, 6358146795ull, 6464843385ull, 6581239665ull, 6610338735ull,
6733882155ull, 6804332535ull, 6980458485ull, 7027425405ull, 7040548515ull,
7066224165ull, 7250518275ull, 7320968655ull, 7357214865ull, 7454211765ull,
7461869415ull, 7473611145ull, 7716103395ull, 7825863045ull, 7919796885ull,
7968295335ull, 8001988995ull, 8142889755ull, 8298084795ull, 8365982625ull,
8433880455ull, 8550276735ull, 8553850305ull, 8579375805ull, 8589075495ull,

8812168365ull, 8870366505ull, 9023519505ull, 9171056895ull, 9258354105ull,
9345651315ull, 9423248835ull, 9481446975ull, 9510546045ull, 9704539845ull,
9801536745ull, 9927632715ull, 9985830855ull, 10073128065ull, 10267121865ull,
10373818455ull, 10587211635ull, 10674508845ull, 10820004195ull, 10975199235ull,
11110994895ull, 11266189935ull, 11489282805ull, 11673576915ull, 11877270405ull,
11935468545ull, 12158561415ull, 12381654285ull, 12478651185ull, 12517449945ull,
12604747155ull, 13080031965ull, 13274025765ull, 13642613985ull, 13720211505ull,
14166397245ull, 14205196005ull, 14389490115ull, 14486487015ull, 14612582985ull,
15049069035ull, 15281861595ull, 15504954465ull, 15611651055ull, 15728047335ull,
15833722905ull, 16174233075ull, 16620418815ull, 16688316645ull, 16843511685ull,

17289697425ull, 17299397115ull, 17696513835ull, 17735883165ull, 18143270145ull,
18492458985ull, 18898314435ull, 19394530155ull, 19743718995ull, 19831016205ull,
20412997605ull, 20636090475ull, 21198672495ull, 21800053275ull, 22071644595ull,
22362635295ull, 22420833435ull, 22643926305ull, 22867019175ull, 23759390655ull,
23904886005ull, 23982483525ull, 24205576395ull, 24428669265ull, 24894254385ull,
26019418425ull, 26213412225ull, 26389538175ull, 26436505095ull, 27775062315ull,
27813861075ull, 28832328525ull, 28890526665ull, 29113619535ull, 29494189725ull,
29559805275ull, 29782898145ull, 30238783575ull, 30520074585ull, 30820764975ull,
31121455365ull, 31497190725ull, 32013826845ull, 32324216925ull, 32460012585ull,
32906198325ull, 33051693675ull, 33575476935ull, 33798569805ull, 34021662675ull,

35583312765ull, 36427185795ull, 36698777115ull, 37144962855ull, 37368055725ull,
38037334335ull, 38939405505ull, 39240095895ull, 39598984425ull, 40365259935ull,
41160634515ull, 41383727385ull, 42053005995ull, 42334297005ull, 43168470345ull,
44953213305ull, 45399399045ull, 46068677655ull, 46834953165ull, 47184142005ull,
47630327745ull, 47960117205ull, 48522699225ull, 50064949935ull, 50307442185ull,
51869092275ull, 51898191345ull, 52315278015ull, 52761463755ull, 53207649495ull,
53653835235ull, 54429810435ull, 54992392455ull, 55438578195ull, 56680138515ull,
57242720535ull, 57669506895ull, 58784971245ull, 59231156985ull, 59454249855ull,
61190494365ull, 61238992815ull, 61462085685ull, 61908271425ull, 62354457165ull,
63431122755ull, 65031571605ull, 65400159825ull, 65700850215ull, 66525323865ull,
67262500305ull, 67650487905ull, 67931778915ull, 68154871785ull, 68601057525ull,

70162707615ull, 71947450575ull, 72393636315ull, 72616729185ull, 74120181135ull,
76855493715ull, 77747865195ull, 78058255275ull, 78417143805ull, 78640236675ull,
80871165375ull, 81540443985ull, 82840202445ull, 85556115645ull, 86448487125ull,
86671579995ull, 86894672865ull, 87340858605ull, 88456322955ull, 88466022645ull,
88679415825ull, 89125601565ull, 91133437395ull, 91560223755ull, 93364366095ull,
96041480535ull, 98718594975ull, 98873790015ull, 99610966455ull, 100280245065ull,
100503337935ull, 100726430805ull, 103180452375ull, 106749938295ull,
109000266375ull, 109281557385ull, 110096331345ull, 111434888565ull,
111881074305ull, 112327260045ull, 113219631525ull, 114112003005ull,
114335095875ull, 115673653095ull, 117720287685ull, 119689324755ull,
119912417625ull, 121027881975ull, 123481903545ull, 124151182155ull,
125043553635ull, 126159017985ull, 127943760945ull, 129282318165ull,
129505411035ull, 130097092125ull, 131067061125ull, 132628711215ull,
134859639915ull, 136198197135ull,

139098404445ull, 141775518885ull, 145568097675ull, 148691397855ull,
152037790905ull, 158507484135ull, 161320394235ull, 162523155795ull,
164977177365ull, 169215941895ull, 169439034765ull, 171446870595ull,
177470378085ull, 177916563825ull, 178760436855ull, 183270792705ull,
184386257055ull, 187480458165ull, 190186671675ull, 190855950285ull,
193979250465ull, 195282582495ull, 196656364905ull, 197325643515ull,
203795336745ull, 204018429615ull, 204920500785ull, 210265029975ull,
210934308585ull, 213359231085ull, 214949980245ull, 216394213035ull,
216734723205ull, 218257003965ull, 218742559035ull, 223766998455ull,
225435345135ull, 226950028305ull, 228073660815ull, 229674109665ull,
231080564715ull, 231681945495ull, 236143802895ull, 239199205245ull,
242390403255ull, 242613496125ull, 243505867605ull, 244582533195ull,
245513703435ull, 247958025315ull, 248061658845ull, 249083189355ull,
251760303795ull, 252429582405ull, 254214325365ull, 255552882585ull,
257240628645ull, 258276963945ull, 261450294105ull, 265058578785ull,
265960649955ull, 268492269045ull, 269830826265ull, 271023888135ull,
273177219315ull,

276523612365ull, 277862169585ull, 278977633935ull, 282993305595ull,
287008977255ull, 287901348735ull, 292120713885ull, 292586299005ull,
293032484745ull, 300840735195ull, 301286920935ull, 306418056945ull,
307756614165ull, 309318264255ull, 309560756505ull, 313780121655ull,
318280777815ull, 319803629145ull, 321365279235ull, 324711672285ull,
326050229505ull, 326719508115ull, 333189201345ull, 333858479955ull,
335420130045ull, 340551266055ull, 342559101885ull, 346128587805ull,
349251887985ull, 350144259465ull, 352598281035ull, 356167766955ull,
359067974265ull, 361298902965ull, 365537667495ull, 367322410455ull,
369999524895ull, 376915403865ull, 383831282835ull, 390747161805ull,
397886133645ull, 400340155215ull, 410825520105ull, 416849027595ull,
417295213335ull, 423764906565ull, 425326556655ull, 439158314595ull,
443173986255ull, 446074193565ull, 449643679485ull, 452990072535ull,
456113372715ull, 473737709445ull, 480653588415ull, 481992145635ull,
483961182705ull, 487569467385ull, 488461838865ull, 501401225325ull,
507647825685ull, 507870918555ull, 514340611785ull, 515232983265ull,
522148862235ull, 527279998245ull, 532411134255ull, 536281310565ull,
540219384705ull, 542896499145ull,

549812378115ull, 562441374495ull, 572567850855ull, 581937751395ull,
589969094715ull, 591976930545ull, 598446623775ull, 612055288845ull,
614761502355ull, 632802925755ull, 639718804725ull, 640077693255ull,
644849940735ull, 650204169615ull, 656227677105ull, 663143556075ull,
676306035405ull, 689022328995ull, 693241694145ull, 695045836485ull,
695492022225ull, 701961715455ull, 708431408685ull, 708877594425ull,
727171209765ull, 730517602815ull, 736541110305ull, 743456989275ull,
755280911385ull, 757288747215ull, 760188954525ull, 762642976095ull,
766658647755ull, 805476807135ull, 806601971175ull, 812615778975ull,
819531657945ull, 837825273285ull, 846079709475ull, 857234352975ull,
861026931765ull, 863704046205ull, 887351890425ull, 893802184275ull,
895606326615ull, 902522205585ull, 916353963525ull, 923269842495ull,
928400978505ull, 937402290825ull, 941340364965ull, 954279751425ull,
969896252325ull, 973688831115ull, 980158524345ull, 983281824525ull,
986628217575ull, 992428632195ull, 999567604035ull, 1006260390135ull,
1020092148075ull, 1024602503925ull, 1038385763415ull, 1040839784985ull,
1047755663955ull, 1054671542925ull, 1057794843105ull, 1064264536335ull,
1066796155425ull, 1068503300865ull, 1074749901225ull, 1077203922795ull,
1083673616025ull, 1093712795175ull, 1096613002485ull,

1103082695715ull, 1129242759645ull, 1137662090565ull, 1148370548325ull,
1151493848505ull, 1158409727475ull, 1180719014475ull, 1184511593265ull,
1193658400935ull, 1200128094165ull, 1242292646595ull, 1251323057985ull,
1251885640005ull, 1275979669965ull, 1277764412925ull, 1282895548935ull,
1303643185845ull, 1316582572305ull, 1335991651995ull, 1338222580695ull,
1357854753255ull, 1368340118145ull, 1376594554335ull, 1381279504605ull,
1407381370395ull, 1428129007305ull, 1451883548115ull, 1458915823365ull,
1462708402155ull, 1476540160095ull, 1504203675975ull, 1517143062435ull,
1522943477055ull, 1530082448895ull, 1543021835355ull, 1555961221815ull,
1559530707735ull, 1597233402765ull, 1607718767655ull, 1608843931695ull,
1621773618465ull, 1635605376405ull, 1649437134345ull, 1653006620265ull,
1663268892285ull, 1672415699955ull, 1704764166105ull, 1717703552565ull,
1724173245795ull, 1745813254185ull, 1769907284145ull, 1774524336585ull,
1775930791635ull, 1787754713745ull, 1788870178095ull, 1795339871325ull,
1808279257785ull, 1836165866535ull, 1843081745505ull, 1861375360845ull,
1885915576545ull, 1898408777265ull, 1905324656235ull, 1919156414175ull,
1932988172115ull, 1950612508845ull, 1952174158935ull, 1966364805405ull,
1976491281765ull, 1989430668225ull, 2015978719755ull, 2034718520835ull,
2067066986985ull, 2085137509455ull, 2086476066675ull, 2097165125055ull,
2099415453135ull, 2105885146365ull, 2112801025335ull, 2125294226055ull,
2126632783275ull, 2133771755115ull, 2163220013955ull, 2170582078665ull,
2175043936065ull,
};


/* All the prime factors that can appear in eulerphi(P) */
const unsigned long phiPfactors[] = {2UL, 3UL, 5UL, 7UL, 11UL, 13UL, 
				     17UL, 19UL};

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
  /* We need B2min >= 2 * max(S_1 + S_2) + (2*m_1 - 1)*P + 1, or
     B2min - 2 * max(S_1 + S_2) - 1 >= (2*m_1)*P - P, or
     (B2min - 2*max(S_1 + S_2) + P - 1)/(2P) >= m_1
     Choose m_1 accordingly */
  
  mpz_init (m);
  mpz_init (Ptmp);
  mpz_init (nrtmp);
  sets_max (m, P);
  mpz_mul_2exp (m, m, 1UL); /* m = 2*max(S_1 + S_2) */

  mpz_sub (m_1, B2min, m);
  mpz_sub_ui (m_1, m_1, 1UL); /* m_1 = B2min - 2*max(S_1 + S_2) - 1 */
  mpz_set_uint64(Ptmp, P);
  mpz_add (m_1, m_1, Ptmp);
  mpz_fdiv_q_2exp (m_1, m_1, 1UL);
  mpz_fdiv_q (m_1, m_1, Ptmp);    /* 2UL*P may overflow */
  
  /* Compute effB2min = 2 * max(S_1 + S_2) + (2*(m_1 - 1) + 1)*P + 1 */
  
  mpz_mul_2exp (effB2min, m_1, 1UL);
  mpz_sub_ui (effB2min, effB2min, 1UL);
  mpz_mul (effB2min, effB2min, Ptmp);
  mpz_add (effB2min, effB2min, m);
  mpz_add_ui (effB2min, effB2min, 1UL);
  ASSERT_ALWAYS (mpz_cmp (effB2min, B2min) <= 0);

  /* Compute the smallest value coprime to P at the high end of the stage 2
     interval that will not be covered: 
     2*(min(S_1 + S_2)) + (2*(m_1 + nr) + 1)*P. 
     We assume min(S_1 + S_2) = -max(S_1 + S_2) */
  mpz_set_uint64 (nrtmp, nr);
  mpz_add (effB2, m_1, nrtmp);
  mpz_mul_2exp (effB2, effB2, 1UL);
  mpz_add_ui (effB2, effB2, 1UL);
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
    
    ASSERT_ALWAYS (phiP > 0UL);

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

     Let F(x) = \prod_{k_1 \in S_1} (x - b_1^{2 k_1}).

     If we evaluate F(b_1^{2 k_2 + (2m + 1)P}) for all k_2 \in S_2 with 
     m_1 <= m < m_1+nr, we test all exponents 2 k_2 + (2m + 1)P - 2 k_1.
     The largest value coprime to P at the low end of the stage 2 interval 
     *not* covered will be 
       2*max(S_2) + (2*(m_1-1) + 1)*P - 2*min(S_1).
     The smallest value at the high end not covered will be
       2*min(S_2) + (2*(m_1 + nr) + 1)*P - 2*max(S_1).
     Assume S_1 and S_2 are symmetric around 0, so that max(S_1) = -min(S_1).
     Then the largest ... is:
       2*(max(S_1) + max(S_2)) + (2*m_1 - 1)*P
     The smallest ... is:
       -2*(max(S_1) + max(S_2)) + (2*m_1 + 2*nr + 1)*P
     The effective B2min = 2*(max(S_1) + max(S_2)) + (2*m_1 - 1)*P + 1
     The effective B2max = -2*(max(S_1) + max(S_2)) + (2*m_1 + 2*nr + 1)*P - 1

     Then the difference effB2max - effB2min =
       -4*(max(S_1) + max(S_2)) + 2P*(nr + 1) - 2

     We obviously require B2max - B2min <= 2*nr*P
     Since nr < lmax, B2max - B2min <= 2*lmax*P or
     P >= ceil((B2max - B2min)/(2*lmax))

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
      
      /* We have nr < l and effB2-effB2min <= 2*nr*P. Hence we need 
	 l >= B2l/P/2 */
      mpz_set_uint64(Ptmp, tryP);
      mpz_cdiv_q (lmin, B2l, Ptmp);
      mpz_cdiv_q_2exp (lmin, lmin, 1UL);
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
