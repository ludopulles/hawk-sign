/*
 * Hawk key pair generation.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2022 Hawk Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Ludo Pulles <ludo.pulles@cwi.nl>
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 */

#include "inner.h"
#include <assert.h>

/* ==================================================================== */
/*
 * Modular arithmetics.
 *
 * We implement a few functions for computing modulo a small integer p.
 *
 * All functions require that 2^30 < p < 2^31. Moreover, operands must
 * be in the 0..p-1 range.
 *
 * Modular addition and subtraction work for all such p.
 *
 * Montgomery multiplication requires that p is odd, and must be provided
 * with an additional value p0i = -1/p mod 2^31. See below for some basics
 * on Montgomery multiplication.
 *
 * Division computes an inverse modulo p by an exponentiation (with
 * exponent p-2): this works only if p is prime. Multiplication
 * requirements also apply, i.e. p must be odd and p0i must be provided.
 *
 * The NTT and inverse NTT need all of the above, and also that
 * p = 1 mod 2048.
 *
 * -----------------------------------------------------------------------
 *
 * We use Montgomery representation with 31-bit values:
 *
 *   Let R = 2^31 mod p. When 2^30 < p < 2^31, R = 2^31 - p.
 *   Montgomery representation of an integer x modulo p is x*R mod p.
 *
 *   Montgomery multiplication computes (x*y)/R mod p for
 *   operands x and y. Therefore:
 *
 *    - if operands are x*R and y*R (Montgomery representations of x and
 *      y), then Montgomery multiplication computes (x*R*y*R)/R = (x*y)*R
 *      mod p, which is the Montgomery representation of the product x*y;
 *
 *    - if operands are x*R and y (or x and y*R), then Montgomery
 *      multiplication returns x*y mod p: mixed-representation
 *      multiplications yield results in normal representation.
 *
 * To convert to Montgomery representation, we multiply by R, which is done
 * by Montgomery-multiplying by R^2. Stand-alone conversion back from
 * Montgomery representation is Montgomery-multiplication by 1.
 */

/*
 * Precomputed small primes. Each element contains the following:
 *
 *  p   The prime itself.
 *
 *  g   A primitive root of phi = X^N+1 (in field Z_p).
 *
 *  s   The inverse of the product of all previous primes in the array,
 *      computed modulo p and in Montgomery representation.
 *
 * All primes are such that p = 1 mod 2048, and are lower than 2^31. They
 * are listed in decreasing order.
 */

typedef struct {
	uint32_t p;
	uint32_t g;
	uint32_t s;
} small_prime;

static const small_prime PRIMES[] = {
	{ 2147473409,  383167813,      10239 },
	{ 2147389441,  211808905,  471403745 },
	{ 2147387393,   37672282, 1329335065 },
	{ 2147377153, 1977035326,  968223422 },
	{ 2147358721, 1067163706,  132460015 },
	{ 2147352577, 1606082042,  598693809 },
	{ 2147346433, 2033915641, 1056257184 },
	{ 2147338241, 1653770625,  421286710 },
	{ 2147309569,  631200819, 1111201074 },
	{ 2147297281, 2038364663, 1042003613 },
	{ 2147295233, 1962540515,   19440033 },
	{ 2147239937, 2100082663,  353296760 },
	{ 2147235841, 1991153006, 1703918027 },
	{ 2147217409,  516405114, 1258919613 },
	{ 2147205121,  409347988, 1089726929 },
	{ 2147196929,  927788991, 1946238668 },
	{ 2147178497, 1136922411, 1347028164 },
	{ 2147100673,  868626236,  701164723 },
	{ 2147082241, 1897279176,  617820870 },
	{ 2147074049, 1888819123,  158382189 },
	{ 2147051521,   25006327,  522758543 },
	{ 2147043329,  327546255,   37227845 },
	{ 2147039233,  766324424, 1133356428 },
	{ 2146988033, 1862817362,   73861329 },
	{ 2146963457,  404622040,  653019435 },
	{ 2146959361, 1936581214,  995143093 },
	{ 2146938881, 1559770096,  634921513 },
	{ 2146908161,  422623708, 1985060172 },
	{ 2146885633, 1751189170,  298238186 },
	{ 2146871297,  578919515,  291810829 },
	{ 2146846721, 1114060353,  915902322 },
	{ 2146834433, 2069565474,   47859524 },
	{ 2146818049, 1552824584,  646281055 },
	{ 2146775041, 1906267847, 1597832891 },
	{ 2146756609, 1847414714, 1228090888 },
	{ 2146744321, 1818792070, 1176377637 },
	{ 2146738177, 1118066398, 1054971214 },
	{ 2146736129,   52057278,  933422153 },
	{ 2146713601,  592259376, 1406621510 },
	{ 2146695169,  263161877, 1514178701 },
	{ 2146656257,  685363115,  384505091 },
	{ 2146650113,  927727032,  537575289 },
	{ 2146646017,   52575506, 1799464037 },
	{ 2146643969, 1276803876, 1348954416 },
	{ 2146603009,  814028633, 1521547704 },
	{ 2146572289, 1846678872, 1310832121 },
	{ 2146547713,  919368090, 1019041349 },
	{ 2146508801,  671847612,   38582496 },
	{ 2146492417,  283911680,  532424562 },
	{ 2146490369, 1780044827,  896447978 },
	{ 2146459649,  327980850, 1327906900 },
	{ 2146447361, 1310561493,  958645253 },
	{ 2146441217,  412148926,  287271128 },
	{ 2146437121,  293186449, 2009822534 },
	{ 2146430977,  179034356, 1359155584 },
	{ 2146418689, 1517345488, 1790248672 },
	{ 2146406401, 1615820390, 1584833571 },
	{ 2146404353,  826651445,  607120498 },
	{ 2146379777,    3816988, 1897049071 },
	{ 2146363393, 1221409784, 1986921567 },
	{ 2146355201, 1388081168,  849968120 },
	{ 2146336769, 1803473237, 1655544036 },
	{ 2146312193, 1023484977,  273671831 },
	{ 2146293761, 1074591448,  467406983 },
	{ 2146283521,  831604668, 1523950494 },
	{ 2146203649,  712865423, 1170834574 },
	{ 2146154497, 1764991362, 1064856763 },
	{ 2146142209,  627386213, 1406840151 },
	{ 2146127873, 1638674429, 2088393537 },
	{ 2146099201, 1516001018,  690673370 },
	{ 2146093057, 1294931393,  315136610 },
	{ 2146091009, 1942399533,  973539425 },
	{ 2146078721, 1843461814, 2132275436 },
	{ 2146060289, 1098740778,  360423481 },
	{ 2146048001, 1617213232, 1951981294 },
	{ 2146041857, 1805783169, 2075683489 },
	{ 2146019329,  272027909, 1753219918 },
	{ 2145986561, 1206530344, 2034028118 },
	{ 2145976321, 1243769360, 1173377644 },
	{ 2145964033,  887200839, 1281344586 },
	{ 2145906689, 1651026455,  906178216 },
	{ 2145875969, 1673238256, 1043521212 },
	{ 2145871873, 1226591210, 1399796492 },
	{ 2145841153, 1465353397, 1324527802 },
	{ 2145832961, 1150638905,  554084759 },
	{ 2145816577,  221601706,  427340863 },
	{ 2145785857,  608896761,  316590738 },
	{ 2145755137, 1712054942, 1684294304 },
	{ 2145742849, 1302302867,  724873116 },
	{ 2145728513,  516717693,  431671476 },
	{ 2145699841,  524575579, 1619722537 },
	{ 2145691649, 1925625239,  982974435 },
	{ 2145687553,  463795662, 1293154300 },
	{ 2145673217,  771716636,  881778029 },
	{ 2145630209, 1509556977,  837364988 },
	{ 2145595393,  229091856,  851648427 },
	{ 2145587201, 1796903241,  635342424 },
	{ 2145525761,  715310882, 1677228081 },
	{ 2145495041, 1040930522,  200685896 },
	{ 2145466369,  949804237, 1809146322 },
	{ 2145445889, 1673903706,   95316881 },
	{ 2145390593,  806941852, 1428671135 },
	{ 2145372161, 1402525292,  159350694 },
	{ 2145361921, 2124760298, 1589134749 },
	{ 2145359873, 1217503067, 1561543010 },
	{ 2145355777,  338341402,   83865711 },
	{ 2145343489, 1381532164,  641430002 },
	{ 2145325057, 1883895478, 1528469895 },
	{ 2145318913, 1335370424,   65809740 },
	{ 2145312769, 2000008042, 1919775760 },
	{ 2145300481,  961450962, 1229540578 },
	{ 2145282049,  910466767, 1964062701 },
	{ 2145232897,  816527501,  450152063 },
	{ 2145218561, 1435128058, 1794509700 },
	{ 2145187841,   33505311, 1272467582 },
	{ 2145181697,  269767433, 1380363849 },
	{ 2145175553,   56386299, 1316870546 },
	{ 2145079297, 2106880293, 1391797340 },
	{ 2145021953, 1347906152,  720510798 },
	{ 2145015809,  206769262, 1651459955 },
	{ 2145003521, 1885513236, 1393381284 },
	{ 2144960513, 1810381315,   31937275 },
	{ 2144944129, 1306487838, 2019419520 },
	{ 2144935937,   37304730, 1841489054 },
	{ 2144894977, 1601434616,  157985831 },
	{ 2144888833,   98749330, 2128592228 },
	{ 2144880641, 1772327002, 2076128344 },
	{ 2144864257, 1404514762, 2029969964 },
	{ 2144827393,  801236594,  406627220 },
	{ 2144806913,  349217443, 1501080290 },
	{ 2144796673, 1542656776, 2084736519 },
	{ 2144778241, 1210734884, 1746416203 },
	{ 2144759809, 1146598851,  716464489 },
	{ 2144757761,  286328400, 1823728177 },
	{ 2144729089, 1347555695, 1836644881 },
	{ 2144727041, 1795703790,  520296412 },
	{ 2144696321, 1302475157,  852964281 },
	{ 2144667649, 1075877614,  504992927 },
	{ 2144573441,  198765808, 1617144982 },
	{ 2144555009,  321528767,  155821259 },
	{ 2144550913,  814139516, 1819937644 },
	{ 2144536577,  571143206,  962942255 },
	{ 2144524289, 1746733766,    2471321 },
	{ 2144512001, 1821415077,  124190939 },
	{ 2144468993,  917871546, 1260072806 },
	{ 2144458753,  378417981, 1569240563 },
	{ 2144421889,  175229668, 1825620763 },
	{ 2144409601, 1699216963,  351648117 },
	{ 2144370689, 1071885991,  958186029 },
	{ 2144348161, 1763151227,  540353574 },
	{ 2144335873, 1060214804,  919598847 },
	{ 2144329729,  663515846, 1448552668 },
	{ 2144327681, 1057776305,  590222840 },
	{ 2144309249, 1705149168, 1459294624 },
	{ 2144296961,  325823721, 1649016934 },
	{ 2144290817,  738775789,  447427206 },
	{ 2144243713,  962347618,  893050215 },
	{ 2144237569, 1655257077,  900860862 },
	{ 2144161793,  242206694, 1567868672 },
	{ 2144155649,  769415308, 1247993134 },
	{ 2144137217,  320492023,  515841070 },
	{ 2144120833, 1639388522,  770877302 },
	{ 2144071681, 1761785233,  964296120 },
	{ 2144065537,  419817825,  204564472 },
	{ 2144028673,  666050597, 2091019760 },
	{ 2144010241, 1413657615, 1518702610 },
	{ 2143952897, 1238327946,  475672271 },
	{ 2143940609,  307063413, 1176750846 },
	{ 2143918081, 2062905559,  786785803 },
	{ 2143899649, 1338112849, 1562292083 },
	{ 2143891457,   68149545,   87166451 },
	{ 2143885313,  921750778,  394460854 },
	{ 2143854593,  719766593,  133877196 },
	{ 2143836161, 1149399850, 1861591875 },
	{ 2143762433, 1848739366, 1335934145 },
	{ 2143756289, 1326674710,  102999236 },
	{ 2143713281,  808061791, 1156900308 },
	{ 2143690753,  388399459, 1926468019 },
	{ 2143670273, 1427891374, 1756689401 },
	{ 2143666177, 1912173949,  986629565 },
	{ 2143645697, 2041160111,  371842865 },
	{ 2143641601, 1279906897, 2023974350 },
	{ 2143635457,  720473174, 1389027526 },
	{ 2143621121, 1298309455, 1732632006 },
	{ 2143598593, 1548762216, 1825417506 },
	{ 2143567873,  620475784, 1073787233 },
	{ 2143561729, 1932954575,  949167309 },
	{ 2143553537,  354315656, 1652037534 },
	{ 2143541249,  577424288, 1097027618 },
	{ 2143531009,  357862822,  478640055 },
	{ 2143522817, 2017706025, 1550531668 },
	{ 2143506433, 2078127419, 1824320165 },
	{ 2143488001,  613475285, 1604011510 },
	{ 2143469569, 1466594987,  502095196 },
	{ 2143426561, 1115430331, 1044637111 },
	{ 2143383553,    9778045, 1902463734 },
	{ 2143377409, 1557401276, 2056861771 },
	{ 2143363073,  652036455, 1965915971 },
	{ 2143260673, 1464581171, 1523257541 },
	{ 2143246337, 1876119649,  764541916 },
	{ 2143209473, 1614992673, 1920672844 },
	{ 2143203329,  981052047, 2049774209 },
	{ 2143160321, 1847355533,  728535665 },
	{ 2143129601,  965558457,  603052992 },
	{ 2143123457, 2140817191,    8348679 },
	{ 2143100929, 1547263683,  694209023 },
	{ 2143092737,  643459066, 1979934533 },
	{ 2143082497,  188603778, 2026175670 },
	{ 2143062017, 1657329695,  377451099 },
	{ 2143051777,  114967950,  979255473 },
	{ 2143025153, 1698431342, 1449196896 },
	{ 2143006721, 1862741675, 1739650365 },
	{ 2142996481,  756660457,  996160050 },
	{ 2142976001,  927864010, 1166847574 },
	{ 2142965761,  905070557,  661974566 },
	{ 2142916609,   40932754, 1787161127 },
	{ 2142892033, 1987985648,  675335382 },
	{ 2142885889,  797497211, 1323096997 },
	{ 2142871553, 2068025830, 1411877159 },
	{ 2142861313, 1217177090, 1438410687 },
	{ 2142830593,  409906375, 1767860634 },
	{ 2142803969, 1197788993,  359782919 },
	{ 2142785537,  643817365,  513932862 },
	{ 2142779393, 1717046338,  218943121 },
	{ 2142724097,   89336830,  416687049 },
	{ 2142707713,    5944581, 1356813523 },
	{ 2142658561,  887942135, 2074011722 },
	{ 2142638081,  151851972, 1647339939 },
	{ 2142564353, 1691505537, 1483107336 },
	{ 2142533633, 1989920200, 1135938817 },
	{ 2142529537,  959263126, 1531961857 },
	{ 2142527489,  453251129, 1725566162 },
	{ 2142502913, 1536028102,  182053257 },
	{ 2142498817,  570138730,  701443447 },
	{ 2142416897,  326965800,  411931819 },
	{ 2142363649, 1675665410, 1517191733 },
	{ 2142351361,  968529566, 1575712703 },
	{ 2142330881, 1384953238, 1769087884 },
	{ 2142314497, 1977173242, 1833745524 },
	{ 2142289921,   95082313, 1714775493 },
	{ 2142283777,  109377615, 1070584533 },
	{ 2142277633,   16960510,  702157145 },
	{ 2142263297,  553850819,  431364395 },
	{ 2142208001,  241466367, 2053967982 },
	{ 2142164993, 1795661326, 1031836848 },
	{ 2142097409, 1212530046,  712772031 },
	{ 2142087169, 1763869720,  822276067 },
	{ 2142078977,  644065713, 1765268066 },
	{ 2142074881,  112671944,  643204925 },
	{ 2142044161, 1387785471, 1297890174 },
	{ 2142025729,  783885537, 1000425730 },
	{ 2142011393,  905662232, 1679401033 },
	{ 2141974529,  799788433,  468119557 },
	{ 2141943809, 1932544124,  449305555 },
	{ 2141933569, 1527403256,  841867925 },
	{ 2141931521, 1247076451,  743823916 },
	{ 2141902849, 1199660531,  401687910 },
	{ 2141890561,  150132350, 1720336972 },
	{ 2141857793, 1287438162,  663880489 },
	{ 2141833217,  618017731, 1819208266 },
	{ 2141820929,  999578638, 1403090096 },
	{ 2141786113,   81834325, 1523542501 },
	{ 2141771777,  120001928,  463556492 },
	{ 2141759489,  122455485, 2124928282 },
	{ 2141749249,  141986041,  940339153 },
	{ 2141685761,  889088734,  477141499 },
	{ 2141673473,  324212681, 1122558298 },
	{ 2141669377, 1175806187, 1373818177 },
	{ 2141655041, 1113654822,  296887082 },
	{ 2141587457,  991103258, 1585913875 },
	{ 2141583361, 1401451409, 1802457360 },
	{ 2141575169, 1571977166,  712760980 },
	{ 2141546497, 1107849376, 1250270109 },
	{ 2141515777,  196544219,  356001130 },
	{ 2141495297, 1733571506, 1060744866 },
	{ 2141483009,  321552363, 1168297026 },
	{ 2141458433,  505818251,  733225819 },
	{ 2141360129, 1026840098,  948342276 },
	{ 2141325313,  945133744, 2129965998 },
	{ 2141317121, 1871100260, 1843844634 },
	{ 2141286401, 1790639498, 1750465696 },
	{ 2141267969, 1376858592,  186160720 },
	{ 2141255681, 2129698296, 1876677959 },
	{ 2141243393, 2138900688, 1340009628 },
	{ 2141214721, 1933049835, 1087819477 },
	{ 2141212673, 1898664939, 1786328049 }
};

/*
 * Reduce a small signed integer modulo a small prime. The source
 * value x MUST be such that -p < x < p.
 */
static inline uint32_t
modp_set(int32_t x, uint32_t p)
{
	uint32_t w;

	w = (uint32_t)x;
	w += p & -(w >> 31);
	return w;
}

/*
 * Normalize a modular integer around 0.
 */
static inline int32_t
modp_norm(uint32_t x, uint32_t p)
{
	return (int32_t)(x - (p & (((x - ((p + 1) >> 1)) >> 31) - 1)));
}

/*
 * Compute -1/p mod 2^31. This works for all odd integers p that fit
 * on 31 bits.
 */
static uint32_t
modp_ninv31(uint32_t p)
{
	uint32_t y;

	y = 2 - p;
	y *= 2 - p * y;
	y *= 2 - p * y;
	y *= 2 - p * y;
	y *= 2 - p * y;
	return (uint32_t)0x7FFFFFFF & -y;
}

/*
 * Compute R = 2^31 mod p.
 */
static inline uint32_t
modp_R(uint32_t p)
{
	/*
	 * Since 2^30 < p < 2^31, we know that 2^31 mod p is simply
	 * 2^31 - p.
	 */
	return ((uint32_t)1 << 31) - p;
}

/*
 * Addition modulo p.
 */
static inline uint32_t
modp_add(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d;

	d = a + b - p;
	d += p & -(d >> 31);
	return d;
}

/*
 * Subtraction modulo p.
 */
static inline uint32_t
modp_sub(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d;

	d = a - b;
	d += p & -(d >> 31);
	return d;
}

/*
 * Montgomery multiplication modulo p. The 'p0i' value is -1/p mod 2^31.
 * It is required that p is an odd integer.
 */
static inline uint32_t
modp_montymul(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i)
{
	uint64_t z, w;
	uint32_t d;

	z = (uint64_t)a * (uint64_t)b;
	w = ((z * p0i) & (uint64_t)0x7FFFFFFF) * p;
	d = (uint32_t)((z + w) >> 31) - p;
	d += p & -(d >> 31);
	return d;
}

/*
 * Compute R2 = 2^62 mod p.
 */
static uint32_t
modp_R2(uint32_t p, uint32_t p0i)
{
	uint32_t z;

	/*
	 * Compute z = 2^31 mod p (this is the value 1 in Montgomery
	 * representation), then double it with an addition.
	 */
	z = modp_R(p);
	z = modp_add(z, z, p);

	/*
	 * Square it five times to obtain 2^32 in Montgomery representation
	 * (i.e. 2^63 mod p).
	 */
	z = modp_montymul(z, z, p, p0i);
	z = modp_montymul(z, z, p, p0i);
	z = modp_montymul(z, z, p, p0i);
	z = modp_montymul(z, z, p, p0i);
	z = modp_montymul(z, z, p, p0i);

	/*
	 * Halve the value mod p to get 2^62.
	 */
	z = (z + (p & -(z & 1))) >> 1;
	return z;
}

/*
 * Compute 2^(31*x) modulo p. This works for integers x up to 2^11.
 * p must be prime such that 2^30 < p < 2^31; p0i must be equal to
 * -1/p mod 2^31; R2 must be equal to 2^62 mod p.
 */
static inline uint32_t
modp_Rx(unsigned x, uint32_t p, uint32_t p0i, uint32_t R2)
{
	int i;
	uint32_t r, z;

	/*
	 * 2^(31*x) = (2^31)*(2^(31*(x-1))); i.e. we want the Montgomery
	 * representation of (2^31)^e mod p, where e = x-1.
	 * R2 is 2^31 in Montgomery representation.
	 */
	x --;
	r = R2;
	z = modp_R(p);
	for (i = 0; (1U << i) <= x; i ++) {
		if ((x & (1U << i)) != 0) {
			z = modp_montymul(z, r, p, p0i);
		}
		r = modp_montymul(r, r, p, p0i);
	}
	return z;
}

/*
 * Division modulo p. If the divisor (b) is 0, then 0 is returned.
 * This function computes proper results only when p is prime.
 * Parameters:
 *   a     dividend
 *   b     divisor
 *   p     odd prime modulus
 *   p0i   -1/p mod 2^31
 *   R     2^31 mod R
 */
static uint32_t
modp_div(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i)
{
	uint32_t z, e;
	int i;

	e = p - 2;
	z = modp_R(p);
	for (i = 30; i >= 0; i --) {
		uint32_t z2;

		z = modp_montymul(z, z, p, p0i);
		z2 = modp_montymul(z, b, p, p0i);
		z ^= (z ^ z2) & -(uint32_t)((e >> i) & 1);
	}

	/*
	 * The loop above just assumed that b was in Montgomery
	 * representation, i.e. really contained b*R; under that
	 * assumption, it returns 1/b in Montgomery representation,
	 * which is R/b. But we gave it b in normal representation,
	 * so the loop really returned R/(b/R) = R^2/b.
	 *
	 * We want a/b, so we need one Montgomery multiplication with a,
	 * which also remove one of the R factors, and another such
	 * multiplication to remove the second R factor.
	 */
	z = modp_montymul(z, 1, p, p0i);
	return modp_montymul(a, z, p, p0i);
}

/*
 * Bit-reversal index table.
 */
static const uint16_t REV10[] = {
	   0,  512,  256,  768,  128,  640,  384,  896,   64,  576,  320,  832,
	 192,  704,  448,  960,   32,  544,  288,  800,  160,  672,  416,  928,
	  96,  608,  352,  864,  224,  736,  480,  992,   16,  528,  272,  784,
	 144,  656,  400,  912,   80,  592,  336,  848,  208,  720,  464,  976,
	  48,  560,  304,  816,  176,  688,  432,  944,  112,  624,  368,  880,
	 240,  752,  496, 1008,    8,  520,  264,  776,  136,  648,  392,  904,
	  72,  584,  328,  840,  200,  712,  456,  968,   40,  552,  296,  808,
	 168,  680,  424,  936,  104,  616,  360,  872,  232,  744,  488, 1000,
	  24,  536,  280,  792,  152,  664,  408,  920,   88,  600,  344,  856,
	 216,  728,  472,  984,   56,  568,  312,  824,  184,  696,  440,  952,
	 120,  632,  376,  888,  248,  760,  504, 1016,    4,  516,  260,  772,
	 132,  644,  388,  900,   68,  580,  324,  836,  196,  708,  452,  964,
	  36,  548,  292,  804,  164,  676,  420,  932,  100,  612,  356,  868,
	 228,  740,  484,  996,   20,  532,  276,  788,  148,  660,  404,  916,
	  84,  596,  340,  852,  212,  724,  468,  980,   52,  564,  308,  820,
	 180,  692,  436,  948,  116,  628,  372,  884,  244,  756,  500, 1012,
	  12,  524,  268,  780,  140,  652,  396,  908,   76,  588,  332,  844,
	 204,  716,  460,  972,   44,  556,  300,  812,  172,  684,  428,  940,
	 108,  620,  364,  876,  236,  748,  492, 1004,   28,  540,  284,  796,
	 156,  668,  412,  924,   92,  604,  348,  860,  220,  732,  476,  988,
	  60,  572,  316,  828,  188,  700,  444,  956,  124,  636,  380,  892,
	 252,  764,  508, 1020,    2,  514,  258,  770,  130,  642,  386,  898,
	  66,  578,  322,  834,  194,  706,  450,  962,   34,  546,  290,  802,
	 162,  674,  418,  930,   98,  610,  354,  866,  226,  738,  482,  994,
	  18,  530,  274,  786,  146,  658,  402,  914,   82,  594,  338,  850,
	 210,  722,  466,  978,   50,  562,  306,  818,  178,  690,  434,  946,
	 114,  626,  370,  882,  242,  754,  498, 1010,   10,  522,  266,  778,
	 138,  650,  394,  906,   74,  586,  330,  842,  202,  714,  458,  970,
	  42,  554,  298,  810,  170,  682,  426,  938,  106,  618,  362,  874,
	 234,  746,  490, 1002,   26,  538,  282,  794,  154,  666,  410,  922,
	  90,  602,  346,  858,  218,  730,  474,  986,   58,  570,  314,  826,
	 186,  698,  442,  954,  122,  634,  378,  890,  250,  762,  506, 1018,
	   6,  518,  262,  774,  134,  646,  390,  902,   70,  582,  326,  838,
	 198,  710,  454,  966,   38,  550,  294,  806,  166,  678,  422,  934,
	 102,  614,  358,  870,  230,  742,  486,  998,   22,  534,  278,  790,
	 150,  662,  406,  918,   86,  598,  342,  854,  214,  726,  470,  982,
	  54,  566,  310,  822,  182,  694,  438,  950,  118,  630,  374,  886,
	 246,  758,  502, 1014,   14,  526,  270,  782,  142,  654,  398,  910,
	  78,  590,  334,  846,  206,  718,  462,  974,   46,  558,  302,  814,
	 174,  686,  430,  942,  110,  622,  366,  878,  238,  750,  494, 1006,
	  30,  542,  286,  798,  158,  670,  414,  926,   94,  606,  350,  862,
	 222,  734,  478,  990,   62,  574,  318,  830,  190,  702,  446,  958,
	 126,  638,  382,  894,  254,  766,  510, 1022,    1,  513,  257,  769,
	 129,  641,  385,  897,   65,  577,  321,  833,  193,  705,  449,  961,
	  33,  545,  289,  801,  161,  673,  417,  929,   97,  609,  353,  865,
	 225,  737,  481,  993,   17,  529,  273,  785,  145,  657,  401,  913,
	  81,  593,  337,  849,  209,  721,  465,  977,   49,  561,  305,  817,
	 177,  689,  433,  945,  113,  625,  369,  881,  241,  753,  497, 1009,
	   9,  521,  265,  777,  137,  649,  393,  905,   73,  585,  329,  841,
	 201,  713,  457,  969,   41,  553,  297,  809,  169,  681,  425,  937,
	 105,  617,  361,  873,  233,  745,  489, 1001,   25,  537,  281,  793,
	 153,  665,  409,  921,   89,  601,  345,  857,  217,  729,  473,  985,
	  57,  569,  313,  825,  185,  697,  441,  953,  121,  633,  377,  889,
	 249,  761,  505, 1017,    5,  517,  261,  773,  133,  645,  389,  901,
	  69,  581,  325,  837,  197,  709,  453,  965,   37,  549,  293,  805,
	 165,  677,  421,  933,  101,  613,  357,  869,  229,  741,  485,  997,
	  21,  533,  277,  789,  149,  661,  405,  917,   85,  597,  341,  853,
	 213,  725,  469,  981,   53,  565,  309,  821,  181,  693,  437,  949,
	 117,  629,  373,  885,  245,  757,  501, 1013,   13,  525,  269,  781,
	 141,  653,  397,  909,   77,  589,  333,  845,  205,  717,  461,  973,
	  45,  557,  301,  813,  173,  685,  429,  941,  109,  621,  365,  877,
	 237,  749,  493, 1005,   29,  541,  285,  797,  157,  669,  413,  925,
	  93,  605,  349,  861,  221,  733,  477,  989,   61,  573,  317,  829,
	 189,  701,  445,  957,  125,  637,  381,  893,  253,  765,  509, 1021,
	   3,  515,  259,  771,  131,  643,  387,  899,   67,  579,  323,  835,
	 195,  707,  451,  963,   35,  547,  291,  803,  163,  675,  419,  931,
	  99,  611,  355,  867,  227,  739,  483,  995,   19,  531,  275,  787,
	 147,  659,  403,  915,   83,  595,  339,  851,  211,  723,  467,  979,
	  51,  563,  307,  819,  179,  691,  435,  947,  115,  627,  371,  883,
	 243,  755,  499, 1011,   11,  523,  267,  779,  139,  651,  395,  907,
	  75,  587,  331,  843,  203,  715,  459,  971,   43,  555,  299,  811,
	 171,  683,  427,  939,  107,  619,  363,  875,  235,  747,  491, 1003,
	  27,  539,  283,  795,  155,  667,  411,  923,   91,  603,  347,  859,
	 219,  731,  475,  987,   59,  571,  315,  827,  187,  699,  443,  955,
	 123,  635,  379,  891,  251,  763,  507, 1019,    7,  519,  263,  775,
	 135,  647,  391,  903,   71,  583,  327,  839,  199,  711,  455,  967,
	  39,  551,  295,  807,  167,  679,  423,  935,  103,  615,  359,  871,
	 231,  743,  487,  999,   23,  535,  279,  791,  151,  663,  407,  919,
	  87,  599,  343,  855,  215,  727,  471,  983,   55,  567,  311,  823,
	 183,  695,  439,  951,  119,  631,  375,  887,  247,  759,  503, 1015,
	  15,  527,  271,  783,  143,  655,  399,  911,   79,  591,  335,  847,
	 207,  719,  463,  975,   47,  559,  303,  815,  175,  687,  431,  943,
	 111,  623,  367,  879,  239,  751,  495, 1007,   31,  543,  287,  799,
	 159,  671,  415,  927,   95,  607,  351,  863,  223,  735,  479,  991,
	  63,  575,  319,  831,  191,  703,  447,  959,  127,  639,  383,  895,
	 255,  767,  511, 1023
};

/*
 * Compute the roots for NTT and inverse NTT (binary case). Input
 * parameter g is a primitive 2048-th root of 1 modulo p (i.e. g^1024 =
 * -1 mod p). This fills gm[] and igm[] with powers of g and 1/g:
 *   gm[rev(i)] = g^i mod p
 *   igm[rev(i)] = (1/g)^i mod p
 * where rev() is the "bit reversal" function over 10 bits. It fills
 * the arrays only up to N = 2^logn values.
 *
 * The values stored in gm[] and igm[] are in Montgomery representation.
 *
 * p must be a prime such that p = 1 mod 2048.
 */
static void
modp_mkgm2(uint32_t *restrict gm, uint32_t *restrict igm, unsigned logn,
	uint32_t g, uint32_t p, uint32_t p0i)
{
	size_t u, n;
	unsigned k;
	uint32_t ig, x1, x2, R2;

	n = MKN(logn);

	/*
	 * We want g such that g^(2N) = 1 mod p, but the provided
	 * generator has order 2048. We must square it a few times.
	 */
	R2 = modp_R2(p, p0i);
	g = modp_montymul(g, R2, p, p0i);
	for (k = logn; k < 10; k ++) {
		g = modp_montymul(g, g, p, p0i);
	}

	ig = modp_div(R2, g, p, p0i);
	k = 10 - logn;
	x1 = x2 = modp_R(p);
	for (u = 0; u < n; u ++) {
		size_t v;

		v = REV10[u << k];
		gm[v] = x1;
		igm[v] = x2;
		x1 = modp_montymul(x1, g, p, p0i);
		x2 = modp_montymul(x2, ig, p, p0i);
	}
}

/*
 * Compute the NTT over a polynomial (binary case). Polynomial elements
 * are a[0], a[stride], a[2 * stride]...
 */
static void
modp_NTT2_ext(uint32_t *a, size_t stride, const uint32_t *gm, unsigned logn,
	uint32_t p, uint32_t p0i)
{
	size_t t, m, n;

	if (logn == 0) {
		return;
	}
	n = MKN(logn);
	t = n;
	for (m = 1; m < n; m <<= 1) {
		size_t ht, u, v1;

		ht = t >> 1;
		for (u = 0, v1 = 0; u < m; u ++, v1 += t) {
			uint32_t s;
			size_t v;
			uint32_t *r1, *r2;

			s = gm[m + u];
			r1 = a + v1 * stride;
			r2 = r1 + ht * stride;
			for (v = 0; v < ht; v ++, r1 += stride, r2 += stride) {
				uint32_t x, y;

				x = *r1;
				y = modp_montymul(*r2, s, p, p0i);
				*r1 = modp_add(x, y, p);
				*r2 = modp_sub(x, y, p);
			}
		}
		t = ht;
	}
}

/*
 * Compute the inverse NTT over a polynomial (binary case).
 */
static void
modp_iNTT2_ext(uint32_t *a, size_t stride, const uint32_t *igm, unsigned logn,
	uint32_t p, uint32_t p0i)
{
	size_t t, m, n, k;
	uint32_t ni;
	uint32_t *r;

	if (logn == 0) {
		return;
	}
	n = MKN(logn);
	t = 1;
	for (m = n; m > 1; m >>= 1) {
		size_t hm, dt, u, v1;

		hm = m >> 1;
		dt = t << 1;
		for (u = 0, v1 = 0; u < hm; u ++, v1 += dt) {
			uint32_t s;
			size_t v;
			uint32_t *r1, *r2;

			s = igm[hm + u];
			r1 = a + v1 * stride;
			r2 = r1 + t * stride;
			for (v = 0; v < t; v ++, r1 += stride, r2 += stride) {
				uint32_t x, y;

				x = *r1;
				y = *r2;
				*r1 = modp_add(x, y, p);
				*r2 = modp_montymul(
					modp_sub(x, y, p), s, p, p0i);
			}
		}
		t = dt;
	}

	/*
	 * We need 1/n in Montgomery representation, i.e. R/n. Since
	 * 1 <= logn <= 10, R/n is an integer; morever, R/n <= 2^30 < p,
	 * thus a simple shift will do.
	 */
	ni = (uint32_t)1 << (31 - logn);
	for (k = 0, r = a; k < n; k ++, r += stride) {
		*r = modp_montymul(*r, ni, p, p0i);
	}
}

/*
 * Simplified macros for NTT and iNTT (binary case) when the elements
 * are consecutive in RAM.
 */
#define modp_NTT2(a, gm, logn, p, p0i)   modp_NTT2_ext(a, 1, gm, logn, p, p0i)
#define modp_iNTT2(a, igm, logn, p, p0i) modp_iNTT2_ext(a, 1, igm, logn, p, p0i)

/*
 * Given polynomial f in NTT representation modulo p, compute f' of degree
 * less than N/2 such that f' = f0^2 - X*f1^2, where f0 and f1 are
 * polynomials of degree less than N/2 such that f = f0(X^2) + X*f1(X^2).
 *
 * The new polynomial is written "in place" over the first N/2 elements
 * of f.
 *
 * If applied logn times successively on a given polynomial, the resulting
 * degree-0 polynomial is the resultant of f and X^N+1 modulo p.
 *
 * This function applies only to the binary case; it is invoked from
 * solve_NTRU_binary_depth1().
 */
static void
modp_poly_rec_res(uint32_t *f, unsigned logn,
	uint32_t p, uint32_t p0i, uint32_t R2)
{
	size_t hn, u;

	hn = MKN(logn - 1);
	for (u = 0; u < hn; u ++) {
		uint32_t w0, w1;

		w0 = f[(u << 1) + 0];
		w1 = f[(u << 1) + 1];
		f[u] = modp_montymul(modp_montymul(w0, w1, p, p0i), R2, p, p0i);
	}
}

/* ==================================================================== */
/*
 * Custom bignum implementation.
 *
 * This is a very reduced set of functionalities. We need to do the
 * following operations:
 *
 *  - Rebuild the resultant and the polynomial coefficients from their
 *    values modulo small primes (of length 31 bits each).
 *
 *  - Compute an extended GCD between the two computed resultants.
 *
 *  - Extract top bits and add scaled values during the successive steps
 *    of Babai rounding.
 *
 * When rebuilding values using CRT, we must also recompute the product
 * of the small prime factors. We always do it one small factor at a
 * time, so the "complicated" operations can be done modulo the small
 * prime with the modp_* functions. CRT coefficients (inverses) are
 * precomputed.
 *
 * All values are positive until the last step: when the polynomial
 * coefficients have been rebuilt, we normalize them around 0. But then,
 * only additions and subtractions on the upper few bits are needed
 * afterwards.
 *
 * We keep big integers as arrays of 31-bit words (in uint32_t values);
 * the top bit of each uint32_t is kept equal to 0. Using 31-bit words
 * makes it easier to keep track of carries. When negative values are
 * used, two's complement is used.
 */

/*
 * Subtract integer b from integer a. Both integers are supposed to have
 * the same size. The carry (0 or 1) is returned. Source arrays a and b
 * MUST be distinct.
 *
 * The operation is performed as described above if ctr = 1. If
 * ctl = 0, the value a[] is unmodified, but all memory accesses are
 * still performed, and the carry is computed and returned.
 */
static uint32_t
zint_sub(uint32_t *restrict a, const uint32_t *restrict b, size_t len,
	uint32_t ctl)
{
	size_t u;
	uint32_t cc, m;

	cc = 0;
	m = -ctl;
	for (u = 0; u < len; u ++) {
		uint32_t aw, w;

		aw = a[u];
		w = aw - b[u] - cc;
		cc = w >> 31;
		aw ^= ((w & 0x7FFFFFFF) ^ aw) & m;
		a[u] = aw;
	}
	return cc;
}

/*
 * Mutiply the provided big integer m with a small value x.
 * This function assumes that x < 2^31. The carry word is returned.
 */
static uint32_t
zint_mul_small(uint32_t *m, size_t mlen, uint32_t x)
{
	size_t u;
	uint32_t cc;

	cc = 0;
	for (u = 0; u < mlen; u ++) {
		uint64_t z;

		z = (uint64_t)m[u] * (uint64_t)x + cc;
		m[u] = (uint32_t)z & 0x7FFFFFFF;
		cc = (uint32_t)(z >> 31);
	}
	return cc;
}

/*
 * Reduce a big integer d modulo a small integer p.
 * Rules:
 *  d is unsigned
 *  p is prime
 *  2^30 < p < 2^31
 *  p0i = -(1/p) mod 2^31
 *  R2 = 2^62 mod p
 */
static uint32_t
zint_mod_small_unsigned(const uint32_t *d, size_t dlen,
	uint32_t p, uint32_t p0i, uint32_t R2)
{
	uint32_t x;
	size_t u;

	/*
	 * Algorithm: we inject words one by one, starting with the high
	 * word. Each step is:
	 *  - multiply x by 2^31
	 *  - add new word
	 */
	x = 0;
	u = dlen;
	while (u -- > 0) {
		uint32_t w;

		x = modp_montymul(x, R2, p, p0i);
		w = d[u] - p;
		w += p & -(w >> 31);
		x = modp_add(x, w, p);
	}
	return x;
}

/*
 * Similar to zint_mod_small_unsigned(), except that d may be signed.
 * Extra parameter is Rx = 2^(31*dlen) mod p.
 */
static uint32_t
zint_mod_small_signed(const uint32_t *d, size_t dlen,
	uint32_t p, uint32_t p0i, uint32_t R2, uint32_t Rx)
{
	uint32_t z;

	if (dlen == 0) {
		return 0;
	}
	z = zint_mod_small_unsigned(d, dlen, p, p0i, R2);
	z = modp_sub(z, Rx & -(d[dlen - 1] >> 30), p);
	return z;
}

/*
 * Add y*s to x. x and y initially have length 'len' words; the new x
 * has length 'len+1' words. 's' must fit on 31 bits. x[] and y[] must
 * not overlap.
 */
static void
zint_add_mul_small(uint32_t *restrict x,
	const uint32_t *restrict y, size_t len, uint32_t s)
{
	size_t u;
	uint32_t cc;

	cc = 0;
	for (u = 0; u < len; u ++) {
		uint32_t xw, yw;
		uint64_t z;

		xw = x[u];
		yw = y[u];
		z = (uint64_t)yw * (uint64_t)s + (uint64_t)xw + (uint64_t)cc;
		x[u] = (uint32_t)z & 0x7FFFFFFF;
		cc = (uint32_t)(z >> 31);
	}
	x[len] = cc;
}

/*
 * Normalize a modular integer around 0: if x > p/2, then x is replaced
 * with x - p (signed encoding with two's complement); otherwise, x is
 * untouched. The two integers x and p are encoded over the same length.
 */
static void
zint_norm_zero(uint32_t *restrict x, const uint32_t *restrict p, size_t len)
{
	size_t u;
	uint32_t r, bb;

	/*
	 * Compare x with p/2. We use the shifted version of p, and p
	 * is odd, so we really compare with (p-1)/2; we want to perform
	 * the subtraction if and only if x > (p-1)/2.
	 */
	r = 0;
	bb = 0;
	u = len;
	while (u -- > 0) {
		uint32_t wx, wp, cc;

		/*
		 * Get the two words to compare in wx and wp (both over
		 * 31 bits exactly).
		 */
		wx = x[u];
		wp = (p[u] >> 1) | (bb << 30);
		bb = p[u] & 1;

		/*
		 * We set cc to -1, 0 or 1, depending on whether wp is
		 * lower than, equal to, or greater than wx.
		 */
		cc = wp - wx;
		cc = ((-cc) >> 31) | -(cc >> 31);

		/*
		 * If r != 0 then it is either 1 or -1, and we keep its
		 * value. Otherwise, if r = 0, then we replace it with cc.
		 */
		r |= cc & ((r & 1) - 1);
	}

	/*
	 * At this point, r = -1, 0 or 1, depending on whether (p-1)/2
	 * is lower than, equal to, or greater than x. We thus want to
	 * do the subtraction only if r = -1.
	 */
	zint_sub(x, p, len, r >> 31);
}

/*
 * Rebuild integers from their RNS representation. There are 'num'
 * integers, and each consists in 'xlen' words. 'xx' points at that
 * first word of the first integer; subsequent integers are accessed
 * by adding 'xstride' repeatedly.
 *
 * The words of an integer are the RNS representation of that integer,
 * using the provided 'primes' are moduli. This function replaces
 * each integer with its multi-word value (little-endian order).
 *
 * If "normalize_signed" is non-zero, then the returned value is
 * normalized to the -m/2..m/2 interval (where m is the product of all
 * small prime moduli); two's complement is used for negative values.
 */
static void
zint_rebuild_CRT(uint32_t *restrict xx, size_t xlen, size_t xstride,
	size_t num, int normalize_signed, uint32_t *restrict tmp)
{
	size_t u;
	uint32_t *x;

	tmp[0] = PRIMES[0].p;
	for (u = 1; u < xlen; u ++) {
		/*
		 * At the entry of each loop iteration:
		 *  - the first u words of each array have been
		 *    reassembled;
		 *  - the first u words of tmp[] contains the
		 * product of the prime moduli processed so far.
		 *
		 * We call 'q' the product of all previous primes.
		 */
		uint32_t p, p0i, s, R2;
		size_t v;

		p = PRIMES[u].p;
		s = PRIMES[u].s;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);

		for (v = 0, x = xx; v < num; v ++, x += xstride) {
			uint32_t xp, xq, xr;
			/*
			 * xp = the integer x modulo the prime p for this
			 *      iteration
			 * xq = (x mod q) mod p
			 */
			xp = x[u];
			xq = zint_mod_small_unsigned(x, u, p, p0i, R2);

			/*
			 * New value is (x mod q) + q * (s * (xp - xq) mod p)
			 */
			xr = modp_montymul(s, modp_sub(xp, xq, p), p, p0i);
			zint_add_mul_small(x, tmp, u, xr);
		}

		/*
		 * Update product of primes in tmp[].
		 */
		tmp[u] = zint_mul_small(tmp, u, p);
	}

	/*
	 * Normalize the reconstructed values around 0.
	 */
	if (normalize_signed) {
		for (u = 0, x = xx; u < num; u ++, x += xstride) {
			zint_norm_zero(x, tmp, xlen);
		}
	}
}

/*
 * Negate a big integer conditionally: value a is replaced with -a if
 * and only if ctl = 1. Control value ctl must be 0 or 1.
 */
static void
zint_negate(uint32_t *a, size_t len, uint32_t ctl)
{
	size_t u;
	uint32_t cc, m;

	/*
	 * If ctl = 1 then we flip the bits of a by XORing with
	 * 0x7FFFFFFF, and we add 1 to the value. If ctl = 0 then we XOR
	 * with 0 and add 0, which leaves the value unchanged.
	 */
	cc = ctl;
	m = -ctl >> 1;
	for (u = 0; u < len; u ++) {
		uint32_t aw;

		aw = a[u];
		aw = (aw ^ m) + cc;
		a[u] = aw & 0x7FFFFFFF;
		cc = aw >> 31;
	}
}

/*
 * Replace a with (a*xa+b*xb)/(2^31) and b with (a*ya+b*yb)/(2^31).
 * The low bits are dropped (the caller should compute the coefficients
 * such that these dropped bits are all zeros). If either or both
 * yields a negative value, then the value is negated.
 *
 * Returned value is:
 *  0  both values were positive
 *  1  new a had to be negated
 *  2  new b had to be negated
 *  3  both new a and new b had to be negated
 *
 * Coefficients xa, xb, ya and yb may use the full signed 32-bit range.
 */
static uint32_t
zint_co_reduce(uint32_t *a, uint32_t *b, size_t len,
	int64_t xa, int64_t xb, int64_t ya, int64_t yb)
{
	size_t u;
	int64_t cca, ccb;
	uint32_t nega, negb;

	cca = 0;
	ccb = 0;
	for (u = 0; u < len; u ++) {
		uint32_t wa, wb;
		uint64_t za, zb;

		wa = a[u];
		wb = b[u];
		za = wa * (uint64_t)xa + wb * (uint64_t)xb + (uint64_t)cca;
		zb = wa * (uint64_t)ya + wb * (uint64_t)yb + (uint64_t)ccb;
		if (u > 0) {
			a[u - 1] = (uint32_t)za & 0x7FFFFFFF;
			b[u - 1] = (uint32_t)zb & 0x7FFFFFFF;
		}
		cca = *(int64_t *)&za >> 31;
		ccb = *(int64_t *)&zb >> 31;
	}
	a[len - 1] = (uint32_t)cca;
	b[len - 1] = (uint32_t)ccb;

	nega = (uint32_t)((uint64_t)cca >> 63);
	negb = (uint32_t)((uint64_t)ccb >> 63);
	zint_negate(a, len, nega);
	zint_negate(b, len, negb);
	return nega | (negb << 1);
}

/*
 * Finish modular reduction. Rules on input parameters:
 *
 *   if neg = 1, then -m <= a < 0
 *   if neg = 0, then 0 <= a < 2*m
 *
 * If neg = 0, then the top word of a[] is allowed to use 32 bits.
 *
 * Modulus m must be odd.
 */
static void
zint_finish_mod(uint32_t *a, size_t len, const uint32_t *m, uint32_t neg)
{
	size_t u;
	uint32_t cc, xm, ym;

	/*
	 * First pass: compare a (assumed nonnegative) with m. Note that
	 * if the top word uses 32 bits, subtracting m must yield a
	 * value less than 2^31 since a < 2*m.
	 */
	cc = 0;
	for (u = 0; u < len; u ++) {
		cc = (a[u] - m[u] - cc) >> 31;
	}

	/*
	 * If neg = 1 then we must add m (regardless of cc)
	 * If neg = 0 and cc = 0 then we must subtract m
	 * If neg = 0 and cc = 1 then we must do nothing
	 *
	 * In the loop below, we conditionally subtract either m or -m
	 * from a. Word xm is a word of m (if neg = 0) or -m (if neg = 1);
	 * but if neg = 0 and cc = 1, then ym = 0 and it forces mw to 0.
	 */
	xm = -neg >> 1;
	ym = -(neg | (1 - cc));
	cc = neg;
	for (u = 0; u < len; u ++) {
		uint32_t aw, mw;

		aw = a[u];
		mw = (m[u] ^ xm) & ym;
		aw = aw - mw - cc;
		a[u] = aw & 0x7FFFFFFF;
		cc = aw >> 31;
	}
}

/*
 * Replace a with (a*xa+b*xb)/(2^31) mod m, and b with
 * (a*ya+b*yb)/(2^31) mod m. Modulus m must be odd; m0i = -1/m[0] mod 2^31.
 */
static void
zint_co_reduce_mod(uint32_t *a, uint32_t *b, const uint32_t *m, size_t len,
	uint32_t m0i, int64_t xa, int64_t xb, int64_t ya, int64_t yb)
{
	size_t u;
	int64_t cca, ccb;
	uint32_t fa, fb;

	/*
	 * These are actually four combined Montgomery multiplications.
	 */
	cca = 0;
	ccb = 0;
	fa = ((a[0] * (uint32_t)xa + b[0] * (uint32_t)xb) * m0i) & 0x7FFFFFFF;
	fb = ((a[0] * (uint32_t)ya + b[0] * (uint32_t)yb) * m0i) & 0x7FFFFFFF;
	for (u = 0; u < len; u ++) {
		uint32_t wa, wb;
		uint64_t za, zb;

		wa = a[u];
		wb = b[u];
		za = wa * (uint64_t)xa + wb * (uint64_t)xb
			+ m[u] * (uint64_t)fa + (uint64_t)cca;
		zb = wa * (uint64_t)ya + wb * (uint64_t)yb
			+ m[u] * (uint64_t)fb + (uint64_t)ccb;
		if (u > 0) {
			a[u - 1] = (uint32_t)za & 0x7FFFFFFF;
			b[u - 1] = (uint32_t)zb & 0x7FFFFFFF;
		}
		cca = *(int64_t *)&za >> 31;
		ccb = *(int64_t *)&zb >> 31;
	}
	a[len - 1] = (uint32_t)cca;
	b[len - 1] = (uint32_t)ccb;

	/*
	 * At this point:
	 *   -m <= a < 2*m
	 *   -m <= b < 2*m
	 * (this is a case of Montgomery reduction)
	 * The top words of 'a' and 'b' may have a 32-th bit set.
	 * We want to add or subtract the modulus, as required.
	 */
	zint_finish_mod(a, len, m, (uint32_t)((uint64_t)cca >> 63));
	zint_finish_mod(b, len, m, (uint32_t)((uint64_t)ccb >> 63));
}

/*
 * Compute a GCD between two positive big integers x and y. The two
 * integers must be odd. Returned value is 1 if the GCD is 1, 0
 * otherwise. When 1 is returned, arrays u and v are filled with values
 * such that:
 *   0 <= u <= y
 *   0 <= v <= x
 *   x*u - y*v = 1
 * x[] and y[] are unmodified. Both input values must have the same
 * encoded length. Temporary array must be large enough to accommodate 4
 * extra values of that length. Arrays u, v and tmp may not overlap with
 * each other, or with either x or y.
 */
static int
zint_bezout(uint32_t *restrict u, uint32_t *restrict v,
	const uint32_t *restrict x, const uint32_t *restrict y,
	size_t len, uint32_t *restrict tmp)
{
	/*
	 * Algorithm is an extended binary GCD. We maintain 6 values
	 * a, b, u0, u1, v0 and v1 with the following invariants:
	 *
	 *  a = x*u0 - y*v0
	 *  b = x*u1 - y*v1
	 *  0 <= a <= x
	 *  0 <= b <= y
	 *  0 <= u0 < y
	 *  0 <= v0 < x
	 *  0 <= u1 <= y
	 *  0 <= v1 < x
	 *
	 * Initial values are:
	 *
	 *  a = x   u0 = 1   v0 = 0
	 *  b = y   u1 = y   v1 = x-1
	 *
	 * Each iteration reduces either a or b, and maintains the
	 * invariants. Algorithm stops when a = b, at which point their
	 * common value is GCD(a,b) and (u0,v0) (or (u1,v1)) contains
	 * the values (u,v) we want to return.
	 *
	 * The formal definition of the algorithm is a sequence of steps:
	 *
	 *  - If a is even, then:
	 *        a <- a/2
	 *        u0 <- u0/2 mod y
	 *        v0 <- v0/2 mod x
	 *
	 *  - Otherwise, if b is even, then:
	 *        b <- b/2
	 *        u1 <- u1/2 mod y
	 *        v1 <- v1/2 mod x
	 *
	 *  - Otherwise, if a > b, then:
	 *        a <- (a-b)/2
	 *        u0 <- (u0-u1)/2 mod y
	 *        v0 <- (v0-v1)/2 mod x
	 *
	 *  - Otherwise:
	 *        b <- (b-a)/2
	 *        u1 <- (u1-u0)/2 mod y
	 *        v1 <- (v1-v0)/2 mod y
	 *
	 * We can show that the operations above preserve the invariants:
	 *
	 *  - If a is even, then u0 and v0 are either both even or both
	 *    odd (since a = x*u0 - y*v0, and x and y are both odd).
	 *    If u0 and v0 are both even, then (u0,v0) <- (u0/2,v0/2).
	 *    Otherwise, (u0,v0) <- ((u0+y)/2,(v0+x)/2). Either way,
	 *    the a = x*u0 - y*v0 invariant is preserved.
	 *
	 *  - The same holds for the case where b is even.
	 *
	 *  - If a and b are odd, and a > b, then:
	 *
	 *      a-b = x*(u0-u1) - y*(v0-v1)
	 *
	 *    In that situation, if u0 < u1, then x*(u0-u1) < 0, but
	 *    a-b > 0; therefore, it must be that v0 < v1, and the
	 *    first part of the update is: (u0,v0) <- (u0-u1+y,v0-v1+x),
	 *    which preserves the invariants. Otherwise, if u0 > u1,
	 *    then u0-u1 >= 1, thus x*(u0-u1) >= x. But a <= x and
	 *    b >= 0, hence a-b <= x. It follows that, in that case,
	 *    v0-v1 >= 0. The first part of the update is then:
	 *    (u0,v0) <- (u0-u1,v0-v1), which again preserves the
	 *    invariants.
	 *
	 *    Either way, once the subtraction is done, the new value of
	 *    a, which is the difference of two odd values, is even,
	 *    and the remaining of this step is a subcase of the
	 *    first algorithm case (i.e. when a is even).
	 *
	 *  - If a and b are odd, and b > a, then the a similar
	 *    argument holds.
	 *
	 * The values a and b start at x and y, respectively. Since x
	 * and y are odd, their GCD is odd, and it is easily seen that
	 * all steps conserve the GCD (GCD(a-b,b) = GCD(a, b);
	 * GCD(a/2,b) = GCD(a,b) if GCD(a,b) is odd). Moreover, either a
	 * or b is reduced by at least one bit at each iteration, so
	 * the algorithm necessarily converges on the case a = b, at
	 * which point the common value is the GCD.
	 *
	 * In the algorithm expressed above, when a = b, the fourth case
	 * applies, and sets b = 0. Since a contains the GCD of x and y,
	 * which are both odd, a must be odd, and subsequent iterations
	 * (if any) will simply divide b by 2 repeatedly, which has no
	 * consequence. Thus, the algorithm can run for more iterations
	 * than necessary; the final GCD will be in a, and the (u,v)
	 * coefficients will be (u0,v0).
	 *
	 *
	 * The presentation above is bit-by-bit. It can be sped up by
	 * noticing that all decisions are taken based on the low bits
	 * and high bits of a and b. We can extract the two top words
	 * and low word of each of a and b, and compute reduction
	 * parameters pa, pb, qa and qb such that the new values for
	 * a and b are:
	 *    a' = (a*pa + b*pb) / (2^31)
	 *    b' = (a*qa + b*qb) / (2^31)
	 * the two divisions being exact. The coefficients are obtained
	 * just from the extracted words, and may be slightly off, requiring
	 * an optional correction: if a' < 0, then we replace pa with -pa
	 * and pb with -pb. Each such step will reduce the total length
	 * (sum of lengths of a and b) by at least 30 bits at each
	 * iteration.
	 */
	uint32_t *u0, *u1, *v0, *v1, *a, *b;
	uint32_t x0i, y0i;
	uint32_t num, rc;
	size_t j;

	if (len == 0) {
		return 0;
	}

	/*
	 * u0 and v0 are the u and v result buffers; the four other
	 * values (u1, v1, a and b) are taken from tmp[].
	 */
	u0 = u;
	v0 = v;
	u1 = tmp;
	v1 = u1 + len;
	a = v1 + len;
	b = a + len;

	/*
	 * We'll need the Montgomery reduction coefficients.
	 */
	x0i = modp_ninv31(x[0]);
	y0i = modp_ninv31(y[0]);

	/*
	 * Initialize a, b, u0, u1, v0 and v1.
	 *  a = x   u0 = 1   v0 = 0
	 *  b = y   u1 = y   v1 = x-1
	 * Note that x is odd, so computing x-1 is easy.
	 */
	memcpy(a, x, len * sizeof *x);
	memcpy(b, y, len * sizeof *y);
	u0[0] = 1;
	memset(u0 + 1, 0, (len - 1) * sizeof *u0);
	memset(v0, 0, len * sizeof *v0);
	memcpy(u1, y, len * sizeof *u1);
	memcpy(v1, x, len * sizeof *v1);
	v1[0] --;

	/*
	 * Each input operand may be as large as 31*len bits, and we
	 * reduce the total length by at least 30 bits at each iteration.
	 */
	for (num = 62 * (uint32_t)len + 30; num >= 30; num -= 30) {
		uint32_t c0, c1;
		uint32_t a0, a1, b0, b1;
		uint64_t a_hi, b_hi;
		uint32_t a_lo, b_lo;
		int64_t pa, pb, qa, qb;
		int i;
		uint32_t r;

		/*
		 * Extract the top words of a and b. If j is the highest
		 * index >= 1 such that a[j] != 0 or b[j] != 0, then we
		 * want (a[j] << 31) + a[j-1] and (b[j] << 31) + b[j-1].
		 * If a and b are down to one word each, then we use
		 * a[0] and b[0].
		 */
		c0 = (uint32_t)-1;
		c1 = (uint32_t)-1;
		a0 = 0;
		a1 = 0;
		b0 = 0;
		b1 = 0;
		j = len;
		while (j -- > 0) {
			uint32_t aw, bw;

			aw = a[j];
			bw = b[j];
			a0 ^= (a0 ^ aw) & c0;
			a1 ^= (a1 ^ aw) & c1;
			b0 ^= (b0 ^ bw) & c0;
			b1 ^= (b1 ^ bw) & c1;
			c1 = c0;
			c0 &= (((aw | bw) + 0x7FFFFFFF) >> 31) - (uint32_t)1;
		}

		/*
		 * If c1 = 0, then we grabbed two words for a and b.
		 * If c1 != 0 but c0 = 0, then we grabbed one word. It
		 * is not possible that c1 != 0 and c0 != 0, because that
		 * would mean that both integers are zero.
		 */
		a1 |= a0 & c1;
		a0 &= ~c1;
		b1 |= b0 & c1;
		b0 &= ~c1;
		a_hi = ((uint64_t)a0 << 31) + a1;
		b_hi = ((uint64_t)b0 << 31) + b1;
		a_lo = a[0];
		b_lo = b[0];

		/*
		 * Compute reduction factors:
		 *
		 *   a' = a*pa + b*pb
		 *   b' = a*qa + b*qb
		 *
		 * such that a' and b' are both multiple of 2^31, but are
		 * only marginally larger than a and b.
		 */
		pa = 1;
		pb = 0;
		qa = 0;
		qb = 1;
		for (i = 0; i < 31; i ++) {
			/*
			 * At each iteration:
			 *
			 *   a <- (a-b)/2 if: a is odd, b is odd, a_hi > b_hi
			 *   b <- (b-a)/2 if: a is odd, b is odd, a_hi <= b_hi
			 *   a <- a/2 if: a is even
			 *   b <- b/2 if: a is odd, b is even
			 *
			 * We multiply a_lo and b_lo by 2 at each
			 * iteration, thus a division by 2 really is a
			 * non-multiplication by 2.
			 */
			uint32_t rt, oa, ob, cAB, cBA, cA;
			uint64_t rz;

			/*
			 * rt = 1 if a_hi > b_hi, 0 otherwise.
			 */
			rz = b_hi - a_hi;
			rt = (uint32_t)((rz ^ ((a_hi ^ b_hi)
				& (a_hi ^ rz))) >> 63);

			/*
			 * cAB = 1 if b must be subtracted from a
			 * cBA = 1 if a must be subtracted from b
			 * cA = 1 if a must be divided by 2
			 *
			 * Rules:
			 *
			 *   cAB and cBA cannot both be 1.
			 *   If a is not divided by 2, b is.
			 */
			oa = (a_lo >> i) & 1;
			ob = (b_lo >> i) & 1;
			cAB = oa & ob & rt;
			cBA = oa & ob & ~rt;
			cA = cAB | (oa ^ 1);

			/*
			 * Conditional subtractions.
			 */
			a_lo -= b_lo & -cAB;
			a_hi -= b_hi & -(uint64_t)cAB;
			pa -= qa & -(int64_t)cAB;
			pb -= qb & -(int64_t)cAB;
			b_lo -= a_lo & -cBA;
			b_hi -= a_hi & -(uint64_t)cBA;
			qa -= pa & -(int64_t)cBA;
			qb -= pb & -(int64_t)cBA;

			/*
			 * Shifting.
			 */
			a_lo += a_lo & (cA - 1);
			pa += pa & ((int64_t)cA - 1);
			pb += pb & ((int64_t)cA - 1);
			a_hi ^= (a_hi ^ (a_hi >> 1)) & -(uint64_t)cA;
			b_lo += b_lo & -cA;
			qa += qa & -(int64_t)cA;
			qb += qb & -(int64_t)cA;
			b_hi ^= (b_hi ^ (b_hi >> 1)) & ((uint64_t)cA - 1);
		}

		/*
		 * Apply the computed parameters to our values. We
		 * may have to correct pa and pb depending on the
		 * returned value of zint_co_reduce() (when a and/or b
		 * had to be negated).
		 */
		r = zint_co_reduce(a, b, len, pa, pb, qa, qb);
		pa -= (pa + pa) & -(int64_t)(r & 1);
		pb -= (pb + pb) & -(int64_t)(r & 1);
		qa -= (qa + qa) & -(int64_t)(r >> 1);
		qb -= (qb + qb) & -(int64_t)(r >> 1);
		zint_co_reduce_mod(u0, u1, y, len, y0i, pa, pb, qa, qb);
		zint_co_reduce_mod(v0, v1, x, len, x0i, pa, pb, qa, qb);
	}

	/*
	 * At that point, array a[] should contain the GCD, and the
	 * results (u,v) should already be set. We check that the GCD
	 * is indeed 1. We also check that the two operands x and y
	 * are odd.
	 */
	rc = a[0] ^ 1;
	for (j = 1; j < len; j ++) {
		rc |= a[j];
	}
	return (int)((1 - ((rc | -rc) >> 31)) & x[0] & y[0]);
}

/*
 * Add k*y*2^sc to x. The result is assumed to fit in the array of
 * size xlen (truncation is applied if necessary).
 * Scale factor 'sc' is provided as sch and scl, such that:
 *   sch = sc / 31
 *   scl = sc % 31
 * xlen MUST NOT be lower than ylen.
 *
 * x[] and y[] are both signed integers, using two's complement for
 * negative values.
 */
static void
zint_add_scaled_mul_small(uint32_t *restrict x, size_t xlen,
	const uint32_t *restrict y, size_t ylen, int32_t k,
	uint32_t sch, uint32_t scl)
{
	size_t u;
	uint32_t ysign, tw;
	int32_t cc;

	if (ylen == 0) {
		return;
	}

	ysign = -(y[ylen - 1] >> 30) >> 1;
	tw = 0;
	cc = 0;
	for (u = sch; u < xlen; u ++) {
		size_t v;
		uint32_t wy, wys, ccu;
		uint64_t z;

		/*
		 * Get the next word of y (scaled).
		 */
		v = u - sch;
		wy = v < ylen ? y[v] : ysign;
		wys = ((wy << scl) & 0x7FFFFFFF) | tw;
		tw = wy >> (31 - scl);

		/*
		 * The expression below does not overflow.
		 */
		z = (uint64_t)((int64_t)wys * (int64_t)k + (int64_t)x[u] + cc);
		x[u] = (uint32_t)z & 0x7FFFFFFF;

		/*
		 * Right-shifting the signed value z would yield
		 * implementation-defined results (arithmetic shift is
		 * not guaranteed). However, we can cast to unsigned,
		 * and get the next carry as an unsigned word. We can
		 * then convert it back to signed by using the guaranteed
		 * fact that 'int32_t' uses two's complement with no
		 * trap representation or padding bit, and with a layout
		 * compatible with that of 'uint32_t'.
		 */
		ccu = (uint32_t)(z >> 31);
		cc = *(int32_t *)&ccu;
	}
}

/*
 * Subtract y*2^sc from x. The result is assumed to fit in the array of
 * size xlen (truncation is applied if necessary).
 * Scale factor 'sc' is provided as sch and scl, such that:
 *   sch = sc / 31
 *   scl = sc % 31
 * xlen MUST NOT be lower than ylen.
 *
 * x[] and y[] are both signed integers, using two's complement for
 * negative values.
 */
static void
zint_sub_scaled(uint32_t *restrict x, size_t xlen,
	const uint32_t *restrict y, size_t ylen, uint32_t sch, uint32_t scl)
{
	size_t u;
	uint32_t ysign, tw;
	uint32_t cc;

	if (ylen == 0) {
		return;
	}

	ysign = -(y[ylen - 1] >> 30) >> 1;
	tw = 0;
	cc = 0;
	for (u = sch; u < xlen; u ++) {
		size_t v;
		uint32_t w, wy, wys;

		/*
		 * Get the next word of y (scaled).
		 */
		v = u - sch;
		wy = v < ylen ? y[v] : ysign;
		wys = ((wy << scl) & 0x7FFFFFFF) | tw;
		tw = wy >> (31 - scl);

		w = x[u] - wys - cc;
		x[u] = w & 0x7FFFFFFF;
		cc = w >> 31;
	}
}

/*
 * Convert a one-word signed big integer into a signed value.
 */
static inline int32_t
zint_one_to_plain(const uint32_t *x)
{
	uint32_t w;

	w = x[0];
	w |= (w & 0x40000000) << 1;
	return *(int32_t *)&w;
}

/* ==================================================================== */

/*
 * Convert a polynomial to floating-point values.
 *
 * Each coefficient has length flen words, and starts fstride words after
 * the previous.
 *
 * IEEE-754 binary64 values can represent values in a finite range,
 * roughly 2^(-1023) to 2^(+1023); thus, if coefficients are too large,
 * they should be "trimmed" by pointing not to the lowest word of each,
 * but upper.
 */
static void
poly_big_to_fp(fpr *d, const uint32_t *f, size_t flen, size_t fstride,
	unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	if (flen == 0) {
		for (u = 0; u < n; u ++) {
			d[u] = fpr_zero;
		}
		return;
	}
	for (u = 0; u < n; u ++, f += fstride) {
		size_t v;
		uint32_t neg, cc, xm;
		fpr x, fsc;

		/*
		 * Get sign of the integer; if it is negative, then we
		 * will load its absolute value instead, and negate the
		 * result.
		 */
		neg = -(f[flen - 1] >> 30);
		xm = neg >> 1;
		cc = neg & 1;
		x = fpr_zero;
		fsc = fpr_one;
		for (v = 0; v < flen; v ++, fsc = fpr_mul(fsc, fpr_ptwo31)) {
			uint32_t w;

			w = (f[v] ^ xm) + cc;
			cc = w >> 31;
			w &= 0x7FFFFFFF;
			w -= (w << 1) & neg;
			x = fpr_add(x, fpr_mul(fpr_of(*(int32_t *)&w), fsc));
		}
		d[u] = x;
	}
}

/*
 * Convert a polynomial to small integers. Source values are supposed
 * to be one-word integers, signed over 31 bits. Returned value is 0
 * if any of the coefficients exceeds the provided limit (in absolute
 * value), or 1 on success.
 *
 * This is not constant-time; this is not a problem here, because on
 * any failure, the NTRU-solving process will be deemed to have failed
 * and the (f,g) polynomials will be discarded.
 */
static int
poly_big_to_small(int8_t *d, const uint32_t *s, int lim, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		int32_t z;

		z = zint_one_to_plain(s + u);
		if (z < -lim || z > lim) {
			return 0;
		}
		d[u] = (int8_t)z;
	}
	return 1;
}

/*
 * Subtract k*f from F, where F, f and k are polynomials modulo X^N+1.
 * Coefficients of polynomial k are small integers (signed values in the
 * -2^31..2^31 range) scaled by 2^sc. Value sc is provided as sch = sc / 31
 * and scl = sc % 31.
 *
 * This function implements the basic quadratic multiplication algorithm,
 * which is efficient in space (no extra buffer needed) but slow at
 * high degree.
 */
static void
poly_sub_scaled(uint32_t *restrict F, size_t Flen, size_t Fstride,
	const uint32_t *restrict f, size_t flen, size_t fstride,
	const int32_t *restrict k, uint32_t sch, uint32_t scl, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		int32_t kf;
		size_t v;
		uint32_t *x;
		const uint32_t *y;

		kf = -k[u];
		x = F + u * Fstride;
		y = f;
		for (v = 0; v < n; v ++) {
			zint_add_scaled_mul_small(
				x, Flen, y, flen, kf, sch, scl);
			if (u + v == n - 1) {
				x = F;
				kf = -kf;
			} else {
				x += Fstride;
			}
			y += fstride;
		}
	}
}

/*
 * Subtract k*f from F. Coefficients of polynomial k are small integers
 * (signed values in the -2^31..2^31 range) scaled by 2^sc. This function
 * assumes that the degree is large, and integers relatively small.
 * The value sc is provided as sch = sc / 31 and scl = sc % 31.
 */
static void
poly_sub_scaled_ntt(uint32_t *restrict F, size_t Flen, size_t Fstride,
	const uint32_t *restrict f, size_t flen, size_t fstride,
	const int32_t *restrict k, uint32_t sch, uint32_t scl, unsigned logn,
	uint32_t *restrict tmp)
{
	uint32_t *gm, *igm, *fk, *t1, *x;
	const uint32_t *y;
	size_t n, u, tlen;

	n = MKN(logn);
	tlen = flen + 1;
	gm = tmp;
	igm = gm + n;
	fk = igm + n;
	t1 = fk + n * tlen;

	/*
	 * Compute k*f in fk[], in RNS notation.
	 */
	for (u = 0; u < tlen; u ++) {
		uint32_t p, p0i, R2, Rx;
		size_t v;

		p = PRIMES[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);
		Rx = modp_Rx((unsigned)flen, p, p0i, R2);
		modp_mkgm2(gm, igm, logn, PRIMES[u].g, p, p0i);

		for (v = 0; v < n; v ++) {
			t1[v] = modp_set(k[v], p);
		}
		modp_NTT2(t1, gm, logn, p, p0i);
		for (v = 0, y = f, x = fk + u;
			v < n; v ++, y += fstride, x += tlen)
		{
			*x = zint_mod_small_signed(y, flen, p, p0i, R2, Rx);
		}
		modp_NTT2_ext(fk + u, tlen, gm, logn, p, p0i);
		for (v = 0, x = fk + u; v < n; v ++, x += tlen) {
			*x = modp_montymul(
				modp_montymul(t1[v], *x, p, p0i), R2, p, p0i);
		}
		modp_iNTT2_ext(fk + u, tlen, igm, logn, p, p0i);
	}

	/*
	 * Rebuild k*f.
	 */
	zint_rebuild_CRT(fk, tlen, tlen, n, 1, t1);

	/*
	 * Subtract k*f, scaled, from F.
	 */
	for (u = 0, x = F, y = fk; u < n; u ++, x += Fstride, y += tlen) {
		zint_sub_scaled(x, Flen, y, tlen, sch, scl);
	}
}

/* ==================================================================== */

/*
 * The following average lengths, in bits, have been measured on thousands
 * of random keys (fg = max length of the absolute value of coefficients
 * of f and g at that depth; FG = idem for the unreduced F and G; for the
 * maximum depth, F and G are the output of binary GCD;
 * for each value, the average and standard deviation are provided).
 *
 * Binary case:
 *    depth:  9    fg: 2369.72 (24.19)    FG:  2367.83 (24.18)
 *    depth:  8    fg: 1184.95 (12.09)    FG:  3535.44 (32.19)
 *    depth:  7    fg:  598.30 ( 7.09)    FG:  1772.94 (16.40)
 *    depth:  6    fg:  302.97 ( 4.17)    FG:  893.58 (9.30)
 *    depth:  5    fg:  153.62 ( 2.40)    FG:  451.73 (5.47)
 *    depth:  4    fg:   77.63 ( 1.31)    FG:  228.97 (3.18)
 *    depth:  3    fg:   38.68 ( 0.67)    FG:  116.22 (1.76)
 *    depth:  2    fg:   18.62 ( 0.49)    FG:  58.89 (0.95)
 *    depth:  1    fg:    8.05 ( 0.21)    FG:  29.77 (0.52)
 *    depth:  0    fg:    3.00 ( 0.04)    FG:  15.00 (0.14)
 *
 * Integers are actually represented either in binary notation over
 * 31-bit words (signed, using two's complement), or in RNS, modulo
 * many small primes. These small primes are close to, but slightly
 * lower than, 2^31. Use of RNS loses less than two bits, even for
 * the largest values.
 *
 * IMPORTANT: if these values are modified, then the temporary buffer
 * sizes (HAWK_KEYGEN_TEMP_*, in inner.h) must be recomputed
 * accordingly.
 *
 * Note, to get the values in MAX_BL_SMALL and MAX_BL_LARGE,
 * take ceil( (avg + 6 stddev) / 31 ) to arrive at the number of ints used to
 * represent a number at a certain depth.
*/

// NIST-1 parameter set
static const size_t MAX_BL_SMALL_512[10] = { 1, 1, 1, 2, 3,  6, 11, 21,  42, 83 };
static const size_t MAX_BL_LARGE_512[ 9] = { 2, 2, 3, 5, 9, 16, 31, 62, 123 };
// NIST-5 parameter set
static const size_t MAX_BL_SMALL_1024[11] = { 1, 1, 1, 2,  4,  7, 13, 25,  49,  96, 192 };
static const size_t MAX_BL_LARGE_1024[10] = { 2, 2, 3, 5, 10, 19, 37, 73, 144, 286 };

/*
 * Average and standard deviation for the maximum size (in bits) of
 * coefficients of (f,g), depending on depth. These values are used
 * to compute bounds for Babai's reduction.
 */

static const struct {
        int avg, std;
} BITLENGTH_512[10] = {
	// NIST-1 parameter set
	{ 3, 0 },
	{ 8, 1 },
	{ 19, 1 },
	{ 39, 1 },
	{ 78, 2 },
	{ 155, 3 },
	{ 307, 4 },
	{ 606, 8 },
	{ 1201, 13 },
	{ 2402, 25 }
}, BITLENGTH_1024[11] = {
	// NIST-5 parameter set
	{ 3, 1 },
	{ 10, 1 },
	{ 22, 1 },
	{ 46, 1 },
	{ 93, 2 },
	{ 186, 3 },
	{ 368, 4 },
	{ 728, 6 },
	{ 1444, 10 },
	{ 2870, 16 },
	{ 5741, 32 }
};



/*
 * Minimal recursion depth at which we rebuild intermediate values
 * when reconstructing f and g.
 */
#define DEPTH_INT_FG   4

/*
 * Align (upwards) the provided 'data' pointer with regards to 'base'
 * so that the offset is a multiple of the size of 'fpr'.
 */
static fpr *
align_fpr(void *base, void *data)
{
	uint8_t *cb, *cd;
	size_t k, km;

	cb = (uint8_t *)base;
	cd = (uint8_t *)data;
	k = (size_t)(cd - cb);
	km = k % sizeof(fpr);
	if (km) {
		k += (sizeof(fpr)) - km;
	}
	return (fpr *)(cb + k);
}

/*
 * Align (upwards) the provided 'data' pointer with regards to 'base'
 * so that the offset is a multiple of the size of 'uint32_t'.
 */
static uint32_t *
align_u32(void *base, void *data)
{
	uint8_t *cb, *cd;
	size_t k, km;

	cb = (uint8_t *)base;
	cd = (uint8_t *)data;
	k = (size_t)(cd - cb);
	km = k % sizeof(uint32_t);
	if (km) {
		k += (sizeof(uint32_t)) - km;
	}
	return (uint32_t *)(cb + k);
}

// =============================================================================

/*
 * Table below incarnates a discrete Gaussian distribution:
 *    D(x) = exp(-(x^2)/(2*sigma^2))
 * where sigma = 1.500 for HAWK-512 and sigma = 2.000 for HAWK-1024.
 * Element k (k >= 0) contains P(|X| >= k+1) scaled up by 2^63.
 * To generate the values in the table below, run `sage code/renyi.sage`.
 */

/*
 * RD_{513}(keygen) = 1 + 6.859024E-20 < 1 + 2^-63
 */
static const uint64_t gauss_keygen_512[13] = {
	6770309987939008324u, 2841792919453817158u, 824825004081786282u,
	160853309707784581u, 20707417942076380u, 1740733985516594u,
	94912702842187u, 3342501151111u, 75826385177u, 1106214542u, 10367241u,
	62372u, 240u
};

/*
 * RD_{513}(keygen) = 1 + 6.189417E-20 < 1 + 2^-63
 */
static const uint64_t gauss_keygen_1024[18] = {
	7383575500167950194u, 4136346010143971104u, 1904559995876614063u,
	709971025731164996u, 211992254950663721u, 50322218324695908u,
	9445631610447647u, 1396554443915520u, 162188481909715u, 14764491139957u,
	1051923941110u, 58588904910u, 2548776934u, 86545550u, 2292626u, 47361u,
	762u, 9u
};

static int8_t mkgauss_keygen_512(prng *rng)
{
	uint64_t r;
	uint8_t v, k, neg;

	/*
	 * Generate a 64 bit value.
	 */
	r = prng_get_u64(rng);

	/*
	 * Get the sign bit, and unset this bit in r.
	 */
	neg = (uint8_t)(r >> 63);
	r &= ~((uint64_t)1u << 63);

	v = 0;
	for (k = 0; k < 13; k++) {
		/*
		 * Add 1 iff r < gauss_keygen[k].
		 */
		v += (uint8_t)((uint64_t)(r - gauss_keygen_512[k]) >> 63);
	}

	/*
	 * Apply the sign ('neg' flag). If neg = 0, this has no effect.
	 * However, if neg = 1, this changes v into -v = (~v) + 1.
	 */
	v = (v ^ -neg) + neg;
	return *(int8_t *)&v;
}

static int8_t mkgauss_keygen_1024(prng *rng)
{
	uint64_t r;
	uint8_t v, k, neg;

	/*
	 * Generate a 64 bit value.
	 */
	r = prng_get_u64(rng);

	/*
	 * Get the sign bit, and unset this bit in r.
	 */
	neg = (uint8_t)(r >> 63);
	r &= ~((uint64_t)1u << 63);

	v = 0;
	for (k = 0; k < 18; k++) {
		/*
		 * Add 1 iff r < gauss_keygen[k].
		 */
		v += (uint8_t)((uint64_t)(r - gauss_keygen_1024[k]) >> 63);
	}

	/*
	 * Apply the sign ('neg' flag). If neg = 0, this has no effect.
	 * However, if neg = 1, this changes v into -v = (~v) + 1.
	 */
	v = (v ^ -neg) + neg;
	return *(int8_t *)&v;
}


/*
 * Generate a random polynomial with a Gaussian distribution and return the XOR
 * of coefficients modulo 2. The algebraic norm of the polynomial is odd if and
 * only if the sum is odd.
 */
static uint8_t
poly_small_mkgauss(prng *rng, int8_t *f, unsigned logn)
{
	size_t n, u;
	uint8_t mod2;

	n = MKN(logn);
	mod2 = 0u;

	if (logn == 10) {
		for (u = 0; u < n; u++) {
			f[u] = mkgauss_keygen_1024(rng);
			mod2 ^= (uint8_t)f[u];
		}
	} else {
		for (u = 0; u < n; u++) {
			f[u] = mkgauss_keygen_512(rng);
			mod2 ^= (uint8_t)f[u];
		}
	}
	return mod2;
}

/*
 * Input: f,g of degree N = 2^logn; 'depth' is used only to get their
 * individual length.
 *
 * Output: f',g' of degree N/2, with the length for 'depth+1'.
 *
 * Values are in RNS; input and/or output may also be in NTT.
 */
static void
make_fg_step(uint32_t *data, unsigned logn, unsigned depth,
	int in_ntt, int out_ntt)
{
	size_t n, hn, u;
	size_t slen, tlen;
	uint32_t *fd, *gd, *fs, *gs, *gm, *igm, *t1;

	n = MKN(logn);
	hn = n >> 1;

	if (logn + depth == 10) {
		slen = MAX_BL_SMALL_1024[depth];
		tlen = MAX_BL_SMALL_1024[depth + 1];
	} else {
		slen = MAX_BL_SMALL_512[depth];
		tlen = MAX_BL_SMALL_512[depth + 1];
	}

	/*
	 * Prepare room for the result.
	 */
	fd = data;
	gd = fd + hn * tlen;
	fs = gd + hn * tlen;
	gs = fs + n * slen;
	gm = gs + n * slen;
	igm = gm + n;
	t1 = igm + n;
	memmove(fs, data, 2 * n * slen * sizeof *data);

	/*
	 * First slen words: we use the input values directly, and apply
	 * inverse NTT as we go.
	 */
	for (u = 0; u < slen; u ++) {
		uint32_t p, p0i, R2;
		size_t v;
		uint32_t *x;

		p = PRIMES[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);
		modp_mkgm2(gm, igm, logn, PRIMES[u].g, p, p0i);

		for (v = 0, x = fs + u; v < n; v ++, x += slen) {
			t1[v] = *x;
		}
		if (!in_ntt) {
			modp_NTT2(t1, gm, logn, p, p0i);
		}
		for (v = 0, x = fd + u; v < hn; v ++, x += tlen) {
			uint32_t w0, w1;

			w0 = t1[(v << 1) + 0];
			w1 = t1[(v << 1) + 1];
			*x = modp_montymul(
				modp_montymul(w0, w1, p, p0i), R2, p, p0i);
		}
		if (in_ntt) {
			modp_iNTT2_ext(fs + u, slen, igm, logn, p, p0i);
		}

		for (v = 0, x = gs + u; v < n; v ++, x += slen) {
			t1[v] = *x;
		}
		if (!in_ntt) {
			modp_NTT2(t1, gm, logn, p, p0i);
		}
		for (v = 0, x = gd + u; v < hn; v ++, x += tlen) {
			uint32_t w0, w1;

			w0 = t1[(v << 1) + 0];
			w1 = t1[(v << 1) + 1];
			*x = modp_montymul(
				modp_montymul(w0, w1, p, p0i), R2, p, p0i);
		}
		if (in_ntt) {
			modp_iNTT2_ext(gs + u, slen, igm, logn, p, p0i);
		}

		if (!out_ntt) {
			modp_iNTT2_ext(fd + u, tlen, igm, logn - 1, p, p0i);
			modp_iNTT2_ext(gd + u, tlen, igm, logn - 1, p, p0i);
		}
	}

	/*
	 * Since the fs and gs words have been de-NTTized, we can use the
	 * CRT to rebuild the values.
	 */
	zint_rebuild_CRT(fs, slen, slen, n, 1, gm);
	zint_rebuild_CRT(gs, slen, slen, n, 1, gm);

	/*
	 * Remaining words: use modular reductions to extract the values.
	 */
	for (u = slen; u < tlen; u ++) {
		uint32_t p, p0i, R2, Rx;
		size_t v;
		uint32_t *x;

		p = PRIMES[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);
		Rx = modp_Rx((unsigned)slen, p, p0i, R2);
		modp_mkgm2(gm, igm, logn, PRIMES[u].g, p, p0i);
		for (v = 0, x = fs; v < n; v ++, x += slen) {
			t1[v] = zint_mod_small_signed(x, slen, p, p0i, R2, Rx);
		}
		modp_NTT2(t1, gm, logn, p, p0i);
		for (v = 0, x = fd + u; v < hn; v ++, x += tlen) {
			uint32_t w0, w1;

			w0 = t1[(v << 1) + 0];
			w1 = t1[(v << 1) + 1];
			*x = modp_montymul(
				modp_montymul(w0, w1, p, p0i), R2, p, p0i);
		}
		for (v = 0, x = gs; v < n; v ++, x += slen) {
			t1[v] = zint_mod_small_signed(x, slen, p, p0i, R2, Rx);
		}
		modp_NTT2(t1, gm, logn, p, p0i);
		for (v = 0, x = gd + u; v < hn; v ++, x += tlen) {
			uint32_t w0, w1;

			w0 = t1[(v << 1) + 0];
			w1 = t1[(v << 1) + 1];
			*x = modp_montymul(
				modp_montymul(w0, w1, p, p0i), R2, p, p0i);
		}

		if (!out_ntt) {
			modp_iNTT2_ext(fd + u, tlen, igm, logn - 1, p, p0i);
			modp_iNTT2_ext(gd + u, tlen, igm, logn - 1, p, p0i);
		}
	}
}

/*
 * Compute f and g at a specific depth, in RNS notation.
 *
 * Returned values are stored in the data[] array, at slen words per integer.
 *
 * Conditions:
 *   0 <= depth <= logn
 *
 * Space use in data[]: enough room for any two successive values (f', g',
 * f and g).
 */
static void
make_fg(uint32_t *data, const int8_t *f, const int8_t *g,
	unsigned logn, unsigned depth, int out_ntt)
{
	size_t n, u;
	uint32_t *ft, *gt, p;
	unsigned d;

	n = MKN(logn);
	ft = data;
	gt = ft + n;
	p = PRIMES[0].p;
	for (u = 0; u < n; u ++) {
		ft[u] = modp_set(f[u], p);
		gt[u] = modp_set(g[u], p);
	}

	if (depth == 0 && out_ntt) {
		uint32_t *gm, *igm, p0i;

		p0i = modp_ninv31(p);
		gm = gt + n;
		igm = gm + MKN(logn);
		modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);
		modp_NTT2(ft, gm, logn, p, p0i);
		modp_NTT2(gt, gm, logn, p, p0i);
		return;
	}

	for (d = 0; d < depth; d ++) {
		make_fg_step(data, logn - d, d, d != 0, (d + 1) < depth || out_ntt);
	}
}

/*
 * Solving the NTRU equation for q = 1, deepest level: compute the resultants
 * of f and g with X^N+1, and use binary GCD. The F and G values are returned
 * in tmp[].
 *
 * Returned value: 1 on success, 0 on error.
 */
static int
solve_NTRU_deepest(unsigned logn_top,
	const int8_t *f, const int8_t *g, uint32_t *tmp)
{
	size_t len;
	uint32_t *Fp, *Gp, *fp, *gp, *t1;

	if (logn_top == 10) {
		len = MAX_BL_SMALL_1024[logn_top];
	} else {
		len = MAX_BL_SMALL_512[logn_top];
	}

	Fp = tmp;
	Gp = Fp + len;
	fp = Gp + len;
	gp = fp + len;
	t1 = gp + len;

	make_fg(fp, f, g, logn_top, logn_top, 0);

	/*
	 * We use the CRT to rebuild the resultants as big integers.
	 * There are two such big integers. The resultants are always
	 * nonnegative.
	 */
	zint_rebuild_CRT(fp, len, len, 2, 0, t1);

	/*
	 * Apply the binary GCD. The zint_bezout() function works only
	 * if both inputs are odd.
	 *
	 * We can test on the result and return 0 because that would
	 * imply failure of the NTRU solving equation, and the (f,g)
	 * values will be abandoned in that case.
	 */
	return zint_bezout(Gp, Fp, fp, gp, len, t1);
}

static int
solve_NTRU_intermediate(unsigned logn_top,
	const int8_t *f, const int8_t *g, unsigned depth, uint32_t *tmp)
{
	/*
	 * In this function, 'logn' is the log2 of the degree for
	 * this step. If N = 2^logn, then:
	 *  - the F and G values already in fk->tmp (from the deeper
	 *    levels) have degree N/2;
	 *  - this function should return F and G of degree N.
	 */
	unsigned logn;
	size_t n, hn, slen, dlen, llen, rlen, FGlen, u;
	uint32_t *Fd, *Gd, *Ft, *Gt, *ft, *gt, *t1;
	fpr *rt1, *rt2, *rt3, *rt4, *rt5;
	int scale_fg, minbl_fg, maxbl_fg, maxbl_FG, scale_k;
	uint32_t *x, *y;
	int32_t *k;

	logn = logn_top - depth;
	n = MKN(logn);
	hn = n >> 1;

	/*
	 * slen = size for our input f and g; also size of the reduced
	 *        F and G we return (degree N)
	 *
	 * dlen = size of the F and G obtained from the deeper level
	 *        (degree N/2 or N/3)
	 *
	 * llen = size for intermediary F and G before reduction (degree N)
	 *
	 * We build our non-reduced F and G as two independent halves each,
	 * of degree N/2 (F = F0 + X*F1, G = G0 + X*G1).
	 */
	if (logn_top == 10) {
		slen = MAX_BL_SMALL_1024[depth];
		dlen = MAX_BL_SMALL_1024[depth + 1];
		llen = MAX_BL_LARGE_1024[depth];
	} else {
		slen = MAX_BL_SMALL_512[depth];
		dlen = MAX_BL_SMALL_512[depth + 1];
		llen = MAX_BL_LARGE_512[depth];
	}

	/*
	 * Fd and Gd are the F and G from the deeper level.
	 */
	Fd = tmp;
	Gd = Fd + dlen * hn;

	/*
	 * Compute the input f and g for this level. Note that we get f
	 * and g in RNS + NTT representation.
	 */
	ft = Gd + dlen * hn;
	make_fg(ft, f, g, logn_top, depth, 1);

	/*
	 * Move the newly computed f and g to make room for our candidate
	 * F and G (unreduced).
	 */
	Ft = tmp;
	Gt = Ft + n * llen;
	t1 = Gt + n * llen;
	memmove(t1, ft, 2 * n * slen * sizeof *ft);
	ft = t1;
	gt = ft + slen * n;
	t1 = gt + slen * n;

	/*
	 * Move Fd and Gd _after_ f and g.
	 */
	memmove(t1, Fd, 2 * hn * dlen * sizeof *Fd);
	Fd = t1;
	Gd = Fd + hn * dlen;

	/*
	 * We reduce Fd and Gd modulo all the small primes we will need,
	 * and store the values in Ft and Gt (only n/2 values in each).
	 */
	for (u = 0; u < llen; u ++) {
		uint32_t p, p0i, R2p, Rx;
		size_t v;
		uint32_t *xs, *ys, *xd, *yd;

		p = PRIMES[u].p;
		p0i = modp_ninv31(p);
		R2p = modp_R2(p, p0i);
		Rx = modp_Rx((unsigned)dlen, p, p0i, R2p);
		for (v = 0, xs = Fd, ys = Gd, xd = Ft + u, yd = Gt + u;
			v < hn;
			v ++, xs += dlen, ys += dlen, xd += llen, yd += llen)
		{
			*xd = zint_mod_small_signed(xs, dlen, p, p0i, R2p, Rx);
			*yd = zint_mod_small_signed(ys, dlen, p, p0i, R2p, Rx);
		}
	}

	/*
	 * We do not need Fd and Gd after that point.
	 */

	/*
	 * Compute our F and G modulo sufficiently many small primes.
	 */
	for (u = 0; u < llen; u ++) {
		uint32_t p, p0i, R2p;
		uint32_t *gm, *igm, *fx, *gx, *Fp, *Gp;
		size_t v;

		/*
		 * All computations are done modulo p.
		 */
		p = PRIMES[u].p;
		p0i = modp_ninv31(p);
		R2p = modp_R2(p, p0i);

		/*
		 * If we processed slen words, then f and g have been
		 * de-NTTized, and are in RNS; we can rebuild them.
		 */
		if (u == slen) {
			zint_rebuild_CRT(ft, slen, slen, n, 1, t1);
			zint_rebuild_CRT(gt, slen, slen, n, 1, t1);
		}

		gm = t1;
		igm = gm + n;
		fx = igm + n;
		gx = fx + n;

		modp_mkgm2(gm, igm, logn, PRIMES[u].g, p, p0i);

		if (u < slen) {
			for (v = 0, x = ft + u, y = gt + u;
				v < n; v ++, x += slen, y += slen)
			{
				fx[v] = *x;
				gx[v] = *y;
			}
			modp_iNTT2_ext(ft + u, slen, igm, logn, p, p0i);
			modp_iNTT2_ext(gt + u, slen, igm, logn, p, p0i);
		} else {
			uint32_t Rx;

			Rx = modp_Rx((unsigned)slen, p, p0i, R2p);
			for (v = 0, x = ft, y = gt;
				v < n; v ++, x += slen, y += slen)
			{
				fx[v] = zint_mod_small_signed(x, slen,
					p, p0i, R2p, Rx);
				gx[v] = zint_mod_small_signed(y, slen,
					p, p0i, R2p, Rx);
			}
			modp_NTT2(fx, gm, logn, p, p0i);
			modp_NTT2(gx, gm, logn, p, p0i);
		}

		/*
		 * Get F' and G' modulo p and in NTT representation
		 * (they have degree n/2). These values were computed in
		 * a previous step, and stored in Ft and Gt.
		 */
		Fp = gx + n;
		Gp = Fp + hn;
		for (v = 0, x = Ft + u, y = Gt + u;
			v < hn; v ++, x += llen, y += llen)
		{
			Fp[v] = *x;
			Gp[v] = *y;
		}
		modp_NTT2(Fp, gm, logn - 1, p, p0i);
		modp_NTT2(Gp, gm, logn - 1, p, p0i);

		/*
		 * Compute our F and G modulo p.
		 *
		 * General case:
		 *
		 *   we divide degree by d = 2 or 3
		 *   f'(x^d) = N(f)(x^d) = f * adj(f)
		 *   g'(x^d) = N(g)(x^d) = g * adj(g)
		 *   f'*G' - g'*F' = 1
		 *   F = F'(x^d) * adj(g)
		 *   G = G'(x^d) * adj(f)
		 *
		 * We compute things in the NTT. We group roots of phi
		 * such that all roots x in a group share the same x^d.
		 * If the roots in a group are x_1, x_2... x_d, then:
		 *
		 *   N(f)(x_1^d) = f(x_1)*f(x_2)*...*f(x_d)
		 *
		 * Thus, we have:
		 *
		 *   G(x_1) = f(x_2)*f(x_3)*...*f(x_d)*G'(x_1^d)
		 *   G(x_2) = f(x_1)*f(x_3)*...*f(x_d)*G'(x_1^d)
		 *   ...
		 *   G(x_d) = f(x_1)*f(x_2)*...*f(x_{d-1})*G'(x_1^d)
		 *
		 * In all cases, we can thus compute F and G in NTT
		 * representation by a few simple multiplications.
		 * Moreover, in our chosen NTT representation, roots
		 * from the same group are consecutive in RAM.
		 */
		for (v = 0, x = Ft + u, y = Gt + u; v < hn;
			v ++, x += (llen << 1), y += (llen << 1))
		{
			uint32_t ftA, ftB, gtA, gtB;
			uint32_t mFp, mGp;

			ftA = fx[(v << 1) + 0];
			ftB = fx[(v << 1) + 1];
			gtA = gx[(v << 1) + 0];
			gtB = gx[(v << 1) + 1];
			mFp = modp_montymul(Fp[v], R2p, p, p0i);
			mGp = modp_montymul(Gp[v], R2p, p, p0i);
			x[0] = modp_montymul(gtB, mFp, p, p0i);
			x[llen] = modp_montymul(gtA, mFp, p, p0i);
			y[0] = modp_montymul(ftB, mGp, p, p0i);
			y[llen] = modp_montymul(ftA, mGp, p, p0i);
		}
		modp_iNTT2_ext(Ft + u, llen, igm, logn, p, p0i);
		modp_iNTT2_ext(Gt + u, llen, igm, logn, p, p0i);
	}

	/*
	 * Rebuild F and G with the CRT.
	 */
	zint_rebuild_CRT(Ft, llen, llen, n, 1, t1);
	zint_rebuild_CRT(Gt, llen, llen, n, 1, t1);

	/*
	 * At that point, Ft, Gt, ft and gt are consecutive in RAM (in that
	 * order).
	 */

	/*
	 * Apply Babai reduction to bring back F and G to size slen.
	 *
	 * We use the FFT to compute successive approximations of the
	 * reduction coefficient. We first isolate the top bits of
	 * the coefficients of f and g, and convert them to floating
	 * point; with the FFT, we compute adj(f), adj(g), and
	 * 1/(f*adj(f)+g*adj(g)).
	 *
	 * Then, we repeatedly apply the following:
	 *
	 *   - Get the top bits of the coefficients of F and G into
	 *     floating point, and use the FFT to compute:
	 *        (F*adj(f)+G*adj(g))/(f*adj(f)+g*adj(g))
	 *
	 *   - Convert back that value into normal representation, and
	 *     round it to the nearest integers, yielding a polynomial k.
	 *     Proper scaling is applied to f, g, F and G so that the
	 *     coefficients fit on 32 bits (signed).
	 *
	 *   - Subtract k*f from F and k*g from G.
	 *
	 * Under normal conditions, this process reduces the size of F
	 * and G by some bits at each iteration. For constant-time
	 * operation, we do not want to measure the actual length of
	 * F and G; instead, we do the following:
	 *
	 *   - f and g are converted to floating-point, with some scaling
	 *     if necessary to keep values in the representable range.
	 *
	 *   - For each iteration, we _assume_ a maximum size for F and G,
	 *     and use the values at that size. If we overreach, then
	 *     we get zeros, which is harmless: the resulting coefficients
	 *     of k will be 0 and the value won't be reduced.
	 *
	 *   - We conservatively assume that F and G will be reduced by
	 *     at least 25 bits at each iteration.
	 *
	 * Even when reaching the bottom of the reduction, reduction
	 * coefficient will remain low. If it goes out-of-range, then
	 * something wrong occurred and the whole NTRU solving fails.
	 */

	/*
	 * Memory layout:
	 *  - We need to compute and keep adj(f), adj(g), and
	 *    1/(f*adj(f)+g*adj(g)) (sizes N, N and N/2 fp numbers,
	 *    respectively).
	 *  - At each iteration we need two extra fp buffer (N fp values),
	 *    and produce a k (N 32-bit words). k will be shared with one
	 *    of the fp buffers.
	 *  - To compute k*f and k*g efficiently (with the NTT), we need
	 *    some extra room; we reuse the space of the temporary buffers.
	 *
	 * Arrays of 'fpr' are obtained from the temporary array itself.
	 * We ensure that the base is at a properly aligned offset (the
	 * source array tmp[] is supposed to be already aligned).
	 */

	rt3 = align_fpr(tmp, t1);
	rt4 = rt3 + n;
	rt5 = rt4 + n;
	rt1 = rt5 + (n >> 1);
	k = (int32_t *)align_u32(tmp, rt1);
	rt2 = align_fpr(tmp, k + n);
	if (rt2 < (rt1 + n)) {
		rt2 = rt1 + n;
	}
	t1 = (uint32_t *)k + n;

	/*
	 * Get f and g into rt3 and rt4 as floating-point approximations.
	 *
	 * We need to "scale down" the floating-point representation of
	 * coefficients when they are too big. We want to keep the value
	 * below 2^310 or so. Thus, when values are larger than 10 words,
	 * we consider only the top 10 words. Array lengths have been
	 * computed so that average maximum length will fall in the
	 * middle or the upper half of these top 10 words.
	 */

	rlen = (slen > 10) ? 10 : slen;
	poly_big_to_fp(rt3, ft + slen - rlen, rlen, slen, logn);
	poly_big_to_fp(rt4, gt + slen - rlen, rlen, slen, logn);

	/*
	 * Values in rt3 and rt4 are downscaled by 2^(scale_fg).
	 */
	scale_fg = 31 * (int)(slen - rlen);

	/*
	 * Estimated boundaries for the maximum size (in bits) of the
	 * coefficients of (f,g). We use the measured average, and
	 * allow for a deviation of at most six times the standard
	 * deviation.
	 */
	if (logn_top == 10) {
		minbl_fg = BITLENGTH_1024[depth].avg - 6 * BITLENGTH_1024[depth].std;
		maxbl_fg = BITLENGTH_1024[depth].avg + 6 * BITLENGTH_1024[depth].std;
	} else {
		minbl_fg = BITLENGTH_512[depth].avg - 6 * BITLENGTH_512[depth].std;
		maxbl_fg = BITLENGTH_512[depth].avg + 6 * BITLENGTH_512[depth].std;
	}

	/*
	 * Compute 1/(f*adj(f)+g*adj(g)) in rt5. We also keep adj(f)
	 * and adj(g) in rt3 and rt4, respectively.
	 */

	Zf(FFT)(rt3, logn);
	Zf(FFT)(rt4, logn);
	Zf(poly_invnorm2_fft)(rt5, rt3, rt4, logn);
	Zf(poly_adj_fft)(rt3, logn);
	Zf(poly_adj_fft)(rt4, logn);

	/*
	 * Reduce F and G repeatedly.
	 *
	 * The expected maximum bit length of coefficients of F and G
	 * is kept in maxbl_FG, with the corresponding word length in
	 * FGlen.
	 */
	FGlen = llen;
	maxbl_FG = 31 * (int)llen;

	/*
	 * Each reduction operation computes the reduction polynomial
	 * "k". We need that polynomial to have coefficients that fit
	 * on 32-bit signed integers, with some scaling; thus, we use
	 * a descending sequence of scaling values, down to zero.
	 *
	 * The size of the coefficients of k is (roughly) the difference
	 * between the size of the coefficients of (F,G) and the size
	 * of the coefficients of (f,g). Thus, the maximum size of the
	 * coefficients of k is, at the start, maxbl_FG - minbl_fg;
	 * this is our starting scale value for k.
	 *
	 * We need to estimate the size of (F,G) during the execution of
	 * the algorithm; we are allowed some overestimation but not too
	 * much (poly_big_to_fp() uses a 310-bit window). Generally
	 * speaking, after applying a reduction with k scaled to
	 * scale_k, the size of (F,G) will be size(f,g) + scale_k + dd,
	 * where 'dd' is a few bits to account for the fact that the
	 * reduction is never perfect (intuitively, dd is on the order
	 * of sqrt(N), so at most 5 bits; we here allow for 10 extra
	 * bits).
	 *
	 * The size of (f,g) is not known exactly, but maxbl_fg is an
	 * upper bound.
	 */
	scale_k = maxbl_FG - minbl_fg;

	for (;;) {
		int scale_FG, dc, new_maxbl_FG;
		uint32_t scl, sch;
		fpr pdc, pt;

		/*
		 * Convert current F and G into floating-point. We apply
		 * scaling if the current length is more than 10 words.
		 */
		rlen = (FGlen > 10) ? 10 : FGlen;
		poly_big_to_fp(rt1, Ft + FGlen - rlen, rlen, llen, logn);
		poly_big_to_fp(rt2, Gt + FGlen - rlen, rlen, llen, logn);

		/*
		 * Values in rt1 and rt2 are downscaled by 2^(scale_FG).
		 */
		scale_FG = 31 * (int)(FGlen - rlen);

		/*
		 * Compute (F*adj(f)+G*adj(g))/(f*adj(f)+g*adj(g)) in rt2.
		 */
		Zf(FFT)(rt1, logn);
		Zf(FFT)(rt2, logn);
		Zf(poly_mul_fft)(rt1, rt3, logn);
		Zf(poly_mul_fft)(rt2, rt4, logn);
		Zf(poly_add)(rt2, rt1, logn);
		Zf(poly_mul_autoadj_fft)(rt2, rt5, logn);
		Zf(iFFT)(rt2, logn);

		/*
		 * (f,g) are scaled by 'scale_fg', meaning that the
		 * numbers in rt3/rt4 should be multiplied by 2^(scale_fg)
		 * to have their true mathematical value.
		 *
		 * (F,G) are similarly scaled by 'scale_FG'. Therefore,
		 * the value we computed in rt2 is scaled by
		 * 'scale_FG-scale_fg'.
		 *
		 * We want that value to be scaled by 'scale_k', hence we
		 * apply a corrective scaling. After scaling, the values
		 * should fit in -2^31-1..+2^31-1.
		 */
		dc = scale_k - scale_FG + scale_fg;

		/*
		 * We will need to multiply values by 2^(-dc). The value
		 * 'dc' is not secret, so we can compute 2^(-dc) with a
		 * non-constant-time process.
		 * (We could use ldexp(), but we prefer to avoid any
		 * dependency on libm. When using FP emulation, we could
		 * use our fpr_ldexp(), which is constant-time.)
		 */
		if (dc < 0) {
			dc = -dc;
			pt = fpr_two;
		} else {
			pt = fpr_onehalf;
		}
		pdc = fpr_one;
		while (dc != 0) {
			if ((dc & 1) != 0) {
				pdc = fpr_mul(pdc, pt);
			}
			dc >>= 1;
			pt = fpr_sqr(pt);
		}

		for (u = 0; u < n; u ++) {
			fpr xv;

			xv = fpr_mul(rt2[u], pdc);

			/*
			 * Sometimes the values can be out-of-bounds if
			 * the algorithm fails; we must not call
			 * fpr_rint() (and cast to int32_t) if the value
			 * is not in-bounds. Note that the test does not
			 * break constant-time discipline, since any
			 * failure here implies that we discard the current
			 * secret key (f,g).
			 */
			if (!fpr_lt(fpr_mtwo31m1, xv)
				|| !fpr_lt(xv, fpr_ptwo31m1))
			{
				return 0;
			}
			k[u] = (int32_t)fpr_rint(xv);
		}

		/*
		 * Values in k[] are integers. They really are scaled
		 * down by maxbl_FG - minbl_fg bits.
		 *
		 * If we are at low depth, then we use the NTT to
		 * compute k*f and k*g.
		 */
		sch = (uint32_t)(scale_k / 31);
		scl = (uint32_t)(scale_k % 31);
		if (depth <= DEPTH_INT_FG) {
			poly_sub_scaled_ntt(Ft, FGlen, llen, ft, slen, slen, k, sch, scl, logn, t1);
			poly_sub_scaled_ntt(Gt, FGlen, llen, gt, slen, slen, k, sch, scl, logn, t1);
		} else {
			poly_sub_scaled(Ft, FGlen, llen, ft, slen, slen, k, sch, scl, logn);
			poly_sub_scaled(Gt, FGlen, llen, gt, slen, slen, k, sch, scl, logn);
		}

		/*
		 * We compute the new maximum size of (F,G), assuming that
		 * (f,g) has _maximal_ length (i.e. that reduction is
		 * "late" instead of "early". We also adjust FGlen
		 * accordingly.
		 */
		new_maxbl_FG = scale_k + maxbl_fg + 10;
		if (new_maxbl_FG < maxbl_FG) {
			maxbl_FG = new_maxbl_FG;
			if ((int)FGlen * 31 >= maxbl_FG + 31) {
				FGlen --;
			}
		}

		/*
		 * We suppose that scaling down achieves a reduction by
		 * at least 25 bits per iteration. We stop when we have
		 * done the loop with an unscaled k.
		 */
		if (scale_k <= 0) {
			break;
		}
		scale_k -= 25;
		if (scale_k < 0) {
			scale_k = 0;
		}
	}

	/*
	 * If (F,G) length was lowered below 'slen', then we must take
	 * care to re-extend the sign.
	 */
	if (FGlen < slen) {
		for (u = 0; u < n; u ++, Ft += llen, Gt += llen) {
			size_t v;
			uint32_t sw;

			sw = -(Ft[FGlen - 1] >> 30) >> 1;
			for (v = FGlen; v < slen; v ++) {
				Ft[v] = sw;
			}
			sw = -(Gt[FGlen - 1] >> 30) >> 1;
			for (v = FGlen; v < slen; v ++) {
				Gt[v] = sw;
			}
		}
	}

	/*
	 * Compress encoding of all values to 'slen' words (this is the
	 * expected output format).
	 */
	for (u = 0, x = tmp, y = tmp;
		u < (n << 1); u ++, x += slen, y += llen)
	{
		memmove(x, y, slen * sizeof *y);
	}
	return 1;
}

/*
 * Solving the NTRU equation, binary case, depth = 1. Upon entry, the
 * F and G from the previous level should be in the tmp[] array.
 *
 * Returned value: 1 on success, 0 on error.
 */
static int
solve_NTRU_binary_depth1(unsigned logn_top,
	const int8_t *f, const int8_t *g, uint32_t *tmp)
{
	/*
	 * The first half of this function is a copy of the corresponding
	 * part in solve_NTRU_intermediate(), for the reconstruction of
	 * the unreduced F and G. The second half (Babai reduction) is
	 * done differently, because the unreduced F and G fit in 53 bits
	 * of precision, allowing a much simpler process with lower RAM
	 * usage.
	 */
	unsigned depth, logn;
	size_t n_top, n, hn, slen, dlen, llen, u;
	uint32_t *Fd, *Gd, *Ft, *Gt, *ft, *gt, *t1;
	fpr *rt1, *rt2, *rt3, *rt4, *rt5, *rt6;
	uint32_t *x, *y;

	depth = 1;
	n_top = MKN(logn_top);
	logn = logn_top - depth;
	n = MKN(logn);
	hn = n >> 1;

	/*
	 * Equations are:
	 *
	 *   f' = f0^2 - X^2*f1^2
	 *   g' = g0^2 - X^2*g1^2
	 *   F' and G' are a solution to f'G' - g'F' = 1 (from deeper levels)
	 *   F = F'*(g0 - X*g1)
	 *   G = G'*(f0 - X*f1)
	 *
	 * f0, f1, g0, g1, f', g', F' and G' are all "compressed" to
	 * degree N/2 (their odd-indexed coefficients are all zero).
	 */

	/*
	 * slen = size for our input f and g; also size of the reduced
	 *        F and G we return (degree N)
	 *
	 * dlen = size of the F and G obtained from the deeper level
	 *        (degree N/2)
	 *
	 * llen = size for intermediary F and G before reduction (degree N)
	 *
	 * We build our non-reduced F and G as two independent halves each,
	 * of degree N/2 (F = F0 + X*F1, G = G0 + X*G1).
	 */
	if (logn_top == 10) {
		slen = MAX_BL_SMALL_1024[depth];
		dlen = MAX_BL_SMALL_1024[depth + 1];
		llen = MAX_BL_LARGE_1024[depth];
	} else {
		slen = MAX_BL_SMALL_512[depth];
		dlen = MAX_BL_SMALL_512[depth + 1];
		llen = MAX_BL_LARGE_512[depth];
	}

	/*
	 * Fd and Gd are the F and G from the deeper level. Ft and Gt
	 * are the destination arrays for the unreduced F and G.
	 */
	Fd = tmp;
	Gd = Fd + dlen * hn;
	Ft = Gd + dlen * hn;
	Gt = Ft + llen * n;

	/*
	 * We reduce Fd and Gd modulo all the small primes we will need,
	 * and store the values in Ft and Gt.
	 */
	for (u = 0; u < llen; u ++) {
		uint32_t p, p0i, R2, Rx;
		size_t v;
		uint32_t *xs, *ys, *xd, *yd;

		p = PRIMES[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);
		Rx = modp_Rx((unsigned)dlen, p, p0i, R2);
		for (v = 0, xs = Fd, ys = Gd, xd = Ft + u, yd = Gt + u;
			v < hn;
			v ++, xs += dlen, ys += dlen, xd += llen, yd += llen)
		{
			*xd = zint_mod_small_signed(xs, dlen, p, p0i, R2, Rx);
			*yd = zint_mod_small_signed(ys, dlen, p, p0i, R2, Rx);
		}
	}

	/*
	 * Now Fd and Gd are not needed anymore; we can squeeze them out.
	 */
	memmove(tmp, Ft, llen * n * sizeof(uint32_t));
	Ft = tmp;
	memmove(Ft + llen * n, Gt, llen * n * sizeof(uint32_t));
	Gt = Ft + llen * n;
	ft = Gt + llen * n;
	gt = ft + slen * n;

	t1 = gt + slen * n;

	/*
	 * Compute our F and G modulo sufficiently many small primes.
	 */
	for (u = 0; u < llen; u ++) {
		uint32_t p, p0i, R2;
		uint32_t *gm, *igm, *fx, *gx, *Fp, *Gp;
		unsigned e;
		size_t v;

		/*
		 * All computations are done modulo p.
		 */
		p = PRIMES[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);

		/*
		 * We recompute things from the source f and g, of full
		 * degree. However, we will need only the n first elements
		 * of the inverse NTT table (igm); the call to modp_mkgm()
		 * below will fill n_top elements in igm[] (thus overflowing
		 * into fx[]) but later code will overwrite these extra
		 * elements.
		 */
		gm = t1;
		igm = gm + n_top;
		fx = igm + n;
		gx = fx + n_top;
		modp_mkgm2(gm, igm, logn_top, PRIMES[u].g, p, p0i);

		/*
		 * Set ft and gt to f and g modulo p, respectively.
		 */
		for (v = 0; v < n_top; v ++) {
			fx[v] = modp_set(f[v], p);
			gx[v] = modp_set(g[v], p);
		}

		/*
		 * Convert to NTT and compute our f and g.
		 */
		modp_NTT2(fx, gm, logn_top, p, p0i);
		modp_NTT2(gx, gm, logn_top, p, p0i);
		for (e = logn_top; e > logn; e --) {
			modp_poly_rec_res(fx, e, p, p0i, R2);
			modp_poly_rec_res(gx, e, p, p0i, R2);
		}

		/*
		 * From that point onward, we only need tables for
		 * degree n, so we can save some space.
		 */
		if (depth > 0) { /* always true */
			memmove(gm + n, igm, n * sizeof *igm);
			igm = gm + n;
			memmove(igm + n, fx, n * sizeof *ft);
			fx = igm + n;
			memmove(fx + n, gx, n * sizeof *gt);
			gx = fx + n;
		}

		/*
		 * Get F' and G' modulo p and in NTT representation
		 * (they have degree n/2). These values were computed
		 * in a previous step, and stored in Ft and Gt.
		 */
		Fp = gx + n;
		Gp = Fp + hn;
		for (v = 0, x = Ft + u, y = Gt + u;
			v < hn; v ++, x += llen, y += llen)
		{
			Fp[v] = *x;
			Gp[v] = *y;
		}
		modp_NTT2(Fp, gm, logn - 1, p, p0i);
		modp_NTT2(Gp, gm, logn - 1, p, p0i);

		/*
		 * Compute our F and G modulo p.
		 *
		 * Equations are:
		 *
		 *   f'(x^2) = N(f)(x^2) = f * adj(f)
		 *   g'(x^2) = N(g)(x^2) = g * adj(g)
		 *
		 *   f'*G' - g'*F' = 1
		 *
		 *   F = F'(x^2) * adj(g)
		 *   G = G'(x^2) * adj(f)
		 *
		 * The NTT representation of f is f(w) for all w which
		 * are roots of phi. In the binary case, as well as in
		 * the ternary case for all depth except the deepest,
		 * these roots can be grouped in pairs (w,-w), and we
		 * then have:
		 *
		 *   f(w) = adj(f)(-w)
		 *   f(-w) = adj(f)(w)
		 *
		 * and w^2 is then a root for phi at the half-degree.
		 *
		 * At the deepest level in the ternary case, this still
		 * holds, in the following sense: the roots of x^2-x+1
		 * are (w,-w^2) (for w^3 = -1, and w != -1), and we
		 * have:
		 *
		 *   f(w) = adj(f)(-w^2)
		 *   f(-w^2) = adj(f)(w)
		 *
		 * In all case, we can thus compute F and G in NTT
		 * representation by a few simple multiplications.
		 * Moreover, the two roots for each pair are consecutive
		 * in our bit-reversal encoding.
		 */
		for (v = 0, x = Ft + u, y = Gt + u;
			v < hn; v ++, x += (llen << 1), y += (llen << 1))
		{
			uint32_t ftA, ftB, gtA, gtB;
			uint32_t mFp, mGp;

			ftA = fx[(v << 1) + 0];
			ftB = fx[(v << 1) + 1];
			gtA = gx[(v << 1) + 0];
			gtB = gx[(v << 1) + 1];
			mFp = modp_montymul(Fp[v], R2, p, p0i);
			mGp = modp_montymul(Gp[v], R2, p, p0i);
			x[0] = modp_montymul(gtB, mFp, p, p0i);
			x[llen] = modp_montymul(gtA, mFp, p, p0i);
			y[0] = modp_montymul(ftB, mGp, p, p0i);
			y[llen] = modp_montymul(ftA, mGp, p, p0i);
		}
		modp_iNTT2_ext(Ft + u, llen, igm, logn, p, p0i);
		modp_iNTT2_ext(Gt + u, llen, igm, logn, p, p0i);

		/*
		 * Also save ft and gt (only up to size slen).
		 */
		if (u < slen) {
			modp_iNTT2(fx, igm, logn, p, p0i);
			modp_iNTT2(gx, igm, logn, p, p0i);
			for (v = 0, x = ft + u, y = gt + u;
				v < n; v ++, x += slen, y += slen)
			{
				*x = fx[v];
				*y = gx[v];
			}
		}
	}

	/*
	 * Rebuild f, g, F and G with the CRT. Note that the elements of F
	 * and G are consecutive, and thus can be rebuilt in a single
	 * loop; similarly, the elements of f and g are consecutive.
	 */
	zint_rebuild_CRT(Ft, llen, llen, n << 1, 1, t1);
	zint_rebuild_CRT(ft, slen, slen, n << 1, 1, t1);

	/*
	 * Here starts the Babai reduction, specialized for depth = 1.
	 *
	 * Candidates F and G (from Ft and Gt), and base f and g (ft and gt),
	 * are converted to floating point. There is no scaling, and a
	 * single pass is sufficient.
	 */

	/*
	 * Convert F and G into floating point (rt1 and rt2).
	 */
	rt1 = align_fpr(tmp, gt + slen * n);
	rt2 = rt1 + n;
	poly_big_to_fp(rt1, Ft, llen, llen, logn);
	poly_big_to_fp(rt2, Gt, llen, llen, logn);

	/*
	 * Integer representation of F and G is no longer needed, we
	 * can remove it.
	 */
	memmove(tmp, ft, 2 * slen * n * sizeof *ft);
	ft = tmp;
	gt = ft + slen * n;
	rt3 = align_fpr(tmp, gt + slen * n);
	memmove(rt3, rt1, 2 * n * sizeof *rt1);
	rt1 = rt3;
	rt2 = rt1 + n;
	rt3 = rt2 + n;
	rt4 = rt3 + n;

	/*
	 * Convert f and g into floating point (rt3 and rt4).
	 */
	poly_big_to_fp(rt3, ft, slen, slen, logn);
	poly_big_to_fp(rt4, gt, slen, slen, logn);

	/*
	 * Remove unneeded ft and gt.
	 */
	memmove(tmp, rt1, 4 * n * sizeof *rt1);
	rt1 = (fpr *)tmp;
	rt2 = rt1 + n;
	rt3 = rt2 + n;
	rt4 = rt3 + n;

	/*
	 * We now have:
	 *   rt1 = F
	 *   rt2 = G
	 *   rt3 = f
	 *   rt4 = g
	 * in that order in RAM. We convert all of them to FFT.
	 */
	Zf(FFT)(rt1, logn);
	Zf(FFT)(rt2, logn);
	Zf(FFT)(rt3, logn);
	Zf(FFT)(rt4, logn);

	/*
	 * Compute:
	 *   rt5 = F*adj(f) + G*adj(g)
	 *   rt6 = 1 / (f*adj(f) + g*adj(g))
	 * (Note that rt6 is half-length.)
	 */
	rt5 = rt4 + n;
	rt6 = rt5 + n;
	Zf(poly_add_muladj_fft)(rt5, rt1, rt2, rt3, rt4, logn);
	Zf(poly_invnorm2_fft)(rt6, rt3, rt4, logn);

	/*
	 * Compute:
	 *   rt5 = (F*adj(f)+G*adj(g)) / (f*adj(f)+g*adj(g))
	 */
	Zf(poly_mul_autoadj_fft)(rt5, rt6, logn);

	/*
	 * Compute k as the rounded version of rt5. Check that none of
	 * the values is larger than 2^63-1 (in absolute value)
	 * because that would make the fpr_rint() do something undefined;
	 * note that any out-of-bounds value here implies a failure and
	 * (f,g) will be discarded, so we can make a simple test.
	 */
	Zf(iFFT)(rt5, logn);
	for (u = 0; u < n; u ++) {
		fpr z;

		z = rt5[u];
		if (!fpr_lt(z, fpr_ptwo63m1) || !fpr_lt(fpr_mtwo63m1, z)) {
			return 0;
		}
		rt5[u] = fpr_of(fpr_rint(z));
	}
	Zf(FFT)(rt5, logn);

	/*
	 * Subtract k*f from F, and k*g from G.
	 */
	Zf(poly_mul_fft)(rt3, rt5, logn);
	Zf(poly_mul_fft)(rt4, rt5, logn);
	Zf(poly_sub)(rt1, rt3, logn);
	Zf(poly_sub)(rt2, rt4, logn);
	Zf(iFFT)(rt1, logn);
	Zf(iFFT)(rt2, logn);

	/*
	 * Convert back F and G to integers, and return.
	 */
	Ft = tmp;
	Gt = Ft + n;
	rt3 = align_fpr(tmp, Gt + n);
	memmove(rt3, rt1, 2 * n * sizeof *rt1);
	rt1 = rt3;
	rt2 = rt1 + n;
	for (u = 0; u < n; u ++) {
		Ft[u] = (uint32_t)fpr_rint(rt1[u]);
		Gt[u] = (uint32_t)fpr_rint(rt2[u]);
	}

	return 1;
}

/*
 * Solving the NTRU equation, top level. Upon entry, the F and G from the
 * previous level should be in the tmp[] array.
 *
 * Returned value: 1 on success, 0 on error.
 */
static int
solve_NTRU_binary_depth0(unsigned logn,
	const int8_t *f, const int8_t *g, uint32_t *tmp)
{
	size_t n, hn, u;
	uint32_t p, p0i, R2;
	uint32_t *Fp, *Gp, *t1, *t2, *t3, *t4, *t5;
	uint32_t *gm, *igm, *ft, *gt;
	fpr *rt1, *rt2, *rt3;

	n = MKN(logn);
	hn = n >> 1;

	/*
	 * Equations are:
	 *
	 *   f' = f0^2 - X^2*f1^2
	 *   g' = g0^2 - X^2*g1^2
	 *   F' and G' are a solution to f'G' - g'F' = 1 (from deeper levels)
	 *   F = F'*(g0 - X*g1)
	 *   G = G'*(f0 - X*f1)
	 *
	 * f0, f1, g0, g1, f', g', F' and G' are all "compressed" to
	 * degree N/2 (their odd-indexed coefficients are all zero).
	 *
	 * Everything should fit in 31-bit integers, hence we can just use
	 * the first small prime p = 2147473409.
	 */
	p = PRIMES[0].p;
	p0i = modp_ninv31(p);
	R2 = modp_R2(p, p0i);

	Fp = tmp;
	Gp = Fp + hn;
	ft = Gp + hn;
	gt = ft + n;
	gm = gt + n;
	igm = gm + n;

	modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);

	/*
	 * Convert F' anf G' in NTT representation.
	 */
	for (u = 0; u < hn; u ++) {
		Fp[u] = modp_set(zint_one_to_plain(Fp + u), p);
		Gp[u] = modp_set(zint_one_to_plain(Gp + u), p);
	}
	modp_NTT2(Fp, gm, logn - 1, p, p0i);
	modp_NTT2(Gp, gm, logn - 1, p, p0i);

	/*
	 * Load f and g and convert them to NTT representation.
	 */
	for (u = 0; u < n; u ++) {
		ft[u] = modp_set(f[u], p);
		gt[u] = modp_set(g[u], p);
	}
	modp_NTT2(ft, gm, logn, p, p0i);
	modp_NTT2(gt, gm, logn, p, p0i);

	/*
	 * Build the unreduced F,G in ft and gt.
	 */
	for (u = 0; u < n; u += 2) {
		uint32_t ftA, ftB, gtA, gtB;
		uint32_t mFp, mGp;

		ftA = ft[u + 0];
		ftB = ft[u + 1];
		gtA = gt[u + 0];
		gtB = gt[u + 1];
		mFp = modp_montymul(Fp[u >> 1], R2, p, p0i);
		mGp = modp_montymul(Gp[u >> 1], R2, p, p0i);
		ft[u + 0] = modp_montymul(gtB, mFp, p, p0i);
		ft[u + 1] = modp_montymul(gtA, mFp, p, p0i);
		gt[u + 0] = modp_montymul(ftB, mGp, p, p0i);
		gt[u + 1] = modp_montymul(ftA, mGp, p, p0i);
	}
	modp_iNTT2(ft, igm, logn, p, p0i);
	modp_iNTT2(gt, igm, logn, p, p0i);

	Gp = Fp + n;
	t1 = Gp + n;
	memmove(Fp, ft, 2 * n * sizeof *ft);

	/*
	 * We now need to apply the Babai reduction. At that point,
	 * we have F and G in two n-word arrays.
	 *
	 * We can compute F*adj(f)+G*adj(g) and f*adj(f)+g*adj(g)
	 * modulo p, using the NTT. We still move memory around in
	 * order to save RAM.
	 */
	t2 = t1 + n;
	t3 = t2 + n;
	t4 = t3 + n;
	t5 = t4 + n;

	/*
	 * Compute the NTT tables in t1 and t5.
	 */
	memmove(t1, gm, n * sizeof *gm);
	memmove(t5, igm, n * sizeof *igm);
	// modp_mkgm2(t1, t5, logn, PRIMES[0].g, p, p0i);

	/*
	 * Convert F and G to NTT.
	 */
	modp_NTT2(Fp, t1, logn, p, p0i);
	modp_NTT2(Gp, t1, logn, p, p0i);

	/*
	 * Load f in t4 and convert it to NTT representation.
	 */
	for (u = 0; u < n; u ++) {
		t4[u] = modp_set(f[u], p);
	}
	modp_NTT2(t4, t1, logn, p, p0i);

	/*
	 * Compute F*adj(f) in t2, and f*adj(f) in t3.
	 */
	for (u = 0; u < n; u ++) {
		uint32_t w;

		w = modp_montymul(t4[n - 1 - u], R2, p, p0i);
		t2[u] = modp_montymul(w, Fp[u], p, p0i);
		t3[u] = modp_montymul(w, t4[u], p, p0i);
	}

	/*
	 * Load g in t4, and convert it to NTT representation.
	 */
	for (u = 0; u < n; u ++) {
		t4[u] = modp_set(g[u], p);
	}
	modp_NTT2(t4, t1, logn, p, p0i);

	/*
	 * Add G*adj(g) to t2, and g*adj(g) to t3.
	 */
	for (u = 0; u < n; u ++) {
		uint32_t w;

		w = modp_montymul(t4[n - 1 - u], R2, p, p0i);
		t2[u] = modp_add(t2[u],
			modp_montymul(w, Gp[u], p, p0i), p);
		t3[u] = modp_add(t3[u],
			modp_montymul(w, t4[u], p, p0i), p);
	}

	/*
	 * Convert back t2 and t3 to normal representation (normalized
	 * around 0), and then move them to t1 and t2.
	 */
	modp_iNTT2(t2, t5, logn, p, p0i);
	modp_iNTT2(t3, t5, logn, p, p0i);
	for (u = 0; u < n; u ++) {
		t1[u] = (uint32_t)modp_norm(t2[u], p);
		t2[u] = (uint32_t)modp_norm(t3[u], p);
	}

	/*
	 * At that point, array contents are:
	 *
	 *   F (NTT representation) (Fp)
	 *   G (NTT representation) (Gp)
	 *   F*adj(f)+G*adj(g) (t1)
	 *   f*adj(f)+g*adj(g) (t2)
	 *
	 * We want to divide t1 by t2. The result is not integral; it
	 * must be rounded. We thus need to use the FFT.
	 */

	/*
	 * Get F*adj(f)+G*adj(g) and f*adj(f)+g*adj(g) in FFT representation.
	 */
	rt1 = align_fpr(tmp, t1);
	rt2 = rt1 + n;
	rt3 = rt2 + n;
	for (u = 0; u < n; u ++) {
		rt3[u] = fpr_of(((int32_t *)t2)[u]);
		rt2[u] = fpr_of(((int32_t *)t1)[u]);
	}
	Zf(FFT)(rt2, logn);
	Zf(FFT)(rt3, logn);
	memmove(rt1, rt3, n * sizeof *rt3);

	/*
	 * Compute (F*adj(f)+G*adj(g))/(f*adj(f)+g*adj(g)) and get
	 * its rounded normal representation in t1.
	 */
	Zf(poly_div_autoadj_fft)(rt2, rt1, logn);

	/*
	 * Babai (F, G) with respect to the GSO of (f, g) using Ducas--Prest's Fast
	 * Fourier Orthogonalization ffNP routine.
	 * Comment out if you don't want to perform ffNP.
	 */
	Zf(ffNearestPlane_dyn)(rt2, rt1, logn, rt3);

	Zf(iFFT)(rt2, logn);
	for (u = 0; u < n; u ++) {
		t1[u] = modp_set((int32_t)fpr_rint(rt2[u]), p);
	}

	/*
	 * RAM contents are now:
	 *
	 *   F (NTT representation) (Fp)
	 *   G (NTT representation) (Gp)
	 *   k (t1)
	 *
	 * We want to compute F-k*f, and G-k*g.
	 */
	modp_mkgm2(t2, t3, logn, PRIMES[0].g, p, p0i);
	for (u = 0; u < n; u ++) {
		t4[u] = modp_set(f[u], p);
		t5[u] = modp_set(g[u], p);
	}
	modp_NTT2(t1, t2, logn, p, p0i);
	modp_NTT2(t4, t2, logn, p, p0i);
	modp_NTT2(t5, t2, logn, p, p0i);
	for (u = 0; u < n; u ++) {
		uint32_t kw;

		kw = modp_montymul(t1[u], R2, p, p0i);
		Fp[u] = modp_sub(Fp[u],
			modp_montymul(kw, t4[u], p, p0i), p);
		Gp[u] = modp_sub(Gp[u],
			modp_montymul(kw, t5[u], p, p0i), p);
	}
	modp_iNTT2(Fp, t3, logn, p, p0i);
	modp_iNTT2(Gp, t3, logn, p, p0i);
	for (u = 0; u < n; u ++) {
		Fp[u] = (uint32_t)modp_norm(Fp[u], p);
		Gp[u] = (uint32_t)modp_norm(Gp[u], p);
	}

	return 1;
}


/* ==================================================================== */

/*
 * Solve the NTRU equation, but now for q = 1. Returned value is 1 on success,
 * 0 on error. G can be NULL, in which case that value is computed but not
 * returned. If any of the coefficients of F and G exceeds lim (in absolute
 * value), then 0 is returned.
 */
static int
solve_NTRU(unsigned logn, int8_t *F, int8_t *G,
	const int8_t *f, const int8_t *g, int lim, uint32_t *tmp)
{
	size_t n, u, depth;
	uint32_t *ft, *gt, *Ft, *Gt, *gm;
	uint32_t p, p0i, r, z;

	n = MKN(logn);

	if (!solve_NTRU_deepest(logn, f, g, tmp)) {
		return 0;
	}

	for (depth = logn - 1; depth >= 2; depth --) {
		if (!solve_NTRU_intermediate(logn, f, g, depth, tmp)) {
			return 0;
		}
	}

	if (logn >= 2 && !solve_NTRU_binary_depth1(logn, f, g, tmp)) {
		return 0;
	}

	/*
	 * Solve the top level of the NTRU equation, assuming the F, G from
	 * previous level are in tmp. Note that this method also does FFO Babai
	 * reduction on (F, G).
	 */
	if (!solve_NTRU_binary_depth0(logn, f, g, tmp)) {
		return 0;
	}

	/*
	 * Final F and G are in fk->tmp, one word per coefficient
	 * (signed value over 31 bits).
	 */
	if (!poly_big_to_small(F, tmp, lim, logn)
		|| !poly_big_to_small(G, tmp + n, lim, logn))
	{
		return 0;
	}

	/*
	 * Verify that the NTRU equation fG - gF = 1 is fulfilled. Since all
	 * elements have short lengths, verifying modulo a small prime p works and
	 * allows using the NTT.
	 */
	ft = tmp;
	gt = ft + n;
	Ft = gt + n;
	Gt = Ft + n;
	gm = Gt + n;

	p = PRIMES[0].p;
	p0i = modp_ninv31(p);
	modp_mkgm2(gm, tmp, logn, PRIMES[0].g, p, p0i);
	for (u = 0; u < n; u ++) {
		ft[u] = modp_set(f[u], p);
		gt[u] = modp_set(g[u], p);
		Ft[u] = modp_set(F[u], p);
		Gt[u] = modp_set(G[u], p);
	}
	modp_NTT2(ft, gm, logn, p, p0i);
	modp_NTT2(gt, gm, logn, p, p0i);
	modp_NTT2(Ft, gm, logn, p, p0i);
	modp_NTT2(Gt, gm, logn, p, p0i);

	r = modp_montymul(1, 1, p, p0i);
	for (u = 0; u < n; u ++) {
		z = modp_sub(modp_montymul(ft[u], Gt[u], p, p0i),
			modp_montymul(gt[u], Ft[u], p, p0i), p);
		if (z != r) {
			return 0;
		}
	}

	return 1;
}

/* ========================================================================= */

/* see inner.h */
void
Zf(make_public)(const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	int16_t *restrict iq00, int16_t *restrict iq01, int16_t *restrict iq11,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u, v;
	uint32_t p, p0i, R2, *gm, *igm, *t0, *t1, *t2, *t3;

	n = MKN(logn);
	p = PRIMES[0].p;
	p0i = modp_ninv31(p);
	R2 = modp_R2(p, p0i);

	gm = (uint32_t *)tmp;
	igm = gm + n;
	t0 = igm + n;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;

	modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);

	for (u = 0; u < n; u++) {
		t0[u] = modp_set((int32_t)f[u], p);
		t1[u] = modp_set((int32_t)g[u], p);
		t2[u] = modp_set((int32_t)F[u], p);
	}
	modp_NTT2(t0, gm, logn, p, p0i);
	modp_NTT2(t1, gm, logn, p, p0i);
	modp_NTT2(t2, gm, logn, p, p0i);

	if (G == NULL) {
		/*
		 * Compute (in NTT representation) G = (1 + gF) / f.
		 */
		for (u = 0; u < n; u++) {
			t3[u] = modp_montymul(t1[u], R2, p, p0i);
			t3[u] = modp_add(1U, modp_montymul(t3[u], t2[u], p, p0i), p);
			t3[u] = modp_div(t3[u], t0[u], p, p0i);
		}
	} else {
		for (u = 0; u < n; u++) {
			t3[u] = modp_set((int32_t)G[u], p);
		}
		modp_NTT2(t3, gm, logn, p, p0i);
	}

	for (u = 0, v = n - 1; u < v; u++, v--) {
		uint32_t q00, q01a, q01b;

		/*
		 * Calculate q00 = adj(f) * f + adj(g) * g.
		 */
		q00 = modp_add(
			modp_montymul(t0[u], t0[v], p, p0i),
			modp_montymul(t1[u], t1[v], p, p0i), p);
		q00 = modp_montymul(q00, R2, p, p0i);

		/*
		 * Calculate q01 = adj(f) * F + adj(g) * G.
		 */
		q01a = modp_add(
			modp_montymul(t0[v], t2[u], p, p0i),
			modp_montymul(t1[v], t3[u], p, p0i), p);
		q01b = modp_add(
			modp_montymul(t0[u], t2[v], p, p0i),
			modp_montymul(t1[u], t3[v], p, p0i), p);

		/*
		 * Now store the result in t0, t1.
		 */
		t0[u] = t0[v] = q00;
		t1[u] = modp_montymul(q01a, R2, p, p0i);
		t1[v] = modp_montymul(q01b, R2, p, p0i);

		if (iq11 != NULL) {
			/*
			 * Calculate q11 = adj(f) * f + adj(g) * g and store result in t2.
			 */
			t2[u] = modp_add(
				modp_montymul(t2[u], t2[v], p, p0i),
				modp_montymul(t3[u], t3[v], p, p0i), p);
			t2[u] = t2[v] = modp_montymul(t2[u], R2, p, p0i);
		}
	}

	modp_iNTT2(t0, igm, logn, p, p0i);
	modp_iNTT2(t1, igm, logn, p, p0i);

	for (u = 0; u < n; u++) {
		iq00[u] = modp_norm(t0[u], p);
		iq01[u] = modp_norm(t1[u], p);
	}

	if (iq11 != NULL) {
		modp_iNTT2(t2, igm, logn, p, p0i);
		for (u = 0; u < n; u++) {
			iq11[u] = modp_norm(t2[u], p);
		}
	}
}

static const int32_t l2bound_ssec_512[10] = {
	0u /* unused */, 8, 16, 32, 64, 129, 259, 519, 1039, 2079
};
static const int32_t l2bound_ssec_1024[11] = {
	0u /* unused */, 15, 31, 62, 124, 249, 498, 997, 1995, 3990, 7980
};


/* see inner.h */
void
Zf(keygen)(inner_shake256_context *sc,
	int8_t *restrict f, int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G,
	int16_t *restrict iq00, int16_t *restrict iq01,
	unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * Algorithm is the following:
	 *
	 *  - Generate f and g with the Gaussian distribution.
	 *
	 *  - If either N(f) or N(g) is even, try again.
	 *
	 *  - Solve the NTRU equation fG - gF = 1; if the solving fails, try again.
	 *    Usual failure condition is when N(f) and N(g) are not coprime.
	 *
	 *  - Use Babai Reduction on F, G.
	 *
	 *  - Calculate the Gram matrix of the basis [[f, g], [F, G]].
	 */
	size_t n, hn, u;
	uint8_t fg_okay;
	int16_t *iq11, *atmp;
	int32_t norm, bound, x;
	uint32_t p, p0i, *gm, *igm, *rf1, *rf2;
	prng rng;
	fpr *rt1, *rt2, *rt3;

	n = MKN(logn);
	hn = n >> 1;

	iq11 = (int16_t *)tmp;
	atmp = iq11 + n;

	rt1 = (fpr *)tmp;
	rt2 = rt1 + n;
	rt3 = rt2 + n;

	gm = (uint32_t *)tmp;
	igm = gm + n;
	rf1 = igm + n;
	rf2 = rf1 + n;

	p = PRIMES[0].p;
	p0i = modp_ninv31(p);

	for (;;) {
		/*
		 * The coefficients of polynomials f and g will be generated from a
		 * discrete gaussian that draws random numbers from a fast PRNG that is
		 * seeded from a SHAKE context ('sc').
		 */
		Zf(prng_init)(&rng, sc);

		/*
		 * The coefficients of f and g are generated independently of each
		 * other, with a discrete Gaussian distribution of standard deviation
		 * 1.500. The expected l2-norm of (f, g) is 2n sigma^2.
		 *
		 * We require that N(f) and N(g) are both odd (the binary GCD in the
		 * NTRU solver requires it), so we require (fg_okay & 1) == 1.
		 */
		fg_okay = poly_small_mkgauss(&rng, f, logn)
			& poly_small_mkgauss(&rng, g, logn) & 1u;

		fg_okay &= Zf(mf_is_invertible)(f, logn, tmp);

		Zf(int8_to_fft)(rt2, f, logn);
		Zf(int8_to_fft)(rt3, g, logn);
		Zf(poly_invnorm2_fft)(rt1, rt2, rt3, logn);
		Zf(iFFT)(rt1, logn);

		if (logn == 9) {
			/*
			 * For n = 512, we reject a key pair if cst(1/q00) >= 1/1000, as
			 * the failure probability of decompressing a signature is bounded
			 * from above by 1.9e-32 < 2^{-105}. Experimentally this fails
			 * with a probability of 9%.
			 */
			fg_okay &= fpr_lt(rt1[0], fpr_inv(fpr_of(1000)));
		} else if (logn == 10) {
			/*
			 * For n = 1024, we reject a key pair if cst(1/q00) >= 1/3000, as
			 * the failure probability of decompressing a signature is bounded
			 * from above by 1.2e-95 < 2^{-315}. Experimentally this fails
			 * with a probability of 0.9%.
			 */
			fg_okay &= fpr_lt(rt1[0], fpr_inv(fpr_of(3000)));
		}

		/*
		 * If NTT_p(q00) has a zero, we cannot invert it in
		 * Zf(uncompressed_verify_NTT), so remove this key.
		 */
		modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);
		for (u = 0; u < n; u++) {
			rf1[u] = modp_set((int32_t)f[u], p);
			rf2[u] = modp_set((int32_t)g[u], p);
		}
		modp_NTT2(rf1, gm, logn, p, p0i);
		modp_NTT2(rf2, gm, logn, p, p0i);
		for (u = 0; u < hn; u++) {
			uint32_t prod;

			prod = modp_add(modp_montymul(rf1[u], rf1[n - 1 - u], p, p0i),
							modp_montymul(rf2[u], rf2[n - 1 - u], p, p0i), p);
			fg_okay &= (1U - prod) >> 31;
		}

		/*
		 * If the l2-norm of (f, g) is shorter than sigma_ver^2 * 2n, BKZ may
		 * return a shortest vector when given the public key much faster than
		 * other instances, so this secret key is not secure to use.
		 * Thus, set fg_okay to 0 when ||(f, g)||^2 < l2bound_ssec[logn].
		 */
		norm = 0;
		for (u = 0; u < n; u++) {
			norm += (int32_t)f[u] * (int32_t)f[u];
			norm += (int32_t)g[u] * (int32_t)g[u];
		}

		norm -= (logn == 10 ? l2bound_ssec_1024[logn] : l2bound_ssec_512[logn]);
		fg_okay &= ((uint32_t) -norm) >> 31;

		if (fg_okay == 0) {
			/*
			 * Generation of (f, g) failed because:
			 * 1) N(f) or N(g) was even,
			 * 2) NTT(f) had zero coefficient,
			 * 3) NTT(q00) had zero coefficient,
			 * 4) cst(q00) = ||(f,g)||^2 < l2bound(logn)/4, or
			 * 5) cst(1/q00) was too large.
			 *
			 * Thus, resample f and g.
			 */
			continue;
		}

		/*
		 * Try to solve the NTRU equation for polynomials f and g, i.e. find
		 * polynomials F, G that satisfy
		 *
		 *     f * G - g * F = 1 (mod X^n + 1).
		 */
		if (!solve_NTRU(logn, F, G, f, g, 127, (uint32_t *)tmp)) {
			continue;
		}

		Zf(make_public)(f, g, F, G, iq00, iq01, iq11, logn, (uint8_t *)atmp);

		/*
		 * Check the bounds on q00 and q11.
		 */
		bound = (int32_t)(1U << Zf(bits_q00)[logn]);
		for (u = 1; u < hn; u++) {
			x = iq00[u];
			fg_okay &= (+x - bound) >> 31;
			fg_okay &= (-x - bound) >> 31;
			fg_okay &= x == -iq00[n - u];
		}

		bound = (int32_t)(1U << Zf(bits_q11)[logn]);
		for (u = 1; u < hn; u++) {
			x = iq11[u];
			fg_okay &= (+x - bound) >> 31;
			fg_okay &= (-x - bound) >> 31;
			fg_okay &= x == -iq11[n - u];
		}

		bound = (int32_t)(1U << Zf(bits_q01)[logn]);
		for (u = 0; u < n; u++) {
			x = iq01[u];
			fg_okay &= (+x - bound) >> 31;
			fg_okay &= (-x - bound) >> 31;
		}

		if (fg_okay == 0) {
			/*
			 * There was a coefficient that was too large, i.e. there
			 * was an x in the above with
			 *   x <= -bound or x >= bound,
			 * or q00 and q11 failed to be selfadjoint (should not happen).
			 */
			continue;
		}

		/*
		 * A valid key pair is generated.
		 */
		break;
	}
}
