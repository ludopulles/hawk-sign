/*
 * Hawk signature generation.
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
 * @author   Ludo Pulles <ludo.pulles@cwi.nl>
 */

#include "inner.h"

/*
 * Constants for NTT.
 *
 *   n = 2^logn  (2 <= n <= 1024)
 *   phi = X^n + 1
 *   q == 1 (mod 2048), q prime and < 2^16
 *   q0i = -1/q mod 2^16
 *   R = 2^16 mod q
 *   R2 = 2^32 mod q
 *
 * Note that the options for q < 2^16, == 1 (mod 2048) are:
 *   12289, 18433, 40961, 59393, 61441.
 */

#define Q   18433
#define Q0I 18431
#define R   10237
#define R2   4564

/*
 * Table for NTT, binary case:
 *   GMb[x] = R*(g^rev(x)) mod q
 * where g = 19 (it is a 2048-th primitive root of 1 modulo q)
 * and rev() is the bit-reversal function over 10 bits.
 */
static const uint16_t GMb[1024] = {
	10237, 17077,  6948,  4658,  9087,  7063,  9974,  2028,
	 5729,  2891, 14386, 16468,  7338,  1322,  6113,  1875,
	 3627, 16901, 16006, 16790,   479,  5261,  1753, 16483,
	  872,   765, 10272,  9688, 15631, 14326, 14682,   324,
	13141,   177,  9533,  6651,  1954, 12495,  1879,  4629,
	 3942,  5699, 15658,  3886, 13968, 18342,  4055,  5016,
	13596, 14718,  2886,  8493,   957, 17053,  3156, 14691,
	 7553, 16498,  7631,  4771,  9256,  9304, 13825, 12192,
	  193, 11394, 12251,  6372, 17686, 12345, 15045,  7428,
	 7568, 10698, 12881,  2401,  3754, 16949,  5157, 15157,
	 7613, 11731, 10198, 13724,  5681,  3018, 16019,  5619,
	12574, 16654, 13846,  4072,  3480,     9,  1422,  3150,
	 3849,  4793,  1541,   147, 14901,  7809, 17244,  5066,
	 1459,  1132, 12959,  9107,  6769, 12428,  9726, 18045,
	15441,  1772,  3481, 11911,  6232, 17305,  6106, 10726,
	  724,  8837, 13771, 14639, 18305,  6483, 10499,  1791,
	 7702,  1895,  4482, 18095, 17275,  5368,   226, 17067,
	17210,  5924, 14342,  8904,  9891,  9544, 14879,  4027,
	15281, 14484,  2780,   325,  9621,  3346, 12544,  9821,
	  467,  9901, 15986, 18379, 16721, 10674,  9089, 12434,
	11963,  7134,  2759,  8445, 13772,  8108,  9187, 17551,
	16105, 15376, 14685, 17597,  9679, 11641, 14411,   657,
	 9057,   230, 17907,  6768, 17952,  7801, 15980,  2266,
	10746, 10738,   768, 16401, 14089,  2277,  9539,  4331,
	 9341,  7159,  6709, 17195, 13319, 17371, 16534, 15393,
	11663, 12536,  8357,   546, 13214,  2672, 16650, 13550,
	14414, 17930, 12691,  8280, 10589,  3857,  1117,  4341,
	17853,  9215, 18196, 17908,  9981, 11610,  9513,  8240,
	 6733,  7915, 15559,  5300, 15104,  9192, 14562,  9858,
	 1944,  4073, 16812,  6209, 13201, 13843, 12100, 15604,
	11250,   188, 11271, 10501,  2492,  1087,  5849, 11790,
	12310,  8136, 13611,  8918,   777, 12921, 13888,  6265,
	 8957,  8175,  1340,  4135, 17766,  5989,  6179, 13221,
	 3673, 11403, 13673,  9522, 13099, 16417, 13266, 13287,
	  604, 18371,  8637, 15166,  7022,   622,  6111, 14937,
	 1431, 18103,  3159, 13531,  2293, 10446,  9931,  6366,
	 6620,  8598, 12875,  4721,  4940, 13043, 14731, 12099,
	 3721, 11276, 12040,  1938, 10239,  4015,  7648,  4342,
	 4978,  4494,  9598,  6095, 13683, 17944, 14903, 13180,
	16982,  1919,  8274,  8062,  5043,  3938, 13915, 14258,
	14374,  2675, 17124, 14600,  6828, 14192, 11943,  8723,
	 3964,  9481,  4925,   410,  9585, 17266, 18377, 15509,
	 1776, 13734, 13311, 14320,  4778,  1951, 13330,   829,
	 4696,  2936,  3063, 13785,  2836,  3249, 15651, 12737,
	 7014, 16004,  3311, 16201, 15258, 17233, 13163,  3959,
	10575,  3126, 14650,  6553, 10453,  7289,  8816,  7396,
	 7435, 12970,  3197,  4982,   926, 16751, 10739,  1156,
	 3302,  1248, 12854, 12841,  9091, 17805, 11374,  1396,
	15472,  2074, 14331,  7013,  5921,  2383,  7854,  4565,
	  889,   336, 16222,  7002, 13082, 16846,  7316, 15973,
	 4974, 12185,  8198,  6727,  7777,  9761, 12299,  6245,
	 2690, 16692,  1417, 17372,  8690,   817,    55,  9455,
	 5321, 13187,   617,  7200, 13215, 14574, 17000, 13392,
	 7510,  2403, 11014, 11565, 10282, 18110,  4265, 15981,
	 6936,  9298, 12877, 10092, 10689, 14345, 17684,  6974,
	 8376,  5488,   753,  3768, 17054, 10945, 15041, 15119,
	17295,  3779,  7226, 13907, 15579,  3711, 14915,  8540,
	 7619,  9411, 12298, 12776, 14540,  6076,  1492,  6805,
	 5348,  2747, 10067,  2934,  6998,  9902, 16144,   296,
	11816,  8675,  6608, 13238,  8706,  6919,  5655,  6927,
	15890,   200, 13167, 14701, 14809,   372,  3477,  1169,
	 1330, 14146,  4675, 11056,  9847,  1980, 17912, 10979,
	 5990, 12569, 13571, 12096, 14828,  5314, 10127, 16600,
	 1557,  6249, 10393, 12056,  4002,   932, 18225, 12839,
	10173, 11102,  2981, 14770,  6756,  5166,  5176,  1666,
	16686, 18063, 15272, 17964, 10391,  6685,  5549, 17192,
	13614,  7758,  9186,  5649,  9101,  7794, 14874, 18249,
	16568, 14535, 10838, 18175,  2061, 14132,  2463,  6156,
	10050,  3363, 15230, 15771,   260, 16209, 17268, 14219,
	 1166, 16116,  2574,   102,  7330, 16704,  3313,  3139,
	  262,  3147, 17968, 13903, 18183, 10646,  4665,  2634,
	14476,   101, 15958, 16917,  9967, 10879,  4613, 10452,
	 3667, 13723, 11573, 10470,  4240, 13359,  9360, 12101,
	14761,   499,  5110,  8753, 16027,  8670,  5818, 11488,
	15616,  1693,  9432,  2694, 15774,  2043,  9433, 14596,
	17710,  3065,  5012,  3636, 10821,   171,  8585,  4551,
	17832, 17335, 10846,  2793,  6624,   907, 14275,  4089,
	 9288,  3075,  6592,  7136, 18013, 14936,   464, 11061,
	16884, 15235, 10840,  5113,  7810, 15434,  5416,  1031,
	13756,  2006,  3587,  1646, 16001, 12579, 15151, 15596,
	17307, 17572, 11426, 12011, 14864,  9827,  4294, 10912,
	13629,  1958, 14436,  3279,  3599, 15439,  6206,  2781,
	13844, 17134, 15954,  6175, 16902,  8275, 17140,  2269,
	 8873,  3789,  8806, 17407,  4338,    43,  6794, 15050,
	 6101,  6515, 15555, 12991,  3606,  6588,  8656,  1675,
	11067, 15649,  2520,  2549, 18004, 18416, 15747, 12483,
	 6186,  4370,  8439, 17994,  9294,   755,  8692,  6188,
	 1411,  1259, 14592, 16691,  9629,  6397, 15344,  8557,
	11582,  6990, 16873, 13344, 13432, 16688,   785, 15972,
	  401, 16988, 11319, 10374, 11437, 13902,  2989, 17821,
	15804,  8876,  1500,  9856, 16861, 17984,  2790,  8747,
	 7413,  9188, 13930,  8458,  5309, 17827, 14850,  9096,
	17329,  2921,   693,  8535, 10481,  8751,   183,  2972,
	   70,  3655,  6067,  7373, 11190,  4955,  8704,  1548,
	10987,  3572, 11386, 15189, 10482,  2220,   533,  2814,
	12694,  7120,   547,  3545, 14763,  5870,  5810,  8437,
	 4286,  7861,  7027,  4833,  5760,  3193,  6803, 11570,
	14488, 13894,  1725, 15021,  9252, 16995, 12425, 12824,
	11476, 17255, 16639, 11659,  4387, 11818,  5511,  7308,
	 8756, 12163,  4722, 17460,  6701, 14144,  4359, 10356,
	15182, 15898,  4996, 15967,  1695,  8188,  3394,  8685,
	15400, 11481,  7564, 18389, 10211,  2553, 16281,  8766,
	 2417, 11654, 16465,  5207,  1915,  9142,  6662, 10791,
	 9297, 18028,  9742,  5714,  3652,  1090,  6323, 12840,
	15044, 13959, 11995,   905,   701, 11586,  5721, 18273,
	 1584, 14242,  1410,  7790, 16218, 14693, 17369, 18176,
	15311,  2884, 13280, 14018, 17050,   203, 13641, 15751,
	15492,   485,  2898,  3853, 17018,  6432,  2441,  2374,
	 4235,  9148,  7610, 12891, 13407, 14066, 10468,  1489,
	16595,  4095,  1855, 13909, 14277,  9460,  1607, 11493,
	12234,  6801,  5444,  2493, 17594,  4908,  1278,  3531,
	 7439,  5279,  4597,  4350,  6832,  6501, 13343,  8091,
	17473,  2540, 14227,  4216,  1901,  8411,  1762, 13003,
	16891,  6384, 13290,  4007,  8929,  6713,  9973,  8559,
	 2341, 10319,  8298, 17215,   299,  1129, 12485,  8057,
	14244,  3787,  8490, 16707, 17646, 15523,  1045, 13748,
	 8934, 10924, 11723,  7769, 11456,   411,  9639, 14819,
	13659,  8791,  6503, 16972, 11028, 12296,  7303,  8711,
	 2753, 10765,  5034,  7418,   328, 14493,  4202,  3475,
	11680, 12107, 14307, 16293, 10665,  5192,  9284, 10766,
	15244, 16502,  8263,  6171,  1073, 15210,  6890, 14796,
	15730, 12912, 12466,  3115, 18198,  4846,  9915,   264,
	 9447, 15327,  6943,   447,  3931,  3808, 11808,  5624,
	 3308, 17361, 14954, 11893, 17950,  2430, 15280,  2582,
	 6982,  3800, 10544,  2824,  4876,  7068, 10764,  3778,
	 6837, 10712, 15093,  7301,  2763,   754,  8534,  5838,
	 3212, 17615, 18220,  8628,  5237,  8801,  8083,  2039,
	11150,  8133, 13137,  7868,  2306, 17708, 14481,  4312,
};

/*
 * Table for inverse NTT, binary case:
 *   iGMb[x] = R*((1/g)^rev(x)) mod q
 * Since g = 19, 1/g = 5821 mod 18433.
 */
static const uint16_t iGMb[1024] = {
	10237,  1356, 13775, 11485, 16405,  8459, 11370,  9346,
	16558, 12320, 17111, 11095,  1965,  4047, 15542, 12704,
	18109,  3751,  4107,  2802,  8745,  8161, 17668, 17561,
	 1950, 16680, 13172, 17954,  1643,  2427,  1532, 14806,
	 6241,  4608,  9129,  9177, 13662, 10802,  1935, 10880,
	 3742, 15277,  1380, 17476,  9940, 15547,  3715,  4837,
	13417, 14378,    91,  4465, 14547,  2775, 12734, 14491,
	13804, 16554,  5938, 16479, 11782,  8900, 18256,  5292,
	16642,  7934, 11950,   128,  3794,  4662,  9596, 17709,
	 7707, 12327,  1128, 12201,  6522, 14952, 16661,  2992,
	  388,  8707,  6005, 11664,  9326,  5474, 17301, 16974,
	13367,  1189, 10624,  3532, 18286, 16892, 13640, 14584,
	15283, 17011, 18424, 14953, 14361,  4587,  1779,  5859,
	12814,  2414, 15415, 12752,  4709,  8235,  6702, 10820,
	 3276, 13276,  1484, 14679, 16032,  5552,  7735, 10865,
	11005,  3388,  6088,   747, 12061,  6182,  7039, 18240,
	12168,  4545,  5512, 17656,  9515,  4822, 10297,  6123,
	 6643, 12584, 17346, 15941,  7932,  7162, 18245,  7183,
	 2829,  6333,  4590,  5232, 12224,  1621, 14360, 16489,
	 8575,  3871,  9241,  3329, 13133,  2874, 10518, 11700,
	10193,  8920,  6823,  8452,   525,   237,  9218,   580,
	14092, 17316, 14576,  7844, 10153,  5742,   503,  4019,
	 4883,  1783, 15761,  5219, 17887, 10076,  5897,  6770,
	 3040,  1899,  1062,  5114,  1238, 11724, 11274,  9092,
	14102,  8894, 16156,  4344,  2032, 17665,  7695,  7687,
	16167,  2453, 10632,   481, 11665,   526, 18203,  9376,
	17776,  4022,  6792,  8754,   836,  3748,  3057,  2328,
	  882,  9246, 10325,  4661,  9988, 15674, 11299,  6470,
	 5999,  9344,  7759,  1712,    54,  2447,  8532, 17966,
	 8612,  5889, 15087,  8812, 18108, 15653,  3949,  3152,
	14406,  3554,  8889,  8542,  9529,  4091, 12509,  1223,
	 1366, 18207, 13065,  1158,   338, 13951, 16538, 10731,
	 5594,   208, 17501, 14431,  6377,  8040, 12184, 16876,
	 1833,  8306, 13119,  3605,  6337,  4862,  5864, 12443,
	 7454,   521, 16453,  8586,  7377, 13758,  4287, 17103,
	17264, 14956, 18061,  3624,  3732,  5266, 18233,  2543,
	11506, 12778, 11514,  9727,  5195, 11825,  9758,  6617,
	18137,  2289,  8531, 11435, 15499,  8366, 15686, 13085,
	11628, 16941, 12357,  3893,  5657,  6135,  9022, 10814,
	 9893,  3518, 14722,  2854,  4526, 11207, 14654,  1138,
	 3314,  3392,  7488,  1379, 14665, 17680, 12945, 10057,
	11459,   749,  4088,  7744,  8341,  5556,  9135, 11497,
	 2452, 14168,   323,  8151,  6868,  7419, 16030, 10923,
	 5041,  1433,  3859,  5218, 11233, 17816,  5246, 13112,
	 8978, 18378, 17616,  9743,  1061, 17016,  1741, 15743,
	12188,  6134,  8672, 10656, 11706, 10235,  6248, 13459,
	 2460, 11117,  1587,  5351, 11431,  2211, 18097, 17544,
	13868, 10579, 16050, 12512, 11420,  4102, 16359,  2961,
	17037,  7059,   628,  9342,  5592,  5579, 17185, 15131,
	17277,  7694,  1682, 17507, 13451, 15236,  5463, 10998,
	11037,  9617, 11144,  7980, 11880,  3783, 15307,  7858,
	14474,  5270,  1200,  3175,  2232, 15122,  2429, 11419,
	 5696,  2782, 15184, 15597,  4648, 15370, 15497, 13737,
	17604,  5103, 16482, 13655,  4113,  5122,  4699, 16657,
	 2924,    56,  1167,  8848, 18023, 13508,  8952, 14469,
	 9710,  6490,  4241, 11605,  3833,  1309, 15758,  4059,
	 4175,  4518, 14495, 13390, 10371, 10159, 16514,  1451,
	 5253,  3530,   489,  4750, 12338,  8835, 13939, 13455,
	14091, 10785, 14418,  8194, 16495,  6393,  7157, 14712,
	 6334,  3702,  5390, 13493, 13712,  5558,  9835, 11813,
	12067,  8502,  7987, 16140,  4902, 15274,   330, 17002,
	 3496, 12322, 17811, 11411,  3267,  9796,    62, 17829,
	 5146,  5167,  2016,  5334,  8911,  4760,  7030, 14760,
	 5212, 12254, 12444,   667, 14298, 17093, 10258,  9476,
	14121,  3952,   725, 16127, 10565,  5296, 10300,  7283,
	16394, 10350,  9632, 13196,  9805,   213,   818, 15221,
	12595,  9899, 17679, 15670, 11132,  3340,  7721, 11596,
	14655,  7669, 11365, 13557, 15609,  7889, 14633, 11451,
	15851,  3153, 16003,   483,  6540,  3479,  1072, 15125,
	12809,  6625, 14625, 14502, 17986, 11490,  3106,  8986,
	18169,  8518, 13587,   235, 15318,  5967,  5521,  2703,
	 3637, 11543,  3223, 17360, 12262, 10170,  1931,  3189,
	 7667,  9149, 13241,  7768,  2140,  4126,  6326,  6753,
	14958, 14231,  3940, 18105, 11015, 13399,  7668, 15680,
	 9722, 11130,  6137,  7405,  1461, 11930,  9642,  4774,
	 3614,  8794, 18022,  6977, 10664,  6710,  7509,  9499,
	 4685, 17388,  2910,   787,  1726,  9943, 14646,  4189,
	10376,  5948, 17304, 18134,  1218, 10135,  8114, 16092,
	 9874,  8460, 11720,  9504, 14426,  5143, 12049,  1542,
	 5430, 16671, 10022, 16532, 14217,  4206, 15893,   960,
	10342,  5090, 11932, 11601, 14083, 13836, 13154, 10994,
	14902, 17155, 13525,   839, 15940, 12989, 11632,  6199,
	 6940, 16826,  8973,  4156,  4524, 16578, 14338,  1838,
	16944,  7965,  4367,  5026,  5542, 10823,  9285, 14198,
	16059, 15992, 12001,  1415, 14580, 15535, 17948,  2941,
	 2682,  4792, 18230,  1383,  4415,  5153, 15549,  3122,
	  257,  1064,  3740,  2215, 10643, 17023,  4191, 16849,
	  160, 12712,  6847, 17732, 17528,  6438,  4474,  3389,
	 5593, 12110, 17343, 14781, 12719,  8691,   405,  9136,
	 7642, 11771,  9291, 16518, 13226,  1968,  6779, 16016,
	 9667,  2152, 15880,  8222,    44, 10869,  6952,  3033,
	 9748, 15039, 10245, 16738,  2466, 13437,  2535,  3251,
	 8077, 14074,  4289, 11732,   973, 13711,  6270,  9677,
	11125, 12922,  6615, 14046,  6774,  1794,  1178,  6957,
	 5609,  6008,  1438,  9181,  3412, 16708,  4539,  3945,
	 6863, 11630, 15240, 12673, 13600, 11406, 10572, 14147,
	 9996, 12623, 12563,  3670, 14888, 17886, 11313,  5739,
	15619, 17900, 16213,  7951,  3244,  7047, 14861,  7446,
	16885,  9729, 13478,  7243, 11060, 12366, 14778, 18363,
	15461, 18250,  9682,  7952,  9898, 17740, 15512,  1104,
	 9337,  3583,   606, 13124,  9975,  4503,  9245, 11020,
	 9686, 15643,   449,  1572,  8577, 16933,  9557,  2629,
	  612, 15444,  4531,  6996,  8059,  7114,  1445, 18032,
	 2461, 17648,  1745,  5001,  5089,  1560, 11443,  6851,
	 9876,  3089, 12036,  8804,  1742,  3841, 17174, 17022,
	12245,  9741, 17678,  9139,   439,  9994, 14063, 12247,
	 5950,  2686,    17,   429, 15884, 15913,  2784,  7366,
	16758,  9777, 11845, 14827,  5442,  2878, 11918, 12332,
	 3383, 11639, 18390, 14095,  1026,  9627, 14644,  9560,
	16164,  1293, 10158,  1531, 12258,  2479,  1299,  4589,
	15652, 12227,  2994, 14834, 15154,  3997, 16475,  4804,
	 7521, 14139,  8606,  3569,  6422,  7007,   861,  1126,
	 2837,  3282,  5854,  2432, 16787, 14846, 16427,  4677,
	17402, 13017,  2999, 10623, 13320,  7593,  3198,  1549,
	 7372, 17969,  3497,   420, 11297, 11841, 15358,  9145,
	14344,  4158, 17526, 11809, 15640,  7587,  1098,   601,
	13882,  9848, 18262,  7612, 14797, 13421, 15368,   723,
	 3837,  9000, 16390,  2659, 15739,  9001, 16740,  2817,
	 6945, 12615,  9763,  2406,  9680, 13323, 17934,  3672,
	 6332,  9073,  5074, 14193,  7963,  6860,  4710, 14766,
	 7981, 13820,  7554,  8466,  1516,  2475, 18332,  3957,
	15799, 13768,  7787,   250,  4530,   465, 15286, 18171,
	15294, 15120,  1729, 11103, 18331, 15859,  2317, 17267,
	 4214,  1165,  2224, 18173,  2662,  3203, 15070,  8383,
	12277, 15970,  4301, 16372,   258,  7595,  3898,  1865,
	  184,  3559, 10639,  9332, 12784,  9247, 10675,  4819,
	 1241, 12884, 11748,  8042,   469,  3161,   370,  1747,
	16767, 13257, 13267, 11677,  3663, 15452,  7331,  8260,
};


/* see inner.h */
uint32_t
Zf(mq_conv_small)(int x)
{
	/*
	 * If x < 0, the cast to uint32_t will set the high bit to 1.
	 */
	uint32_t y;

	y = (uint32_t)x;
	y += Q & -(y >> 31);
	return y;
}

/* see inner.h */
int32_t
Zf(mq_conv_signed)(uint32_t x)
{
	/*
	 * If x > Q/2, subtract Q.
	 */
	int32_t y;

	y = (int32_t)x;
	y -= (int32_t)(Q & -(((Q >> 1) - x) >> 31));
	return y;
}


/* see inner.h */
uint32_t
Zf(mq_add)(uint32_t x, uint32_t y)
{
	/*
	 * We compute x + y - q. If the result is negative, then the
	 * high bit will be set, and 'd >> 31' will be equal to 1;
	 * thus '-(d >> 31)' will be an all-one pattern. Otherwise,
	 * it will be an all-zero pattern. In other words, this
	 * implements a conditional addition of q.
	 */
	uint32_t d;

	d = x + y - Q;
	d += Q & -(d >> 31);
	return d;
}

/* see inner.h */
uint32_t
Zf(mq_sub)(uint32_t x, uint32_t y)
{
	/*
	 * As in mq_add(), we use a conditional addition to ensure the
	 * result is in the 0..q-1 range.
	 */
	uint32_t d;

	d = x - y;
	d += Q & -(d >> 31);
	return d;
}

/*
 * Division by 2 modulo q. Operand must be in the 0..q-1 range.
 */
static inline uint32_t
mq_rshift1(uint32_t x)
{
	x += Q & -(x & 1);
	return (x >> 1);
}

/* see inner.h */
uint32_t
Zf(mq_montymul)(uint32_t x, uint32_t y)
{
	uint32_t z, w;

	/*
	 * We compute x*y + k*q with a value of k chosen so that the 16
	 * low bits of the result are 0. We can then shift the value.
	 * After the shift, result may still be larger than q, but it
	 * will be lower than 2*q, so a conditional subtraction works.
	 */

	z = x * y;
	w = ((z * Q0I) & 0xFFFF) * Q;

	/*
	 * When adding z and w, the result will have its low 16 bits
	 * equal to 0. Since x, y and z are lower than q, the sum will
	 * be no more than (2^15 - 1) * q + (q - 1)^2, which will
	 * fit on 29 bits.
	 */
	z = (z + w) >> 16;

	/*
	 * After the shift, analysis shows that the value will be less
	 * than 2q. We do a subtraction then conditional subtraction to
	 * ensure the result is in the expected range.
	 */
	z -= Q;
	z += Q & -(z >> 31);
	return z;
}

/* see inner.h */
uint32_t
Zf(mq_mul)(uint32_t x, uint32_t y)
{
	x = Zf(mq_montymul)(x, R2);
	return Zf(mq_montymul)(x, y);
}

/*
 * Montgomery squaring (computes (x^2)/R).
 */
static inline uint32_t
mq_montysqr(uint32_t x)
{
	return Zf(mq_montymul)(x, x);
}

/* see inner.h */
static uint32_t
mq_div(uint32_t x, uint32_t y)
{
	/*
	 * We invert y by computing y^(q-2) mod q.
	 *
	 * We use the following addition chain for exponent e = 18431:
	 *
	 *   e0 = 1
	 *   e1 = 2 * e0 = 2
	 *   e2 = e1 + e0 = 3
	 *   e3 = 2 * e2 = 6
	 *   e4 = e3 + e0 = 7
	 *   e5 = 2 * e4 = 14
	 *   e6 = 2 * e5 = 28
	 *   e7 = e6 + e4 = 35
	 *   e8 = e7 + e6 = 63
	 *   e9 = 2 * e8 = 126
	 *   e10 = 2 * e9 = 252
	 *   e11 = e10 + e7 = 287
	 *   e12 = 2 * e11 = 574
	 *   e13 = 2 * e12 = 1148
	 *   e14 = 2 * e13 = 2296
	 *   e15 = 2 * e14 = 4592
	 *   e16 = 2 * e15 = 9184
	 *   e17 = 2 * e16 = 18368
	 *   e18 = e17 + e8 = 18431
	 *
	 * Additions on exponents are converted to Montgomery
	 * multiplications. We define all intermediate results as so
	 * many local variables, and let the C compiler work out which
	 * must be kept around.
	 */
	uint32_t y0, y1, y2, y3, y4, y5, y6, y7, y8, y9;
	uint32_t y10, y11, y12, y13, y14, y15, y16, y17, y18;

	y0 = Zf(mq_montymul)(y, R2);
	y1 = mq_montysqr(y0);
	y2 = Zf(mq_montymul)(y1, y0);
	y3 = mq_montysqr(y2);
	y4 = Zf(mq_montymul)(y3, y0);
	y5 = mq_montysqr(y4);
	y6 = mq_montysqr(y5);
	y7 = Zf(mq_montymul)(y6, y4);
	y8 = Zf(mq_montymul)(y7, y6);
	y9 = mq_montysqr(y8);
	y10 = mq_montysqr(y9);
	y11 = Zf(mq_montymul)(y10, y7);
	y12 = mq_montysqr(y11);
	y13 = mq_montysqr(y12);
	y14 = mq_montysqr(y13);
	y15 = mq_montysqr(y14);
	y16 = mq_montysqr(y15);
	y17 = mq_montysqr(y16);
	y18 = Zf(mq_montymul)(y17, y8);

	/*
	 * Final multiplication with x, which is not in Montgomery
	 * representation, computes the correct division result.
	 */
	return Zf(mq_montymul)(y18, x);
}

/* see inner.h */
void
Zf(mq_NTT)(uint16_t *a, unsigned logn)
{
	size_t n, t, m;

	n = MKN(logn);
	t = n;
	for (m = 1; m < n; m <<= 1) {
		size_t ht, i, j1;

		ht = t >> 1;
		for (i = 0, j1 = 0; i < m; i ++, j1 += t) {
			size_t j, j2;
			uint32_t s;

			s = GMb[m + i];
			j2 = j1 + ht;
			for (j = j1; j < j2; j ++) {
				uint32_t u, v;

				u = a[j];
				v = Zf(mq_montymul)(a[j + ht], s);
				a[j] = (uint16_t)Zf(mq_add)(u, v);
				a[j + ht] = (uint16_t)Zf(mq_sub)(u, v);
			}
		}
		t = ht;
	}
}

/* see inner.h */
void
Zf(mq_iNTT)(uint16_t *a, unsigned logn)
{
	size_t n, t, m;
	uint32_t ni;

	n = MKN(logn);
	t = 1;
	m = n;
	while (m > 1) {
		size_t hm, dt, i, j1;

		hm = m >> 1;
		dt = t << 1;
		for (i = 0, j1 = 0; i < hm; i ++, j1 += dt) {
			size_t j, j2;
			uint32_t s;

			j2 = j1 + t;
			s = iGMb[hm + i];
			for (j = j1; j < j2; j ++) {
				uint32_t u, v, w;

				u = a[j];
				v = a[j + t];
				a[j] = (uint16_t)Zf(mq_add)(u, v);
				w = Zf(mq_sub)(u, v);
				a[j + t] = (uint16_t) Zf(mq_montymul)(w, s);
			}
		}
		t = dt;
		m = hm;
	}

	/*
	 * To complete the inverse NTT, we must now divide all values by
	 * n (the vector size). We thus need the inverse of n, i.e. we
	 * need to divide 1 by 2 logn times. But we also want it in
	 * Montgomery representation, i.e. we also want to multiply it
	 * by R = 2^16. In the common case, this should be a simple right
	 * shift. The loop below is generic and works also in corner cases;
	 * its computation time is negligible.
	 */
	ni = R;
	for (m = n; m > 1; m >>= 1) {
		ni = mq_rshift1(ni);
	}
	for (m = 0; m < n; m ++) {
		a[m] = (uint16_t)Zf(mq_montymul)(a[m], ni);
	}
}

/* see inner.h */
void
Zf(mq_int8_to_NTT)(uint16_t *restrict p, const int8_t *restrict f, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u++) {
		p[u] = Zf(mq_conv_small)(f[u]);
	}
	Zf(mq_NTT)(p, logn);
}

/* see inner.h */
void
Zf(NTT_NTRU)(
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	uint16_t *restrict bf, uint16_t *restrict bg,
	uint16_t *restrict bF, uint16_t *restrict bG, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);

	Zf(mq_int8_to_NTT)(bf, f, logn);
	Zf(mq_int8_to_NTT)(bg, g, logn);
	Zf(mq_int8_to_NTT)(bF, F, logn);

	if (G == NULL) {
		/*
		 * Compute (in NTT representation) G = (1 + gF) / f.
		 */
		for (u = 0; u < n; u++) {
			bG[u] = Zf(mq_mul)(bg[u], bF[u]);
			bG[u] = Zf(mq_add)(1, bG[u]);
		}
		Zf(mq_poly_div)(bG, bf, logn);
	} else {
		Zf(mq_int8_to_NTT)(bG, G, logn);
	}
}

/* see inner.h */
void
Zf(mq_poly_muladj)(uint16_t *restrict a, const uint16_t *restrict b,
	unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u++) {
		a[u] = (uint16_t) Zf(mq_mul)(a[u], b[n - 1 - u]);
	}
}

/* see inner.h */
void
Zf(mq_poly_mulselfadj)(uint16_t *a, unsigned logn)
{
	size_t n, hn, u;

	n = MKN(logn);
	hn = n >> 1;
	for (u = 0; u < hn; u++) {
		a[u] = (uint16_t) Zf(mq_mul)(a[u], a[n - 1 - u]);
		a[n - 1 - u] = a[u];
	}
}

/* see inner.h */
void
Zf(mq_poly_tomonty)(uint16_t *f, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		f[u] = (uint16_t)Zf(mq_montymul)(f[u], R2);
	}
}

/* see inner.h */
void
Zf(mq_poly_div)(uint16_t *f, uint16_t *g, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		f[u] = (uint16_t)mq_div(f[u], g[u]);
	}
}

/* see inner.h */
int
Zf(mq_is_invertible)(int8_t *f, unsigned logn, uint8_t *restrict tmp)
{
	uint16_t *p;
	size_t n, u;
	int res;

	n = MKN(logn);
	p = (uint16_t*)tmp;
	res = 1;
	Zf(mq_int8_to_NTT)(p, f, logn);
	for (u = 0; u < n; u++) {
		// res &= (p[u] > 0)
		res &= 1U - ((p[u] - 1U) >> 15);
	}
	return res;
}

