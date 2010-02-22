#include <cctbx/eltbx/xray_scattering.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {

  namespace {

    /* International Tables for Crystallography
        Volume C
        Mathematical, Physical and Chemical Tables
        Edited by A.J.C. Wilson
        Kluwer Academic Publishers
        Dordrecht/Boston/London
        1992

      Table 6.1.1.4 (pp. 500-502)
        Coefficients for analytical approximation to the scattering factors
        of Tables 6.1.1.1 and 6.1.1.3

      [ Table 6.1.1.4 is a reprint of Table 2.2B, pp. 99-101,
        International Tables for X-ray Crystallography, Volume IV,
        The Kynoch Press: Birmingham, England, 1974.
        There is just one difference, see "Tl3+".
      ]
     */

    static const detail::raw_table_entry<4> it1992_raw_table[] =
    {
// BEGIN_COMPILED_IN_REFERENCE_DATA
      { "H",     { 0.489918, 0.262003, 0.196767, 0.049879 },
                 { 20.6593, 7.74039, 49.5519, 2.20159 },
                 0.001305 },
      /*
      { "H",     { 0.493002, 0.322912, 0.140191, 0.040810 },
                 { 10.5109, 26.1257, 3.14236, 57.7997 },
                 0.003038 },
      */
      { "H1-",   { 0.897661, 0.565616, 0.415815, 0.116973 },
                 { 53.1368, 15.1870, 186.576, 3.56709 },
                 0.002389 },
      { "He",    { 0.873400, 0.630900, 0.311200, 0.178000 },
                 { 9.10370, 3.35680, 22.9276, 0.982100 },
                 0.006400 },
      { "Li",    { 1.12820, 0.750800, 0.617500, 0.465300 },
                 { 3.95460, 1.05240, 85.3905, 168.261 },
                 0.037700 },
      { "Li1+",  { 0.696800, 0.788800, 0.341400, 0.156300 },
                 { 4.62370, 1.95570, 0.631600, 10.0953 },
                 0.016700 },
      { "Be",    { 1.59190, 1.12780, 0.539100, 0.702900 },
                 { 43.6427, 1.86230, 103.483, 0.542000 },
                 0.038500 },
      { "Be2+",  { 6.26030, 0.884900, 0.799300, 0.164700 },
                 { 0.002700, 0.831300, 2.27580, 5.11460 },
                 -6.1092 },
      { "B",     { 2.05450, 1.33260, 1.09790, 0.706800 },
                 { 23.2185, 1.02100, 60.3498, 0.140300 },
                 -0.19320 },
      { "C",     { 2.31000, 1.02000, 1.58860, 0.865000 },
                 { 20.8439, 10.2075, 0.568700, 51.6512 },
                 0.215600 },
      { "Cval",  { 2.26069, 1.56165, 1.05075, 0.839259 },
                 { 22.6907, 0.656665, 9.75618, 55.5949 },
                 0.286977 },
      { "N",     { 12.2126, 3.13220, 2.01250, 1.16630 },
                 { 0.005700, 9.89330, 28.9975, 0.582600 },
                 -11.529 },
      { "O",     { 3.04850, 2.28680, 1.54630, 0.867000 },
                 { 13.2771, 5.70110, 0.323900, 32.9089 },
                 0.250800 },
      { "O1-",   { 4.19160, 1.63969, 1.52673, -20.307 },
                 { 12.8573, 4.17236, 47.0179, -0.01404 },
                 21.9412 },
      { "O2-",   { 3.75040, 2.84294, 1.54298, 1.62091 },
                 { 16.5151, 6.59203, 0.319201, 43.3486 },
                 0.242060 },
                 /* Hovestreydt, Acta Cryst. (1983) A39, 268-269 */
      { "F",     { 3.53920, 2.64120, 1.51700, 1.02430 },
                 { 10.2825, 4.29440, 0.261500, 26.1476 },
                 0.277600 },
      { "F1-",   { 3.63220, 3.51057, 1.26064, 0.940706 },
                 { 5.27756, 14.7353, 0.442258, 47.3437 },
                 0.653396 },
      { "Ne",    { 3.95530, 3.11250, 1.45460, 1.12510 },
                 { 8.40420, 3.42620, 0.230600, 21.7184 },
                 0.351500 },
      { "Na",    { 4.76260, 3.17360, 1.26740, 1.11280 },
                 { 3.28500, 8.84220, 0.313600, 129.424 },
                 0.676000 },
      { "Na1+",  { 3.25650, 3.93620, 1.39980, 1.00320 },
                 { 2.66710, 6.11530, 0.200100, 14.0390 },
                 0.404000 },
      { "Mg",    { 5.42040, 2.17350, 1.22690, 2.30730 },
                 { 2.82750, 79.2611, 0.380800, 7.19370 },
                 0.858400 },
      { "Mg2+",  { 3.49880, 3.83780, 1.32840, 0.849700 },
                 { 2.16760, 4.75420, 0.185000, 10.1411 },
                 0.485300 },
      { "Al",    { 6.42020, 1.90020, 1.59360, 1.96460 },
                 { 3.03870, 0.742600, 31.5472, 85.0886 },
                 1.11510 },
      { "Al3+",  { 4.17448, 3.38760, 1.20296, 0.528137 },
                 { 1.93816, 4.14553, 0.228753, 8.28524 },
                 0.706786 },
      { "Si",    { 6.29150, 3.03530, 1.98910, 1.54100 }, /* "Siv" */
                 { 2.43860, 32.3337, 0.678500, 81.6937 },
                 1.14070 },
      { "Sival", { 5.66269, 3.07164, 2.62446, 1.39320 },
                 { 2.66520, 38.6634, 0.916946, 93.5458 },
                 1.24707 },
      { "Si4+",  { 4.43918, 3.20345, 1.19453, 0.416530 },
                 { 1.64167, 3.43757, 0.214900, 6.65365 },
                 0.746297 },
      { "P",     { 6.43450, 4.17910, 1.78000, 1.49080 },
                 { 1.90670, 27.1570, 0.526000, 68.1645 },
                 1.11490 },
      { "S",     { 6.90530, 5.20340, 1.43790, 1.58630 },
                 { 1.46790, 22.2151, 0.253600, 56.1720 },
                 0.866900 },
      { "Cl",    { 11.4604, 7.19640, 6.25560, 1.64550 },
                 { 0.010400, 1.16620, 18.5194, 47.7784 },
                 -9.5574 },
      { "Cl1-",  { 18.2915, 7.20840, 6.53370, 2.33860 },
                 { 0.006600, 1.17170, 19.5424, 60.4486 },
                 -16.378 },
      { "Ar",    { 7.48450, 6.77230, 0.653900, 1.64420 },
                 { 0.907200, 14.8407, 43.8983, 33.3929 },
                 1.44450 },
      { "K",     { 8.21860, 7.43980, 1.05190, 0.865900 },
                 { 12.7949, 0.774800, 213.187, 41.6841 },
                 1.42280 },
      { "K1+",   { 7.95780, 7.49170, 6.35900, 1.19150 },
                 { 12.6331, 0.767400, -0.00200, 31.9128 },
                 -4.9978 },
      { "Ca",    { 8.62660, 7.38730, 1.58990, 1.02110 },
                 { 10.4421, 0.659900, 85.7484, 178.437 },
                 1.37510 },
      { "Ca2+",  { 15.6348, 7.95180, 8.43720, 0.853700 },
                 { -0.00740, 0.608900, 10.3116, 25.9905 },
                 -14.875 },
      { "Sc",    { 9.18900, 7.36790, 1.64090, 1.46800 },
                 { 9.02130, 0.572900, 136.108, 51.3531 },
                 1.33290 },
      { "Sc3+",  { 13.4008, 8.02730, 1.65943, 1.57936 },
                 { 0.298540, 7.96290, -0.28604, 16.0662 },
                 -6.6667 },
      { "Ti",    { 9.75950, 7.35580, 1.69910, 1.90210 },
                 { 7.85080, 0.500000, 35.6338, 116.105 },
                 1.28070 },
      { "Ti2+",  { 9.11423, 7.62174, 2.27930, 0.087899 },
                 { 7.52430, 0.457585, 19.5361, 61.6558 },
                 0.897155 },
      { "Ti3+",  { 17.7344, 8.73816, 5.25691, 1.92134 },
                 { 0.220610, 7.04716, -0.15762, 15.9768 },
                 -14.652 },
      { "Ti4+",  { 19.5114, 8.23473, 2.01341, 1.52080 },
                 { 0.178847, 6.67018, -0.29263, 12.9464 },
                 -13.280 },
      { "V",     { 10.2971, 7.35110, 2.07030, 2.05710 },
                 { 6.86570, 0.438500, 26.8938, 102.478 },
                 1.21990 },
      { "V2+",   { 10.1060, 7.35410, 2.28840, 0.022300 },
                 { 6.88180, 0.440900, 20.3004, 115.122 },
                 1.22980 },
      { "V3+",   { 9.43141, 7.74190, 2.15343, 0.016865 },
                 { 6.39535, 0.383349, 15.1908, 63.9690 },
                 0.656565 },
      { "V5+",   { 15.6887, 8.14208, 2.03081, -9.5760 },
                 { 0.679003, 5.40135, 9.97278, 0.940464 },
                 1.71430 },
      { "Cr",    { 10.6406, 7.35370, 3.32400, 1.49220 },
                 { 6.10380, 0.392000, 20.2626, 98.7399 },
                 1.18320 },
      { "Cr2+",  { 9.54034, 7.75090, 3.58274, 0.509107 },
                 { 5.66078, 0.344261, 13.3075, 32.4224 },
                 0.616898 },
      { "Cr3+",  { 9.68090, 7.81136, 2.87603, 0.113575 },
                 { 5.59463, 0.334393, 12.8288, 32.8761 },
                 0.518275 },
      { "Mn",    { 11.2819, 7.35730, 3.01930, 2.24410 },
                 { 5.34090, 0.343200, 17.8674, 83.7543 },
                 1.08960 },
      { "Mn2+",  { 10.8061, 7.36200, 3.52680, 0.218400 },
                 { 5.27960, 0.343500, 14.3430, 41.3235 },
                 1.08740 },
      { "Mn3+",  { 9.84521, 7.87194, 3.56531, 0.323613 },
                 { 4.91797, 0.294393, 10.8171, 24.1281 },
                 0.393974 },
      { "Mn4+",  { 9.96253, 7.97057, 2.76067, 0.054447 },
                 { 4.84850, 0.283303, 10.4852, 27.5730 },
                 0.251877 },
      { "Fe",    { 11.7695, 7.35730, 3.52220, 2.30450 },
                 { 4.76110, 0.307200, 15.3535, 76.8805 },
                 1.03690 },
      { "Fe2+",  { 11.0424, 7.37400, 4.13460, 0.439900 },
                 { 4.65380, 0.305300, 12.0546, 31.2809 },
                 1.00970 },
      { "Fe3+",  { 11.1764, 7.38630, 3.39480, 0.072400 },
                 { 4.61470, 0.300500, 11.6729, 38.5566 },
                 0.970700 },
      { "Co",    { 12.2841, 7.34090, 4.00340, 2.34880 },
                 { 4.27910, 0.278400, 13.5359, 71.1692 },
                 1.01180 },
      { "Co2+",  { 11.2296, 7.38830, 4.73930, 0.710800 },
                 { 4.12310, 0.272600, 10.2443, 25.6466 },
                 0.932400 },
      { "Co3+",  { 10.3380, 7.88173, 4.76795, 0.725591 },
                 { 3.90969, 0.238668, 8.35583, 18.3491 },
                 0.286667 },
      { "Ni",    { 12.8376, 7.29200, 4.44380, 2.38000 },
                 { 3.87850, 0.256500, 12.1763, 66.3421 },
                 1.03410 },
      { "Ni2+",  { 11.4166, 7.40050, 5.34420, 0.977300 },
                 { 3.67660, 0.244900, 8.87300, 22.1626 },
                 0.861400 },
      { "Ni3+",  { 10.7806, 7.75868, 5.22746, 0.847114 },
                 { 3.54770, 0.223140, 7.64468, 16.9673 },
                 0.386044 },
      { "Cu",    { 13.3380, 7.16760, 5.61580, 1.67350 },
                 { 3.58280, 0.247000, 11.3966, 64.8126 },
                 1.19100 },
      { "Cu1+",  { 11.9475, 7.35730, 6.24550, 1.55780 },
                 { 3.36690, 0.227400, 8.66250, 25.8487 },
                 0.89000 },
      { "Cu2+",  { 11.8168, 7.11181, 5.78135, 1.14523 },
                 { 3.37484, 0.244078, 7.98760, 19.8970 },
                 1.14431 },
      { "Zn",    { 14.0743, 7.03180, 5.16520, 2.41000 },
                 { 3.26550, 0.233300, 10.3163, 58.7097 },
                 1.30410 },
      { "Zn2+",  { 11.9719, 7.38620, 6.46680, 1.39400 },
                 { 2.99460, 0.203100, 7.08260, 18.0995 },
                 0.780700 },
      { "Ga",    { 15.2354, 6.70060, 4.35910, 2.96230 },
                 { 3.06690, 0.241200, 10.7805, 61.4135 },
                 1.71890 },
      { "Ga3+",  { 12.6920, 6.69883, 6.06692, 1.00660 },
                 { 2.81262, 0.227890, 6.36441, 14.4122 },
                 1.53545 },
      { "Ge",    { 16.0816, 6.37470, 3.70680, 3.68300 },
                 { 2.85090, 0.251600, 11.4468, 54.7625 },
                 2.13130 },
      { "Ge4+",  { 12.9172, 6.70003, 6.06791, 0.859041 },
                 { 2.53718, 0.205855, 5.47913, 11.6030 },
                 1.45572 },
      { "As",    { 16.6723, 6.07010, 3.43130, 4.27790 },
                 { 2.63450, 0.264700, 12.9479, 47.7972 },
                 2.53100 },
      { "Se",    { 17.0006, 5.81960, 3.97310, 4.35430 },
                 { 2.40980, 0.272600, 15.2372, 43.8163 },
                 2.84090 },
      { "Br",    { 17.1789, 5.23580, 5.63770, 3.98510 },
                 { 2.17230, 16.5796, 0.260900, 41.4328 },
                 2.95570 },
      { "Br1-",  { 17.1718, 6.33380, 5.57540, 3.72720 },
                 { 2.20590, 19.3345, 0.287100, 58.1535 },
                 3.17760 },
      { "Kr",    { 17.3555, 6.72860, 5.54930, 3.53750 },
                 { 1.93840, 16.5623, 0.226100, 39.3972 },
                 2.82500 },
      { "Rb",    { 17.1784, 9.64350, 5.13990, 1.52920 },
                 { 1.78880, 17.3151, 0.274800, 164.934 },
                 3.48730 },
      { "Rb1+",  { 17.5816, 7.65980, 5.89810, 2.78170 },
                 { 1.71390, 14.7957, 0.160300, 31.2087 },
                 2.07820 },
      { "Sr",    { 17.5663, 9.81840, 5.42200, 2.66940 },
                 { 1.55640, 14.0988, 0.166400, 132.376 },
                 2.50640 },
      { "Sr2+",  { 18.0874, 8.13730, 2.56540, -34.193 },
                 { 1.49070, 12.6963, 24.5651, -0.01380 },
                 41.4025 },
      { "Y",     { 17.7760, 10.2946, 5.72629, 3.26588 },
                 { 1.40290, 12.8006, 0.125599, 104.354 },
                 1.91213 },
      { "Y3+",   { 17.9268, 9.15310, 1.76795, -33.108 },
                 { 1.35417, 11.2145, 22.6599, -0.01319 },
                 40.2602 },
      { "Zr",    { 17.8765, 10.9480, 5.41732, 3.65721 },
                 { 1.27618, 11.9160, 0.117622, 87.6627 },
                 2.06929 },
      { "Zr4+",  { 18.1668, 10.0562, 1.01118, -2.6479 },
                 { 1.21480, 10.1483, 21.6054, -0.10276 },
                 9.41454 },
      { "Nb",    { 17.6142, 12.0144, 4.04183, 3.53346 },
                 { 1.18865, 11.7660, 0.204785, 69.7957 },
                 3.75591 },
      { "Nb3+",  { 19.8812, 18.0653, 11.0177, 1.94715 },
                 { 0.019175, 1.13305, 10.1621, 28.3389 },
                 -12.912 },
      { "Nb5+",  { 17.9163, 13.3417, 10.7990, 0.337905 },
                 { 1.12446, 0.028781, 9.28206, 25.7228 },
                 -6.3934 },
      { "Mo",    { 3.70250, 17.2356, 12.8876, 3.74290 },
                 { 0.277200, 1.09580, 11.0040, 61.6584 },
                 4.38750 },
      { "Mo3+",  { 21.1664, 18.2017, 11.7423, 2.30951 },
                 { 0.014734, 1.03031, 9.53659, 26.6307 },
                 -14.421 },
      { "Mo5+",  { 21.0149, 18.0992, 11.4632, 0.740625 },
                 { 0.014345, 1.02238, 8.78809, 23.3452 },
                 -14.316 },
      { "Mo6+",  { 17.8871, 11.1750, 6.57891, 0.000000 },
                 { 1.03649, 8.48061, 0.058881, 0.000000 },
                 0.344941 },
      { "Tc",    { 19.1301, 11.0948, 4.64901, 2.71263 },
                 { 0.864132, 8.14487, 21.5707, 86.8472 },
                 5.40428 },
      { "Ru",    { 19.2674, 12.9182, 4.86337, 1.56756 },
                 { 0.808520, 8.43467, 24.7997, 94.2928 },
                 5.37874 },
      { "Ru3+",  { 18.5638, 13.2885, 9.32602, 3.00964 },
                 { 0.847329, 8.37164, 0.017662, 22.8870 },
                 -3.1892 },
      { "Ru4+",  { 18.5003, 13.1787, 4.71304, 2.18535 },
                 { 0.844582, 8.12534, 0.36495, 20.8504 },
                 1.42357 },
      { "Rh",    { 19.2957, 14.3501, 4.73425, 1.28918 },
                 { 0.751536, 8.21758, 25.8749, 98.6062 },
                 5.32800 },
      { "Rh3+",  { 18.8785, 14.1259, 3.32515, -6.1989 },
                 { 0.764252, 7.84438, 21.2487, -0.01036 },
                 11.8678 },
      { "Rh4+",  { 18.8545, 13.9806, 2.53464, -5.6526 },
                 { 0.760825, 7.62436, 19.3317, -0.01020 },
                 11.2835 },
      { "Pd",    { 19.3319, 15.5017, 5.29537, 0.605844 },
                 { 0.698655, 7.98929, 25.2052, 76.8986 },
                 5.26593 },
      { "Pd2+",  { 19.1701, 15.2096, 4.32234, 0.000000 },
                 { 0.696219, 7.55573, 22.5057, 0.000000 },
                 5.29160 },
      { "Pd4+",  { 19.2493, 14.7900, 2.89289, -7.9492 },
                 { 0.683839, 7.14833, 17.9144, 0.005127 },
                 13.0174 },
      { "Ag",    { 19.2808, 16.6885, 4.80450, 1.04630 },
                 { 0.644600, 7.47260, 24.6605, 99.8156 },
                 5.17900 },
      { "Ag1+",  { 19.1812, 15.9719, 5.27475, 0.357534 },
                 { 0.646179, 7.19123, 21.7326, 66.1147 },
                 5.21572 },
      { "Ag2+",  { 19.1643, 16.2456, 4.37090, 0.000000 },
                 { 0.645643, 7.18544, 21.4072, 0.000000 },
                 5.21404 },
      { "Cd",    { 19.2214, 17.6444, 4.46100, 1.60290 },
                 { 0.594600, 6.90890, 24.7008, 87.4825 },
                 5.06940 },
      { "Cd2+",  { 19.1514, 17.2535, 4.47128, 0.000000 },
                 { 0.597922, 6.80639, 20.2521, 0.000000 },
                 5.11937 },
      { "In",    { 19.1624, 18.5596, 4.29480, 2.03960 },
                 { 0.547600, 6.37760, 25.8499, 92.8029 },
                 4.93910 },
      { "In3+",  { 19.1045, 18.1108, 3.78897, 0.000000 },
                 { 0.551522, 6.32470, 17.3595, 0.000000 },
                 4.99635 },
      { "Sn",    { 19.1889, 19.1005, 4.45850, 2.46630 },
                 { 5.83030, 0.503100, 26.8909, 83.9571 },
                 4.78210 },
      { "Sn2+",  { 19.1094, 19.0548, 4.56480, 0.487000 },
                 { 0.503600, 5.83780, 23.3752, 62.2061 },
                 4.78610 },
      { "Sn4+",  { 18.9333, 19.7131, 3.41820, 0.019300 },
                 { 5.76400, 0.465500, 14.0049, -0.75830 },
                 3.91820 },
      { "Sb",    { 19.6418, 19.0455, 5.03710, 2.68270 },
                 { 5.30340, 0.460700, 27.9074, 75.2825 },
                 4.59090 },
      { "Sb3+",  { 18.9755, 18.9330, 5.10789, 0.288753 },
                 { 0.467196, 5.22126, 19.5902, 55.5113 },
                 4.69626 },
      { "Sb5+",  { 19.8685, 19.0302, 2.41253, 0.000000 },
                 { 5.44853, 0.467973, 14.1259, 0.000000 },
                 4.69263 },
      { "Te",    { 19.9644, 19.0138, 6.14487, 2.52390 },
                 { 4.81742, 0.420885, 28.5284, 70.8403 },
                 4.35200 },
      { "I",     { 20.1472, 18.9949, 7.51380, 2.27350 },
                 { 4.34700, 0.381400, 27.7660, 66.8776 },
                 4.07120 },
      { "I1-",   { 20.2332, 18.9970, 7.80690, 2.88680 },
                 { 4.35790, 0.381500, 29.5259, 84.9304 },
                 4.07140 },
      { "Xe",    { 20.2933, 19.0298, 8.97670, 1.99000 },
                 { 3.92820, 0.344000, 26.4659, 64.2658 },
                 3.71180 },
      { "Cs",    { 20.3892, 19.1062, 10.6620, 1.49530 },
                 { 3.56900, 0.310700, 24.3879, 213.904 },
                 3.33520 },
      { "Cs1+",  { 20.3524, 19.1278, 10.2821, 0.961500 },
                 { 3.55200, 0.308600, 23.7128, 59.4565 },
                 3.27910 },
      { "Ba",    { 20.3361, 19.2970, 10.8880, 2.69590 },
                 { 3.21600, 0.275600, 20.2073, 167.202 },
                 2.77310 },
      { "Ba2+",  { 20.1807, 19.1136, 10.9054, 0.77634 },
                 { 3.21367, 0.283310, 20.0558, 51.7460 },
                 3.02902 },
      { "La",    { 20.5780, 19.5990, 11.3727, 3.28719 },
                 { 2.94817, 0.244475, 18.7726, 133.124 },
                 2.14678 },
      { "La3+",  { 20.2489, 19.3763, 11.6323, 0.336048 },
                 { 2.92070, 0.250698, 17.8211, 54.9453 },
                 2.40860 },
      { "Ce",    { 21.1671, 19.7695, 11.8513, 3.33049 },
                 { 2.81219, 0.226836, 17.6083, 127.113 },
                 1.86264 },
      { "Ce3+",  { 20.8036, 19.5590, 11.9369, 0.612376 },
                 { 2.77691, 0.231540, 16.5408, 43.1692 },
                 2.09013 },
      { "Ce4+",  { 20.3235, 19.8186, 12.1233, 0.144583 },
                 { 2.65941, 0.218850, 15.7992, 62.2355 },
                 1.59180 },
      { "Pr",    { 22.0440, 19.6697, 12.3856, 2.82428 },
                 { 2.77393, 0.222087, 16.7669, 143.644 },
                 2.05830 },
      { "Pr3+",  { 21.3727, 19.7491, 12.1329, 0.975180 },
                 { 2.64520, 0.214299, 15.3230, 36.4065 },
                 1.77132 },
      { "Pr4+",  { 20.9413, 20.0539, 12.4668, 0.296689 },
                 { 2.54467, 0.202481, 14.8137, 45.4643 },
                 1.24285 },
      { "Nd",    { 22.6845, 19.6847, 12.7740, 2.85137 },
                 { 2.66248, 0.210628, 15.8850, 137.903 },
                 1.98486 },
      { "Nd3+",  { 21.9610, 19.9339, 12.1200, 1.51031 },
                 { 2.52722, 0.199237, 14.1783, 30.8717 },
                 1.47588 },
      { "Pm",    { 23.3405, 19.6095, 13.1235, 2.87516 },
                 { 2.56270, 0.202088, 15.1009, 132.721 },
                 2.02876 },
      { "Pm3+",  { 22.5527, 20.1108, 12.0671, 2.07492 },
                 { 2.41740, 0.185769, 13.1275, 27.4491 },
                 1.19499 },
      { "Sm",    { 24.0042, 19.4258, 13.4396, 2.89604 },
                 { 2.47274, 0.196451, 14.3996, 128.007 },
                 2.20963 },
      { "Sm3+",  { 23.1504, 20.2599, 11.9202, 2.71488 },
                 { 2.31641, 0.174081, 12.1571, 24.8242 },
                 0.954586 },
      { "Eu",    { 24.6274, 19.0886, 13.7603, 2.92270 },
                 { 2.38790, 0.194200, 13.7546, 123.174 },
                 2.57450 },
      { "Eu2+",  { 24.0063, 19.9504, 11.8034, 3.87243 },
                 { 2.27783, 0.173530, 11.6096, 26.5156 },
                 1.36389 },
      { "Eu3+",  { 23.7497, 20.3745, 11.8509, 3.26503 },
                 { 2.22258, 0.163940, 11.3110, 22.9966 },
                 0.759344 },
      { "Gd",    { 25.0709, 19.0798, 13.8518, 3.54545 },
                 { 2.25341, 0.181951, 12.9331, 101.398 },
                 2.41960 },
      { "Gd3+",  { 24.3466, 20.4208, 11.8708, 3.71490 },
                 { 2.13553, 0.155525, 10.5782, 21.7029 },
                 0.645089 },
      { "Tb",    { 25.8976, 18.2185, 14.3167, 2.95354 },
                 { 2.24256, 0.196143, 12.6648, 115.362 },
                 3.58324 },
      { "Tb3+",  { 24.9559, 20.3271, 12.2471, 3.77300 },
                 { 2.05601, 0.149525, 10.0499, 21.2773 },
                 0.691967 },
      { "Dy",    { 26.5070, 17.6383, 14.5596, 2.96577 },
                 { 2.18020, 0.202172, 12.1899, 111.874 },
                 4.29728 },
      { "Dy3+",  { 25.5395, 20.2861, 11.9812, 4.50073 },
                 { 1.98040, 0.143384, 9.34972, 19.5810 },
                 0.689690 },
      { "Ho",    { 26.9049, 17.2940, 14.5583, 3.63837 },
                 { 2.07051, 0.197940, 11.4407, 92.6566 },
                 4.56796 },
      { "Ho3+",  { 26.1296, 20.0994, 11.9788, 4.93676 },
                 { 1.91072, 0.139358, 8.80018, 18.5908 },
                 0.852795 },
      { "Er",    { 27.6563, 16.4285, 14.9779, 2.98233 },
                 { 2.07356, 0.223545, 11.3604, 105.703 },
                 5.92046 },
      { "Er3+",  { 26.7220, 19.7748, 12.1506, 5.17379 },
                 { 1.84659, 0.137290, 8.36225, 17.8974 },
                 1.17613 },
      { "Tm",    { 28.1819, 15.8851, 15.1542, 2.98706 },
                 { 2.02859, 0.238849, 10.9975, 102.961 },
                 6.75621 },
      { "Tm3+",  { 27.3083, 19.3320, 12.3339, 5.38348 },
                 { 1.78711, 0.136974, 7.96778, 17.2922 },
                 1.63929 },
      { "Yb",    { 28.6641, 15.4345, 15.3087, 2.98963 },
                 { 1.98890, 0.257119, 10.6647, 100.417 },
                 7.56672 },
      { "Yb2+",  { 28.1209, 17.6817, 13.3335, 5.14657 },
                 { 1.78503, 0.159970, 8.18304, 20.3900 },
                 3.70983 },
      { "Yb3+",  { 27.8917, 18.7614, 12.6072, 5.47647 },
                 { 1.73272, 0.138790, 7.64412, 16.8153 },
                 2.26001 },
      { "Lu",    { 28.9476, 15.2208, 15.1000, 3.71601 },
                 { 1.90182, 9.98519, 0.261033, 84.3298 },
                 7.97628 },
      { "Lu3+",  { 28.4628, 18.1210, 12.8429, 5.59415 },
                 { 1.68216, 0.142292, 7.33727, 16.3535 },
                 2.97573 },
      { "Hf",    { 29.1440, 15.1726, 14.7586, 4.30013 },
                 { 1.83262, 9.59990, 0.275116, 72.0290 },
                 8.58154 },
      { "Hf4+",  { 28.8131, 18.4601, 12.7285, 5.59927 },
                 { 1.59136, 0.128903, 6.76232, 14.0366 },
                 2.39699 },
      { "Ta",    { 29.2024, 15.2293, 14.5135, 4.76492 },
                 { 1.77333, 9.37046, 0.295977, 63.3644 },
                 9.24354 },
      { "Ta5+",  { 29.1587, 18.8407, 12.8268, 5.38695 },
                 { 1.50711, 0.116741, 6.31524, 12.4244 },
                 1.78555 },
      { "W",     { 29.0818, 15.4300, 14.4327, 5.11982 },
                 { 1.72029, 9.22590, 0.321703, 57.0560 },
                 9.88750 },
      { "W6+",   { 29.4936, 19.3763, 13.0544, 5.06412 },
                 { 1.42755, 0.104621, 5.93667, 11.1972 },
                 1.01074 },
      { "Re",    { 28.7621, 15.7189, 14.5564, 5.44174 },
                 { 1.67191, 9.09227, 0.350500, 52.0861 },
                 10.4720 },
      { "Os",    { 28.1894, 16.1550, 14.9305, 5.67589 },
                 { 1.62903, 8.97948, 0.382661, 48.1647 },
                 11.0005 },
      { "Os4+",  { 30.4190, 15.2637, 14.7458, 5.06795 },
                 { 1.37113, 6.84706, 0.165191, 18.0030 },
                 6.49804 },
      { "Ir",    { 27.3049, 16.7296, 15.6115, 5.83377 },
                 { 1.59279, 8.86553, 0.417916, 45.0011 },
                 11.4722 },
      { "Ir3+",  { 30.4156, 15.8620, 13.6145, 5.82008 },
                 { 1.34323, 7.10909, 0.204633, 20.3254 },
                 8.27903 },
      { "Ir4+",  { 30.7058, 15.5512, 14.2326, 5.53672 },
                 { 1.30923, 6.71983, 0.167252, 17.4911 },
                 6.96824 },
      { "Pt",    { 27.0059, 17.7639, 15.7131, 5.78370 },
                 { 1.51293, 8.81174, 0.424593, 38.6103 },
                 11.6883 },
      { "Pt2+",  { 29.8429, 16.7224, 13.2153, 6.35234 },
                 { 1.32927, 7.38979, 0.263297, 22.9426 },
                 9.85329 },
      { "Pt4+",  { 30.9612, 15.9829, 13.7348, 5.92034 },
                 { 1.24813, 6.60834, 0.168640, 16.9392 },
                 7.39534 },
      { "Au",    { 16.8819, 18.5913, 25.5582, 5.86000 },
                 { 0.461100, 8.62160, 1.48260, 36.3956 },
                 12.0658 },
      { "Au1+",  { 28.0109, 17.8204, 14.3359, 6.58077 },
                 { 1.35321, 7.73950, 0.356752, 26.4043 },
                 11.2299 },
      { "Au3+",  { 30.6886, 16.9029, 12.7801, 6.52354 },
                 { 1.21990, 6.82872, 0.212867, 18.6590 },
                 9.09680 },
      { "Hg",    { 20.6809, 19.0417, 21.6575, 5.96760 },
                 { 0.545000, 8.44840, 1.57290, 38.3246 },
                 12.6089 },
      { "Hg1+",  { 25.0853, 18.4973, 16.8883, 6.48216 },
                 { 1.39507, 7.65105, 0.443378, 28.2262 },
                 12.0205 },
      { "Hg2+",  { 29.5641, 18.0600, 12.8374, 6.89912 },
                 { 1.21152, 7.05639, 0.284738, 20.7482 },
                 10.6268 },
      { "Tl",    { 27.5446, 19.1584, 15.5380, 5.52593 },
                 { 0.655150, 8.70751, 1.96347, 45.8149 },
                 13.1746 },
      { "Tl1+",  { 21.3985, 20.4723, 18.7478, 6.82847 },
                 { 1.47110, 0.517394, 7.43463, 28.8482 },
                 12.5258 },
      { "Tl3+",  { 30.8695, 18.3841, 11.9328, 7.00574 },
                 { 1.10080, 6.53852, 0.219074, 17.2114 },
                 9.80270 },
                 /* IT Vol IV 1974: a2 = 18.3841
                    IT Vol C  1992: a2 = 18.3481 */
      { "Pb",    { 31.0617, 13.0637, 18.4420, 5.96960 },
                 { 0.690200, 2.35760, 8.61800, 47.2579 },
                 13.4118 },
      { "Pb2+",  { 21.7886, 19.5682, 19.1406, 7.01107 },
                 { 1.33660, 0.488383, 6.77270, 23.8132 },
                 12.4734 },
      { "Pb4+",  { 32.1244, 18.8003, 12.0175, 6.96886 },
                 { 1.00566, 6.10926, 0.147041, 14.7140 },
                 8.08428 },
      { "Bi",    { 33.3689, 12.9510, 16.5877, 6.46920 },
                 { 0.704000, 2.92380, 8.79370, 48.0093 },
                 13.5782 },
      { "Bi3+",  { 21.8053, 19.5026, 19.1053, 7.10295 },
                 { 1.23560, 6.24149, 0.469999, 20.3185 },
                 12.4711 },
      { "Bi5+",  { 33.5364, 25.0946, 19.2497, 6.91555 },
                 { 0.916540, 0.39042, 5.71414, 12.8285 },
                 -6.7994 },
      { "Po",    { 34.6726, 15.4733, 13.1138, 7.02588 },
                 { 0.700999, 3.55078, 9.55642, 47.0045 },
                 13.6770 },
      { "At",    { 35.3163, 19.0211, 9.49887, 7.42518 },
                 { 0.685870, 3.97458, 11.3824, 45.4715 },
                 13.7108 },
      { "Rn",    { 35.5631, 21.2816, 8.00370, 7.44330 },
                 { 0.663100, 4.06910, 14.0422, 44.2473 },
                 13.6905 },
      { "Fr",    { 35.9299, 23.0547, 12.1439, 2.11253 },
                 { 0.646453, 4.17619, 23.1052, 150.645 },
                 13.7247 },
      { "Ra",    { 35.7630, 22.9064, 12.4739, 3.21097 },
                 { 0.616341, 3.87135, 19.9887, 142.325 },
                 13.6211 },
      { "Ra2+",  { 35.2150, 21.6700, 7.91342, 7.65078 },
                 { 0.604909, 3.57670, 12.6010, 29.8436 },
                 13.5431 },
      { "Ac",    { 35.6597, 23.1032, 12.5977, 4.08655 },
                 { 0.589092, 3.65155, 18.5990, 117.020 },
                 13.5266 },
      { "Ac3+",  { 35.1736, 22.1112, 8.19216, 7.05545 },
                 { 0.579689, 3.41437, 12.9187, 25.9443 },
                 13.4637 },
      { "Th",    { 35.5645, 23.4219, 12.7473, 4.80703 },
                 { 0.563359, 3.46204, 17.8309, 99.1722 },
                 13.4314 },
      { "Th4+",  { 35.1007, 22.4418, 9.78554, 5.29444 },
                 { 0.555054, 3.24498, 13.4661, 23.9533 },
                 13.3760 },
      { "Pa",    { 35.8847, 23.2948, 14.1891, 4.17287 },
                 { 0.547751, 3.41519, 16.9235, 105.251 },
                 13.4287 },
      { "U",     { 36.0228, 23.4128, 14.9491, 4.18800 },
                 { 0.529300, 3.32530, 16.0927, 100.613 },
                 13.3966 },
      { "U3+",   { 35.5747, 22.5259, 12.2165, 5.37073 },
                 { 0.520480, 3.12293, 12.7148, 26.3394 },
                 13.3092 },
      { "U4+",   { 35.3715, 22.5326, 12.0291, 4.79840 },
                 { 0.516598, 3.05053, 12.5723, 23.4582 },
                 13.2671 },
      { "U6+",   { 34.8509, 22.7584, 14.0099, 1.21457 },
                 { 0.507079, 2.89030, 13.1767, 25.2017 },
                 13.1665 },
      { "Np",    { 36.1874, 23.5964, 15.6402, 4.18550 },
                 { 0.511929, 3.25396, 15.3622, 97.4908 },
                 13.3573 },
      { "Np3+",  { 35.7074, 22.6130, 12.9898, 5.43227 },
                 { 0.502322, 3.03807, 12.1449, 25.4928 },
                 13.2544 },
      { "Np4+",  { 35.5103, 22.5787, 12.7766, 4.92159 },
                 { 0.498626, 2.96627, 11.9484, 22.7502 },
                 13.2116 },
      { "Np6+",  { 35.0136, 22.7286, 14.3884, 1.75669 },
                 { 0.489810, 2.81099, 12.3300, 22.6581 },
                 13.1130 },
      { "Pu",    { 36.5254, 23.8083, 16.7707, 3.47947 },
                 { 0.499384, 3.26371, 14.9455, 105.980 },
                 13.3812 },
      { "Pu3+",  { 35.8400, 22.7169, 13.5807, 5.66016 },
                 { 0.484938, 2.96118, 11.5331, 24.3992 },
                 13.1991 },
      { "Pu4+",  { 35.6493, 22.6460, 13.3595, 5.18831 },
                 { 0.481422, 2.89020, 11.3160, 21.8301 },
                 13.1555 },
      { "Pu6+",  { 35.1736, 22.7181, 14.7635, 2.28678 },
                 { 0.473204, 2.73848, 11.5530, 20.9303 },
                 13.0582 },
      { "Am",    { 36.6706, 24.0992, 17.3415, 3.49331 },
                 { 0.483629, 3.20647, 14.3136, 102.273 },
                 13.3592 },
      { "Cm",    { 36.6488, 24.4096, 17.3990, 4.21665 },
                 { 0.465154, 3.08997, 13.4346, 88.4834 },
                 13.2887 },
      { "Bk",    { 36.7881, 24.7736, 17.8919, 4.23284 },
                 { 0.451018, 3.04619, 12.8946, 86.0030 },
                 13.2754 },
      { "Cf",    { 36.9185, 25.1995, 18.3317, 4.24391 },
                 { 0.437533, 3.00775, 12.4044, 83.7881 },
                 13.2674 },
      { 0,       { 0., 0., 0., 0. },
                 { 0., 0., 0., 0. },
                 0. }
// END_COMPILED_IN_REFERENCE_DATA
    };

  } // namespace <anonymous>

  it1992::it1992(std::string const& label, bool exact)
  :
    base<4>(it1992_raw_table, "IT1992", label, exact)
  {}

  it1992_iterator::it1992_iterator()
  :
    current_("H", true)
  {}

  it1992
  it1992_iterator::next()
  {
    it1992 result = current_;
    current_.next_entry();
    return result;
  }

}}} // namespace cctbx::eltbx::xray_scattering
