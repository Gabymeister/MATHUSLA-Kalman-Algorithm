#!/usr/bin/env python
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import event
import ROOT as root
import analyzer
import sys
import util
import warnings
import numpy as np
import visualization

if __name__ == "__main__":

    warnings.filterwarnings("ignore")

#    directory = '/scratch/liaojiah/MATHUSLA/tracker/MATHUSLA-MLTracker/data/13_07_21/11_27_04/trees'
#    file = directory+'/'+'stat_2_0.root' #'10_48_42/stat_1_1.root' #'10_44_51/stat_2_1.root'

#    ev = event.Event(file, 192)
    # ev = event.Event(sys.argv[1], 12)

    '''
    ev.unused = True
    ev.used = True
    ev.vertex = True
    '''
    opt = 12

    if opt == 0: # visualise

        bracket = [34, 47, 70, 81, 96, 147, 176, 215, 234, 
		251, 268, 283, 295, 319, 395, 398, 404, 412, 
                422, 432, 448, 487, 494, 569, 570, 581, 587, 
                599, 605, 608, 629, 633, 642, 643, 674, 678, 
                735, 898, 920, 926, 951, 990]

        for evnum in bracket:
            ev = event.Event(sys.argv[1], evnum)
            ev.ExtractTruthPhysics()
            ev.Print()

            ev.kalman = True
            ev.scattered = True
            ev.merged = False

            ev.unused = True
            ev.used = True
            ev.vertex = True
            ev.merged = True
            ev.DrawRecoPoint()

            ev.merged = True

            ev.DrawRecoPoint()

    elif opt == 1:
        ev.PlotMomentum(num=1000, cut=1000)

    elif opt == 2:
        ev.PlotHitAccuracy(num=100)

    elif opt == 3:
        ev.kalman = True

        ev.DrawRecoPoint()

    elif opt == 4:

        ev = event.Event(sys.argv[1], 0)
        missed_inds_k, missed_inds_lin, high_mult = ev.find_missed_tracks(num=(0,1000), mult=3)

        # print(missed_inds_k)
        # print("Missed by kalman: ", len(missed_inds_k))
        # print(missed_inds_lin)
        # print("Missed by lin: ", len(missed_inds_lin))

        print(high_mult)


    elif opt == 5: # for automator

        ev = event.Event(sys.argv[1], 0)
        missed_inds_k, missed_inds_lin, high_mult = ev.find_missed_tracks(num=(0,1000), mult=10)

        '''
        print(missed_inds_k)
        print("Missed by kalman: ", len(missed_inds_k))
        print(missed_inds_lin)
        print("Missed by lin: ", len(missed_inds_lin))
        '''
        print("high multiplicity events: \n",high_mult)

        for evn in high_mult:

            ev = event.Event(sys.argv[1], evn)
#            ev = event.Event(file, evn)
            ev.writeDirectory = str(sys.argv[2])

            ev.unused = True
            ev.used = True
            ev.vertex = True

            ev.ExtractTruthPhysics()
            ev.Print()

            ev.kalman = True

            ev.DrawRecoPoint()




#            ev.kalman = False

#            ev.DrawRecoPoint()

    elif opt == 6: # for automator
        # missed_inds_k, missed_inds_lin, high_mult = ev.find_missed_tracks(num=(0,5000), mult=1)
        ev = event.Event(sys.argv[1], 0)
        bracket = ev.find_bracket(num=(0,5000))

        print(bracket)
        print(len(bracket))
        # bracket = [542, 867]
        ev.writeDirectory = str(sys.argv[2])



        for evn in bracket:
            print("Current Event:", evn)

            ev = event.Event(sys.argv[1], evn)

            ev.unused = True
            ev.used = True
            ev.vertex = True

            ev.ExtractTruthPhysics()
            ev.Print()
            ev.kalman = True
            ev.merged = False

            ev.DrawRecoPoint()
            ev.merged = True


            ev.DrawRecoPoint()

    elif opt == 7:

        num_truths = []

        for i in range(5000):
            ev = event.Event(sys.argv[1], i)

            ev.writeDirectory = str(sys.argv[2])

            ev.ExtractTruthPhysics()
            ev.Print()

            num_truths.append(ev.num_truth)

        num_truths = np.array(num_truths)

        num_truths[num_truths < 2] = 0
        num_truths = np.count_nonzero(num_truths)

        print("number of truth tracks \n",num_truths)

    elif opt == 8:

        inds = []

        for ind in inds:
            ev = event.Event(sys.argv[1],ind)
            ev.writeDirectory = str(sys.argv[2])

            ev.ExtractTruthPhysics()
            ev.Print()


            if ev.num_truth < 2:
                ev.kalman = True
                ev.merged = True
                ev.vertex = True
                ev.DrawRecoPoint()

    elif opt == 9: # plot vertex info

        def bool_func(tree, op=0):

            tree.SetBranchStatus("Vertex_k_m_x",1)
            tree.SetBranchStatus("Vertex_x",1)

            tree.SetBranchStatus("NumTracks_k_m",1)
            tree.SetBranchStatus("NumTracks_k_m",1)

            tree.SetBranchStatus("king_move_inds",1)


#            if len(util.unzip(tree.king_move_inds)) > 0: # event selection criteria
#                print(util.unzip(tree.king_move_inds))


            if op == 0 and len(tree.Vertex_k_m_x) >= 1: # event selection criteria
                return True

            if op == 1 and tree.NumTracks_k_m == 1:
                return True

            if op == 2 and tree.NumTracks_k_m >= 2 and len(tree.Vertex_k_m_x) == 0:
                return True

            else:
                return False

        ev = event.Event(sys.argv[1], 0)
        inds = ev.find_with_bool(bool_func, Op=0)

        print(inds)

#        inds = set(inds)

#        linear = set([0, 3, 7, 26, 32, 34, 36, 40, 41, 43, 44, 47, 48, 60, 64, 70, 87, 94, 106, 108, 110, 112, 114, 120, 123, 124, 127, 131, 134, 144, 145, 150, 153, 154, 155, 159, 161, 164, 165, 174, 180, 187, 188, 190, 192, 193, 195, 197, 209, 212, 222, 226, 227, 228, 230, 231, 239, 241, 245, 246, 256, 262, 270, 272, 280, 282, 283, 290, 291, 292, 293, 300, 301, 304, 318, 324, 328, 332, 336, 337, 338, 341, 343, 349, 352, 364, 369, 372, 373, 375, 377, 379, 384, 386, 388, 390, 393, 395, 403, 412, 423, 426, 432, 433, 438, 443, 445, 446, 449, 451, 453, 456, 458, 462, 463, 468, 473, 475, 482, 487, 489, 492, 495, 497, 500, 508, 510, 514, 526, 528, 534, 537, 542, 546, 554, 556, 560, 585, 591, 593, 597, 603, 607, 608, 617, 618, 619, 620, 634, 640, 641, 642, 644, 645, 647, 650, 651, 658, 665, 671, 677, 678, 679, 689, 693, 694, 697, 707, 709, 712, 714, 718, 719, 721, 726, 727, 732, 733, 734, 736, 737, 739, 744, 751, 755, 758, 760, 765, 766, 770, 778, 779, 780, 788, 813, 817, 818, 823, 826, 833, 834, 844, 850, 851, 856, 859, 860, 865, 866, 867, 868, 870, 875, 877, 881, 885, 886, 895, 899, 903, 906, 909, 926, 927, 928, 931, 936, 941, 944, 945, 949, 955, 960, 964, 971, 972, 973, 974, 979, 983, 985, 990, 993, 994, 995, 1000, 1006, 1018, 1022, 1024, 1026, 1039, 1040, 1045, 1048, 1057, 1058, 1063, 1077, 1082, 1090, 1095, 1100, 1109, 1112, 1117, 1118, 1125, 1127, 1130, 1133, 1134, 1138, 1141, 1147, 1149, 1150, 1154, 1155, 1158, 1159, 1160, 1166, 1167, 1172, 1183, 1190, 1191, 1201, 1204, 1211, 1212, 1213, 1217, 1218, 1223, 1226, 1228, 1239, 1240, 1242, 1252, 1253, 1254, 1256, 1259, 1261, 1262, 1263, 1268, 1273, 1281, 1285, 1292, 1304, 1305, 1311, 1313, 1316, 1321, 1324, 1332, 1336, 1346, 1347, 1348, 1352, 1356, 1362, 1364, 1365, 1370, 1372, 1375, 1377, 1380, 1382, 1390, 1391, 1394, 1397, 1398, 1399, 1403, 1417, 1420, 1421, 1425, 1426, 1432, 1436, 1440, 1442, 1443, 1445, 1451, 1457, 1462, 1464, 1466, 1471, 1476, 1478, 1479, 1482, 1483, 1493, 1495, 1497, 1500, 1501, 1505, 1506, 1510, 1515, 1517, 1524, 1531, 1534, 1540, 1556, 1559, 1562, 1565, 1568, 1570, 1578, 1579, 1580, 1583, 1584, 1588, 1592, 1593, 1594, 1600, 1608, 1613, 1614, 1615, 1627, 1628, 1629, 1634, 1636, 1639, 1643, 1645, 1646, 1655, 1660, 1661, 1663, 1670, 1676, 1681, 1683, 1687, 1696, 1698, 1700, 1704, 1705, 1711, 1713, 1715, 1717, 1725, 1727, 1735, 1743, 1744, 1751, 1762, 1764, 1773, 1775, 1778, 1781, 1786, 1791, 1794, 1803, 1805, 1806, 1808, 1811, 1816, 1821, 1823, 1826, 1833, 1841, 1844, 1854, 1855, 1856, 1866, 1869, 1873, 1876, 1877, 1886, 1887, 1891, 1896, 1900, 1904, 1905, 1911, 1919, 1924, 1931, 1933, 1937, 1939, 1942, 1946, 1952, 1954, 1960, 1961, 1962, 1963, 1966, 1970, 1973, 1977, 1979, 1980, 1987, 1990, 1994, 1995, 2011, 2014, 2015, 2022, 2030, 2032, 2033, 2034, 2056, 2060, 2074, 2075, 2086, 2087, 2094, 2097, 2102, 2103, 2110, 2116, 2118, 2119, 2123, 2125, 2130, 2131, 2132, 2133, 2134, 2135, 2140, 2142, 2145, 2148, 2149, 2154, 2155, 2156, 2157, 2162, 2178, 2180, 2181, 2184, 2195, 2196, 2206, 2209, 2213, 2214, 2222, 2232, 2237, 2240, 2241, 2243, 2246, 2249, 2255, 2260, 2263, 2266, 2272, 2274, 2276, 2279, 2280, 2281, 2283, 2285, 2287, 2289, 2291, 2302, 2306, 2307, 2310, 2314, 2316, 2322, 2323, 2324, 2327, 2331, 2335, 2341, 2342, 2343, 2344, 2347, 2349, 2359, 2364, 2365, 2371, 2373, 2374, 2381, 2382, 2391, 2394, 2395, 2396, 2398, 2401, 2402, 2406, 2408, 2419, 2421, 2423, 2426, 2428, 2438, 2439, 2444, 2447, 2453, 2461, 2462, 2463, 2480, 2491, 2496, 2499, 2501, 2503, 2508, 2512, 2516, 2523, 2531, 2532, 2536, 2539, 2546, 2548, 2554, 2557, 2561, 2566, 2568, 2570, 2575, 2583, 2592, 2599, 2606, 2608, 2614, 2615, 2616, 2618, 2620, 2625, 2631, 2639, 2641, 2648, 2660, 2664, 2669, 2674, 2684, 2686, 2691, 2692, 2694, 2698, 2702, 2710, 2713, 2716, 2718, 2720, 2722, 2726, 2733, 2735, 2741, 2742, 2743, 2748, 2751, 2763, 2764, 2768, 2769, 2773, 2775, 2778, 2789, 2790, 2797, 2798, 2799, 2801, 2802, 2804, 2806, 2807, 2810, 2811, 2823, 2830, 2833, 2836, 2838, 2840, 2842, 2844, 2848, 2857, 2859, 2861, 2865, 2868, 2870, 2876, 2881, 2885, 2891, 2897, 2910, 2917, 2919, 2921, 2927, 2932, 2933, 2937, 2941, 2943, 2948, 2954, 2960, 2965, 2969, 2974, 2976, 2981, 2985, 2992, 2994, 2996, 3000, 3001, 3002, 3003, 3004, 3015, 3017, 3026, 3030, 3041, 3043, 3047, 3052, 3055, 3056, 3062, 3065, 3068, 3069, 3071, 3072, 3076, 3077, 3079, 3080, 3082, 3083, 3085, 3086, 3090, 3091, 3092, 3096, 3104, 3108, 3115, 3116, 3117, 3124, 3130, 3133, 3134, 3136, 3137, 3143, 3147, 3150, 3151, 3154, 3155, 3162, 3165, 3168, 3170, 3172, 3183, 3184, 3185, 3187, 3189, 3190, 3191, 3193, 3200, 3205, 3207, 3211, 3215, 3216, 3219, 3224, 3226, 3237, 3241, 3243, 3244, 3246, 3255, 3268, 3270, 3278, 3284, 3290, 3291, 3299, 3307, 3309, 3320, 3329, 3330, 3331, 3332, 3333, 3340, 3342, 3348, 3355, 3367, 3369, 3376, 3377, 3381, 3382, 3386, 3389, 3390, 3392, 3393, 3398, 3404, 3405, 3411, 3414, 3417, 3419, 3422, 3438, 3440, 3447, 3448, 3452, 3455, 3459, 3460, 3470, 3473, 3475, 3477, 3479, 3481, 3486, 3487, 3492, 3495, 3499, 3503, 3505, 3506])
#        kalman = set([0, 2, 3, 12, 26, 34, 36, 40, 48, 56, 60, 64, 77, 81, 87, 98, 104, 108, 112, 120, 122, 124, 127, 132, 134, 137, 145, 146, 150, 151, 154, 155, 158, 159, 162, 164, 165, 169, 171, 173, 175, 187, 188, 192, 195, 196, 209, 210, 226, 227, 228, 231, 239, 262, 267, 270, 272, 280, 282, 290, 291, 292, 293, 295, 301, 302, 304, 306, 318, 319, 320, 328, 332, 333, 336, 337, 339, 341, 343, 349, 352, 364, 372, 373, 377, 379, 380, 387, 403, 410, 414, 422, 426, 431, 438, 445, 446, 451, 453, 456, 460, 476, 486, 488, 493, 495, 497, 514, 515, 516, 526, 528, 537, 538, 546, 552, 556, 557, 563, 580, 585, 586, 593, 598, 601, 603, 606, 615, 617, 618, 619, 622, 634, 635, 640, 641, 650, 662, 665, 671, 677, 678, 689, 690, 697, 709, 716, 727, 731, 737, 739, 743, 746, 752, 758, 763, 766, 774, 777, 778, 790, 799, 801, 802, 813, 816, 819, 823, 824, 826, 833, 851, 856, 858, 863, 866, 873, 875, 879, 883, 888, 899, 903, 909, 910, 917, 919, 926, 928, 929, 936, 944, 949, 955, 964, 972, 973, 983, 987, 990, 991, 995, 996, 1000, 1004, 1017, 1021, 1024, 1031, 1032, 1048, 1058, 1060, 1063, 1072, 1093, 1095, 1117, 1118, 1123, 1125, 1133, 1134, 1141, 1147, 1149, 1154, 1155, 1160, 1171, 1172, 1177, 1178, 1183, 1190, 1226, 1235, 1239, 1242, 1243, 1252, 1253, 1259, 1261, 1263, 1267, 1273, 1278, 1285, 1304, 1313, 1314, 1316, 1321, 1324, 1336, 1338, 1339, 1343, 1351, 1353, 1364, 1365, 1375, 1391, 1394, 1398, 1399, 1401, 1420, 1421, 1426, 1432, 1436, 1440, 1441, 1442, 1445, 1447, 1448, 1451, 1455, 1456, 1457, 1462, 1464, 1469, 1471, 1478, 1482, 1487, 1493, 1495, 1497, 1500, 1505, 1506, 1510, 1515, 1517, 1522, 1523, 1531, 1536, 1540, 1546, 1550, 1552, 1555, 1559, 1562, 1568, 1570, 1580, 1583, 1588, 1591, 1593, 1594, 1595, 1600, 1606, 1607, 1613, 1623, 1628, 1629, 1634, 1639, 1642, 1643, 1645, 1655, 1660, 1663, 1671, 1676, 1681, 1687, 1689, 1690, 1691, 1696, 1698, 1705, 1710, 1711, 1715, 1717, 1725, 1748, 1762, 1764, 1773, 1778, 1780, 1786, 1791, 1799, 1803, 1805, 1806, 1811, 1813, 1816, 1821, 1833, 1844, 1848, 1849, 1855, 1856, 1859, 1886, 1887, 1889, 1891, 1900, 1904, 1919, 1924, 1927, 1937, 1946, 1954, 1962, 1973, 1986, 1987, 1994, 1995, 1999, 2013, 2014, 2015, 2022, 2032, 2033, 2034, 2041, 2043, 2044, 2063, 2075, 2077, 2081, 2082, 2086, 2087, 2089, 2090, 2094, 2096, 2099, 2104, 2110, 2118, 2119, 2123, 2131, 2134, 2138, 2148, 2149, 2154, 2156, 2157, 2161, 2162, 2164, 2166, 2185, 2191, 2196, 2201, 2204, 2209, 2214, 2215, 2222, 2229, 2237, 2239, 2241, 2243, 2247, 2249, 2253, 2260, 2263, 2266, 2279, 2289, 2294, 2306, 2307, 2314, 2323, 2331, 2335, 2338, 2341, 2344, 2346, 2347, 2348, 2349, 2355, 2371, 2374, 2375, 2382, 2386, 2394, 2395, 2396, 2398, 2401, 2402, 2406, 2407, 2411, 2412, 2414, 2421, 2425, 2426, 2433, 2434, 2438, 2441, 2444, 2447, 2452, 2461, 2463, 2465, 2471, 2479, 2480, 2482, 2486, 2491, 2492, 2499, 2502, 2503, 2512, 2516, 2531, 2543, 2548, 2561, 2566, 2568, 2570, 2573, 2592, 2594, 2597, 2599, 2608, 2618, 2620, 2623, 2625, 2633, 2639, 2641, 2648, 2653, 2656, 2667, 2675, 2684, 2686, 2690, 2691, 2692, 2701, 2702, 2707, 2713, 2716, 2719, 2722, 2726, 2733, 2742, 2745, 2748, 2751, 2756, 2763, 2764, 2769, 2775, 2797, 2799, 2802, 2806, 2807, 2810, 2816, 2820, 2830, 2833, 2835, 2836, 2838, 2842, 2844, 2846, 2848, 2857, 2862, 2868, 2891, 2892, 2895, 2897, 2902, 2919, 2921, 2929, 2930, 2932, 2933, 2937, 2938, 2941, 2944, 2950, 2953, 2956, 2961, 2965, 2967, 2968, 2969, 2985, 2991, 2994, 2996, 3002, 3004, 3016, 3017, 3023, 3030, 3034, 3039, 3041, 3043, 3056, 3066, 3069, 3071, 3072, 3076, 3079, 3082, 3083, 3092, 3096, 3106, 3108, 3116, 3130, 3134, 3135, 3136, 3137, 3142, 3143, 3147, 3149, 3151, 3159, 3162, 3164, 3168, 3170, 3183, 3185, 3187, 3189, 3190, 3191, 3195, 3211, 3215, 3216, 3219, 3226, 3229, 3232, 3237, 3246, 3255, 3259, 3270, 3278, 3288, 3307, 3309, 3320, 3323, 3329, 3332, 3333, 3352, 3355, 3356, 3367, 3369, 3376, 3377, 3382, 3389, 3390, 3393, 3399, 3402, 3403, 3405, 3411, 3419, 3422, 3440, 3448, 3452, 3455, 3460, 3465, 3473, 3490, 3503, 3504, 3505, 3506, 3511])
#        print('linear - kalman \n',linear.difference(kalman))
#        print('kalman - linear \n',kalman.difference(linear))

#        inds = (linear.difference(kalman)).union(kalman.difference(linear))

#        inds = set()

        for ind in inds:

            print("Event number: ",ind)
            ev = event.Event(sys.argv[1],ind)
            ev.writeDirectory = str(sys.argv[2])

            ev.used = True
            ev.unused = True
            ev.vert_vel = True
            ev.kalman = True

            ev.ExtractTruthPhysics()
            ev.Print()

            ev.RecoVertex()

            ev.vert_vel = False
            ev.vertex = True
            ev.kalman = False
            ev.DrawRecoPoint()

            ev.VertexAccuracy()

    elif opt == 10:

        ev = event.Event(sys.argv[1], 0)
        ev.PlotBeta()

    elif opt == 11:
        ev = event.Event(sys.argv[1], 0)
        ev.PlotVertexBeta()

    elif opt == 12:

        ev = event.Event(sys.argv[1], 0)

        drs = []

        above = 0

        for ev_num in range(ev.Tree.GetEntries()):

            print("event {}".format(ev_num)) if ev_num % 100 == 0 else None

            ev = event.Event(sys.argv[1], ev_num)

            dr = ev.VertexAccuracy()

            drs.extend(dr)

            dr = np.array(dr)
            dr[dr < 1e4] = 0
#            print(np.count_nonzero(dr)," were above range")
            above += np.count_nonzero(dr)


        print(above," were above range")
        vis_engine = visualization.Histogram(drs, rng=(0,1e4))

#    input('pause')
    # ev.GetRecoInfo()
#    ev.DrawReco()
#    ev.DrawReco_k()ev.DrawRecoPoint()

    #h_mumu.StudyPassedEvents(0)
