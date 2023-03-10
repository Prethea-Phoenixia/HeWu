"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

Implements the test cases found in the Appendix A of:
And Table 2. for quantities for 1-kT free-air blast wave:

PSR Report 1419-3
"AIRBLAST FROM NUCLEAR BURSTS—ANALYTIC APPROXIMATIONS"
Harold L. Brode, July 1987
Technical Report
CONTRACT No. DNA 001-85-C-0089  

"""


"""
test cases for airburst, in the format of:
    X: scaled ground range, ft/kT^(1/3)
    Y: scaled burst height, ft/kT^(1/3)
    sigma-tau: scaled time - (scaled) time of arrival, ms/kT^(1/3) 
    DeltaP: overpressure, psi.
"""

abtests = (
    (18.2, 0.0, 0.042, 107926.6),
    (28.1, 0.0, 0.363, 10113.7),
    (37.9, 0.0, 1.22, 2069.0),
    (47.1, 0.0, 0.0045, 27696.0),
    (66.3, 0.0, 6.39, 306.5),
    (80.3, 0.0, 19.8, 78.34),
    (0.0, 25.0, 0.137, 55580.0),
    (19.7, 25.0, 25.3, 9.284),
    (32.8, 25.0, 0.00417, 255273.1),
    (49.2, 25.0, 0.046, 10462.5),
    (64.7, 25.0, 0.189, 15346.0),
    (83.0, 25.0, 57.8, 7.799),
    (97.6, 25.0, 124.0, 1.265),
    (0.0, 50.0, 0.0416, 110106.0),
    (6.72, 50.0, 0.116, 43739.3),
    (24.2, 50.0, 3.07, 270.9),
    (41.5, 50.0, 56.1, 0.1174),
    (55.9, 50.0, 0.137, 21674.8),
    (72.6, 50.0, 0.251, 15168.9),
    (90.6, 50.0, 11.7, 70.46),
    (14.9, 22.5, 0.0139, 201704.4),
    (33.1, 18.9, 1.35, 1397.0),
    (37.7, 8.22, 0.214, 16365.0),
    (82.8, 1410.0, 6.95, 8.374),
    (157.0, 257.0, 0.448, 375.8),
    (393.0, 682.0, 10.9, 24.23),
    (1250.0, 220.0, 22.1, 6.593),
    (1160.0, 2520.0, 0.538, 3.1185),
    (3140.0, 4260.0, 13.4, 1.233),
)
"""
test caes for 1 kT free-air bursts, in the format of:
    ΔP_s  peak overpressure in psi
    r    shock radius in ft / kT^(1/3) # original table annotated sr, should be r if scaled
    T     time of arrival in ms/kT^(1/3)
    D^+_p overpressure positive phase duration in ms/kT^(1/3)
    I^+_p overpressure impulse in psi-ms/kT^(1/3)
    D^+_u dynamic pressure duration ms/kT^(1/3)
    I^+_u dynamic impulse in psi-ms/kT^(1/3)
    θ_m   maximum fireball temperature in 10^3 deg C


"""
fatest = (
    (0.1, 26310, 22760, 332, 45.8, 315, 0.0274, 0.0005),
    (0.15, 18160, 15550, 327, 56.1, 314, 0.0579, 0.0008),
    (0.2, 14020, 11890, 322, 64.7, 313, 0.0977, 0.0011),
    (0.3, 9816, 8161, 315, 79.2, 311, 0.203, 0.0016),
    (0.4, 7667, 6261, 308, 91.5, 310, 0.338, 0.002),
    (0.6, 5467, 4317, 297, 112, 306, 0.692, 0.003),
    (0.8, 4335, 3320, 288, 129, 302, 1.15, 0.004),
    (1.0, 3639, 2709, 279, 144, 298, 1.70, 0.005),
    (1.5, 2680, 1872, 262, 177, 289, 3.46, 0.007),
    (2, 2178, 1440, 248, 204, 280, 5.75, 0.010),
    (3, 1650, 994.0, 226, 249, 264, 11.8, 0.015),
    (4, 1370, 763.4, 209, 288, 249, 19.4, 0.019),
    (6, 1069, 525.9, 183, 352, 226, 38.8, 0.028),
    (7, 977.3, 456.3, 173, 380, 217, 49.9, 0.032),
    (10, 801.0, 328.9, 151, 453, 200, 87.4, 0.045),
    (15, 648.0, 227.1, 126, 553, 191, 156, 0.064),
    (20, 562.2, 174.9, 110, 637, 195, 226, 0.082),
    (30, 465.0, 121.5, 91.5, 778, 211, 358, 0.117),
    (40, 409.0, 94.05, 81.4, 895, 224, 475, 0.15),
    (50, 371.5, 77.22, 75.5, 998, 230, 579, 0.18),
    (70, 322.9, 57.50, 69.8, 1175, 231, 754, 0.24),
    (100, 279.8, 42.17, 67.7, 1396, 222, 961, 0.56),
    (150, 239.1, 29.72, 69.1, 1700, 209, 1210, 1.88),
    (200, 214.5, 23.21, 72.0, 1940, 202, 1390, 4.02),
    (300, 184.6, 16.39, 77.4, 2350, 195, 1630, 9.36),
    (450, 159.5, 11.57, 83.5, 2840, 192, 1790, 16.3),
    (700, 136.4, 7.903, 89.9, 3480, 191, 1890, 23.7),
    (1000, 120.4, 5.798, 94.5, 4090, 191, 1950, 29.8),
    (1500, 104.7, 4.064, 99.0, 4890, 191, 2060, 37.3),
    (2000, 94.85, 3.152, 102, 5530, 191, 2200, 43.4),
    (3000, 83.67, 2.106, 105, 6560, 191, 2460, 53.5),
    (4000, 75.06, 1.696, 107, 7370, 191, 2690, 62.0),
    (6000, 65.64, 1.177, 110, 8650, 192, 3060, 76.3),
    (8000, 59.78, 0.9091, 111, 9650, 192, 3350, 88.4),
    (10000, 55.68, 0.7460, 112, 10500, 192, 3570, 99.1),
    (15000, 49.17, 0.5270, 114, 12100, 192, 3970, 121),
    (20000, 45.26, 0.4180, 115, 13300, 192, 4230, 141),
    (30000, 40.67, 0.3100, 117, 15100, 192, 4530, 174),
    (40000, 37.92, 0.2570, 118, 16400, 192, 4700, 201),
    (60000, 34.12, 0.1936, 120, 18300, 192, 4880, 247),
    (80000, 31.00, 0.1510, 121, 19600, 192, 4960, 286),
    (100000, 28.40, 0.1206, 122, 20700, 192, 4940, 321),
    (150000, 23.91, 0.07897, 124, 22500, 192, 4660, 394),
    (200000, 21.22, 0.05719, 125, 23800, 192, 4300, 456),
    (300000, 18.11, 0.03659, 127, 25500, 192, 3690, 561),
    (400000, 16.29, 0.02583, 128, 26700, 192, 3260, 650),
    (600000, 14.10, 0.01438, 129, 28200, 192, 2710, 799),
    (800000, 12.76, 0.008623, 131, 29200, 192, 2390, 926),
    (1000000, 11.82, 0.005415, 132, 29900, 192, 2180, 1037),
    (1500000, 10.30, 0.002133, 133, 31100, 192, 1900, 1275),
    (2000000, 9.352, 0.001029, 135, 31800, 192, 1760, 1477),
    (3000000, 8.163, 0.0003522, 136, 32750, 192, 1670, 1816),
    (4000000, 7.414, 0.0001623, 138, 33300, 192, 1660, 2103),
    (6000000, 6.475, 0.0000514, 140, 34100, 192, 1750, 2587),
    (8000000, 5.882, 0.00002485, 141, 34500, 192, 1920, 2996),
    (10000000, 5.460, 0.00001360, 142, 34800, 192, 2140, 3357),
)


from HeWu.uc import _uc_ft2m, _uc_m2ft, _uc_psi2pa, _uc_pa2psi
from random import randint


def runABtest(airburst_to_op):
    """
    given a (lambda?) function to return the overpressure at time for

    """
    print(
        "{:^10}|{:^10},{:^10},{:^10}|{:^10}-{:^10}:{:^10}".format(
            "Y",
            "g.r.",
            "h.o.b",
            "sigma-tau",
            "Δp (calc)",
            "Δp (ref)",
            "Δ",
        )
    )
    print(
        "{:^10}|{:^10},{:^10},{:^10}|{:^10}-{:^10}:{:^10}".format(
            "kT",
            "m/kT^\u2153",
            "m/kT^\u2153",  # this bloody thing -> ⅓
            "ms/kT^\u2153",
            "Pa",
            "Pa",
            "1",
        )
    )
    print(
        "{:-^10}+{:-^10}-{:-^10}-{:-^10}+{:-^10}-{:-^10}-{:-^10}".format(
            "", "", "", "", "", "", ""
        )
    )
    for testPoint in abtests:
        sgr_ft, sbh_ft, sigma_tau, op_psi = testPoint
        op_pa = _uc_psi2pa(op_psi)
        Y = randint(1, 25000)  # 1kt to 25Mt
        Y3 = Y ** (1 / 3)

        sgr_m = _uc_ft2m(sgr_ft)
        sbh_m = _uc_ft2m(sbh_ft)

        gr_m = sgr_m * Y3
        hob_m = sbh_m * Y3

        time_elapsed = sigma_tau * Y3 / 1000  # ms -> s

        dp = airburst_to_op(gr_m, hob_m, Y, time_elapsed)
        delta = (dp - op_pa) / op_pa
        print(
            "{:^10}|{:^10.3g},{:^10.3g},{:^10.3g}|{:^10.3g}-{:^10.3g}:{:^15.1%}".format(
                Y, sgr_m, sbh_m, sigma_tau, dp, op_pa, delta
            )
        )


def runFAtest(freeair_from_r):
    """
    runs a test against the Table.2

    Arguments:
    - function freeair_from_r takes:
        - range in meter
        - yield in kiloton
        and should return:
        -peak overpressure in Pa
        -time of arrival in s
        -op positive phase in s
        -op impulse in pa-s
        -dp positive phase in s
        -dp impulse in pa-s
        -max. fireball impulse in C
    """
    Y = 1
    print(
        "{:^12} {:^12} {:^12} {:^12} {:^12} {:^12} {:^12} {:^12}".format(
            "ΔP_s - Pa",
            "r - m",
            "t.o.a - ms",
            "D_p^+ - ms",
            "I_p^+ - Pa-s",
            "D_u^+ - ms",
            "I_u^+ - Pa-s",
            "θ - deg C",
        )
    )
    print(
        "{:-^12}+{:-^12}+{:-^12}+{:-^12}+{:-^12}+{:-^12}+{:-^12}+{:-^12}".format(
            "", "", "", "", "", "", "", ""
        )
    )
    for testPoint in fatest:

        m = Y ** (1 / 3)

        (
            dP_s_ref,
            r_ref,
            sT_ref,
            sD_p_ref,
            sI_p_ref,
            sD_u_ref,
            sI_u_ref,
            theta_ref,
        ) = testPoint

        dP_s_ref = _uc_psi2pa(dP_s_ref)

        r = _uc_ft2m(r_ref * m)

        T_ref = sT_ref * m / 1000

        D_p_ref = sD_p_ref * m / 1000
        I_p_ref = _uc_psi2pa(sI_p_ref * m / 1000)

        D_u_ref = sD_u_ref * m / 1000
        I_u_ref = _uc_psi2pa(sI_u_ref * m / 1000)

        theta_ref *= 1000

        dP_s, T, D_p, I_p, D_u, I_u, theta = freeair_from_r(r, Y)
        # units expected: pa, s, s, pa-s, s, pa-s, deg C

        delta1 = (dP_s - dP_s_ref) / dP_s_ref
        delta2 = (T - T_ref) / T_ref
        delta3 = (D_p - D_p_ref) / D_p_ref
        delta4 = (I_p - I_p_ref) / I_p_ref
        delta5 = (D_u - D_u_ref) / D_u_ref
        delta6 = (I_u - I_u_ref) / I_u_ref
        delta7 = (theta - theta_ref) / theta_ref
        print(
            "{:^12.4g} {:^12.4g} {:^12.4g} {:^12.4g} {:^12.4g} {:^12.4g} {:^12.4g} {:^12.4g}".format(
                dP_s,
                r,
                T,
                D_p,
                I_p,
                D_u,
                I_u,
                theta,
            )
        )
        print(
            "{:^12.1%} {:^12} {:^12.1%} {:^12.1%} {:^12.1%} {:^12.1%} {:^12.1%} {:^12.1%}".format(
                delta1,
                "",
                delta2,
                delta3,
                delta4,
                delta5,
                delta6,
                delta7,
            )
        )
        print(
            "{:-^12}+{:-^12}+{:-^12}+{:-^12}+{:-^12}+{:-^12}+{:-^12}+{:-^12}".format(
                "", "", "", "", "", "", "", ""
            )
        )
