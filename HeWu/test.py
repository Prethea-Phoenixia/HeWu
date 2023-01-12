"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

Implements the test cases found in the Appendix A of:
PSR Report 1419-3
"AIRBLAST FROM NUCLEAR BURSTS—ANALYTIC APPROXIMATIONS"
Harold L. Brode, July 1987
Technical Report
CONTRACT No. DNA 001-85-C-0089  

"""


"""
test cases, in the format of:
    X: scaled ground range, ft/kT^(1/3)
    Y: scaled burst height, ft/kT^(1/3)
    sigma-tau: scaled time - (scaled) time of arrival, ms/kT^(1/3) 
    DeltaP: overpressure, psi.
"""

tests = (
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


from HeWu.uc import _uc_ft2m, _uc_m2ft, _uc_psi2pa, _uc_pa2psi
from random import randint


def runtest(airburst_to_op):
    """
    given a (lambda?) function to return the overpressure at time for

    """
    print(
        "{:^10}|{:^10},{:^10},{:^10}|{:^10}-{:^10}:{:^10}".format(
            "Y",
            "g.r.",
            "h.o.b",
            "sigma-tau",
            "Δp (ref)",
            "Δp (calc)",
            "Δ",
        )
    )
    print(
        "{:^10}|{:^10},{:^10},{:^10}|{:^10}-{:^10}:{:^10}".format(
            "kT",
            "m",
            "m",
            "ms",
            "kPa",
            "kPa",
            "1",
        )
    )
    print(
        "{:-^10}+{:-^10}-{:-^10}-{:-^10}+{:-^10}-{:-^10}-{:-^10}".format(
            "", "", "", "", "", "", ""
        )
    )
    for testPoint in tests:
        sgr, sbh, sigma_tau, op_psi = testPoint
        op_pa = _uc_psi2pa(op_psi)
        Y = randint(1, 25000)  # 1kt to 10Mt
        Y3 = Y ** (1 / 3)
        gr_ft = sgr * Y3
        hob_ft = sbh * Y3
        gr_m = _uc_ft2m(gr_ft)
        hob_m = _uc_ft2m(hob_ft)
        time_elapsed = sigma_tau * Y3 / 1000  # ms -> s

        dp = airburst_to_op(gr_m, hob_m, Y, time_elapsed)
        delta = (dp - op_pa) / op_pa
        print(
            "{:^10d}|{:^10.3g},{:^10.3g},{:^10.3g}|{:^10.3g}-{:^10.3g}:{:^15.1%}".format(
                Y, gr_m, hob_m, time_elapsed * 1000, op_pa, dp, delta
            )
        )
