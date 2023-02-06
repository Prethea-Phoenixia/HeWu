"""
Zhai Jinpeng, 翟锦鹏, 2023
Contact: 914962409@qq.com

implements the code in Appendix B from the book:
PSR Report 1419-3
"AIRBLAST FROM NUCLEAR BURSTS—ANALYTIC APPROXIMATIONS"
Harold L. Brode, July 1987
Technical Report
CONTRACT No. DNA 001-85-C-0089  

"""


def PPEAK(X, Y, CAPR, Z, DELTPS=None):
    # THIS IS A FORTRAN IMPLEMENTATION OF THE ANALYTIC EXPRESSION FOR
    # PEAK OVERPRESSURE (BRODE AND SPEICHER, MAY 1986)
    #
    # THE PARAMETERS ARE:
    #       X = SCALED GR          (FT/KT**(1/3))
    #       Y = SCALED HOB         (FT/KT**(1/3))
    #    CAPR = (X*X+Y*Y)**(1/2)   (FT/KT**(1/3))
    #       Z = Y/X
    # DELTAPS = PEAK PRESSURE (PSI)

    R = CAPR / 1000
    A = 1.22 - (3.908 * Z * Z) / (1 + 810.2 * Z**5)
    B = (
        2.321
        + (6.195 * (Z**18) / (1 + 1.113 * (Z**18)))
        - (0.03831 * (Z**17) / (1 + 0.02415 * (Z**17)))
        + (0.6692 / (1 + 4164 * (Z**8)))
    )
    BB = (
        0.0629
        * ((X / Y) ** 8.34)
        / (1 + 0.00509 * ((X / Y) ** 13.05))
        * 0.05
        * Y
        / (1 + 2.56e-8 * (Y**5))
    )

    C = (
        4.153
        - (1.149 * (Z**18)) / (1 + 1.641 * (Z**18))
        - (1.1 / (1 + 2.771 * (Z**2.5)))
    )

    D = (
        (-4.166)
        + 25.76 * (Z**1.75) / (1 + 1.382 * (Z**18))
        + 8.257 * Z / (1 + 3.219 * Z)
    )

    E = 1 - 0.004642 * (Z**18) / (1 + 0.003886 * (Z**18))

    F = (
        0.6096
        + 2.879 * (Z**9.25) / (1 + 2.359 * (Z**14.5))
        - 17.15 * Z * Z / (1 + 71.66 * (Z**3))
    )

    G = 1.83 + 5.361 * Z * Z / (1 + 0.3139 * (Z**6))
    H = (
        -(0.2905 + 64.67 * (Z**5)) / (1 + 441.5 * (Z**5))
        - 1.389 * Z / (1 + 49.03 * (Z**5))
        + 8.808 * (Z**1.5) / (1 + 154.5 * (Z**3.5))
        + 1.094
        * (CAPR**2)
        / (
            (0.7813e9 - 1.234e5 * CAPR + 1201 * (CAPR**1.5) + (CAPR**2))
            * (1 + 2 * Y)
        )
    )
    P = 1.8008e-7 * (Y**4) / (1 + 0.0002863 * (Y**4)) - 2.121 * Y * Y / (
        794300 + (Y**4.3)
    )
    Q = 5.18 + 8.864 * (Y**3.5) / (3.788e6 + (Y**4))
    DELTPS = (
        (10.47) / (R**A)
        + (B - BB) / (R**C)
        + (D * E) / (1 + F * (R**G))
        + H
        + P / (R**Q)
    )
    return DELTPS


def PT(Y, X, SIGMA, DELTAP=None):
    # THIS IS A FORTRAN IMPLEMENTATION OF THE ANALYTIC EXPRESSION FOR
    # PRESSURE TIME HISTORY (BRODE AND SPEICHER, MAY 1986)
    # THE PARAMETERS ARE:
    # X      = SCALED GROUND RANGE   (FT/KT**(1/3))
    # Y      = SCALED HOB            (FT/KT**(1/3))
    # SIGMA  = SCALED TIME           (MSEC/KT**(1/3))
    #        ** NOTE THAT SIGMA IS THE TIME AFTER BURST IN THE           **
    #        ** ANALYTIC EXPRESSION. FOR PURPOSE OF CALCULATION          **
    #        ** SIGMA IN THIS ROUTINE IS TIME AFTER TIME OF ARRIVAL,     **
    #        ** WHICH IS THEN ADDED TO THE CALCULATED TIME OF ARRIVAL    **
    # DELTAP = PRESSURE              (PSI)
    #
    # IF PRESSURE IS DESIRED IN KPA THEN MULTIPLY RESULT BY 100/14.504
    #
    # NOTE THAT LIMITS ARE PLACED O VALUES SUCH AS X, Y, Z ETC...
    # THIS IS DONE TO AVOID OVERFLOWS AND THE VALUES ARE MACHINE DEPENDENT.

    X = max(X, 1e-9)
    Y = max(Y, 1e-9)
    CAPR = (X * X + Y * Y) ** 0.5
    R = CAPR / 1000
    Z = Y / X

    Z = min(Z, 100)

    DELTPS = PPEAK(X, Y, CAPR, Z)
    XM = 170 * Y / (1.0 + 60.0 * (Y**0.25)) + 2.89 * ((Y / 100) ** 2.5)
    U = (
        (0.543 - 21.8 * R + 386 * (R**2) + 2383 * (R**3))
        * (R**8)
        / (
            2.99e-14
            - 1.91e-10 * (R**2)
            + 1.032e-6 * (R**4)
            - 4.43e-6 * (R**6)
            + (1.028 + 2.087 * R + 2.69 * (R**2)) * (R**8)
        )
    )

    if X < XM:
        TAU = U
    else:
        W = (
            (1.086 - 34.605 * R + 486.3 * (R**2) + 2383 * (R**3))
            * (R**8)
            / (
                3.0137e-13
                - 1.2128e-9 * (R**2)
                + 4.128e-6 * (R**4)
                - 1.116e-5 * (R**6)
                + (1.632 + 2.629 * R + 2.69 * (R**2)) * (R**8)
            )
        )
        TAU = U * XM / X + W * (1 - XM / X)

    SIGMA = SIGMA + TAU

    S2 = (
        1
        - 15.18 * ((Y / 100) ** 3.5) / (1 + 15.18 * ((Y / 100) ** 3.5))
        - (0.02441 * ((Y / 1.0e6) ** 2) / (1 + 9000 * ((Y / 100)) ** 7))
        * (1.0e10 / (0.441 + ((X / 100) ** 10)))
    )
    CAPD = (
        (1640700 + 24629 * TAU + 416.15 * (TAU**2))
        / (10880 + 619.76 * TAU + (TAU**2))
    ) * (
        0.4
        + 0.001204 * (TAU**1.5) / (1 + 0.001559 * (TAU**1.5))
        + (
            0.6126
            + 0.5486 * (TAU**0.25) / (1 + 0.00357 * (TAU**1.5))
            - 3.47 * (TAU**0.637) / (1 + 5.696 * (TAU**0.645))
        )
        * S2
    )
    S = (
        1
        - 1100 * ((Y / 100) ** 7) / (1 + 1100 * ((Y / 100) ** 7))
        - (2.441e-14 * Y * Y / (1 + 9000 * ((Y / 100) ** 7)))
        * (1.0e10 / (0.441 + ((X / 100) ** 10)))
    )
    F2 = (
        (
            0.445
            - 5.44 * (R**1.02) / (1 + 100000 * (R**5.84))
            + 7.571
            * (Z**7.15)
            / (
                1 + 5.135 * (Z**12.9)
            )  # change from code: - to + before coefficient 5.135
            - 8.07 * (Z**7.31) / (1 + 5.583 * (Z**12.23))
        )
        * (0.435 * ((Y / 10) ** 1.26) / (1 + 0.03096 * ((Y / 10) ** 3.12)))
        * (1 - 0.000019 * (TAU**8) / (1 + 0.000019 * (TAU**8)))
    )
    F = (
        (
            0.01477 * (TAU**0.75) / (1 + 0.005836 * TAU)
            + 7.402e-5 * (TAU**2.5) / (1 + 1.429e-8 * (TAU**4.75))
            - 0.216
        )
        * S
        + 0.7076
        - 3.077e-5 * (TAU**3) / (1 + 4.367e-5 * (TAU**3))
        + F2
        - (0.452 - 9.94e-7 * (X**4.13) / (1 + 2.1868e-6 * (X**4.13)))
        * (1 - 1.5397e-4 * (Y**4.3) / (1 + 1.5397e-4 * (Y**4.3)))
    )

    G = 10 + (77.58 - 64.99 * (TAU**0.125) / (1 + 0.04348 * (TAU**5))) * S
    H = (
        3.003
        + 0.05601 * TAU / (1 + 1.473e-9 * (TAU**5))
        + (
            0.01769 * TAU / (1 + 3.207e-10 * (TAU**4.25))
            - 0.03209 * (TAU**1.25) / (1 + 9.914e-8 * (TAU**4))
            - 1.6
        )
        * S
        - 0.1966 * (TAU**1.22) / (1 + 0.767 * (TAU**1.22))
    )
    B = (F * ((TAU / SIGMA) ** G) + (1 - F) * ((TAU / SIGMA) ** H)) * (
        1 - (SIGMA - TAU) / CAPD
    )
    if X < XM or Y > 380:
        DELTAP = DELTPS * B
        # print(B)
    else:
        XE = 3.039 * Y / (1 + 0.0067 * Y)
        AK = abs((X - XM) / (XE - XM))
        AK = min(AK, 50)
        D2 = 2.99 + 31240 * ((Y / 100) ** 9.86) / (1 + 15530 * ((Y / 100) ** 9.87))
        D = (
            0.23
            + 0.583 * Y * Y / (26667 + Y * Y)
            + 0.27 * AK
            + (0.5 - 0.583 * Y * Y / (26667 + Y * Y)) * (AK**D2)
        )
        A = (D - 1) * (1 - (AK**20) / (1 + AK**20))
        AJ = 11860 * (SIGMA - TAU) / (Y * ((X - XM) ** 1.25))
        AJ = min(AJ, 200)
        V = 1 + (
            0.003744 * ((Y / 10) ** 5.185) / (1 + 0.004684 * ((Y / 10) ** 4.189))
            + 0.004755 * ((Y / 10) ** 8.049) / (1 + 0.003444 * ((Y / 10) ** 7.497))
            - 0.04852 * ((Y / 10) ** 3.423) / (1 + 0.03038 * ((Y / 10) ** 2.538))
        ) * (AJ**3) / (6.13 + (AJ**3)) * (1 / (1 + 9.23 * (AK**2)))

        C3 = 1 + (
            1.094
            * (AK**0.738)
            / (1 + 3.687 * (AK**2.63))
            * (
                1
                - 83.01 * ((Y / 100) ** 6.5) / (1 + 172.3 * ((Y / 100) ** 6.04))
                - 0.15
            )
        ) * (1 / (1 + 0.5089 * (AK**13)))
        C2 = 23000 * ((Y / 100) ** 9) / (1 + 23000 * ((Y / 100) ** 9))
        TEMP = (X / 100) ** 4
        C = (
            (1.04 - 0.02409 * TEMP / (1 + 0.02317 * TEMP))
            * (AJ**7)
            / ((1 + A) * (1 + 0.923 * (AJ**8.5)))
            * (C2 + (1 - C2) * (1 - 0.09 * (AK**2.5) / (1 + 0.09 * (AK**2.5))))
            * C3
            * (1 - (((SIGMA - TAU) / CAPD) ** 8))
        )

        DELTAP = DELTPS * (1 + A) * (B * V + C)
        print(B, F)
        print(S, F2, F)
    return DELTAP


if __name__ == "__main__":
    from uc import _uc_m2ft, _uc_psi2pa

    def _airburst_to_op(gr, h, Y, t):
        m = Y ** (1 / 3)
        return _uc_psi2pa(PT(_uc_m2ft(h / m), _uc_m2ft(gr / m), t * 1000 / m))

    from HeWu.test import runABtest

    runABtest(_airburst_to_op)
