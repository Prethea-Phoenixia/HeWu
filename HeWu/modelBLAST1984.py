"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

a faithful (with correction) reproduction of the algorithm used
in BLAST.EXE from Horizons Technology for the Defense Nuclear Agency
dating to 1984. 
"""


from math import sqrt, log, exp, atan, pi, sin


def clamp(x, a, b):

    llim = min(a, b)
    hlim = max(a, b)

    return max(min(x, hlim), llim)


def _dPdna(x):
    return (
        3.04e11 / x**3
        + 1.13e9 / x**2
        + 7.9e6 / (x * sqrt(log(x / 445.42 + 3 * exp(-1 / 3 * sqrt(x / 445.42)))))
    )


def _n(x):
    xi = x / 101325 + 1
    t = 1e-12 * (xi) ** 6
    z = log(xi) - 0.47 * t / (100 + t)
    gs = 1.402 - 3.4e-4 * z**4 / (1 + 2.22e-5 * z**6)

    mus = (gs + 1) / (gs - 1)
    n = (1 + mus * xi) / (5.975 + xi)
    return n, gs


def _ta(x):
    return x**2 * (6.7 + x) / (7.12e6 + 7.32e4 * x + 340.5 * x**2)


def freeair(Y, ALT, RANGE):

    Y3 = Y ** (1 / 3)

    if ALT < 0:
        raise ValueError(
            "free air model is not applicable in case of a underground burst"
        )
    if 0 <= ALT < 11000:
        T = 1 - ALT / sqrt(2e9)
        P = T**5.3

    elif ALT <= 20000:
        T = 0.7537 * (1 + 2.09e-7 * ALT)
        P = sqrt(1.6) * (1 + 2.09e-7 * ALT) ** (-754)

    else:
        T = 0.684 * (1 + 5.16e-6)
        P = 1.4762 * (1 + 5.16e-6) * (-33.6)

    SP = P
    SD = SP ** (-1 / 3)
    ST = SD / sqrt(T)
    C = 340.5 * SD / ST

    R = RANGE / (SD * Y3)

    dPfree = _dPdna(R)

    PFREE = dPfree * SP

    n, _ = _n(dPfree)

    qfree = 0.5 * dPfree * (n - 1)

    QFREE = qfree * SP

    tafree = _ta(R)

    TAFREE = tafree * ST * Y3

    withinLimit = True

    if ALT < 0 or ALT > 32000:
        withinLimit = False

    if Y < 0.1 or Y > 25000:
        withinLimit = False

    if RANGE < 16 * Y3 * SD or RANGE > 4000 * Y3 * SD:
        withinLimit = False

    return PFREE, QFREE, TAFREE, withinLimit


def airburst(Y, HOB, GR):
    """
    inputs:
    Y: yield, kt
    ALT: burst height, m
    RANGE: ground range, m

    return:
    PAIR: peak overpressure, pa
    QAIR: peak dynamic pressure, pa
    TAAIR: time of arrival, s
    IPTOTAL: overpressure impulse, pa-s
    DPP: over pressure positive phase duration, s
    limit1: whether model limit is exceeded for PAIR to DPP

    XM: mach stem formation range, m
    HTP: height of triple point, m
    limit2: whether XM and HTP estimates are within model limit

    IQTOTAL: total dynamic pressure impulse, pa-s
    DPQ: dynamic pressure positive phase duration, s
    limit3: whether the IQTOTAL and DPQ estimations are within limit.
    (this model is only applicable in the mach reflection region for these
    parameters)
    """
    if HOB < 0:
        raise ValueError("model is not applicable to underground bursts")

    N = 100

    Y3 = Y ** (1 / 3)

    SGR = GR / Y3
    SHOB = HOB / Y3
    SR = sqrt(SGR**2 + SHOB**2)

    alpha = atan(SHOB / SGR)

    dPfree = _dPdna(SR)

    T = 340 / dPfree**0.55
    U = 1 / (7782 / dPfree**0.7 + 0.9)
    W = 1 / (7473 / dPfree**0.5 + 6.6)
    V = 1 / (647 / dPfree**0.8 + W)

    alphaMach = atan(1 / (T + U))
    beta = atan(1 / (T + V))

    s = (alpha - alphaMach) / beta

    so = clamp(s, 1, -1)

    sigma = 0.5 * (sin(pi * so * 0.5) + 1)

    inMach, inReg = True, True

    dPreg, dPmach = 0, 0

    if sigma == 0:
        inReg = False

    elif sigma == 1:
        inMach = False

    if inMach:
        A = min(3.7 - 0.94 * log(SGR), 0.7)
        B = 0.77 * log(SGR) - 3.8 - 18 / SGR
        C = max(A, B)

        dPmach += _dPdna((SGR / 2 ** (1 / 3))) / (1 - C * sin(alpha))

    if inReg:
        n, gs = _n(dPfree)
        Rn = 2 + 0.5 * (gs + 1) * (n - 1)
        f = dPfree / 75842

        D = f**6 * (1.2 + 0.07 * sqrt(f)) / (f**6 + 1)

        dPreg += dPfree * ((Rn - 2) * sin(alpha) ** D + 2)

    PAIR = dPreg * sigma + dPmach * (1 - sigma)

    nq, _ = _n(PAIR)
    QAIR = 0.5 * PAIR * (nq - 1) * (1 - sigma * sin(alpha) ** 2)

    xm = SHOB**2.5 / 5822 + 2.09 * SHOB**0.75
    XM = xm * Y3

    S = 1 / (5.98e-5 * SHOB**2 + 3.8e-3 * SHOB + 0.766)
    h = 0.9 * xm - 3.6 * SHOB
    try:
        HTP = S * (h + sqrt(h**2 + (SGR - 0.9 * xm) ** 2 - xm**2 / 100)) * Y3
    except ValueError:
        HTP = None

    if SGR <= xm:
        v = 1
    else:
        v = 1.26 - 0.26 * (xm / SGR)

    R = SR / v

    taair = _ta(R)

    TAAIR = taair * Y3 * v

    """total pressure impulse"""

    SGR = max(SGR, 1e-7)
    SHOB = max(SHOB, 1e-7)
    SR = sqrt(SGR**2 + SHOB**2)

    s = (
        1
        - 1 / (1 + 1 / (4.5e-8 * SHOB**7))
        - (5.958e-3 * SHOB**2)
        / (1 + 3.682e-7 * SHOB**7)
        / (1 + SGR**10 / 3.052e14)
    )

    f = (
        s
        * (
            2.627 * taair**0.75 / (1 + 5.836 * taair)
            + 2341 * taair**2.5 / (1 + 2.541e6 * taair**4.75)
            - 0.216
        )
        + 0.7076
        - 3.077 / (1e-4 * taair ** (-3) + 4.367)
    )

    g = 10 + s * (77.58 - 154 * taair**0.125 / (1 + 1.375 * taair**0.5))

    h = (
        s
        * (
            17.69 * taair / (1 + 1803 * taair**4.25)
            - 180.5 * taair**1.25 / (1 + 99140 * taair**4)
            - 1.6
        )
        + 2.753
        + 56 * taair / (1 + 1.473e6 * taair**5)
    )

    to = log(1000 * taair) / 3.77

    dpsurf = 1e-3 * (155 * exp(-20.8 * taair) + exp(-(to**2) + 4.86 * to + 0.25))

    dpunmod = dpsurf * (
        1
        - (1 - 1 / (1 + 4.5e-8 * SHOB**7))
        * (0.04 + 0.61 / (1 + taair**1.5 / 0.027))
    )

    dpDp = dpunmod * (1.16 * exp(-abs(SHOB / 0.3048 - 156) / 1062))

    DPP = dpDp * (Y3)

    singlePeak = False
    if SGR < xm or SHOB > 116:  # singlepeak
        singlePeak = True

    if singlePeak:
        pass
    else:
        xe = 138.3 / (1 + 45.5 / SHOB)
        e = clamp(abs((SGR - xm) / (xe - SGR)), 50, 0.02)
        w = 0.583 / (1 + 2477 / SHOB**2)

        d = 0.23 + w + 0.27 * e + e**5 * (0.5 - w)
        a = (d - 1) * (1 - 1 / (1 + e ** (-20)))

        dt = max(SHOB * (SGR - xm) ** 1.25 / 8.186e5, 1e-12)

        vo = SHOB**6 / (2445 * (1 + SHOB**6.75 / 3.9e4) * (1 + 9.23 * e**2))
        co = (1.04 - 1.04 / (1 + 3.725e7 / SGR**4)) / (
            (a + 1) * (1 + 9.872e8 / SHOB**9)
        )

    dp = dpDp

    accumulator = 0

    if singlePeak:
        for i in range(N):
            # simple newton-raphson quadrature using midpoint rule
            t = taair + dp * (i + 0.5) / N
            b = (f * (taair / t) ** g + (1 - f) * (taair / t) ** h) * (
                1 - (t - taair) / dp
            )

            accumulator += b  # dpt
    else:
        for i in range(N):
            t = taair + dp * (i + 0.5) / N
            b = (f * (taair / t) ** g + (1 - f) * (taair / t) ** h) * (
                1 - (t - taair) / dp
            )
            ga = clamp((t - taair) / dt, 400, 0.0001)

            v = 1 + vo * ga**3 / (ga**3 + 6.13)
            c = co / (ga ** (-7) + 0.923 * ga**1.5) * (1 - ((t - taair) / dp) ** 8)

            accumulator += (1 + a) * (b * v + c)

    IPTOTAL = Y3 * PAIR * accumulator * dp / N

    SHOBo = SHOB / 0.3048
    SGRo = SGR / 0.3048

    SHOBx = abs(SHOBo - 200) + 200
    SGRx = SGRo - 1000

    dpo = 0.3 + 0.42 * exp(-SHOBx / 131)

    if SGRx > 0:
        dpx = dpo + 4.4e-5 * SGRx
    else:
        dpx = dpo + SGRx * (1 / 2361 - (SHOBx - 533) ** 2 / 7.88e7)

    if SHOBo >= 200:
        dpq = dpx
    else:
        dpq = dpx * (1 + 0.2 * sin(pi / 200 * SHOBo))

    DPQ = dpq * Y3

    deltao = max(SHOBo**1.52 / 16330 - 0.29, 0)

    delta = 2.38 * exp(-7e-7 * abs(SHOBo - 750) ** 2.7 - 4e-7 * SGRo**2) + deltao
    qo = 0.5 * (1 - sigma * sin(alpha) ** 2)

    dp = dpq

    accumulator = 0

    if singlePeak:
        for i in range(N):
            # simple newton-raphson quadrature using midpoint rule
            t = taair + dp * (i + 0.5) / N
            b = (f * (taair / t) ** g + (1 - f) * (taair / t) ** h) * (
                1 - (t - taair) / dp
            )
            dpt = b * PAIR  # dpt
            nq, _ = _n(dpt)
            qt = 0.5 * dpt * (nq - 1) * (dpt / PAIR) ** delta
            accumulator += qt

    else:
        for i in range(N):
            t = taair + dp * (i + 0.5) / N
            b = (f * (taair / t) ** g + (1 - f) * (taair / t) ** h) * (
                1 - (t - taair) / dp
            )
            ga = clamp((t - taair) / dt, 400, 0.0001)
            # ga = (t - taair) / dt
            v = 1 + vo * ga**3 / (ga**3 + 6.13)
            c = co / (ga ** (-7) + 0.923 * ga**1.5) * (1 - ((t - taair) / dp) ** 8)
            dpt = (1 + a) * (b * v + c) * PAIR
            nq, _ = _n(dpt)
            qt = 0.5 * dpt * (nq - 1) * (dpt / PAIR) ** delta
            accumulator += qt

    IQTOTAL = Y3 * accumulator * dp / N

    """overpressure, dynamic pressure, time of arrival, overpressure impulse"""
    limit1 = True

    if Y > 25e3 or Y < 0.1:
        limit1 = False

    if HOB < 0 or HOB > 4000 * Y3:
        limit1 = False

    if HOB >= 25 * Y3:
        LM = 0
    else:
        LM = 20 * Y3

    if GR < LM or GR > 4000 * Y3:
        limit1 = False

    """mach stem and triple point height"""
    limit2 = True

    if Y > 25e3 or Y < 0.1:
        limit2 = False

    if HOB < 0 or HOB > 800 * Y3:
        limit2 = False

    LM = max(XM, 20 * Y3)

    if GR < LM or GR > 4000 * Y3:
        limit2 = False

    """ dynamic pressure total impulse"""
    limit3 = True

    if Y > 25e3 or Y < 0.1:
        limit3 = False

    if HOB < 0 or HOB > 750 * Y3:
        limit3 = False

    LM = max(1.3 * XM, 80 * Y3)

    if GR < LM or GR > 4000 * Y3:
        limit3 = False

    """ sanitize the output to remove misleading elements"""
    if GR < XM:
        HTP = None  # triple point cannot exist inside of mach stem formation range
        IQTOTAL = None  # there is still dynamic pressure yet we cannot model it.
        DPQ = None

    return (
        PAIR,
        QAIR,
        TAAIR,
        IPTOTAL,
        DPP,
        limit1,
        XM,
        HTP,
        limit2,
        IQTOTAL,
        DPQ,
        limit3,
    )


if __name__ == "__main__":
    print(*airburst(40, 738.3, 446.4), sep="\n")

    print("----------------------------")
