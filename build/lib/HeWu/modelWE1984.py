"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

a faithful (with correction) reproduction of the algorithm used
in WE.EXE from Horizons Technology for the Defense Nuclear Agency
dating to 1984. 

"""


from math import sqrt, log, log10, exp, pi, sin


def clamp(x, a, b):

    llim = min(a, b)
    hlim = max(a, b)

    return max(min(x, hlim), llim)


def sign(x):
    if x > 0:
        return 1
    elif x == 0:
        return 0
    else:
        return -1


def phi(z):
    # probability density function of the standard normal distribution

    return exp(-(z**2) / 2) / sqrt(2 * pi)


def gamma(u):
    return 0.994 - 0.446 * (u - 1) + 0.455 * (u - 1) ** 2


cDn = {
    1: (253, 71, 138),
    2: (98, 49, 225),
    3: (253, 71, 138),
    4: (394, 31, 261),
    5: (347, 52, 147),
    6: (177, 67, 152),
    7: (394, 31, 261),
    8: (450, 24, 323),
    9: (753, -13, 492),
    10: (272, 13, 933),
    11: (394, 31, 261),
    12: (300, -13, 492),
    13: (1431, -2, 63),
}

cDg = {
    1: (357, 32, -188),
    2: (605, 0, 793),
    3: (433, 27, 68.8),
    4: (542, 14, 439),
    5: (542, 14, 439),
    6: (433, 27, 68.8),
    7: (542, 14, 439),
    8: (490, 18, 297),
    9: (542, 14, 439),
    11: (542, 14, 439),
    12: (604, 7, 633),
    13: (722, 8, 651),
}

cCg = {
    1: (59, 9, 8, 5, 98, 108),
    2: (61, 9, 34, 13, 108, 142),
    3: (59, 9, 8, 5, 98, 108),
    4: (59, 9, 8, 5, 98, 108),
    5: (50, 11, -25, 1, 90, 92),
    6: (61, 9, 34, 13, 108, 142),
    7: (42, 12, -4, 2, 100, 93),
    8: (44, 11, -43, 0, 86, 87),
    9: (44, 11, -43, 0, 86, 87),
    10: (61, 9, 34, 13, 108, 142),
    11: (44, 11, -43, 0, 86, 87),
    12: (42, 12, -4, 2, 100, 93),
    13: (-87, 12, 0, 0, 62, 54),
}
""" last one is 54 in the original implementation. but 84 in accompanying documentation"""


def iniRad(Y, AIR, H, GR, FF, WT):
    """
    rad = exp(-(H**1.04) / 13800)
    AIR = (AIR + rad) / 2
    """

    if H < 0:
        raise ValueError(
            "underground burst is not yet considered for initial radiation model"
        )

    if isinstance(WT, int) and 13 >= WT >= 1:
        pass
    else:
        raise ValueError("invalid weapon type")

    if FF < 0 or FF > 1:
        raise ValueError("impossible fission fragment")

    SR = sqrt(GR**2 + H**2)
    SRo = AIR * SR
    Ho = AIR * H
    s = clamp((Ho - 277) / 50, 1, -1)
    sigma = 0.5 * (1 + sin(s * pi * 0.5))

    """Dn"""

    ap, bp, cp = cDn[WT]

    a = 1e6 * ap
    b = -(500 + bp) / 1e5
    c = 1 + cp / 1000
    Dn = a / (SRo) ** c * exp(b * SRo)

    """Dg"""

    if WT == 10:
        a = 13.5
        b = -0.344
        c = -1.537
        d = 0.5173

    else:
        ap, bp, cp = cDg[WT]
        a = 10 ** (ap / 100)
        b = -(193 + bp) / 1e4
        c = cp / 1000
        d = 0.8

    Dg = a / (SRo) ** c * exp(b * SRo**d)

    """Dff"""

    p = AIR * (0.264 - AIR / 12.6) + 0.815
    SRx = sign(SR - 140) * abs(SR - 140) ** p + 140
    Dff = 1.42e7 / SRx**0.9516 * exp(-(SRx**0.774) / 32.7)

    """Cg"""

    ap, app, bp, bpp, cp, cpp = cCg[WT]

    Hx = min(1000, Ho)

    if WT == 13:
        Cgmax = 1.1
        ao = 0.002
        b = -2.12 + exp(Hx**0.903 / 204) + 0.977 * (1 - sigma)

    else:
        Cgmax = 1
        ao = 0.0035
        b = bp / 100 + exp(Hx**0.88 / (bpp + 192)) + (1 - sigma) * pi / 6

    x = 0.9 - ap / 1000
    t1 = ao * Ho**x
    t2 = ao * 277**x + 0.00011 * Ho ** (1 + app / 100)
    a = min(Cgmax, 0.31 + (1 - sigma) * t1 + sigma * t2)
    c = cp + sigma / cpp * max(0, Hx - 277) ** 1.4

    Cg = a + exp(b + c * (1 - exp(4e-5 * SRo)))

    """Cf """

    """ issue: is log e or log 10 ?
    log10 is numerically closer"""
    SH = min(H / Y ** (1 / 3), 250)
    sf = sin(1.16 * log(Y, 10) - 1.39)
    z = (log(Y, 10) - 1.4) / 3.5

    Cf = (
        0.5 * (z + sqrt(z**2 + 0.044)) * AIR * (1 - SH / 125)
        + 0.038 * log(Y, 10)
        - 0.22
        + (GR / 1000 - 0.65 * log(Y, 10) - 0.4)
        * 0.075
        * (sf + 2 * abs(sf) * (AIR - 0.9))
    )

    Cf = 10**Cf

    """ He Hydrodynamic enhancement factor"""

    if Y < 1:
        He = 1
    else:
        a = 0.1455 * AIR - 0.0077
        b = 2.55 - 0.35 * AIR

        AH = 10 ** (a * log10(Y) ** b)

        c = 0.05875 * AIR + 0.004
        d = 0.04 * AIR - 0.03 * sign(AIR - 0.6) * abs(AIR - 0.6) ** 1.3

        BY = 1 - (1 + c * log10(Y) ** 2.6) * exp(-d * log10(Y) ** 2.6)

        He = min(AH * exp(BY * SR / 1000), exp(SR * (-0.26 + 2.563 * AIR) / 1000))

    """Cn"""

    t1 = 0.205 + 2.2e-3 * Ho**0.839
    """ issue: where does ^2 operate on ? resolved, (/637)^2"""
    t2 = 0.4514 + Ho**1.636 / 637**2
    a = min(1, (1 - sigma) * t1 + sigma * t2)
    t3 = 0.388 + 0.116 * Ho**0.27
    t4 = 0.9176 + Ho**3.726 / 1.78e10
    b = (1 - sigma) * t3 + sigma * t4
    """ issue: 0.0728 or 0.728? resovled: 0.0728"""
    c = 8e-4 + 0.0728 / (Ho**1.23 + 23.6)
    d = 0.9
    if Ho > 1:
        d += log(Ho) / 25

    Cn = a + b * exp(-c * SRo**d)

    """end"""
    N = Dn * Cn * Y * AIR**2
    SG = Dg * Cg * Y * AIR**2
    FFG = Dff * Cf * He * Y * FF
    TD = N + SG + FFG
    NtG = N / (SG + FFG)

    if WT == 10:
        fn = exp(SR / 800) / 250
    elif WT == 13:
        fn = exp(-SR / 234.7) / 20 + exp(-SR / 4329) / 25
    else:
        fn = 0.015

    DS = FFG + SG + fn * N

    """limits"""

    x1 = sqrt(max(0, 1e4 - H**2))
    x2 = sqrt(max(0, (150 * max(1, Y ** (1 / 3))) ** 2 - H**2))
    y = sqrt(1e8 - H**2)

    NSGLim = True
    OTHLim = True

    if Y < 0.01 or Y > 25000:
        NSGLim = False
        OTHLim = False

    if AIR < 0.6 or AIR > 1:
        NSGLim = False
        OTHLim = False

    if H < 1.5 / AIR or H > 10000 / AIR:
        NSGLim = False
        OTHLim = False

    if GR < x1 / AIR or GR > y / AIR:
        NSGLim = False
    if GR < x2 / AIR or GR > y / AIR:
        OTHLim = False

    return (N, SG, NSGLim, FFG, TD, DS, NtG, OTHLim)


def therm(Y, H, GR, VIS):
    if H < 0:
        raise ValueError(
            "underground and surface burst are not yet considered for thermal models"
        )
    SR = sqrt(H**2 + GR**2)
    HT = 4 * Y ** (1 / 3)
    A1 = 0.32 * (1 - exp(-12 * Y ** (-VIS / 17e3)))
    B1 = -log10(Y) ** 2 / 275 + 0.0186 * log10(Y) - 0.025
    A2 = ((30 * Y**-0.26) ** 4 + 1350) ** (-1 / 4)
    B2 = -(1.457 / VIS + 9.3e-6)

    FS = A1 * exp(B1 * SR) + A2 * exp(B2 * SR) + 0.006
    A3 = H ** (3 / 2) / 5e7 + 97 / (281 + sqrt(Y))

    if H == 0:
        B3 = -1.112 / VIS
    else:
        B3 = 0.139 / H * (exp(-8 * H / VIS) - 1)
    FA = A3 * exp(B3 * SR)

    if H >= HT:
        F = FA
    else:
        F = FA * (H / HT) + FS * (1 - H / HT)

    return 8e6 * F * Y / SR**2


Mn = {
    1: "dry soil",
    2: "wet soil",
    3: "dry soft rock",
    4: "wet soft rock",
    5: "hard rock",
}


def crater(Y, H, M1, T1, M2, T2, M3, GR, HIGHRAD=False):
    if isinstance(M1, int) and 0 < M1 < 6:
        pass
    else:
        raise ValueError("invalid material for layer 1")
    if isinstance(M2, int) and 0 < M2 < 6:
        pass
    else:
        raise ValueError("invalid material code for layer 2")
    if isinstance(M3, int) and 0 < M3 < 6:
        pass
    else:
        raise ValueError("invalid material code for layer 3")

    if T1 < 0:
        raise ValueError("invalid thickness for layer 1")
    if T2 < 0:
        raise ValueError("invalid thickness for layer 2")

    if H > 3 * Y ** (1 / 3):
        raise ValueError("cratering is insignificant")

    isWithinLimit = True

    S1 = -H / Y ** (1 / 3)
    if -3 <= S1 < 0.15:
        alpha1 = 1 / 3
    elif 0.15 <= S1 < 5:
        alpha1 = 0.2946 + exp(-S1 * log10(583)) / sqrt(305)
    elif S1 >= 5:
        alpha1 = 1 / 3.4

    def V(M):
        if M == 1:
            if -3 <= S1 < 0:
                alpha20 = 1 / 3.1
            elif 0 <= S1 < 5:
                alpha20 = 1 / (3.4 - 0.3 * exp(-2.2 * S1))
            elif S1 >= 5:
                alpha20 = 1 / 3.4
        else:
            alpha20 = 1 / 3.4

        if Y > 20:
            alpha = alpha20
        else:
            if HIGHRAD:
                g = 0
            else:
                g = 1 - clamp(log(Y, 20), 0, 1)

            alpha = alpha20 + g * (alpha1 - alpha20)

        S2 = -H / Y**alpha

        if M == 1:
            J = 4
            P, Q, R = 9.7, 0.103, -0.00143
        elif M == 2:
            J = 1
            P, Q, R = 12.54, 0.029, -0.00078
        elif M == 3:
            J = 5
            P, Q, R = 9.34, 0.131, -0.00231
        elif M == 4:
            J = 2
            P, Q, R = 10.45, 0.089, -0.00134
        elif M == 5:
            J = 8
            P, Q, R = 8.72, 0.1634, -1 / 370

        def SV(y):

            if S2 >= 5:
                if M == 2:
                    S = 503**2 * exp(-S2 / 30)
                else:
                    S = 0
                if S2 >= 40:
                    isWithinLimit = False

                SV = exp(P + Q * S2 + R * S2**2) - S
            else:
                if M == 1 and y >= 20:
                    if -3 <= S2 < 0:
                        SV = 354 * 10 ** (
                            0.506 * (exp(2.6 * S2 + 0.486 * S2**2) - 1) + 2 * S2 / 9
                        )

                    elif 0 <= S2 < 5:
                        SV = 354 * exp(
                            (1 - exp(-3.967 * S2**1.139))
                            * (4.283 - 0.0515 * (5 - S2) ** 2.068)
                        )
                    elif S2 < -3:
                        raise ValueError("cratering is insignificant")
                else:
                    if y <= 1 and H <= 0:
                        F, G, D, K, L = -1.05, -0.105, 0.0573, -0.5, 16989
                    elif y <= 1 and H > 0:
                        F, G, D, K, L = 0.258, 0.01, 0.1, 1.9, 16989
                    elif y >= 20 and H <= 0:
                        F, G, D, K, L = -2, -0.3044, 0.0707, -0.9059, 5663
                    elif y >= 20 and H > 0:
                        F, G, D, K, L = 0.53, 0.028, -1 / 46, 1.74, 5663

                    SV = L / J * 10 ** (K * (exp(F * S2 + G * S2**2) - 1) + D * S2)

            return SV

        if Y > 20:
            return SV(Y) * Y ** (3 * alpha)
        # or Y < 1
        else:
            return SV(20) * (SV(1) / SV(20)) ** g * Y ** (3 * alpha)

    VM1, VM2, VM3 = V(M1), V(M2), V(M3)

    def hat(VU, VL, T):
        Vi = sqrt(VU * VL)
        for _ in range(5):
            Vi = (VL - VU) * exp(-5.4 * T / Vi ** (1 / 3)) + VU
        r = 1 - exp(-5.4 * T / Vi ** (1 / 3))
        return Vi, r

    VA, ra = hat(VM2, VM3, T2)
    VA, rb = hat(VM1, VA, T1)

    r3 = ra * rb
    r2 = rb - r3
    r1 = 1 - rb

    CR = 1.2 * VA ** (1 / 3)
    CD = 0.5 * VA ** (1 / 3)

    if GR is not None:

        if GR > 10000 or GR < 1.8 * CR:
            reliable = False
        else:
            reliable = True
        EJ = 0
        for M, V, r in zip((M1, M2, M3), (VM1, VM2, VM3), (r1, r2, r3)):
            if M == 1 or M == 2:
                k = 0.9
            else:
                k = 1.17
            EJ += r * k * V**1.62 * GR**-3.86  # stacking of ejecta
    else:
        EJ = None
        reliable = False

    if Y < 0.1 or Y > 25000:
        isWithinLimit = False
    if H < -40 * Y ** (1 / 3):
        isWithinLimit = False
    if T1 > 1000:
        isWithinLimit = False
    if T2 > 1000:
        isWithinLimit = False

    return VA, CR, CD, isWithinLimit, EJ, reliable


def fallout(Y, H, DW, CW, W, SY, FF, T, TI, TEXP):
    if H < 0:
        raise ValueError("invalid burst height")
    if DW < 0 or CW < 0:
        raise ValueError("invalid wind ground range")
    if W < 0:
        raise ValueError("invalid wind speed")
    if SY < 0:
        raise ValueError("invalid crosswind shear")
    if FF < 0 or FF > 1:
        raise ValueError("impossible fission fraction")
    if T <= 0 or TI <= 0 or TEXP <= 0:
        raise ValueError("invalid time")

    Hhat = H / 0.3048

    if Hhat == 0:
        AF = 1
    else:
        z = 0.01 * Hhat / Y**0.4

        if z > 1:
            AF = 0
        else:
            AF = 0.5 * (1 - z) ** 2 * (2 + z) + 0.001 * z

    Ym = Y / 1000

    # cloud radius in nm
    if Y >= 1:
        sigma0 = Ym ** (1 / 3) * exp(0.56 - 3.25 / (4 + (log(Ym) + 5.4) ** 2))
    else:
        sigma0 = 0.1 * Y**0.2665

    # cloud height in kft

    if Y >= 1:
        h0 = 44 + 6.1 * log(Ym) - 0.205 * (log(Ym) + 2.42) * abs(log(Ym) + 2.42)
    else:
        h0 = 6 * Y**0.25
    # cloud duration (hours)
    Ti = 1.057 * h0 * (0.2 - h0 / 1440) * (1 - exp(-(h0**2) / 625) / 2)
    # cloud thickness (nm)
    sigmah = 0.18 * h0
    # effective particle distance (nm)
    L0 = W * Ti
    # change in fallout distribution
    sigmax = sigma0 * sqrt((L0**2 + 8 * sigma0**2) / (L0**2 + 2 * sigma0**2))
    # modified L0
    L = sqrt(L0**2 + 2 * sigmax**2)
    # constant for symmetry
    N = (L0**2 + sigmax**2) / (L0**2 + 0.5 * sigmax**2)
    # downwind distance
    d = DW / 1853
    c = CW / 1853
    # area reduction
    alpha1 = (1 + 0.001 * h0 * W / sigma0) ** (-1)
    alpha2 = (1 + 0.001 * h0 * W / sigma0 * (1 - phi(2 * d / W))) ** (-1)
    # crosswind spread parameter:
    P = SY * Ti * sigmah / L
    Q = abs(d + 2 * sigmax) / L
    R = min(4, 1 + 8 * Q)
    sigmay = sqrt(sigma0**2 * R + 2 * (P * sigmax) ** 2 + (P * Q * L0) ** 2)
    # crosswind transport function:
    F1 = exp(-((c / (alpha2 * sigmay)) ** 2) / 2) / (sigmay * sqrt(2 * pi))
    # downwind transport function
    F2 = phi(L0 * d / (L * alpha1 * sigmax))
    # deposition function
    F3 = exp(-abs(d / L) ** N) / (L * gamma(1 + 1 / N))
    F = F1 * F2 * F3
    # one hour dose rate
    DHP1 = 1510 * Y * FF * AF * F
    # debris arrival time
    T0 = sqrt(
        0.25 + ((L0 * Q * Ti) ** 2 + 2 * sigmax**2) / (L0**2 + 0.5 * sigmax**2)
    )
    z0 = log(T0)
    MBD = DHP1 * (2.737 - 0.7809 * z0 + 2 * z0**2 / 29 - z0**3 / 617)
    DT = DHP1 * T ** (-1.2)
    FD = 5 * DHP1 * (TI**-0.2 - (TI + TEXP) ** -0.2)

    withinLimit = True
    if DW > 1e6:
        withinLimit = False
    if CW > 4e4:
        withinLimit = False
    if W < 1 or W > 40:
        withinLimit = False
    if SY < 0 or SY > 10:
        withinLimit = False
    if T < 0.1 or T > 5000:
        withinLimit = False
    if TI < 0.1 or TI > 5000:
        withinLimit = False
    if TEXP < 0 or TEXP > 5000 - TI:
        withinLimit = False

    return DHP1, DT, T0, MBD, FD, withinLimit


if __name__ == "__main__":
    print("----------")
    print(*iniRad(1, 0.975, 10, 1000, 0.85, 10), sep="\n")

    print("----------")
    print(therm(1, 0, 500, 50000), sep="\n")

    print("----------")
    print(*crater(300, -10, 1, 1, 1, 1, 1, None, False), sep="\n")

    print("----------")
    print(*fallout(1, 20, 0, 0, 5, 0, 0.5, 1, 1, 1), sep="\n")
