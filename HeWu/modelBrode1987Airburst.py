"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

implements the Brode, 1978 model from the book:
PSR Report 1419-3
"AIRBLAST FROM NUCLEAR BURSTS—ANALYTIC APPROXIMATIONS"
Harold L. Brode, July 1987
Technical Report
CONTRACT No. DNA 001-85-C-0089  

This section mainly reflects equation from section IV, falling back to using the 
equations in III when in the free-air region.


"""
from math import log10, exp
from HeWu.uc import _uc_m2ft, _uc_psi2pa, _uc_ft2m

from HeWu.intg import intg


def _DeltaP_s(GR, H, W):
    """
    Peak overpressure in psi

    GR: ground range, in ft
    H: height, in ft
    W: yield, kT

    """
    m = W ** (1 / 3)

    z = H / GR

    a = 1.22 - 3.908 * z**2 / (1 + 810.2 * z**5)
    b = (
        2.321
        + 6.195 * z**18 / (1 + 1.113 * z**18)
        - 0.03831 * z**17 / (1 + 0.02415 * z**17)
        + 0.6692 / (1 + 4164 * z**8)
    )

    c = 4.153 - 1.149 * z**18 / (1 + 1.641 * z**18) - 1.1 / (1 + 2.771 * z**2.5)
    d = -4.166 + 25.76 * z**1.75 / (1 + 1.382 * z**18) + 8.257 * z / (1 + 3.219 * z)
    e = 1 - 0.004642 * z**18 / (1 + 0.003886 * z**18)
    f = (
        0.6096
        + 2.879 * z**9.25 / (1 + 2.359 * z**14.5)
        - 17.15 * z**2 / (1 + 71.66 * z**3)
    )
    g = 1.83 + 5.361 * z**2 / (1 + 0.3139 * z**6)

    # In-text version

    x = GR / m / 1000
    y = H / m / 1000

    r = (x**2 + y**2) ** 0.5

    h = (
        8.808 * z**1.5 / (1 + 154.5 * z**3.5)
        - (0.2905 + 64.67 * z**5) / (1 + 441.5 * z**5)
        - 1.389 * z / (1 + 49.03 * z**5)
        + 1.094
        * r**2
        / ((781.2 - 123.4 * r + 37.98 * r**1.5 + r**2) * (1 + 2 * y))
    )

    j = 0.000629 * y**4 / (3.493e-9 + y**4) - 2.67 * y**2 / (1 + 1e7 * y**4.3)
    k = 5.18 + 0.2803 * y**3.5 / (3.788e-6 + y * 4)
    DeltaP_s = 10.47 / r**a + b / r**c + d * e / (1 + f * r**g) + h + j / r**k

    return DeltaP_s


def _u(r):
    """
    r: scaled range in kilofeet per cube-root kiloton

    """

    return (
        (0.543 - 21.8 * r + 386 * r**2 + 2383 * r**3)
        * r**8
        / (
            2.99e-14
            - 1.91e-10 * r**2
            + 1.032e-6 * r**4
            - 4.43e-6 * r**6
            + (1.028 + 2.087 * r + 2.69 * r**2) * r**8
        )
    )


def _w(r):
    """
    r: scaled range in kilofeet per cube-root kiloton

    """
    return (
        (1.086 - 34.605 * r + 486.3 * r**2 + 2383 * r**3)
        * r**8
        / (
            3.0137e-13
            - 1.2128e-9 * r**2
            + 4.128e-6 * r**4
            - 1.116e-5 * r**6
            + (1.632 + 2.629 * r + 2.69 * r**2) * r**8
        )
    )


def _Xm(GR, H, W):
    """
    scaled range at which Mach reflection begins for a given burst height,
    onset of Mach reflection locus,in ft per kT**(1/3)   

                                        __ path of triple point      
                                   ____/
                            ______/    \\
                    _______/            \\Mach stem (merged front)
    _______________/---------------------|----
    reg.ref.region | mach reflection region
                   Xm

    GR: ground range in feet
    H: height of burst in feet
    W: yield in kiloton
    """
    m = W ** (1 / 3)
    X = GR / m  # ft/kT^(1/3)
    Y = H / m  # ft/kT^(1/3)

    Xm = 170 * Y / (1 + 60 * Y**0.25) + 2.89 * (Y / 100) ** 2.5
    return Xm


def _tau(GR, H, W, Xm=None):
    """
    scaled time of arrival for GR and H, in ms/kT^(1/3), based on Eq. (41).

    GR: ground range in feet
    H: height of burst in feet
    W: yield in kiloton

    Xm: optional, when supplied skips the call to _Xm
    """
    m = W ** (1 / 3)
    x = GR / m  # ft/kT^(1/3)
    y = H / m  # ft/kT^(1/3)

    if Xm is None:
        Xm = _Xm(GR, H, W)

    rm = (Xm**2 + y**2) ** 0.5 / 1000
    r = (x**2 + y**2) ** 0.5 / 1000
    if x <= Xm:
        tau = _u(r)
    else:
        tau = _u(rm) + _w(r) - _w(rm)

    return tau


def _D(GR, H, W, tau=None, Xm=None):
    """
    overpressure duration of positive phase in milliseconds

    GR: ground range in feet
    H: burst height in feet
    W: yield in kiloton

    tau: optional, supplied to skip a call to _tau
    """

    m = W ** (1 / 3)
    X = GR / m  # ft/kT^(1/3)
    Y = H / m  # ft/kT^(1/3)

    if tau is None:
        if Xm is None:
            Xm = _Xm(GR, H, W)
        tau = _tau(GR, H, W, Xm)

    s2 = (
        1
        - 15.18 * (Y / 100) ** 3.5 / (1 + 15.18 * (Y / 100) ** 3.5)
        - 0.02441
        * (Y / 1e6) ** 2
        / (1 + 9e3 * (Y / 100) ** 7)
        * 1e10
        / (0.441 + (X / 100) ** 10)
    )

    """ a very careful reading of the work suggest that the duration time must
    also scale with yield, but the original formulation does not change. We fix
    it by multiplying m onto D."""
    D = (
        (1640700 + 24629 * tau + 416.15 * tau**2)
        / (10880 + 619.76 * tau + tau**2)
        * (
            0.4
            + 0.001204 * tau**1.5 / (1 + 0.001559 * tau**1.5)
            + (
                0.6126
                + 0.5486 * tau**0.25 / (1 + 0.00357 * tau**1.5)
                - 3.47 * tau**0.637 / (1 + 5.696 * tau**0.645)
            )
            * s2
        )
    ) * m  # duration of positive phase in milliseconds,

    return D


def _D_u(GR, H, W, D=None, Xm=None, DeltaP_s=None, tau=None):
    """
    positive phase duration for dynamic pressure in milliseconds,
    duration of outward blast wind

    The below use the more complex form seen in Eqn.63). However, it is realized
    that the below form can go below positive. Further, this behaviour is mostly
    seen in case where H>GR, i.e. close to the situation of a free-air burst.

    Therefore, we transition to a simpler fit in free air, Eqn.51)

    GR: ground range, feet
    H: height of burst, feet
    W: yield, kiloton
    """

    m = W ** (1 / 3)
    x = GR / m / 1000  # kft/kT^(1/3)
    y = H / m / 1000

    if Xm is None:
        Xm = _Xm(GR, H, W)

    if tau is None:
        tau = _tau(GR, H, W, Xm)

    if x * 1000 < Xm:  # this is fitted to DeltaP_s so should work by first principle
        if DeltaP_s is None:
            _pi = _DeltaP_s(GR, H, W) / 1000  # in ksi
        else:
            _pi = DeltaP_s / 1000

        """ see equation 52) """

        D_u_pos = m * (
            317 / (1 + 85 * _pi + 7500 * _pi**2)
            + 6110 * _pi / (1 + 420 * _pi**2)
            + 2113 * _pi / (1 + 11 * _pi)
        )  # this correctly scales with yield.
        return D_u_pos

    else:  # in surface-burst
        if D is None:
            D = _D(GR, H, W, tau, Xm)

        C = (
            89.6 * y**5.2 / (1 + 20.5 * y**5.4)
            + 4.51 / (1 + 130.7 * y**8.6)
            + 2.466 * y**0.5 / (1 + 99 * y**2.5)
            - 12.8 * (x**2 + y**2) ** 1.25 / (1 + 3.63 * (x**2 + y**2) ** 1.25)
        )
        return C * D


def _DeltaP(GR, H, W, t, integrate=True):
    """
    Overpressure over time in psi.
    Optionally, integrates it from arrival time to specified time.

        "This time-dependent formulation is fitted to the BM-3 cratering calculation
        impulses [Pyatt, 1983] in the close-in range (scaled ground range between
        25 and 100ft, or peak pressure between 1.5-150 ksi)......
        ......This may not be correct, although further calculation are being made
        to verify that region of the BM-3 calculation. It has been suggested that
        the lower impulse may be due ot the quenching effect of cratering ejected
        into the fireball"

    input:
        GR: ground range in feet
        H: burst height in feet
        W: yield in kiloton
        t: time after detonation in milliseconds.

        integrate: boolean value, controls whether an integration is done over time
        from the time of arrival to the supplied time.

    return:
        integrate = False:
            overpressure in psi
        integrate = True:
            overpressure in psi, overpressure impulse in psi-ms

    """
    m = W ** (1 / 3)
    X = GR / m  # ft/kT^(1/3)
    Y = H / m  # ft/kT^(1/3)

    z = H / GR

    Xm = _Xm(GR, H, W)  # onset of Mach reflection locus, scaled, in ft per kT**(1/3)

    Xe = (
        3.039 * Y / (1 + 0.0067 * Y)
    )  # locus of points where second peak equals first peak, scaled in ft per kT**(1/3)

    """

    y  x = Xm _____
    |      __/  /
    |    _/    / 
    |r.r/    _/ x = Xe
    |  /2>1_/ 
    | / _/
    |/_/   1>2
    //--------------x

    r.r.: regular reflection region
    2>1 : mach reflection region, second peak > first peak
    1>2 : mach reflection region, second peak < first peak

    """

    K = abs((X - Xm) / (Xe - Xm))

    d2 = 2.99 + 31240 * (Y / 100) ** 9.86 / (1 + 15530 * (Y / 100) ** 9.87)
    d = (
        0.23
        + 0.583 * Y**2 / (26667 + Y**2)
        + 0.27 * K
        + (0.5 - 0.583 * Y**2 / (26667 + Y**2)) * K**d2
    )

    a = (d - 1) * (1 - K**20 / (1 + K**20))

    r = (X**2 + Y**2) ** 0.5 / 1000

    tau = _tau(GR, H, W, Xm)

    sd = _D(GR, H, W, tau, Xm) / m  # scaled time of arrival

    s = (
        1
        - 1100 * (Y / 100) ** 7 / (1 + 1100 * (Y / 100) ** 7)
        - 2.441e-14
        * Y**2
        / (1 + 9000 * (Y / 100) ** 7)
        * 1e10
        / (0.441 + (X / 100) ** 10)
    )

    f2 = (
        (
            0.445
            - 5.44 * r**1.02 / (1 + 1e5 * r**5.84)
            + 7.571 * z**7.15 / (1 + 5.135 * z**12.9)
            - 8.07 * z**7.31 / (1 + 5.583 * z**12.23)
        )
        * 0.4530  # 0.4530 in paper, 0.435 in FORTAN
        * (Y / 10) ** 1.26
        / (1 + 0.03096 * (Y / 10) ** 3.12)
        * (1 - 0.000019 * tau**8 / (1 + 0.000019 * tau**8))
    )

    f = (
        (
            0.01477 * tau**0.75 / (1 + 0.005836 * tau)
            + 7.402e-5 * tau**2.5 / (1 + 1.429e-8 * tau**4.75)
            - 0.216
        )
        * s
        + 0.7076
        - 3.077e-5 * tau**3 / (1 + 4.367e-5 * tau**3)
        + f2
        - (0.452 - 9.94e-7 * X**4.13 / (1 + 2.1868e-6 * X**4.13))
        * (1 - 1.5397e-4 * Y**4.3 / (1 + 1.5397e-4 * Y**4.3))
    )

    g = (
        10 + (77.58 - 64.99 * tau**0.125 / (1 + 0.04348 * tau**0.5)) * s
    )  # early time decay power, tau raised to 0.5 in work, 5 in FORTRAN
    # which ever one is better is not established at this moment

    h = (
        3.003
        + 0.05601 * tau / (1 + 1.473e-9 * tau**5)
        + (
            0.01769 * tau / (1 + 3.207e-10 * tau**4.25)
            - 0.03209 * tau**1.25 / (1 + 9.914e-8 * tau**4)
            - 1.6
        )
        * s
        - 0.1966 * tau**1.22 / (1 + 0.767 * tau**1.22)
    )

    # supporting values for calculatng c. Time-independent

    c2 = 23000 * (Y / 100) ** 9 / (1 + 23000 * (Y / 100) ** 9)
    c3 = 1 + (
        1.094
        * K**0.738
        / (1 + 3.687 * K**2.63)
        * (1 - 83.01 * (Y / 100) ** 6.5 / (1 + 172.3 * (Y / 100) ** 6.04))
        - 0.15
    ) / (1 + 0.5089 * K**13)

    """
    we split the calculation of v into v = v0 * vj to simplify per iteration
    calculation
    """
    v0 = (
        0.003744 * (Y / 10) ** 5.185 / (1 + 0.004684 * (Y / 10) ** 4.189)
        + 0.004755 * (Y / 10) ** 8.049 / (1 + 0.003444 * (Y / 10) ** 7.497)
        - 0.04852 * (Y / 10) ** 3.423 / (1 + 0.03038 * (Y / 10) ** 2.538)
    ) / (1 + 9.23 * K**2)

    def atSigma(sigma):
        """partial function, implicitly inheriting the non-time-dependent factors
        sigma: scaled time in question

        returns the ratio: DeltaP @ scaled time sigma / DeltaP_s
        """

        b = (f * (tau / sigma) ** g + (1 - f) * (tau / sigma) ** h) * (
            1 - (sigma - tau) / sd
        )

        if X >= Xm and Y <= 380:
            j = min(
                11860 * (sigma - tau) / (Y * (X - Xm) ** 1.25), 200
            )  # ratio of time after TOA to time to second peak after TOA,

            v = v0 * j**3 / (6.13 + j**3) + 1

            c = (
                (
                    (1.04 - 0.02409 * (X / 100) ** 4 / (1 + 0.02317 * (X / 100) ** 4))
                    * j**7
                    / ((1 + a) * (1 + 0.923 * j**8.5))
                )
                * (c2 + (1 - c2) * (1 - 0.09 * K**2.5 / (1 + 0.09 * K**2.5)))
                * c3
                * (1 - ((sigma - tau) / sd) ** 8)
            )  # this is harder to split into time and non-time dependent part and left as is

            return (1 + a) * (b * v + c)
        else:
            return b

    DeltaP_s = _DeltaP_s(GR, H, W)

    sigma = t / m  # sigma is the scaled time in ms/kT^(1/3)

    start = tau
    end = tau + sd

    if start - sigma > 1e-9:
        raise ValueError(
            "blast wave hasn't arrived at the specified time ( {} ms < {} ms )".format(
                t, start * m
            )
        )
    if sigma - end > 1e-9:

        raise ValueError(
            "positive phase for over pressure is over ( {} ms > {} ms )".format(
                t, end * m
            )
        )

    if integrate:
        """the integral DeltaP_s * intg(atSigma, tau, sigma)[0]
        is the scaled overpressure total impulse, psi-ms/kT^(1/3)
        """
        return (
            DeltaP_s * atSigma(sigma),
            DeltaP_s * intg(atSigma, tau, sigma)[0] * m,
        )
    else:
        return DeltaP_s * atSigma(sigma)


def _Q(GR, H, W, t, integrate=True):
    """
    Overpressure or dynamic pressure horizontal component over time in psi.
    Optionally, integrates it from arrival time to specified time.

        "The peak dynamic pressure and the dynamic impulse derived from this
        quick fix fit are not as accurate as those given by Eqs. (64) and (65).
        To illustrate the limited usefulness of this approximation [Eqs. (66) and (67)],
        the peak dynamic pressures are compared with those from the KA
        calculations in Figs. Ill and 112, all at the scaled burst height of 200 ft/KT^(1/3).
        While the peaka from this time-history fit are appreciably low at the innermost
        range plotted in Fig. Ill (24.5 psi versus 30 psi), the agreement gets better
        at larger ranges (Fig. 112). The impulses, shown in Figs. 113 and 114, are
        similarly poor at the closest range shown, but they are in good agreement far¬
        ther out at intermediate ranges. This fit falls about 30 percent low on peak
        and impulse at the 100 psi overpressure range. This quick fix time-history fit
         was provided as an analytic expression useful in dynamic analyses in a limited
        (low) overpressure range, and should not be used at high overpressures.
        The expressions for peak dynamic pressure and dynamic impulse [Eqs. (64)
        and (65)] are more accurate for peaks and total impulse, but they do not
        provide the time- dependent or transient behavior necessary for calculations
        of response of structures or vehicles. It is anticipated that this quick fix
        time-history fit will be improved in the near future."

    input:
        GR: ground range in feet
        H: burst height in feet
        W: yield in kiloton
        t: time or times in milliseconds.

        integrate: boolean value, controls whether an integration is done over time
        from the time of arrival to the supplied time.

    return:
        integrate = False:
            dynm.press.hz.component in psi
        integrate = True:
            dynm.press.hz.component in psi, dynm.press.hz.component impulse in psi-ms

    """
    m = W ** (1 / 3)
    X = GR / m  # ft/kT^(1/3)
    Y = H / m  # ft/kT^(1/3)
    z = H / GR

    Xm = _Xm(GR, H, W)  # onset of Mach reflection locus, scaled, in ft per kT**(1/3)

    Xe = (
        3.039 * Y / (1 + 0.0067 * Y)
    )  # locus of points where second peak equals first peak, scaled in ft per kT**(1/3)

    K = abs((X - Xm) / (Xe - Xm))

    d2 = 2.99 + 31240 * (Y / 100) ** 9.86 / (1 + 15530 * (Y / 100) ** 9.87)
    d = (
        0.23
        + 0.583 * Y**2 / (26667 + Y**2)
        + 0.27 * K
        + (0.5 - 0.583 * Y**2 / (26667 + Y**2)) * K**d2
    )

    a = (d - 1) * (1 - K**20 / (1 + K**20))

    r = (X**2 + Y**2) ** 0.5 / 1000

    tau = _tau(GR, H, W, Xm)

    DeltaP_s = _DeltaP_s(GR, H, W)

    sD_u = _D_u(GR, H, W, None, Xm, DeltaP_s, tau) / m  # scaled time

    s = (
        1
        - 1100 * (Y / 100) ** 7 / (1 + 1100 * (Y / 100) ** 7)
        - 2.441e-14
        * Y**2
        / (1 + 9000 * (Y / 100) ** 7)
        * 1e10
        / (0.441 + (X / 100) ** 10)
    )

    f2 = (
        (
            0.445
            - 5.44 * r**1.02 / (1 + 1e5 * r**5.84)
            + 7.571 * z**7.15 / (1 + 5.135 * z**12.9)
            - 8.07 * z**7.31 / (1 + 5.583 * z**12.23)
        )
        * 0.4530  # 0.4530 in paper, 0.435 in FORTAN
        * (Y / 10) ** 1.26
        / (1 + 0.03096 * (Y / 10) ** 3.12)
        * (1 - 0.000019 * tau**8 / (1 + 0.000019 * tau**8))
    )

    f = (
        (
            0.01477 * tau**0.75 / (1 + 0.005836 * tau)
            + 7.402e-5 * tau**2.5 / (1 + 1.429e-8 * tau**4.75)
            - 0.216
        )
        * s
        + 0.7076
        - 3.077e-5 * tau**3 / (1 + 4.367e-5 * tau**3)
        + f2
        - (0.452 - 9.94e-7 * X**4.13 / (1 + 2.1868e-6 * X**4.13))
        * (1 - 0.00015397 * Y**4.3 / (1 + 0.00015397 * Y**4.3))
    )

    g = (
        10 + (77.58 - 64.99 * tau**0.125 / (1 + 0.04348 * tau**0.5)) * s
    )  # early time decay power

    h = (
        3.003
        + 0.05601 * tau / (1 + 1.473e-9 * tau**5)
        + (
            0.01769 * tau / (1 + 3.207e-10 * tau**4.25)
            - 0.03209 * tau**1.25 / (1 + 9.914e-8 * tau**4)
            - 1.6
        )
        * s
        - 0.1966 * tau**1.22 / (1 + 0.767 * tau**1.22)
    )

    # supporting values for calculatng c. Time-independent

    c2 = 23000 * (Y / 100) ** 9 / (1 + 23000 * (Y / 100) ** 9)
    c3 = 1 + (
        1.094
        * K**0.738
        / (1 + 3.687 * K**2.63)
        * (1 - 83.01 * (Y / 100) ** 6.5 / (1 + 172.3 * (Y / 100) ** 6.04))
        - 0.15
    ) / (1 + 0.5089 * K**13)

    """
    we split the calculation of v into v = v0 * vj to simplify per iteration
    calculation
    """
    v0 = (
        0.003744 * (Y / 10) ** 5.185 / (1 + 0.004684 * (Y / 10) ** 4.189)
        + 0.004755 * (Y / 10) ** 8.049 / (1 + 0.003444 * (Y / 10) ** 7.497)
        - 0.04852 * (Y / 10) ** 3.423 / (1 + 0.03038 * (Y / 10) ** 2.538)
    ) / (1 + 9.23 * K**2)

    """ time independent scaling factors for Q calculation"""

    y = H / m / 1000
    a1 = 2 - 2 / (1 + 3817 * y**9)  # a in the original
    b1 = 2 + 1.011 / (1 + 33660 * y**15)  # b in the original

    def atSigma(sigma):
        """partial function, implicitly inheriting the non-time-dependent factors
        sigma: scaled time in question

        returns the ratio: Q, dynamic pressure hz.component in psi
        """

        b = (f * (tau / sigma) ** g + (1 - f) * (tau / sigma) ** h) * (
            1 - (sigma - tau) / sD_u
        )

        if X >= Xm and Y <= 380:
            j = min(
                11860 * (sigma - tau) / (Y * (X - Xm) ** 1.25), 200
            )  # ratio of time after TOA to time to second peak after TOA,

            v = v0 * j**3 / (6.13 + j**3) + 1

            c = (
                (
                    (1.04 - 0.02409 * (X / 100) ** 4 / (1 + 0.02317 * (X / 100) ** 4))
                    * j**7
                    / ((1 + a) * (1 + 0.923 * j**8.5))
                )
                * (c2 + (1 - c2) * (1 - 0.09 * K**2.5 / (1 + 0.09 * K**2.5)))
                * c3
                * (1 - ((sigma - tau) / sD_u) ** 8)
            )  # this is harder to split into time and non-time dependent part and left as is

            ratio = (1 + a) * (b * v + c)

        else:
            ratio = b

        DeltaP = ratio * DeltaP_s

        if X >= 1.3 * Xm:
            return ratio**a1 * 2.5 * DeltaP**2 / (102.9 + DeltaP)
        else:
            return (
                ratio**a1
                * 2.5
                * DeltaP**2
                / (102.9 + DeltaP)
                * (X / (1.3 * Xm)) ** b1
            )

    sigma = t / m  # sigma is the scaled time in ms/kT^(1/3)
    start = tau
    end = tau + sD_u

    if start - sigma > 1e-9:
        raise ValueError(
            "blast wave hasn't arrived at the specified time ( {} ms < {} ms )".format(
                t, start * m
            )
        )
    if sigma - end > 1e-9:
        raise ValueError(
            "positive phase for dynamic pressure is over ( {} ms > {} ms )".format(
                t, end * m
            )
        )

    if integrate:
        """the integral intg(atSigma, tau, sigma)[0]
        is the scaled dynamic hz. impulse, psi-ms/kT^(1/3)
        """
        return atSigma(sigma), intg(atSigma, tau, sigma)[0] * m
    else:
        return atSigma(sigma)


def _Q_1(x, y, xq):
    """
    helper function for _Q_s
    """

    r = (x**2 + y**2) ** 0.5
    M = xq / x

    A = -236.1 + 17.72 * M**0.593 / (1 + 10.4 * M**3.124)
    B = 12.27 - 21.69 * M**2.24 / (1 + 6.976 * M**0.484)
    C = 20.26 + 14.7 * M**2 / (1 + 0.08747 * M**3.05)
    D = -1.137 - 0.5606 * M**0.895 / (1 + 3.046 * M**7.48)
    E = 1.731 + 10.84 * M**1.12 / (1 + 12.26 * M**0.0014)
    F = 2.84 + 0.855 * M**0.9 / (1 + 1.05 * M**2.84)

    return A * r**D / (1 + B * r**E) + C / r**F


def _Q_s(GR, H, W):
    """
    Peak (horizontal) dynamic pressure in psi

    GR: ground range, feet
    H: burst height, feet
    W: yield, kiloton
    """
    m = W ** (1 / 3)
    x = GR / 1000 / m  # scaled ground range in kft/kt^1/3
    y = H / 1000 / m  # scaled ground range in kft/kft^1/3

    r = (x**2 + y**2) ** 0.5

    xq = (
        63.5 * y**7.26 / (1 + 67.11 * y**4.746) + 0.6953 * y**0.808
    )  # approximate interface between regular and Mach reflection in kft/kt^(1/3)

    Q1 = _Q_1(x, y, xq)
    Qm = _Q_1(xq, y, xq)  # Q_1 evaluated at x=xq, and M = 1

    G = 50 - 1843 * y**2.153 / (1 + 3.95 * y**5.08)
    H = 0.294 + 71.56 * y**8.7 / (1 + 115.3 * y**6.186)
    I = abs(-3.324 + 987.5 * y**4.77 / (1 + 211.8 * y**5.166))
    J = 1.955 + 169.7 * y**9.317 / (1 + 97.36 * y**6.513)
    K = 8.123e-6 + 0.001613 * y**6.428 / (1 + 60.26 * y**7.358)
    L = log10(xq / x)

    if x >= xq:
        return Q1
    else:
        return Qm * exp(
            G * L**I / (1 + 649 * L**I)
            - 4.01 * L**J / (1 + H * L**J)
            + 7.67e-6 * (1 / (K + L**3.22) - 1 / K)
        )


def _I_u_pos(GR, H, W, DeltaP_s=None):
    """
    simple fit for integral of dynamic pressure, psi-ms, with time over the
    positive (outward flow) phase. This is only valid in the Mach reflection
    region but is computationally simpler than integration routine above.

    GR: ground range, ft
    H: height of burst, ft
    W: yield, kiloton

    """
    m = W ** (1 / 3)
    x = GR / 1000 / m
    y = H / 1000 / m

    psi = y + 0.09

    E = 183 * (y**2 + 0.00182) / (y**2 + 0.00222)
    F = 0.00058 * exp(9.5 * y) + 0.0117 * exp(-22 * y)
    G = 2.3 + 29 * y / (1 + 1760 * y**5) + 25 * y**4 / (1 + 3.76 * y**6)

    if x > 170 * psi / (1 + 337 * psi**0.25) + 0.914 * psi**2.5:  # x> Xi
        return (E * x / (F + x**3.61) + G / (1 + 0.22 * x**2)) * m
    else:
        return None  # we were unable to source a good enough estimation for this


def _I_p_pos(GR, H, W, DeltaP_s=None, Xm=None):
    """
    positive phase impulse in psi-ms, using a simple fit against overpressure peak
    this should work for airburst as well as free-air.

    GR: ground range, feet
    H: height of burst, feet
    W: yield, kiloton

    """
    m = W ** (1 / 3)
    X = GR / m  # ft/kT^(1/3)
    Y = H / m  # ft/kT^(1/3)

    if DeltaP_s is None:
        DeltaP_s = _DeltaP_s(GR, H, W)

    if Xm is None:
        Xm = _Xm(GR, H, W)

    if X <= Xm:  # airburst
        return 145 * DeltaP_s**0.5 / (1 + 0.00385 * DeltaP_s**0.5) * m
    else:  # surface burst
        return 183 * DeltaP_s**0.5 / (1 + 0.00385 * DeltaP_s**0.5) * m


def airburst(GR_m, H_m, W, t=None, prettyPrint=True):
    """
    Calculate various air-burst parameters, using the Brode 1987 model and adapting
    to SI unit system.

    The total integral for overpressure and dynamic pressure horizontal component
    are provided for illustrative purpose only. In terms of total impulse they are
    less accurate than using the estimates. However, the estimations behave poorly
    outside of their indicated range, and especially poorly at high pressures.

    input:
        GR: ground range, meter
        H : height of burst, meter
        W : yield, kiloton
        t : partial time after arrival, second, default to None
        prettyPrint: boolean, whether to produce a pretty print output using print()

    output:
        TAAIR  : time of arrival, second
        PAAIR  : maximum overpressure, Pa
        DPP    : overpressure positive phase duration, second
        IPTOTAL: overpressure total positive impulse, Pa-s
        IPEST  : estimation of the above, Pa-s
        PPART  : overpressure calculated at partial time, Pa
        IPPART : overpressure impulse integrated to partial time, Pa-s
        QAAIR  : maximum dynamic pressure horizontal component, Pa
        DPQ    : dynamic pressure positive (outward flow) phase duration, s
        IQTOTAL: dynamic pressure horizontal component total impulse, Pa-s
        IQEST  : estimation of the above, Pa-s
        QPART  : dynamic pressure horizontal component at partial time, Pa
        IQPART : dynamic pressure horizontal impulse integrated to partial time, Pa-s
        XM     : "onset of Mach reflection locus", range at which Mach reflection starts, m

    """

    m = W ** (1 / 3)

    GR = max(_uc_m2ft(GR_m), 1e-9 * m)  # clamp the value to > 0.001 meter
    H = max(_uc_m2ft(H_m), 1e-9 * m)  # clamp the value to > 0.001 meter

    Xm = _Xm(GR, H, W)
    tau = _tau(GR, H, W, Xm)

    DeltaP_s = _DeltaP_s(GR, H, W)
    Q_s = _Q_s(GR, H, W)

    D = _D(GR, H, W, tau, Xm)
    D_u = _D_u(GR, H, W, D, Xm, DeltaP_s, tau)

    XM = _uc_ft2m(Xm * m)
    TAAIR = tau * m / 1000  # ms to s

    DPP = D / 1000  # ms to s
    DPQ = D_u / 1000  # ms to s

    p, I_p_pos = _DeltaP(GR, H, W, (tau * m + D), integrate=True)
    q, I_u_pos = _Q(GR, H, W, (tau * m + D_u), integrate=True)

    I_p_est = _I_p_pos(GR, H, W, DeltaP_s, Xm)
    I_u_est = _I_u_pos(GR, H, W, DeltaP_s)

    IPEST = _uc_psi2pa(I_p_est / 1000)
    if I_u_est is not None:
        IQEST = _uc_psi2pa(I_u_est / 1000)
    else:
        IQEST = None

    IPTOTAL = _uc_psi2pa(I_p_pos / 1000)
    IQTOTAL = _uc_psi2pa(I_u_pos / 1000)

    PAAIR = _uc_psi2pa(DeltaP_s)
    QAAIR = _uc_psi2pa(Q_s)

    if t is None:
        PPART, IPPART, QPART, IQPART = None, None, None, None
    else:
        try:
            pt, ipt = _DeltaP(GR, H, W, tau * m + t * 1000, integrate=True)
            PPART = _uc_psi2pa(pt)
            IPPART = _uc_psi2pa(ipt / 1000)
        except ValueError:
            PPART, IPPART = None, None
        try:
            qt, iqt = _Q(GR, H, W, tau * m + t * 1000, integrate=True)
            QPART = _uc_psi2pa(qt)
            IQPART = _uc_psi2pa(iqt / 1000)
        except ValueError:
            QPART, IQPART = None, None

    if prettyPrint:
        print("{:^49}".format("INPUT"))
        print("Ground Range:\n{:.>20,.6g} m".format(GR_m))
        print("Burst Height:\n{:.>20,.6g} m".format(H_m))
        print("Yield:\n{:.>20,.6g} kT".format(W))
        print("")
        print(
            "Time To:\n{:.>20} s".format(
                "######" if t is None else "{:,.6g}".format(t + TAAIR)
            )
        )
        print("")
        print("{:^49}".format("OUTPUT"))
        print("Time of Arrival:\n{:.>20,.6g} s".format(TAAIR))
        print("Range Onset of Mach Stem:\n{:.>20,.6g} m".format(XM))
        print("")
        print("{:>24}{:>24} ".format("Overpressure", "Dyn.Pres.Hz.Comp."))
        print("{:<24}{:<24} ".format("Peak:", ""))
        print("{:.>21,.6g} Pa{:.>21,.6g} Pa ".format(PAAIR, QAAIR))
        print("{:<24}{:<24} ".format("Duration:", ""))
        print("{:.>22,.6g} s{:.>22,.6g} s ".format(DPP, DPQ))
        print("{:<24}{:<24} ".format("Impulse:", ""))
        print(".integral.{:.>9,.6g} Pa-s{:.>19,.6g} Pa-s".format(IPTOTAL, IQTOTAL))
        print(
            ".estimate.{:.>9,.6g} Pa-s{:.>19} Pa-s".format(
                IPEST, "######" if IQEST is None else "{:,.6g}".format(IQEST)
            )
        )

        print("")

        print("{:<24}{:<24} ".format("Pressure at Time:", ""))
        print(
            "{:.>21} Pa{:.>21} Pa ".format(
                "######" if PPART is None else "{:,.6g}".format(PPART),
                "######" if QPART is None else "{:,.6g}".format(QPART),
            )
        )
        print("{:<24}{:<24} ".format("Partial Impulse to Time:", ""))
        print(
            ".integral.{:.>9} Pa-s{:.>19} Pa-s ".format(
                "######" if IPPART is None else "{:,.6g}".format(IPPART),
                "######" if IQPART is None else "{:,.6g}".format(IQPART),
            )
        )
        print("")

    return (
        TAAIR,
        PAAIR,
        DPP,
        IPTOTAL,
        IPEST,
        PPART,
        IPPART,
        QAAIR,
        DPQ,
        IQTOTAL,
        IQEST,
        QPART,
        IQPART,
        XM,
    )


if __name__ == "__main__":
    """
    by default, runs a test
    """
    from HeWu.test import runABtest

    def _airburst_to_op(gr, h, Y, t):
        _, _, _, _, _, p, _, _, _, _, _, _, _, _ = airburst(gr, h, Y, t, False)
        return p

    runABtest(_airburst_to_op)
