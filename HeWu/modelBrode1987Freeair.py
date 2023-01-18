"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

implements the Brode, 1978 model from the book:
PSR Report 1419-3
"AIRBLAST FROM NUCLEAR BURSTS—ANALYTIC APPROXIMATIONS"
Harold L. Brode, July 1987
Technical Report
CONTRACT No. DNA 001-85-C-0089  

This part mainly reflects Section III of the original work


sr = m x r
T = m x tau

X to Y means "fitting X to Y" or "fitting X against Y"
"""

from math import exp
from HeWu.uc import _uc_psi2pa, _uc_m2ft


def _DeltaP_s_to_r(r):
    """
    peak overpressure versus range in psi, as in equation 33)

    the ideal ground burst equation in 34) can be reproduced by supplying an
    argument of 2W instead of W in the scaling for the s.g.r.

    overpressure from near-surface burst observed early atmospheric nuclear test
    can be reproduced by scaling using 1.8W instead of 2W

    Arguments:
    - r: scaled range in kft/kT^(1/3)
    Returns:
    - DeltaP_s: peak overpressure in psi

    "
    An expression that is appropriate for overpressures in the range
    0.07 < ΔP_s < 400,000 psi.

    In general, the agreement with a detailed one-dimensional calculation
    [Brode, 1959b] is within 5 percent. Extension to pressures below 2 psi is
    in agreement with the "Combined Airborne Polynomial" to within 5 percent,
    as well. The term in brackets is unnecessary for overpressures below 200 psi
    (viz., r > 0.22 kft/KT^(1/3)). This fit is only 10 percent high at 600,000 psi.
    At the lower extreme, it is 20 percent high at 0.01 psi in comparison with
    the suggested curve, which is an extension of the Airborne curve and other
    (theoretical) extrapolations. The fit is in agreement with the 1-KT standard
    [Needham and Crepeau, 1981] to within 10 percent (beyond 20-m scaled range).

    The accuracy of the fit should not mislead the user into believing that such
    overpressures are repeatable or reproducible with compparable accuracy.
    A review of the original atmospheric nuclear test measurements suggests a
    much larger data scatter even for the controlled test-site conditions and
    surfaces. The scatter in test data far exceeds the few percent difference
    between "best-guess" curves and fits. Figure 11 is a graph of approximate
    data scatter, as derived from Brode [1981], plotted against peak overpressure.
    At the higher overpressure, the plot suggests a minimum uncertainty of 70 percent
    (from minimum to maximum datum), which translates into a pressure variation
    of a factor of 5.
    "
    """
    return (
        2 / r
        + 3 / r**1.5
        + 1.6 / r**3 * (1 + 105 * (10 * r) ** 4 / (1 + 12920 * (10 * r) ** 8))
    )  # Eqn. (33)


def _r(DeltaP_s):
    """
    shock radius versus peak overpressure for free air burst, in kft/kT^(1/3)

    Arguments:
    - DeltaP_s: peak overpressure in psi

    Returns:
    - r: scaled ground range in kft per cube root kiloton

    "
    The following fit can be used to search for the range to a given
    shock overpressure. Employing this fit saves one from the tedium of
    inverting the approximation of Eq. (33) for overpressure as a function of
    range. The fit is accurate to within 2 percent of Eq. (33):

    The comparable range for an ideal-surface burst (2W) is 2^1/3
    (~1.26 times larger). For bursts over the (empirical) test-site surface (1.8W),
    the range is ~1.22 times larger than that for a free-air burst.
    "
    """
    _pi = DeltaP_s / 1000

    return 2.463 / (DeltaP_s**0.9846) + 1.176 / (DeltaP_s ** (1 / 3)) * (
        1 + 0.0004726 * _pi**1.5 / (1 + 2.952e-6 * _pi**3)
    )  # Eqn. (36)


def _DeltaP_s_to_tau(tau):
    """
    peak overpressure versus scaled time of arrival, equation 37)

    apply the same procedure when scaling for ideal or test-surface bursts, as in
    equation 33) or _DeltaP_s_to_r, see above.

    Arguments:
    - tau: scaled time of arrival, T/m, in milliseconds per cube root kiloton
    Returns:
    - DeltaP_s: peak overpressure in psi

    "
    This fit to the time of arrival and peak overpressure is based on
    early detailed calculations [Brode, 1959b, 1966] scaled to 1 KT for
    pressures between 2 <= ΔP_s <= 500,000 psi. The fit is good to better than
    4 percent below 50,000 psi, and to better than 3 percent below 10,000 psi.
    "
    """

    return (7000 + 5.18 * tau**0.75) / (
        0.007643 + tau**1.1334
    ) + 7.11e8 * tau**2.973 / (
        1 + 430800 * tau**4.656 + 3052e3 * tau**6.535
    )  # Eqn. (37)


def _tau_to_DeltaP_s(DeltaP_s):
    """
    scaled time of arrival versus peak overpressure

    Arguemnts:
    - DeltaP_s, peak overpressure in psi
    Returns:
    - tau: scaled time of arrival in millisecond per cube root kiloton

    "
    Conversely, the time of arrival T can be expressed as a function of the peak
    overpressure (for a free-air burst) as Eqn.39).
    For a surface burst, one uses 2W. The above form is accurate to ±2 percent
    for 2 <= ΔP_s <= 100,000 psi

    As is evident in Fig. 12, the time of arrival is nearly inversely proportional
    to the overpressure. Actually, the product of time of arrival and peak
    overpressure to the fractional power 0.875 is a slowly varying function,
    and can generate a more readable curve. That relation is plotted in Fig. 13.
    The approximate proportionality means that the product of overpressure times
    arrival time to the 1.14 power is roughly constant (APS x T^1.14 ~ constant).
    "

    """

    return 0.03394 + 893 / DeltaP_s**0.80424 + 2015 / DeltaP_s  # Eqn.(39)


def _tau_to_r(r):
    """

    Arguemnts:
    - r, scaled ground range in kilofeet per cube-root kiloton
    Returns:
    - tau: scaled time of arrival in millisecond per cube root kiloton

    "
    This shock arrival-time form is limited to the range of the "Empirical 59"
    data [Moulton, I960], namely, 620 μs to 26 s at 1 KT, which spans a range of
    free-air peak overpressures from 17,200 to 0.07 psi.
    (Use 2W in m for a surface burst; i.e., replace m with m x 2^(1/3).)
    "
    """

    return (0.54291 - 21.185 * r + 361.8 * r**2 + 2383 * r**3) / (
        1 + 2.048 * r + 2.6872 * r**2
    )  # Eqn. (40), scaled


def _ci_tau_to_r(r):
    """
    close-in time of arrival versus shock radius for higher pressure or earlier
    time of arrival T,  in millisecond per square root kiloton. We use the al-
    ternative formualtion of equation 42) dividing out the scaling factor m

    Arguments:
    - r, scaled ground range in kft per cube-root kiloton

    Returns:
    -tau: scaled time of arrival in milliseconds per cube root kiloton


    "
    For higher pressure (or earlier time of arrival T), the following fit is
    appropriate [Brode, 1970]: Eq.41)
    Equation (41) is alternatively written in terms of the scaled range r in
    kilofeet per cube-root kiloton: Eq.42)
    and is valid for 10^(-3) < T < 26,000 ms at 1 KT. This more complex fit is
    advisable for pressures above 10,000 psi, or scaled times less than
    0.6 ms/KT^(1/3).
    For a surface burst, use m = (2W)^3.
    "
    """

    a = (0.543 - 21.8 * r + 386 * r**2 + 2383 * r**3) * r**8
    b = (2.99e-8 - 1.91e-4 * r**2 + 1.032 * r**4 - 4.43 * r**6) * 1e-6
    c = (1.028 + 2.087 * r + 2.69 * r**2) * r**8

    return a / (b + c)  # Eqn.(41)


def _s_D_p_pos_to_tau(tau):
    """
    scaled positive overpressure duration versus time of arrival

    Arguments:
    - tau, scaled time of arrival in milliseconds per square root kiloton
    Returns:
    - scaled D_p_pos, time of positive overpressure duration in ms per cube
    root killton

    "
    The duration of the positive phase for overpressure can be expressed
    variously as a function of overpressure, shock range, or time of arrival
    in a fit to detailed one-dimensional blast calculations [Brode, 1959b]
    "
    """

    return (813140 + 11412 * tau + 313 * tau**2) / (6780 + 444.7 * tau + tau**2)


def _s_D_p_pos_to_DeltaP_s(DeltaP_s):
    """
    scaled positive overpressure duration versus peak overpressure

    Arguments:
    - DeltaP_s, peak overpressure
    Returns:
    - D_p_pos, scaled time of positive overpressure duration in ms per square
    root killton

    "
    The duration of the positive phase for overpressure can be expressed
    variously as a function of overpressure, shock range, or time of arrival
    in a fit to detailed one-dimensional blast calculations [Brode, 1959b]

    This form, in combination with Eq. (33), is compared with both the
    calculation to which it is a fit [Brode, 1959b] and the DNA 1-KT standard
    values [Needham and Crepeau, 1981]. The difference between those two
    predictions is typical of the uncertainty and expected variation
    (particularly at large distances) in blast parameters such as positive phase
    durations. For a surface burst, this form should be multiplied by
    2^(1/3) = 1.26 (i.e., use 2W in m).
    "
    """

    _pi = DeltaP_s / 1000  # in ksi

    return (
        -148.6
        + 497.3 / (1 + 18.68 * _pi**0.6783)
        + 1629 * _pi**0.8711 / (1 + 6.477 * _pi**0.8555)
    )  # Eqn.(45), in scaled


def _s_D_p_pos_to_r(r):
    """
    scaled positive overpressure duration versus scaled range

    Arguments:
    - r, scaled ground range for kilofeet per cube-root kiloton
    Returns:
    - scaled D_p_pos, scaled time of positive overpressure duration in ms per square
    root killton

    "
    The duration of the positive phase for overpressure can be expressed
    variously as a function of overpressure, shock range, or time of arrival
    in a fit to detailed one-dimensional blast calculations [Brode, 1959b]

    Again, for a surface burst, Dp should be increased by a factor of 2^(1/3),
    i.e., the free-air form for twice the yield.
    "
    """
    return (
        69.12
        + 46.19 / (1 + 3e6 * r**7.217)
        + 4043 * r**6.329 / (1 + 37.16 * r**5.621)
    )  # Eqn.(46), in scaled


def _s_I_p_pos(DeltaP_s):
    """
    overpressure impulse in positive phase versus peak overpressure

    Arguments:
    - DeltaP_s: peak overpressure in psi
    Returns:
    - scaled I_p_pos, scaled impulse in psi-ms per cube-root kiloton

    "
    This fit, when used with Eq. (33) leads to the values compared in Fig. 16
    versus radius. The impulse approaches a constant at small ranges and decays
    approximately as the inverse of the range elsewhere. For a surface burst,
    that expression should be multiplied by 1.26 (i.e., by 2^(1/3)), which leads
    to replacing the coefficient 145 by 183. This form is good to better than
    10 percent for 2 < APQ < 100 000 psi. Comparison between the approximation
    as a function of peak overpressure and the detailed numerical results
    [Brode, 1964] is made in Fig. 17. In that plot, it is evident that impulse
    increases approximately as the cube root of overpressure, but tends toward
    a constant at about a few thousand pounds per square inch.
    "
    """

    return (
        145 * DeltaP_s**0.5 / (1 + 0.00385 * DeltaP_s**0.5)
    )  # Eqn. (48) in scaled formulation


def _DeltaP(t, tau, m, s_D_p_pos=None):
    """
    Overpressure versus time

    Arguments:
    - t: arbitrary time after detonation, t >= T, in milliseconds
    - tau: scaled time of arrival in milliseconds per square root kiloton
    - m: yield scaling factor, cube-root kiloton
    - s_D_p_pos: scaled overpressure positive phase duration in ms/kT^(1/3)

    Returns:
    - scaled DeltaP: overpressure at time t in psi per cube-root kiloton
    "
    The following analytic expression is valid for overpressures less than
    about 15,000 psi. It is an approximate form, modified from earlier fits
    [Brode, 1970, 1978] for the overpressure in the positive phase as a function
     of time. In these fits, time is zero at the instant of burst: Eqn.(49)
    "
    """

    sigma = t / m
    DeltaP_s = _DeltaP_s_to_tau(sigma)  # yes

    if s_D_p_pos is None:  # any of the Eqn.(43) to (46) could be used here
        s_D_p_pos = _s_D_p_pos_to_tau(tau)

    return (
        DeltaP_s
        * (
            0.417
            + 0.583
            * (tau / sigma) ** 6
            * (40 * (tau / sigma) ** 6 + tau**2)
            / (40 + tau**2)
        )
        * (1 - (sigma - tau) / s_D_p_pos)
    )  # in Eqn.(49)


"""
Equation (50) is a simplified version of Eqn (63), see that file
"""


def _s_D_u_pos(DeltaP_s):
    """
    scaled duration of the outward dynamic pressure versus peak overpressure

    Arguments:
    - DeltaP_s, peak overpressure

    Returns:
    - scaled D_u_pos, in psi-ms per cube root kiloton
    "
    For a surface burst, substitute 2W for W, i.e., m ~ 1.26m.
    The fit is in good agreement with earlier calculations [Brode, 1959b, 1966],
    but it differs significantly from the 1-KT standard [Needham and Crepeau, 1981].
    The difference is illustrated in Fig. 18 (scaled to 1-MT surface burst).
    Such a disparity between results of detailed numerical (one-dimensional)
    calculations is a measure of the differences introduced by dissimilarities
    in boundary and initial conditions, equations of state and opacities, and by
    various treatments of radiation transport, thermal radiation losses, and
    accumulated numerical errors in detailed computer calculations.
    "
    """

    _pi = DeltaP_s / 1000  # ksi

    return (
        317 / (1 + 85 * _pi + 7500 * _pi**2)
        + 6100 * _pi / (1 + 420 * _pi**2)
        + 2113 * _pi / (1 + 11 * _pi)
    )  # Eqn. (52)


def _s_I_u_pos_to_DeltaP_s(DeltaP_s):
    """
    scaled dynamic pressure positive phase impulse impulse

    Arguments:
    - DeltaP_s, peak overpressure

    Returns:
    - scaled I_u_pos, in psi-ms per cube root kiloton

    "
    That form is within 10 percent of the scaled values from the detailed
    calculations [Brode, 1959b, 1966] for 3 < ΔP_s < 10,000 psi. It is high by
    nearly 20 percent at ΔP_s = 100,000 psi.
    "
    """

    return 2.14 * DeltaP_s**1.637 / (1 + 0.00434 * DeltaP_s**1.431)  # Eqn. (54)


def _s_I_u_pos_to_r(r):
    """
    scaled dynamic pressure positive phase impulse impulse

    Arguments:
    - r, scaled range in kft per cube root kiloton

    Returns:
    - scaled I_u_pos, in psi-ms per cube root kiloton

    "
    A fit to the dynamic impulse versus range for the early calculations
    [Brode, 1959b] agrees to better than 10 percent for 0.0025 <= r <= 2 kft/KT^(1/3).
    The relation, when scaled to a 1-KT free-air burst, is:

    This expression is compared with the detailed calculation results to which
    it was fit. The fit is good to a few percent over the entire range. For
    a surface burst, m = (2W)^(1/3).
    "
    """

    return (
        18.8 * r**2 / (1e-6 + 0.06896 * r**3 + 0.5963 * r**5.652)
        + 92.64 / (100 * r) ** 5
        + 2935
        * (r - 0.00597)
        * (0.01 - r)
        * (0.0003552 - r**4)
        / (1e-10 + 0.003377 * r**2.5 + 155.8 * r**8)
    )  # Eqn.(55)


def _Q_s_over_Delta_P_s(DeltaP_s):
    """
    In an ideal atmosphere, this would be the hugoniot equation

    Arguments:
    - DeltaP_s, peak overpressure

    Returns:
    - Q_s, peak dynamic pressure
    "
    The above expression is accurate to better than 3 percent for all
    overpressures less than 7,000,000 psi.
    "
    """

    _pi = DeltaP_s / 1000  # in ksi
    zeta = _pi / 1000  # lol wut

    return (
        _pi
        * (1 + 0.241 * _pi + 0.4376 * _pi**2)
        / (0.041 + 0.4 * _pi + 0.02891 * _pi**2 + 0.1015 * _pi**3)
        + 0.01251 * _pi**2 / (1 + 9.649e-7 * _pi**5)
        + 7.29e-8 * _pi**4 / (1 + 2.61e-21 * _pi**12)
        + 9.763e-10 * _pi**4 / (1 + 6.957e-28 * _pi**12)
        - 5.052e-8 * _pi**6 / (1 + 1.368e-14 * _pi**12)
        - 6.021e5 * zeta**6 / (1 + 3.541e12 * zeta**12)
        - 2.17e8 * zeta**14 / (1 + 1.62e9 * zeta**15)
        - 0.7670 * zeta**2.839 / (1 + 0.1646 * zeta**3.678)
    )  # Eqn.(18)


def _Q(t, tau, m, DeltaP_s=None, Q_s=None, D_u_pos=None):
    """
    dynamic pressure over time

    Arguments:
    - t, time, t>=T
    - tau, scaled time of arrival
    - m, scaling factor
    - (optionally) DeltaP_s, peak overpressure
    - (optionally) Q_s, peak dynamic pressure
    - (optionally) D_u_pos, dynamic pressure duration

    Returns:
    - Q, dynamic pressure in psi

    "An older approximate analytic expression for dynamic pressure versus time
    covers the range 2 <= ΔP_x (??) <= 1000 psi (0.1 < Q_s < 3000 psi) [Brode, 1964]
    "
    """
    if DeltaP_s is None:
        DeltaP_s = _DeltaP_s_to_tau(tau)

    if Q_s is None:
        Q_s = _Q_s_over_Delta_P_s(DeltaP_s) * DeltaP_s

    if D_u_pos is None:
        D_u_pos = _s_D_u_pos(DeltaP_s) * m

    _pi = DeltaP_s / 1000

    T = tau * m  # scale it back

    w = (t - T) / D_u_pos

    d = 1.06 * _pi**0.035 / (1 + 147 * _pi**3) + 2.13 * _pi**3 / (
        1 + 67.9 * _pi**3.5
    )
    a = 0.38 * DeltaP_s**0.8605
    b = 5.4 * DeltaP_s**0.604

    return Q_s * (1 - w) ** 2 * (d * exp(-a * w) + (1 - d) * exp(-b * w))  # Eqn.(56)


def _theta_and_t_m(DeltaP_s):
    """

    Arguments:
    - DeltaP_s, peak overpressure
    Returns:
    - theta_m, maximum temperature, in 10^3 deg C
    - t_m, scaled time to maximum temperature, in milliseconds in square root kiloton

    "
    Another parameter of interest is the maximum temperature at a given range.
    The maximum temperature at overpressure levels above 100 psi occurs after
    shock arrival, as hotter fireball air expands past that point. That maximum
    temperature is somewhat dependent on yield, since the fireball of a megaton
    burst cools more slowly than that of a kiloton burst (even if cube-root
    scaling is applied). The rough fit of Eq. (59) is appropriate for a megaton
    and is based on earlier calculations by Brode [1959b]: Eqn. (59)

    More recently, improved opacities for air have led to slightly higher fireball
    temperatures at corresponding times, so that temperatures somewhat above
    those predicted by Eq. (59) are likely.

    The time to maximum temperature cannot be rigorously scaled; but assuming
    cube-root scaling, that time is roughly inversely proportional to the peak
    overpressure (free-air burst): Eqn. (60)
    "

    It is observed that calculating the maximum of the two values results in a
    closer fit in all range (see validation)
    """

    _pi = DeltaP_s / 1000  # in ksi

    theta_m = 1090 * _pi**3.26 / (1 + 35.6 * _pi**2.75)  # Eqs. (59)
    # this formulation works at high pressure, drops off too low at lower ones

    theta_s = (
        5.310
        * _pi
        * (1 + 34 * _pi)
        / (1 + 58.3 * _pi + 6.53 * _pi**3.5 / (1 + 0.2027 * _pi**2))
    )  # Eqs. (23) modified by dividing by 1000
    # this formulation works at lower pressure, drops off too low at high pressure

    t_m = 19 / _pi

    return max(theta_m, theta_s), t_m


def _D_p_neg(m):
    return 1051 * m  # duration of negative phase in ms


def _DeltaP_neg(t, m, DeltaP_s, D_p_neg):
    """
    time history for negative overpressure, or underpressure
    Arguments:
    - t, some time after the beginning of negative pressure, in milliseconds
    - m, yiled scaling factor in square root kiloton

    Returns:
    - DeltaP_neg, underpressure in psi

    """

    t_n = (
        151.4 + 2844 / DeltaP_s**0.9638
    ) * m  #  time of beginning of negative phase

    An = 0.2532 * DeltaP_s / (1 + 0.1262 * DeltaP_s) + 413.2 * (DeltaP_s / 100) ** 4 / (
        1 + 668.1 * (DeltaP_s / 100) ** 5
    )
    Bn = 2.481 * DeltaP_s / (1 + 0.004272 * DeltaP_s**1.7)
    Cn = 18.55 * (DeltaP_s / 100) ** 8 / (1 + 2.75 * (DeltaP_s / 100) ** 7.335)

    tau = (t - t_n) / D_p_neg

    Po = 14.7  # we use a fixed value for the ambient here

    return -Po * An * tau * (1 - tau) / (1 + Bn * tau**2 + Cn * tau**3)


def freeair(R, Y, t=None, prettyPrint=True):

    m = Y ** (1 / 3)
    r = _uc_m2ft(R) / m / 1000
    DeltaP_s = _DeltaP_s_to_r(r)
    Q_s = _Q_s_over_Delta_P_s(DeltaP_s) * DeltaP_s
    tau = _ci_tau_to_r(r)

    D_p_pos = _s_D_p_pos_to_DeltaP_s(DeltaP_s) * m  # in ms
    I_p_pos = _s_I_p_pos(DeltaP_s) * m  # in ms
    D_u_pos = _s_D_u_pos(DeltaP_s) * m  # in ms
    I_u_pos = _s_I_u_pos_to_r(r) * m  # in ms

    theta_m, t_m = _theta_and_t_m(DeltaP_s)

    D_p_neg = _D_p_neg(m)

    T = tau * m / 1000  # in seconds

    if t is not None:

        if t < T:
            OP = "Not Yet Arrived"
            DeltaP_t = 0

        elif t < T + D_p_pos / 1000:
            OP = "Positive"
            DeltaP_t = _DeltaP(t * 1000, tau, m, D_p_pos / m)

        elif t < T + D_p_pos / 1000 + D_p_neg / 1000:  # in seconds
            OP = "Negative"
            DeltaP_t = _DeltaP_neg(t * 1000, m, DeltaP_s, D_p_neg)

        else:
            OP = "Passed"
            DeltaP_t = 0

        if t < T:
            DP = "Not Yet Arrived"
            Q_t = 0

        elif t < T + D_p_pos / 1000:  # in seconds
            DP = "Positive"
            Q_t = _Q(t * 1000, tau, m, DeltaP_s, Q_s, D_u_pos)
            # Q_t = _Q(t, tau, m)
        else:
            DP = "Negative/Passed"
            Q_t = None

    else:
        DeltaP_t = None
        Q_t = None

    if DeltaP_t is not None:
        DeltaP_t = _uc_psi2pa(DeltaP_t)
    if Q_t is not None:
        Q_t = _uc_psi2pa(Q_t)

    t_m *= m / 1000  # scaling as well as ms-> s
    theta_m *= 1e3  # convert to deg C

    DeltaP_s = _uc_psi2pa(DeltaP_s)
    Q_s = _uc_psi2pa(Q_s)
    I_p_pos = _uc_psi2pa(I_p_pos / 1000)
    I_u_pos = _uc_psi2pa(I_u_pos / 1000)

    D_p_pos /= 1000  # ms to s
    D_u_pos /= 1000  # ms to s
    D_p_neg /= 1000  # ms to s

    if prettyPrint:
        print("{:^49}".format("INPUT"))
        print("Range{:.>16,.6g} m Yield{:.>16,.6g} m ".format(R, Y))

        print(
            "Time To:{:.>37} s".format("{:,.6g}".format(t) if t is not None else "N/A")
        )

        print("{:^49}".format("OUTPUT"))
        print("Time of Arrival:{:.>29,.6g} s".format(T))

        print("")
        print("TOTAL{:>19}{:>24} ".format("Overpressure", "Dynamic Pressure"))
        print("Peak:{:.>16,.6g} Pa{:.>21,.6g} Pa ".format(DeltaP_s, Q_s))
        print("Duration:{:>12,.6g} s{:>22,.6g} s ".format(D_p_pos, D_u_pos))
        print("Impulse:{:.>13,.6g} Pa-s{:.>19,.6g} Pa-s".format(I_p_pos, I_u_pos))

        if t is not None:
            print("")
            print("PARTIAL{:>17}{:>24} ".format("Overpressure", "Dynamic Pressure"))
            print("Phase:{:.>18}{:.>24}".format(OP, DP))
            print(
                "Value:{:.>15,.6g} Pa{:.>21} Pa".format(
                    DeltaP_t, "{:,.6g}".format(Q_t) if Q_t is not None else "N/A"
                )
            )

    return (DeltaP_s, DeltaP_t, Q_t, T, D_p_pos, I_p_pos, D_u_pos, I_u_pos, theta_m)


if __name__ == "__main__":
    from HeWu.test import runFAtest

    def _freeair_test(R, Y):
        (DeltaP_s, _, _, T, D_p_pos, I_p_pos, D_u_pos, I_u_pos, theta_m) = freeair(
            R, Y, None, False
        )
        return DeltaP_s, T, D_p_pos, I_p_pos, D_u_pos, I_u_pos, theta_m

    runFAtest(_freeair_test)

    pass
