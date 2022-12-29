"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

implements the Brode, 1970 model from:
Stephen J. Speicher & Harold L. Brode (1980-1981?)
"Analytical Approximation for Dynamic Pressure Versus Time"

Since it is not apparent how this model made provisions for pressure v.s. time,
with only vague mentions of a 9th power decay that does not lend itself to a 
straightforward implementation, and since the dynamic pressure model is rather
basic anyway, only the following caluclations are implemented:

free air models requires range
air burst models require range and burst height 

time of arrival (free air) (air burst)
max overpressure (free air) (air burst)
max dynamic pressure given max overpressure
max dynamic pressure horizontal component (air burst)

Although provisions for P_0 (atmospheric pressure) is included, its not very 
useful as the overpressure curves themselves must be corrected if the result
for different atmospheric condition is required. The correction factors can 
be sourced from other models, too.
"""


def _t_fa(r, W):
    """
    free-air-burst time-of-arrival in milliseconds
    W: yield, kilotons
    r: slant range, kft
    """

    m = W ** (1 / 3)

    return (
        0.5429 * m**3 - 21.185 * r * m**2 + 361.8 * r**2 * m + 2383 * r**3
    ) / (m**2 + 2.048 * r * m + 2.6872 * r**2)


def _t_a(x, y, W):
    """
    arrival time for bursts near surface
    (surface burst is approximated by supplying 2W, interpolated for cases
     between)
    x: ground range, kilofeet
    y: height of burst, kilofeet
    W: yield, kiloton
    """
    r = (x**2 + y**2) ** 0.5
    if x < y:
        return _t_fa(r, W)
    else:
        return _t_fa(r, W) * y / x + _t_fa(r, 2 * W) * (1 - y / x)


def _DeltaP_fs(t_a, W):
    """
    free-air-burst peak overpressure in psi
    t_a: time of arrival, milliseconds?
    W: yield, kiloton
    """
    m = W ** (1 / 3)
    # return 1.05e6 / (1 + 130 * (t_a / m) ** 1.14)  # alternative version

    return (
        (14843 * m)
        / (0.0135 * m + t_a)
        * (m**2 + 0.6715 * m * t_a + 0.00481 * t_a**2)
        / (m**2 + 1.8836 * m * t_a + 0.02161 * t_a**2)
    )


def _DeltaP_s(x, y, W):
    """
    Peak overpressure in psi
    t_a: time of arrival in milliseconds
    x: ground range in kilofeet
    y: height of burst in kilofeet
    W: yield in kiloton
    """
    # print(t_a)
    t_a = _t_a(x, y, W)

    r = (x**2 + y**2) ** 0.5

    z = y / (x * t_a)
    P = 1.58 * W / r**3 + 5.3 * (W / r) ** 0.5 / r + 0.0215

    A = 0.743 * (1.136 - z) * z**2 / (1.544 + z**6) - 0.0257 * z**6 / (
        0.004435 + z**12
    )
    B = z * (20.42 + 35.5 * z) / (3.57 + z**2) + 2500 * z**4 / (29.3 + z**14)

    C = (
        1
        + z * (2.23 * z - 0.225) / (0.148 + z**2)
        + (28.4 * z**7) / (0.905 + z**7)
    ) ** 3

    E = (
        1
        + 0.002655 * P / (1 + 0.0001728 * P + 1.921e-9 * P**2)
        + (0.004218 + 0.04824 * P + 6.856e-6 * P**2)
        / (1 + 0.008 * P + 3.844e-6 * P**2)
    )

    F = 2.07 * z**2 / (0.00125 + 0.0146 * z**2 + z**8) + 221.25 * z**8 / (
        1 + z**20
    )

    I = 40000 - 17650 * z**2 / (0.235 + z**6)

    H = 1 + A + B * P ** (3 / 2) / (C + P**3) + F * P / (I + P**2)

    a = z**2 * (1 + 2 * z**4) / (1 + 2 * z**6)

    return (
        H
        * (1 + E / (1 + 0.4 / z**4))
        * (a * _DeltaP_fs(t_a, W) + (1 - a) * _DeltaP_fs(t_a, 2 * W))
    )


def _Q_s(DeltaP_s, P_0=14.7):
    """
    Rankine-Hugoniot relation for dynamic pressure
    DeltaP_s: overpressure in psi
    P_0: ambient pressure in psi
    """
    gamma = 1.4  # for OP <= 300 psi
    return DeltaP_s**2 / (2 * gamma * P_0 + (gamma - 1) * DeltaP_s)


def _Q_H(x, y, W, P_0=14.7):
    """
    Horizontal component of dynamic pressure
    x: ground range, kft
    y: burst height, kft
    DeltaP_s: overpressure in psi
    P_0: ambient pressure in psi
    """
    if x < y:
        """in regular reflection region, flow is constrained to horizontal
        velocity only"""
        return _Q_s(_DeltaP_s(x, y, W), P_0) * x / y
    else:
        """in mach reflection region, flow is turned parallel to surface and
        the horizontal component is the total"""
        return _Q_s(_DeltaP_s(x, y, W), P_0)


if __name__ == "__main__":
    y = 0.684
    x = 2.977
    W = 40
    t_a = _t_a(x, y, W)
    print(t_a)
    DeltaP_s = _DeltaP_s(x, y, W)
    print(DeltaP_s)
    Q_H = _Q_H(x, y, W)
    print(Q_H)
