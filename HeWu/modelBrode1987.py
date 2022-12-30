"""
Zhai Jinpeng, 翟锦鹏, 2022
Contact: 914962409@qq.com

implements the Brode, 1978 model from the book:
PSR Report 1419-3
"AIRBLAST FROM NUCLEAR BURSTS—ANALYTIC APPROXIMATIONS"
Harold L. Brode, July 1987
Technical Report
CONTRACT No. DNA 001-85-C-0089  

Specifically, the below came from SECTION 4.
"""


def _DeltaP_s(GR, H, W):
    """
    Peak overpressure in psi

    GR: ground range, in kft
    H: height, in kft
    W: yield, kT

    """
    m = W ** (1 / 3)
    x = GR / m / 1000
    y = H / m / 1000

    r = (x**2 + y**2) ** 0.5

    z = y / x

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

    return 10.47 / r**a + b / r**c + d * e / (1 + f * r**g) + h + j / r**k
