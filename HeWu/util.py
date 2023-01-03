def intg(f, l, u, tol=1e-3):
    """
    To apply the quadrature procedure, first the problem is transformed on interval to:

    u              1                        given:
    ∫ f(x) dx -> a ∫ f(ax+b) dx             a = (u - l) / 2
    l             -1                        b = (u + l) / 2

    another transformation on the variable of integration eliminates the need to sample
    at either end points, which makes it possible to evaluate improper integrals if the
    asymptotes is at either end point.

    1                                        1
    ∫ f(u) du -- let u = 1.5v-0.5v**3 -> 1.5 ∫ f(1.5v-0.5v^3)*(1-v^2) dv
    -1                                      -1

    as the weight (1-v^2) is exactly 0 on both end points. We then sample evenly along
    v, take quadrature using the mid-point rule and doubling the number of nodes taken
    for each pass. This helps with suppressing harmonics if the function integrated is
    periodic. In addition, all of the previously calcualted quadratures can be reused
    in the next round, after dividing by half. This is especially important when funct-
    ion calls are expensive. Specifically, for pass k (k>=1 & integer) we consider 2^k-1
    points (besides the end points):

    v(i) = -1 + 2^(1-k) * i as i increments from 1 to 2^k-1 (inclusive).

                                   2^k-1
    let I(k) =  2^(1-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2)
                                   i=1

                                     2^k+1
    then I(k+1) = 2^(-k) * 1.5 * a * Σ f(1.5v-0.5v^3)*(1-v^2) for every odd i + I(k)/2
                                     i=1

    as a rough approximation, the error is simply taken to be the change in estimated value
    between two successive evaluations:

    ΔI(k) = I(k) - I(k-1)

    if the quadrature procedure results in a converging result, then the error should dec-
    rease faster than the increment in the result, speaking in absolute terms. Although
    this is no-way guaranteed, it is convenient to take the increment as an upper bound
    on error. Therefore we check for three consecutive increments smaller than the specified
    tolerance before submitting the result as a good enough estimate for the integral.
    """

    tol = abs(tol)  # ensure a positive tolerance

    a = (u - l) / 2
    b = (u + l) / 2

    k = 1  # iteration counter
    I = 0  # integral counter
    c = 0  # trend counter, No. of iterations with reducing delta.

    while c < 3:
        dI = 0  # change to integral
        for i in range(1, 2**k, 2):
            v = -1 + 2 ** (1 - k) * i
            u = 1.5 * v - 0.5 * v**3
            dI += f(a * u + b) * (1 - v**2)

        dI *= 1.5 * a * 2 ** (1 - k)
        d = abs(I * 0.5 - dI)  # delta, change per iteration
        I = I * 0.5 + dI
        k += 1

        if d < tol:
            c += 1
        else:
            c = 0

    return I, d


if __name__ == "__main__":

    def f(x):
        return x**2

    # x^2 <- 1/3 * x^3

    print(intg(f, 0, 1))
