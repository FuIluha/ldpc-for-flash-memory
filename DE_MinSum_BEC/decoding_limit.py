import numpy as np

def decoding_limit() -> list[np.ndarray, np.ndarray]:
    def H2(p):
        """Binary entropy function in bits (base 2), with special handling of 0·log(0) = 0."""
        p = np.asarray(p)
        term1, term2 = np.zeros(p.shape), np.zeros(p.shape)
        term1[p > 0] = -p * np.log2(p[p > 0])
        term2[p < 1] = -(1 - p[p > 0]) * np.log2(1 - p[p > 0])
        return term1 + term2

    def invert_entropy(target, tol=1e-12, max_iter=200):
        """
        Invert H2(p) = target for p in [0, 0.5].
        For target <= 0 -> p=0; for target >= 1 -> p=0.5.
        Uses bisection (H2 is monotone increasing on [0, 0.5]).
        """
        if target <= 0:
            return 0.0
        if target >= 1:
            return 0.5
        a, b = 0.0, 0.5
        for _ in range(max_iter):
            m = 0.5 * (a + b)
            fm = H2(m)
            if abs(b - a) < tol:
                return m
            if fm < target:
                a = m
            else:
                b = m
        return 0.5 * (a + b)

    # Domain where solutions exist: p0 in [0, 0.5]
    p0 = np.linspace(0.0, 0.5, 501)
    # Target entropy values from the equation: (1 - H2(p_-)) * (1 - p0) = 0.5  =>  H2(p_-) = 1 - 0.5 / (1 - p0)
    targets = 1.0 - 0.5 / (1.0 - p0)

    # Invert entropy to get the lower-branch p_- in [0, 0.5], then map to p_minus
    p_minus_lower = np.empty_like(p0)
    for i, t in enumerate(targets):
        p_ = invert_entropy(float(t))
        p_minus_lower[i] = (1.0 - p0[i]) * p_

    # Upper branch is symmetric: p_- -> 1 - p_  => p_minus -> (1 - p0) * (1 - p_)
    return p0, p_minus_lower