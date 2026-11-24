import math
import numpy as np

def bh_fdr(p_values):
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.ones_like(p_values)
    finite = np.isfinite(p_values)
    if not np.any(finite):
        return adjusted
    indices = np.flatnonzero(finite)
    values = p_values[finite]
    order = np.argsort(values)
    ranked = values[order]
    count = len(ranked)
    corrected = ranked * count / np.arange(1, count + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    adjusted[indices[order]] = np.clip(corrected, 0.0, 1.0)
    return adjusted

def beta_continued_fraction(a, b, x, max_iterations=200, epsilon=3e-14):
    tiny = np.finfo(float).tiny
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = 1.0 / (tiny if abs(d) < tiny else d)
    result = d
    for iteration in range(1, max_iterations + 1):
        twice = 2 * iteration
        coefficient = iteration * (b - iteration) * x / (
            (qam + twice) * (a + twice)
        )
        d = 1.0 + coefficient * d
        d = tiny if abs(d) < tiny else d
        c = 1.0 + coefficient / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        result *= d * c
        coefficient = -(a + iteration) * (qab + iteration) * x / (
            (a + twice) * (qap + twice)
        )
        d = 1.0 + coefficient * d
        d = tiny if abs(d) < tiny else d
        c = 1.0 + coefficient / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < epsilon:
            break
    return result

def regularized_incomplete_beta(x, a, b):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_factor = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    factor = math.exp(log_factor)
    if x < (a + 1.0) / (a + b + 2.0):
        return factor * beta_continued_fraction(a, b, x) / a
    return 1.0 - factor * beta_continued_fraction(b, a, 1.0 - x) / b

def pearson_p_values(correlations, n):
    p_values = np.ones_like(correlations, dtype=float)
    if n < 4:
        return p_values
    a = (n - 2) / 2.0
    for index, correlation in enumerate(correlations):
        if np.isfinite(correlation):
            x = np.clip(1.0 - float(correlation) ** 2, 0.0, 1.0)
            p_values[index] = regularized_incomplete_beta(x, a, 0.5)
    return np.clip(p_values, 0.0, 1.0)

def vectorized_pearson(x, y, calculate_p=True):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.shape[0]
    if n < 4 or np.nanstd(y) == 0:
        return np.full(x.shape[1], np.nan), np.ones(x.shape[1])
    x_centered = x - np.mean(x, axis=0)
    y_centered = y - np.mean(y)
    numerator = np.sum(x_centered * y_centered[:, None], axis=0)
    denominator = np.sqrt(np.sum(x_centered**2, axis=0) * np.sum(y_centered**2))
    correlations = np.divide(
        numerator,
        denominator,
        out=np.full(x.shape[1], np.nan),
        where=denominator > 0,
    )
    correlations = np.clip(correlations, -1.0, 1.0)
    p_values = pearson_p_values(correlations, n) if calculate_p else None
    return correlations, p_values
