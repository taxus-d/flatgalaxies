import numpy as np
import scipy.stats


def linconfband(xx, yy, a, b, conf=0.95, x_grid=None, need_eps=False):
    """
    Confidence bands for linear regression line

    adopted from
    1. https://gist.github.com/rsnemmen/f2c03beb391db809c90f
    2. https://gist.github.com/rsnemmen/0eb32832c657c19e4d39
    3. https://en.wikipedia.org/wiki/Simple_linear_regression
    """
    alpha = 1.0 - conf  # significance
    N = xx.size      # data sample size

    x_grid = np.linspace(xx.min(), xx.max(), 100) if x_grid is None else x_grid
    y_grid = a*x_grid + b

    eps = yy - (a*xx + b)
    sd = np.sqrt(1/(N-2) * np.sum(eps**2))

    sxx = np.sum((xx - xx.mean())**2)
    sx  = (x_grid - xx.mean())**2

    # Quantile of Student's t distribution for p=1-alpha/2
    alpha = alpha/N  # bonferroni correction?
    q = scipy.stats.t.ppf(1.-alpha/2., N-2)
    dy = q * sd * np.sqrt(1./N + sx/sxx)

    ucb = y_grid + dy  # Upper confidence band
    lcb = y_grid - dy  # Lower confidence band
    res = [lcb, ucb]
    if need_eps:
        res.append(eps)
    return res
