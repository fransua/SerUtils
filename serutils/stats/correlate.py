# References:
# - Statistics in Geography by David Ebdon (ISBN: 978-0631136880)
# - Reliability Engineering Resource Website:
# - http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
# - University of Glascow, Department of Statistics:
# - http://www.stats.gla.ac.uk/steps/glossary/confidence_intervals.html#conflim

import re

import numpy as np

from scipy.stats import t as t_distr
from scipy.optimize import curve_fit


def fit_with_uncertainty(x, y, func_string='A*x+B', df=None, conf=0.95,
                         x_range=(None, None), precision=100):
    """
    Calculates the confidence and prediction band of the polynomial regression
    model at the desired confidence level.
    The 2 sigma confidence interval is 95% sure to contain the best-fit
    regression line.
    The 2 sigma prediction interval contains 95% of the data points.

    :params .95 conf: desired confidence level, by default 0.95 (2 sigma)
    :params x: data array or list
    :params y: data array or list
    :params "A*x+B" func_string: a string representing the function to optimize.
      E.g.: "A*x^B+C*x" constants in upper-case.
    :params None df: polynomial degree for the fit
    :param None x_range: range in X in which to compute the fit.
    :param 100 precision: number of points used for curve fitting.

    :returns: an array with the confidence values to add/subtract, another
       with prediction values, an array with the fitted x values, an array with
       the fitted y values and the scipy.poly1d object, the R-square value of
       the fit, and a latex formatted string representing the formula with
       fitted values.

    Example:
    --------

    ```
    from random import random

    from matplotlib.patches import Rectangle
    from matplotlib import pyplot as plt
    import numpy as np

    from serutils.stats.correlate import fit_with_uncertainty


    # fake data
    size = 50
    x = [i+(2*random()+.5)**4 for i in xrange(size)]
    y = [np.log((1.+i))+np.log(1+i**(float(i+1)/(random()*20+i))) for i in xrange(size)]
    x = np.array(x)
    y = np.array(y)

    func_string = "A^(x*B)+C"
    # or
    func_string = "A*x^3+ B*x^2+C*x+D"

    confs, preds, p_x, p_y, z, r2, formula = fit_with_uncertainty(x, y, func_string)
    # plot sample data
    dots = plt.plot(x, y, 'o', mec='none', color='red', label='Sample observations', alpha=.8)

    # plot line of best fit
    fit_line = plt.plot(p_x, p_y,color= 'darkgreen', lw=2, label='Regression line')

    # plot confidence limits
    plt.fill_between(p_x, p_y + confs, p_y - confs, color='darkgreen', alpha=0.3)
    plt.fill_between(p_x, p_y - preds, p_y + preds, color='orangered', alpha=0.2)

    p1 = Rectangle((0, 0), 1, 1, fc="darkgreen", alpha=.3)
    p2 = Rectangle((0, 0), 1, 1, fc="orangered", alpha=.2)
    plt.legend(dots + fit_line + [p1, p2],
               ['almost random dots',
                'Fit ($R^2=%.3f$):\ny = $%s$' % (r2, formula),
                '95% Confidence band',
                '95% Prediction band'],
               loc='upper left', frameon=False, bbox_to_anchor=[1,1])

    plt.xlim((min(x), max(x)))
    plt.subplots_adjust(right=0.65)
    plt.show()
    ```

    """
    recomp = re.compile("[A-Z]")
    func_restring = re.sub(recomp, "%s", func_string)

    df = df or len(re.findall(recomp, func_string))

    def func(x, *args):
        cmd = "zzz = " + func_restring.replace('^', '**') % (args)
        exec(cmd) in globals(), locals()
        #print cmd
        try:
            return np.lib.asarray_chkfinite(zzz)
        except:
            # avoid the creation of NaNs when invalid values for power or log
            return x

    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if not isinstance(y, np.ndarray):
        y = np.array(y)

    # fit a curve to the data using a least squares of "df" order polynomial fit
    z, covariance = curve_fit(func, x, y, [1. for _ in xrange(df)])
    # predict y values of original data using the fit
    p_y = func(x, *z)

    # create series of new test x values to predict for
    p_x = np.linspace(min(x) if x_range[0] is None else x_range[0],
                      max(x) if x_range[1] is None else x_range[1], precision)

    # number of samples in original fit
    n = x.size
    # alpha 1 minus the wanted probability
    alpha = 1. - conf
    # t distribution with n-2 degrees of freedom
    t = t_distr.ppf(1. - alpha / 2., n - df - 1)

    # mean of x
    mean_x = np.mean(x)
    # Error sum of squares
    sse = sum((y - p_y)**2)
    # Error mean square (estimate of the variance)
    mse = sse / (n - df - 1)
    # Square individual deviation
    sdi = (p_x - mean_x)**2
    # standard deviation
    sd = sum((x - mean_x)**2)
    # relative individual deviation
    sdi_sd = sdi / sd

    confs = t * np.sqrt(mse * (      1.0 / n + sdi_sd))
    preds = t * np.sqrt(mse * (1.0 + 1.0 / n + sdi_sd))

    # calculate R-square
    sstot = sum((y - np.mean(y))**2)
    sserr = sum((y - func(x, *z))**2)
    r2 = 1 - sserr / sstot

    # now predict y based on test x-values
    p_y = func(p_x, *z)

    formula = re.sub('\^(\([^)(]+\))', '^{\\1}',
                     re.sub('\^([0-9.]+)', '^{\\1}',
                            func_restring.replace('np.', '') % tuple(
                                ["%.2g" % i for i in z])))
    formula = re.sub('\+ *-', '-', formula.replace('*',''))
    formula = re.sub('e([-+][0-9]+)', 'e^{\\1}', formula)

    return confs, preds, p_x, p_y, z, r2, formula
