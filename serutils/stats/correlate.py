# References:
# - Statistics in Geography by David Ebdon (ISBN: 978-0631136880)
# - Reliability Engineering Resource Website:
# - http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
# - University of Glascow, Department of Statistics:
# - http://www.stats.gla.ac.uk/steps/glossary/confidence_intervals.html#conflim

import re

import numpy as np

from scipy import odr
from scipy.stats import t as t_distr
from scipy.optimize import curve_fit


def get_confidence(x, y, p_y, new_x, df, conf=0.95):
    """
    Computes prediction band and confidence band from distribution and it's fit

    :params x: data array or list
    :params y: data array or list
    :params p_y: predicted data array corresponding to x
    :params new_x: data array to be predicted
    :params None df: number of degrees of freedom needed to fit
    :params .95 conf: desired confidence level, by default 0.95 (2 sigma)

    :returns: an array with the confidence values to add/subtract, another
       with prediction values, and the R-square value of the fit.
    """
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if not isinstance(y, np.ndarray):
        y = np.array(y)
    if not isinstance(p_y, np.ndarray):
        p_y = np.array(p_y)
    if not isinstance(new_x, np.ndarray):
        new_x = np.array(new_x)
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
    sdi = (new_x - mean_x)**2
    # standard deviation
    sd = sum((x - mean_x)**2)
    # relative individual deviation
    sdi_sd = sdi / sd

    confs = t * np.sqrt(mse * (      1.0 / n + sdi_sd))
    preds = t * np.sqrt(mse * (1.0 + 1.0 / n + sdi_sd))

    # calculate R-square
    mean_y = np.mean(y)
    sst = sum((y - mean_y)**2)
    r2 = 1 - sse / sst

    return confs, preds, r2


def build_func(func_string, df=None):
    recomp = re.compile("[A-Z]")
    func_restring = re.sub(recomp, "%s", func_string)

    df = len(re.findall(recomp, func_string)) if df is None else df

    def func(x, *args):
        if USE_ODR:
            args, x = tuple(x), np.array(args)
        cmd = "zzz = " + func_restring.replace('^', '**') % (args)

        exec(cmd) in globals(), locals()

        try:
            return np.lib.asarray_chkfinite(zzz)
        except:
            # avoid the creation of NaNs when invalid values for power or log
            return x
    return func, df


def fit_function(x, y, func, df=None, **kwargs):
    """
    Fit a given function using input data points

    :params x: data array or list
    :params y: data array or list
    :params func: a function.
    :params None df: polynomial degree for the fit

    :returns: fitted parameters
    """

    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if not isinstance(y, np.ndarray):
        y = np.array(y)

    if USE_ODR:
        # Fit using orthogonal distance
        quad_model = odr.Model(func)
        data = odr.RealData(x, y)
        odr_obj = odr.ODR(data, quad_model, beta0=[1. for _ in xrange(df)])
        out = odr_obj.run()
        z = out.beta
        _ = out.sd_beta
    else:
        # fit a curve to the data using a least squares of "df" order polynomial fit
        z, _ = curve_fit(func, x, y, [1. for _ in xrange(df)], **kwargs)

    return z


def fit_with_uncertainty(x, y, func_string='A*x+B', df=None, conf=0.95,
                         use_odr=False, x_range=(None, None), precision=100,
                         **kwargs):
    """
    Calculates the confidence and prediction band of the polynomial regression
    model at the desired confidence level.
    The 2 sigma confidence interval is 95% sure to contain the best-fit
    regression line.
    The 2 sigma prediction interval contains 95% of the data points.

    :params x: data array or list
    :params y: data array or list
    :params "A*x+B" func_string: a string representing the function to optimize.
      E.g.: "A*x^B+C*x" constants in upper-case.
    :params None df: polynomial degree for the fit
    :params .95 conf: desired confidence level, by default 0.95 (2 sigma)
    :param None x_range: range in X in which to compute the fit.
    :param 100 precision: number of points used for curve fitting.
    :param False use_odr: if True, uses orthogonal distance regression
       (scipy.ODR), to fit the input function. Returned confidence interval is
       based on the prediction errors of betas, and is to be applied on X and Y\
       axis.

    :returns: the number of degrees of freedom, and a latex formatted string
       representing the formula with fitted values, an array with the fitted x
       values, an array with the fitted y values and the scipy.poly1d object,
       the R-square value of the fit.

    Example:
    --------

    ```
    from random import random

    from matplotlib.patches import Rectangle
    from matplotlib import pyplot as plt
    import numpy as np

    from serutils.stats.correlate import fit_with_uncertainty, latex_formula


    # fake data
    size = 50
    x = [i+(2*random()+.5)**4 for i in xrange(size)]
    y = [np.log((1.+i))+np.log(1+i**(float(i+1)/(random()*20+i))) for i in xrange(size)]
    x = np.array(x)
    y = np.array(y)

    func_string = "A^(x*B)+C"
    # or
    func_string = "A*x^3+ B*x^2+C*x+D"

    p_x, p_y, z, confs, preds, r2 = fit_with_uncertainty(x, y, func_string)
    formula = latex_formula(func_string, z)

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
                '''Fit ($R^2=%.3f$):\ny = $%s$''' % (r2, formula),
                '95% Confidence band',
                '95% Prediction band'],
               loc='upper left', frameon=False, bbox_to_anchor=[1,1])

    plt.xlim((min(x), max(x)))
    plt.subplots_adjust(right=0.65)
    plt.show()
    ```

    """
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if not isinstance(y, np.ndarray):
        y = np.array(y)

    global USE_ODR
    USE_ODR = use_odr
    # build function
    func, df = build_func(func_string, df=df)

    # find optimal parameters of the funciton
    z = fit_function(x, y, func=func, df=df, USE_ODR=USE_ODR, **kwargs)

    # create series of new test x values to predict for
    new_x = np.linspace(min(x) if x_range[0] is None else x_range[0],
                      max(x) if x_range[1] is None else x_range[1], precision)

    # now predict y based on test x-values
    USE_ODR=False
    new_y = func(new_x, *z)

    # estimate confidence and prediction bands
    p_y = func(x, *z)
    confs, preds, r2 = get_confidence(x, y, p_y, new_x, df, conf=conf)

    return new_x, new_y, z, confs, preds, r2


def latex_formula(func_string, z):
    recomp = re.compile("[A-Z]")
    func_restring = re.sub(recomp, "%s", func_string)
    formula = re.sub('\^(\([^)(]+\))', '^{\\1}',
                     re.sub('\^([0-9.]+)', '^{\\1}',
                            func_restring.replace('np.', '')))
    formula = formula % tuple(["%.2g" % i for i in z])
    formula = re.sub('\+ *-', '-', formula.replace('*',''))
    formula = re.sub('e([-+][0-9]+)', 'e^{\\1}', formula)
    return formula
