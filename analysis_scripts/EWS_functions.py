import numpy as np
import statsmodels.api as sm
import scipy.stats as st
from scipy.optimize import curve_fit


def runmean(x, w):
    ## running mean of timeseries x with window size w
    n = x.shape[0]
    xs = np.zeros_like(x)
    for i in range(w // 2):
        xs[i] = np.nanmean(x[: i + w // 2 + 1])
    for i in range(n - w // 2, n):
        xs[i] = np.nanmean(x[i - w // 2 + 1:])

    for i in range(w // 2, n - w // 2):
        xs[i] = np.nanmean(x[i - w // 2 : i + w // 2 + 1])
    return xs

## EWS functions for timeseries without missing values ##

def runstd(x, w): 
    ## calculate standard deviation of timeseries x in running windows of length w
    n = x.shape[0]
    xs = np.zeros_like(x)
    for i in range(w // 2): # for beginning edge bit, the window centre starts at w/2
        xw = x[: i + w // 2 + 1] # window always starts at 0 and more and more covers into range
        xw = xw - xw.mean() # make into variations around zero
        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]
            xw = xw - p0 * np.arange(xw.shape[0]) - p1 # remove linear trend

            xs[i] = np.std(xw) # calculate standard deviation
        else:
            xs[i] = np.nan
    for i in range(n - w // 2, n):
        xw = x[i - w // 2 + 1:] # for end bit window always ends at 0 and covers less and less
        xw = xw - xw.mean()
        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]

            xw = xw - p0 * np.arange(xw.shape[0]) - p1
            xs[i] = np.std(xw)
        else:
            xs[i] = np.nan

    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2 : i + w // 2 + 1] # standard window
        xw = xw - xw.mean()
        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]
            xw = xw - p0 * np.arange(xw.shape[0]) - p1

            xs[i] = np.std(xw)
        else:
            xs[i] = np.nan

    return xs

def runac2(x, w):
    ## calcualte the autocorrelation of timeseries x in running window size w
    n = x.shape[0]
    xs = np.zeros_like(x)
    for i in range(w // 2):
        xw = x[: i + w // 2 + 1]
        xw = xw - xw.mean()
        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]
            xw = xw - p0 * np.arange(xw.shape[0]) - p1 # remove linear trend
            xs[i] = sm.tsa.acf(xw)[1]
        else:
            xs[i] = np.nan

    for i in range(n - w // 2, n):
        xw = x[i - w // 2 + 1:]

        xw = xw - xw.mean()
        if np.std(xw) > 0:
            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]

            xw = xw - p0 * np.arange(xw.shape[0]) - p1

            xs[i] = sm.tsa.acf(xw)[1]
        else:
            xs[i] = np.nan

    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2 : i + w // 2 + 1]

        xw = xw - xw.mean()
        if np.std(xw) > 0:

            lg = st.linregress(np.arange(xw.shape[0]), xw)[:]
            p0 = lg[0]
            p1 = lg[1]

            xw = xw - p0 * np.arange(xw.shape[0]) - p1

            xs[i] = sm.tsa.acf(xw)[1]
        else:
            xs[i] = np.nan

    return xs

def run_fit_a_ar1(x, w):
    ## calculate the restoring rate of timeseries x in running windows w
    n = x.shape[0]
    xs = np.zeros_like(x)

    for i in range(w // 2):
        xs[i] = np.nan

    for i in range(n - w // 2, n):
        xs[i] = np.nan

    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2 : i + w // 2 + 1]
        xw = xw - xw.mean() # variations in the window

        p0, p1 = np.polyfit(np.arange(xw.shape[0]), xw, 1)

        xw = xw - p0 * np.arange(xw.shape[0]) - p1 # remove linear trend


        dxw = xw[1:] - xw[:-1]

        xw = sm.add_constant(xw)
        model = sm.GLSAR(dxw, xw[:-1], rho=1)
        results = model.iterative_fit(maxiter=10)

        a = results.params[1]

        xs[i] = a

## EWS functions for timeseries with missing values

def runstd_mv(x, w): 
    n = x.shape[0]
    xs = np.zeros_like(x)
    for i in range(w // 2): 
        xw = x[: i + w // 2 + 1] # window always starts at 0 and more and more covers into range
        xw = xw - np.nanmean(xw) # make into variations around zero

        if np.nanstd(xw) > 0:
            xw = xw - runmean(xw,w) # remove large window running mean instead of linear trend
            xs[i] = np.nanstd(xw) # calculate standard deviation
        else:
            xs[i] = np.nan

    for i in range(n - w // 2, n):
        xw = x[i - w // 2 + 1:] # for end bit window always ends at 0 and covers less and less
        xw = xw - np.nanmean(xw)
        if np.nanstd(xw) > 0:
            xw = xw - runmean(xw,w)
            xs[i] = np.nanstd(xw)
        else:
            xs[i] = np.nan

    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2 : i + w // 2 + 1] # standard window
        xw = xw - np.nanmean(xw)
        if np.nanstd(xw) > 0:
            xw = xw - runmean(xw,w)
            xs[i] = np.nanstd(xw)
    else:
        xs[i] = np.nan

    return xs

def runac_mv2(x, w):
    n = x.shape[0]
    xs = np.zeros_like(x)
    for i in range(w // 2):
        xw = x[: i + w // 2 + 1]
        xw = xw - np.nanmean(xw)
        if np.count_nonzero(np.isnan(xw))<=20:
            if np.nanstd(xw) > 0:
                xw = xw - runmean(xw,w)
                xw1, xw2 = xw[1:], xw[:-1]

                idc = np.where((~np.isnan(xw1))&(~np.isnan(xw2))) # get indices where both x_i and x_(i-1) are not nan
                xs[i] = sm.tsa.acf(xw[idc])[1]
            else:
                xs[i] = np.nan
        else:
            xs[i]=np.nan

    for i in range(n - w // 2, n):
        xw = x[i - w // 2 + 1:]

        xw = xw - np.nanmean(xw)
        if np.count_nonzero(np.isnan(xw))<=20:
            if np.nanstd(xw) > 0:
                xw = xw - runmean(xw,w)
                xw1, xw2 = xw[1:], xw[:-1]

                idc = np.where((~np.isnan(xw1))&(~np.isnan(xw2))) # get indices where both x_i and x_(i-1) are not nan
                xs[i] = sm.tsa.acf(xw[idc])[1]

            else:
                xs[i] = np.nan
        else:
            xs[i]=np.nan

    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2 : i + w // 2 + 1]

        xw = xw - np.nanmean(xw)
        if np.count_nonzero(np.isnan(xw))<=20:
            if np.nanstd(xw) > 0:
                xw = xw - runmean(xw,w)
                xw1, xw2 = xw[1:], xw[:-1]

                idc = np.where((~np.isnan(xw1))&(~np.isnan(xw2))) # get indices where both x_i and x_(i-1) are not nan

                xs[i] = sm.tsa.acf(xw[idc])[1]

            else:
                xs[i] = np.nan
        else:
            xs[i]=np.nan

    return xs

def run_fit_a_ar1_mv(x, w):
    n = x.shape[0]
    xs = np.zeros_like(x)

    for i in range(w // 2):
        xs[i] = np.nan

    for i in range(n - w // 2, n):
        xs[i] = np.nan

    for i in range(w // 2, n - w // 2):
        xw = x[i - w // 2 : i + w // 2 + 1]
        xw = xw - np.nanmean(xw)
        xw = xw - runmean(xw,w)

        dxw = xw[1:] - xw[:-1]

        idc = np.where(~np.isnan(dxw)) # get non-nan indices in difference

        # use these for regression:
        xw2 = xw[idc]
        dxw2 = dxw[idc] 

        xw2 = sm.add_constant(xw2)
        model = sm.GLSAR(dxw2, xw2, rho=1)
        results = model.iterative_fit(maxiter=10)

        a = results.params[1]

        xs[i] = a
    return xs