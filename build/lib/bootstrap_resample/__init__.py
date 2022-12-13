import pandas as pd
from numpy.linalg import norm
import numpy as np
from warnings import warn
from bootstrap_resample.pmap import pmap
from itertools import combinations_with_replacement

def possible_boots(n):
    """Calculate size of full bootstrap distribution"""
    from scipy.special import comb
    return comb(n+n - 1, n)

def _rs(df): 
    """generate random sample"""
    return df.sample(len(df), replace=True)

def Storey(pvalues, m, interpolate_pi0=True):
    """
    Estimates q-values from p-values using Storey's FDR approach
    Args
    =====
    interpolate_pi0 (default: True): Estimated P_i0 as suggested in Storey and 
    Tibshirani, 2003. For most GWAS this is unnecessary, since pi0 is ~1
    """
    from scipy import interpolate
    assert type(pvalues) == pd.Series, "pvalues must be indexed Series"
    assert(pvalues.min() >= 0 and pvalues.max() <= 1), "p-values should be between 0 and 1"

    if interpolate_pi0 and m > 100:
            # evaluate pi0 for different lambdas
        lam = np.arange(0, 0.9, 0.01)
        pi0 = np.array([(pvalues > lam_i).sum()/(m*(1-lam_i)) for lam_i in lam]) 
            # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = np.clip(interpolate.splev(lam[-1], tck), 0, 1)
    else:
        pi0 = 1.0

    pvalues = pvalues.sort_values()
    qv = pi0 * pvalues
    qv[-1] = min(qv[-1], 1.0)

    for i in range(m-2, -1, -1):
        qv[i] = min(pi0*m*pvalues[i]/(i+1.0), qv[i+1])
    return qv

def FWER(P, method='Bonferroni'):
    """Map array of P-values (P) to Family-Wise Error Rates. 

Number of hypotheses is determined by length of P. `method` determines 
correction approach. Can be either 'Bonferroni' (default), 'Sidak', 'Hochberg', 
'Storey', or None. See Wikipedia for details of each. 
"""
    if method is None:
        return P
    
    methods = dict(
        Bonferroni  = lambda P, m: P*m,
        Sidak       = lambda P, m: 1 - np.exp(np.log(1 - P)*m),
        Hochberg    = lambda P, m: P*(m - P.sort_values()),
        Storey      = Storey)
    try: 
        return methods[method](P, len(P)).clip(upper=1.)
    except KeyError:
        raise ValueError(f"Unknown FWER method: {method}")

@np.vectorize
def pstars(p_value, pstar_values=np.array([1e-4, 1e-3, 1e-2, 0.05])):
    """pstars(p_value) -> str of `*` 
    
Statistical significance thresholds exceeded by p_value (can be scalar or vector).

Parameters:
-----------
pstar_values (default:APA/NEJM thresholds, i.e. [0.05, 0.01, 0.001, 0.0001]) : list or numpy.array
    p-value thresholds defining each increasing `*`. 
"""
    return (len(pstar_values) - np.searchsorted(pstar_values, p_value))*'*'

class sample(object):
    permissible_uncertainty = 0.21

    def bs_estimate(self, *_args):
        """Bootstrap data once and calculate estimator"""
        return self.estimator(_rs(self.df) if self.constrain_level is None else self.gb.apply(_rs), **self.kwargs)

    def __init__(self, df, estimator, N=4000, map=pmap, constrain_level=None, min_pvalue=None, **kwargs):
        """sample(df, estimator) -> bootstrap object

Creates bootstrap sampling distribution of df.apply(estimator). 

Input:
    df  : pandas.DataFrame
        Contains data to analysis with estimator
    
    estimator   : function(1-D array-like) -> float
        Operates on each column of df. Can be any estimator e.g. mean, 
        std, median
    
    N (default:2000) : int
        Number of bootstrap samples to create for sampling distribution. In
        general, 1,825 replications are needed to obtain 1% accuracy of 95% CI.
    
    constrain_level (default: None) : Index identifier (int or str)
        Samples an equal number of rows from each category in *constrain_level*.

    min_pvalue (default: None) : If not None, will set `N` to appropriate value 
        for rejecting null hypotheses up to a chosen minimum P-value. I.e. if 
        you want to reject a null hypothesis beyond 1e-4, you need to draw more 
        bootstraps than you'd need to reject a null hypothesis beyond 0.05.
    
    **kargs will be passed to 'estimator' function.  

Properties of the Bootstrap Method:
    -Assumes independence of samples (a limitation of nearly all CI estimators)
    -Simple
    -Applicable to any estimator
    -More accurate than CI estimators that assume normality when N is 
    asymptotically large. However, N may need to be very large to converge to 
    true CI, rendering the method infeasible. 
"""
        if min_pvalue is not None:
            N = int(np.ceil(4/(min_pvalue*self.permissible_uncertainty**2)))
            self.min_pvalue = min_pvalue
        else:
            self.min_pvalue = 4/(N*self.permissible_uncertainty**2)
        assert N >= 1, "must have at least one bootstrap sample"
        self.df         = df
        self.estimator  = estimator 
        self.constrain_level = constrain_level
        self.map = map

        self.kwargs = kwargs
        self.true_estimate = estimator(self.df, **kwargs)
        
        if constrain_level is not None:
            self.gb = df.groupby(level=constrain_level)
        
        if possible_boots(len(df)) < N:
            print("N > len(df)! --> Calculating *entire* bootstrap distribution")
            self.bootstrap_samples = pd.DataFrame([self.estimator(df.loc[list(idxs)], **kwargs) for idxs in combinations_with_replacement(df.index.tolist(), len(df))]) 
        else:
            self.bootstrap_samples = pd.DataFrame(list(map(self.bs_estimate, range(N))))

    def Jackknife(self):
        """Calculate Jackknife distribution (see Wikipedia)."""
        return np.array([self.estimator(self.df.drop(ix), **self.kwargs) for ix in self.df.index])

    def percentileofscore(self, scores):
        """Calculate percentile of `scores` in the bootstrap distribution."""
        from scipy.stats import percentileofscore
        return (self.bootstrap_samples - scores).apply(percentileofscore, args=(0, 'mean'))
    
    def CI(self, alpha=0.05, bias_corrected=False, for_matplotlib=False):
        """sample.CI(alpha=0.05) -> pandas.df

Estimates Confidence Interval (CI) of any estimator using Bootstrap Method.

Parameters:
    alpha (default:0.05) : float or array-like
        Level of significance or quantiles of sampling distribution to report:
            if float        -> 1 - alpha CI
            if array-like   -> quantiles of sampling distribution
    
    bias_corrected (default: False) : bool
        Compensates for skewness in bootstrap distribution. See Efron, B.
        (1987) "Better Bootstrap Confidence Intervals" for details. 
    
        There are a number of sanity-checks and conditions for applying this correction. 

    for_matplotlib (default:False) : bool
        Output can be fed directly into matplotlib errorbar assignment. 
"""
        quantiles = [alpha/2, 1-alpha/2] 
        
        if quantiles[0] < self.min_pvalue:
            warn("There is insufficient sampling for the {:.2%} Confidence Interval.".format(1 - alpha))

        samples = self.bootstrap_samples
        if bias_corrected == 'auto':
            bias_corrected = len(self.df) < len(samples)    # Not useful when the Jackknife takes longer to calculate than the bootstrap distribution (bias scales with len(df)**-0.5, anyways). 
            if bias_corrected:
                print("Correcting bootstrap distribution bias...")
        
        if bias_corrected and not self.constrain_level:    # Only apply the bias correction if bootstrapping was on individual samples 

            Jackknife = self.Jackknife()
            Theta = Jackknife - Jackknife.mean()
            a = (Theta**3).sum(axis=0) / ( 6*((Theta**2).sum(axis=0)**1.5) )
            z_true = norm.ppf(self.percentileofscore(self.true_estimate)*1e-2)[:, np.newaxis]
            z_q = z_true + norm.ppf(quantiles)
            adjusted_z = z_true + z_q/(1 - a[:, np.newaxis]*z_q)
            adjusted_quantiles = norm.cdf(adjusted_z) 
            
            assert (adjusted_quantiles[:,0] <= 0.5).all() and (adjusted_quantiles[:,1] >= 0.5).all(), "Bias Correction gave absurd adjusted_quantiles."
            CIs = pd.DataFrame([S.quantile(q=quantiles).values for quantiles, (_col, S) in zip(adjusted_quantiles, samples.iteritems()) ], columns=['low', 'high'], index=samples.columns)
        else:
            CIs = samples.quantile(q=quantiles).T
            CIs.columns = 'low', 'high'

        above_true = CIs['low'] > self.true_estimate
        if above_true.any():
            warn("Low CI above True Estimate.")
            CIs.loc[above_true, 'low'] = np.nan
        
        below_true = CIs['high'] < self.true_estimate
        if below_true.any():
            warn("High CI below True Estimate.")
            CIs.loc[below_true, 'high'] = np.nan
        
        if for_matplotlib:
            return np.vstack((self.true_estimate - CIs['low'], CIs['high'] - self.true_estimate)) 
        
        CIs.insert(len(alpha)/2 if hasattr(alpha, '__len__') else 1, 'true', self.true_estimate)
        CIs.columns.name = 'Estimates' 
        return CIs 

    def print_CI(self, alpha=0.05, formatter='.2f', **kargs):
        """Pretty string-format of Confidence Interval."""
        CIs = self.CI(alpha, **kargs)
        return '\n'.join(["{name:<10}: {true: FORMAT} {0:.0%} CI [{low: FORMAT}, {high: FORMAT}]".replace('FORMAT', formatter).format(
                        1 - alpha, **locals()) for name, (low, true, high) in CIs.iterrows()])

    def pscores(self, null_hypothesis=0, two_sided=True):
        """Calculate P-value(s) of `null_hypothesis` (default: 0).

two_sided [def: True] : If False, reports the probability that estimator exceeds the 
    null_hypothesis; otherwise symmetric. 
"""
        pscores = self.percentileofscore(null_hypothesis)*1e-2
        return pd.DataFrame([pscores, 1 - pscores]).min()*2 if two_sided else pscores

    def pTable(self, method='Hochberg', **kwargs):
        """Create a table of P-values using pscores method.

FWER_method [def: 'Hochberg']: Method for calculating Family-Wise Error Rate.
    See FWER function for other options. 

**kwargs: Keyword arguments for pscores method. 
"""
        return self.CI().assign(**{
            'P-value':self.pscores(**kwargs),
            'FWER': lambda df: FWER(df['P-value'], method=method),
            'pstars' : lambda df: pstars(df['FWER'])})

def describe(df, estimator, **kargs):
    """ describe(df, estimator_string) -> Readable CI for statstic_str
"""
    obj = sample(df, estimator, **kargs) 
    return "Bootstrap estimates of {:}:".format(estimator.__name__ if hasattr(estimator, '__name__') else 'estimator') + '\n' + obj.print_CI()

def seaborn_hacker(S, axis=None):
    if axis is None:
        return S[1]
    if axis == 1:
        return S[:, 1]
    if axis == 0:
        return S[1, :]
   
def rsquared(df, levels=0, agg_level=None):
    from statsmodels.api import OLS
    expanded = df.reset_index(levels)
    expanded.index = df.index
    expanded = expanded.dropna()
    condensed = expanded.groupby(level=agg_level).agg(np.mean) if agg_level else expanded
    Y = condensed.pop(condensed.columns[-1])
    result = OLS(Y, condensed).fit()
    return result.rsquared
