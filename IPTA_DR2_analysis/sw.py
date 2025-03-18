import numpy as np
import argparse, json, os
from h5pulsar import FilePulsar

from dr3_noise.models import model_singlepsr_noise
#from dr3_noise.model_utils import get_freqs, get_flag_groups_by_PTA
#from dr3_noise.pre_processing_utils import make_pulsar_object
from enterprise_extensions.sampler import setup_sampler, group_from_params, get_parameter_groups
import enterprise.constants as const
from enterprise.signals import signal_base, parameter, utils, selections, gp_signals
from enterprise_extensions.blocks import red_noise_block
from enterprise_extensions import chromatic as chrom

@signal_base.function
def sw_dm_triangular_basis(toas, planetssb, sunssb, pos_t, freqs, fref=1400):
    """
    Construct SWGP basis using triangular basis from Nitu et al, 2024
    :param toas: vector of time series in seconds
    :param planetssb: solar system barycenter positions
    :param pos_t: pulsar position as 3-vector
    :param freqs: radio frequencies of observations [MHz]
    :param fref: reference frequency [MHz]

    :return: V: Nc x Ntoa design matrix
    :return: Tc: SW conjunctions
    """

    #First get SW conjunctions
    theta, R_earth, _, _ = chrom.solar_wind.theta_impact(planetssb, sunssb, pos_t)
    # Estimate conjunction from closest approach
    toa_min_theta = toas[np.argmin(theta)]
    #if np.max(theta) < 3:
    #    raise ValueError('No close SW conjunctions. Revise the code!')
    Tc = toa_min_theta + np.arange(100)*const.yr - 50*const.yr
    Tc = Tc[(Tc > np.min(toas))*(Tc < np.max(toas))]

    # Set up triangular basis matrix functions
    Nc = len(Tc)
    Nt = len(toas)
    Lambda = np.max(np.array([1 - np.abs(toas[:,None] - Tc[None,:])/const.yr, np.zeros((Nt,Nc))]), axis=0)

    # Geometric factor (units of DM)
    S_theta = chrom.solar_wind.dm_solar(1.0, theta, R_earth)

    # Convert to a time delay
    S_theta *= 1e-12/const.DM_K/freqs**2

    # Apply DM scaling
    # V = (fref/freqs[:,None])**2*S_theta[:,None]*Lambda

    # return basis and conjunctions
    V = S_theta[:,None]*Lambda
    return V, Tc

@parameter.function
def sw_dm_wn_prior(Tc, log10_sigma_ne=-7):
    """
    Gaussian prior on variance of n_earth at SW conjunctions using Nitu+2024 model
    :param Tc: vector of SW conjunctions
    """
    ne_sigma = 10**log10_sigma_ne
    return ne_sigma**2 * np.ones_like(Tc)

def sw_gp_signal(name='sw_gp'):
    log10_sigma_ne = parameter.Uniform(-4, 2)
    sw_basis = sw_dm_triangular_basis()
    sw_prior = sw_dm_wn_prior(log10_sigma_ne=log10_sigma_ne)
    sw_gp = gp_signals.BasisGP(sw_prior, sw_basis, name=name)
    return sw_gp


elats = {'J0023+0923': 6.309102101039059,
         'J0030+0451': 1.445689005048072,
         'J0034-0534': -8.528452356973881,
         'J0218+4232': 27.0116061250998,
         'J0340+4130': 21.334470118644838,
         'J0437-4715': -67.87296926107412,
         'J0610-2100': -44.417381862577855,
         'J0613-0200': -25.407135074664332,
         'J0621+1002': -13.295007182863769,
         'J0645+5158': 28.85264313780267,
         'J0711-6830': -82.88861704859559,
         'J0751+1807': -2.8075786456756546,
         'J0900-3144': -46.14670241989052,
         'J0931-1902': -31.77671744063506,
         'J1012+5307': 38.75532205517446,
         'J1022+1001': -0.0638937216937713,
         'J1024-0719': -16.044667986871342,
         'J1045-4509': -47.71477467865307,
         'J1455-3330': -16.04479370795446,
         'J1600-3053': -10.071830090669334,
         'J1603-7202': -49.96301036916258,
         'J1614-2230': -1.2567905048204142,
         'J1640+2224': 44.05852415171132,
         'J1643-1224': 9.778334281390634,
         'J1713+0747': 30.700364706239444,
         'J1721-2457': -1.8096246728910548,
         'J1730-2304': 0.18886277506806667,
         'J1732-5049': -27.491601827167283,
         'J1738+0333': 26.88423486352859,
         'J1741+1351': 37.21120358084678,
         'J1744-1134': 11.80520548647091,
         'J1747-4036': -17.201537231988922,
         'J1751-2857': -5.537272946827567,
         'J1801-1417': 9.145627488099496,
         'J1802-2124': 2.0373843795917432,
         'J1804-2717': -3.8564165186072152,
         'J1824-2452A': -1.5487992298720248,
         'J1832-0836': 14.590723616438343,
         'J1843-1113': 11.799994284128172,
         'J1853+1303': 35.743354319738636,
         'J1857+0943': 32.321487755158984,
         'J1903+0327': 25.937987418597157,
         'J1909-3744': -15.155503527773728,
         'J1910+1256': 35.10722753384467,
         'J1911+1347': 35.88643590521716,
         'J1911-1114': 11.088111810821761,
         'J1918-0642': 15.351062753500955,
         'J1923+2515': 46.69620498968441,
         'J1939+2134': 42.296750941090686,
         'J1944+0907': 29.891047892332423,
         'J1949+3106': 50.93091311129789,
         'J1955+2908': 48.684547744242266,
         'J2010-1323': 6.4909513014911076,
         'J2017+0603': 25.044492410900347,
         'J2019+2425': 42.565378323921465,
         'J2033+1734': 35.06285612145178,
         'J2043+1711': 33.96432244869372,
         'J2124-3358': -17.81881083853565,
         'J2129-5721': -39.89996264073268,
         'J2145-0750': 5.3130527751715455,
         'J2214+3000': 37.713147857069536,
         'J2229+2643': 33.29017543516259,
         'J2302+4442': 45.66543169489286,
         'J2317+1439': 17.680224954317637,
         'J2322+2057': 22.878371719047973}