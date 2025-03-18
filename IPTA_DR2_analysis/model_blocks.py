import numpy as np
from enterprise.signals import parameter, utils, selections, gp_signals
from IPTA_DR2_analysis.sw import elats, sw_gp_signal
from dr3_noise.models import model_singlepsr_noise


# timespan of IPTA DR2
full_Tspan=952746385.6296968
# timespan of DR2 lite
lite_Tspan=473032183.90760994
# timespan of EDR2 with no legacy data
filtered_Tspan=478180995.23838425


def std_noise_block(psr, full_pta_analysis=False, dataset_Tspan=952746385.6296968, psr_model=True,
                    tm_marg=False, tm_svd=True):
    '''
    Std noise models from IPTA DR2 GWB analysis paper.
    
    The Fourier bases are setup as follows:
        - 1/Tspan ... 30/Tspan freqs RN
        - 1/psr_Tspan ... 30/psr_Tspan freqs for DM
    As such the red noise frequencies will be different for a full PTA run vs a
    single pulsar noise analysis. The full_pta_analysis kwarg
    is used to determine this.
    '''
    psr_Tspan = psr.toas.max() - psr.toas.min()
    if full_pta_analysis:
        Tspan = dataset_Tspan
    else:
        Tspan = psr_Tspan
    
    # WN groups
    if 'NANOGrav' in psr.flags['pta']:
        inc_ecorr = True
        ecorr_groups_by_PTA = {'NANOGrav':'group'}
    else:
        inc_ecorr = False
        ecorr_groups_by_PTA = None
    efeq_groups_by_PTA = {}
    for pta in np.unique(psr.flags['pta']):
        efeq_groups_by_PTA[pta] = 'group'
        
    # custom signals - dip
    if psr.name == 'J1713+0747':
        dm_expdip=True
    else:
        dm_expdip=False
    
    noise = model_singlepsr_noise(psr, Tspan=Tspan, psr_model=psr_model,
                                  # timing
                                  tm_svd=tm_svd, tm_marg=tm_marg,
                                  # red noise
                                  log_A_min=-20, log_A_max=-11, components=30,
                                  # white noise
                                  white_vary=False, tnequad=True, inc_ecorr=inc_ecorr,
                                  efeq_groups=efeq_groups_by_PTA, ecorr_groups=ecorr_groups_by_PTA,
                                  log_equad_min=-10, log_equad_max=-4,
                                  # DM
                                  dm_var=True, dm_type='gp',
                                  dmgp_kernel='diag', dm_psd='powerlaw',
                                  dm_Nfreqs=30, dm_Tspan=psr_Tspan,
                                  # solar wind
                                  dm_sw_deter=False,
                                  # dm dip
                                  dm_expdip=dm_expdip, dm_expdip_basename='exp', dm_expdip_idx=2,
                                  dm_expdip_tau_min=np.log10(5), dm_expdip_tau_max=np.log10(500),
                                  dm_expdip_tmin=54742, dm_expdip_tmax=54768)
    
    return noise


def adv_noise_block(psr, full_pta_analysis=False, dataset=None, psr_model=True,
                    tm_marg=False, tm_svd=True, extra_sigs=False, extra_swgp=False,
                    use_pta_Tspan_for_chrom=False, gp_ecorr=False):
    '''
    Custom noise models. Differences here includes treatment of Fourier bases and 
    addition of custom signals.
    
    The Fourier bases are setup as follows:
        - 1/psr_Tspan ... 30/full_Tspan freqs RN (single pulsar noise)
        - 1/dataset_Tspan ... 30/full_Tspan freqs RN (full PTA/GWB analysis)
        - 1/psr_Tspan ... 150/full_Tspan freqs for DM, chrom, HF RN
    So the Fourier basis will have a different number of frequencies in a single pulsar vs a full
    PTA run, but will always be consistent in the high frequency cutoff. The full_pta_analysis kwarg
    Will use data_Tspan for the lowest bin, otherwise the single psr Tspan is used. A FL analysis
    should use full_pta_analysis=True. This also toggles whether to vary white noise
    
    Additional signals are included for select pulsars:
        - sw_gp for low elat pulsars (using white noise kernel and triangular basis)
        - chrom_gp for select pulsars (using Gaussian prior on chromatic index)
        - exp dip for J1713 (using uniform prior on chromatic index)
        - HF red noise for J1012

    args/kwargs:
        - full_pta_analysis: leave False for single pulsar noise, True for FL or GWB analysis
        - dataset is a required input to set freqs if full_pta_analysis: must be 'full' or 'lite' or 'litec_filtered'
        - psr_model, tm_marg, tm_svd: kwargs for model_singlepsr_noise
        - extra_sigs, extra_sw: Add even more signals to select pulsars (insignificant)
        - use_pta_Tspan_for_chrom: old setting, leave False
    '''
    psr_Tspan = psr.toas.max() - psr.toas.min()
    if full_pta_analysis:
        # set red noise to use full dataset Tspan (extending a bit past ind Tspans)
        if dataset == 'lite':
            Tspan = lite_Tspan
        elif dataset == 'full' or dataset == 'lite_unfiltered':
            Tspan = full_Tspan
        elif dataset == 'litec_filtered':
            Tspan = filtered_Tspan
        else:
            raise KeyError(f'Must specify "full" or "lite" or "litec_filtered" for var dataset (currently {dataset})')
        white_vary = False
        chrom_idx = parameter.Constant()
        dmexp_idx = parameter.Constant()
    else:
        Tspan = psr_Tspan # for red noise and GWB
        white_vary = True
        chrom_idx = parameter.TruncNormal(mu=4, sigma=0.5, pmin=2.5, pmax=10)
        dmexp_idx = parameter.Uniform(pmin=1, pmax=6)
    Nfreqs_low = int(np.ceil(Tspan*30/full_Tspan)) # for red noise and GWB. Is 30 if Tspan = full_Tspan
    Nfreqs_high = int(np.ceil(psr_Tspan*150/full_Tspan)) # for HF (chromatic) noise
        
    chrom_psrs = ['J0437-4715', 'J0613-0200', 'J1600-3053', 'J1643-1224',
                  'J1713+0747', 'J1903+0327', 'J1939+2134']
    hf_psrs = ['J1012+5307']
    sw_psrs = ['J0034-0534', 'J1022+1001', 'J1600-3053', 'J1614-2230', 'J1721-2457',
               'J1751-2857', 'J1843-1113', 'J1918-0642', 'J2124-3358', 'J2145-0750']
    if extra_sigs:
        chrom_psrs += ['J1744-1134']
        hf_psrs += ['J0900-3144']
    
    # WN groups
    if 'NANOGrav' in psr.flags['pta']:
        inc_ecorr = True
        ecorr_groups_by_PTA = {'NANOGrav':'group'}
    else:
        inc_ecorr = False
        ecorr_groups_by_PTA = None
    efeq_groups_by_PTA = {}
    for pta in np.unique(psr.flags['pta']):
        efeq_groups_by_PTA[pta] = 'group'
        
    # custom signals - dip
    if psr.name == 'J1713+0747':
        dm_expdip=True
    else:
        dm_expdip=False
    # custom signals - chrom gp
    if psr.name in chrom_psrs:
        chrom_gp=True
    else:
        chrom_gp=False
    # custom signals - J1012 RN
    if psr.name in hf_psrs:
        log10_A = parameter.Uniform(-20, -11)
        gamma = parameter.Uniform(0, 7)
        pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
        selection = selections.Selection(selections.no_selection)
        hfrn = gp_signals.FourierBasisGP(pl, components=Nfreqs_high,
                                         combine=True,
                                         selection=selection,
                                         name='hf_red_noise')
        extra_sigs = hfrn
    # custom signals - SWGP
    elif (np.abs(elats[psr.name]) < 20 and extra_swgp) or psr.name in sw_psrs:
        extra_sigs = sw_gp_signal()
    else:
        extra_sigs = None
    
    # below Tspan used for red noise and GWB
    noise = model_singlepsr_noise(psr, Tspan=Tspan, psr_model=psr_model,
                                  # timing
                                  tm_svd=tm_svd, tm_marg=tm_marg,
                                  # red noise
                                  log_A_min=-20, log_A_max=-11, components=Nfreqs_low,
                                  # white noise
                                  white_vary=white_vary, tnequad=True, inc_ecorr=inc_ecorr,
                                  efeq_groups=efeq_groups_by_PTA, ecorr_groups=ecorr_groups_by_PTA,
                                  log_equad_min=-10, log_equad_max=-4, gp_ecorr=gp_ecorr,
                                  # DM - uses psr TSpan
                                  dm_var=True, dm_type='gp',
                                  dmgp_kernel='diag', dm_psd='powerlaw',
                                  dm_Nfreqs=Nfreqs_high,
                                  # chrom_gp
                                  chrom_gp=chrom_gp,
                                  chrom_gp_kernel='diag', chrom_Nfreqs=Nfreqs_high,
                                  chrom_idx=chrom_idx,
                                  # solar wind
                                  dm_sw_deter=False,
                                  # dm dip
                                  dm_expdip=dm_expdip, dm_expdip_basename='exp',
                                  dm_expdip_idx=dmexp_idx,
                                  dm_expdip_tau_min=np.log10(5), dm_expdip_tau_max=np.log10(500),
                                  dm_expdip_tmin=54742, dm_expdip_tmax=54768,
                                  # extra sigs
                                  extra_sigs=extra_sigs)
    
    return noise



