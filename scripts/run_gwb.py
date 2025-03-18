#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:02:10 2022

Std noise run script - for use w/ parallel tempering
"""

import numpy as np
import argparse, json, os, glob
from h5pulsar import FilePulsar

from dr3_noise.models import model_singlepsr_noise
#from dr3_noise.model_utils import get_freqs, get_flag_groups_by_PTA
#from dr3_noise.pre_processing_utils import make_pulsar_object
from enterprise_extensions.blocks import common_red_noise_block
from enterprise_extensions.sampler import setup_sampler, group_from_params
import enterprise.constants as const
from targeted_cws_ng15.Dists_Parameters import get_parameter_groups_CAW_target
from enterprise.signals import signal_base, deterministic_signals

parser = argparse.ArgumentParser()

parser.add_argument('--outdir', type=str, default=None,
                    help='Path to store output from the PTMCMCSampler run')

parser.add_argument('--dataset', type=str, default=None,
                    help='Path to directory of premade PSR hdf5 files comprising the dataset')

parser.add_argument('--noisedict_path', type=str, default=None,
                    help='Path to noise dictionary json')

parser.add_argument('--emp_dist_path', type=str, default=None,
                    help='Path to empirical distribution file',
                    nargs='?', const=None)

parser.add_argument("--Niter", default=2_000_000, type=int, 
                    help='Number of PTMCMC iterations')

parser.add_argument('--fixed_gamma', action='store_true',
                    dest='fixed_gamma', default=False,
                    help='Run with fixed gamma=13/3 (bool type)')

parser.add_argument('--bayesephem', action='store_true',
                    dest='bayesephem', default=False,
                    help='Run with BayesEphem')

parser.add_argument('--human', type=str, default='Mr. Meeseeks',
                    help='Name of human who ran this (saved to runtime info)')

args = parser.parse_args()


print('kicking off parallel tempering')
if args.outdir == None:
    raise NameError('Missing output directory')
elif args.dataset == None:
    raise NameError('Missing dataset')
elif args.noisedict_path == None:
    raise NameError('Missing noise dictionary')
else:
    print(f'Dataset from: {args.dataset}')
    print(f'Output to: {args.outdir}')

# setup directory for chains
if not os.path.isdir(args.outdir):
    try:
        os.mkdir(args.outdir) # sometimes this breaks for no apparent reason
    except:
        pass
# support for chain parallelization using job arrays
if "SLURM_ARRAY_TASK_ID" in list(os.environ.keys()):
    task_id = os.environ["SLURM_ARRAY_TASK_ID"]
    args.outdir += f'/{task_id}'
    if not os.path.isdir(args.outdir):
        try:
            os.mkdir(args.outdir) 
        except:
            pass

# check number of samples
try:
    with open(f'{args.outdir}/chain_1.0.txt', 'r') as f:
        n_samples = len(f.readlines())
except:
    n_samples = 0
if n_samples < 100:
    print(f'only {n_samples} samples so far, setting resume = False!')
    resume = False
else:
    print(f'{n_samples} samples so far')
    resume = True

# load noise params
with open(args.noisedict_path,'r') as f:
    noise_params = json.load(f)
    
# get PSRs
print(f'Loading Pulsars from file')
psrs = []
psr_paths = np.sort(glob.glob(f'{args.dataset}/J*'))
for psr_path in psr_paths:
    with open(psr_path,'rb') as f:
        psrs.append(FilePulsar(f))
print(f'Loaded {len(psrs)} pulsars')

# set up model
max_mjd = np.max([np.max(psr.toas) for psr in psrs])
min_mjd = np.min([np.min(psr.toas) for psr in psrs])
Tspan = max_mjd - min_mjd
if args.fixed_gamma:
    gamma_val = 13/3
else:
    gamma_val = None
crn = common_red_noise_block(psd='powerlaw', prior='log-uniform', Tspan=Tspan,
                             components=13, gamma_val=gamma_val,
                             logmin=-18, logmax=-12,
                             orf=None, name='crn')
be = deterministic_signals.PhysicalEphemerisSignal(use_epoch_toas=True, model='orbel')

# loop through pulsars, add noise
models = []
for psr in psrs:
    psr_Tspan = np.max(psr.toas) - np.min(psr.toas)
    if 'NANOGrav' in psr.flags['pta']:
        inc_ecorr = True
        ecorr_groups_by_PTA = {'NANOGrav':'group'}
    else:
        inc_ecorr = False
        ecorr_groups_by_PTA = None
    efeq_groups_by_PTA = {}
    for pta in np.unique(psr.flags['pta']):
        efeq_groups_by_PTA[pta] = 'group'
    if psr.name == 'J1713+0747':
        dm_expdip=True
    else:
        dm_expdip=False
    print(f'{psr.name} DM Tspan: {psr_Tspan}')
    noise = model_singlepsr_noise(psr, Tspan=Tspan, psr_model=True,
                                  # timing
                                  tm_svd=True, tm_marg=True,
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
                                  dm_expdip=dm_expdip, dm_expdip_basename='exp',
                                  dm_expdip_tau_min=np.log10(5), dm_expdip_tau_max=np.log10(500), 
                                  # red noise
                                  log_A_min=-20, log_A_max=-11)
    if args.bayesephem:
        signals = crn + noise + be
    else:
        signals = crn + noise
    models.append(signals(psr))
pta = signal_base.PTA(models)
pta.set_default_params(noise_params)

# set groups for adaptive metropolis
gr = get_parameter_groups_CAW_target(pta)
with open(f'{args.outdir}/groups.txt', 'w') as fi:
    for group in gr:
        line = np.array(pta.param_names)[np.array(group)]
        fi.write("[" + " ".join(line) + "]\n")

# initial sample
x0 = np.hstack([p.sample() for p in pta.params])
print("lnlike: ", pta.get_lnlikelihood(x0))



################ SAMPLER SETUP #####################
if args.emp_dist_path is not None:
    print(f"Using empirical distributions file path: {args.emp_dist_path}")
else:
    print("No empirical distributions file given")
sampler = setup_sampler(pta, outdir=args.outdir, groups=gr, resume=resume,
                        human=args.human, empirical_distr=args.emp_dist_path)

'''def draw_from_sw_prior(self, x, iter, beta):

    q = x.copy()
    lqxy = 0

    signal_name = 'sw_r2'

    # draw parameter from signal model
    param = np.random.choice(self.snames[signal_name])
    if param.size:
        idx2 = np.random.randint(0, param.size)
        q[self.pmap[str(param)]][idx2] = param.sample()[idx2]

    # scalar parameter
    else:
        q[self.pmap[str(param)]] = param.sample()

    # forward-backward jump probability
    lqxy = (param.get_logpdf(x[self.pmap[str(param)]]) -
            param.get_logpdf(q[self.pmap[str(param)]]))

    return q, float(lqxy)

sampler.jp.draw_from_sw_prior = draw_from_sw_prior
sampler.addProposalToCycle(sampler.jp.draw_from_sw_prior, 10)'''

print("Beginning to sample ...")
sampler.sample(x0, args.Niter, burn=3_000, thin=10,
               SCAMweight=200, AMweight=100, DEweight=200,
               writeHotChains=True)


