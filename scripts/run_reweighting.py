#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:02:10 2022

@author: bjornlarsen
"""

import numpy as np
import argparse, json, os, glob, pickle, cloudpickle
from h5pulsar import Pulsar, FilePulsar
from mpi4py import MPI

from enterprise.signals import parameter
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import selections
from enterprise.signals.selections import Selection
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals
import enterprise.constants as const
from enterprise.signals import utils
from enterprise_extensions.blocks import common_red_noise_block
from IPTA_DR2_analysis.model_blocks import adv_noise_block

from dr3_noise.models import model_singlepsr_noise
from dr3_noise.model_utils import get_freqs
from dr3_noise.selections import CustomSelections
import la_forge.core as co
import time


# reweight a CRN run to HD

PARSER = argparse.ArgumentParser()

PARSER.add_argument('--coredir', type=str, default=None,
                    help='Path to directory where la_forge core should be loaded and results output')

PARSER.add_argument('--dataset', type=str, default=None,
                    help='Path to directory of premade PSR hdf5 files comprising the dataset')

PARSER.add_argument('--noisedict_path', type=str, default=None,
                    help='Path to noise dictionary json')

PARSER.add_argument("-N", dest="Niter", type=int, default=100_000,
                    help='Target number of reweighted samples')

PARSER.add_argument('--fixed_gamma', action='store_true',
                    dest='fixed_gamma', default=False,
                    help='Run with fixed gamma=13/3 (bool type)')

PARSER.add_argument('--freespec', action='store_true', default=False,
                    help='Flag to run a GWB free spectrum')

PARSER.add_argument('--gwb_Nfreqs', type=int, default=13,
                    help='Number of frequencies for GWB model',
                    nargs='?', const=None)

PARSER.add_argument('--human', type=str, default='Mr. Meeseeks',
                    help='Name of human who ran this (saved to runtime info)')

PARSER.add_argument('--restart', action='store_true', default=False,
                    help='Flag to restart the likelihood reweighting run')

args = PARSER.parse_args()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print(f'kicking off likelihood reweighting')
if args.dataset == None:
    raise NameError('Missing dataset')
elif args.coredir == None:
    raise NameError('Missing chains for reweighting')
else:
    print(f'Dataset: {args.dataset}')
    print(f'Corepath: {args.coredir}/core.h5')

# setup directory for chains
# head directory (dataset info)
outdir = f'{args.coredir}/HD_reweighting'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
# check if using job array to perform task in parallel
if "SLURM_ARRAY_TASK_ID" in list(os.environ.keys()):
    task_id = os.environ["SLURM_ARRAY_TASK_ID"]
    outdir += f'/{task_id}'
    if not os.path.isdir(outdir):
        try:
            os.mkdir(outdir)
        except Exception as e:
            print('Exception when setting up directories:', e)
            print('Directory probably exists already!')
else:
    task_id = 0
print(f'Output to: {outdir}')
if 'litec_filtered' in args.dataset:
    dataset = 'litec_filtered'
elif 'full' in args.dataset or 'litec' in args.dataset or 'unfiltered' in args.dataset or 'edr2' in args.dataset:
    dataset = 'full'
elif 'lite' in args.dataset:
    dataset = 'lite'
else:
    raise KeyError(f'Could not determine Tspan for dataset from outdir {outdir}')

# load noise params
with open(args.noisedict_path,'r') as f:
    noise_params = json.load(f)

# let's get our pulsars
psrnames = []
psrs = []
psr_paths = np.sort(glob.glob(f'{args.dataset}/J*'))
try:
    for psr_path in psr_paths:
        with open(psr_path,'rb') as f:
            psrs.append(FilePulsar(f))
except:
    for psr_path in psr_paths:
        with open(psr_path,'rb') as f:
            psrs.append(pickle.load(f))
print(f'Loaded {len(psrs)} pulsars from {args.dataset}')

# set up model
max_mjd = np.max([np.max(psr.toas) for psr in psrs])
min_mjd = np.min([np.min(psr.toas) for psr in psrs])
Tspan = max_mjd - min_mjd
print(f'Tspan for GWB: {Tspan/const.yr} yr')
print(f'Nfreqs for GWB: {args.gwb_Nfreqs}')
print(f'Fixed gamma = {args.fixed_gamma}')
# setup GWB priors/PSD
gamma_val = None
psd = 'powerlaw'
logmin = -18
logmax = -12
if args.fixed_gamma:
    gamma_val = 13/3
if args.freespec:
    psd = 'spectrum'
    logmin = -10
    logmax = -4
crn = common_red_noise_block(psd=psd, prior='log-uniform', Tspan=Tspan,
                             components=args.gwb_Nfreqs, gamma_val=gamma_val,
                             logmin=logmin, logmax=logmax, orf='hd', name='gwb')
models = []
for psr in psrs:
    noise = adv_noise_block(psr, full_pta_analysis=True, dataset=dataset, psr_model=True,
                            tm_marg=True, tm_svd=True)
    signals = crn + noise
    models.append(signals(psr))
pta = signal_base.PTA(models)
pta.set_default_params(noise_params)

with open(f'{outdir}/pta_summary.txt','w') as f:
    f.write(pta.summary())

# time initial likelihood
params = {p.name:p.sample() for p in pta.params}
t0 = time.time()
lnlike = pta.get_lnlikelihood(params)
print(f'lnlike = {lnlike}')
print(f'time = {time.time() - t0} s')

# ----------------------------
# REWEIGHTING SETUP
# ----------------------------

# load samples
core_crn = co.Core(corepath=f'{args.coredir}/core.h5', burn=0)

# determine thinning
Ns = len(core_crn('lnpost'))
thin = Ns//args.Niter
print(f'Aiming for {args.Niter} reweighted samples')
print(f'Starting from {Ns} reweighted samples')
print(f'Thinning by {thin}')

# resuming run
if os.path.isfile(f'{outdir}/ln_weights.txt') and not args.restart:
    # load in idxs + values so far
    chain_idxs = np.loadtxt(f'{outdir}/chain_idxs.txt') # which values of the chain to reference
    ln_weights = np.loadtxt(f'{outdir}/ln_weights.txt') # computed ln_weights for each sample in the chain
    new_array_idxs = np.arange(len(ln_weights))#reference to keep track of where in the chain_idxs and ln_weights array we are at
    # print stuff
    if not np.any(np.isnan(ln_weights)):
        print('You have already finished processing this run! (no more nans in output ln_weights file)')
        raise ValueError('You have already finished processing this run! (no more nans in output ln_weights file)')
    print(f'resuming run ({len(chain_idxs[~np.isnan(ln_weights)])}/{len(chain_idxs)} processed)')
    
# starting new run
else:
    # setup arrays of chain idxs + ln weight values + ln weight idxs
    # uses different idxs for different jobs in array
    chain_idxs = np.arange(int(task_id),Ns,thin,dtype=int) # which values of the chain to reference
    ln_weights = np.zeros(len(chain_idxs))*np.nan # computed ln_weights for each sample in the chain
    new_array_idxs = np.arange(len(chain_idxs))#reference to keep track of where in the chain_idxs and ln_weights array we are at
    # save idxs to process
    np.savetxt(f'{outdir}/chain_idxs.txt', chain_idxs)
    # print stuff
    print(f'Starting new run (0/{len(chain_idxs)} processed)')
    if os.path.isfile(f'{outdir}/ln_weights.txt'):
        print("Overwriting the existing file!!!")
    
# setup mask to filter out processed values (only process samples which are still nan in the ln_weights array)
mask = np.isnan(ln_weights)

start_time = time.time()
last_time = time.time()

ct = 0
crn_params = [p for p in core_crn.params if 'crn' in p]
for chain_idx, new_array_idx in zip(chain_idxs, new_array_idxs):
    # get sample
    chain_idx = int(chain_idx)
    params = {p:core_crn(p)[chain_idx] for p in core_crn.params}
    # duplicate GWB params for HD model
    if args.freespec:
        params['gwb_log10_rho'] = np.array([params[p] for p in crn_params])
    else:
        for p in crn_params:
            params[p.replace('crn','gwb')] = params[p]
    # get ln_weight for given sample
    ln_weights[new_array_idx] = pta.get_lnlikelihood(params) - core_crn('lnlike')[chain_idx]
    
    if time.time() - last_time > 60:
        # compute stats
        Nvals = np.count_nonzero(~np.isnan(ln_weights))
        if Nvals > 0:
            mean_wgt = np.nanmean(np.exp(ln_weights))
            sigma_wgt = np.nanstd(np.exp(ln_weights))
            neff = Nvals/(1 + (sigma_wgt/mean_wgt)**2)
            eff = neff/Nvals
            sigma_B = sigma_wgt/np.sqrt(neff)

            # print progress
            postfix_str = (f'{Nvals}/{len(chain_idxs)} | B = {np.round(mean_wgt,decimals=3)} +- {np.round(sigma_B,decimals=3)} |'
                           f' neff = {np.round(neff,decimals=3)} | efficiency = {np.round(eff,decimals=3)}')
            print(postfix_str)

            # update array
            np.savetxt(f'{outdir}/ln_weights.txt', ln_weights)
        else:
            print('No non-nan values in array yet')
        last_time = time.time()

# Save final results
print('FINISHED!!')
Nvals = np.count_nonzero(~np.isnan(ln_weights))
mean_wgt = np.nanmean(np.exp(ln_weights))
sigma_wgt = np.nanstd(np.exp(ln_weights))
neff = Nvals/(1 + (sigma_wgt/mean_wgt)**2)
eff = neff/Nvals
sigma_B = sigma_wgt/np.sqrt(neff)
postfix_str = (f'{Nvals}/{len(chain_idxs)} | B = {np.round(mean_wgt,decimals=3)} +- {np.round(sigma_B,decimals=3)} | '
               f'neff = {np.round(neff,decimals=3)} | efficiency = {np.round(eff,decimals=3)}')
print(postfix_str)
np.savetxt(f'{outdir}/ln_weights.txt', ln_weights)
print(f'Total time: {time.time() - start_time}')
