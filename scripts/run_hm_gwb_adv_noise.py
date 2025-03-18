#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:02:10 2022

Std noise run script - for use w/ parallel tempering
"""

import numpy as np
import argparse, json, os, glob, pickle
from h5pulsar import FilePulsar

from dr3_noise.models import model_singlepsr_noise
#from dr3_noise.model_utils import get_freqs, get_flag_groups_by_PTA
#from dr3_noise.pre_processing_utils import make_pulsar_object
from enterprise_extensions.blocks import common_red_noise_block, red_noise_block
from enterprise_extensions.sampler import setup_sampler, group_from_params, get_parameter_groups
from enterprise_extensions import chromatic as chrom
from enterprise_extensions import hypermodel
import enterprise.constants as const
from targeted_cws_ng15.Dists_Parameters import get_parameter_groups_CAW_target
from enterprise.signals import signal_base, parameter, utils, selections, gp_signals, deterministic_signals
from IPTA_DR2_analysis.sw import elats, sw_gp_signal
from IPTA_DR2_analysis.model_blocks import adv_noise_block

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

parser.add_argument('--freespec', action='store_true', default=False,
                    help='Flag to run a GWB free spectrum')

parser.add_argument('--gwb_Nfreqs', type=int, default=13,
                    help='Number of frequencies for GWB model',
                    nargs='?', const=None)

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
if 'full' in args.dataset or 'litec' in args.dataset or 'unfiltered' in args.dataset or 'edr2' in args.dataset:
    dataset = 'full'
elif 'lite' in args.dataset:
    dataset = 'lite'
else:
    raise KeyError(f'Could not determine Tspan for dataset from outdir {outdir}')

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
try:
    for psr_path in psr_paths:
        with open(psr_path,'rb') as f:
            psrs.append(FilePulsar(f))
except:
    for psr_path in psr_paths:
        with open(psr_path,'rb') as f:
            psrs.append(pickle.load(f))
print(f'Loaded {len(psrs)} pulsars')

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
crn = common_red_noise_block(psd='powerlaw', prior='log-uniform', Tspan=Tspan,
                             components=args.gwb_Nfreqs, gamma_val=gamma_val,
                             logmin=logmin, logmax=logmax, orf=None, name='crn')
gwb = common_red_noise_block(psd='powerlaw', prior='log-uniform', Tspan=Tspan,
                             components=args.gwb_Nfreqs, gamma_val=gamma_val,
                             logmin=logmin, logmax=logmax, orf='hd', name='gwb')
be = deterministic_signals.PhysicalEphemerisSignal(use_epoch_toas=True, model='orbel')

# loop through pulsars, add noise
crn_models = []
gwb_models = []
for psr in psrs:
    noise = adv_noise_block(psr, full_pta_analysis=True, dataset=dataset, psr_model=True,
                            tm_marg=True, tm_svd=True)
    # apply crn model
    signals = crn + noise
    crn_models.append(signals(psr))
    # apply gwb model
    signals = gwb + noise
    gwb_models.append(signals(psr))
pta_crn = signal_base.PTA(crn_models)
pta_gwb = signal_base.PTA(gwb_models)
pta_crn.set_default_params(noise_params)
pta_gwb.set_default_params(noise_params)

model_dictionary = {0:pta_crn, 1:pta_gwb}

# make a dictionary in the file for model labels and model params
model_params = {}
for ky, pta in model_dictionary.items():
    model_params.update({str(ky) : pta.param_names})
with open(args.outdir+'/model_params.json' , 'w') as fout:
    json.dump(model_params, fout, sort_keys=True,
          indent=4, separators=(',', ': '))
with open(args.outdir+'/model_labels.json' , 'w') as fout:
    json.dump([[0,'CRN'],[1,'HD']], fout, sort_keys=True,
          indent=4, separators=(',', ': '))  

hyper_model = hypermodel.HyperModel(model_dictionary)
    
# set groups for adaptive metropolis
# groups - start with the default groups
pnames = hyper_model.param_names
gr = hyper_model.get_parameter_groups()
# custom groups don't work yet because there are indices for multiple models in the hypermodel
'''# add custom groups
gr.extend(get_parameter_groups_CAW_target(pta_crn))
gr.extend(get_parameter_groups_CAW_target(pta_gwb, gw_only=True))
gr.extend(['nmodel'])
gr.extend(['nmodel'])
gr.extend(['nmodel'])
with open(f'{args.outdir}/groups.txt', 'w') as fi:
    for group in gr:
        line = np.array(pta.param_names)[np.array(group)]
        fi.write("[" + " ".join(line) + "]\n")'''

# initial sample
x0 = hyper_model.initial_sample()
print("lnlike: ", pta.get_lnlikelihood(x0))



################ SAMPLER SETUP #####################
if args.emp_dist_path is not None:
    print(f"Using empirical distributions file path: {args.emp_dist_path}")
else:
    print("No empirical distributions file given")
sampler = hyper_model.setup_sampler(outdir=args.outdir, resume=True,
                                    empirical_distr=args.emp_dist_path, 
                                    groups=gr)

print("Beginning to sample ...")
sampler.sample(x0, args.Niter, burn=3_000, thin=10,
               SCAMweight=150, AMweight=100, DEweight=200,
               writeHotChains=True)


