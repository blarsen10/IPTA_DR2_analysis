#!/usr/bin/env python
# coding: utf-8

import numpy as np
import argparse, sys, os, glob, json, pickle, copy
import cloudpickle
import logging

from nautilus import Prior
from ipta_gwb_analysis.models import pta_model
from h5pulsar import FilePulsar
import numpy as np
from scipy.stats import multivariate_normal, uniform
from nautilus import Sampler
import time

from dr3_noise.models import model_singlepsr_noise
#from dr3_noise.model_utils import get_freqs, get_flag_groups_by_PTA
#from dr3_noise.pre_processing_utils import make_pulsar_object
from enterprise_extensions.blocks import common_red_noise_block, red_noise_block
from enterprise_extensions.sampler import setup_sampler, group_from_params, get_parameter_groups
from enterprise_extensions import chromatic as chrom
import enterprise.constants as const
from targeted_cws_ng15.Dists_Parameters import get_parameter_groups_CAW_target
from enterprise.signals import signal_base, parameter, utils, selections, gp_signals, deterministic_signals
from IPTA_DR2_analysis.sw import elats, sw_gp_signal
from IPTA_DR2_analysis.model_blocks import adv_noise_block, full_Tspan, lite_Tspan

from mpi4py.futures import MPIPoolExecutor

parser = argparse.ArgumentParser()

parser.add_argument('--outdir', type=str, default=None,
                    help='Path to store sampler/output from the nautilus run')

parser.add_argument('--psr', type=str, default=None,
                    help='Path to premade enterprise PSR hdf5 file',
                    nargs='?', const=None)

parser.add_argument('--noisepath', type=str, default=None,
                    help='Path to noise dictionary json')

parser.add_argument('--gwb_Nfreqs', type=int, default=13,
                    help='Number of frequencies for GWB model',
                    nargs='?', const=None)

parser.add_argument("--Nlive", default=1000, type=int, 
                    help='Number of live points')

parser.add_argument("--Nnet", default=4, type=int, 
                    help='Number of networks')

parser.add_argument("--Ncpus", default=1, type=int, 
                    help='Number of cores for multiprocessing')

parser.add_argument('--MPI', dest='MPI', action='store_true',
                    help="Whether to use MPI")

parser.add_argument('--restart', dest='restart', action='store_true',
                    help="Whether to restart the run")

args = parser.parse_args()



print('setting up nested sampler')
if args.psr == None:
    raise NameError('Missing pulsar file')
elif args.outdir == None:
    raise NameError('Missing output directory')
else:
    # identify dataset
    if 'full' in args.psr or 'litec' in args.psr or 'unfiltered' in args.psr:
        dataset = 'full'
        Tspan = full_Tspan
    elif 'lite' in args.psr:
        dataset = 'lite'
        Tspan = lite_Tspan
    else:
        raise NameError(f'Could not determine Tspan for dataset from psrfile {args.psr}')
    psrname = args.psr.split('/')[-1].split('.')[0]
    print(f'Pulsar: {psrname}')
    print(f'Output to: {args.outdir}')
    print(f'Dataset: {dataset}')

    
# get PSR pickle
if args.psr is None:
    raise NameError('Missing data (psr set to None)')
else:
    print(f'Loading Pulsar object from file: {args.psr}')
    try:
        with open(args.psr,'rb') as f:
            psr = FilePulsar(f)
    except:
        with open(args.psr,'rb') as f:
            psr = pickle.load(f)
    
    
if not os.path.isdir(args.outdir):
    try:
        os.mkdir(args.outdir) # sometimes this breaks for no apparent reason
    except:
        pass
if args.MPI:
    filepath=f'{args.outdir}/{psrname}_sampler_MPI.hdf5'
else:
    filepath=f'{args.outdir}/{psrname}_sampler.hdf5'
print(f'output to {filepath}')

# check if need to restart
if args.restart:
    print('Restarting run!')
    resume = False
else:
    resume = True
        
# check pta flags
if np.any(np.unique(psr.flags['pta']) == ''):
    raise ValueError('ERROR: Blank -pta flags found in the parfile. '
                     'This will break the noise modeling!')

# print fitpars
print(f'Free TM params from PSR pickle: {psr.fitpars}')

# load noise params
with open(args.noisepath,'r') as f:
    noise_params = json.load(f)

# set up model
print(f'Tspan for FL GWB: {Tspan/const.yr} yr')
print(f'Nfreqs for FL GWB: {args.gwb_Nfreqs}')
crn = common_red_noise_block(psd='powerlaw', prior='log-uniform', Tspan=Tspan,
                             components=args.gwb_Nfreqs, gamma_val=13/3,
                             logmin=-18, logmax=-12, orf=None, name='crn')
noise = adv_noise_block(psr, full_pta_analysis=True, dataset=dataset, psr_model=True,
                        tm_marg=True, tm_svd=True)
signals = crn + noise
pta = signal_base.PTA([signals(psr)])
pta.set_default_params(noise_params)

# save PTA
#with open(args.outdir+f'/{psrname}_pta.pkl','wb') as f:
#    cloudpickle.dump(pta, f)

# saving model params in chain directory
with open(args.outdir+f'/{psrname}_model_params.json' , 'w') as fout:
    json.dump(pta.param_names, fout, sort_keys=True,
              indent=4, separators=(',', ': '))
# saving runtime info
with open(args.outdir+f"/{psrname}_runtime_info.txt", "w") as fout:
    fout.write(pta.summary())
    
def get_prior_distr(param):
    pline = str(param)
    pmin = float(pline[pline.index('pmin')+5:pline.index(', pmax')])
    pmax = float(pline[pline.index('pmax')+5:-1])
    prior_dist = uniform(loc=pmin, scale=pmax-pmin)
    return prior_dist

# set up prior for nautilus
prior = Prior()
for i in range(len(pta.params)):
    prior.add_parameter(pta.param_names[i], dist=get_prior_distr(pta.params[i]))
print(f'dim = {prior.dimensionality()}')

# set up sampler
if args.MPI:
    print('using MPI')
    sampler = Sampler(prior, pta.get_lnlikelihood, n_live=args.Nlive, filepath=filepath,
                      n_networks=args.Nnet, pass_dict=False, pool=MPIPoolExecutor(), resume=resume)
else:
    sampler = Sampler(prior, pta.get_lnlikelihood, n_live=args.Nlive, filepath=filepath,
                      n_networks=args.Nnet, pass_dict=False, pool=args.Ncpus, resume=resume)

# run
print('Starting sampler')
t0 = time.time()
sampler.run(verbose=True, discard_exploration=True)
print(f'Ran in {time.time()-t0} s')

# indicate sampler has converged
with open(args.outdir+f"/{psrname}_converged.txt", "w") as fout:
    fout.write('')

