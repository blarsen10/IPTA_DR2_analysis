#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:02:10 2022

Std noise run script - for use w/ parallel tempering
"""

import numpy as np
import argparse, json, os, pickle
from h5pulsar import FilePulsar

from dr3_noise.models import model_singlepsr_noise
#from dr3_noise.model_utils import get_freqs, get_flag_groups_by_PTA
#from dr3_noise.pre_processing_utils import make_pulsar_object
from enterprise_extensions.sampler import setup_sampler, group_from_params, get_parameter_groups
import enterprise.constants as const
from enterprise.signals import signal_base, parameter, utils, selections, gp_signals
from enterprise_extensions.blocks import red_noise_block
from enterprise_extensions import chromatic as chrom
from IPTA_DR2_analysis.sw import sw_gp_signal
from IPTA_DR2_analysis.model_blocks import adv_noise_block

parser = argparse.ArgumentParser()

parser.add_argument('--psrname', type=str, default=None,
                    help='Name of the pulsar')

parser.add_argument('--outdir', type=str, default=None,
                    help='Path to store output from the PTMCMCSampler run')

# TODO: Add support for input from hdf5 file instead of pickle
parser.add_argument('--psr', type=str, default=None,
                    help='Path to premade enterprise PSR hdf5 file',
                    nargs='?', const=None)

parser.add_argument('--dataset', type=str, default=None,
                    help='Used to set frequency basis. Should be "lite" or "full"',
                    nargs='?', const=None)

parser.add_argument('--emp_dist_path', type=str, default=None,
                    help='Path to empirical distribution file',
                    nargs='?', const=None)

parser.add_argument("--Niter", default=2_000_000, type=int, 
                    help='Number of PTMCMC iterations')

parser.add_argument('--human', type=str, default='Mr. Meeseeks',
                    help='Name of human who ran this (saved to runtime info)')

args = parser.parse_args()

print('kicking off parallel tempering')
if args.psrname == None:
    raise NameError('Missing pulsar name')
elif args.dataset == None:
    raise NameError('Missing dataset name')
elif args.outdir == None:
    raise NameError('Missing output directory')
else:
    print(f'Pulsar: {args.psrname}')
    print(f'Dataset: {args.dataset}')
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
print(f'Output to: {args.outdir}')

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
    
# make pta in one line!!
pta = adv_noise_block(psr, full_pta_analysis=False, psr_model=False,
                      tm_marg=False, tm_svd=True)

# saving model params in chain directory
with open(args.outdir+'/model_params.json' , 'w') as fout:
    json.dump(pta.param_names, fout, sort_keys=True,
              indent=4, separators=(',', ': '))

# groups - start with the default groups
pnames = pta.param_names
groups = get_parameter_groups(pta)
# add per-backend white noise
backends = np.unique([p[p.index('_')+1:p.index('efac')-1] for p in pnames
                      if 'efac' in p])
print(f'Backends w/ unique noise params: {backends}')
for be in backends:
    groups.append(group_from_params(pta,[f'{args.psrname}_{be}_']))
# group red noise + chromatic parameters
exclude = ['linear_timing_model','measurement_noise','tnequad',
           'ecorr_sherman-morrison','ecorr_fast-sherman-morrison']
red_signals = [p[p.index('_')+1:] for p in list(pta.signals.keys())
               if not p[p.index('_')+1:] in exclude]
rn_ct = 0
for rs in red_signals:
    if len(group_from_params(pta,[rs])) > 0:
        rn_ct += 1
        groups.append(group_from_params(pta,[rs]))
if rn_ct > 1:
    groups.append(group_from_params(pta,red_signals))
# add cross chromatic groups
if np.any(['sw_gp' in p for p in pnames]):
    # cross SW and chrom groups
    dmgp_sw = [idx for idx, nm in enumerate(pnames)
               if any([flag in nm for flag in ['dm_gp','sw_gp']])]
    groups.append(dmgp_sw)
    if np.any(['chrom' in param for param in pnames]):
        chromgp_sw = [idx for idx, nm in enumerate(pnames)
                      if any([flag in nm for flag in ['chrom','sw_gp']])]
        dmgp_chromgp_sw = [idx for idx, nm in enumerate(pnames)
                           if any([flag in nm for flag in ['dm_gp','chrom','sw_gp']])]
        groups.append(chromgp_sw)
        groups.append(dmgp_chromgp_sw)
if np.any(['chrom' in param for param in pnames]):
    # cross dmgp and chromgp group
    dmgp_chromgp = [idx for idx, nm in enumerate(pnames)
                    if any([flag in nm for flag in ['dm_gp','chrom']])]
    groups.append(dmgp_chromgp)
# group with all params
groups.append([i for i in range(len(pta.params))])
# save list of params corresponding to groups
with open(f'{args.outdir}/groups.txt', 'w') as fi:
    for group in groups:
        line = np.array(pnames)[np.array(group)]
        fi.write("[" + " ".join(line) + "]\n")
        
        
        
        
# initial sample
x0 = np.hstack([p.sample() for p in pta.params])
print("lnlike: ", pta.get_lnlikelihood(x0))



################ SAMPLER SETUP #####################
if args.emp_dist_path is not None:
    print(f"Using empirical distributions file path: {args.emp_dist_path}")
else:
    print("No empirical distributions file given")
sampler = setup_sampler(pta, outdir=args.outdir, groups=groups, resume=resume,
                        human=args.human, empirical_distr=args.emp_dist_path)

def draw_from_sw_gp_prior(self, x, iter, beta):

    q = x.copy()
    lqxy = 0

    signal_name = 'sw_gp'

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

def draw_from_exp_prior(self, x, iter, beta):

    q = x.copy()
    lqxy = 0

    signal_name = 'exp'

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

def draw_from_hf_prior(self, x, iter, beta):

    q = x.copy()
    lqxy = 0

    signal_name = 'hf_red_noise'

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

# this does not ork for some reason
# something about missing 1 required positional argument: 'beta'
'''if np.any(['sw_gp' in p for p in pnames]):
    print('Adding SW GP prior draw...')
    sampler.jp.draw_from_sw_gp_prior = draw_from_sw_gp_prior
    sampler.addProposalToCycle(sampler.jp.draw_from_sw_gp_prior, 25)
if np.any(['exp' in p for p in pnames]):
    print('Adding exp dip prior draw...')
    sampler.jp.draw_from_exp_prior = draw_from_exp_prior
    sampler.addProposalToCycle(sampler.jp.draw_from_exp_prior, 25)
if np.any(['hf_red_noise' in p for p in pnames]):
    print('Adding HF red noise prior draw...')
    sampler.jp.draw_from_hf_prior = draw_from_hf_prior
    sampler.addProposalToCycle(sampler.jp.draw_from_hf_prior, 25)'''

print("Beginning to sample ...")
sampler.sample(x0, args.Niter, burn=3_000, thin=10,
               SCAMweight=200, AMweight=100, DEweight=200,
               writeHotChains=True)


