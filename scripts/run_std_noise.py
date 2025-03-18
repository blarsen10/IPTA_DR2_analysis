#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:02:10 2022

Std noise run script - for use w/ parallel tempering
"""

import numpy as np
import argparse, json, os
from h5pulsar import FilePulsar

from dr3_noise.models import model_singlepsr_noise
#from dr3_noise.model_utils import get_freqs, get_flag_groups_by_PTA
#from dr3_noise.pre_processing_utils import make_pulsar_object
from enterprise_extensions.sampler import setup_sampler, group_from_params
import enterprise.constants as const

parser = argparse.ArgumentParser()

parser.add_argument('--psrname', type=str, default=None,
                    help='Name of the pulsar')

parser.add_argument('--outdir', type=str, default=None,
                    help='Path to store output from the PTMCMCSampler run')

# TODO: Add support for input from hdf5 file instead of pickle
parser.add_argument('--psr', type=str, default=None,
                    help='Path to premade enterprise PSR hdf5 file',
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
elif args.outdir == None:
    raise NameError('Missing output directory')
else:
    print(f'Pulsar: {args.psrname}')
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
    with open(args.psr,'rb') as f:
        psr = FilePulsar(f)

# check pta flags
if np.any(np.unique(psr.flags['pta']) == ''):
    raise ValueError('ERROR: Blank -pta flags found in the parfile. '
                     'This will break the noise modeling!')

# print fitpars
print(f'Free TM params from PSR pickle: {psr.fitpars}')

# set up noise flags + selections
if 'NANOGrav' in psr.flags['pta']:
    inc_ecorr = True
    ecorr_groups_by_PTA = {'NANOGrav':'group'}
else:
    inc_ecorr = False
    ecorr_groups_by_PTA = None
efeq_groups_by_PTA = {}
for pta in np.unique(psr.flags['pta']):
    efeq_groups_by_PTA[pta] = 'group'
    
# set up J1713 dip
if psr.name == 'J1713+0747':
    dm_expdip=True
else:
    dm_expdip=False

# set up DM Nfreqs
pta = model_singlepsr_noise(psr, #Tspan=952746385.6296968,
                            # timing
                            tm_svd=True,
                            # white noise
                            tnequad=True, inc_ecorr=inc_ecorr,
                            efeq_groups=efeq_groups_by_PTA, ecorr_groups=ecorr_groups_by_PTA,
                            log_equad_min=-10, log_equad_max=-4,
                            # DM
                            dm_var=True, dm_type='gp',
                            dmgp_kernel='diag', dm_psd='powerlaw',
                            dm_Nfreqs=30,
                            # solar wind
                            dm_sw_deter=False,
                            # dm dip
                            dm_expdip=dm_expdip, dm_expdip_basename='exp',
                            dm_expdip_tau_min=np.log10(5), dm_expdip_tau_max=np.log10(500), 
                            # red noise
                            log_A_min=-20, log_A_max=-11)

# saving model params in chain directory
with open(args.outdir+'/model_params.json' , 'w') as fout:
    json.dump(pta.param_names, fout, sort_keys=True,
              indent=4, separators=(',', ': '))

# set groups for adaptive metropolis
gr = []
# group white noise parameters
pnames = pta.param_names
backends = np.unique([p[p.index('_')+1:p.index('efac')-1] for p in pnames
                      if 'efac' in p])
print(f'Backends w/ unique noise params: {backends}')
for be in backends:
    gr.append(group_from_params(pta,[f'{args.psrname}_{be}_']))
# group correlated noise parameters
# group red noise parameters
exclude = ['linear_timing_model','measurement_noise','tnequad',
           'ecorr_sherman-morrison','ecorr_fast-sherman-morrison']
red_signals = [p[p.index('_')+1:] for p in list(pta.signals.keys())
               if not p[p.index('_')+1:] in exclude]
rn_ct = 0
for rs in red_signals:
    if len(group_from_params(pta,[rs])) > 0:
        rn_ct += 1
        gr.append(group_from_params(pta,[rs]))
if rn_ct > 1:
    gr.append(group_from_params(pta,red_signals))
# group with all params
gr.append([i for i in range(len(pta.params))])
# save list of params corresponding to groups
with open(f'{args.outdir}/groups.txt', 'w') as fi:
    for group in gr:
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


