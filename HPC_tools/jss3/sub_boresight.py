import numpy as np
import subprocess
import sys
import os
import uuid
import time

#------ User dependent parameters ------
jss_account = ""
user_email  = ""
coderoot    = ""
conda_base  = ""
bizcode     = ""

base_path = ""
julia_path = ""
coderoot  = ""
imo_path = ""

#--------- HPC setting -------------
job_name    = "ScnFldCal"
vnode       = 1
vnode_core  = 2   # Maxmum of the RURI: 36 cores
vnode_mem   = 256   # Unit: GiB
#mode        = "debug"
mode        = "default"

elapse      = "24:00:00" # When you use the `debug` mode you should requesgt <= 1800 == "00:30:00"
if mode == "debug":
    elapse  = "00:10:00"

#------------------------------------

logdir    = os.path.join(coderoot, 'log')
ancillary = os.path.join(coderoot, 'ancillary')
if not os.path.exists(base_path):
    os.makedirs(base_path)
if not os.path.exists(logdir):
    os.makedirs(logdir, exist_ok=True)
if not os.path.exists(ancillary):
    os.makedirs(ancillary, exist_ok=True)
toml_filename   = str(uuid.uuid4())


#---------- Simulation setting ---------
# general
imo_version = 'v2'
telescope   = "boresight"
channel     = "boresight"
det         = "boresight"

# simulation
name            = 'Scan field calculation'
nside           = 256
n_year          = 3
one_year        = 3600*24*365
duration_s      = n_year * one_year
division        = 320
sampling_rate   = 19.
hwp_rpm         = 61.0
coord           = "G"
gamma           = 0

jobid = 0
jobname = "scnfld_"+str(nside)+"_"+channel+"_"+"id_"+str(jobid)
step = 0.5
spin_max = 20
spin_n = np.arange(0, spin_max+step, step) #[0,1,2,3,4,5,6,7,8,9,10]
spin_m = np.arange(-spin_max, spin_max+step, step) #[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]

tomlfile_path = os.path.join(ancillary, toml_filename+'.toml')
tomlfile_data = f"""
[general]
imo_path = '{imo_path}'
imo_version = '{imo_version}'
telescope = '{telescope}'
channel = '{channel}'
det_name = '{det}'
sampling_rate = {sampling_rate}
nside = {nside}

[simulation]
name = '{name}'
base_path = '{base_path}'
gamma = {gamma}
spin_n = [{', '.join(map(str, spin_n))}]
spin_m = [{', '.join(map(str, spin_m))}]
duration_s = '{duration_s}'
division = '{division}'
hwp_rpm = '{hwp_rpm}'
coord = '{coord}'
"""
with open(tomlfile_path, 'w') as f:
    f.write(tomlfile_data)

jobscript_path = os.path.join(ancillary, toml_filename+".jx")
jobscript_data = f"""#!/bin/zsh
#JX --bizcode {bizcode}
#JX -L rscunit=RURI
#JX -L rscgrp={mode}
#JX -L elapse={elapse}
#JX -L vnode={vnode}
#JX -L vnode-core={vnode_core}
#JX -L vnode-mem={vnode_mem}Gi

#JX -o {coderoot}/log/%n_%j.out
#JX -e {coderoot}/log/%n_%j.err
#JX --spath {coderoot}/log/%n_%j.stats
#JX -N {job_name}
#JX -m e
#JX --mail-list {user_email}
#JX -S
#export OMP_NUM_THREADS=1

module load intel
source {conda_base}

export PATH="{julia_path}"
cd {coderoot}
julia -e 'using Falcons; sim_det_scanfields("{tomlfile_path}")'
"""

with open(jobscript_path, 'w') as f:
    f.write(jobscript_data)

process = subprocess.Popen("jxsub "+jobscript_path, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(stdout_data, stderr_data) = process.communicate()

#print useful information
print("out: "+str(stdout_data).split('b\'')[1][:-3])
print("err: "+str(stderr_data).split('b\'')[1][:-3])
os.remove(jobscript_path)
jobid += 1
