import litebird_sim as lbs
import numpy as np
import subprocess
import sys
import os
import uuid
import time

#------ User dependent parameters ------
base_path = "<enter your path>"
julia_path = "<enter your path>"
coderoot  = "<enter your path>"
imo_path = "<enter your path>"
#----------------------------------------

logdir    = os.path.join(coderoot, 'log')
ancillary = os.path.join(coderoot, 'ancillary')
if not os.path.exists(base_path):
    os.makedirs(base_path)
if not os.path.exists(logdir):
    os.makedirs(logdir, exist_ok=True)
if not os.path.exists(ancillary):
    os.makedirs(ancillary, exist_ok=True)

# LSF job settings
jobq  = "l"

# general
imo_version = 'v2'

# simulation
name            = 'Scan field calculation'
nside           = 128
n_year          = 1
one_year        = 3600*24*365
duration_s      = n_year * one_year
division        = 320
sampling_rate   = 1.
coord           = "G"
hwp_rpm         = "IMO" # in None, rot. rate will be chosen from IMO
gamma           = 0

spin_n = [0,1,2,3,4,5,6,7,8,9,10]
spin_m = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]

imo = lbs.Imo(flatfile_location=imo_path)
telescopes = ["LFT", "MFT", "HFT"]

jobid = 0
for telescope in telescopes:
    inst = lbs.InstrumentInfo.from_imo(imo, f"/releases/{imo_version}/satellite/{telescope}/instrument_info",)
    for channel in inst.channel_names:
        channel_info = lbs.FreqChannelInfo.from_imo(imo, f"/releases/{imo_version}/satellite/{telescope}/{channel}/channel_info",)
        for det in channel_info.detector_names:
            if telescope == "LFT":
                gamma           = 270. # LFT has gamma=270. because IMO has wrong orientation. We have to modify here.
            else:
                gamma           = 0.
            if det.split("_")[-1] == "T":
                jobname = "scnfld_"+str(nside)+"_"+channel+"_"+"id_"+str(jobid)
                toml_filename   = str(uuid.uuid4())
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
spin_n = {spin_n}
spin_m = {spin_m}
duration_s = '{duration_s}'
division = '{division}'
hwp_rpm = '{hwp_rpm}'
coord = '{coord}'
                """
                with open(tomlfile_path, 'w') as f:
                    f.write(tomlfile_data)

                jobscript_path = os.path.join(ancillary, toml_filename+".lsf")
                jobscript_data = f"""#!/bin/zsh
#BSUB -q {jobq}
#BSUB -J {jobname}
#BSUB -o {logdir}/{jobname}.out
#BSUB -e {logdir}/{jobname}.err

eval "$(conda shell.bash hook)"
export PATH="{julia_path}"
cd {coderoot}
julia -e 'using Falcons; sim_det_scanfields("{tomlfile_path}")'
                """

                with open(jobscript_path, 'w') as f:
                    f.write(jobscript_data)

                process = subprocess.Popen("bsub < " + jobscript_path, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                (stdout_data, stderr_data) = process.communicate()
                #print useful information
                print("out: "+str(stdout_data).split('b\'')[1][:-3])
                print("err: "+str(stderr_data).split('b\'')[1][:-3])
                os.remove(jobscript_path)
                jobid += 1
