from pyuvdata import UVData
import hera_pspec as hp
import numpy as np
import matplotlib.pyplot as plt
import copy, os, itertools, inspect
from hera_pspec.data import DATA_PATH

# select the data file to load
"""NEED TO CHANGE"""
dfile = os.path.join(DATA_PATH, 'zen.all.xx.LST.1.06964.uvA')

# Load into UVData objects
uvd = UVData()
uvd.read_miriad(dfile)

# # Check which baselines are available
# print(uvd.get_antpairs())

cosmo = hp.conversions.Cosmo_Conversions()
# print(cosmo)
beamfile = os.path.join(DATA_PATH, 'HERA_NF_dipole_power.beamfits')

# intantiate beam and pass cosmology, if not fed, a default Planck cosmology will be assumed
uvb = hp.pspecbeam.PSpecBeamUV(beamfile, cosmo=cosmo)

# find conversion factor from Jy to mK
Jy_to_mK = uvb.Jy_to_mK(np.unique(uvd.freq_array), pol='XX')

# reshape to appropriately match a UVData.data_array object and multiply in!
uvd.data_array *= Jy_to_mK[None, None, :, None]

uvd1 = uvd.select(times=np.unique(uvd.time_array)[:-1:2], inplace=False)
uvd2 = uvd.select(times=np.unique(uvd.time_array)[1::2], inplace=False)

ds = hp.PSpecData(dsets=[uvd1, uvd2], wgts=[None, None], beam=uvb)

ds.rephase_to_dset(0)

ds.dsets[0].vis_units = 'mK'
ds.dsets[1].vis_units = 'mK'

"""NEED TO CHANGE"""
baselines = [(24, 25), (37, 38), (38, 39)]

# Power spectra object
uvp = ds.pspec(baselines, baselines, (0, 1), [('xx', 'xx')], spw_ranges=[(600, 721)],
               input_data_weight='identity',
               norm='I', taper='blackman-harris', verbose=True)

# fig, ax = plt.subplots(figsize=(12, 8))

# spw = 1
# blp = ((24, 25), (24, 25))
"""NEED TO CHANGE"""
"""TRY TO DO IT FOR ALL THE KEYS"""
# print(list_of_keys.length())

# key = (spw, blp, 'xx')

""" Get all keys from UVP as a list and process each individual one. From each attained key parse out the spw values 
baselines and polarization pairs. This will be done through a for loop. 

Spw values will need to be parse out in order to attain the delays from each corresponding key. Once all delays 
are attained high delays can be parsed out. High delays are wanted because they have the least influence from 
systematic noise and other noise from equipment. Although thermal noise is still in the data set, because it is 
thermal is can be averaged out and thus trivial. 

From each key get the corresponding local sidereal times as that will be used to plot with the mean attained from high
delays. To get local sidereal times use the uvp.lst_avg_array method """

list_of_keys = uvp.get_all_keys()

list_of_spw = []
list_of_baselines = []
list_of_avg_powers = []

for key in list_of_keys:
    # Attains spectral window selection to each corresponding key
    list_of_spw.append(key[0])
    # Attains base line to each corresponding key
    list_of_baselines.append(key[1])
    # Gets avg power value to respective key
    power = np.abs(np.real(uvp.get_data(key)))
    # We want to average high delays for each time and each key therefore will
    # need to iterate through power variable which is a 2d array containing, time in the first
    # index and delay values in the second
    for delay_value in power:
        list_of_avg_powers.append(np.mean(delay_value))

# Want to reshape into a 3x3 as there are 3 averages that correspond to 3 times
# may want to revisit as shape is hard coded
# make it a double diminsion array n key x n time
list_of_avg_powers = np.array(list_of_avg_powers).reshape(3, 3)

# Attains time indices and sidereal time with each corresponding key
# only need to do it once because data from all baselines (keys) are recorded
# at the same time. Therefore in principal sidereal times should all be the same
time_index = uvp.key_to_indices(list_of_keys[0])[1]
sidereal_time = uvp.lst_avg_array[time_index]

# get the list of delays from corresponding spw values
list_of_delays = []
for spw in list_of_spw:
    delays = uvp.get_dlys(spw) * 1e9
    list_of_delays.append(delays)

# from the list of delays parse out the high delays
high_delay_indices = []
high_delays_list = []
for delays in list_of_delays:
    for delays_index, delays_value in np.ndenumerate(delays):
        if np.abs(delays_value) >= 2000:
            high_delay_indices.append(delays_index)
    high_delays = np.squeeze(delays[high_delay_indices])
    high_delays_list.append(high_delays)

# print(sidereal_time_list)
# print(list_of_avg_powers)
# print(len(sidereal_time_list))
# print(len(list_of_avg_powers))
# p1_big = ax.plot(high_delays, power_spectrum_high_delays.T)
# ax.set_yscale('log')
# ax.grid()
# ax.set_xlabel("delay [ns]", fontsize=14)
# ax.set_ylabel(r"$P(k)\ \rm [mK^2\ h^{-3}\ Mpc^3]$", fontsize=14)
# ax.set_title("spw : {}, blpair : {}, pol : {}".format(*key), fontsize=14)
#
fig, ax = plt.subplots(figsize=(12, 8))
p1_big_mean = ax.plot(sidereal_time, list_of_avg_powers)
# ax.set_yscale('log')
ax.grid()
ax.set_xlabel("Local Sidereal Time", fontsize=14)
ax.set_ylabel("Mean of High Delays", fontsize=14)
# ax.set_title("spw : {}, blpair : {}, pol : {}".format(*key), fontsize=14)
# plt.scatter(sidereal_time, list_of_avg_powers, marker='o')
# plt.xlabel("Local Sidereal Time", fontsize=14)
# plt.ylabel("Mean of High Delays", fontsize=14)
plt.show()
