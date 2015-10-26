#!/bin/env/python

from msmbuilder.utils import verboseload,verbosedump
import mdtraj as mdt 
import numpy as np

import os 
from data_loader import load_frame, load_mutant_data
from data_transformer import *

from IPython import embed

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.nanargmin(np.abs(a - a0))
    return a.flat[idx]



def main():
    #load some of the basic data
    mutant="mut_wt"
    base_dir,wt_dir,msm_mdl,tica_mdl,tica_data,kmeans_mdl,assignments = load_mutant_data(mutant,sanity=False)
    #make sure we are
    os.chdir(wt_dir)
    key_mapping,assignment_matrix = create_assignment_matrix(assignments)
    tics_array = create_tics_array(assignments, kmeans_mdl, tica_data)
    
    #a few nice things to pre calculate 
    n_traj = assignment_matrix.shape[0]
    max_traj_len =assignment_matrix.shape[1]
    n_tics = kmeans_mdl.cluster_centers_.shape[1]
    tic_index = 0
    
    #figure out how much trajectories move through tic space in general 
    tic_ms_list = []
    for i in range(n_traj):
        len_of_current_traj = np.count_nonzero(np.isnan(tics_array[i,:,tic_index]))
        if len_of_current_traj:
            tic_ms_list.append(np.mean(tics_array[i,:-len_of_current_traj,tic_index]))
        else:
            tic_ms_list.append(np.mean(tics_array[i,:,tic_index])) 
    
    #get some statistics about the data 
    mean_tic_movement = np.mean(np.abs(tic_ms_list))
    std_tic_movement = np.std(np.abs(tic_ms_list))
    
    #figure out what bin width would work 
    bin_width = mean_tic_movement -2*std_tic_movement
    
    #get min and max of tic
    tic_min = np.nanmin(tics_array[:,:,tic_index])
    tic_max = np.nanmax(tics_array[:,:,tic_index])
    
    #get number of bins 
    n_bins = np.ceil((tic_max-tic_min)/bin_width)
    
    #get histogram centers    
    lin_place_points = np.linspace(tic_min,tic_max,n_bins)
    traj_list = []
    actual_tic_val_list=[]
    for i in lin_place_points:
        actual_tic_val = find_nearest(tics_array[:,:,tic_index],i)
        actual_tic_val_list.append(actual_tic_val)
    
        traj_index,frame_index = np.where(tics_array[:,:,tic_index]==actual_tic_val)
        traj_name = key_mapping[traj_index[0]]
       
        print i, actual_tic_val, traj_name,frame_index[0]
        traj_list.append(load_frame(mutant,traj_name,frame_index[0]))
    
    trj = traj_list[0]
    for i in traj_list[1:]:
        trj += i
    
    save_dir="/home/msultan/research/methods/tica_umbrella/"
    trj.save_xtc(save_dir+"test_case.xtc")
    trj[0].save_pdb(save_dir+"prot.pdb")
    verbosedump(tica_mdl,save_dir+"tica_obj.pkl")
    verbosedump(actual_tic_val_list,save_dir+"tic_val_list.pkl")
    return 

if __name__ == '__main__':
     main()
