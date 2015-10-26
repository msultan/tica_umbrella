#!/bin/env/python

import os 
from msmbuilder.utils import verboseload,verbosedump
import mdtraj as mdt
import numpy as np
'''
script to load pertinent data for a given mutant
and perform some sanity checks 
'''
def load_frame(mutant,filename,frame_index):
    filename = os.path.splitext(filename)[0]
    return mdt.load("/nobackup/msultan/research/kinase/her_kinase/fah_data/%s/protein_traj/%s.hdf5"%(mutant,filename))[frame_index]


def sanity_test(msm_mdl,tica_data,kmeans_mdl,assignments):
    tics_to_use = kmeans_mdl.cluster_centers_.shape[1]
    for i,v in enumerate(tica_data.keys()[:20]):
        print v
        #skip 
        if not np.isnan(assignments[v]).any():

            assert((msm_mdl.transform(kmeans_mdl.transform([tica_data[v][:,:tics_to_use]]))==assignments[v]).all())
            trj = mdt.load("../../protein_traj/%s.hdf5"%v.split(".h5")[0])
            #print assignments[v],
            print trj.n_frames,assignments[v].shape[0]
            assert(trj.n_frames == assignments[v].shape[0])
    print "successful"



def load_mutant_data(mutant,sanity=True):
    mutant= "mut_wt"
    if os.getcwd()!="/nobackup/msultan/research/kinase/her_kinase/fah_data/%s/features/common_basis/"%mutant:
        os.chdir("/nobackup/msultan/research/kinase/her_kinase/fah_data/%s/features/common_basis"%mutant)
    base_dir = os.getcwd()
    wt_dir = os.path.join("/nobackup/msultan/research/kinase/her_kinase/fah_data/mut_wt/features/common_basis/")
    tica_mdl = verboseload("tica_obj.pkl")
    tica_data = verboseload("tica_data.pkl")
    kmeans_mdl = verboseload(wt_dir+"kmeans_obj.pkl")

    #need the fixed assignments because otherwise we will have issues 
    assignments = verboseload("fixed_assignments.pkl")
    msm_mdl = verboseload("msm_mdl.pkl")
    #some sanity tests
    if sanity:
        sanity_test(msm_mdl,tica_data,kmeans_mdl,assignments)
    return base_dir,wt_dir,msm_mdl,tica_mdl,tica_data,kmeans_mdl,assignments
    
 
