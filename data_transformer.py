#!/bin/env/python
import numpy as np
def create_assignment_matrix(assignments):
    n_traj = np.shape(assignments.keys())[0]
    max_length = np.max([np.shape(assignments[i]) for i in assignments.keys()])
    assignment_matrix = np.zeros((n_traj,max_length))-1
    key_mapping={}
    for i,v in enumerate(assignments.keys()):
        current_traj_length = np.shape(assignments[v])[0]        
        assignment_matrix[i,:current_traj_length] = assignments[v]
        key_mapping[i]=v
    
    return key_mapping,assignment_matrix


def create_tics_array(assignments,kmeans_mdl,tica_data):
    n_traj = np.shape(assignments.keys())[0]
    max_length = np.max([np.shape(assignments[i]) for i in assignments.keys()])
    
    tics_to_use = kmeans_mdl.cluster_centers_.shape[1]
    tics_array = np.empty((n_traj,max_length,tics_to_use))
    tics_array[:,:,:] = np.NAN
    for index,trjname in enumerate(assignments.keys()):
        current_traj_length = np.shape(assignments[trjname])[0]
        for tic_index in range(tics_to_use):
            tics_array[index,:current_traj_length,tic_index] = tica_data[trjname][:,tic_index]
    return tics_array
