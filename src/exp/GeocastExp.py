"""
Contains all experimental methods
"""
import random

import numpy as np
import time
import datetime
import logging
import math
import sys
from datetime import datetime

#sys.path.append('.')
sys.path.append('../geocast')
sys.path.append('../icde12')
sys.path.append('../localness')
sys.path.append('../htree')

from Params import Params   
from PSDExp import data_readin
from multiprocessing import Process
 
from Grid_adaptive import Grid_adaptive
from Grid_adaptive_localness import Grid_adaptive_localness

from Geocast import geocast, post_geocast
from GeocastLog import geocast_log
from GeocastInfo import GeocastInfo

from GeocastNaive import geocast_naive
from GeocastKNN import geocast_knn
from Utils import is_rect_cover, performed_tasks
import os.path

eps_list = [.1,.4,.7, 1.0]
#eps_list = [0.1]

seed_list = [9110, 4064, 6903, 7509, 5342, 3230, 3584, 7019, 3564, 6456]

def tasks_gen(data, taskNo, x1=-124.8193, y1=31.3322, x2=-103.0020, y2=49.0025):
    """
    Generate random points within dataset
    """
    all_points = []
    if os.path.isfile(Params.TASKPATH):
        with open(Params.TASKPATH) as f:
            content = f.readlines()
        for i in range(len(seed_list)):
            ran_points = []
            for j in range(taskNo):
                ran_points.append(map(float, content[i*taskNo + j].split()))
            all_points.append(ran_points)
    else:
        tasks = ""
        logging.debug('tasks_gen: generating tasks...')
        
        boundary = np.array([[x1,y1],[x2,y2]])
        for seed in seed_list:
            ran_points = []
            np.random.seed(seed)
            count = 0
            while count < taskNo:
                idx = np.random.randint(0,data.shape[1])
                _ran_point = data[:,idx]
                if is_rect_cover(boundary, _ran_point):
                    ran_points.append(_ran_point)
                    count = count + 1
            all_points.append(ran_points)
            for item in ran_points:
                tasks = tasks + ("%s\n" % " " . join(map(str, item)))
        outfile = open(Params.TASKPATH, "w")
        outfile.write(tasks)
        outfile.close()
    return all_points

def sort_data(data):
    """
    Get true answer by linear search along each dimension
    """
    ndim = data.shape[0]
    for dim in range(ndim):
        if data.shape[1] == 0:
            break
        idx = np.argsort(data[dim,:],kind='mergesort')
        data[:,:] = data[:,idx]

def evalGeocast_Utility(data, all_tasks):
    """
    Varying expect utility U and privacy budget eps
    """
    logging.info("evalGeocast_Utility")
    exp_name = "Geocast_Utility"
    temp_U = Params.U
    U_list = [0.6, 0.7, 0.8, 0.9]
    
    res_cube_anw = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    res_cube_atd = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    res_cube_atd_fcfs = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    res_cube_appt = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    res_cube_cell = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    res_cube_cmp = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    res_cube_hop = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    res_cube_hop2 = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    res_cube_cov = np.zeros((len(eps_list),len(seed_list),len(U_list)))
    
    for j in range(len(seed_list)):
        for i in range(len(eps_list)):
            log_str = ""

            Params.CUSTOMIZED_GRANULARITY = True
            Params.PARTIAL_CELL_SELECTION = True
            p = Params(seed_list[j])
            p.Eps = eps_list[i]
            tree = Grid_adaptive(data, p)
            tree.buildIndex()
            
            for k in range(len(U_list)):
                Params.U = U_list[k]

                totalANW, totalATD, totalATD_FCFS, totalCell = 0, 0, 0, 0
                totalPerformedTasks, totalCompactness = 0, 0
                totalHop_GDY, totalHop2_GDY, totalCov_GDY = 0, 0, 0
                
                tasks = all_tasks[j]
                for l in range(len(tasks)):
                    if (l+1)%Params.LOGGING_STEPS == 0:
                        print ">> " + str(l+1) + " tasks completed"                    
                    t = tasks[l]
                    q, q_log = geocast(tree, t, p.Eps)
                    no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                    performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
                    if performed:
                        totalPerformedTasks += 1
                        totalANW += no_workers
                        totalATD += dist
                        totalCell += len(Cells)
                        totalCompactness += q_log[-1][3]
                        totalHop_GDY += no_hops
                        totalHop2_GDY += no_hops2
                        totalCov_GDY += coverage 
                        
                        # logging
                        if Params.GEOCAST_LOG:
                            info = GeocastInfo(1, t, Cells)
                            log_str = log_str + str(info.logging()) + "\n"
		    else:   # the task is not performed
                        if Params.GEOCAST_LOG:
                            info = GeocastInfo(0, t, Cells)
                            log_str = log_str + str(info.logging()) + "\n"
                        
                    # FCFS
                    performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
                    if performed:
                        totalATD_FCFS += dist

                ANW = (totalANW + 0.0)/totalPerformedTasks
                ATD = totalATD/totalPerformedTasks
                ATD_FCFS = totalATD_FCFS/totalPerformedTasks
                ASC = (totalCell + 0.0)/totalPerformedTasks
                APPT = 100*float(totalPerformedTasks)/Params.TASK_NO
                CMP = totalCompactness/totalPerformedTasks
                HOP_GDY = (totalHop_GDY + 0.0)/Params.TASK_NO
                HOP2_GDY = (totalHop2_GDY + 0.0)/Params.TASK_NO
                COV_GDY = 100*(totalCov_GDY + 0.0)/Params.TASK_NO
                
                res_cube_anw[i,j,k] = ANW
                res_cube_atd[i,j,k] = ATD
                res_cube_atd_fcfs[i,j,k] = ATD_FCFS
                res_cube_appt[i,j,k] = APPT
                res_cube_cell[i,j,k] = ASC
                res_cube_cmp[i,j,k] = CMP                
                res_cube_hop[i,j,k] = HOP_GDY
                res_cube_hop2[i,j,k] = HOP2_GDY
                res_cube_cov[i,j,k] = COV_GDY   
                
            if Params.GEOCAST_LOG:
                geocast_log(exp_name, log_str, p.Eps)
                
    res_summary_anw = np.average(res_cube_anw,axis=1)
    np.savetxt(Params.resdir+exp_name + '_anw_' + `Params.TASK_NO`, res_summary_anw, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd_fcfs,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_fcfs_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_appt = np.average(res_cube_appt,axis=1)
    np.savetxt(Params.resdir+exp_name + '_appt_' + `Params.TASK_NO`, res_summary_appt, fmt='%.4f\t')
    res_summary_cell = np.average(res_cube_cell,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cell_' + `Params.TASK_NO`, res_summary_cell, fmt='%.4f\t')
    res_summary_cmp = np.average(res_cube_cmp,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cmp_' + `Params.TASK_NO`, res_summary_cmp, fmt='%.4f\t')
    res_summary_hop = np.average(res_cube_hop,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop_' + `Params.TASK_NO`, res_summary_hop, fmt='%.4f\t')
    res_summary_hop2 = np.average(res_cube_hop2,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop2_' + `Params.TASK_NO`, res_summary_hop2, fmt='%.4f\t')    
    res_summary_cov = np.average(res_cube_cov,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cov_' + `Params.TASK_NO`, res_summary_cov, fmt='%.4f\t')
    
    Params.U = temp_U

def evalGeocast_MAR(data, all_tasks):
    """
    Varying max acceptance rate MAR and privacy budget eps
    """
    logging.info("evalGeocast_MAR")
    exp_name = "Geocast_MAR"
    temp_MAR = Params.MAR
    MAR_list = [0.1, 0.4, 0.7, 1.0]
#    MAR_list = [0.5]
       
    res_cube_anw = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    res_cube_atd = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    res_cube_atd_fcfs = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    res_cube_appt = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    res_cube_cell = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    res_cube_cmp = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    res_cube_hop = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    res_cube_hop2 = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    res_cube_cov = np.zeros((len(eps_list),len(seed_list),len(MAR_list)))
    
    for j in range(len(seed_list)):
        for i in range(len(eps_list)):
            log_str = ""
            Params.CUSTOMIZED_GRANULARITY = True
            Params.PARTIAL_CELL_SELECTION = True
            p = Params(seed_list[j])
            p.Eps = eps_list[i]
            tree = Grid_adaptive(data, p)
            tree.buildIndex()
            for k in range(len(MAR_list)):
                Params.MAR = MAR_list[k]
                totalANW, totalATD, totalATD_FCFS, totalCell = 0, 0, 0, 0
                totalPerformedTasks, totalCompactness = 0, 0
                totalHop_GDY, totalHop2_GDY, totalCov_GDY = 0, 0, 0

                tasks = all_tasks[j]
                for l in range(len(tasks)):
                    if (l+1)%Params.LOGGING_STEPS == 0:
                        print ">> " + str(l+1) + " tasks completed"
                    t = tasks[l]
                    q, q_log = geocast(tree, t, p.Eps)
                    no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                    performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
                    if performed:
                        totalPerformedTasks += 1
                        totalANW += no_workers
                        totalATD += dist
                        totalCell += len(Cells)
                        totalCompactness += q_log[-1][3]
                        totalHop_GDY += no_hops
                        totalHop2_GDY += no_hops2
                        totalCov_GDY += coverage 
                        
                        # logging
                        if Params.GEOCAST_LOG:
                            info = GeocastInfo(1, t, Cells)
                            log_str = log_str + str(info.logging()) + "\n"
		    else:   # the task is not performed
                        if Params.GEOCAST_LOG:
                            info = GeocastInfo(0, t, Cells)
                            log_str = log_str + str(info.logging()) + "\n"
                        
                    # FCFS
                    performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
                    if performed:
                        totalATD_FCFS += dist

                ANW = (totalANW + 0.0)/totalPerformedTasks
                ATD = totalATD/totalPerformedTasks
                ATD_FCFS = totalATD_FCFS/totalPerformedTasks
                ASC = (totalCell + 0.0)/totalPerformedTasks
                APPT = 100*float(totalPerformedTasks)/Params.TASK_NO
                CMP = totalCompactness/totalPerformedTasks
                HOP_GDY = (totalHop_GDY + 0.0)/Params.TASK_NO
                HOP2_GDY = (totalHop2_GDY + 0.0)/Params.TASK_NO
                COV_GDY = 100*(totalCov_GDY + 0.0)/Params.TASK_NO
            
                res_cube_anw[i,j,k] = ANW
                res_cube_atd[i,j,k] = ATD
                res_cube_atd_fcfs[i,j,k] = ATD_FCFS
                res_cube_appt[i,j,k] = APPT
                res_cube_cell[i,j,k] = ASC
                res_cube_cmp[i,j,k] = CMP
                res_cube_hop[i,j,k] = HOP_GDY
                res_cube_hop2[i,j,k] = HOP2_GDY
                res_cube_cov[i,j,k] = COV_GDY            
            
            if Params.GEOCAST_LOG:
                geocast_log(exp_name, log_str, p.Eps)
            
    res_summary_anw = np.average(res_cube_anw,axis=1)
    np.savetxt(Params.resdir+exp_name + '_anw_' + `Params.TASK_NO`, res_summary_anw, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd_fcfs,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_fcfs_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_appt = np.average(res_cube_appt,axis=1)
    np.savetxt(Params.resdir+exp_name + '_appt_' + `Params.TASK_NO`, res_summary_appt, fmt='%.4f\t')
    res_summary_cell = np.average(res_cube_cell,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cell_' + `Params.TASK_NO`, res_summary_cell, fmt='%.4f\t')
    res_summary_cmp = np.average(res_cube_cmp,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cmp_' + `Params.TASK_NO`, res_summary_cmp, fmt='%.4f\t')
    res_summary_hop = np.average(res_cube_hop,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop_' + `Params.TASK_NO`, res_summary_hop, fmt='%.4f\t')
    res_summary_hop2 = np.average(res_cube_hop2,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop2_' + `Params.TASK_NO`, res_summary_hop2, fmt='%.4f\t')    
    res_summary_cov = np.average(res_cube_cov,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cov_' + `Params.TASK_NO`, res_summary_cov, fmt='%.4f\t')
    
    Params.MAR = temp_MAR

def evalGeocast_GRA_PAR(data, all_tasks):
    """
    Evaluate the baseline approach (GREEDY) with improved approaches, Customized Granularity (GRA) and Partial Cell Selection (PAR)
    """
    logging.info("evalGeocast_GRA_PAR")
    exp_name = "Geocast_GRA_PAR"
    methodList = ["GDY", "GRA", "PAR", "GRA_PAR"]
    
    res_cube_anw = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_atd = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_atd_fcfs = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_appt = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cell = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cmp = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_hop = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_hop2 = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cov = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    
    for j in range(len(seed_list)):
        for i in range(len(eps_list)):
            # GDY tree
            Params.CUSTOMIZED_GRANULARITY = False
            p = Params(seed_list[j])
            p.Eps = eps_list[i]
            tree_GDY = Grid_adaptive(data, p)
            tree_GDY.buildIndex()
            
            # GRA tree
            Params.CUSTOMIZED_GRANULARITY = True
            p = Params(seed_list[j])
            p.Eps = eps_list[i]
            tree_GRA = Grid_adaptive(data, p)
            tree_GRA.buildIndex()
            
            totalANW_GDY, totalANW_GRA, totalANW_PAR, totalANW_GRA_PAR= 0, 0, 0, 0
            totalATD_GDY, totalATD_FCFS_GDY, totalATD_GRA, totalATD_FCFS_GRA, totalATD_PAR, totalATD_FCFS_PAR, totalATD_GRA_PAR, totalATD_FCFS_GRA_PAR = 0, 0, 0, 0, 0, 0, 0, 0
            totalCell_GDY, totalCell_GRA, totalCell_PAR, totalCell_GRA_PAR = 0, 0, 0, 0
	    totalCompactness_GDY, totalCompactness_GRA, totalCompactness_PAR, totalCompactness_GRA_PAR = 0, 0, 0, 0
            totalPerformedTasks_GDY, totalPerformedTasks_GRA, totalPerformedTasks_PAR, totalPerformedTasks_GRA_PAR = 0, 0, 0, 0
            totalHop_GDY, totalHop_GRA, totalHop_PAR, totalHop_GRA_PAR = 0, 0, 0, 0
            totalHop2_GDY, totalHop2_GRA, totalHop2_PAR, totalHop2_GRA_PAR = 0, 0, 0, 0
            totalCov_GDY, totalCov_GRA, totalCov_PAR, totalCov_GRA_PAR = 0, 0, 0, 0
            
            tasks = all_tasks[j]
            for l in range(len(tasks)):
                if (l+1)%Params.LOGGING_STEPS == 0:
                    print ">> " + str(l+1) + " tasks completed"
                    
                t = tasks[l]
                
                # GDY
                Params.PARTIAL_CELL_SELECTION = False
		q, q_log = geocast(tree_GDY, t, p.Eps)
		no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
		if performed:
                    totalPerformedTasks_GDY += 1
                    totalANW_GDY += no_workers
                    totalATD_GDY += dist
                    totalCell_GDY += len(Cells)
		    totalCompactness_GDY += q_log[-1][3]
                    totalHop_GDY += no_hops
                    totalHop2_GDY += no_hops2
                    totalCov_GDY += coverage                    
		performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
		if performed:
		    totalATD_FCFS_GDY += dist


                # GRA
                Params.PARTIAL_CELL_SELECTION = False               
                q, q_log = geocast(tree_GRA, t, p.Eps)
		no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
		if performed:
                    totalPerformedTasks_GRA += 1
                    totalANW_GRA += no_workers
                    totalATD_GRA += dist
                    totalCell_GRA += len(Cells)
		    totalCompactness_GRA += q_log[-1][3]
                    totalHop_GRA += no_hops
                    totalHop2_GRA += no_hops2
                    totalCov_GRA += coverage                        
		performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
		if performed:
		    totalATD_FCFS_GRA += dist

                # PAR
                Params.PARTIAL_CELL_SELECTION = True
                q, q_log = geocast(tree_GDY, t, p.Eps)
		no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
		if performed:
                    totalPerformedTasks_PAR += 1
                    totalANW_PAR += no_workers
                    totalATD_PAR += dist
                    totalCell_PAR += len(Cells)
		    totalCompactness_PAR += q_log[-1][3]
                    totalHop_PAR += no_hops
                    totalHop2_PAR += no_hops2
                    totalCov_PAR += coverage                      
		performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
		if performed:
		    totalATD_FCFS_PAR += dist
                    
                # GRA_PAR
                Params.PARTIAL_CELL_SELECTION = True              
                q, q_log = geocast(tree_GRA, t, p.Eps)
		no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
		if performed:
                    totalPerformedTasks_GRA_PAR += 1
                    totalANW_GRA_PAR += no_workers
                    totalATD_GRA_PAR += dist
                    totalCell_GRA_PAR += len(Cells)
		    totalCompactness_GRA_PAR += q_log[-1][3]
                    totalHop_GRA_PAR += no_hops
                    totalHop2_GRA_PAR += no_hops2
                    totalCov_GRA_PAR += coverage                     
		performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
		if performed:
		    totalATD_FCFS_GRA_PAR += dist
                    
            # GDY 
            ANW_GDY = (totalANW_GDY + 0.0)/totalPerformedTasks_GDY
            ATD_GDY = totalATD_GDY/totalPerformedTasks_GDY
	    ATD_FCFS_GDY = totalATD_FCFS_GDY/totalPerformedTasks_GDY
            ASC_GDY = (totalCell_GDY + 0.0)/totalPerformedTasks_GDY
	    CMP_GDY = totalCompactness_GDY/totalPerformedTasks_GDY
            APPT_GDY = 100*float(totalPerformedTasks_GDY)/Params.TASK_NO
            HOP_GDY = (totalHop_GDY + 0.0)/Params.TASK_NO
            HOP2_GDY = (totalHop2_GDY + 0.0)/Params.TASK_NO
            COV_GDY = 100*(totalCov_GDY + 0.0)/Params.TASK_NO

            # GRA
            ANW_GRA = (totalANW_GRA + 0.0)/totalPerformedTasks_GRA
            ATD_GRA = totalATD_GRA/totalPerformedTasks_GRA
	    ATD_FCFS_GRA = totalATD_FCFS_GRA/totalPerformedTasks_GRA
            ASC_GRA = (totalCell_GRA + 0.0)/totalPerformedTasks_GRA
	    CMP_GRA = totalCompactness_GRA/totalPerformedTasks_GRA
            APPT_GRA = 100*float(totalPerformedTasks_GRA)/Params.TASK_NO
            HOP_GRA = (totalHop_GRA + 0.0)/Params.TASK_NO
            HOP2_GRA = (totalHop2_GRA + 0.0)/Params.TASK_NO
            COV_GRA = 100*(totalCov_GRA + 0.0)/Params.TASK_NO
            	    
            # PAR
            ANW_PAR = (totalANW_PAR + 0.0)/totalPerformedTasks_PAR
            ATD_PAR = totalATD_PAR/totalPerformedTasks_PAR
	    ATD_FCFS_PAR = totalATD_FCFS_PAR/totalPerformedTasks_PAR
            ASC_PAR = (totalCell_PAR + 0.0)/totalPerformedTasks_PAR
	    CMP_PAR = totalCompactness_PAR/totalPerformedTasks_PAR
            APPT_PAR = 100*float(totalPerformedTasks_PAR)/Params.TASK_NO
            HOP_PAR = (totalHop_PAR + 0.0)/Params.TASK_NO
            HOP2_PAR = (totalHop2_PAR + 0.0)/Params.TASK_NO
            COV_PAR = 100*(totalCov_PAR + 0.0)/Params.TASK_NO
            
            # GRA_PAR
            ANW_GRA_PAR = (totalANW_GRA_PAR + 0.0)/totalPerformedTasks_GRA_PAR
            ATD_GRA_PAR = totalATD_GRA_PAR/totalPerformedTasks_GRA_PAR
	    ATD_FCFS_GRA_PAR = totalATD_FCFS_GRA_PAR/totalPerformedTasks_GRA_PAR
            ASC_GRA_PAR = (totalCell_GRA_PAR + 0.0)/totalPerformedTasks_GRA_PAR
	    CMP_GRA_PAR = totalCompactness_GRA_PAR/totalPerformedTasks_GRA_PAR
            APPT_GRA_PAR = 100*float(totalPerformedTasks_GRA_PAR)/Params.TASK_NO
            HOP_GRA_PAR = (totalHop_GRA_PAR + 0.0)/Params.TASK_NO
            HOP2_GRA_PAR = (totalHop2_GRA_PAR + 0.0)/Params.TASK_NO
            COV_GRA_PAR = 100*(totalCov_GRA_PAR + 0.0)/Params.TASK_NO
            
	    res_cube_anw[i,j,0] = ANW_GDY
            res_cube_atd[i,j,0] = ATD_GDY
            res_cube_atd_fcfs[i,j,0] = ATD_FCFS_GDY
            res_cube_appt[i,j,0] = APPT_GDY
            res_cube_cell[i,j,0] = ASC_GDY
            res_cube_cmp[i,j,0] = CMP_GDY
            res_cube_hop[i,j,0] = HOP_GDY
            res_cube_hop2[i,j,0] = HOP2_GDY
            res_cube_cov[i,j,0] = COV_GDY
            
	    res_cube_anw[i,j,1] = ANW_GRA
            res_cube_atd[i,j,1] = ATD_GRA
            res_cube_atd_fcfs[i,j,1] = ATD_FCFS_GRA
            res_cube_appt[i,j,1] = APPT_GRA
            res_cube_cell[i,j,1] = ASC_GRA
            res_cube_cmp[i,j,1] = CMP_GRA
            res_cube_hop[i,j,1] = HOP_GRA
            res_cube_hop2[i,j,1] = HOP2_GRA
            res_cube_cov[i,j,1] = COV_GRA

	    res_cube_anw[i,j,2] = ANW_PAR
            res_cube_atd[i,j,2] = ATD_PAR
            res_cube_atd_fcfs[i,j,2] = ATD_FCFS_PAR
            res_cube_appt[i,j,2] = APPT_PAR
            res_cube_cell[i,j,2] = ASC_PAR
            res_cube_cmp[i,j,2] = CMP_PAR
            res_cube_hop[i,j,2] = HOP_PAR
            res_cube_hop2[i,j,2] = HOP2_PAR
            res_cube_cov[i,j,2] = COV_PAR
            
	    res_cube_anw[i,j,3] = ANW_GRA_PAR
            res_cube_atd[i,j,3] = ATD_GRA_PAR
            res_cube_atd_fcfs[i,j,3] = ATD_FCFS_GRA_PAR
            res_cube_appt[i,j,3] = APPT_GRA_PAR
            res_cube_cell[i,j,3] = ASC_GRA_PAR
            res_cube_cmp[i,j,3] = CMP_GRA_PAR
            res_cube_hop[i,j,3] = HOP_GRA_PAR
            res_cube_hop2[i,j,3] = HOP2_GRA_PAR
            res_cube_cov[i,j,3] = COV_GRA_PAR        
            
    res_summary_anw = np.average(res_cube_anw,axis=1)
    np.savetxt(Params.resdir+exp_name + '_anw_' + `Params.TASK_NO`, res_summary_anw, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd_fcfs,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_fcfs_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_appt = np.average(res_cube_appt,axis=1)
    np.savetxt(Params.resdir+exp_name + '_appt_' + `Params.TASK_NO`, res_summary_appt, fmt='%.4f\t')
    res_summary_cell = np.average(res_cube_cell,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cell_' + `Params.TASK_NO`, res_summary_cell, fmt='%.4f\t')
    res_summary_cmp = np.average(res_cube_cmp,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cmp_' + `Params.TASK_NO`, res_summary_cmp, fmt='%.4f\t')
    res_summary_hop = np.average(res_cube_hop,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop_' + `Params.TASK_NO`, res_summary_hop, fmt='%.4f\t')
    res_summary_hop2 = np.average(res_cube_hop2,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop2_' + `Params.TASK_NO`, res_summary_hop2, fmt='%.4f\t')    
    res_summary_cov = np.average(res_cube_cov,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cov_' + `Params.TASK_NO`, res_summary_cov, fmt='%.4f\t') 

def evalGeocast_Baseline(data, all_tasks):
    """
    Evaluate geocast algorithm in privacy mode and in non-privacy mode
    """
    logging.info("evalGeocast_Baseline")
    exp_name = "Geocast_Baseline"
    methodList = ["Geocast", "Baseline"]
    
    res_cube_anw = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_atd = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_atd_fcfs = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_appt = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cell = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cmp = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_hop = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_hop2 = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cov = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    
    for j in range(len(seed_list)):
        for i in range(len(eps_list)):
            # Geocast tree
            Params.CUSTOMIZED_GRANULARITY = True
	    Params.PARTIAL_CELL_SELECTION = True  
	    Params.COST_FUNCTION = "hybrid"	    
            p = Params(seed_list[j])
            p.Eps = eps_list[i]
            tree = Grid_adaptive(data, p)
            tree.buildIndex()
            
            totalANW_Naive, totalANW_Geocast, totalANW_Knn = 0, 0, 0
            totalATD_Naive, totalATD_Naive_FCFS, totalATD_Geocast, totalATD_FCFS_Geocast, totalATD_Knn, totalATD_Knn_FCFS = 0, 0, 0, 0, 0, 0
            totalCell_Geocast = 0
	    totalCompactness_Geocast = 0
            totalPerformedTasks_Naive, totalPerformedTasks_Geocast, totalPerformedTasks_Knn = 0, 0, 0
            totalHop_Geocast, totalHop_Knn = 0, 0
            totalHop2_Geocast, totalHop2_Knn = 0, 0
            totalCov_Geocast, totalCov_Knn = 0, 0
            
            tasks = all_tasks[j]
            for l in range(len(tasks)):
                if (l+1)%Params.LOGGING_STEPS == 0:
                    print ">> " + str(l+1) + " tasks completed"
                t = tasks[l]
                                        
                # Naive approach
#                workers_naive, performed, dist_naive = geocast_naive(tree, data, t, False, Params.U)
#                if performed:
#                    totalPerformedTasks_Naive +=1
#                    totalANW_Naive += no_workers_naive
#                    totalATD_Naive += dist_naive
#                workers_naive, performed, dist_naive_FCFS = geocast_naive(tree, data, t, True, Params.U)
#                if performed:
#                    totalATD_Naive_FCFS += dist_naive_FCFS
                    
                # Geocast
                q, q_log = geocast(tree, t, p.Eps)
		no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
                if performed:
                    totalPerformedTasks_Geocast += 1
                    totalANW_Geocast += no_workers
                    totalATD_Geocast += dist
                    totalCell_Geocast += len(Cells)
		    totalCompactness_Geocast += q_log[-1][3]
                    totalHop_Geocast += no_hops
                    totalHop2_Geocast += no_hops2
                    totalCov_Geocast += coverage
		performed, worker, dist_fcfs = performed_tasks(workers, Params.MTD, t, True)
		if performed:
		    totalATD_FCFS_Geocast += dist_fcfs
		    
                # Baseline (no privacy)
                no_workers_knn, performed, dist_knn, dist_knn_FCFS, no_hops, coverage, no_hops2 = geocast_knn(data, t)
                if performed:
                    totalPerformedTasks_Knn +=1
                    totalANW_Knn += no_workers_knn
                    totalATD_Knn += dist_knn
                    totalATD_Knn_FCFS += dist_knn_FCFS
                    totalHop_Knn += no_hops
                    totalHop2_Knn += no_hops2
                    totalCov_Knn += coverage                    

            
            # Naive
#	    ANW_naive = (totalANW_Naive + 0.0)/totalPerformedTasks_Naive
#            ATD_naive = totalATD_Naive/totalPerformedTasks_Naive
#	    ATD_naive_FCFS = totalATD_Naive_FCFS/totalPerformedTasks_Naive
#            APPT_naive = 100*float(totalPerformedTasks_Naive)/Params.TASK_NO
            
            # Geocast
            ANW_Geocast = (totalANW_Geocast + 0.0)/totalPerformedTasks_Geocast
            ATD_Geocast = totalATD_Geocast/totalPerformedTasks_Geocast
	    ATD_FCFS_Geocast = totalATD_FCFS_Geocast/totalPerformedTasks_Geocast
            ASC_Geocast = (totalCell_Geocast + 0.0)/totalPerformedTasks_Geocast
	    CMP_Geocast = totalCompactness_Geocast/totalPerformedTasks_Geocast
            APPT_Geocast = 100*float(totalPerformedTasks_Geocast)/Params.TASK_NO
            HOP_Geocast = float(totalHop_Geocast)/Params.TASK_NO
            HOP2_Geocast = float(totalHop2_Geocast)/Params.TASK_NO
            COV_Geocast = 100*float(totalCov_Geocast)/Params.TASK_NO

            # Baseline
            ANW_Knn = (totalANW_Knn + 0.0)/totalPerformedTasks_Knn
            ATD_Knn = totalATD_Knn/totalPerformedTasks_Knn
            ATD_FCFS_Knn = totalATD_Knn_FCFS/totalPerformedTasks_Knn
            APPT_Knn = 100*float(totalPerformedTasks_Knn)/Params.TASK_NO
            HOP_Knn = float(totalHop_Knn)/Params.TASK_NO
            HOP2_Knn = float(totalHop2_Knn)/Params.TASK_NO
            COV_Knn = 100*float(totalCov_Knn)/Params.TASK_NO            
            
#	    res_cube_anw[i,j,0] = ANW_naive
#            res_cube_atd[i,j,0] = ATD_naive
#            res_cube_atd_fcfs[i,j,0] = ATD_naive_FCFS
#            res_cube_appt[i,j,0] = APPT_naive
#            res_cube_cell[i,j,0] = 0
#            res_cube_cmp[i,j,0] = 0
            
	    res_cube_anw[i,j,0] = ANW_Geocast
            res_cube_atd[i,j,0] = ATD_Geocast
            res_cube_atd_fcfs[i,j,0] = ATD_FCFS_Geocast
            res_cube_appt[i,j,0] = APPT_Geocast
            res_cube_cell[i,j,0] = ASC_Geocast
            res_cube_cmp[i,j,0] = CMP_Geocast
            res_cube_hop[i,j,0] = HOP_Geocast
            res_cube_hop2[i,j,0] = HOP2_Geocast
            res_cube_cov[i,j,0] = COV_Geocast
            
            res_cube_anw[i,j,1] = ANW_Knn
            res_cube_atd[i,j,1] = ATD_Knn
            res_cube_atd_fcfs[i,j,1] = ATD_FCFS_Knn
            res_cube_appt[i,j,1] = APPT_Knn
            res_cube_cell[i,j,1] = 0
            res_cube_cmp[i,j,1] = 0
            res_cube_hop[i,j,1] = HOP_Knn
            res_cube_hop2[i,j,1] = HOP2_Knn
            res_cube_cov[i,j,1] = COV_Knn        
            
    res_summary_anw = np.average(res_cube_anw,axis=1)
    np.savetxt(Params.resdir+exp_name + '_anw_' + `Params.TASK_NO`, res_summary_anw, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd_fcfs,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_fcfs_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_appt = np.average(res_cube_appt,axis=1)
    np.savetxt(Params.resdir+exp_name + '_appt_' + `Params.TASK_NO`, res_summary_appt, fmt='%.4f\t')
    res_summary_cell = np.average(res_cube_cell,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cell_' + `Params.TASK_NO`, res_summary_cell, fmt='%.4f\t')
    res_summary_cmp = np.average(res_cube_cmp,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cmp_' + `Params.TASK_NO`, res_summary_cmp, fmt='%.4f\t')
    res_summary_hop = np.average(res_cube_hop,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop_' + `Params.TASK_NO`, res_summary_hop, fmt='%.4f\t')
    res_summary_hop2 = np.average(res_cube_hop2,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop2_' + `Params.TASK_NO`, res_summary_hop2, fmt='%.4f\t')    
    res_summary_cov = np.average(res_cube_cov,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cov_' + `Params.TASK_NO`, res_summary_cov, fmt='%.4f\t')    


def evalGeocast_Compactness(data, all_tasks):
    """
    Evaluate the geocast algorith under different cost functions, 
    utility-based heuristic, compactness-based heuristic and hybrid heuristic
    """
    logging.info("evalGeocast_Compactness")
    exp_name = "Geocast_Compactness"
    methodList = ["utility", "hybrid", "compactness"]
    
    res_cube_anw = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_atd = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_atd_fcfs = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_appt = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cell = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cmp = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_hop = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_hop2 = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    res_cube_cov = np.zeros((len(eps_list),len(seed_list),len(methodList)))
    
    Params.CUSTOMIZED_GRANULARITY = True
    Params.PARTIAL_CELL_SELECTION = True
            
    for j in range(len(seed_list)):
        for i in range(len(eps_list)):
            p = Params(seed_list[j])
            p.Eps = eps_list[i]
            tree = Grid_adaptive(data, p)
            tree.buildIndex()
            
            totalANW_Utility, totalANW_Hybrid, totalANW_Compactness = 0, 0, 0
            totalATD_Utility, totalATD_FCFS_Utility, totalATD_Hybrid, totalATD_FCFS_Hybrid, totalATD_Compactness, totalATD_FCFS_Compactness = 0, 0, 0, 0, 0, 0
            totalCell_Utility, totalCell_Hybrid, totalCell_Compactness = 0, 0, 0
	    totalCompactness_Utility, totalCompactness_Hybrid, totalCompactness_Compactness = 0, 0, 0
            totalPerformedTasks_Utility, totalPerformedTasks_Hybrid, totalPerformedTasks_Compactness= 0, 0, 0
            totalHop_Utility, totalHop_Hybrid, totalHop_Compactness = 0, 0, 0
            totalHop2_Utility, totalHop2_Hybrid, totalHop2_Compactness = 0, 0, 0
            totalCov_Utility, totalCov_Hybrid, totalCov_Compactness = 0, 0, 0
            
            tasks = all_tasks[j]
            for l in range(len(tasks)):
                if (l+1)%Params.LOGGING_STEPS == 0:
                    print ">> " + str(l+1) + " tasks completed"
                t = tasks[l]
                
                # Utility
                Params.COST_FUNCTION = "utility"
		q, q_log = geocast(tree, t, p.Eps)
		no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
		if performed:
                    totalPerformedTasks_Utility += 1
                    totalANW_Utility += no_workers
                    totalATD_Utility += dist
                    totalCell_Utility += len(Cells)
		    totalCompactness_Utility += q_log[-1][3]
                    totalHop_Utility += no_hops
                    totalHop2_Utility += no_hops2
                    totalCov_Utility += coverage
		performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
		if performed:
		    totalATD_FCFS_Utility += dist

                # Hybrid
                Params.COST_FUNCTION = "hybrid"           
                q, q_log = geocast(tree, t, p.Eps)
		no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
		if performed:
                    totalPerformedTasks_Hybrid += 1
                    totalANW_Hybrid += no_workers
                    totalATD_Hybrid += dist
                    totalCell_Hybrid += len(Cells)
		    totalCompactness_Hybrid += q_log[-1][3]
                    totalHop_Hybrid += no_hops
                    totalHop2_Hybrid += no_hops2
                    totalCov_Hybrid += coverage
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
		if performed:
		    totalATD_FCFS_Hybrid += dist

                # Compactness
                Params.COST_FUNCTION = "compactness"
                q, q_log = geocast(tree, t, p.Eps)
		no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
		if performed:
                    totalPerformedTasks_Compactness += 1
                    totalANW_Compactness += no_workers
                    totalATD_Compactness += dist
                    totalCell_Compactness += len(Cells)
		    totalCompactness_Compactness += q_log[-1][3]
                    totalHop_Compactness += no_hops
                    totalHop2_Compactness += no_hops2
                    totalCov_Compactness += coverage
		performed, worker, dist = performed_tasks(workers, Params.MTD, t, True)
		if performed:
		    totalATD_FCFS_Compactness += dist

            # Utility
            ANW_Utility = (totalANW_Utility + 0.0)/totalPerformedTasks_Utility
            ATD_Utility = totalATD_Utility/totalPerformedTasks_Utility
	    ATD_FCFS_Utility = totalATD_FCFS_Utility/totalPerformedTasks_Utility
            ASC_Utility = (totalCell_Utility + 0.0)/totalPerformedTasks_Utility
	    CMP_Utility = totalCompactness_Utility/totalPerformedTasks_Utility
            APPT_Utility = 100*float(totalPerformedTasks_Utility)/Params.TASK_NO
            HOP_Utility = (totalHop_Utility + 0.0)/Params.TASK_NO
            HOP2_Utility = (totalHop2_Utility + 0.0)/Params.TASK_NO
            COV_Utility = 100*(totalCov_Utility + 0.0)/Params.TASK_NO


            # Hybrid
            ANW_Hybrid = (totalANW_Hybrid + 0.0)/totalPerformedTasks_Hybrid
            ATD_Hybrid = totalATD_Hybrid/totalPerformedTasks_Hybrid
	    ATD_FCFS_Hybrid = totalATD_FCFS_Hybrid/totalPerformedTasks_Hybrid
            ASC_Hybrid = (totalCell_Hybrid + 0.0)/totalPerformedTasks_Hybrid
	    CMP_Hybrid = totalCompactness_Hybrid/totalPerformedTasks_Hybrid
            APPT_Hybrid = 100*float(totalPerformedTasks_Hybrid)/Params.TASK_NO
	    HOP_Hybrid = (totalHop_Hybrid + 0.0)/Params.TASK_NO
            HOP2_Hybrid = (totalHop2_Hybrid + 0.0)/Params.TASK_NO
            COV_Hybrid = 100*(totalCov_Hybrid + 0.0)/Params.TASK_NO
            
            # Compactness
            ANW_Compactness = (totalANW_Compactness + 0.0)/totalPerformedTasks_Compactness
            ATD_Compactness = totalATD_Compactness/totalPerformedTasks_Compactness
	    ATD_FCFS_Compactness = totalATD_FCFS_Compactness/totalPerformedTasks_Compactness
            ASC_Compactness = (totalCell_Compactness + 0.0)/totalPerformedTasks_Compactness
	    CMP_Compactness = totalCompactness_Compactness/totalPerformedTasks_Compactness
            APPT_Compactness = 100*float(totalPerformedTasks_Compactness)/Params.TASK_NO            
            HOP_Compactness = (totalHop_Compactness + 0.0)/Params.TASK_NO
            HOP2_Compactness = (totalHop2_Compactness + 0.0)/Params.TASK_NO
            COV_Compactness = 100*(totalCov_Compactness + 0.0)/Params.TASK_NO
             
	    res_cube_anw[i,j,0] = ANW_Utility
            res_cube_atd[i,j,0] = ATD_Utility
            res_cube_atd_fcfs[i,j,0] = ATD_FCFS_Utility
            res_cube_appt[i,j,0] = APPT_Utility
            res_cube_cell[i,j,0] = ASC_Utility
            res_cube_cmp[i,j,0] = CMP_Utility
            res_cube_hop[i,j,0] = HOP_Utility
            res_cube_hop2[i,j,0] = HOP2_Utility
            res_cube_cov[i,j,0] = COV_Utility

	    res_cube_anw[i,j,1] = ANW_Hybrid
            res_cube_atd[i,j,1] = ATD_Hybrid
            res_cube_atd_fcfs[i,j,1] = ATD_FCFS_Hybrid
            res_cube_appt[i,j,1] = APPT_Hybrid
            res_cube_cell[i,j,1] = ASC_Hybrid
            res_cube_cmp[i,j,1] = CMP_Hybrid
            res_cube_hop[i,j,1] = HOP_Hybrid
            res_cube_hop2[i,j,1] = HOP2_Hybrid
            res_cube_cov[i,j,1] = COV_Hybrid
            
            res_cube_anw[i,j,2] = ANW_Compactness
            res_cube_atd[i,j,2] = ATD_Compactness
            res_cube_atd_fcfs[i,j,2] = ATD_FCFS_Compactness
            res_cube_appt[i,j,2] = APPT_Compactness
            res_cube_cell[i,j,2] = ASC_Compactness
            res_cube_cmp[i,j,2] = CMP_Compactness
            res_cube_hop[i,j,2] = HOP_Compactness
            res_cube_hop2[i,j,2] = HOP2_Compactness
            res_cube_cov[i,j,2] = COV_Compactness
            
    res_summary_anw = np.average(res_cube_anw,axis=1)
    np.savetxt(Params.resdir+exp_name + '_anw_' + `Params.TASK_NO`, res_summary_anw, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_atd = np.average(res_cube_atd_fcfs,axis=1)
    np.savetxt(Params.resdir+exp_name + '_atd_fcfs_' + `Params.TASK_NO`, res_summary_atd, fmt='%.4f\t')
    res_summary_appt = np.average(res_cube_appt,axis=1)
    np.savetxt(Params.resdir+exp_name + '_appt_' + `Params.TASK_NO`, res_summary_appt, fmt='%.4f\t')
    res_summary_cell = np.average(res_cube_cell,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cell_' + `Params.TASK_NO`, res_summary_cell, fmt='%.4f\t')
    res_summary_cmp = np.average(res_cube_cmp,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cmp_' + `Params.TASK_NO`, res_summary_cmp, fmt='%.4f\t')
    res_summary_hop = np.average(res_cube_hop,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop_' + `Params.TASK_NO`, res_summary_hop, fmt='%.4f\t')
    res_summary_hop2 = np.average(res_cube_hop2,axis=1)
    np.savetxt(Params.resdir+exp_name + '_hop2_' + `Params.TASK_NO`, res_summary_hop2, fmt='%.4f\t')    
    res_summary_cov = np.average(res_cube_cov,axis=1)
    np.savetxt(Params.resdir+exp_name + '_cov_' + `Params.TASK_NO`, res_summary_cov, fmt='%.4f\t')    
    
def _evalGeocast_Parameter_Fitting(data, all_tasks):
    logging.info("evalGeocast_Parameter_Fitting")
    exp_name = "Geocast_Parameter_Fitting"
    
    alphaList = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    res_cube = np.zeros((len(alphaList), len(eps_list), 3))
        
    for j in range(len(eps_list)):
        Params.CUSTOMIZED_GRANULARITY = True
        Params.PARTIAL_CELL_SELECTION = True
        p = Params(seed_list[j])
        p.Eps = eps_list[j]
        tree = Grid_adaptive(data, p)
        tree.buildIndex()
        
        for i in range(len(alphaList)):
            Params.ALPHA = alphaList[i]
            totalANW_Utility, totalANW_Hybrid, totalANW_Compactness = 0, 0, 0
	    totalCompactness_Utility, totalCompactness_Hybrid, totalCompactness_Compactness = 0, 0, 0
            totalPerformedTasks_Utility, totalPerformedTasks_Hybrid, totalPerformedTasks_Compactness= 0, 0, 0
            
#            tasks = tasks_gen(data, Params.TASK_NO, random.randint(0,1000),Params.x_min,Params.y_min,Params.x_max,Params.y_max)
            for l in range(len(tasks)):
                if (l+1)%Params.LOGGING_STEPS == 0:
                    print ">> " + str(l+1) + " tasks completed"
                t = tasks[l]
                
                # Utility
                Params.COST_FUNCTION = "utility"
		q, q_log = geocast(tree, t, p.Eps)
		no_workers, workers, performed, worker, dist, Cells, no_hops, coverage = post_geocast(t, q, q_log)
                
		if performed:
                    totalPerformedTasks_Utility += 1
                    totalANW_Utility += no_workers
#		    totalCompactness_Utility += q_log[-1][3]

                # Hybrid
                Params.COST_FUNCTION = "hybrid"           
                q, q_log = geocast(tree, t, p.Eps)
		no_workers, workers, performed, worker, dist, Cells, no_hops, coverage = post_geocast(t, q, q_log)
		if performed:
                    totalPerformedTasks_Hybrid += 1
                    totalANW_Hybrid += no_workers
		    totalCompactness_Hybrid += q_log[-1][3]

                # Compactness
                Params.COST_FUNCTION = "compactness"  
                q, q_log = geocast(tree, t, p.Eps)
		no_workers, workers, performed, worker, dist, Cells, no_hops, coverage = post_geocast(t, q, q_log)
		if performed:
                    totalPerformedTasks_Compactness += 1
#                    totalANW_Compactness += no_workers
		    totalCompactness_Compactness += q_log[-1][3]

            # Utility
            ANW_Utility = (totalANW_Utility + 0.0)/totalPerformedTasks_Utility
#	    CMP_Utility = totalCompactness_Utility/totalPerformedTasks_Utility


            # Hybrid
            ANW_Hybrid = totalANW_Hybrid/totalPerformedTasks_Hybrid
	    CMP_Hybrid = totalCompactness_Hybrid/totalPerformedTasks_Hybrid
	    
            
            # Compactness
#            ANW_Compactness = totalANW_Compactness/totalPerformedTasks_Compactness
	    CMP_Compactness = totalCompactness_Compactness/totalPerformedTasks_Compactness            

            
            res_cube[i,j,0] = (ANW_Hybrid - ANW_Utility)/ANW_Utility
            print res_cube[i,j,0]
            res_cube[i,j,1] = (CMP_Compactness - CMP_Hybrid)/CMP_Compactness
            res_cube[i,j,2] = (res_cube[i,j,0] + res_cube[i,j,1])/2
#            print res_cube[i,j,0], res_cube[i,j,1], res_cube[i,j,2]
            
    res_summary = np.sum(res_cube,axis=1)
    np.savetxt(Params.resdir+exp_name + `Params.TASK_NO`, res_summary, fmt='%.4f\t')

def _evalGeocast_Test(data, all_tasks):
    Params.CUSTOMIZED_GRANULARITY = True
    Params.PARTIAL_CELL_SELECTION = True
    Params.CONSTRAINT_INFERENCE = False
    p = Params(1000)
    p.Eps = 1
    
    prev_time = datetime.now()
    tree = Grid_adaptive(data, p)
    tree.buildIndex()
    curr_time = datetime.now()
    print curr_time - prev_time
    
    Params.CONSTRAINT_INFERENCE = True
    prev_time = datetime.now()
    tree = Grid_adaptive(data, p)
    tree.buildIndex()
    curr_time = datetime.now()
    print curr_time - prev_time

    
def evalGeocast_Granularity_Check(data, all_tasks):
    """
    Evaluate APPT of geocast algorithm by varying parameter c2
    """
    logging.info("evalGeocast_Granularity_Check")
    exp_name = "Geocast_Granularity_Check_"
    Params.TASK_NO = 200
    temp_c2_c = Params.c2_c
    c2List = [0.1,0.2,0.4,0.8,1.414,1.6,3.2,6.4,12.8,25.6]
    res_cube = np.zeros((len(eps_list), len(seed_list), len(c2List)))

    for j in range(len(seed_list)):
        for i in range(len(eps_list)):
            for k in range(len(c2List)):
                Params.CUSTOMIZED_GRANULARITY = True
                Params.PARTIAL_CELL_SELECTION = False
                p = Params(seed_list[j])
                p.Eps = eps_list[i]
                Params.c2_c = c2List[k]
                tree = Grid_adaptive(data, p)
                tree.buildIndex()

                totalPerformedTasks_Hybrid = 0
                tasks = all_tasks[j]
                for l in range(len(tasks)):
                    if (l+1)%Params.LOGGING_STEPS == 0:
                        print ">> " + str(l+1) + " tasks completed"
                    t = tasks[l]

                    Params.COST_FUNCTION = "utility"
                    q, q_log = geocast(tree, t, p.Eps)
                    no_workers, workers, Cells, no_hops, coverage, no_hops2 = post_geocast(t, q, q_log)
                    performed, worker, dist = performed_tasks(workers, Params.MTD, t, False)
                    if performed:
                        totalPerformedTasks_Hybrid += 1

                # Hybrid 
                APPT_Hybrid = 100*float(totalPerformedTasks_Hybrid)/Params.TASK_NO
                res_cube[i,j,k] = APPT_Hybrid
            
    res_summary = np.sum(res_cube,axis=1)
    np.savetxt(Params.resdir+exp_name + `Params.TASK_NO`, res_summary, fmt='%.4f\t')
    Params.c2_c = temp_c2_c

    
def read_tasks():
    p = Params(0)
    p.select_dataset()
    data = np.genfromtxt(Params.dataset_task,unpack = True)
    return data

if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG, filename='../log/debug.log')
    logging.info(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "  START")
    worker_data = data_readin()
    task_data = read_tasks()
    all_tasks = tasks_gen(task_data, Params.TASK_NO, Params.x_min,Params.y_min,Params.x_max,Params.y_max)
    
    # Experiment: geocast query
#    evalGeocast_Test(data, all_tasks)

    evalGeocast_GRA_PAR(worker_data, all_tasks)
#    evalGeocast_Compactness(worker_data, all_tasks)
#    evalGeocast_Baseline(worker_data, all_tasks)
#    evalGeocast_MAR(worker_data, all_tasks)
#    evalGeocast_Utility(worker_data, all_tasks)
#    evalGeocast_Granularity_Check(worker_data, all_tasks)

#    _evalGeocast_Test(worker_data, all_tasks)

#    evalGeocast_Parameter_Fitting(data)
#    p1 = Process(target = evalGeocast_GRA_PAR, args=(data))
#    p1.start()
#    p2 = Process(target = evalGeocast_Compactness, args=(data))
#    p2.start()
#    p3 = Process(target = evalMTD, args=(data))
#    p3.start()
#    p4 = Process(target = evalGeocastAlgo, args=(data))
#    p4.start()
#    p1.join()
#    p2.join()
#    p3.join()
#    p4.join()