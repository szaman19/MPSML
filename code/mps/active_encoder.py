#!python -m pip install torch
import torch
import torch.nn as nn 
import torch.nn.functional as f
import numpy as np 
from torch.utils.data import TensorDataset, DataLoader, random_split
import matplotlib.pyplot as plt
#!python -m pip install tensorflow
import auto_encoder
import math
import bisect
import pickle
import copy
from mpl_toolkits.mplot3d import Axes3D


def seed(num_samples):
    pts = []
    i = 0
    it1 = int(math.sqrt(num_samples)/10)
    it2 = int((num_samples/10) - math.sqrt(num_samples))
    for j in range (10):
        for k in range(10):
            pts.append(i)
            i+=it1
        i+=it2
    return pts


def seq_to_magnetization(arr_seq, num_qubits):
    
    mag_vec = []
    for elem in arr_seq:
        magnetization = 0
        for char in elem:
            temp = (int(char)*-2)+1
            magnetization += temp 
        mag_vec.append(magnetization)
    mag_vec = np.array(mag_vec)
    mag_vec = mag_vec / (num_qubits)
    return mag_vec   


def seq_gen(num_q):
    if num_q == 2:
        return ['00','01', '10','11']
    else:
        temp = []
        
        smaller_vals = seq_gen(num_q-1)
        for i in ['0','1']:
            for each in smaller_vals:
                temp.append(i+each)
        return temp 
    

def error_points(predicted, truth, num):
    diffs = abs((predicted) - (truth))/(truth)
    points = (-np.asarray(diffs)).argsort()[:num]
    return points


def expand_pool(old_pool, error_points, max_pool_size):
    for i in error_points:
        num = int(i)
        if(i not in old_pool):
            bisect.insort(old_pool,i)
        else:
            Bz_low_num = num-1
            Bz_high_num= num+1
            Bz_low_flag=True
            Bz_high_flag=True
            
            Bx_low_num = int(num-math.sqrt(max_pool_size))
            Bx_high_num= int(num+math.sqrt(max_pool_size))
            Bx_low_flag=True
            Bx_high_flag=True
            
            while(Bz_low_flag or Bz_high_flag or Bx_low_flag or Bx_high_flag):
                if(Bz_low_flag):
                    if(Bz_low_num<=0 or Bz_low_num % math.sqrt(max_pool_size)==math.sqrt(max_pool_size)-1):
                        Bz_low_flag=False
                    elif(Bz_low_num not in old_pool):
                        bisect.insort(old_pool,Bz_low_num)
                        Bz_low_flag=False
                    else:
                        Bz_low_num-=1
                if(Bz_high_flag):
                    if(Bz_high_num>=max_pool_size or Bz_high_num % math.sqrt(max_pool_size)==0):
                        Bz_high_flag=False
                    elif(Bz_high_num not in old_pool):
                        bisect.insort(old_pool, Bz_high_num)
                        Bz_high_flag=False
                    else:
                        Bz_high_num+=1
                if(Bx_low_flag):
                    if(Bx_low_num<=0):
                        Bx_low_flag=False
                    elif(Bx_low_num not in old_pool):
                        bisect.insort(old_pool, Bx_low_num)
                        Bx_low_flag=False
                    else:
                        Bx_low_num-=int(math.sqrt(max_pool_size))
                if(Bx_high_flag):
                    if(Bx_high_num>=max_pool_size or Bx_high_num % math.sqrt(max_pool_size)==0):
                        Bx_high_flag=False
                    elif(Bx_high_num not in old_pool):
                        bisect.insort(old_pool, Bx_high_num)
                        Bx_high_flag=False
                    else:
                        Bx_high_num+=int(math.sqrt(max_pool_size))    
                    
    return old_pool



def main():
    get_dataset = auto_encoder.get_dataset
    data_2 = '2_qubit_crit_data.npz'
    data_4 = '4_qubit_crit_data.npz'
    data_6 = '6_qubit_crit_data.npz'
    #data_7 = '7_qubit_crit_data.npz'
    data_8 = '8_qubit_crit_data.npz'
    # data_9 = '9_qubit_crit_data.npz'

    training_n_sizes = [2,4,8]
    validation_n_sizes = [6]

    pool_size = 2500

    pts = seed(pool_size)

    stop = 1

    error_p = []

    magnetization_6 = []
    wave_func_6 = []
    
    
    for learning_set in range (stop):

        training_data_2 = auto_encoder.get_dataset_active(data_2,2,pool_size, pts)
        training_data_4 = auto_encoder.get_dataset_active(data_4,4,pool_size, pts)
        training_data_8 = auto_encoder.get_dataset_active(data_8,8,pool_size, pts) 


        split1 = int(len(pts)*0.9)
        split2 = int(len(pts)/10)
        split1 += len(pts)-(split1+split2)

        training_data_2, val_data_2 = random_split(training_data_2, [split1,split2])
        training_data_4, val_data_4 = random_split(training_data_4, [split1,split2])
        training_data_8, val_data_8 = random_split(training_data_8, [split1,split2])

        datasets = [training_data_2,
                    training_data_4,
                    training_data_8]

        training_loaders = [DataLoader(x, batch_size = 32,  shuffle=True, num_workers=20) for x in datasets]

        val_data_6 = get_dataset(data_6, 6, pool_size)
    #     val_data_8 = get_dataset(data_8, 8, pool_size)
    #     val_data_9 = get_dataset(data_9, 9, pool_size)

        val_datasets = [val_data_6, val_data_2, val_data_4, val_data_8] # val_data_8, val_data_9

        val_loaders = [DataLoader(x, batch_size = 10, num_workers=20) for x in val_datasets]

        warmup_2 = next(iter(training_loaders[0]))
        warmup_4 = next(iter(training_loaders[1]))
        warmup_8 = next(iter(training_loaders[2]))

        mps_size = 5
        model = auto_encoder.MPS_autoencoder(mps_size = mps_size)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
        loss_func = nn.MSELoss(reduction='sum')

        warmup_data = [(warmup_2,2), (warmup_4,4), (warmup_8,8)]

        print("WARMUP TRAINING ITERATION: ", learning_set)
        for j in range(10):
            for i in range(3):
                for epoch in range(10):
                    fields,wf = warmup_data[i][0]
                    gs = model(fields, warmup_data[i][1])            
                    loss = loss_func(gs, wf)
                    if (epoch % 10 == 0):
                        current_loss = loss.item() *(2**warmup_data[i][1])
                        print(warmup_data[i][1],"\t", current_loss)
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                if(j==9):
                    final_warmup_loss[learning_set][i] = current_loss
        print("__________________________________________________")
        print()

        print("Training Validation ", learning_set)
        train_sizes, train_loader = enumerate(warmup_data)
        val_6 = next(iter(val_loaders[0]))
        #val_8 = next(iter(val_loaders[1]))
        #val_9 = next(iter(val_loaders[2]))
        val_2 = next(iter(val_loaders[1]))
        val_4 = next(iter(val_loaders[2]))
        val_8 = next(iter(val_loaders[3]))

        val_data = [(val_6,6),(val_2,2),(val_4,4),(val_8,8)]
        val_loaders, val_sizes = val_data

        model, tot_err, val_err, t_errs, val_errs = auto_encoder.mps_fit(mps_size, train_loaders, train_sizes, val_loaders, val_sizes)


main()