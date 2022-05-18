import torch
import torch.nn as nn 
import torch.nn.functional as f
import numpy as np 
from torch.utils.data import TensorDataset, DataLoader, random_split
import matplotlib.pyplot as plt
import auto_encoder
import math
import bisect
import pickle
import copy
from mpl_toolkits.mplot3d import Axes3D
import ray
from ray import tune
from ray.tune.schedulers import ASHAScheduler


def seed(num_samples):
    sampler = qmc.LatinHypercube(d=2)
    sample = sampler.random(n=100)
    bx, bz = np.squeeze(np.split(sample,2, axis = 1))
    bx = [int(i*num) for i in bx]
    bz = [int(i*num) for i in bz]
    array = []
    for i in range (100):
        array.append(bx[i]*num+bz[i])
    return array


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
    pred = np.sum(predicted,axis=1)
    tru = np.sum(truth,axis=1)
    diffs = abs(((pred+1) - (tru+1))/(tru+1))
    points = (-np.asarray(diffs)).argsort()[:num]
    return points

def mean_error(predicted, truth):
    pred = np.sum(predicted,axis=1)
    tru = np.sum(truth,axis=1)
    diffs = abs(((pred+1) - (tru+1))/(tru+1))
    return np.mean(diffs) 

def max_error(predicted, truth):
    pred = np.sum(predicted,axis=1)
    tru = np.sum(truth,axis=1)
    diffs = abs((pred+1) - (tru+1))/(tru+1)
    return np.sort(diffs)[np.size(diffs)-1]

def error_data(predicted, truth):
    pred = np.sum(predicted,axis=1)
    tru = np.sum(truth,axis=1)
    diffs = abs(((pred+1) - (tru+1))/(tru+1))
    return np.mean(diffs, axis=1)

def print_points(data1, qubits, error_points):    
    data = np.load(data1)
    Bx = data['fields'].T[qubits]
    Bz = data['fields'].T[2*qubits]
    points = []
    for i in error_points:
        points.append((Bx[i],Bz[i]))
    print(points)
    

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
    ############################################## SETUP ########################################################
    get_dataset = auto_encoder.get_dataset
    data_2 = '2_qubit_crit_data_1.npz'
    data_4 = '4_qubit_crit_data_1.npz'
    data_6 = '6_qubit_crit_data_1.npz'
    data_7 = '7_qubit_crit_data_1.npz'
    data_8 = '8_qubit_crit_data.npz'

    v_data_2 = '2_qubit_val_data.npz'
    v_data_4 = '4_qubit_val_data.npz'
    v_data_7 = '7_qubit_val_data.npz'
    v_data_6 = '6_qubit_val_data.npz'

    data_3 = '3_qubit_test_data.npz'
    data_5 = '5_qubit_test_data.npz'
    data_10 = '10_qubit_test_data.npz'

    t_data_2 = '2_qubit_test_data.npz'
    t_data_4 = '4_qubit_test_data.npz'
    t_data_7 = '7_qubit_test_data.npz'
    t_data_6 = '6_qubit_test_data.npz'
    t_data_8 = '8_qubit_test_data.npz'


    #     val_data_9 = get_dataset(data_9, 9, pool_size)

    validation_n_sizes = [6,8]
    mps_size = 5
    device = torch.device("cpu")
    model = auto_encoder.MPS_autoencoder(mps_size = mps_size).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    train_sizes = [2,4,7]
    val_sizes = [6,8,2,4,7]

    pool_size = 10000

    val_size = 1089

    pts = seed2(100)
    error_p = []
    magnetization_6 = []
    wave_func_6 = []
    final_weights = []

    val_data_2 = get_dataset(v_data_2, 2, val_size)
    val_data_4 = get_dataset(v_data_4, 4, val_size)
    val_data_7 = get_dataset(v_data_7, 7, val_size)
    val_data_6 = get_dataset(v_data_6, 6, 2500)
    val_data_8 = get_dataset(data_8, 8, 2500)

    active_val_6 = get_dataset(data_6, 6, 10000)


    test_data_3 = get_dataset(data_3, 3, pool_size)
    test_data_5 = get_dataset(data_5, 5, 2500)
    test_data_10 = get_dataset(data_10, 10, 100)
    test_data_2 = get_dataset(t_data_2, 2, 1681)
    test_data_4 = get_dataset(t_data_4, 4, 1681)
    test_data_7 = get_dataset(t_data_7, 7, 1681)
    test_data_6 = get_dataset(t_data_6, 6, 1681)
    test_data_8 = get_dataset(t_data_8, 8, 1681)


    val_datasets = [val_data_6, val_data_8, val_data_2, val_data_4, val_data_7] # val_data_9

    val_loaders = [DataLoader(x, batch_size = 10, num_workers=5) for x in val_datasets]


    val_6 = next(iter(val_loaders[0]))
    val_8 = next(iter(val_loaders[1]))
    val_2 = next(iter(val_loaders[2]))
    val_4 = next(iter(val_loaders[3]))
    val_7 = next(iter(val_loaders[4]))

    stop = 500
    runs = 0
    increase_points = 5
    sizes = []
    
    ###################################################Active learning########################################################
    
    while(True):
        runs = runs+1
        increase_points = 4 + len(pts)/100   

        training_data_2 = auto_encoder.get_dataset_active(data_2,2,pool_size, pts)
        training_data_4 = auto_encoder.get_dataset_active(data_4,4,pool_size, pts)
        training_data_7 = auto_encoder.get_dataset_active(data_7,7,pool_size, pts) 


        datasets = [training_data_2,
                    training_data_4,
                    training_data_7]

        training_loaders = [DataLoader(x, batch_size = 32,  shuffle=True, num_workers=5) for x in datasets]

        print("Training Validation ", runs)

        model, tot_err, val_err, t_errs, val_errs = auto_encoder.mps_fit(device, model, optimizer, mps_size, training_loaders, train_sizes, val_loaders, val_sizes) 

        sizes.append(len(pts))    

        print("VALIDATION LOSS ITERATION: ", runs)
        val_data = [(val_6,6),(val_8,8),(val_2,2),(val_4,4),(val_7,7)]
        #val_data = [(val_6,6),(val_8,8),(val_9,9),(val_2,2),(val_4,4),(val_7,7)]
        loss_func = nn.MSELoss()
        count=0
        for data, size in val_data:
            with torch.no_grad():
                fields, wf = data
                gs = model(fields, size)
                loss = loss_func(gs,wf)
                current = loss.item() * (2**size)
                print(size,"\t" ,current)

        print("__________________________________________________")
        print()

        f = open("test_data_1.p", 'wb')
        pickle.dump(val_data, f)

        for N, train_loader in enumerate(training_loaders):            
            temp = 0
            sys_size = train_sizes[N]
            for i, (fields,wf) in enumerate(train_loader):
                fields = fields.to(device)
                gs = model(fields, sys_size)            
                loss = loss_func(gs, wf.to(device))
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                temp += loss.item()
            temp = temp / (len(train_loader)) 

        print(model)    

        data_sizes = [2,3,4,5,6,7,8,10]

        mag_dat = [test_data_2, test_data_3, test_data_4, test_data_5, active_val_6,test_data_7, test_data_8, test_data_10]
        mag_loaders = [DataLoader(x, batch_size = 10000, num_workers=5) for x in mag_dat]

        model.eval()
        with torch.no_grad():
            n_systems = {}
            for j,loader in enumerate((mag_loaders)):
                sys_size = data_sizes[j]
                wave_functions = []
                true_wave = []
                for i, (fields,wf) in enumerate(loader):
                    fields = fields
                    gs = model(fields, sys_size)

                    wave_functions.append(gs)
                    true_wave.append(wf)

                n_systems[sys_size] = (wave_functions,true_wave)


        data_y_6 = n_systems[6][0][0].numpy()        
        data_y_6_t = n_systems[6][1][0].numpy()
        vec_6 = seq_to_magnetization(seq_gen(6),6).reshape((64,1))
        mag_6 = np.squeeze((np.power(data_y_6,2) @ vec_6))
        mag_6_t = np.squeeze((np.power(data_y_6_t,2) @ vec_6))

        magnetization_6.append(mag_6)

        final_weights.append(model.encoder[6].weight)

        error_pts = error_points(data_y_6, data_y_6_t, increase_points)
        print("ERROR POINTS:")
        print_points(data_6, 6, error_pts)

        error_p.append(error_data(data_y_6,data_y_6_t))

        if(len(pts)<stop):
            pts = expand_pool(pts, error_pts, pool_size)
            print("NEW SET: ")
            print(pts)
            print()
            print("_________________________________")
            prev_err = new_err
            prev_diff = new_diff
        else:
            break


##################################################Save Results of training###############################################
    from tempfile import TemporaryFile

    with open('magnetization_6.npy', 'wb') as f:
        np.save(f, magnetization_6)

    with open('pts.npy', 'wb') as f:
        np.save(f, pts)   

    with open('error_p.npy', 'wb') as f:
        np.save(f,error_p)

    with open('final_weights.npy','wb') as f:
        np.save(f, final_weights)
        
    with open('sizes.npy', 'wb') as f:
        np.save(f, sizes)

    torch.save(mag_loaders, 'mag_loaders.pth')    
    torch.save(model.state_dict(), "Active.pt")
        
main()
        
    

    