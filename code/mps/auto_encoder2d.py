import torch 
import torch.nn as nn 
import torch.nn.functional as f
import numpy as np 
from torch.utils.data import TensorDataset, DataLoader, random_split
import matplotlib.pyplot as plt 

import pickle
import tensor_net_func

def seq_gen(num_q):
    if num_q == 2:
        return ['0000','0001', '0010','0011', '0100','0101', '0110','0111', '1000','1001', '1010','1011','1100','1101', '1110','1111']
    else:
        temp = []
        
        smaller_vals = seq_gen(num_q-1)
        for i in ['0','1']:
            for each in smaller_vals:
                temp.append(i+each)
        return temp   
    
class MPS_autoencoder2d(nn.Module):
    def __init__(self, mps_size):
        super(MPS_autoencoder2d, self).__init__()        
        self.mps_size = mps_size
        
        mat_size = mps_size ** 4
        self.encoder = nn.Sequential(nn.Linear(4, 1024),
                                     nn.ReLU(),
                                     nn.Dropout(.2),
                                     nn.Linear(1024, 32),
                                     nn.ReLU(),
                                     nn.Dropout(.2),
                                     nn.Linear(32, 2 * (mat_size))
                                     )
    def encode(self, x):
        encoded = self.encoder(x)
        temp = encoded.view(-1, 2,self.mps_size,self.mps_size)
        spin_up, spin_down = torch.split(temp, 1, dim = 1)
        return spin_up, spin_down
     
    def decode(self, spin_up, spin_down, num_qubits):
        return wave_func(spin_up, spin_down, num_qubits, 1000)
    
    def forward(self, x, num_qubits):
        spin_up, spin_down = self.encode(x)
        gs = self.decode(spin_up, spin_down, num_qubits)
        gs = gs / torch.norm(gs, dim = 1).view(-1,1)
        return gs    
    
def check_converged(prev_losses, cur_loss):
    rolling_avg = (sum(prev_losses) / len(prev_losses) )
        
    return 1.05 * rolling_avg < cur_loss

def get_dataset(fname, num_qubits, num_samples):
    data = np.load(fname)
    
    _y = data['ground_state']
    _x = data['fields'][:,[0, num_qubits, 2*num_qubits]]

    _num_qubits_column = num_qubits * (np.ones((num_samples,1)))
    _data_x = np.hstack((_x, _num_qubits_column))
    
    dataset = TensorDataset(torch.Tensor(_data_x),torch.Tensor(_y))

    return dataset 

def print_errors(dic, epoch):
    print_string = ""
    
    template = "{}_qubit_loss: {:.5f} \t \t |"
    for each in dic.keys():
        error = dic[each][epoch] 
        print_string += template.format(each, error)

    print(print_string)

def mps_fit(mps_size,
            train_loaders,
            train_sizes,
            val_loaders,
            val_sizes,
            avg_threshold = 10,
            print_interval = 2):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = MPS_autoencoder2d(mps_size = mps_size).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_func = nn.MSELoss()

    train_errors = {}
    val_errors = {}

    for size in train_sizes:
        train_errors[size] = []

    for size in val_sizes:
        val_errors[size] = []
    
    total_err = []

    val_err = [] 

    counter = 0
    while True:
        counter +=1
        epoch_training_error = 0 
        
        model.train()
        for N, train_loader in enumerate(train_loaders):
            
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
            temp = (temp * (2 ** sys_size)) / (len(train_loader)) 
            train_errors[sys_size].append(temp)
            epoch_training_error += temp 
        
    
        epoch_val_error = 0
        model.eval()
        with torch.no_grad():
            for N, val_loader in enumerate(val_loaders):
                temp = 0
                sys_size = val_sizes[N]

                for i, (fields,wf) in enumerate(val_loader):
                    fields = fields.to(device)
                    gs = model(fields.to(device), sys_size)
                    loss = loss_func(gs, wf.to(device))
                    temp += loss.item()
                    if (sys_size > 8):
                        print(gs[0])
                        print(gs[-1])
                temp = (temp * 2 **(sys_size) )/ (len(val_loader))
                val_errors[sys_size].append(temp)
                epoch_val_error += temp 
        
        total_err.append(epoch_training_error)
        val_err.append(epoch_val_error)

        
        if(counter > avg_threshold):
            if(check_converged(val_err[-avg_threshold:], epoch_val_error)):
                print("Model converged!")

                break

        if (counter % print_interval == 0):
            print("Epoch {} : \t Training_Error: {:.5f} \t Val Error: {:.5f}".format(counter, epoch_training_error, epoch_val_error))
            print('*' *80)
            print_errors(train_errors, counter -1)
            print('*' * 80)
            print_errors(val_errors, counter - 1)
            print('\n\n')
    return model, total_err, val_err, train_errors, val_errors 
                

    