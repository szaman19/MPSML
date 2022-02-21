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
        return ['00','01', '10','11']
    else:
        temp = []
        
        smaller_vals = seq_gen(num_q-1)
        for i in ['0','1']:
            for each in smaller_vals:
                temp.append(i+each)
        return temp   
    
class MPS_autoencoder(nn.Module):
    def __init__(self, mps_size):
        super(MPS_autoencoder, self).__init__()        
        self.mps_size = mps_size
        
        mat_size = mps_size ** 2
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
        # spin_up and spin_down (1,mps.size,mps.size)
        mps = {'0':spin_up, '1':spin_down}
                
        coeffs = [] 
        
        states = seq_gen(num_qubits)
        for state in states:
            mat = mps[state[0]]
            
            for site in state[1:]:
                mat = torch.matmul(mat, mps[site])
            diagonal = torch.diagonal(mat, dim1=-1,dim2=-2)
            coeffs.append(torch.sum(diagonal, dim = -1, keepdim = True))                
        
        c_i = coeffs[0]
        for i in coeffs[1:]:
            c_i = torch.cat((c_i, i), dim = 1)
        return c_i.squeeze()
    
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

def get_dataset_active(fname, num_qubits, num_samples, data_pts):
    data = np.load(fname, allow_pickle=True)
    
    _x = data['fields'][:,[0, num_qubits, 2*num_qubits]]
    _num_qubits_column = num_qubits * (np.ones((num_samples,1)))
    _data_x = np.hstack((_x, _num_qubits_column))
    
    y=[]
    x=[]
    for i in data_pts:
        y.append(data['ground_state'][i])
        x.append(_data_x[i])
    dataset = TensorDataset(torch.Tensor(np.array(x)),torch.Tensor(np.array(y)))
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
            print_interval = 1):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = MPS_autoencoder(mps_size = mps_size).to(device)
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
                


def main():
    data_2 = '2_qubit_crit_data.npz'
    data_4 = '4_qubit_crit_data.npz'
    data_6 = '6_qubit_crit_data.npz'
    data_7 = '7_qubit_crit_data.npz'
    data_10 = '10_qubit_crit_data.npz'
    
    data_11 = '11_qubit_crit_data.npz'
    data_12 = '12_qubit_crit_data.npz'



    training_n_sizes = [2,4,7]
    validation_n_sizes = [6,2,4,7]


    training_data_2 = get_dataset(data_2, 2, 10000)
    training_data_4 = get_dataset(data_4,4,10000)
    training_data_7 = get_dataset(data_7,7,10000) 
    
    training_data_2, test_data_2 = random_split(training_data_2, [9000,1000])
    training_data_4, test_data_4 = random_split(training_data_4, [9000,1000])
    training_data_7, test_data_7 = random_split(training_data_7, [9000,1000])
    
    training_data_2, val_data_2 = random_split(training_data_2, [8000,1000])
    training_data_4, val_data_4 = random_split(training_data_4, [8000,1000])
    training_data_7, val_data_7 = random_split(training_data_7, [8000,1000])
    training_data_10 = get_dataset(data_10,10,10000)


    datasets = [training_data_2,
                training_data_4,
                training_data_7]

    training_loaders = [DataLoader(x, batch_size = 32,  shuffle=True, num_workers=20) for x in datasets]

    val_data_6 = get_dataset(data_6, 6, 10000)
    val_data_6, test_data_6 = random_split(val_data_6, [9000,1000])
    #val_data_11 = get_dataset(data_11,11,100)
    #val_data_12 = get_dataset(data_12,12,100)

    val_datasets = [val_data_6, val_data_2, val_data_4, val_data_7]
    val_loaders = [DataLoader(x, batch_size=10000, num_workers=20) for x in val_datasets]

    test_loader_10 = DataLoader(training_data_10, batch_size = 10000, num_workers = 20)
    test_datasets = [test_data_2, test_data_4, test_data_6,test_data_7]
    test_loaders = [DataLoader(x, batch_size=1000, num_workers=20) for x in test_datasets]
    test_loaders.append(test_loader_10)
    
    for mps_size in range(2,11):
        err_name = "{}_dump_errors.p".format(mps_size)
        model_name = "{}_site_model.pt".format(mps_size)

        model, tot_err, val_err, t_errs, val_errs = mps_fit(mps_size,
                                                            training_loaders,
                                                            training_n_sizes,
                                                            val_loaders,
                                                            validation_n_sizes)
        loss_func = nn.MSELoss()
        model.eval()
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        test_sizes = [2, 4, 6, 7, 10]
        test_errors = {}
        with torch.no_grad():
            for j, loader in enumerate(test_loaders):
                sys_size = test_sizes[j]
                temp = 0
                for i, (fields, wf) in enumerate(loader):
                    fields = fields.to(device)
                    gs = model(fields.to(device), sys_size)
                    loss = loss_func(gs, wf.to(device))
                    print(loss.item() * (2**sys_size))
                    print(gs[0])
                    print(gs[-1])
                    temp +=loss.item()

                test_errors[sys_size] = (loss/len(loader))
        test_error_file = open("{}_test_err_.p".format(mps_size),'wb')
        pickle.dump(test_errors, test_error_file)
        test_error_file.close()

        errs = [tot_err, val_err, t_errs, val_errs]
        torch.save(model.state_dict(), model_name)
        f = open(err_name, 'wb')
        pickle.dump(errs,f)
        f.close()

    '''
    errors = []
    for i in range(2,7):
        print("MPS SIZE: {}x{}".format(i,i))
        print("*" * 110)
        errors.append(tune_mps(i, 100)[1])
        print("*" * 110)
        
    print("*"*110)
    
    f = open('logging_errors.p', 'wb')
    pickle.dump(errors, f)
    f.close()
    
    for mps_size in range(11):
        val_error = errors[mps_size][-1]
        min_err , idx = min((val, idx) for (idx, val) in enumerate(val_error))
        print("Min Validation {}x{}, error: {:.5f} \t epoch: {} ".format(mps_size + 2,mps_size + 2, min_err, idx ) )

    '''
if __name__ == '__main__':
    main()

