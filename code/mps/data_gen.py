import struct
import os 
import os.path as osp 
import numpy as np 
from tqdm import tqdm

cur_dir = osp.dirname(osp.realpath(__file__))
data_dir = osp.join(cur_dir, "release")


#print(data_dir)

def parser(num_qubits):
    data_file = open(osp.join(data_dir, str(num_qubits)+"-qubits_critical.bin"), 'rb')

    num_fields = 3 * num_qubits
    num_energies = 2 ** (num_qubits)
    num_coefficients = 2 ** (2 * num_qubits)

    num_doubles = num_fields + num_energies + num_coefficients
    struct_size = (num_fields + num_energies + num_coefficients) * 8 

    format_specificer = "".join(['d' for i in range(num_doubles)])
    format_specificer = "@"+format_specificer

    data = data_file.read()

    data = struct.iter_unpack(format_specificer, data)


    ground_states = []

    fields = [] 
    for values in data:
        start = num_fields + num_energies 

        end = start + (2 ** num_qubits)
        ground_states.append(values[start:end])
        fields.append(values[0:num_fields])

    ground_states = np.array(ground_states)
    fields = np.array(fields)
    print(ground_states.shape)
    print(fields.shape)
    np.savez(str(num_qubits)+"_qubit_crit_data",ground_state = ground_states, fields = fields)

def main():
    sample_sites = [2,4,6,7,10]
    for each in tqdm(sample_sites):
        parser(each)
main()
