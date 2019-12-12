#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
from skopt.space import Integer, Real, Categorical
from skopt.utils import use_named_args, point_asdict
from skopt import gp_minimize, dump, load
from os import system
from subprocess import call
from skopt.plots import plot_convergence
import sys

cmd = "phynet.x"
model = sys.argv[1]
qubits = sys.argv[2]
sampling = sys.argv[3]
algorithm = sys.argv[4]

def base_output_dir(**kwargs):
    p = "/home/csingh5/Documents/Guided-Machine-Learning/code/chris/data/output/"
    p += str(kwargs["model"]) + "/"
    p += str(kwargs["qubits"]) + "-qubits/"
    p += str(kwargs["instances"]) + "-instances/"
    p += str(kwargs["sampling"]) + "-phase/"
    return p

def input_file_handle(**kwargs):
    return base_output_dir(**kwargs) + "phynet-" + algorithm + ".in"

def write_input_file(**kwargs):
    with open(input_file_handle(**kwargs), 'w') as file:
        for key,value in kwargs.items():
            line = str(key) + " = " + str(value) + "\n"
            file.write(line)


def final_validation_mse(**kwargs):
    p = base_output_dir(**kwargs) + "mse-" + algorithm + ".dat"
    with open(p, 'r') as file:
        for line in file:
            pass
        return float([i for i in line.split()][-1])

def res_to_dict(res):
    d = {}
    d.update({'algorithm'                   : res.x[0]})
    d.update({'batch_size'                  : res.x[1]})
    d.update({'hidden_activation'           : res.x[2]})
    d.update({'hidden_layer_size'           : res.x[3]})
    d.update({'input'                       : res.x[4]})
    d.update({'instances'                   : res.x[5]})
    d.update({'lambda'                      : res.x[6]})
    d.update({'learning_rate'               : res.x[7]})
    d.update({'model'                       : res.x[8]})
    d.update({'num_hidden_layers'           : res.x[9]})
    d.update({'qubits'                      : res.x[10]})
    d.update({'sampling'                    : res.x[11]})
    d.update({'save_model'                  : "true"})
    d.update({'write_average_magnetization' : "true"})
    return d

space = [Categorical([algorithm],   name='algorithm'),
         Integer(low=10, high=1000, name='batch_size'),
         Categorical(['tanh'],      name='hidden_activation'),
         Integer(low=10, high=1000, name='hidden_layer_size'),
         Categorical(['fields'],    name='input'),
         Categorical(['20k','40k','60k','80k','100k'], name='instances'),
         Real(low=0.0001, high=2,   name='lambda', prior='log-uniform'),
         Real(low=0.0001, high=2,   name='learning_rate', prior='log-uniform'),
         Categorical([model],       name='model'),
         Integer(low=1, high=3,     name='num_hidden_layers'),
         Categorical([qubits],      name='qubits'),
         Categorical([sampling],    name='sampling'),
         Categorical(['true'],      name='save_model'),
         Categorical(['true'],      name='write_average_magnetization')
        ]

@use_named_args(space)
def objective(**kwargs):
    write_input_file(**kwargs)
    state = call([cmd, input_file_handle(**kwargs)])
    if (state == 0):
        return final_validation_mse(**kwargs)
    else:
        print("ERROR IN PHYNET PROBABLY")
        exit(-10)

res = gp_minimize(objective, space, n_calls=100) 

dump_location = base_output_dir(**res_to_dict(res)) + "opt-" + algorithm + ".pkl"
save_location = base_output_dir(**res_to_dict(res)) + "con-" + algorithm + ".png"

dump(res, dump_location, store_objective=False)
write_input_file(**res_to_dict(res))
plot_convergence(res)
plt.savefig(save_location)

call([cmd, input_file_handle(**res_to_dict(res))])
