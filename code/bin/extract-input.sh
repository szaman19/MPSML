#!/bin/bash 

module load psxe-2019
source activate qm
export KMP_INIT_AT_FORK=FALSE

cat > tmp.py << EOF
from skopt import load

def print_input_file(**kwargs):
	for key,value in kwargs.items():
		print(str(key) + " = " + str(value))

def res_to_dict(res):
    d = {}
    d.update({'algorithm'                   : res.x[0]})
    d.update({'batch_size'                  : res.x[1]})
    d.update({'cutout_option'               : res.x[2]})
    d.update({'epochs'                      : res.x[3]})
    d.update({'hidden_activation'           : res.x[4]})
    d.update({'hidden_layer_size'           : res.x[5]})
    d.update({'input'                       : res.x[6]})
    d.update({'instances'                   : res.x[7]})
    d.update({'lambda'                      : res.x[8]})
    d.update({'learning_rate'               : res.x[9]})
    d.update({'model'                       : res.x[10]})
    d.update({'num_hidden_layers'           : res.x[11]})
    d.update({'qubits'                      : res.x[12]})
    d.update({'sampling'                    : res.x[13]})
    d.update({'save_model'                  : res.x[14]})
    d.update({'seed'                        : res.x[15]})
    d.update({'target_activation'           : res.x[16]})
    d.update({'validation_history_size'     : res.x[17]})
    d.update({'validation_threshold'        : res.x[18]})
    d.update({'write_average_magnetization' : "true"})
    d.update({'write_coefficients'          : "true"})
    d.update({'write_lyapunov_estimate'     : "true"})
    d.update({'write_schrodinger_error'     : "true"})
    return d


res = load("$1")
print_input_file(**res_to_dict(res))
EOF

python tmp.py 
rm tmp.py
