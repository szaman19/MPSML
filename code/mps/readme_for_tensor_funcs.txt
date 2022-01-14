uniform_random: 
gives a stratified random sampling of numbers
Takes n: the number of qubits, and samples: the number of coefficients desired
Used in the net functions 4 qubits and above

indices:
returns a list of all spin configurations for an input
takes num_indices: Number of expected indices (length of number in binary format) 
takes _range: the number to go up until
Used in the net functions for 2 and 3 qubits

instant_indices:
returns the spin configuration (binary format)
takes num_independent_indices: Number of expected indices (length of number in binary format) 
takes i: the index of the desired coefficient to be translated to binary format
Used in all the single net functions

net_NxN:

Returns a tensor of coefficients
For 2x2 and 3x3, returns all coefficients
takes A: the spin up tensor
takes B: the spin down tensor
For 4x4, 5x5, and 6x6, returns a stratified sample of all coefficients
takes samples: the number of desired coefficients

single_net_NxN:
Returns one single coefficient
Takes A the spin up tensor and B the spin down tensor
Takes CI the spin configuration (proper format is given from instant_indices)
Used in Monte_Carlo

batched_single_net_NxN:
Returns a 1 x n batch of a single coefficient
Takes A, a batch of spin up tensors
Takes B, a batch of spin down tensors
Assumes A and B are aligned such that the first spin up tensor in A corresponds to the first spin down tensor in B
Takes CI the spin configuration (proper format is given from instant_indices)

wave_func:
puts the net_NxN functions into one
Takes A, B, and samples to be used the same way as in net_NxN
Also takes net for which net to use; 2, 3, 4, 5, or 6
ignores samples when net is 2 or 3

single_wave_func:
puts the single_net_NxN functions into one
Takes A, B to be used the same way as in single_net_NxN
sample is the same as CI in single_net_NxN
Also takes net for which net to use; 2, 3, 4, 5, or 6

batched_single_wave_func:
puts the batched_single_net_NxN functions into one
Takes A, B to be used the same way as in batched_single_net_NxN
sample is the same as CI in batched_single_net_NxN
Also takes net for which net to use; 2, 3, 4, 5, or 6

Monte_Carlo:
Completes the first 4 steps of Monte_Carlo: Takes a random configuration and goes through each index of it to flip the configuration
Returns a new spin configuration
Takes A and B, the spin up and down tensors
Takes net for which net to use; 2, 3, 4, 5, or 6
