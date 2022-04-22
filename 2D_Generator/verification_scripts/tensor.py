import numpy as np
import math
import copy
import scipy.sparse as sps
class tensor:

    """
    Tensor class to hold sparse data
    Built on top of a dictionary. The keys are a indices (pair of ints in the case of matrics) and value
    Assuming square matrices, the stride holds the number of columnds in the matrix
    """

    def __init__(self, num_rows, _data):
        self._data = _data
        self.stride = num_rows

    def tensor_prod(self, _tensor):
        """ Creates a tensor using a tensor or Kronecker product of it's self and another tensor. 
            Performs the operation A \otimescross B. Where A is itself. 

           _tensor (tensor): Tensor to do the tensor product. The B in the aforementioned equaiton 
        """
        new_data = {}
        new_stride = self.stride * _tensor.stride
        for key, value in self._data.items():
            for other_key, other_value in _tensor._data.items():
                _i, _j = key
                _k, _l = other_key
                new_i = _i * _tensor.stride + _k
                new_j = _j * _tensor.stride + _l
                new_data[(new_i, new_j)] = value * other_value
        return tensor(new_stride, new_data)

    def numpy(self):
        _tensor = np.zeros((self.stride, self.stride))
        for indices, val in self._data.items():
            i, j = indices
            _tensor[i][j] = val
        return _tensor

    def get_as_scipy_sparse(self):
        _tensor = sps.dok_matrix((self.stride,self.stride))
        for indices, val in self._data.items():
            i, j = indices
            _tensor[i,j] = val
        return _tensor
    
    def __add__(self, other_tensor):
        new_data = copy.deepcopy(self._data)

        for key, value in other_tensor._data.items():
            i,j = key
            if((i,j) in new_data):
                new_data[(i,j)] += value
            else:
                new_data[(i,j)] = value
        return tensor(self.stride, new_data)
    
    def scalar_mul(self, scalar):
        new_data = copy.deepcopy(self._data)
        for key, value in new_data.items():
            i,j = key
            new_data[(i,j)] *= scalar

        return tensor(self.stride, new_data)

    @staticmethod
    def join_tensors(t1, t2):
        _data = t1._data
        for _indices, val in t2._data.items():
            if _indices in _data.keys():
                _data[_indices] += val
            else:
                _data[_indices] = val
        return tensor(t1.stride, _data)





