import tensorflow as tf
import numpy as np 
import torch
import random
import math
import copy

#For sampling from all coefficients
def uniform_random(n, samples):
    indis = [0 for i in range(samples)]
    for i in range (samples):
        indis[i] = random.randint(int(i*2**(n*n)/samples), int((i+1)*2**(n*n)/samples))
    return indis

#For generating all coefficients
def indices ( num_indices , _range ):
    _list = [[]]
    return ( _index_gen_helper ( num_indices , _range , _list ))
def _index_gen_helper ( num_indices , _range , _list ):
    _ret = []

    for i in _list :
        for j in range ( _range ):
            _temp = i. copy ()
            _temp . extend ([ j ])
            _ret . append ( _temp )
    if num_indices == 1:
        return _ret
    else :
        return _index_gen_helper ( num_indices -1 , _range , _ret )
    
#Generate the spin configuration given the index of the coefficient    
def instant_indices(num_independent_indices, i):
    inds=[0 for i in range (num_independent_indices)]
    for j in range (0, num_independent_indices):
        inds[j] = i%2
        i//=2
    return inds


##################################################### 2 X 2 ##########################################################################

#returns all the coefficients for a 2 X 2 network
def net_2x2(A,B):
    AB=[A, B]
    C = [0 for i in range (2**4)]
    count=0
    for i ,j ,k ,l in indices(4,2):  
        C[count]=torch.einsum('zijab, zjicd, zklba, zlkdc->z', AB[i], AB[j], AB[k], AB[l])
        count+=1  
    ans=torch.stack(C)
    return tf.transpose(ans)

#returns the coefficients for a batch of 2 X 2 network given a spin configuration
def batched_single_net_2x2(A,B, i):
    A.requires_grad = True
    B.requires_grad = True
    AB=[A, B]
    ans=torch.einsum('zijab, zjicd, zklba, zlkdc->z', AB[i[0]], AB[i[1]], AB[i[2]], AB[i[3]])
    return ans.reshape(A.shape[0],1)

#returns the coefficients for a single 2 X 2 network given a spin configuration
def single_net_2x2(A,B, i):
    AB=[A, B]
    ans=torch.einsum('ijab, jicd, klba, lkdc', AB[i[0]], AB[i[1]], AB[i[2]], AB[i[3]])
    return ans


##################################################### 3 X 3 ##########################################################################

def net_3x3(A,B):
    AB=[A,B]
    C = [0 for i in range (2**9)]
    count=0
    for i ,j ,k ,l ,m ,n ,o ,p ,q in indices(9, 2):    
        C[count]=torch.einsum('zabjk, zcamn, zbcpq, zdelj, zfdom, zefrp, zghkl, zigno, zhiqr->z', AB[i], AB[j], AB[k], AB[l], AB[m], AB[n], AB[o], AB[p], AB[q])
        count+=1
    ans=torch.stack(C)
    return tf.transpose(ans)

def batched_single_net_3x3(A,B, i):
    A.requires_grad = True
    B.requires_grad = True
    AB=[A,B]  
    ans=torch.einsum('zabjk, zcamn, zbcpq, zdelj, zfdom, zefrp, zghkl, zigno, zhiqr->z', AB[i[0]], AB[i[1]], AB[i[2]], AB[i[3]], AB[i[4]], AB[i[5]], AB[i[6]], AB[i[7]], AB[i[8]])
    return ans.reshape(A.shape[0],1)

def single_net_3x3(A,B, i):
    AB=[A,B]  
    ans=torch.einsum('abjk, camn, bcpq, delj, fdom, efrp, ghkl, igno, hiqr', AB[i[0]], AB[i[1]], AB[i[2]], AB[i[3]], AB[i[4]], AB[i[5]], AB[i[6]], AB[i[7]], AB[i[8]])
    return ans



##################################################### 4 X 4 ##########################################################################

def net_4x4(A,B, samples):
    AB=[A,B]
    C = [0 for i in range (samples)]
    count=0
    indis = uniform_random(4, samples)
    for i in (indis):
        CI = instant_indices(16,i) #current index
        X = torch.einsum('zabeh, zcail, zdcmp, zbdqt-> zehilmpqt', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]])
        Y = torch.einsum('zabfe, zcaji, zdcnm, zbdrq-> zfejinmrq', AB[CI[4]], AB[CI[5]], AB[CI[6]], AB[CI[7]])
        Z = torch.einsum('zabgf, zcakj, zdcon, zbdsr-> zgfkjonsr', AB[CI[8]], AB[CI[9]], AB[CI[10]], AB[CI[11]])
        W = torch.einsum('zabhg, zcalk, zdcpo, zbdts-> zhglkpots', AB[CI[12]], AB[CI[13]], AB[CI[14]], AB[CI[15]])
        C[count] = torch.einsum('zehilmpqt, zfejinmrq, zgfkjonsr, zhglkpots->z', X, Y, Z, W)
        count+=1
    ans=torch.stack(C)
    return tf.transpose(ans)

def batched_single_net_4x4(A,B, CI):
    A.requires_grad = True
    B.requires_grad = True
    AB=[A,B]
    X = torch.einsum('zabeh, zcail, zdcmp, zbdqt-> zehilmpqt', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]])
    Y = torch.einsum('zabfe, zcaji, zdcnm, zbdrq-> zfejinmrq', AB[CI[4]], AB[CI[5]], AB[CI[6]], AB[CI[7]])
    Z = torch.einsum('zabgf, zcakj, zdcon, zbdsr-> zgfkjonsr', AB[CI[8]], AB[CI[9]], AB[CI[10]], AB[CI[11]])
    W = torch.einsum('zabhg, zcalk, zdcpo, zbdts-> zhglkpots', AB[CI[12]], AB[CI[13]], AB[CI[14]], AB[CI[15]])
    ans = torch.einsum('zehilmpqt, zfejinmrq, zgfkjonsr, zhglkpots->z', X, Y, Z, W)
    return ans.reshape(A.shape[0],1)

def single_net_4x4(A,B, CI):
    AB=[A,B]
    X = torch.einsum('abeh, cail, dcmp, bdqt-> ehilmpqt', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]])
    Y = torch.einsum('abfe, caji, dcnm, bdrq-> fejinmrq', AB[CI[4]], AB[CI[5]], AB[CI[6]], AB[CI[7]])
    Z = torch.einsum('abgf, cakj, dcon, bdsr-> gfkjonsr', AB[CI[8]], AB[CI[9]], AB[CI[10]], AB[CI[11]])
    W = torch.einsum('abhg, calk, dcpo, bdts-> hglkpots', AB[CI[12]], AB[CI[13]], AB[CI[14]], AB[CI[15]])
    ans = torch.einsum('ehilmpqt, fejinmrq, gfkjonsr, hglkpots', X, Y, Z, W)
    return ans


##################################################### 5 X 5 ##########################################################################

def net_5x5(A,B, samples):
    AB=[A,B]
    C = [0 for i in range (samples)]
    count=0
    indis = uniform_random(5, samples)
    for i in (indis):
        CI = instant_indices(25,i) #current index
        X = torch.einsum('zabfg, zcahi, zdcjk, zedlm, zbeno -> zfghijklmno', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]], AB[CI[4]]) 
        Y = torch.einsum('zabpf, zcaqh, zdcrj, zedsl, zbetn -> zpfqhrjsltn', AB[CI[5]], AB[CI[6]], AB[CI[7]], AB[CI[8]], AB[CI[9]])
        XY = torch.einsum('zfghijklmno, zpfqhrjsltn -> zgpiqkrmson', X, Y)
        Z = torch.einsum('zabfp, zcahq, zdcjr, zedls, zbent -> zfphqjrlsnt' , AB[CI[10]], AB[CI[11]], AB[CI[12]], AB[CI[13]], AB[CI[14]])
        XYZ = torch.einsum('zgpiqkrmson, zfphqjrlsnt -> zgfihkjmlon', XY, Z)
        U = torch.einsum('zabpf, zcaqh, zdcrj, zedsl, zbetn -> zpfqhrjsltn', AB[CI[15]], AB[CI[16]], AB[CI[17]], AB[CI[18]], AB[CI[19]])
        XYZU = torch.einsum('zgfihkjmlon, zpfqhrjsltn -> zgpiqkrmsot', XYZ, U)
        V = torch.einsum('zabgp, zcaiq, zdckr, zedms, zbeot -> zgpiqkrmsot', AB[CI[20]], AB[CI[21]], AB[CI[22]], AB[CI[23]], AB[CI[24]])
        C[count] = torch.einsum('zgpiqkrmsot, zgpiqkrmsot->z', XYZU, V)
        count+=1
    ans=torch.stack(C)
    return tf.transpose(ans)

def batched_single_net_5x5(A,B, CI):
    A.requires_grad = True
    B.requires_grad = True
    AB=[A,B]
    X = torch.einsum('zabfg, zcahi, zdcjk, zedlm, zbeno -> zfghijklmno', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]], AB[CI[4]]) 
    Y = torch.einsum('zabpf, zcaqh, zdcrj, zedsl, zbetn -> zpfqhrjsltn', AB[CI[5]], AB[CI[6]], AB[CI[7]], AB[CI[8]], AB[CI[9]])
    XY = torch.einsum('zfghijklmno, zpfqhrjsltn -> zgpiqkrmson', X, Y)
    Z = torch.einsum('zabfp, zcahq, zdcjr, zedls, zbent -> zfphqjrlsnt' , AB[CI[10]], AB[CI[11]], AB[CI[12]], AB[CI[13]], AB[CI[14]])
    XYZ = torch.einsum('zgpiqkrmson, zfphqjrlsnt -> zgfihkjmlon', XY, Z)
    U = torch.einsum('zabpf, zcaqh, zdcrj, zedsl, zbetn -> zpfqhrjsltn', AB[CI[15]], AB[CI[16]], AB[CI[17]], AB[CI[18]], AB[CI[19]])
    XYZU = torch.einsum('zgfihkjmlon, zpfqhrjsltn -> zgpiqkrmsot', XYZ, U)
    V = torch.einsum('zabgp, zcaiq, zdckr, zedms, zbeot -> zgpiqkrmsot', AB[CI[20]], AB[CI[21]], AB[CI[22]], AB[CI[23]], AB[CI[24]])
    ans = torch.einsum('zgpiqkrmsot, zgpiqkrmsot->z', XYZU, V)
    return ans.reshape(A.shape[0],1)

def single_net_5x5(A,B, CI):
#     A.requires_grad = True
#     B.requires_grad = True
    AB=[A,B]
    X = torch.einsum('abfg, cahi, dcjk, edlm, beno -> fghijklmno', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]], AB[CI[4]]) 
    Y = torch.einsum('abpf, caqh, dcrj, edsl, betn -> pfqhrjsltn', AB[CI[5]], AB[CI[6]], AB[CI[7]], AB[CI[8]], AB[CI[9]])
    XY = torch.einsum('fghijklmno, pfqhrjsltn -> gpiqkrmson', X, Y)
    Z = torch.einsum('abfp, cahq, dcjr, edls, bent -> fphqjrlsnt' , AB[CI[10]], AB[CI[11]], AB[CI[12]], AB[CI[13]], AB[CI[14]])
    XYZ = torch.einsum('gpiqkrmson, fphqjrlsnt -> gfihkjmlon', XY, Z)
    U = torch.einsum('abpf, caqh, dcrj, edsl, betn -> pfqhrjsltn', AB[CI[15]], AB[CI[16]], AB[CI[17]], AB[CI[18]], AB[CI[19]])
    XYZU = torch.einsum('gfihkjmlon, pfqhrjsltn -> gpiqkrmsot', XYZ, U)
    V = torch.einsum('abgp, caiq, dckr, edms, beot -> gpiqkrmsot', AB[CI[20]], AB[CI[21]], AB[CI[22]], AB[CI[23]], AB[CI[24]])
    ans = torch.einsum('gpiqkrmsot, gpiqkrmsot', XYZU, V)
    return ans


##################################################### 6 X 6 ##########################################################################

def net_6x6(A, B, samples):
    AB=[A,B]
    C = [0 for i in range (samples)]
    count=0
    indis = uniform_random(6, samples)
    for i in (indis):
        CI = instant_indices(36,i) #current index
        X = torch.einsum('zabnh, zcaoi, zdcpj, zedqk, zferl, zbfsm -> znhoipjqkrlsm', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]], AB[CI[4]], AB[CI[5]]) 
        Y = torch.einsum('zabtn, zcauo, zdcvp, zedwq, zfexr, zbfys -> ztnuovpwqxrys', AB[CI[6]], AB[CI[7]], AB[CI[8]], AB[CI[9]], AB[CI[10]], AB[CI[11]])
        XY = torch.einsum('znhoipjqkrlsm, ztnuovpwqxrys -> zhtiujvkwlxmy', X, Y)
        Z = torch.einsum('zabnt, zcaou, zdcpv, zedqw, zferx, zbfsy -> zntoupvqwrxsy', AB[CI[12]], AB[CI[13]], AB[CI[14]], AB[CI[15]], AB[CI[16]], AB[CI[17]])
        XYZ = torch.einsum('zhtiujvkwlxmy, zntoupvqwrxsy -> zhniojpkqlrms', XY, Z)
        U = torch.einsum('zabtn, zcauo, zdcvp, zedwq, zfexr, zbfys -> ztnuovpwqxrys', AB[CI[18]], AB[CI[19]], AB[CI[20]], AB[CI[21]], AB[CI[22]], AB[CI[23]])
        XYZU = torch.einsum('zhniojpkqlrms, ztnuovpwqxrys -> zhtiujvkwlxmy', XYZ, U)
        V = torch.einsum('zabnt, zcaou, zdcpv, zedqw, zferx, zbfsy -> zntoupvqwrxsy', AB[CI[24]], AB[CI[25]], AB[CI[26]], AB[CI[27]], AB[CI[28]], AB[CI[29]])
        XYZUV = torch.einsum('zhtiujvkwlxmy, zntoupvqwrxsy -> zhniojpkqlrms', XYZU, V)
        W = torch.einsum('zabhn, zcaio, zdcjp, zedkq, zfelr, zbfms -> zhniojpkqlrms', AB[CI[30]], AB[CI[31]], AB[CI[32]], AB[CI[33]], AB[CI[34]], AB[CI[35]])
        C[count] = torch.einsum('zhniojpkqlrms, zhniojpkqlrms->z', XYZUV, W)
        count+=1
    ans=torch.stack(C)
    return tf.transpose(ans)

def batched_single_net_6x6(A, B, CI):
    A.requires_grad = True
    B.requires_grad = True
    AB=[A,B] 
    X = torch.einsum('zabnh, zcaoi, zdcpj, zedqk, zferl, zbfsm -> znhoipjqkrlsm', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]], AB[CI[4]], AB[CI[5]]) 
    Y = torch.einsum('zabtn, zcauo, zdcvp, zedwq, zfexr, zbfys -> ztnuovpwqxrys', AB[CI[6]], AB[CI[7]], AB[CI[8]], AB[CI[9]], AB[CI[10]], AB[CI[11]])
    XY = torch.einsum('znhoipjqkrlsm, ztnuovpwqxrys -> zhtiujvkwlxmy', X, Y)
    Z = torch.einsum('zabnt, zcaou, zdcpv, zedqw, zferx, zbfsy -> zntoupvqwrxsy', AB[CI[12]], AB[CI[13]], AB[CI[14]], AB[CI[15]], AB[CI[16]], AB[CI[17]])
    XYZ = torch.einsum('zhtiujvkwlxmy, zntoupvqwrxsy -> zhniojpkqlrms', XY, Z)
    U = torch.einsum('zabtn, zcauo, zdcvp, zedwq, zfexr, zbfys -> ztnuovpwqxrys', AB[CI[18]], AB[CI[19]], AB[CI[20]], AB[CI[21]], AB[CI[22]], AB[CI[23]])
    XYZU = torch.einsum('zhniojpkqlrms, ztnuovpwqxrys -> zhtiujvkwlxmy', XYZ, U)
    V = torch.einsum('zabnt, zcaou, zdcpv, zedqw, zferx, zbfsy -> zntoupvqwrxsy', AB[CI[24]], AB[CI[25]], AB[CI[26]], AB[CI[27]], AB[CI[28]], AB[CI[29]])
    XYZUV = torch.einsum('zhtiujvkwlxmy, zntoupvqwrxsy -> zhniojpkqlrms', XYZU, V)
    W = torch.einsum('zabhn, zcaio, zdcjp, zedkq, zfelr, zbfms -> zhniojpkqlrms', AB[CI[30]], AB[CI[31]], AB[CI[32]], AB[CI[33]], AB[CI[34]], AB[CI[35]])
    ans = torch.einsum('zhniojpkqlrms, zhniojpkqlrms->z', XYZUV, W)
    return ans.reshape(A.shape[0],1)

def single_net_6x6(A, B, CI):
#     A.requires_grad = True
#     B.requires_grad = True
    AB=[A,B] 
    X = torch.einsum('abnh, caoi, dcpj, edqk, ferl, bfsm -> nhoipjqkrlsm', AB[CI[0]], AB[CI[1]], AB[CI[2]], AB[CI[3]], AB[CI[4]], AB[CI[5]]) 
    Y = torch.einsum('abtn, cauo, dcvp, edwq, fexr, bfys -> tnuovpwqxrys', AB[CI[6]], AB[CI[7]], AB[CI[8]], AB[CI[9]], AB[CI[10]], AB[CI[11]])
    XY = torch.einsum('nhoipjqkrlsm, tnuovpwqxrys -> htiujvkwlxmy', X, Y)
    Z = torch.einsum('abnt, caou, dcpv, edqw, ferx, bfsy -> ntoupvqwrxsy', AB[CI[12]], AB[CI[13]], AB[CI[14]], AB[CI[15]], AB[CI[16]], AB[CI[17]])
    XYZ = torch.einsum('htiujvkwlxmy, ntoupvqwrxsy -> hniojpkqlrms', XY, Z)
    U = torch.einsum('abtn, cauo, dcvp, edwq, fexr, bfys -> tnuovpwqxrys', AB[CI[18]], AB[CI[19]], AB[CI[20]], AB[CI[21]], AB[CI[22]], AB[CI[23]])
    XYZU = torch.einsum('hniojpkqlrms, tnuovpwqxrys -> htiujvkwlxmy', XYZ, U)
    V = torch.einsum('abnt, caou, dcpv, edqw, ferx, bfsy -> ntoupvqwrxsy', AB[CI[24]], AB[CI[25]], AB[CI[26]], AB[CI[27]], AB[CI[28]], AB[CI[29]])
    XYZUV = torch.einsum('htiujvkwlxmy, ntoupvqwrxsy -> hniojpkqlrms', XYZU, V)
    W = torch.einsum('abhn, caio, dcjp, edkq, felr, bfms -> hniojpkqlrms', AB[CI[30]], AB[CI[31]], AB[CI[32]], AB[CI[33]], AB[CI[34]], AB[CI[35]])
    ans = torch.einsum('hniojpkqlrms, hniojpkqlrms', XYZUV, W)
    return ans


##################################################### Single function for 2 through 6 net ####################################################################

def wave_func (A, B, net, samples):
    if(net==2):
        return net_2x2(A, B)
    elif(net==3):
        return net_3x3(A, B)
    elif(net==4):
        return net_4x4(A, B, samples)
    elif(net==5):
        return net_5x5(A, B, samples)
    elif(net==6):
        return net_6x6(A, B, samples)
    else:
        return 0
        
def batched_single_wave_func (A, B, net, sample):
    if(net==2):
        return batched_single_net_2x2(A, B, sample)
    elif(net==3):
        return batched_single_net_3x3(A, B, sample)
    elif(net==4):
        return batched_single_net_4x4(A, B, sample)
    elif(net==5):
        return batched_single_net_5x5(A, B, sample)
    elif(net==6):
        return batched_single_net_6x6(A, B, sample)
    else:
        return 0
        
def single_wave_func (A, B, net, sample):
    if(net==2):
        return single_net_2x2(A, B, sample)
    elif(net==3):
        return single_net_3x3(A, B, sample)
    elif(net==4):
        return single_net_4x4(A, B, sample)
    elif(net==5):
        return single_net_5x5(A, B, sample)
    elif(net==6):
        return single_net_6x6(A, B, sample)
    else:
        return 0
    
##################################################### Monte Carlo ####################################################################    
    
def Monte_Carlo(A, B, net):
    num = random.randint(0, 2**(net**2))
    config = instant_indices(net**2, num)
    print(config)
    Ca = single_wave_func (A, B, net, config)
    unnorma = Ca**2
    
    for i in range (net**2):
        temp_config = copy.copy(config)
        temp_config[i] = (1-temp_config[i])
        Cb = single_wave_func(A, B, net, temp_config)
        unnormb = Cb**2
        r = unnormb/unnorma
        x = random.uniform(0,1)
        
        if (r > x):

            config = copy.copy(temp_config)
            Ca = Cb
        else:
 
    return config