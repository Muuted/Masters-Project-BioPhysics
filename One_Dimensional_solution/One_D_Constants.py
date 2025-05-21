import numpy as np
import random

def gamma(i):
    gam = 1
    if i==0:
        return gam/2
    else:
        return gam

def One_D_Constants():

    """------ constants ---------"""
    L = 10e-6 # micrometers  :  Total length of line
    r0 = 0.5e-6 # micrometer  :   radius of hole
    N = 19 + 1 # Number of chain links
    m = 1e-6 # grams  :   Mass of each chain link
    gamma = 1#e-6 # units=   : The drag coefficient
    ds =  L/(N-1) # micrometers  :  Length of each chain
    T = 10 # s  : total time simulated
    dt = 1e-5 # s time step.
    k = 10 #    :  Mean curvature modulus
    kG = 10 #   :  Guassian curvature modulus

    
    """------ variables list ---------"""
    # list of variables
    # we are gonna assume for now that the membrane is just initially flat.
    x_list = np.zeros(shape=(T,N)) # the x-list has increasing value from the start
    z_list = np.zeros(shape=(T,N)) # the z-list is just flat
    psi_list = np.zeros(shape=(T,N)) # all the angles are just flat
    
    for i in range(len(psi_list[0])):
        psi_list[0][i] = random.uniform(-3.14/4 , 3.14/4)

    lambda_list = np.zeros(N)
    nu_list = np.zeros(N)


    args = [
        L,r0,N,ds,T,dt
        ,x_list ,z_list ,psi_list 
        ,lambda_list ,nu_list
        ]


    return args

if __name__ == "__main__":   
    One_D_Constants()
