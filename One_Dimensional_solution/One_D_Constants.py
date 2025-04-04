import numpy as np

def One_D_Constants():

    """------ constants ---------"""
    L = 10e-6 # micrometers  :  Total length of line
    r0 = 0.5e-6 # micrometer  :   radius of hole
    N = 10 + 1 # Number of chain links
    m = 1e-6 # grams  :   Mass of each chain link
    ds =  L/(N-1) # micrometers  :  Length of each chain
    T = 10 # s  : total time simulated
    dt = 1e-5 # s time step.
    k = 10 #    :  Mean curvature modulus
    kG = 10 #   :  Guassian curvature modulus

    
    """------ variables list ---------"""
    # list of variables
    # we are gonna assume for now that the membrane is just initially flat.
    x_list = np.zeros(shape=(T,N+1)) # the x-list has increasing value from the start
    z_list = np.zeros(shape=(T,N+1)) # the z-list is just flat
    psi_list = np.zeros(shape=(T,N+1)) # all the angles are just flat
    m_list = np.zeros(N+1) # list of all the masses
    

    for i in range(0,N+1):
        if i == 0:
            m_list[i] = m/2
        else:
            m_list[i] = m

    args = [
        L,r0,N,ds,T,dt
        ,x_list,z_list,psi_list,m_list
        ]


    return args

if __name__ == "__main__":
    
    One_D_Constants()
