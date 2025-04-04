import numpy as np
from One_D_Constants import One_D_Constants

def xi_pos(i,t,N,L,r0,psi_list,ds):
    cos_sums = 0
    for j in range(i,N):
        cos_sums += np.cos(psi_list[t][j])
    
    xi = L + r0 - ds*cos_sums

    return xi

def zi_pos(i,t,N,L,r0,psi_list,ds):
    sin_sums = 0
    for j in range(i,N):
        sin_sums += np.sin(psi_list[t][j])
    
    zi = - ds*sin_sums

    return zi


def d_xi_pos_dt(i,t,N,L,r0,psi_list,ds,dt):
    sin_sums=0
    for j in range(i,N):        
        sin_sums += np.sin(psi_list[t][j])*(psi_list[t+1][j]-psi_list[t][j])
    
    return ds*sin_sums/dt

def d_zi_pos_dt(i,t,N,L,r0,psi_list,ds,dt):
    cos_sums=0
    for j in range(i,N):        
        cos_sums += -np.cos(psi_list[t][j])*(psi_list[t+1][j]-psi_list[t][j])
    
    return ds*cos_sums/dt

def Kinetic_energy(t,N,L,r0,ds,dt,psi_list,m_list):
    E_kin = 0
    for j in range(N):
        dxdt = d_xi_pos_dt(j,t,N,L,r0,psi_list,ds,dt)
        dzdt = d_zi_pos_dt(j,t,N,L,r0,psi_list,ds,dt)
        
        E_kin += (m_list[j]/2)*(  dxdt**2  +  dzdt**2  )

    return E_kin
    

if __name__ == "__main__":
    L,r0,N,ds,T,dt,x_list,z_list,psi_list,m_list = One_D_Constants()
    t= 0
    E_kin = 10
    E_kin = Kinetic_energy(t,N,L,r0,ds,dt,psi_list,m_list)

    
    