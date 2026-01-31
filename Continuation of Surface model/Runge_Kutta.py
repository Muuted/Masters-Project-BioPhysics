from Two_D_functions import dpsidt_func,drdt_func,dzdt_func



def drdt_RK4(
        i:int,N:int,dt:float,k:float,c0:float, sigma:float, kG:float, tau:float, ds:float,eta:float
        ,Area:list,psi:list,r:list, z:list
        ,lamb:list , nu:list
        ):
    
    k1 = drdt_func(
        i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
        ,Area=Area,psi=psi,lamb=lamb , nu=nu, z_list=z
        ,raid=r
        )

    k2 = drdt_func(
        i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
        ,Area=Area,psi=psi,lamb=lamb , nu=nu, z_list=z
        ,raid=[j if j!=r[i] else j + dt*k1/2 for j in r ]
        )
    
    k3 = drdt_func(
        i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
        ,Area=Area,psi=psi,lamb=lamb , nu=nu, z_list=z
        ,raid=[j if j!=r[i] else j + dt*k2/2 for j in r ]
        )

    k4 = drdt_func(
        i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
        ,Area=Area,psi=psi,lamb=lamb , nu=nu, z_list=z
        ,raid=[j if j!=r[i] else j + dt*k3 for j in r ]
        )
    
    return (dt/6)*(k1 + k2 + k3 + k4)


def  dzdt_RK4(i,dt,Area:list,radi:list, nu:list):
    
    dzdt_1 = dzdt_func(i=i,Area=Area,radi=radi,nu=nu)

    dzdt_2 = dzdt_func(i=i,Area=Area,radi=radi,nu=nu)



def dPsidt_RK4(
        i,N,k,c0,sigma,kG,tau,dt
        ,Area,radi,z_list
        ,lambs,nus
        ,psi
        ):
    
    dpdt_1 = dpsidt_func(
        i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi=psi 
                )

    dpdt_2 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi= psi + (dt/2)*dpdt_1
    )

    dpdt_3 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi= psi + (dt/2)*dpdt_2
    )

    dpdt_4 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi= psi + dt*dpdt_3
    )

    Runge_kutta = dpdt_1 + 2*dpdt_2 + 2*dpdt_3 + dpdt_4
    return Runge_kutta
