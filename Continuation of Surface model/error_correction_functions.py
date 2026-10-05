import numpy as np
import matplotlib.pyplot as plt
from Two_D_functions import constraint_f, constraint_g, check_constraints_truth


def c_diff_f(
        i:int,j:int,N:int
        ,r:list[float],psi:list[float],Area:list[float]
        ,diff_var:str =""
        )-> float:
    diff_var_list = ["r","z","psi",0,1,2]
    modification_counter = 0
    df = ""#0
    if diff_var == "r" or diff_var == 0:
        df = 0
        if i + 1 == j:
            df = 2*np.pi*r[i+1]/Area[i] #*Kronecker(i+1,j) #- r[i]*Kronecker(i,j))/Area[i]
            modification_counter += 1
        if i == j :
            df = -2*np.pi*r[i]/Area[i]
            modification_counter += 1
            
    if diff_var == "z" or diff_var == 1:
        df = 0
        modification_counter += 1

    if diff_var == "psi" or diff_var == 0:
        df = 0
        if i == j:
            df = np.sin(psi[i])#*Kronecker(i,j)
            modification_counter += 1
        
    if diff_var not in diff_var_list:
        #Error handling
        print(f"\n wrong diff variable of either (r,z,psi) error in c_diff_f \n")
        exit()
    if df == "" or modification_counter > 1:
        print(f"\n df never took value \n or modification counter ={modification_counter} should < 2")
        exit()
 
    return df


def c_diff_g(
        i:int,j:int,N:int
        ,r:list[float],z:list[float],psi:list[float],Area:list[float]
        ,diff_var:str = ""
        )-> float:

    diff_var_list = ["r","z","psi",0,1,2]
    modification_counter = 0
    dg = ""#0
    if diff_var == "r" or diff_var == 0:
        dg = 0
        if i + 1 == j or i == j:
            dg = np.pi*(z[i+1]-z[i])/Area[i]
            modification_counter += 1
        
    if diff_var == "z" or diff_var == 1: 
        dg = 0
        if i+1 == j:
            dg = np.pi*(r[i+1] + r[i])/Area[i]
            modification_counter += 1
        if i == j:
            dg = -np.pi*(r[i+1] + r[i])/Area[i]
            modification_counter += 1

    if diff_var == "psi" or diff_var == 2:
        dg = 0
        if i == j:
            dg = -np.cos(psi[i])
            modification_counter += 1

    if diff_var not in diff_var_list:
            #Error handling
            print(f"\n wrong diff variable of either (r,z,psi) error in c_diff_g \n")
            exit()
    if dg == ""or modification_counter > 1:
        print(f"dg never took value \n or modification counter ={modification_counter} should < 2")
        exit()

    return dg


def c_diff(
        i:int,j:int,N:int
        ,r:list[float],z:list[float]
        ,psi:list[float],Area:list[float]
        ,diff_var =""
        ):
    
    diff_var_list = ["r","z","psi",0,1,2]

    if diff_var not in diff_var_list:
        print(f"c diff error in diff_var")
        exit()
    c_diff_val = ""

    if 0 <= i < N :
        c_diff_val = c_diff_f(
                i=i%N,j=j%N,N=N
                ,r=r,psi=psi
                ,Area=Area
                ,diff_var=diff_var
                )
        
    elif  N <= i < 2*N :
        c_diff_val = c_diff_g(
                i=i%N,j=j%N,N=N
                ,r=r,z=z,psi=psi
                ,Area=Area
                ,diff_var=diff_var
                )
    else: 
        print(f"Error c_diff_val didnt take a value because i={i} and shouldnt be higher than N={N}")
        exit()
    if c_diff_val == "":
        print(f"Error c_diff_val didnt take a value because i={i} and shouldnt be higher than N={N}")
        exit()
    
    return c_diff_val


def Epsilon_v2(
        N:int,r:list[float],z:list[float],psi:list[float],Area:list[float]
        ,print_matrix:bool = False
        ,testing:bool= False
        )->list:
    A = np.zeros(shape=(2*N,2*N),dtype=float)
    b = np.zeros(2*N,dtype=float)
    vars = ["r","z","psi"]
    
    for alpha in range(2*N):
        for beta in range(2*N):
            a = 0
            for n in range(N):
                for variables in vars:
                    a += c_diff(i=alpha,j=n,N=N,r=r,z=z,psi=psi,Area=Area,diff_var=variables)*c_diff(i=beta,j=n,N=N,r=r,z=z,psi=psi,Area=Area,diff_var=variables)                    
            A[alpha][beta] = a            

        if 0 <= alpha < N :
            b[alpha] = -constraint_f(i=alpha%N,N=N,r=r,psi=psi,Area=Area)
        elif N <= alpha < 2*N :
            b[alpha] = -constraint_g(i=alpha%N,N=N,r=r,z=z,psi=psi,Area=Area)
        else:
            print("\n Error alpha out of range \n")
            exit()
    
    if print_matrix == True:
        print(f"A: {np.shape(A)[0]}x{np.shape(A)[1]}\n ",A)
        print("b:",b)
        x = np.linalg.solve(A,b)
        #epsilon_f = x[0:N]
        #epsilon_g = x[N:2*N]
    else:
        x = np.linalg.solve(A,b)
        #epsilon_f = x[0:N]
        #epsilon_g = x[N:2*N]

    if testing == True:
        epsilon_f = x[0:N]
        epsilon_g = x[N:2*N]
        return epsilon_f,epsilon_g, A, b
    else:
        return x


def Make_variable_corrections(
        N:int
        ,r:list[float],z:list[float],psi:list[float]
        ,Area:list[float], Area_init:float
        ,Tolerence:float = 1e-10
        ,corr_max:int = 20
        ,t=""
    ):

    do_correction = False
    do_correction = check_constraints_truth(N=N,r=r,z=z,psi=psi ,Area=Area,tol=Tolerence)
    correction_count = 0

    while do_correction == True:
        correction_count += 1
        epsilon = Epsilon_v2(
                N=N, r=r, z=z ,psi=psi ,Area=Area
                        )
        scaleing = 1
        for i in range(N):      
            K_r,K_z,K_psi = 0,0,0
            for beta in range(2*N):
                K_r += epsilon[beta]*c_diff(i=beta,j=i,N=N ,r=r ,psi=psi ,z=z ,Area=Area,diff_var="r")
                
                K_z += epsilon[beta]*c_diff(i=beta,j=i,N=N ,r=r ,psi=psi ,z=z ,Area=Area,diff_var="z")
                
                K_psi += epsilon[beta]*c_diff(i=beta,j=i,N=N ,r=r,psi=psi,z=z,Area=Area,diff_var="psi")
                
            r[i] += K_r
            z[i] += K_z
            psi[i] += K_psi

        #do_correction = False
        do_correction = check_constraints_truth(N=N,r=r,z=z,psi=psi,Area=Area,tol=Tolerence)
        
        if correction_count >= corr_max:
            print(f"{corr_max} corrections, is too many corrections, we close the program. ")
            exit()
            break

    return correction_count




def Perturbation_of_inital_state(
        points_perturbed:int
        , ds:float, N:int
        ,r:list[float],z:list[float],psi:list[float]
        ,Area:list[float]
        ,delta_psi:float = -1
        ,flat = False
        ,Tolerence:float = 1e-10
        ,show_initial_condi:bool = False
        ):
    print("delta_psi=",delta_psi)
    r_unperturb = [i for i in r]
    z_unperturb = [i for i in z]

    if points_perturbed > len(psi)-1:
        print(f"Perturbing too many points \n"
              +f"len(psi)={len(psi)} and points perturbed={points_perturbed}")
        exit()
    i_start = points_perturbed#len(psi)-1 - points_perturbed
    i_stop = -1

    for i in range(0,points_perturbed):
            psi[i] += delta_psi*(points_perturbed - i)
            #psi[i] += delta_psi*np.exp(-i)#/abs(delta_psi))
            #psi[i] += delta_psi

    for i in range(i_start,i_stop,-1):
        r[i] = r[i+1] + np.cos(psi[i]+np.pi)*ds
        z[i] = z[i+1] + np.sin(psi[i]+np.pi)*ds
    
    r_perturb = [i for i in r]
    z_perturb = [i for i in z]

    Area_initial = np.sum(Area)
    corr_count = Make_variable_corrections(N=N,r=r,z=z,psi=psi,Area=Area,Area_init=Area_initial)

    print(f"number of corrections need = {corr_count}")
    if show_initial_condi == True:
        plt.figure()
        plt.plot(r_perturb,z_perturb,"-o",label="perturbed not corrected")
        plt.plot(r_unperturb,z_unperturb,"-o",label="unperturbed")
        plt.plot(r,z,"-o",label="perturbed and corrected")
        plt.legend()
        plt.draw()



if __name__ == "__main__":
    pass