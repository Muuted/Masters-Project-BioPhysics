import numpy as np
import matplotlib.pyplot as plt
import time
import os
import pandas as pd
from two_d_continues_integration import find_init_stationary_state
from Two_D_functions import Perturbation_of_inital_state, drdt_func,dzdt_func,dpsidt_func
from Two_D_functions import Langrange_multi, Make_variable_corrections, gamma
from Runge_Kutta import RungeKutta45
from two_d_data_processing import E_kin, E_pot, Xsqaured_test
from Make_movie import Make_frames, Make_video
from two_d_plot_data import plot_Epot_Ekin, plot_tot_area, plot_comparison_of_plus_minus_un_perturbed_results
from Lagrange_multipliers import Lagrange_multipliers
import multiprocessing

np.set_printoptions(legacy='1.25') #Setting the print format


class Surface_membrane:
    def __init__(self
        ,N:int = 20 ,T:float=1.0e-7 ,dt:float=2.5e-11 ,const_index:int = 0
        , dpsi:float = -0.01, dtau:float = 1.05, perturb_var:str = ""
        ,save_path:str = ""
        ):
        super().__init__()
        # Base constants
        self.N:int = N # Number of links
        self.T:float = T # [s] real time simulatied in seconds
        self.dt:float = dt # [s] time step
        self.sim_steps:int = int(self.T/self.dt) #number of steps in the simulation
        self.c0:float = 25.0 # [1/(mu m)]
        self.k:float = 80.0 # [zJ]
        self.kG:float = None # [zJ]
        self.eta:float = 1.0 # [kg/(m*s)]
        self.alpha:float = None #
        self.Area_init:float = None # [(mu m)^2] 
        self.ds:float = 1.5e-2#/(self.N/20) # [mu m]
        self.dpsi_perturb:float = dpsi # [rad]
        self.dtau_perturb:float = dtau # [nN]
        self.num_perturb:int = 10 #int(self.N/2) # Number of perturbed points
        self.r0:float = 5.0e-2 # [mu m] if the init config is just flat, this is the initial radius of the hole.
        self.L = 100 # [mu m] some times used for the total length of the membrane

        # Phase space variables 
        self.const_index:int = const_index
        self.phase_space_names:list = ["triangle\\","plus\\","cross\\"]
        self.tilde_sigma_list:list = [
        0.0253164556962025 
        ,0.29873417721519
        ,0.116455696202532 
        ]
        self.tilde_tau_list:list = [
        1.31578947368421
        ,4.47368421052632
        ,1.0
        ]
        self.psi2_list:list = [
        -2.83260429562395e-7
        ,-2.23344534748962e-6
        ,-2.26474921864332e-8
        ]

        # Constant variables
        self.lc = 1/self.c0 # [mu m]
        self.tilde_sigma:float = self.tilde_sigma_list[self.const_index] # [dimless]
        self.sigma_c:float = self.k*self.c0**2
        self.sigma:float = self.tilde_sigma*self.sigma_c # [zJ/(mu m)^2]

        self.tilde_tau:float = self.tilde_tau_list[self.const_index] # [dimless]
        self.tau_c:float = self.k*self.c0 # [nN]
        self.tau:float = self.tilde_tau*self.tau_c # [nN]

        self.psi2:float = self.psi2_list[self.const_index] # [rad]

        # Lists for variables
        self.r_list:list = np.zeros(shape=(self.sim_steps,self.N+1),dtype=float)
        self.z_list:list = np.zeros(shape=(self.sim_steps,self.N+1),dtype=float)
        self.psi_list:list = np.zeros(shape=(self.sim_steps,self.N),dtype=float)
        self.Area_list = np.zeros(self.N,dtype=float)

        self.r_unperturbed:list = []
        self.z_unperturbed:list = []
        self.psi_unperturbed:list = []     

        self.Potential_E = np.zeros(self.sim_steps,dtype=float)
        self.Potential_E_before_correction = np.zeros(self.sim_steps,dtype=float)
        self.Kinetic_E = np.zeros(self.sim_steps,dtype=float)
        self.Xsqre = np.zeros(self.sim_steps,dtype=float)
        self.correct_count_list = np.zeros(self.sim_steps)

        # Setting up program
        self.var_corr_tol:float = 1e-5 # The tolerence for when to use variable correction
        self.do_correction:bool = True # deciced whether or not to use variable corrections
        self.perturb:bool = False # decieds if we are going perturb the initial state or not   
        self.var_perturb_options:list = ["psi","tau"] # the types of options for variable perturbations    
        self.var_perturb_choice:str = perturb_var  # The variable choosen
        self.start_flat:bool = False #Choses whether or not the simulation starts from a flat state or not.
        self.use_phase_diagram: bool = False # decides if using the phase space diagram to find the simulation values
        self.phasespace_x = ""
        self.phasespace_y = ""
        self.phasespace_fignum = ""
        self.const_length_diff_N_density = True # if True means that the value of ds changes to conserve length for different N
        self.load_external_init_config = False # This allows the user to load an initial configuration, r,z,psi of shape (1,N)
        self.external_args = [] #The args needed from external user.

        # Printing choices and paths saving
        self.integration_method:str = "RK4" #Type of integration scheme
        if self.var_perturb_choice == "psi":
            a = f"2D sim results\\object results\\T={self.T}\\" + self.phase_space_names[self.const_index] + f"(N,T,dt,dpsi)=({self.N},{self.T:0.1e},{self.dt:0.1e},{self.dpsi_perturb:0.1e})\\"
        elif self.var_perturb_choice == "tau":
            a = f"2D sim results\\object results\\T={self.T}\\" + self.phase_space_names[self.const_index] + f"(N,T,dt,dtau)=({self.N},{self.T:0.1e},{self.dt:0.1e},{self.dtau_perturb:0.1e})\\"
        else:
            a = f"2D sim results\\object results\\T={self.T}\\" + self.phase_space_names[self.const_index] + f"(N,T,dt)=({self.N},{self.T:0.1e},{self.dt:0.1e})\\"
        
        if save_path == "":
            self.save_path:str = a
        else:
            self.save_path:str = save_path #+ f"(N,T,dt)=({self.N},{self.T:0.1e},{self.dt:0.1e})\\"
        self.save_figs_path:str = "figures and movie\\"
        self.figs_for_video_path:str = "figures for video\\"
        self.df_name:str = "2D Surface sim.pkl"
        self.print_progress:bool = True
        self.print_scale:int = int((self.sim_steps-2)/1000)
        self.print_constants:bool = True

        # 
        self.save_data:bool = True #choses whether or not the simulation data is saved.
        self.show_stationary_state:bool = True # Showing the initial configuration before simulation start
        self.init_config_show_time:float = 2 # Showing the initial configuration
        self.close_final_plots:bool = False
        self.overflow_err = False

        # Making movies and plots
        self.make_movie:bool = False #True # make_movie 
        self.make_plots:bool = False #True # make_plots
        self.fps_movie:int = 24 # chooses the frames per second in the movie of the dynamics        

    def setup_simulation(self):

        if self.const_length_diff_N_density == True:
            self.ds = self.ds*(20/self.N)
            print("scaled N")        
        
        # scaling parameters        
        rs2 = 20*self.lc 
        zs2 = 0
        s0, sN = 0, 50*self.lc
        #print(self.N)
        #Initiating the inital state of the membrane
        psi,r,z, r_contin, z_contin, alpha = find_init_stationary_state(
                sigma=self.sigma ,k=self.k ,c0=self.c0 ,tau=self.tau ,ds=self.ds
                ,psi_L=self.psi2 ,r_L=rs2 ,z_L=zs2 ,s0=s0 ,sN=sN
                ,total_points = self.N
                ,use_phase_space= self.use_phase_diagram
                ,use_fig_number = self.phasespace_fignum
            )
        self.alpha = (self.c0 - (psi[1] - psi[0])/self.ds )*r[0]/np.sin(psi[0]) - 1
        self.kG = self.k*self.alpha   

        if self.start_flat == True:
            #self.alpha = 0.75
            #self.kG = 0.75*self.k

            for i in range(int(self.N+1)):
                self.r_list[0][i] = self.r0 + i*self.ds
                self.z_list[0][i] = 0
                if i < self.N:
                    self.psi_list[0][i] = 0
        else:

            """------ variables list ---------"""
            for i in range(int(self.N+1)):
                if i < self.N :
                    #print(i,self.N, self.ds,len(psi))
                    self.psi_list[0][i] = psi[i]
                self.r_list[0][i] = r[i]
                self.z_list[0][i] = z[i]
        
        for i in range(int(self.N+1)):
            if i < self.N :
                self.Area_list[i] =  np.pi*(self.r_list[0][i+1] + self.r_list[0][i])*np.sqrt( 
                    (self.r_list[0][i+1] - self.r_list[0][i])**2
                    + (self.z_list[0][i+1] - self.z_list[0][i])**2 
                    )
                if self.Area_list[i] == 0 :
                    print(f"Area[{i}]=0")
                    exit()

        self.r_unperturb = [i for i in self.r_list[0]]
        self.z_unperturb = [i for i in self.z_list[0]]
        self.psi_unperturb = [i for i in self.psi_list[0]]

        if self.perturb == True:
            if self.var_perturb_choice == "psi":#self.var_perturb_options[0]:
                print("psi perturb chosen")
                Perturbation_of_inital_state(
                    points_perturbed = self.num_perturb 
                    ,ds = self.ds
                    ,N = self.N
                    ,r = self.r_list[0]
                    ,z = self.z_list[0]
                    ,psi = self.psi_list[0]
                    ,Area = self.Area_list
                    ,delta_psi = self.dpsi_perturb
                    ,show_initial_condi = True
                )
            elif self.var_perturb_choice == "tau":#self.var_perturb_options[1]:
                print("tau perturb chosen")
                self.tau = self.tau*self.dtau_perturb
            else:
                print("No perturb method was chosen, yet perturbation was set to true. Some mistake was made")
                exit()

        if self.show_stationary_state == True:
            plt.figure()
            font_size = 10
            #plt.plot(r_contin,z_contin,marker="",label="integration")
            plt.plot(self.r_list[0],self.z_list[0],"o-",label="Discreet")
            plt.plot(self.r_list[0][0],self.z_list[0][0],"o",color="k",label="s1")
            plt.plot(self.r_list[0][len(self.r_list[0])-1],self.z_list[0][len(self.r_list[0])-1],"o",color="y",label="s2")
            if self.start_flat == False:
                plt.plot(r_contin,z_contin,linestyle="--",marker="",color="k",label="integration")
            if self.perturb == True:
                plt.plot(self.r_unperturb,self.z_unperturb,"-o",label="unperturbed")
            plt.xlim(min(self.r_list[0])-self.ds, max(self.r_list[0])+self.ds)
            ceil = max(self.r_list[0]) - min(self.r_list[0]) + 2*self.ds
            plt.ylim(-ceil/10, 9*ceil/10)
            #plt.xlim(3,83)
            #plt.ylim(-10,70)
            plt.xlabel("r",fontsize=font_size)
            plt.ylabel("z",fontsize=font_size)
            plt.title(
                f"Quick peak at the neck configuration before dynanic simulation "
                ,fontsize=font_size
                )
            plt.legend()
            plt.grid()
            plt.draw()
            plt.pause(0.5)
            if self.use_phase_diagram == False:
                plt.pause(self.init_config_show_time)
            elif self.use_phase_diagram == True:
                next = input(f"give input, for close figure write: run sim , for exit program write: end, \n input=")
                if str(next) == "end":
                    print("you close the program.")
                    exit()
            plt.close("all")

    def load_external_config(self):
        
        pass

    def print_consts(self):
        if self.print_constants == True:
            print(
            f" \n \n"
            + "------------- Constant used in Simulation -------------- \n "
            + f"    number of chain links N: {self.N} \n " 
            + f"    Total sim time = {self.T:.1e} [s] \n "
            + f"    dt = {self.dt:0.1e} [s] \n "
            + f"    k = {self.k:0.1e}  [zJ] \n "
            + f"    kG = {self.kG:.1e}  [zJ] \n "
            + f"    c0 = {self.c0:.1e}  [1/(mu m)] \n "
            + f"    tau = {self.tau:0.1e}  [nN] \n "
            + f"    sigma = {self.sigma:0.1e} [zJ/(mu m)^2] \n "
            + f"    ds = {self.ds:0.1e} [mu m] \n "
            + f"    eta = {self.eta:0.1e} [(mu g)/(mu m * s)] \n "
            + f"    gamma(i!=0) = {gamma(i=2,ds=self.ds,eta=self.eta):.2e} [(mu g)/s]  \n "
            + f"    Sim steps = {self.sim_steps:0.1e} \n "
            + f"    dpsi = {self.dpsi_perturb:0.1e} [rad] \n "
            + f"    alpha = {self.alpha:0.1e}  \n "
            + f"    var tol = {self.var_corr_tol:0.1e} \n "
            + f"    Integration Scheme = {self.integration_method} \n "
            + f"------------------------------------------------------- \n \n "
        )

    def dynamics(self):
        start_time = time.time()
        print_scale = (self.sim_steps-2)/1000
        self.Area_init = np.sum(self.Area_list)

        integration_options = ["Euler","RK4"]
        if self.integration_method not in integration_options:
            print("No integration method choosen correctly")
            exit()
        

        self.Potential_E_before_correction[0] = E_pot(
                N=self.N,k=self.k,kG=self.kG,tau=self.tau,c0=self.c0
                ,r=self.r_list[0],psi=self.psi_list[0],Area=self.Area_list
                )
        self.Potential_E[0] = E_pot(
                N=self.N,k=self.k,kG=self.kG,tau=self.tau,c0=self.c0
                ,r=self.r_list[0],psi=self.psi_list[0],Area=self.Area_list
                )
        self.Kinetic_E[0] = E_kin(N=self.N,t=0,dt=self.dt,r=self.r_list,z=self.z_list,Area=self.Area_list)
        self.Xsqre[0] = Xsqaured_test(
                N=self.N
                ,r_init=self.r_unperturb,z_init=self.z_unperturb,psi_init=self.psi_unperturb
                ,r=self.r_list[0],z=self.z_list[0],psi=self.psi_list[0]
                )
        

        if self.print_constants == True:
            self.print_consts()
        
        print(f"integration method={self.integration_method}")
        print(" \n \n Simulation progressbar: \n")

        for t in range(self.sim_steps-1):
            if int(t%print_scale) == 0 and self.print_progress == True:
                time1 = time.time()-start_time
                time_left = (time1/(t+1))*(self.sim_steps-t)
                time_h_start, time_m_start, time_s_start = int((time1/60**2)%24), int((time1/60)%60), int(time1%60)
                time_h_end, time_m_end, time_s_end = int((time_left/60**2)%24), int((time_left/60)%60), int(time_left%60)
                print(
                    f"completion : {round(t/(print_scale*10),1)}%       " 
                    +f"Time since start = {time_h_start}h {time_m_start}m {time_s_start}s        "
                    +f" Estimated time left = {time_h_end}h {time_m_end}m {time_s_end}s   "
                    , end="\r"
                )
            #t1,t2 = t%2, (t+1)%2
            
            lambs,nus = Lagrange_multipliers(#Langrange_multi(
                    N=self.N,k=self.k,c0=self.c0,sigma=self.sigma
                    ,kG=self.kG,tau=self.tau,ds=self.ds,eta=self.eta
                    ,Area=self.Area_list
                    ,psi=self.psi_list[t]
                    ,radi=self.r_list[t]
                    ,z_list=self.z_list[t]
                )
            
            if self.integration_method == "Euler":
                for i in range(self.N+1):
                    if i == self.N:
                        self.z_list[t+1][i] = self.z_list[t][i]
                        self.r_list[t+1][i] = self.r_list[t][i]
                        #psi[t+1][i] = psi[t][i]
                        if self.r_list[t+1][i] != self.r_list[t+1][i] or self.z_list[t+1][i] != self.z_list[t+1][i]:
                            if self.overflow_err == False:
                                print("\n ---- Overflow error occurred ---- \n")
                                self.overflow_err = True

                    if i < self.N:
                        self.z_list[t+1][i] = self.z_list[t][i] + self.dt*dzdt_func(i=i,ds=self.ds,eta=self.eta,Area=self.Area_list,radi=self.r_list[t],nu=nus)

                        drdt = drdt_func(
                                    i=i
                                    ,N=self.N,k=self.k,c0=self.c0,sigma=self.sigma,kG=self.kG,tau=self.tau,ds=self.ds,eta=self.eta
                                    ,Area=self.Area_list,psi=self.psi_list[t],radi=self.r_list[t],z_list=self.z_list[t]
                                    ,lamb=lambs,nu=nus
                                    )
                        self.r_list[t+1][i] = self.r_list[t][i] + self.dt*drdt

                        dpsidt = dpsidt_func(
                                        i=i
                                        ,N=self.N,k=self.k,c0=self.c0,sigma=self.sigma,kG=self.kG,tau=self.tau,ds=self.ds,eta=self.eta
                                        ,Area=self.Area_list,psi=self.psi_list[t],radi=self.r_list[t],z_list=self.z_list[t]
                                        ,lamb=lambs,nu=nus
                                        )
                        self.psi_list[t+1][i] = self.psi_list[t][i] + self.dt*dpsidt

                        if self.r_list[t+1][i] != self.r_list[t+1][i] or self.z_list[t+1][i] != self.z_list[t+1][i] or self.psi_list[t+1][i] != self.psi_list[t+1][i]:
                            if self.overflow_err == False:
                                print("\n ---- Overflow error occurred ---- \n")
                                self.overflow_err = True
                                
            elif self.integration_method == "RK4":
                kr,kz,kpsi = RungeKutta45(
                    N=self.N,dt=self.dt,k=self.k,c0=self.c0, sigma=self.sigma
                    ,kG=self.kG ,tau=self.tau, ds=self.ds,eta=self.eta
                    ,Area=self.Area_list,lamb=lambs,nu=nus
                    ,psi_init=self.psi_list[t],r_init=self.r_list[t], z_init=self.z_list[t]
                    )
                for i in range(self.N+1):
                    if i == self.N:
                        self.z_list[t+1][i] = self.z_list[t][i]
                        self.r_list[t+1][i] = self.r_list[t][i]
                        #psi[t+1][i] = psi[t][i]
                        if self.r_list[t+1][i] != self.r_list[t+1][i] or self.z_list[t+1][i] != self.z_list[t+1][i]:
                            if self.overflow_err == False:
                                print("\n ---- Overflow error occurred ---- \n")
                                self.overflow_err = True

                    if i < self.N:
                        self.r_list[t+1][i] = self.r_list[t][i] + (self.dt/6)*(kr[1][i] + 2*kr[2][i] + 2*kr[3][i] + kr[4][i])
                        self.z_list[t+1][i] = self.z_list[t][i] + (self.dt/6)*(kz[1][i] + 2*kz[2][i] + 2*kz[3][i] + kz[4][i])
                        self.psi_list[t+1][i] = self.psi_list[t][i] + (self.dt/6)*(kpsi[1][i] + 2*kpsi[2][i] + 2*kpsi[3][i] + kpsi[4][i])

                        if self.r_list[t+1][i] != self.r_list[t+1][i] or self.z_list[t+1][i] != self.z_list[t+1][i] or self.psi_list[t+1][i] != self.psi_list[t+1][i]:
                            if self.overflow_err == False:
                                print("\n ---- Overflow error occurred ---- \n")
                                self.overflow_err = True

            else:
                print("No integration method was choosen, program terminates")
                exit()

            self.Potential_E_before_correction[t+1] = E_pot(
                N=self.N,k=self.k,kG=self.kG,tau=self.tau,c0=self.c0
                ,r=self.r_list[t+1],psi=self.psi_list[t+1],Area=self.Area_list
                )

            if self.do_correction == True:
                correction_count = Make_variable_corrections(
                    N = self.N
                    ,r = self.r_list[t+1] ,z = self.z_list[t+1], psi = self.psi_list[t+1] 
                    ,Area = self.Area_list
                    ,Area_init = self.Area_init
                    ,Tolerence = self.var_corr_tol
                    ,corr_max = 10
                    ,t = t
                )
                self.correct_count_list[t+1] = correction_count        

            self.Potential_E[t+1] = E_pot(
                N=self.N,k=self.k,kG=self.kG,tau=self.tau,c0=self.c0
                ,r=self.r_list[t+1],psi=self.psi_list[t+1],Area=self.Area_list
                )
            
            self.Kinetic_E[t+1] = E_kin(N=self.N,t=t,dt=self.dt,r=self.r_list,z=self.z_list,Area=self.Area_list)

            self.Xsqre[t+1] = Xsqaured_test(
                N=self.N
                ,r_init=self.r_unperturb,z_init=self.z_unperturb,psi_init=self.psi_unperturb
                ,r=self.r_list[t+1],z=self.z_list[t+1],psi=self.psi_list[t+1]
                )

        
        print("\n")        

        if self.save_data == True:        
            df = pd.DataFrame({
                'psi': [self.psi_list],
                "r": [self.r_list],
                "z": [self.z_list],
                "r unperturbed": [self.r_unperturb],
                "z unperturbed": [self.z_unperturb],
                "psi unperturbed": [self.psi_unperturb],
                "area list": [self.Area_list],
                "Epot":[self.Potential_E],
                "Ekin":[self.Kinetic_E],
                "Epot before correction": [self.Potential_E_before_correction],
                "Chi squared test": [self.Xsqre],
                #'lambs': [lambs_save],
                #'nus': [nus_save],
                "L" : self.L,
                "r0": self.r0,
                "N": self.N,
                "c0": self.c0,
                "k": self.k,
                "kG": self.kG,
                "sigma": self.sigma,
                "tau": self.tau,
                "eta": self.eta,
                "sim_steps": self.sim_steps,
                "dt": self.dt,
                "ds": self.ds,
                "dpsi perturb": self.dpsi_perturb,
                "dtau perturb": self.dtau_perturb,
                "perturb var": self.var_perturb_choice,
                "gam(i=0)": gamma(0,ds=self.ds,eta=self.eta),
                "gam(i>0)": gamma(2,ds=self.ds,eta=self.eta),
                "correction count": [self.correct_count_list],
                "Tolerance":self.var_corr_tol,
                "Overflow err": self.overflow_err,
                "integration method":self.integration_method,
                "simulation time [s]":int(time.time()-start_time)
                            })

            print(self.save_path + self.df_name)
            if not os.path.exists(self.save_path):
                os.makedirs(self.save_path)
            df.to_pickle(self.save_path + self.df_name)
            
            if not os.path.exists(self.save_path +"figures and video"):
                os.makedirs(self.save_path + "figures and video")

    def plotting_n_movie_data(self):
        if self.make_movie == True:
            Make_frames(
                data_path = self.save_path
                ,figs_save_path = self.figs_for_video_path
                ,df_name = self.df_name
                #,tot_frames= 50
            )
            Make_video(
                output_path = self.save_path + self.save_figs_path
                ,input_path = self.figs_for_video_path
                ,video_name = self.df_name
                ,fps = self.fps_movie
            )
        
        if self.make_plots == True:
            plot_Epot_Ekin(
                data_path = self.save_path
                ,df_name = self.df_name
                ,output_path = self.save_path + self.save_figs_path
            )
            plot_tot_area(
                data_path = self.save_path
                ,df_name = self.df_name
                ,output_path = self.save_path + self.save_figs_path
            )

            if self.close_final_plots == False:
                plt.show()
            else:
                plt.draw()
                plt.pause(2)

    def phase_space_choice(self):
        df_name = "Matlab data"
        data_path2 = "2D sim results\\"

        if not os.path.exists(data_path2 + df_name):
            print(f"The file does not exsist")
        elif os.path.exists(data_path2 + df_name):
            df = pd.read_pickle(data_path2 + df_name)

        fig, ax = plt.subplots()
        cmap = plt.cm.coolwarm
        wm = plt.get_current_fig_manager()
        #wm.window.style('zoomed')
        #wm.window.state('zoomed')

        max_tau = -1e5
        min_tau = 1e5

        diff_taus = []
        """------------ Finding the max and min values of Tau for cmap plot-------------------------------"""
        for n in range(len(df["tau_list_Curved"])):
            for i in range(len(df["tau_list_Curved"][n])):
                if df["tau_list_Curved"][n][i] < min_tau:
                    min_tau = df["tau_list_Curved"][n][i]
                
                if df["tau_list_Curved"][n][i] > max_tau:
                    max_tau = df["tau_list_Curved"][n][i]

                if df["tau_list_Curved"][n][i] not in diff_taus:
                    diff_taus.append(df["tau_list_Curved"][n][i])

        for n in range(len(df["tau_list_Flat"])):
            for i in range(len(df["tau_list_Flat"][n])):
                if df["tau_list_Flat"][n][i] < min_tau:
                    min_tau = df["tau_list_Flat"][n][i]
                
                if df["tau_list_Flat"][n][i] > max_tau:
                    max_tau = df["tau_list_Flat"][n][i]
                
                if df["tau_list_Flat"][n][i] not in diff_taus:
                    diff_taus.append(df["tau_list_Flat"][n][i])
        
        """-------------------------Finding the three marks triangle/plus/cross -----------------------------------"""
        """ Triangle """
        pos_A_min_x1,pos_A_max_x1, pos_r1_min_y1, pos_r1_max_y1 = 23.6 ,23.8 ,2.83 ,2.84
        """ plus """
        pos_A_min_x2,pos_A_max_x2, pos_r1_min_y2, pos_r1_max_y2 = 3.4 ,3.6 ,7.35 ,7.45
        """ cross """
        neg_A_min_x3,neg_A_max_x3,neg_r1_min_y3,neg_r1_max_y3 = -1.08 ,-1.06 ,0.64 ,0.66
        n_neg_A ,n_pos_A = [] ,[]
        i_neg ,i_pos = [],[]

        for n in range(len(df["ExcessAreaCurved"])):
            for i in range(len(df["ExcessAreaCurved"][n])):
                if len(df["ExcessAreaCurved"]) > 0:
                    A_curve = df["ExcessAreaCurved"][n][i]
                    r1_curve = df["r1Curved"][n][i]
                    if pos_A_min_x1 < A_curve < pos_A_max_x1 and pos_r1_min_y1 < r1_curve < pos_r1_max_y1:
                        n_pos_A.append(n)
                        i_pos.append(i)
                    if pos_A_min_x2 < A_curve < pos_A_max_x2 and pos_r1_min_y2 < r1_curve < pos_r1_max_y2:
                        n_pos_A.append(n)
                        i_pos.append(i)

        for n in range(len(df["ExcessAreaFlat"])):
            for i in range(len(df["ExcessAreaFlat"][n])):
                if len(df["ExcessAreaFlat"]) > 0:
                    A_Flat = df["ExcessAreaFlat"][n][i]
                    r1_Flat = df["r1Flat"][n][i]
                    if neg_A_min_x3 < A_Flat < neg_A_max_x3 and neg_r1_min_y3 < r1_Flat < neg_r1_max_y3:
                        n_neg_A.append(n)
                        i_neg.append(i)

        i = 0
        plot_tau_ref = []

        for n in range(len(df["ExcessAreaCurved"])):
            if len(df["tau_list_Curved"][n]) > 0 :
                i = (df["tau_list_Curved"][n][0]-1)/(max_tau - 1)
                if df["tau_list_Curved"][n][0] not in plot_tau_ref:
                    ax.plot(
                        [i/self.lc**2 for i in df["ExcessAreaCurved"][n]]
                        ,[i/self.lc for i in df["r1Curved"][n]],".-",color=cmap(i)
                        ,label=r"$\tau \approx$"+f"{(df["tau_list_Curved"][n][0]/self.lc**2)/1000:0.1f} nN"
                        )
                    plot_tau_ref.append(df["tau_list_Curved"][n][0])
                else:
                    ax.plot(
                        [i/self.lc**2 for i in df["ExcessAreaCurved"][n]]
                        , [i/self.lc for i in df["r1Curved"][n]]
                        ,".-",color=cmap(i)
                        )

        for n in range(len(df["ExcessAreaFlat"])):
            if len(df["tau_list_Flat"][n]) > 0:
                i = (df["tau_list_Flat"][n][0]-1)/(max_tau -1)
                if df["tau_list_Flat"][n][0] not in plot_tau_ref:
                    ax.plot(
                        [i/self.lc**2 for i in df["ExcessAreaFlat"][n]]
                        ,[i/self.lc for i in df["r1Flat"][n]]
                        ,".-",color=cmap(i),label=r"$\tau \approx$"+f"{(df["tau_list_Flat"][n][0]/self.lc**2)/1000:0.1f} nN"
                        )
                    plot_tau_ref.append(df["tau_list_Flat"][n][0])
                else:
                    ax.plot(
                        [i/self.lc**2 for i in df["ExcessAreaFlat"][n]]
                        ,[i/self.lc for i in df["r1Flat"][n]]
                        ,".-",color=cmap(i))
            
        markers = ["^","+","x"]
        markers_latex = [r"$\triangle$",r"$\plus$" ,r"$\times$"]
        ax.plot(
            df["ExcessAreaCurved"][n_pos_A[0]][i_pos[0]]/self.lc**2
            ,df["r1Curved"][n_pos_A[0]][i_pos[0]]/self.lc
            ,marker=markers[0] 
            ,color="k"
            ,markersize = 15
            ,mfc = "none"
            )

        ax.plot(
            df["ExcessAreaCurved"][n_pos_A[1]][i_pos[1]]/self.lc**2
            ,df["r1Curved"][n_pos_A[1]][i_pos[1]]/self.lc
            ,marker=markers[1] 
            ,color="k"
            ,markersize = 12
            )
        ax.plot(
            df["ExcessAreaFlat"][n_neg_A[0]][i_neg[0]]/self.lc**2
            ,df["r1Flat"][n_neg_A[0]][i_neg[0]]/self.lc
            ,marker=markers[2] 
            ,color="k"
            ,markersize = 12
            ,mfc = "none"
            )

        plt.subplots_adjust(
            #left=0.05
            #,right=0.5
            wspace=-0.6
        )
        ax.legend(
            bbox_to_anchor=(1.02, 1)
            ,loc='upper left'
            ,borderaxespad=0
            ,fontsize=15
                    )
        ax.vlines(x=0,ymin=-100,ymax=300,colors="k",linestyles="--")
        #plt.xlim(-60, 50)
        xlim_min = -1.2e4
        xlim_max = 2.9e4
        ax.set_xlim(-1.2e4, 2.9e4)
        ylim_max = 230
        ax.set_ylim(0 ,ylim_max)
        ax.set_title(
            "Edge radius (r1) vs Excess Area"
            #+r"[$\tau]=\frac{\mu g\cdot \mu m }{s^2}$"
            ,fontsize=15)
        ax.set_xlabel(r"$\Delta A_{Excess} = A_{membrane} - A_{disc} $  [$\mu m^2$]",fontsize=15)
        ax.set_ylabel(r"Edge radius (r1) [$\mu m$]",fontsize=15)
        ax.grid()

        plt.draw()
        plt.pause(2)
        running = True
        choose_point = True
        calc_dist = False
        while running:
            if choose_point == True:
                print(f"Choose (x,y) coordinates")
                if self.phasespace_x == "" and self.phasespace_y == "":
                    x1 = input("Excess Area=")
                    y1 = input("r1=")
                elif isinstance(self.phasespace_x ,(float, int)) == True and isinstance(self.phasespace_y , (float,int)) == True:
                    x1 = self.phasespace_x
                    y1 = self.phasespace_y
                else:
                    print(
                        "The chosen values for the phase space coordinates is not either float or integer"
                        +f"program will terminate"
                        )
                    exit()

                choose_point = False
                calc_dist = True

            if calc_dist == True:
                x,y = float(x1), float(y1)
                ax.plot(x,y,marker="o",markersize=5)
                plt.draw()

                min_dist = 1e20
                r1,ExcessArea = 1e20,1e20
                data_set = ["Curved","Flat"]
                data_type:str = ""
                n_point,i_point = None,None
                points_checked = 0
                for data in data_set:
                    for n in range(len(df["ExcessArea" + data])):
                        if len(df["ExcessArea" + data][n]) > 0 :
                            for i in range(len(df["ExcessArea" + data][n])):
                                points_checked += 1
                                r1_test = df["r1" + data][n][i]/self.lc
                                ExA_test = df["ExcessArea" + data][n][i]/self.lc**2

                                dist = np.sqrt( 
                                    ( (x - ExA_test)/(xlim_max + abs(xlim_min)) )**2 
                                    + ( (y - r1_test)/ylim_max )**2
                                    ) 
                                if dist < min_dist:
                                    min_dist = dist
                                    r1 = r1_test
                                    ExcessArea = ExA_test
                                    data_type = data
                                    n_point = n
                                    i_point = i
                                    
                ax.plot([x,ExcessArea],[y,r1],color="k",markersize=1,label="closes point")
                plt.draw()
                calc_dist = False
            
            if running==True:
                next = str(input("choose one (end/e, new point/n, correct point/c) :"))
                if next == "end" or next =="e":
                    running = False

                if next == "new point" or next == "n":
                    choose_point = True

                if next == "correct point" or next == "c":
                    self.sigma = df["sigma_list_" + data_type][n_point][i_point]
                    self.tau = df["tau_list_" + data_type][n_point][i_point]
                    self.psi2 = df["psi_L_list_" + data_type][n_point][i_point]
                    a1 =f"2D sim results\\phase space choice\\"
                    a2 = f"(N,T,dt)=({self.N},{self.T:0.1e},{self.dt:0.1e})\\"
                    a3 = f"(sigma,tau,psi2)=({self.sigma:0.2e},{self.tau:0.2e},{self.psi2:0.2e})\\"
                    self.save_path = a1 + a2 + a3
                    plt.close()
                    running = False
                
    def run_sim(self):
            if self.use_phase_diagram == True:
                self.phase_space_choice()

            if self.load_external_init_config == False:
                self.setup_simulation()
            elif self.load_external_init_config == True:
                print("This is not implementet yet. Program termines")
                exit()

            self.dynamics()
            self.plotting_n_movie_data()



def multi_process(Tot_time:float,cpu_cores:int=5, sim_index:int=0,N:int=20,dt:float=2.5e-11
                  ,inte_scheme = "RK4"
                  ):
    perturb_bool_list = [False,True,True,True,True]
    perturb_list_psi = [0 ,0    ,0      ,0.01 ,-0.01]
    perturb_list_tau = [0 ,1.0+0.05 ,1.0-0.05 ,0    ,0    ]
    perturb_var_choice = ["","tau","tau","psi","psi"]

    process = []
    for i in range(cpu_cores):
        membrane = Surface_membrane(
            N = N
            ,T = Tot_time
            ,dt = dt#1.25e-11
            ,const_index=sim_index
            ,dpsi = perturb_list_psi[i]
            ,dtau = perturb_list_tau[i]
            ,perturb_var = perturb_var_choice[i]
        )
        membrane.init_config_show_time = 10
        membrane.perturb = perturb_bool_list[i]
        membrane.print_constants = False
        membrane.integration_method = inte_scheme
        #membrane.dpsi_perturb *= perturb_list_psi[i]
        #membrane.var_perturb_choice = perturb_var_choice[i]
        #membrane.use_phase_diagram = True
        #membrane.phase_space_choice()
        #membrane.setup_simulation()
        process.append(multiprocessing.Process(
            #target=membrane.dynamics#run_sim
            target=membrane.run_sim
        ))

    for p in process:
        p.start()
    
    for p in process:
        p.join()


def plotting_multi_process_results(
        path:str = "2D sim results\\object results compare with thesis data\\"
        ,make_movie= True
        ,make_figures = True
        ,make_comparison_figs = True
    ):
    #path = "2D sim results\\object results T=1e-06\\"
    #path = "2D sim results\\object results T=1e-06\\plus\\(N,T,dt,dtau)=(20,1.0e-06,2.5e-11,-5.0e-02)\\"
    directory_list = list()

    for root, dirs, files in os.walk(path, topdown=False):
        for df_name in files:
            if ".pkl" in df_name:
                data_path = root + "\\"
                directory_list.append(data_path+ files[0])
                #print(data_path + df_name)
                if make_movie == True:
                    Make_frames(
                        data_path = data_path
                        ,figs_save_path = data_path + "figues for video\\"
                        ,df_name = df_name
                        ,tot_frames = 200 #100
                    )
                    Make_video(
                        output_path=data_path + "figures and video\\"
                        ,input_path=data_path + "figues for video\\"
                        ,video_name= "surface video"
                        ,fps = 12 #24
                    )
                if make_figures == True:
                    plot_Epot_Ekin(
                        data_path=data_path
                        ,df_name=df_name
                        ,output_path=data_path + "figures and video\\"
                    )
                    plot_tot_area(
                        data_path=data_path
                        ,df_name=df_name
                        ,output_path=data_path + "figures and video\\"
                    )
                if make_comparison_figs == True:
                    plot_comparison_of_plus_minus_un_perturbed_results(
                        path=data_path
                    )
                plt.draw()
                plt.pause(2)
                plt.close("all")


if __name__ == "__main__":
    #multi_process(cpu_cores=5,Tot_time=1e-7,N=30,sim_index=0,dt=1.7e-11)
    #multi_process(cpu_cores=5,Tot_time=3e-7,N=20,sim_index=1,dt=2.5e-11)
    #multi_process(cpu_cores=5,Tot_time=1e-7,N=20,sim_index=2,dt=2.5e-11)
    plotting_multi_process_results(path="2D sim results\\object results\\RK4\\")

    exit()
    N = 35
    save_path = f"2D sim results\\obj\\plus\\N={N}\\"
    save_path = f"2D sim results\\obj\\plus larger ds\\N={N}\\"
    membrane = Surface_membrane(
        T=1e-7
        ,dt=1e-11
        ,const_index = 1
        ,N=N
        ,save_path= save_path
        )
    membrane.var_corr_tol = 1e-4
    
    #ds_init = membrane.ds
    membrane.ds = 1.5e-2*(20/30)
    membrane.const_length_diff_N_density = False

    membrane.use_phase_diagram = False
    membrane.phasespace_x = 15000
    membrane.phasespace_y = 125
    membrane.phasespace_fignum = 4

    membrane.init_config_show_time = 3
    membrane.perturb = False
    #membrane.var_perturb_choice = membrane.var_perturb_options[1]

    membrane.make_movie = True
    membrane.make_plots = True
    #membrane.phase_space_choice()
    #membrane.setup_simulation()
    
    membrane.run_sim()
    #membrane.plotting_n_movie_data()
    
