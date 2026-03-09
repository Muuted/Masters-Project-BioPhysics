import numpy as np


np.set_printoptions(legacy='1.25') #Setting the print format



class Surface_membrane():
    def __init__(self
        ,N:int = 30 ,T:float=1.0e-7 ,dt:float=2.5e-11
        ):
        super().__init__()
        # Base constants
        self.T:float = T # [s] real time simulatied in seconds
        self.dt:float = dt #[s] time step
        self.sim_steps:int = int(self.T/self.dt) #number of steps in the simulation
        self.N:int = N # Number of links
        self.c0:float = 25.0 # [1/(mu m)]
        self.k:float = 80.0 #[zJ]
        self.kG:float = None # [zJ]
        self.eta:float = 1.0 # [kg/(m*s)]
        self.alpha:float = None #
        self.Area_init:float = None # [(mu m)^2] 
        self.ds:float = 1.5e-2/(self.N/20) # [mu m]
        self.dpsi_perturb:float = -0.02 # [rad]

        # Phase space variables 
        self.phase_space_names:list = ["triangle","plus","cross"]
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
        self.tilde_sigma:float = 0 # [dimless]
        self.sigma:float = 0 # [zJ/(mu m)^2]

        self.tilde_tau:float = 0 # [dimless]
        self.tau:float = 0 # [zJ/(mu m)]

        self.psi2:float = 0 # [rad]

        # Lists for variables
        self.r_list:list = np.zeros(shape=(self.sim_steps,self.N+1),dtype=float)
        self.z_list:list = np.zeros(shape=(self.sim_steps,self.N+1),dtype=float)
        self.psi_list:list = np.zeros(shape=(self.sim_steps,self.N),dtype=float)
        self.Area_list = np.zeros(self.N,dtype=float)

        self.r_unperturbed:list = []
        self.z_unperturbed:list = []
        self.psi_unperturbed:list = []     

        self.Potential_E = np.zeros(self.sim_steps-1,dtype=float)
        self.Potential_E_before_correction = np.zeros(self.sim_steps-1,dtype=float)
        self.Kinetic_E = np.zeros(self.sim_steps-1,dtype=float)
        self.Xsqre = np.zeros(self.sim_steps-1,dtype=float)
        self.correct_count_list = np.zeros(self.sim_steps-1)


        # Setting up program
        self.var_corr_tol:float = 1e-5 # The tolerence for when to use variable correction
        self.print_scale:int = int((self.sim_steps-2)/1000)

        self.integration_method:str = "RK4" #Type of integration scheme
        self.save_path:str
        self.df_name:str
        self.do_correction:bool = True
        self.print_progress:bool = True
        self.start_flat:bool = True



    def run_simulation(self):
        print("Running" +f"print(N,T,dt)={print(self.N,self.T,self.dt)}")



if __name__ == "__main__":
    membrane = Surface_membrane()

    membrane.run_simulation()