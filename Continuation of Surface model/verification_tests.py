from Surface_class import Surface_membrane, plotting_multi_process_results
import multiprocessing

def verification_of_model():

    membrane = Surface_membrane(
        const_index=0
    )
    
    tau = membrane.tau
    c0 = membrane.c0
    dpsi = membrane.dpsi_perturb

    path = "2D sim results\\model verification\\"

    tau_list = [tau ,0 ,0  ,tau ,tau]
    c0_list =  [ 0  ,0 ,c0 ,0   ,0]
    perturb_bool_list = [False ,False ,False ,True ,True ]
    perturb_psi_val_list = [0 ,0 ,0 ,dpsi ,-dpsi]
    sim_index = 0

    process = []
    for i in range(len(tau_list)):
        membrane = Surface_membrane(
            N = 30
            ,T = 5e-8
            ,dt = 2.5e-11
            ,const_index=sim_index
            ,dpsi = perturb_psi_val_list[i]
            ,perturb_var = "psi"
            ,save_path= path + f"(tau,c0,dpsi) = ({tau_list[i]:.1e}  ,{c0_list[i]:.1e}  ,{perturb_psi_val_list[i]:.1e})\\"
        )
        membrane.show_stationary_state = False
        membrane.start_flat = True
        membrane.const_length_diff_N_density = False
        membrane.perturb = perturb_bool_list[i]
        membrane.print_constants = False
        
        membrane.close_final_plots = True
        membrane.var_corr_tol = 1e-4


        membrane.setup_simulation()

        membrane.tau = tau_list[i]
        membrane.c0 = c0_list[i]

        
        process.append(multiprocessing.Process(
            #target=membrane.dynamics#run_sim
            target=membrane.dynamics
        ))

    for p in process:
        p.start()
    
    for p in process:
        p.join()


    plotting_multi_process_results(
        path = path
        ,make_movie=True
        ,make_figures=True
        ,make_comparison_figs=False
    )



def rolling_test():
    path = "2D sim results\\model verification\\rolling test\\"
    integration_method = "Euler"
    membrane = Surface_membrane(
            N = 40
            ,T = 1e-6
            ,dt = 2.5e-11
            ,const_index=0
            ,dpsi = 0
            ,perturb_var = "psi"
            ,save_path= path + integration_method + f"(tau,c0) = ({0:.1e}  ,{25:.1e} )\\"
        )

    membrane.integration_method = integration_method
    membrane.show_stationary_state = True
    membrane.start_flat = True
    membrane.const_length_diff_N_density = False

    membrane.close_final_plots = True
    membrane.var_corr_tol = 1e-3

    membrane.setup_simulation()

    membrane.tau = 0
    #membrane.c0 = c0_list[i]
    
    membrane.make_movie = True
    membrane.make_plots = True

    membrane.print_consts()
    membrane.dynamics()
    membrane.plotting_n_movie_data()



if __name__ == "__main__":
    #verification_of_model()
    rolling_test()