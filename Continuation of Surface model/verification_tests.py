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

    tau_list = [tau ,0 ,0 ,tau ,tau]
    c0_list = [0 ,0 ,c0 ,0 ,0]
    perturb_bool_list = [False ,False ,False ,True ,True ]
    perturb_psi_val_list = [0 ,0 ,0 ,dpsi ,-dpsi]
    sim_index = 0

    process = []
    for i in range(len(tau_list)):
        membrane = Surface_membrane(
            N = 40
            ,T = 5e-7
            ,dt = 1e-11
            ,const_index=sim_index
            ,dpsi = perturb_psi_val_list[i]
            ,perturb_var = "psi"
            ,save_path= path + f"(tau,c0,dpsi) = ({tau_list[i]:.1e}  ,{c0_list[i]:.1e}  ,{perturb_psi_val_list[i]:.1e})\\"
        )

        membrane.perturb = perturb_bool_list[i]
        membrane.print_constants = False
        membrane.start_flat = True
        membrane.close_final_plots = True
        membrane.var_corr_tol = 1e-4
        #membrane.init_config_show_time = 400
        
        process.append(multiprocessing.Process(
            #target=membrane.dynamics#run_sim
            target=membrane.run_sim
        ))

    for p in process:
        p.start()
    
    for p in process:
        p.join()


    plotting_multi_process_results(
        path = path
    )


if __name__ == "__main__":
    verification_of_model()