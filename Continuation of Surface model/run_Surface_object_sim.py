from Surface_class import Surface_membrane, multi_process,plotting_multi_process_results

if __name__ == "__main__":
    
    N = 30
    save_path = f"2D sim results\\testing plus with Euler\\"
    membrane = Surface_membrane(
        T=1e-7
        ,dt=1e-11
        ,const_index = 1
        ,N=N
        ,save_path= save_path
        )
    membrane.var_corr_tol = 1e-3
    membrane.integration_method = "Euler"
    membrane.make_movie = True
    membrane.make_plots = True
    
    membrane.run_sim()