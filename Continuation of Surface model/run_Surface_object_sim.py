from Surface_class import Surface_membrane, multi_process,plotting_multi_process_results









if __name__ == "__main__":
    
    N = 80
    dt = 1e-13
    T = 0.5e-8
    save_path = f"2D sim results\\Testing with zN subtraction\\(T,N,dt)=({T:.1e},{N},{dt:0.1e})\\"

    membrane = Surface_membrane(
        T = T
        ,dt = dt
        ,const_index = 1
        ,N = N
        ,save_path = save_path
        )
    membrane.var_corr_tol = 5e-4
    membrane.integration_method = "RK4"
    membrane.make_movie = True
    membrane.make_plots = True
    
    membrane.run_sim()