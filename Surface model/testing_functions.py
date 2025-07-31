from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_simulation_function import Two_D_simulation
from Make_movie import Make_frames, Make_video
from two_d_data_processing import check_area
from Two_D_functions import Langrange_multi, Epsilon_values


def test_Lagrange_multi():
    const_args = Two_D_Constants(
        print_val=True
        ,init_rand_psi=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    make_movies = False
    """lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi_list[0]
        ,radi=radi_list[0]
        ,z_list=z_list[0]
        ,print_matrix=True
    )"""
    t = 0
    lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi_list[t]
        ,radi=radi_list[t]
        ,z_list=z_list[t]
            )


def test_make_frames():
    const_args = Two_D_Constants(
        print_val=True
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]


    Make_frames(
        data_path=data_path
        ,figs_save_path=video_fig_path
        ,df_name=df_name
    )

def test_make_video():
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    Make_video(
        output_path = video_save_path
        ,input_path = video_fig_path
        ,video_name = df_name
        ,fps=12
    )


def test_check_area():
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    #sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    Make_video(
        output_path = video_save_path
        ,input_path = video_fig_path
        ,video_name = df_name
        ,fps=12
    )
        

    check_area(
        df_name=df_name
        ,data_path=data_path
    )
    

def test_epsilon_value():
    const_args = Two_D_Constants(
        print_val=True
        ,init_rand_psi=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    print(f"r={radi_list[0]}")
    print(f"Area={Area}")
    ef,eg = Epsilon_values(
        N=N,r=radi_list[0],z=z_list[0],psi=psi_list[0],Area=Area
        ,print_matrix=True
    )

    print(f"ef={ef}")
    print(f"eg={eg}")

if __name__ == "__main__":
    #test_Lagrange_multi()
    #test_make_frames()
    #test_make_video()
    #test_check_area()

    test_epsilon_value()
