% run_me
clear all
Figure6CD = 1;
Figure6EF = 1;% run this
Figure6G = 1;

run_rescale_data = 0;
run_sampling = 0;


data_save_file_path = './raw_data/';%_fay_parameter/';
fig_save_path = './SubFigures2023/';

addpath('./src/')
addpath(data_save_file_path)
addpath(fig_save_path)



if Figure6CD
    %sim_example_
    % sim_example_p38
    sim_representative_p38
end

if Figure6G
    cal_plot_corr_NFkB_p38
end

if run_rescale_data
    rescale_data
end

if Figure6EF
    draw_heatmap_experiments % to be finished
    draw_traj_heatmap_diff_order_2023
end


if run_sampling
    input_paras.Num_sample = 10;
    
    % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
    input_paras.proj_num_vec = {[2];
        [3];
        [4];
        [6]};
    input_paras.Num_sample = [500;
        500;
        500;
        500];
    
    input_paras.proj_ligand_vec = {{'TNF'};
        {'LPS'};
        {'CpG'};
        {'Pam3CSK'}};
    input_paras.proj_dose_str_vec = {{'100ng/mL'};
        {'100ng/mL'};
        {'1000nM'};
        {'100ng/mL'}};
    input_paras.proj_dose_val_vec = {{100};
        {100};
        {1000};
        {100}};

    
    cal_codon =1;
    
    save_metric_name = 'Sim_NFKBfitting_alldose_p38_sampling_r2.mat';%r2
    % data_save_file_path
    para_sample_fitting_sim_sti_doses_var_2023(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
    
end

