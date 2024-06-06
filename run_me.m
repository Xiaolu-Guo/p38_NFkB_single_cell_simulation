% run_me
clear all
Figure6CD = 0;
Figure6EF = 0;% run this
Figure6G = 0;
Figure6_denoise = 0;
sens_analysis = 0;

run_rescale_data = 0;
run_sampling = 0;

data_save_file_path = './raw_data/';%_fay_parameter/';
fig_save_path = './SubFigures2023/';

addpath('./src/')
addpath(data_save_file_path)
addpath(fig_save_path)


%% Figure6CD
if Figure6CD
    %sim_example_
    % sim_example_p38
    sim_representative_p38
end

%% Figure6G
if Figure6G
    cal_plot_corr_NFkB_p38
end

%% run_rescale_data
if run_rescale_data
    rescale_data
end

%% Figure6EF
if Figure6EF
    draw_heatmap_experiments % to be finished
    draw_traj_heatmap_diff_order_2023
end

%% run_sampling
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

%% run sensitivity analysis
if sens_analysis 
    
    run_sens_for_revision = 1;
    if run_sens_for_revision
        
        % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
        input_paras.proj_num_vec = {[3]};
        input_paras.Num_sample = [21];%in total 21 parameters
        
        input_paras.proj_ligand_vec = {{'LPS'}};
        input_paras.proj_dose_str_vec = {{'100ng/mL'}};
        input_paras.proj_dose_val_vec = {{100}};
        
        cal_codon =1;
        MKK4_fold_vec = [10.^(linspace(-1,1,21));
            ones(1,21)];
        MKK6_fold_vec = MKK4_fold_vec;
        p38_fold_vec = MKK4_fold_vec;
        sigma_string_vec = {'Changed','Fixed'};
        for i_MKK4 = 1:size(MKK4_fold_vec,1)
            for i_MKK6 = 1:size(MKK6_fold_vec,1)
                for i_p38 = 1:size(p38_fold_vec,1)
                    
                    p38_model_fold.MKK4 = MKK4_fold_vec(i_MKK4,:);
                    p38_model_fold.MKK6= MKK6_fold_vec(i_MKK6,:);
                    p38_model_fold.p38= p38_fold_vec(i_p38,:);
                    
                    save_metric_name = strcat( 'Sim_p38_sens_MKK4S',sigma_string_vec{i_MKK4},...
                        '_MKK6S',sigma_string_vec{i_MKK6},...
                        '_p38S',sigma_string_vec{i_p38},...
                        '_r2.mat');
                    para_sens_2023_12_forrevision(p38_model_fold,data_save_file_path,input_paras,save_metric_name)
                end
            end
        end
    end
    
    %% run sensitivity analysis for revision, sample NFkB, sens p38
run_sample_NFKB_sens_p38_for_revision = 1;
if run_sample_NFKB_sens_p38_for_revision
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    
    % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
    %     input_paras.proj_num_vec = {[3]};
    %     input_paras.Num_sample = [500];%in total 21 parameters
    %
    %     input_paras.proj_ligand_vec = {{'LPS'}};
    %     input_paras.proj_dose_str_vec = {{'100ng/mL'}};
    %     input_paras.proj_dose_val_vec = {{100}};
    
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
    
    cal_codon = 0;
    MKK4_mean_fold_vec = [0.5,1,2];
    MKK6_mean_fold_vec = [0.5,1,2];
    p38_mean_fold_vec = [0.5,1,2];
    sigma_string_vec = {'p5fold','1fold','2fold'};
    for i_MKK4 = 1:length(MKK4_mean_fold_vec)
        for i_MKK6 = 1:length(MKK6_mean_fold_vec)
            for i_p38 = 1:length(p38_mean_fold_vec)
                
                p38_model.mean_fold.MKK4 = MKK4_mean_fold_vec(i_MKK4);
                p38_model.mean_fold.MKK6= MKK6_mean_fold_vec(i_MKK6);
                p38_model.mean_fold.p38= p38_mean_fold_vec(i_p38);
                
                save_metric_name = strcat( 'Sim_NFKB_sample_p38_sens_MKK4S',sigma_string_vec{i_MKK4},...
                    '_MKK6S',sigma_string_vec{i_MKK6},...
                    '_p38S',sigma_string_vec{i_p38},...
                    '_r1.mat');
                para_sample_fitting_sim_sti_doses_var_2023_12_forrevision(p38_model,data_save_file_path,input_paras,save_metric_name)
                
            end
        end
    end
end

run_sample_NFKB_sens_p38_for_revision_server = 1;
if run_sample_NFKB_sens_p38_for_revision_server
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    
    % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
    %     input_paras.proj_num_vec = {[3]};
    %     input_paras.Num_sample = [500];%in total 21 parameters
    %
    %     input_paras.proj_ligand_vec = {{'LPS'}};
    %     input_paras.proj_dose_str_vec = {{'100ng/mL'}};
    %     input_paras.proj_dose_val_vec = {{100}};
    
    input_paras.proj_num_vec = {[3]};
    input_paras.Num_sample = [1000];
    
    input_paras.proj_ligand_vec = {{'LPS'}};
    input_paras.proj_dose_str_vec = {{'100ng/mL'}};
    input_paras.proj_dose_val_vec = {{100}};
    
    cal_codon = 0;
    MKK4_mean_fold_vec = [0.1,1,10];
    MKK6_mean_fold_vec = [0.1,1,10];
    p38_mean_fold_vec = [0.5,1,2];
    sigma_string_vec = {'p1fold','1fold','10fold'};
    p38_model = cell(27,1);
    i_k = 1;
    for i_MKK4 = 1:length(MKK4_mean_fold_vec)
        for i_MKK6 = 1:length(MKK6_mean_fold_vec)
            for i_p38 = 1:length(p38_mean_fold_vec)
                i_MKK4_ind(i_k) = i_MKK4;
                i_MKK6_ind(i_k) = i_MKK6;
                i_p38_ind(i_k) = i_p38;
                i_k = i_k+1;
                
                
            end
        end
    end
    
    ncpu=8;
    pc=parcluster('local');
    pc.NumThreads=2;%
    parpool(pc,ncpu)
    parfor i_k = 1:length(p38_model)
        p38_model{i_k}.mean_fold.MKK4 = MKK4_mean_fold_vec(i_MKK4_ind(i_k));
        p38_model{i_k}.mean_fold.MKK6= MKK6_mean_fold_vec(i_MKK6_ind(i_k));
        p38_model{i_k}.mean_fold.p38= p38_mean_fold_vec(i_p38_ind(i_k));
        
        save_metric_name = strcat( 'Sim_NFKB_sample_p38_large_sens_MKK4S',sigma_string_vec{i_MKK4_ind(i_k)},...
            '_MKK6S',sigma_string_vec{i_MKK6_ind(i_k)},...
            '_p38S',sigma_string_vec{i_p38_ind(i_k)},...
            '_r2.mat');
        para_sample_fitting_sim_sti_doses_var_2023_12_forrevision(p38_model{i_k},data_save_file_path,input_paras,save_metric_name)
        
    end
    delete(gcp)
    
end
    
    %% draw sensitivity analysis
    if 1 % draw the curves
        draw_sens_p38()
    end
    
    if 1 % draw the heatmap
        draw_sens_withnoise()
    end
    
end



%% add noise: v5 02262024
if 0
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    input_paras.Num_sample = 10;
    
    % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
    input_paras.proj_num_vec = {[3]};
    input_paras.Num_sample = [1000];
    
    input_paras.proj_ligand_vec = {{'LPS'}};
    input_paras.proj_dose_str_vec = {{'100ng/mL'}};
    input_paras.proj_dose_val_vec = {{100}};
    
    cal_codon =1;
    MKK4_sigma_vec = [0,log(2)/3,0,log(2)/3,log(2)/3,log(2)/3];
    MKK6_sigma_vec = [0,log(2)/3,0,log(2)/3,log(2)/3,log(2)/3];
    p38_sigma_vec = [0,log(2)/3,0,log(2)/3,log(2)/3,log(2)/3];
    sigma_string_vec = {'p38m_None','p38m_Noise','p38m_None','p38m_Noise','p38m_Noise','p38m_Noise'};
    
    NFkB_string_vec = {'NFkBm_None','NFkBm_None','NFkBm_Noise','NFkBm_Noise','RCP_Cdeg_end_TAK_None','RCP_Cdeg_end_None'};
    paranames_vec = {{'params99','params101'};
        {'params99','params101'};
        {};
        {};
        {'params44','params40','params36','params52n2','params65n2'}
        {'params44','params40','params36'}};
    paravals_vec = {[0.4,0.4]';
        [0.4,0.4]';
        []';
        []';
        [0.012,0.065681,0.065681,1,1]'
        [0.012,0.065681,0.065681]'};
    
    specienames_vec = {{'NFkB'};
        {'NFkB'};
        {};
        {};
        {};
        {}};
    
    specievals_vec = {[0.08];
        [0.08];
        [];
        [];
        [];
        []};
    
    
    i_vers = 1;
    for ver_num = 1:10
        for i_cond = 1:length(MKK4_sigma_vec)
            vers{i_vers} = strcat('d0226_r',num2str(ver_num));
            i_cond_idx(i_vers) = i_cond;
            i_vers = i_vers +1;
        end
        
    end
    
    var_input = cell(1,length(vers));
    p38_model = cell(1,length(vers));
%     ncpu=2;
%     pc=parcluster('local');
%     pc.NumThreads=2;%
%     parpool(pc,ncpu)
    
    ncpu=8;
    pc=parcluster('local');
    pc.NumThreads=2;%
    parpool(pc,ncpu)
        
    parfor i_vers = 1:length(vers)% 1:length(MKK4_sigma_vec)
        
        % i_cond = i_cond_idx(i_vers);
        p38_model{i_vers}.sigma.MKK4 = MKK4_sigma_vec(i_cond_idx(i_vers));
        p38_model{i_vers}.sigma.MKK6= MKK6_sigma_vec(i_cond_idx(i_vers));
        p38_model{i_vers}.sigma.p38= p38_sigma_vec(i_cond_idx(i_vers));
        
        var_input{i_vers}.paranames = paranames_vec{i_cond_idx(i_vers)};%LPS
        var_input{i_vers}.paravals = paravals_vec{i_cond_idx(i_vers)};
        
        var_input{i_vers}.specienames = specienames_vec{i_cond_idx(i_vers)};
        var_input{i_vers}.specievals = specievals_vec{i_cond_idx(i_vers)};
        
        save_metric_name = strcat( 'Sim_NFKB_sample_',NFkB_string_vec{i_cond_idx(i_vers)},'_p38_sample_',sigma_string_vec{i_cond_idx(i_vers)},'_',vers{i_vers},'.mat');
        para_p38_sample_fitting_sim_alter_var_202402(p38_model{i_vers},var_input{i_vers},data_save_file_path,input_paras,save_metric_name)
        
    end
    delete(gcp)
    
end


%% run figure6 denoise
if Figure6_denoise
    
    if 0
        draw_p38_NFkB_denoise_3model_sim_int_0321(fig_save_path)
    end
    
    if 0 % features update figure 6F
        draw_p38_NFkB_denoise_3model_sim_features_Figure6F_0311(fig_save_path)
    end
    
    if 0
        draw_p38_NFkB_denoise_3model_sim_features_0228(fig_save_path)
        %0227
        %draw_p38_NFkB_denoise_features_0226(fig_save_path)
    end
    
end