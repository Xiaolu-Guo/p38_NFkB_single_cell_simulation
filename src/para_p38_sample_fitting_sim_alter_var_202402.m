% TNFo_dual_para_sample_main
% clear all
function [] = para_p38_sample_fitting_sim_alter_var_202402(p38_model,var_input,data_save_file_path,input_paras,save_metric_name)
data_info.save_file_path = data_save_file_path;

if nargin < 4
    cal_codon = 0;
end

if isfield(p38_model,'sigma')
    sigma_MKK4 = p38_model.sigma.MKK4;
    sigma_MKK6 = p38_model.sigma.MKK6;
    sigma_p38 = p38_model.sigma.p38;
else
    sigma_MKK4 = 0;
    sigma_MKK6 = 0;
    sigma_p38 = 0;
end

if isfield(p38_model,'mean_fold')
    mean_fold_MKK4 = p38_model.mean_fold.MKK4;
    mean_fold_MKK6 = p38_model.mean_fold.MKK6;
    mean_fold_p38 = p38_model.mean_fold.p38;
else
    mean_fold_MKK4 = 1;
    mean_fold_MKK6 = 1;
    mean_fold_p38 = 1;
end


paranames = var_input.paranames;
paravals = var_input.paravals;
specienames = var_input.specienames;
specievals = var_input.specievals;

proj_num_vec = input_paras.proj_num_vec;
proj_ligand_vec = input_paras.proj_ligand_vec;
proj_dose_str_vec = input_paras.proj_dose_str_vec;
proj_dose_val_vec = input_paras.proj_dose_val_vec;

if isfield(input_paras,'Num_sample')
    Num_sample = input_paras.Num_sample;
else
    Num_sample = 500;
end

if length(Num_sample) == length(proj_num_vec)
    Num_sample = Num_sample;
else
    Num_sample = Num_sample * ones(size(proj_num_vec));
end

if isfield(input_paras,'var_fold_mat')
    if length(input_paras.var_fold_mat) == length(proj_num_vec)
        var_fold_mat = input_paras.var_fold_mat;
    else
        var_fold_mat = input_paras.var_fold_vec * ones(size(proj_num_vec));
    end
else
    var_fold_mat = ones(size(proj_num_vec));
end


gene_info.gene_type = {'wt'};
gene_info.gene_parameter_value_vec_genotype = cell(0);

%% parameter from the data

if 1
    data_save_file_path_1 = '../raw_data/';%_fay_parameter/';
    
    load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
    
    data_fitting = data;
    clear data
    for i_data = 1:length(data_fitting.pred_mode)
        
        data_fitting.pred_mode_cv{i_data} = std(data_fitting.pred_mode_filter_nan{i_data},[],2)./mean(data_fitting.pred_mode_filter_nan{i_data},2);
        [~,data_fitting.pred_mode_cv_order{i_data}] = sort(data_fitting.pred_mode_cv{i_data},'descend');
        
    end
    
    % input
    thresh_TNF=0.33;
    thresh_field = {'pred_mode_cv_order'};%pred_mode_cv_order osc
    data_fitting.parameters_mode_filter_TNF = cell(2,1);
    for i_data = 1:3
        index = data_fitting.(thresh_field{1}){i_data}(1: ceil(length(data_fitting.(thresh_field{1}){i_data}) * thresh_TNF));
        data_fitting.parameters_mode_filter_TNF{i_data} = data_fitting.parameters_mode_nan{i_data}(index ,:);
    end
    
    
    for i_data = 4:length(data_fitting.parameters_mode_nan)
        index = data_fitting.(thresh_field{1}){i_data}(1: ceil(length(data_fitting.(thresh_field{1}){i_data}) * thresh_TNF));
        data_fitting.parameters_mode_filter_TNF{i_data} = data_fitting.parameters_mode_nan{i_data}(: ,:);
    end
    % change the parameters for TNF only for high CV
    data_fitting.parameters_mode_nan = data_fitting.parameters_mode_filter_TNF;
    
    data_field_names = fieldnames(data_fitting);
    
    data_index = [1:6,10:12,14:19];
    
    for i_data_field = 1:length(data_field_names)
        data_new.(data_field_names{i_data_field}) = data_fitting.(data_field_names{i_data_field}) (data_index);
    end
    
    data_fitting_all = data_fitting;
    data_fitting = data_new;
    
end


%%

old_str = {'p','p','p','p','p','p','p','p','p','p'};
new_str = {'params','params','params','params','params','params','params','params','params','params'};


for i_proj_num_vec = 1:length(proj_num_vec)
    if length(proj_num_vec{i_proj_num_vec})<1
        error('elements of proj_num_vec has to be non-empty!')
    end
    
    sim_info.ligand = proj_ligand_vec{i_proj_num_vec};
    sim_info.dose_str = proj_dose_str_vec{i_proj_num_vec};
    sim_info.dose_val = proj_dose_val_vec{i_proj_num_vec};
    
    
    var_fold = var_fold_mat(i_proj_num_vec);
    switch length(proj_num_vec{i_proj_num_vec})
        case 1
            
            % all dose parameters are used for sampling
            i_data = find(strcmp(data_fitting_all.info_ligand,sim_info.ligand));
            para_val_ele = [];
            for i_i_data= 1:length(i_data)
                para_val_ele = [para_val_ele;data_fitting_all.parameters_mode_nan{i_data(i_i_data)} ] ; % [cellnumber x parameter]
            end
            i_data = find(strcmp(data_fitting.info_ligand,sim_info.ligand),1);
            
            
            rpt_time = ceil(Num_sample(i_proj_num_vec)*5/size(para_val_ele,1));
            para_val_total = [];
            for i_rpt = 1:rpt_time
                para_val_total = [para_val_total;para_val_ele];
            end
            para_val = para_val_total(randperm(size(para_val_total,1),Num_sample(i_proj_num_vec)),:);
            est.name = cellfun(@replace,data_fitting.para_name{i_data},old_str(1:length(data_fitting.para_name{i_data})),new_str(1:length(data_fitting.para_name{i_data})),'UniformOutput',false);
            NFkB_index = find(strcmp(est.name,'NFkB_cyto_init'));
            shift_index = find(strcmp(est.name,'shift'));
            parameter_index = setdiff(setdiff(1:length(est.name),NFkB_index,'stable'),shift_index,'stable');
            
            % non_NFkB_index = setdiff(1:length(est.name),NFkB_index,'stable');
            gene_info.parameter_name_vec = {est.name(parameter_index)};
            gene_info.gene_parameter_value_vec_genotype{1}{1} = para_val(:,parameter_index)';
            gene_info.species_name_vec= {{'NFkB','MKK4_off','MKK6_off','p38_off'}};
            
            % sampling parameters
            rand_mat_MKK4 = randn(1,Num_sample(i_proj_num_vec));
            fold_change_log_MKK4 = log(10^0) + rand_mat_MKK4 * sigma_MKK4; % sigma_MKK4 = (log(2)/3)
            fold_change_MKK4 = exp(fold_change_log_MKK4)*mean_fold_MKK4;
            
            rand_mat_MKK6 = randn(1,Num_sample(i_proj_num_vec));
            fold_change_log_MKK6 = log(10^0) + rand_mat_MKK6 * sigma_MKK6;
            fold_change_MKK6 = exp(fold_change_log_MKK6)*mean_fold_MKK6;
            
            rand_mat_p38 = randn(1,Num_sample(i_proj_num_vec));
            fold_change_log_p38 = log(10^0) + rand_mat_p38 * sigma_p38;
            fold_change_p38 = exp(fold_change_log_p38)*mean_fold_p38;
            
            gene_info.species_value_vec_genotype{1}{1} = [para_val(:,NFkB_index)';
                0.1*fold_change_MKK4;
                0.1*fold_change_MKK6;
                0.1*fold_change_p38];
            
        if ~isempty(paranames)
        gene_info.parameter_name_vec{1} = {gene_info.parameter_name_vec{1}{:},paranames{1:end}};
        gene_info.gene_parameter_value_vec_genotype{1}{1} = [gene_info.gene_parameter_value_vec_genotype{1}{:};
            paravals*ones(1,size(para_val,1))];
        end
        
        if ~isempty(specienames)
        gene_info.species_name_vec{1}= {gene_info.species_name_vec{1}{:},specienames{1:end}};
        gene_info.species_value_vec_genotype{1}{1} = [gene_info.species_value_vec_genotype{1}{:};
            specievals*ones(1,size(para_val,1))];
        end
            
        otherwise
            error('only one ligand stim supported!')
    end
    
    % stimili info
    
    save_filename = strcat('202307_para_sampled_fitting_NFkBinit_p38_',replace(num2str(var_fold),'.','p'),'varfoldchange');
    
    for i_ligand = 1:length(sim_info.ligand)
        save_filename = strcat(save_filename,'_',sim_info.ligand{i_ligand},...
            '_',replace(replace(sim_info.dose_str{i_ligand},'/',''),'.','p'));
    end
    
    % species that will be saved
    % must be r x 1, for each cell i must be ri x 1
    data_info.species_outputname = {'nucNFkB';'ppp38';'IKK';'TAK1';'ppMKK4';'ppMKK6'};
    data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'ppp38'};{'IKK'};{'TAK1'};{'ppMKK4'};{'ppMKK6'}};
    data_info.save_file_name = save_filename; % only beginning, no .mat
    
    sim_data_tbl = p38_sim_save_2023(sim_info,data_info,gene_info);
    
        
        i_ligand = 1;
        ligand_str= proj_ligand_vec{i_proj_num_vec}{i_ligand};
        dose_str = proj_dose_str_vec{i_proj_num_vec}{i_ligand};
        
        for i_ligand = 2:length(proj_ligand_vec{i_proj_num_vec})
            ligand_str = strcat(ligand_str,'_',proj_ligand_vec{i_proj_num_vec}{i_ligand});
            dose_str = strcat(dose_str,'_',proj_ligand_vec{i_proj_num_vec}{i_ligand});
        end
        
        data.info_ligand{i_proj_num_vec} = ligand_str;
        data.model_sim{i_proj_num_vec} = sim_data_tbl.trajectory(:,1:5:end);
        data.info_dose_index{i_proj_num_vec} = 1;
        data.info_dose_str{i_proj_num_vec} = dose_str;
        data.info_num_cells{i_proj_num_vec} = size(sim_data_tbl.trajectory,1);
        data.para_names{i_proj_num_vec} = sim_data_tbl.parameter_name;
        data.para_values{i_proj_num_vec} = sim_data_tbl.parameter_value;
        data.init_species_name{i_proj_num_vec} = sim_data_tbl.init_species_name;
        data.init_species_values{i_proj_num_vec} = sim_data_tbl.init_species_value;
        [~, data.order{i_proj_num_vec}] = sort(max(sim_data_tbl.trajectory,[],2),'descend');
        
    % save(save_metric_name,'data');
end

save(strcat(data_save_file_path,save_metric_name),'data');

