load('p38_NFkB_Trajectories_all.mat')

rescale_data_or_not = 1;
fig_save_path = './SubFigures2023/';

%% rescaling the data and save
if rescale_data_or_not
    data_LPS = AllTraj.nfkb_NonSmoothedTrajectories.LPS;
    
    % the macrophage cell volume is 6 pl, measured by Stefanie Luecke
    NFkB_max_range = [0.25, 1.75]/6;
    % We assume cells in response to LPS 100ng can reach
    % the highest Nucleus NFkB concentration (all the NFkB can enter into nucleus).
    % (Pam3CSK might not be able to reach the highest Nuc NFkB conc)
    rescale_factor_NFkB = data_rescaling_factor_2023_p38_data(NFkB_max_range,data_LPS);
    rescale_factor_p38 = 0.08/0.73;
else
    rescale_factor_NFkB = 1;
    rescale_factor_p38 = 1;
end
NFkB_exp_data = cell(4,1);
p38_exp_data = cell(4,1);
corr_mat = struct();

for i_sti = 1:4
    NFkB_exp_data{i_sti} = AllTraj.nfkb_SmoothedTrajectories.(AllTraj.StimulusOrder{i_sti}) * rescale_factor_NFkB;
    NFkB_exp_data{i_sti} = NFkB_exp_data{i_sti}(:,13:end);
    p38_exp_data{i_sti} = AllTraj.p38_SmoothedTrajectories.(AllTraj.StimulusOrder{i_sti}) * rescale_factor_p38;
    p38_exp_data{i_sti} = p38_exp_data{i_sti}(:,13:end);
    
    [corr_mat.exp.p38{i_sti}] = func_metrics(p38_exp_data{i_sti},0.01);
    [corr_mat.exp.NFkB{i_sti}] = func_metrics(NFkB_exp_data{i_sti},0.05);
    
end
exp_info_ligand = AllTraj.StimulusOrder;



data_filename = 'Sim_NFKBfitting_alldose_p38_sampling_r2.mat';

% load data, sim_data_tbl
% data: nucNFkB, ppp38, IKK, TAK1, ppMKK4, ppMKK6
data_save_file_path_1 = './raw_data/';%_fay_parameter/';
vis_data_field = {'sampling'};

load(strcat(data_save_file_path_1,data_filename))
data_index = [1;2;3;4];
data_to_draw.(vis_data_field{1}) = data.model_sim(data_index);


NFkB_sample_data = cell(4,1);
p38_sample_data = cell(4,1);

for i_data = 1:length(data_to_draw.(vis_data_field{1}))
    
    NFkB_sample_data{i_data} = data_to_draw.(vis_data_field{1}){i_data}(1:6:end,:);
    p38_sample_data{i_data} = data_to_draw.(vis_data_field{1}){i_data}(2:6:end,:);
    
    [corr_mat.sample.p38{i_data}] = func_metrics(p38_sample_data{i_data},0.01);
    [corr_mat.sample.NFkB{i_data}] = func_metrics(NFkB_sample_data{i_data},0.05);
    
end

sample_info_ligand = data.info_ligand(data_index);


spearman_corr_sample = 2 * ones(4,3);
spearman_corr_rho_sample = -ones(4,3);

spearman_corr_exp = 2 * ones(4,3);
spearman_corr_rho_exp = -ones(4,3);

metric_fields = {'time2halfMax','peakVec','totalIntegral'};

for i_metric_fields = 1:length(metric_fields)
    for i_sti = 1:4
        NFkB_metrics = corr_mat.exp.NFkB{i_sti}.(metric_fields{i_metric_fields});
        p38_metrics = corr_mat.exp.p38{i_sti}.(metric_fields{i_metric_fields});
        [spearman_corr_exp(i_sti,i_metric_fields), spearman_corr_rho_exp(i_sti,i_metric_fields)] = corr(NFkB_metrics,p38_metrics,'Type','Spearman');
        
        if 0
            figure(1)
            paperpos = [0,0,100,80];
            papersize = [100,80];
            draw_pos = [20,20,40,20];
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            scatter(NFkB_metrics,p38_metrics,4)
            title({metric_fields{i_metric_fields}})
            xlable('NFkB')
            
        end
    end
end

for i_metric_fields = 1:length(metric_fields)
    for i_sti = 1:4
        NFkB_metrics = corr_mat.sample.NFkB{i_sti}.(metric_fields{i_metric_fields});
        p38_metrics = corr_mat.sample.p38{i_sti}.(metric_fields{i_metric_fields});
        [spearman_corr_sample(i_sti,i_metric_fields), spearman_corr_rho_sample(i_sti,i_metric_fields)] = corr(NFkB_metrics,p38_metrics,'Type','Spearman');
    end
end
for i_sti = 1:4
    
    for i_hour = 1:10
        NFkB_metrics = corr_mat.exp.NFkB{i_sti}.intigral_half_hour(:,i_hour);
        p38_metrics = corr_mat.exp.p38{i_sti}.intigral_half_hour(:,i_hour);
        [spearman_corr_int_hour(i_sti,i_hour), spearman_corr_rho_int_hour(i_sti,i_hour)] = corr(NFkB_metrics,p38_metrics,'Type','Spearman');
        NFkB_metrics = corr_mat.sample.NFkB{i_sti}.intigral_half_hour(:,i_hour);
        p38_metrics = corr_mat.sample.p38{i_sti}.intigral_half_hour(:,i_hour);
        [spearman_corr_int_hour_sample(i_sti,i_hour), spearman_corr_rho_int_hour_sample(i_sti,i_hour)] = corr(NFkB_metrics,p38_metrics,'Type','Spearman');
    end
end


index_sample = [4,3,2,1];
index_exp = [3,4,1,2];

spearman_corr_sample = spearman_corr_sample(index_sample,:);
spearman_corr_rho_sample = spearman_corr_rho_sample(index_sample,:);

spearman_corr_exp = spearman_corr_exp(index_exp,:);
spearman_corr_rho_exp = spearman_corr_rho_exp(index_exp,:);

spearman_corr_int_hour = spearman_corr_int_hour(index_exp,:);
spearman_corr_rho_int_hour = spearman_corr_rho_int_hour(index_exp,:);

spearman_corr_int_hour_sample = spearman_corr_int_hour_sample(index_sample,:);
spearman_corr_rho_int_hour_sample = spearman_corr_rho_int_hour_sample(index_sample,:);

sti_vec = {'Pam3CSK','CpG','LPS','TNF'};
mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
    ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

metric_fields = {'Time2HalfMax','MaxAmp','TotalIntegral'};

if 1
    
    %LPS only version
    figure(1)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,100,50]*2;
    paper_size = [100,50]*2;
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    h = heatmap( metric_fields,sti_vec(3),spearman_corr_sample(3,:),'Colormap',mymap,'CellLabelColor','k');%'none'
    % h.Title = strcat('');
    caxis([-1,1])
    h.FontColor = [0,0,0];
    
    %         for i_index = 1:length(index)
    %             h.XDisplayLabels{i_index} = ['\color[rgb]{0.8,0.8,0.8}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8} {red}
    %         end
    %
    %         for i_index = 1:length(index_non_wide)
    %             h.XDisplayLabels{i_index} = ['\color[rgb]{0.4,0.4,0.4}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
    %         end
    h.CellLabelColor = [0,0,0];
    h.CellLabelFormat = '%0.2g';
    saveas(gcf,strcat(fig_save_path,'Figure6G_corr_sample_Data_LPS'),'epsc')
    %saveas(gcf,strcat(fig_save_path,'corr_sample_Data_LPS'),'svg')
    
    close()
    
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    paper_pos = [0,0,100,50]*2;
    paper_size = [100,50]*2;
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    h = heatmap( metric_fields,sti_vec(3),spearman_corr_exp(3,:),'Colormap',mymap,'CellLabelColor','k');%'none'
    % h.Title = strcat('');
    caxis([-1,1])
    h.FontColor = [0,0,0];
    
    %         for i_index = 1:length(index)
    %             h.XDisplayLabels{i_index} = ['\color[rgb]{0.8,0.8,0.8}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8} {red}
    %         end
    %
    %         for i_index = 1:length(index_non_wide)
    %             h.XDisplayLabels{i_index} = ['\color[rgb]{0.4,0.4,0.4}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
    %         end
    h.CellLabelColor = [0,0,0];
    h.CellLabelFormat = '%0.2g';
    saveas(gcf,strcat(fig_save_path,'Figure6G_corr_exp_Data_LPS'),'epsc')
    %saveas(gcf,strcat(fig_save_path,'corr_exp_Data_LPS'),'svg')
    
    close()
end





function corr_struc = func_metrics(data,duration_thresh)
corr_struc = struct();

corr_struc.peakVec = max(data,[],2);
corr_struc.time2halfMax= -ones(size(data,1),1);

for i_cell = 1:size(data,1)
    if corr_struc.peakVec(i_cell)<= 0
        corr_struc.time2halfMax(i_cell) = -1;
    else
        corr_struc.time2halfMax(i_cell) = (find(data(i_cell, :) >= corr_struc.peakVec(i_cell)/2,1,'first')-1)*1/12;
    end
end
corr_struc.totalIntegral = sum(data,2)*1/12; % units in hours

if nargin ==2
    corr_struc.duration = sum(data >= duration_thresh,2) /12;
end

for i_hour = 1:10
    corr_struc.intigral_half_hour(:,i_hour) = sum(data(:,(i_hour-1)*6+1:i_hour*6),2)*1/12;
end

end