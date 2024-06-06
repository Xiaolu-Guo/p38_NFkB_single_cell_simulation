function [] = draw_p38_NFkB_denoise_3model_sim_int_0321(fig_save_path)

save_metric_name = {'Sim_NFKB_sample_NFkBm_None_p38_sample_p38m_None_d0226_r';
    'Sim_NFKB_sample_NFkBm_None_p38_sample_p38m_Noise_d0226_r';
    'Sim_NFKB_sample_NFkBm_Noise_p38_sample_p38m_None_d0226_r';
    'Sim_NFKB_sample_NFkBm_Noise_p38_sample_p38m_Noise_d0226_r';
    'Sim_NFKB_sample_RCP_Cdeg_end_TAK_None_p38_sample_p38m_Noise_d0226_r';
    'Sim_NFKB_sample_RCP_Cdeg_end_None_p38_sample_p38m_Noise_d0226_r'};

denoise_label = {'NFkBm-None-p38m-None'
    'NFkBm-None';
    'p38m-None';
    'control';
    'RCPall-TAK-none';
    'RCPall-none'}';


%% features

if 1
    int2hr_label = {'1st-2hr','2nd-2hr','3rd-2hr','4th-2hr'};
    clear duration_all integral_all peak_all time2half integral_by_halfhour
    clear duration_all_NFkB integral_all_NFkB peak_all_NFkB time2half_NFkB integral_by_halfhour_NFkB
    
    
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    
    rho = [];
    pval = [];
    p38_CV = [];
    NFkB_CV = [];
    p38_Fano = [];
    NFkB_Fano = [];
    NFkB_BC = [];
    p38_BC = [];
    
    
    int_num = length(int2hr_label);
    for rep = 1 :10
        for i_sig = 1: 4%length(save_metric_name)
            
            load(strcat('../raw_data/',save_metric_name{i_sig},num2str(rep),'.mat'))
            
            %             codon_struc = func_metrics(data.model_sim{1, 1}(2:6:size(data.model_sim{1, 1},1),:),0.01,24);
            %             NFkB_codon_struc = func_metrics(data.model_sim{1, 1}(1:6:size(data.model_sim{1, 1},1),:),0.05,24);
            
            
            clear traj
            traj.NFkB = data.model_sim{1}(1:6:end,:);
            traj.p38 = data.model_sim{1}(2:6:end,:);
            
            codon_struc_NFkB = func_metrics(traj.NFkB,0.05,6);
            codon_struc_p38 = func_metrics(traj.p38,0.01,6);
            
            integral_by_halfhour{i_sig} = codon_struc_p38.intigral_half_hour;
            integral_by_halfhour_NFkB{i_sig} = codon_struc_NFkB.intigral_half_hour;

            
            int_num = 4;
            for i_int = 1:int_num
                
                
                if 0
                    figure(i_int)
                    Set_figure_size_wide_short
                    figure(i_int)
                    scatter(integral_by_halfhour{i_sig}(:,i_int),integral_by_halfhour_NFkB{i_sig}(:,i_int),4,[0,0,1],'filled');hold on
                    xlim([0,0.2]/4)
                    ylim([0,1]/4);
                    figure(i_int)
                    saveas(gcf,strcat(fig_save_path,'R1P2_',save_metric_name{i_sig},'_int_',num2str(i_int),'_ranged' ),'epsc')
                    close
                end
                [rho(i_sig,i_int,rep), pval(i_sig,i_int,rep)] = corr(integral_by_halfhour{i_sig}(:,i_int), integral_by_halfhour_NFkB{i_sig}(:,i_int), 'Type', 'Spearman');
                rho(i_sig,i_int,rep)
                % scatter(integral_by_halfhour{i_sig}(:,i_int),integral_by_halfhour_NFkB{i_sig}(:,i_int))
                
                % Calculate the mean and standard deviation for each vector
                vector1 = integral_by_halfhour{i_sig}(:,i_int);
                mean1 = mean(vector1);
                stdDev1 = std(vector1);
                p38_CV(i_sig,i_int,rep) = stdDev1 / mean1;
                p38_Fano(i_sig,i_int,rep) = stdDev1^2 / mean1;
                
                
                vector2 = integral_by_halfhour_NFkB{i_sig}(:,i_int);
                mean2 = mean(vector2);
                stdDev2 = std(vector2);
                NFkB_CV(i_sig,i_int,rep) = stdDev2 / mean2;
                NFkB_Fano(i_sig,i_int,rep) = stdDev2^2 / mean2;
                
                [~, p38_BC(i_sig,i_int,rep)] = bimodalitycoeff(vector1);
                [~, NFkB_BC(i_sig,i_int,rep)] = bimodalitycoeff(vector2);
                
            end
            
        end
        
        % save('NFkB_p38_feature_10rep.mat','rho','pval','p38_CV','NFkB_CV','denoise_label','NFkB_feature_list','feature_idx')
        
    end
    
    
end


%% CV
if 1
    
    if 1
        for rep =1:3
            for i_int = 1:int_num% size(p38_CV,1)
                figure(i_int)
                Set_figure_size_wide_short
                scatter(1:size(p38_CV,1),p38_CV(4:-1:1,i_int,rep),10,[40/255,164/255,255/255],'filled');hold on
                scatter(1:size(NFkB_CV,1),NFkB_CV(4:-1:1,i_int,rep),10,[255/255,187/255,40/255],'filled');hold on
                xlim([0.5,4.5])
                ylim([0,2])
                xticks(1:4)                
                saveas(gcf,strcat(fig_save_path,'R1P2_int',num2str(i_int),'_rep',num2str(rep),'_features_CV' ),'epsc')
                saveas(gcf,strcat(fig_save_path,'R1P2_int',num2str(i_int),'_rep',num2str(rep),'_features_CV' ),'svg')                
                close
                
                
            end
            
            for i_int = 1:int_num% size(p38_CV,1)
                figure(i_int)
                Set_figure_size_wide_short
                scatter(1:size(p38_BC,1),p38_BC(4:-1:1,i_int,rep),10,[40/255,164/255,255/255],'filled');hold on
                scatter(1:size(NFkB_BC,1),NFkB_BC(4:-1:1,i_int,rep),10,[255/255,187/255,40/255],'filled');hold on
                xlim([0.5,4.5])
                ylim([0,1])
                xticks(1:4)
                
                saveas(gcf,strcat(fig_save_path,'R1P2_int',num2str(i_int),'_rep',num2str(rep),'_features_BC' ),'epsc')
                saveas(gcf,strcat(fig_save_path,'R1P2_int',num2str(i_int),'_rep',num2str(rep),'_features_BC' ),'svg')
                
                close
                
            end
            
            for i_int = 1:int_num% size(p38_CV,1)
                figure(i_int)
                Set_figure_size_wide_short
                scatter(1:size(p38_Fano,1),p38_Fano(4:-1:1,i_int,rep),10,[40/255,164/255,255/255],'filled');hold on
                scatter(1:size(NFkB_Fano,1),NFkB_Fano(4:-1:1,i_int,rep),10,[255/255,187/255,40/255],'filled');hold on
                xlim([0.5,4.5])
                ylim([0,3])
                xticks(1:4)
                saveas(gcf,strcat(fig_save_path,'R1P2_int',num2str(i_int),'_rep',num2str(rep),'_features_Fano' ),'epsc')
                saveas(gcf,strcat(fig_save_path,'R1P2_int',num2str(i_int),'_rep',num2str(rep),'_features_Fano' ),'svg')
                
                close
                
            end
        end
        
    end
    
    
    %     figure(1)
    %     h = heatmap( int_vec,denoise_vec,p38_CV,'Colormap',mymap,'CellLabelColor','k');%'none'
    %     % h.Title = strcat('');
    %     caxis([-1,1])
    %     saveas(gcf,strcat(fig_save_path,'R1P2_p38CV_int1to3' ),'epsc')
    %     close
    %
    %     figure(2)
    %     h = heatmap( int_vec,denoise_vec,NFkB_CV,'Colormap',mymap,'CellLabelColor','k');%'none'
    %     % h.Title = strcat('');
    %     caxis([-1,1])
    %     saveas(gcf,strcat(fig_save_path,'R1P2_NFkBCV_int1to3' ),'epsc')
    %     close
end
end