function [] = draw_p38_NFkB_denoise_3model_sim_features_0228(fig_save_path)

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
            
            codon_struc_NFkB = func_metrics(traj.NFkB);
            codon_struc_p38 = func_metrics(traj.p38);
            
            feature_list = {'time2halfMax','peakVec','totalIntegral'};
            
            for i_feature = 1:length(feature_list)
                
                
                if 0
                    figure(i_feature)
                    Set_figure_size_wide_short
                    figure(i_feature)
                    scatter(integral_by_halfhour{i_sig}(:,i_feature),integral_by_halfhour_NFkB{i_sig}(:,i_feature),4,[0,0,1],'filled');hold on
                    xlim([0,0.2]/4)
                    ylim([0,1]/4);
                    figure(i_feature)
                    saveas(gcf,strcat(fig_save_path,'R1P2_',save_metric_name{i_sig},'_int_',num2str(i_feature),'_ranged' ),'epsc')
                    close
                end
                [rho(i_sig,i_feature,rep), pval(i_sig,i_feature,rep)] = corr(codon_struc_NFkB.(feature_list{i_feature})(:,:),...
                    codon_struc_p38.(feature_list{i_feature})(:,:), 'Type', 'Spearman');
                rho(i_sig,i_feature,rep)
                % scatter(integral_by_halfhour{i_sig}(:,i_int),integral_by_halfhour_NFkB{i_sig}(:,i_int))
                
                % Calculate the mean and standard deviation for each vector
                vector1 = codon_struc_p38.(feature_list{i_feature})(:,:);
                mean1 = mean(vector1);
                stdDev1 = std(vector1);
                p38_CV(i_sig,i_feature,rep) = stdDev1 / mean1;
                p38_Fano(i_sig,i_feature,rep) = stdDev1^2 / mean1;
                
                vector2 = codon_struc_NFkB.(feature_list{i_feature})(:,:);
                mean2 = mean(vector2);
                stdDev2 = std(vector2);
                NFkB_CV(i_sig,i_feature,rep) = stdDev2 / mean2;
                NFkB_Fano(i_sig,i_feature,rep) = stdDev2^2 / mean2;
                
                
                [~, p38_BC(i_sig,i_feature,rep)] = bimodalitycoeff(vector1);
                [~,  NFkB_BC(i_sig,i_feature,rep)] = bimodalitycoeff(vector2);
                a =1;
            end
            
        end
        
        % save('NFkB_p38_feature_10rep.mat','rho','pval','p38_CV','NFkB_CV','denoise_label','NFkB_feature_list','feature_idx')
        
    end
    
    if 1
        
        
        if 0 % 4 different denoise
            figure(1)
            
            paperpos = [0,0,80,80]*1.5;
            papersize = [80,80]*1.5;
            draw_pos = [10,10,60,60]*1.5;
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            hmap_plot = mean(rho,3);
            hmap_plot = hmap_plot([4,3,2,1],:);
            h = heatmap(hmap_plot, 'Colormap', mymap, 'CellLabelColor', 'k', 'CellLabelFormat', '%.2f');
            
            % heatmap( mean(rho,3),'Colormap',mymap,'CellLabelColor','k');%'none' % int_vec,denoise_vec,
            XLabels = 1:size(hmap_plot,2);
            %
            % Convert each number in the array into a string
            CustomXLabels = string(XLabels);
            % Replace all but the fifth elements by spaces
            % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
            CustomXLabels(:) = " ";
            
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            h.XDisplayLabels = CustomXLabels;
            
            YLabels = 1:size(hmap_plot,1);
            % Convert each number in the array into a string
            YCustomXLabels = string(YLabels);
            % Replace all but the fifth elements by spaces
            YCustomXLabels(:) = " ";
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            h.YDisplayLabels = YCustomXLabels;
            
            
            % colorbar('off')
            
            % h.Title = strcat('');
            caxis([-1,1])
            saveas(gcf,strcat(fig_save_path,'R1P2_spearman_corr_feature_allrep' ),'epsc')
            saveas(gcf,strcat(fig_save_path,'R1P2_spearman_corr_feature_allrep',num2str(rep) ),'svg')
            
            close
            
        end
        
        
        if 0
            figure(1)
            
            paperpos = [0,0,80,80]*1.5;
            papersize = [80,80]*1.5;
            draw_pos = [10,10,60,60]*1.5;
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            hmap_plot = mean(rho,3);
            hmap_plot = hmap_plot(4:end,:);
            h = heatmap(hmap_plot, 'Colormap', mymap, 'CellLabelColor', 'k', 'CellLabelFormat', '%.2f');
            
            % heatmap( mean(rho,3),'Colormap',mymap,'CellLabelColor','k');%'none' % int_vec,denoise_vec,
            XLabels = 1:size(hmap_plot,2);
            %
            % Convert each number in the array into a string
            CustomXLabels = string(XLabels);
            % Replace all but the fifth elements by spaces
            % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
            CustomXLabels(:) = " ";
            
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            h.XDisplayLabels = CustomXLabels;
            
            YLabels = 1:size(hmap_plot,1);
            % Convert each number in the array into a string
            YCustomXLabels = string(YLabels);
            % Replace all but the fifth elements by spaces
            YCustomXLabels(:) = " ";
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            h.YDisplayLabels = YCustomXLabels;
            
            
            % colorbar('off')
            
            % h.Title = strcat('');
            caxis([-1,1])
            saveas(gcf,strcat(fig_save_path,'R1P2_spearman_corr_feature_allrep_TRIF_denoise',num2str(rep) ),'epsc')
            saveas(gcf,strcat(fig_save_path,'R1P2_spearman_corr_feature_allrep_TRIF_denoise',num2str(rep) ),'svg')
            
            close
            
        end
        
        % subplot(1,length(vis_data_field),i_data_field)
        
        
        
        
        %     figure(2)
        %     h = heatmap( pval,'Colormap',mymap,'CellLabelColor','k');%'none'
        %     % h.Title = strcat('');
        %     caxis([-1,1])
        %     saveas(gcf,strcat(fig_save_path,'R1P2_spearman_rho_int1to3_r',num2str(rep) ),'epsc')
        %     saveas(gcf,strcat(fig_save_path,'R1P2_spearman_rho_int1to3_r',num2str(rep) ),'svg')
        %
        %     close
    end
    
end


%% CV
if 1
    
    if 1
        for rep =1:3
        for i_feature = 1:3% size(p38_CV,1)
            figure(i_feature)
            Set_figure_size_wide_short
            scatter(1:size(p38_CV,1),p38_CV(4:-1:1,i_feature,rep),10,[40/255,164/255,255/255],'filled');hold on
            scatter(1:size(NFkB_CV,1),NFkB_CV(4:-1:1,i_feature,rep),10,[255/255,187/255,40/255],'filled');hold on
            xlim([0.5,4.5])
            ylim([0,2])
            xticks(1:4)
            
            %saveas(gcf,strcat(fig_save_path,'R1P2_',feature_list{i_feature},num2str(rep),'_features_CV' ),'epsc')
            saveas(gcf,strcat(fig_save_path,'R1P2_',feature_list{i_feature},num2str(rep),'_features_CV' ),'svg')
            
            close
            
            
        end
        
        for i_feature = 1:3% size(p38_CV,1)
            figure(i_feature)
            Set_figure_size_wide_short
            scatter(1:size(p38_BC,1),p38_BC(4:-1:1,i_feature,rep),10,[40/255,164/255,255/255],'filled');hold on
            scatter(1:size(NFkB_BC,1),NFkB_BC(4:-1:1,i_feature,rep),10,[255/255,187/255,40/255],'filled');hold on
            xlim([0.5,4.5])
            ylim([0,1])
            xticks(1:4)
            
            %saveas(gcf,strcat(fig_save_path,'R1P2_',feature_list{i_feature},num2str(rep),'_features_BC' ),'epsc')
            saveas(gcf,strcat(fig_save_path,'R1P2_',feature_list{i_feature},num2str(rep),'_features_BC' ),'svg')
            
            close
            
        end
        
        for i_feature = 1:3% size(p38_CV,1)
            figure(i_feature)
            Set_figure_size_wide_short
            scatter(1:size(p38_Fano,1),p38_Fano(4:-1:1,i_feature,rep),10,[40/255,164/255,255/255],'filled');hold on
            scatter(1:size(NFkB_Fano,1),NFkB_Fano(4:-1:1,i_feature,rep),10,[255/255,187/255,40/255],'filled');hold on
            xlim([0.5,4.5])
            ylim([0,3])
            xticks(1:4)
%             saveas(gcf,strcat(fig_save_path,'R1P2_',feature_list{i_feature},num2str(rep),'_features_Fano' ),'epsc')
%             saveas(gcf,strcat(fig_save_path,'R1P2_',feature_list{i_feature},num2str(rep),'_features_Fano' ),'svg')
            
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