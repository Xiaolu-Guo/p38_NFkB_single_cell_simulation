function [] = draw_sens_p38()
% monolix_data_save_file_path = '../SAEM_proj_2023/';
% input_paras.Num_sample = 10;
%
% % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
% input_paras.proj_num_vec = {[3]};
% input_paras.Num_sample = [21];%in total 21 parameters
%
% input_paras.proj_ligand_vec = {{'LPS'}};
% input_paras.proj_dose_str_vec = {{'100ng/mL'}};
% input_paras.proj_dose_val_vec = {{100}};
%
% cal_codon =1;
fig_save_path = './SubFigures2023/';
length_para = 21;
color_mapping = [ linspace(1,0.5,floor(length_para/2))',zeros(floor(length_para/2),2);
    0,0,0;
    zeros(floor(length_para/2),2),linspace(0.5,1,floor(length_para/2))'];
Line_wid = 0.75 * ones(length_para,1);
Line_wid(floor(length_para/2)+1) = 1.5;

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
            
            save_metric_name = strcat( './raw_data/Sim_p38_sens_MKK4S',sigma_string_vec{i_MKK4},...
                '_MKK6S',sigma_string_vec{i_MKK6},...
                '_p38S',sigma_string_vec{i_p38},...
                '_r2.mat');
            fig_name = strcat( 'Sim_p38_sens_MKK4S',sigma_string_vec{i_MKK4},...
                '_MKK6S',sigma_string_vec{i_MKK6},...
                '_p38S',sigma_string_vec{i_p38});
            load(save_metric_name)
            figure(1)
            Set_figure_size_square_small
            
            i_para = 1;
            for i_traj = 2:6:size(data.model_sim{1, 1},1)
                p38_traj = data.model_sim{1, 1}(i_traj,:);
                
                plot(1:length(p38_traj),p38_traj,'LineWidth',Line_wid(i_para),'Color',color_mapping(i_para,:));hold on
                i_para = i_para+1;
            end
            
            ylim([0,0.1*p38_model_fold.p38(end)])
            yt = yticks;
            ytlable = cell(1,length(yt));
            for i_yt = 1:length(yt)
                ytlable{i_yt} = '';
            end
            ytlable{1} = num2str(0);
            ytlable{i_yt} = num2str(0.1*p38_model_fold.p38(end));
            
            xlim([0,96])
            xt = 0:(24):(96);
            xtlable = cell(1,length(xt));
            for i_xt = 1:length(xt)
                xtlable{i_xt} = '';
            end
            
            set(gca,'FontSize',8)
            set(gca,'YTickLabel',ytlable,'XTick',xt,'XTickLabel',xtlable)
            set(gca,'fontsize',7,'fontweight','b')
            xticks(0:12:96);
            % saveas(gcf,strcat(fig_save_path,fig_name),'epsc')
            saveas(gcf,strcat(fig_save_path,fig_name),'svg')
            
            %saveas(gcf,strcat(fig_save_path,fig_save_name,string(parameter_name_vec(i_parameter_name_vec))),'svg')
            %set(gca,'YScale','log')
            
            %saveas(gcf,strcat(fig_save_path,fig_name,'_logScale'),'epsc')
            
            close
            
            
            codon_struc = func_metrics(data.model_sim{1, 1}(2:6:size(data.model_sim{1, 1},1),:),0.01);
            
            
            figure(2)
            
            scatter(MKK4_fold_vec(1,:),codon_struc.peakVec/codon_struc.peakVec(11),10,color_mapping,'^');hold on %,'filled'
            scatter(MKK4_fold_vec(1,:),codon_struc.time2halfMax/codon_struc.time2halfMax(11),10,color_mapping,'d'); hold on % ,'filled'
            % scatter(MKK4_fold_vec(1,:),codon_struc.duration/codon_struc.duration(11),10,color_mapping,'filled','s'); hold on
            
            scatter(MKK4_fold_vec(1,:),codon_struc.totalIntegral/codon_struc.totalIntegral(11),10,color_mapping,'o'); %,'filled'
            % ,'YTick',[0.01,1,100],'YTickLabel',{'','',''}
            
            set(gca,'XScale','log','XTick',[0.1,1,10],'XTickLabel',{'','',''},'YScale','log')
            
            YL = ylim(); %ylim([0,0.1*p38_model_fold.p38(end)])
            yt = yticks;
            
            paperpos = [0,0,85,50];
            papersize = [85,50];
            draw_pos = [10,10,65,30];
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            
            ylim(YL)

            
            ytlable = cell(1,length(yt));
            for i_yt = 1:length(yt)
                ytlable{i_yt} = '';
            end
            ytlable{1} = strcat('10^{',num2str(log10(YL(1)), '%.1f'),'}');
            ytlable{2} = '10^{0.0}';
            ytlable{3} = strcat('10^{',num2str(log10(YL(2)), '%.1f'),'}');
            yticks_new = [YL(1),1,YL(2)];
            
            set(gca,'YTick',yticks_new,'YTickLabel',ytlable)
            xlim([min(MKK4_fold_vec(1,:)),max(MKK4_fold_vec(1,:))])
            
            
            
            xt = MKK4_fold_vec(1,:);
            xtlable = cell(1,length(xt));
            for i_xt = 1:length(xt)
                xtlable{i_xt} = '';
            end
            
            
            saveas(gcf,strcat(fig_save_path,fig_name,'_codon'),'epsc')
            saveas(gcf,strcat(fig_save_path,fig_name),'svg')
            
            close
            
        end
    end
end
end

