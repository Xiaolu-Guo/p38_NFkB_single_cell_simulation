MKK4_mean_fold_vec = [0.5,1,2];
MKK6_mean_fold_vec = [0.5,1,2];
p38_mean_fold_vec = [0.5,1,2];
sigma_string_vec = {'p5fold','1fold','2fold'};
savepath = './SubFigures2023/';


paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;

for i_MKK4 = 1:length(MKK4_mean_fold_vec)
    for i_MKK6 = 1:length(MKK6_mean_fold_vec)
        for i_p38 = 1:length(p38_mean_fold_vec)
            
            
            
            save_metric_name = strcat( 'Sim_NFKB_sample_p38_sens_MKK4S',sigma_string_vec{i_MKK4},...
                '_MKK6S',sigma_string_vec{i_MKK6},...
                '_p38S',sigma_string_vec{i_p38},...
                '_r1.mat');
            load(strcat(data_save_file_path,save_metric_name))
            data_plot = data.model_sim{2}(2:6:end,:);
            [~,order_ind] = sort(sum(data_plot,2),'descend');
            
            figure(1)
            
            cell_num=length(order_ind);
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            
            % subplot(1,length(vis_data_field),i_data_field)
            h=heatmap(data_plot(order_ind,:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.1]);%[-0.001,0.2] for TNF
            %
            XLabels = 0:5:((size(data_plot,2)-1)*5);
            % Convert each number in the array into a string
            CustomXLabels = string(XLabels/60);
            % Replace all but the fifth elements by spaces
            % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
            CustomXLabels(:) = " ";
            
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            h.XDisplayLabels = CustomXLabels;
            
            YLabels = 1:cell_num;
            % Convert each number in the array into a string
            YCustomXLabels = string(YLabels);
            % Replace all but the fifth elements by spaces
            YCustomXLabels(:) = " ";
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            h.YDisplayLabels = YCustomXLabels;
            
            % xlabel('Time (hours)');
            % ylabel(vis_data_field{i_data_field});
            % clb=colorbar;
            % clb.Label.String = 'NFkB(A.U.)';
            colorbar('off')
            
            set(gca,'fontsize',14,'fontname','Arial');
            saveas(gcf,strcat(savepath,...
                'Heatmap_Sim_NFKB_sample_p38_sens_MKK4S',sigma_string_vec{i_MKK4},...
                '_MKK6S',sigma_string_vec{i_MKK6},...
                '_p38S',sigma_string_vec{i_p38},...
                '_r1'),'epsc');
            close
        end
    end
end