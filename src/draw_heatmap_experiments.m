load('./raw_data/p38_NFkB_Trajectories_all.mat')

rescale_data_or_not = 1;

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

%% set fig options
paperpos=[0,0,100,130]*3;
papersize=[100 130]*3;
draw_pos=[10,10,90,120]*3;
savepath = './SubFigures2023/';



%% heatmap plot ordered by p38 peak
if 1
for i_sti = 1
    NFkB_data_plot = AllTraj.nfkb_SmoothedTrajectories.(AllTraj.StimulusOrder{i_sti})* rescale_factor_NFkB;
    NFkB_data_plot = NFkB_data_plot(:,14:end);
    p38_data_plot = AllTraj.p38_SmoothedTrajectories.(AllTraj.StimulusOrder{i_sti}) * rescale_factor_p38;
    p38_data_plot = p38_data_plot(:,14:end);
    ligand_name = AllTraj.StimulusOrder{i_sti};
    [~,data_order{i_sti}] = sort(max(p38_data_plot,[],2),'descend');
    figure(1)
    
    cell_num=length(data_order{i_sti});
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    
    % subplot(1,length(vis_data_field),i_data_field)
    h=heatmap(p38_data_plot(data_order{i_sti},:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.05]);%[-0.001,0.2] for TNF
    %
    XLabels = 0:5:((size(p38_data_plot,2)-1)*5);
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
    saveas(gcf,strcat(savepath,'Figure6F_',ligand_name,'_p38_','ordered_by_p38PeakAmp'),'epsc');% svg
    close
     
    
    figure(1)
    
    cell_num=length(data_order{i_sti});
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    
    % subplot(1,length(vis_data_field),i_data_field)
    h=heatmap(NFkB_data_plot(data_order{i_sti},:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
    %
    XLabels = 0:5:((size(NFkB_data_plot,2)-1)*5);
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
    saveas(gcf,strcat(savepath,'Figure6F_',ligand_name,'_NFkB_','ordered_by_p38PeakAmp'),'epsc');%epsc svg
    close
    
    
    
end
end