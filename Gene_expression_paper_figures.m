
clear all
data_save_file_path = './raw_data/';%_fay_parameter/';
fig_save_path = './SubFigures2023/';

addpath('./src/')
addpath(data_save_file_path)
addpath(fig_save_path)


%% dose curve

savepath = './SubFigures2023/';



%% figure S7A
if 1
    figure(1)
    paperpos=[0,0,150,100];
    papersize=[150 100];
    draw_pos=[10,10,130,90];
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    
    figure(2)
    
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    
    
    for i_para = 1
        
        
        params.ksynmax = 1; % Set this value
        params.n = 3; % Set this value
        params.kdegmax = log(2)/30; % Set this value 30 min e^(-lambda 30) = 1/2 lambda = ln2/30
        params.n_p38 = 3; % Set this value
        params.Kd = 0.1;
        params.Kdp38 = 0.0025;
        
        
        
        if 1
            NFkB_dose = 0.01:0.01:0.3;
            
            figure(1)
            hill_curve(i_para,:) = (params.ksynmax * NFkB_dose.^params.n) ./ (params.Kd^params.n + NFkB_dose.^params.n);
            plot(NFkB_dose,hill_curve(i_para,:),'LineWidth',1);hold on
            % label_kd{i_kd_nfkb} =strcat('Kd_{nfkb} = ', num2str(Kd_nfkb_vec(i_kd_nfkb)));
            % l = legend(label_kd);
            % set(l,'box','off')
            %xlabel('NFkB')
            %ylabel('gene synthesis rate')
            
        end
        
        if 1
            
            p38_dose = 0.001:0.001:0.035;
            
            
            figure(2)
            p38_hill_curve(i_para,:) = (params.kdegmax * params.Kdp38^params.n_p38) ./ (params.Kdp38^params.n_p38 + p38_dose.^params.n_p38);
            plot(p38_dose,p38_hill_curve(i_para,:),'LineWidth',1);hold on
            xlim([0,0.03])
            %xlabel('p38')
            %ylabel('gene degradation rate')
            
        end
    end
    
    
    figure(1)
    %l = legend('Kd = 0.5','Kd = 0.25','Kd = 0.1');%,'para4');
    %set(l,'box','off')
    set(gca,'fontsize',7)
    saveas(gcf,strcat(savepath,'syn_NFkB'),'svg');%epsc
    close
    
    figure(2)
    %l = legend('Kd = 0.05','Kd = 0.025','Kd = 0.01');%,'para4');
    %set(l,'box','off')
    set(gca,'fontsize',7)
    saveas(gcf,strcat(savepath,'deg_p38'),'svg');%epsc
    close
end


%% figure 7B and supp
if 1
    
    kd_NFkB_vec = [0.1];
    kd_p38_vec = [0.0025];
    kd_real_val_vec = [0.1];
    kd_p38_real_val_vec = [0.0025];
    
    ylim_vec = [0,6];
    
    load('./p38_NFkB_Trajectories_all.mat')
    
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
    
    rescale_dominator = [4000,2500,2500,2500];
    
    for i_sti = 1:4
        NFkB_data_plot = AllTraj.nfkb_SmoothedTrajectories.(AllTraj.StimulusOrder{i_sti})* rescale_factor_NFkB;
        NFkB_data_plot = NFkB_data_plot(:,13:end);
        p38_data_plot = AllTraj.p38_SmoothedTrajectories.(AllTraj.StimulusOrder{i_sti}) * rescale_factor_p38;
        p38_data_plot = p38_data_plot(:,13:end);
        ligand_name = AllTraj.StimulusOrder{i_sti};
        
        for i_NFkB_kd = 1:length(kd_NFkB_vec)
            for i_p38_kd = 1:length(kd_p38_vec)
                
                % Define your time vector, assuming it starts at 0 with intervals matching your data
                ts = 0:5:470; % for example, if your measurements are in minutes
                
                % Define your parameters
                params.ksynmax = 1; % Set this value
                params.n = 3; % Set this value
                params.kdegmax = log(2)/30; % Set this value 30 min e^(-lambda 30) = 1/2 lambda = ln2/30
                params.n_p38 = 3; % Set this value
                params.Kdp38 = 0.0025;%0.03;%0.01;%0.05;%0.05; % Set this value
                params.Kd = 0.1;%0.05;% 0.25;%0.25;%0.25; % Set this value

%                 params.Kdp38 = kd_p38_vec(i_p38_kd);%0.03;%0.01;%0.05;%0.05; % Set this value
%                 params.Kd = kd_NFkB_vec(i_NFkB_kd);%0.05;% 0.25;%0.25;%0.25; % Set this value
                                
                file_name = strcat(ligand_name,'_exp_NFkB_kd_',replace(num2str(params.Kd),'.','p'),'p38_kd_',replace(num2str(params.Kdp38),'.','p'));
                % load(strcat('../raw_data/',file_name,'.mat'),'mRNA_NFkB_p38_all','mRNA_NFKB_all','params')
                
                if 1
                    % Define NFKB and p38 as vectors
                    mRNA_NFkB_p38_all = [];
                    mRNA_NFKB_all = [];
                    i_cell_num = 1;
                    if 0
                        for i_cell = 1:size(NFkB_data_plot,1)
                            
                            NFKB = NFkB_data_plot(i_cell,:);%sin(ts); % Your 1x481 vector for NFKB
                            p38 = p38_data_plot(i_cell,:); % Your 1x481 vector for p38
                            NFKB = (NFKB>0).*NFKB;
                            p38 = (p38>0).*p38;
                            p38_0 = p38_data_plot(i_cell,:);
                            p38_0(:) = 0;
                            
                            [T,mRNA_NFkB_p38] = solveODE(NFKB, p38, ts, params);
                            [T,mRNA_NFKB] = solveODE(NFKB, p38_0, ts, params);
                            
                            mRNA_NFkB_p38_all(i_cell_num,:) = mRNA_NFkB_p38';
                            mRNA_NFKB_all(i_cell_num,:) = mRNA_NFKB';
                            i_cell_num = i_cell_num +1;
                        end
                        
                        save(strcat('./raw_data/',file_name,'.mat'),'mRNA_NFkB_p38_all','mRNA_NFKB_all','params')
                    else
                        load(strcat('./raw_data/',file_name,'.mat'))
                        
                    end
                    
                end
                
                if 1
                    figure
                    paperpos=[0,0,100,130]*3;
                    papersize=[100 130]*3;
                    draw_pos=[10,10,90,120]*3;
                    
                    cell_num=size(mRNA_NFkB_p38_all,1);
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    [~,order_cells] = sort(max(mRNA_NFkB_p38_all,[],2),'descend');
                    % subplot(1,length(vis_data_field),i_data_field)
                    h=heatmap(log(mRNA_NFkB_p38_all(order_cells,:)+1),'ColorMap',parula,'GridVisible','off','ColorLimits',[0,8]);
                    %
                    XLabels = 0:5:((size(mRNA_NFkB_p38_all,2)-1)*5);
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
                    saveas(gcf,strcat(savepath,'ge_NFkB_p38_exp_',file_name),'epsc');
                    saveas(gcf,strcat(savepath,'ge_NFkB_p38_exp_',file_name),'svg');
                    
                    close
                    
                    
                    figure
                    paperpos=[0,0,100,130]*3;
                    papersize=[100 130]*3;
                    draw_pos=[10,10,90,120]*3;
                    
                    cell_num=size(mRNA_NFKB_all,1);
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    [~,order_cells] = sort(max(mRNA_NFKB_all,[],2),'descend');
                    % subplot(1,length(vis_data_field),i_data_field)
                    h=heatmap(log(mRNA_NFKB_all(order_cells,:)+1),'ColorMap',parula,'GridVisible','off','ColorLimits',[0,8]);
                    %
                    XLabels = 0:5:((size(mRNA_NFKB_all,2)-1)*5);
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
                    saveas(gcf,strcat(savepath,'ge_NFkB_exp_',file_name),'epsc');
                    saveas(gcf,strcat(savepath,'ge_NFkB_exp_',file_name),'svg');
                    
                    close
                    
                end
                
                if 1
                    clear gene_express_pseudo_counts_NFKB_p38_timepoints gene_express_pseudo_counts_NFKB_timepoints
                    time_index = [4,13,37,size(mRNA_NFkB_p38_all,2)];
                    for i_time = 1:4
                        gene_express_pseudo_counts_NFKB_p38_timepoints(:,i_time) = log(mRNA_NFkB_p38_all(:,time_index(i_time))+1);
                        gene_express_pseudo_counts_NFKB_timepoints(:,i_time) = log(mRNA_NFKB_all(:,time_index(i_time))+1);
                    end
                    
                    figure(1)
                    paperpos=[0,0,130,100]*1.5;
                    papersize=[130 100]*1.5;
                    draw_pos=[10,10,120,90]*1.5;
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    
                    y = gene_express_pseudo_counts_NFKB_p38_timepoints;
                    z = gene_express_pseudo_counts_NFKB_timepoints;
                    % subplot(1,length(vis_data_field),i_data_field)
                    
                    al_goodplot_p38(y,[],0.5,ones(size(gene_express_pseudo_counts_NFKB_p38_timepoints,2),1)*[119 172 48]/255 ,'right',[],std(y(:))/rescale_dominator(i_sti));
                    al_goodplot_p38(z,[],0.5,ones(size(gene_express_pseudo_counts_NFKB_timepoints,2),1)*[237 177 32]/255,'left',[],std(z(:))/rescale_dominator(i_sti));
                    
                    xlim([0.4 4.6])
                    
                    xticks([1 2 3 4])
                    xticklabels({'0.25hr','1hr','3hr','8hr'})
                    %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})
                    
                    ylim(ylim_vec(i_NFkB_kd,:));
                    for i_x = 1:4
                        plot([i_x,i_x],ylim_vec(i_NFkB_kd,:),'--','Color','k');hold on
                    end
                    set(gca,'fontsize',14,'fontname','Arial');
                    saveas(gcf,strcat(savepath,'gene_pseudo_timept_distrib_exp_',ligand_name,'_',file_name),'epsc');
                    saveas(gcf,strcat(savepath,'gene_pseudo_timept_distrib_exp_',ligand_name,'_',file_name),'svg');
                    
                    close
                end
                
                if 0
                    clear gene_express_pseudo_counts_NFKB_p38_timepoints gene_express_pseudo_counts_NFKB_timepoints
                    time_index = [4,13,37,size(mRNA_NFkB_p38_all,2)];
                    for i_time = 1:4
                        gene_express_pseudo_counts_NFKB_p38_timepoints(:,i_time) = mRNA_NFkB_p38_all(:,time_index(i_time))+1;
                        gene_express_pseudo_counts_NFKB_timepoints(:,i_time) = mRNA_NFKB_all(:,time_index(i_time))+1;
                    end
                    
                    figure(1)
                    paperpos=[0,0,130,100]*1.5;
                    papersize=[130 100]*1.5;
                    draw_pos=[10,10,120,90]*1.5;
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    
                    y = gene_express_pseudo_counts_NFKB_p38_timepoints;
                    z = gene_express_pseudo_counts_NFKB_timepoints;
                    % subplot(1,length(vis_data_field),i_data_field)
                    
                    al_goodplot_p38(y,[],0.5,ones(size(gene_express_pseudo_counts_NFKB_p38_timepoints,2),1)*[119 172 48]/255 ,'right',[],std(y(:))/rescale_dominator(i_sti));
                    al_goodplot_p38(z,[],0.5,ones(size(gene_express_pseudo_counts_NFKB_timepoints,2),1)*[237 177 32]/255,'left',[],std(z(:))/rescale_dominator(i_sti));
                    
                    xlim([0.4 4.6])
                    
                    xticks([1 2 3 4])
                    xticklabels({'0.25hr','1hr','3hr','8hr'})
                    %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})
                    
                    ylim(ylim_vec(i_NFkB_kd,:));
                    for i_x = 1:4
                        plot([i_x,i_x],ylim_vec(i_NFkB_kd,:),'--','Color','k');hold on
                    end
                    set(gca,'fontsize',14,'fontname','Arial');
                    saveas(gcf,strcat(savepath,'gene_exp_timept_distrib_exp_',ligand_name,'_',file_name),'epsc');
                    saveas(gcf,strcat(savepath,'gene_exp_timept_distrib_exp_',ligand_name,'_',file_name),'svg');
                    
                    close
                end
                
            end
        end
    end
    
end


%% gene expression model function


function [T,mRNA] = solveODE(NFKB, p38, ts, params)
% Interpolate NFKB and p38 values to evaluate at ODE solver time points
NFKB_interp = @(t) interp1(ts, NFKB, t, 'linear', 'extrap');
p38_interp = @(t) interp1(ts, p38, t, 'linear', 'extrap');

% Define the ODE function
    function dmdt = odeFun(t, mRNA)
        NFkB_t = NFKB_interp(t);
        p38_t = p38_interp(t);
        dmdt = (params.ksynmax * NFkB_t^params.n) / (params.Kd^params.n + NFkB_t^params.n) ...
            - (params.kdegmax * mRNA * params.Kdp38^params.n_p38) / (params.Kdp38^params.n_p38 + p38_t^params.n_p38);
    end

% Set initial condition for mRNA
% mRNA0 = 0; % assuming mRNA starts at 0, but you should set this accordingly

NFkB_0 = NFKB_interp(0);
p38_0 = p38_interp(0);

mRNA0 = ((params.ksynmax * NFkB_0^params.n) / (params.Kd^params.n + NFkB_0^params.n))/...
    ((params.kdegmax * params.Kdp38^params.n_p38) / (params.Kdp38^params.n_p38 + p38_0^params.n_p38));

% odeFun(0, mRNA0)
% Solve the ODE using ode45
[T, mRNA] = ode15s(@odeFun, ts, mRNA0);

% Return only the mRNA values if you don't need the time points
end