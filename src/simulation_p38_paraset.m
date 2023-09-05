function [] = simulation_p38_paraset(p0,paraset_name,fig_save_path)

names = {'NFkBn','IkBaNFkBn','ppp38','ppMKK4','ppMKK6'};%

paper_pos = [0,0,100,80]*1.5;
paper_size = [100,80]*1.5;
font_size = 7;
SimTraj = struct();
load('./raw_data/p38_NFkB_Trajectories_representative.mat')
load('./raw_data/p38_NFkB_Trajectories_representative_new.mat');
sti_vec = {'LPS','TNF','Pam3CSK','CpG'};
%% set the ligand, dose val, dose label , and dose scale.
for i_sti = 1:4
    sti = sti_vec{i_sti};
    switch sti
        case 'TNF' % ligand
            doses = [100];% dose val
            dose_str = {'100ng/mL'};% dose label
            dose_scale = 1/5200.0; % dose scale, DO NOT CHANGE for any ligand. This is the scale factor to change the ligand into uM.
        case 'LPS'
            doses = [100];
            dose_str = {'100ng/mL'};
            dose_scale = 1/24000.0;
            % % % % % % % %
        case 'Pam3CSK'
            doses = 100;
            dose_str = {'100ng/mL'};
            dose_scale = 1/1500.0;
            % % % % % % % % %
        case 'CpG'
            doses = 1000;
            dose_str = {'1000nM'};
            dose_scale = 1/1000.0;
    end
    % sti = 'polyIC';
    % doses = 100000;
    % dose_str = {'100ug/mL'};
    % dose_scale = 1/5000000.0;
    
    
    %% initialization
    options = struct;
    options.DEBUG = 1;
    options.SIM_TIME = 8*60;
    [v0.PARAMS, v0.SPECIES] = nfkbInitialize_p38();%_TNF_TRAF2_paraset2
    
    v0.PARAMS(102,1) =  p0(1);%0.03;0.05 4e-5; % single phosphorilation of MKK6 Vmax
    v0.PARAMS(102,2) = 1; % single phosphorilation of MKK6 Hill coefficient n
    v0.PARAMS(102,3) =  p0(2); % single phosphorilation of MKK6 Kd
    v0.PARAMS(103,1) =  p0(3); %0.14 dephosporilation of single pMKK6
    
    v0.PARAMS(104,1) = v0.PARAMS(102,1); % double phosphorilation of MKK6 Vmax
    v0.PARAMS(104,2) = v0.PARAMS(102,2); % double phosphorilation of MKK6 Hill coefficient n
    v0.PARAMS(104,3) = v0.PARAMS(102,3); % double phosphorilation of MKK6 Kd
    v0.PARAMS(105,1) = v0.PARAMS(103,1); % dephosporilation of double ppMKK6
    
    v0.PARAMS(106,1) =  p0(4); %0.005 % 2.7e-4; % single phosphorilation of MKK4 Vmax
    v0.PARAMS(106,2) = 1; % single phosphorilation of MKK4 Hill coefficient n
    v0.PARAMS(106,3) =  p0(5); %0.0001 single phosphorilation of MKK4 Kd
    v0.PARAMS(107,1) =  p0(6); % 0.8 dephosporilation of single pMKK4
    
    v0.PARAMS(108,1) = v0.PARAMS(106,1); % double phosphorilation of MKK4 Vmax
    v0.PARAMS(108,2) = v0.PARAMS(106,2); % double phosphorilation of MKK4 Hill coefficient n
    v0.PARAMS(108,3) = v0.PARAMS(106,3); % double phosphorilation of MKK4 Kd
    v0.PARAMS(109,1) = v0.PARAMS(107,1); % dephosporilation of double pMKK4
    
    v0.PARAMS(110,1) =  p0(7);%1.6e5; % single phosphorilation of p38 through MKK6 Vmax
    v0.PARAMS(110,2) = 1; % single phosphorilation of p38 through MKK6 Hill coefficient n
    v0.PARAMS(110,3) =  p0(8); % single phosphorilation of p38 through MKK6 Kd
    v0.PARAMS(111,1) =  p0(9);% 345.7; % single phosphorilation of p38 through MKK4 Vmax
    v0.PARAMS(111,2) = 1; % single phosphorilation of p38 through MKK4 Hill coefficient n
    v0.PARAMS(111,3) =  p0(10); % single phosphorilation of p38 through MKK4 Kd
    v0.PARAMS(112,1) =  p0(11); % dephosporilation of single pp38
    
    v0.PARAMS(113,1) = v0.PARAMS(110,1); % single phosphorilation of p38 through MKK6 Vmax
    v0.PARAMS(113,2) = v0.PARAMS(110,2); % single phosphorilation of p38 through MKK6 Hill coefficient n
    v0.PARAMS(113,3) = v0.PARAMS(110,3); % single phosphorilation of p38 through MKK6 Kd
    v0.PARAMS(114,1) = v0.PARAMS(111,1); % single phosphorilation of p38 through MKK4 Vmax
    v0.PARAMS(114,2) = v0.PARAMS(111,2); % single phosphorilation of p38 through MKK4 Hill coefficient n
    v0.PARAMS(114,3) = v0.PARAMS(111,3); % single phosphorilation of p38 through MKK4 Kd
    v0.PARAMS(115,1) = v0.PARAMS(112,1); % dephosporilation of single pp38
    
    options.v.PARAMS = v0.PARAMS;
    options.v.SPECIES = v0.SPECIES;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate all doses (only need to equilibrate on first iteration)
    output = [];
    
    for i = 1:length(doses)
        if isempty(output)
            [t,x,simdata] = nfkbSimulate_p38({sti,doses(i)*dose_scale},names, [], {},options);% _TNF_TRAF2
        else
            options.STEADY_STATE = simdata.STEADY_STATE;
            [t,x] = nfkbSimulate_p38({sti,doses(i)*dose_scale},names, [], {},options); % _TNF_TRAF2
        end
        output = cat(3,output,x);
    end
    
    p38_curves = squeeze(output(:,strcmp(names,'ppp38'),:));
    ppMKK6_curves = squeeze(output(:,strcmp(names,'ppMKK6'),:));
    ppMKK4_curves = squeeze(output(:,strcmp(names,'ppMKK4'),:));
    NFkB_cureves = squeeze(output(:,strcmp(names,'NFkBn'),:)) + squeeze(output(:,strcmp(names,'IkBaNFkBn'),:)) ;
    
    if 1
        figure(1)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        plot(0:480,ppMKK6_curves,'LineWidth',1); hold on
        % title(sti)
        set(gca,'FontSize',font_size,'FontName','Arial')
        ylabel('ppMKK6','FontWeight','b','FontSize',9)
        xlabel('Time','FontWeight','b','FontSize',9)
        xlim([0,480])
        saveas(gcf,strcat(fig_save_path,'Figure6D_',paraset_name,'_sim_ppMKK6',sti),'epsc')
%         saveas(gcf,strcat(fig_save_path,paraset_name,'_sim_ppMKK6',sti),'svg')
        
        close
        
        figure(2)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        plot(0:480,ppMKK4_curves,'LineWidth',1); hold on
        % title(sti)
        set(gca,'FontSize',font_size,'FontName','Arial')
        ylabel('ppMKK4','FontWeight','b','FontSize',9)
        xlabel('Time','FontWeight','b','FontSize',9)
        xlim([0,480])
        saveas(gcf,strcat(fig_save_path,'Figure6D_',paraset_name,'_sim_ppMKK4',sti),'epsc')
%         saveas(gcf,strcat(fig_save_path,paraset_name,'_sim_ppMKK4',sti),'svg')
        
        close
        
    end
    
    if 1
        
        
        figure(i_sti+4)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        plot(0:480,p38_curves,'LineWidth',1); hold on
        
        title(sti)
        set(gca,'FontSize',font_size,'FontName','Arial')
        ylabel('ppp38','FontWeight','b','FontSize',9)
        xlabel('Time','FontWeight','b','FontSize',9)
        xlim([0,480])
        ylim([-0.01,0.04])
        
        saveas(gcf,strcat(fig_save_path,'Figure6C_',paraset_name,'_sim_p38_',sti),'epsc')
%         saveas(gcf,strcat(fig_save_path,paraset_name,'_sim_representative_',sti),'svg')
        
        close
        
        figure(i_sti+4)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        plot(0:480,NFkB_cureves,'LineWidth',1); hold on
        %         plot((0:length(RepTraj.p38.SmoothedTrajectories(i_sti,13:end) )-1)*5, ...
        %             RepTraj.p38.SmoothedTrajectories(i_sti,13:end)/0.73*0.08,'LineWidth',1); hold on
        % legend('sim','exp')
        title(sti)
        set(gca,'FontSize',font_size,'FontName','Arial')
        ylabel('NFkB','FontWeight','b','FontSize',9)
        xlabel('Time','FontWeight','b','FontSize',9)
        xlim([0,480])
        ylim([-0.01,0.3])
        saveas(gcf,strcat(fig_save_path,'Figure6C_',paraset_name,'_sim_NFkB_',sti),'epsc')
        
%         saveas(gcf,strcat(fig_save_path,paraset_name,'_sim_NFkB_',sti),'svg')
        close
        
    end
    SimTraj.p38.(sti_vec{i_sti}) = p38_curves(1:5:end)';
    SimTraj.NFkB.(sti_vec{i_sti}) = NFkB_cureves(1:5:end)';
    SimTraj.MKK4.(sti_vec{i_sti}) = ppMKK4_curves(1:5:end)';
    SimTraj.MKK6.(sti_vec{i_sti}) = ppMKK6_curves(1:5:end)';
    
end
SimTraj.TimePtsInMin = 0:5:480;

save('SimRepresentativeTrajecotories.mat','SimTraj')