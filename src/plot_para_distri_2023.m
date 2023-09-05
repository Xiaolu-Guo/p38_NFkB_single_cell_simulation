function [] = plot_para_distri_2023(params,save_path)

mu = log( 0.1);
sigma = (log(2)/3);

% Generate values for x
x = linspace(0, 0.25, 1000);  % Adjust the range of x values if needed

% Compute the log-normal PDF values
pdf_values = lognpdf(x, mu, sigma);


for i_ligand = 1:length(params.info_ligand)
%     for i_params = 1:size(params.params{i_ligand},2)
%         para_name = params.params_name{i_ligand}(1,i_params);
%         figure(1)
%         h1 = histogram(params.params{i_ligand}(:,i_params),20);hold on
%         
%         h1.FaceColor = [0.5,0,0];
%         % h1.Normalization = 'countdensity';
%         h1.Normalization = 'pdf';
%         
%         plot(x, pdf_values, 'b','LineWidth', 1.5);
%         xlim([0,0.25])
%         
%         xlabel(para_name);
%         set(gca,'fontsize',7)
%         Set_figure_size_wide_short
%         
%         saveas(gcf, strcat(save_path,'Figure6E_params_distrib_',params.info_ligand{i_ligand},'_',string(para_name)),'epsc')
%         % saveas(gcf, strcat(save_path,'params_distrib_',params.info_ligand{i_ligand},'_',string(para_name)),'svg')
%         close
%     end
    
    for i_init_species = 1:size(params.init_species{i_ligand},2)
        para_name = params.init_species_name{i_ligand}(1,i_init_species);
        figure(1)
        h1 = histogram(params.init_species{i_ligand}(:,i_init_species),20);hold on        
        h1.FaceColor = [0.5,0,0];
        h1.Normalization = 'countdensity';
        % h1.Normalization = 'pdf';
        
        xlim([0,0.25])        
        xlabel(replace(string(para_name),'_','-'));
        set(gca,'fontsize',7)
        Set_figure_size_wide_short
        
        saveas(gcf, strcat(save_path,'Figure6E_params_sampling_distrib_',params.info_ligand{i_ligand},'_',string(para_name)),'epsc')
        saveas(gcf, strcat(save_path,'Figure6E_params_sampling_distrib_',params.info_ligand{i_ligand},'_',string(para_name)),'svg')
        close
        
        figure(1)
        
        plot(x, pdf_values, 'b','LineWidth', 1.5);
        xlim([0,0.25])
        
        xlabel(replace(string(para_name),'_','-'));
        set(gca,'fontsize',7)
        Set_figure_size_wide_short
        
        saveas(gcf, strcat(save_path,'Figure6E_params_theoretical_distrib_',params.info_ligand{i_ligand},'_',string(para_name)),'epsc')
        saveas(gcf, strcat(save_path,'Figure6E_params_theoretical_distrib_',params.info_ligand{i_ligand},'_',string(para_name)),'svg')
        close
    end
end



