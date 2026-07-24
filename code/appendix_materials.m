%APPENDIX_MATERIALS
%The script only reads existing baseline and IRF output files. It does not
%modify ss_one_industry.m, ss_two_firms_both_paths.m, or test.m.

clear; clc; close all;

root = fileparts(mfilename('fullpath'));
results_dir = fullfile(root,'Results');
appendix_dir = fullfile(results_dir,'Appendix');

if exist(appendix_dir,'dir') ~= 7
    mkdir(appendix_dir);
end

baseline_paths = { ...
    fullfile(root,'German_Baseline_OneType_TAUCHEN.mat'), ...
    fullfile(root,'German_Baseline_OneType_ROUWENHORST.mat'), ...
    fullfile(root,'German_Baseline_TwoTypes_TAUCHEN.mat'), ...
    fullfile(root,'German_Baseline_TwoTypes_ROUWENHORST.mat')};

for i = 1:numel(baseline_paths)
    assert(exist(baseline_paths{i},'file') == 2, ...
        'Missing baseline file: %s',baseline_paths{i});
end

one_tau = load(baseline_paths{1});
one_rou = load(baseline_paths{2});
two_tau = load(baseline_paths{3});
two_rou = load(baseline_paths{4});

rho_path = latest_result_file(results_dir, ...
    'IRF_two_type_tauchen_marginexit_rho_robustness_*.mat');
shock_path = latest_result_file(results_dir, ...
    'IRF_two_type_tauchen_marginexit_shock_size_robustness_*.mat');
main_irf_path = latest_result_file(results_dir, ...
    'IRF_two_type_tauchen_marginexit_yearly_*.mat');

rho_data = load(rho_path);
shock_data = load(shock_path);
main_irf = load(main_irf_path);

assert(isfield(rho_data,'rho_robust_irfs'), ...
    'The persistence-robustness file does not contain rho_robust_irfs.');
assert(isfield(shock_data,'shock_size_robust_irfs'), ...
    'The shock-size file does not contain shock_size_robust_irfs.');

%% B. Productivity discretization

specification = ["One industry"; "One industry"; "Two type"; "Two type"];
method = ["Tauchen"; "Rouwenhorst"; "Tauchen"; "Rouwenhorst"];
n_states = [one_tau.n_theta; one_rou.n_theta; two_tau.n_theta; two_rou.n_theta];

rho_theta = [0.9500; 0.9500; two_tau.rho_theta; two_rou.rho_theta];
sigma_theta = [0.1200; 0.1200; two_tau.sigma_theta; two_rou.sigma_theta];

theta_min = [min(one_tau.theta_grid); min(one_rou.theta_grid); ...
    min(two_tau.theta_grid); min(two_rou.theta_grid)];
theta_max = [max(one_tau.theta_grid); max(one_rou.theta_grid); ...
    max(two_tau.theta_grid); max(two_rou.theta_grid)];
stationarity_error = [stationarity_residual(one_tau); ...
    stationarity_residual(one_rou); stationarity_residual(two_tau); ...
    stationarity_residual(two_rou)];

discretization_table = table(specification,method,n_states,rho_theta, ...
    sigma_theta,theta_min,theta_max,stationarity_error, ...
    'VariableNames',{'specification','method','n_theta','rho_theta', ...
    'sigma_theta','theta_min','theta_max','stationarity_error'});

discretization_csv = fullfile(appendix_dir, ...
    'Table_A1_discretization_diagnostics.csv');
discretization_tex = fullfile(appendix_dir, ...
    'Table_A1_discretization_diagnostics.tex');
writetable(discretization_table,discretization_csv);
write_discretization_latex(discretization_tex,discretization_table);

%% Steady state comparison table

[~,tau_high] = max(two_tau.gamma_type);
tau_low = 3 - tau_high;
[~,rou_high] = max(two_rou.gamma_type);
rou_low = 3 - rou_high;

output = [one_tau.y_agg; one_rou.y_agg; two_tau.y_agg; two_rou.y_agg];
energy = [one_tau.e_agg; one_rou.e_agg; two_tau.e_agg; two_rou.e_agg];
firm_mass = [one_tau.M; one_rou.M; two_tau.M_total; two_rou.M_total];
exit_aggregate_pct = 100*[one_tau.exit_rate; one_rou.exit_rate; ...
    two_tau.exit_rate_agg; two_rou.exit_rate_agg];
exit_low_pct = [NaN; NaN; 100*two_tau.exit_rate_type(tau_low); ...
    100*two_rou.exit_rate_type(rou_low)];
exit_high_pct = [NaN; NaN; 100*two_tau.exit_rate_type(tau_high); ...
    100*two_rou.exit_rate_type(rou_high)];

steady_state_table = table(specification,method,output,energy,firm_mass, ...
    exit_aggregate_pct,exit_low_pct,exit_high_pct, ...
    'VariableNames',{'specification','method','output','energy_use', ...
    'firm_mass','exit_aggregate_pct','exit_low_pct','exit_high_pct'});

steady_state_csv = fullfile(appendix_dir, ...
    'Table_A2_steady_state_method_comparison.csv');
steady_state_tex = fullfile(appendix_dir, ...
    'Table_A2_steady_state_method_comparison.tex');
writetable(steady_state_table,steady_state_csv);
write_steady_state_latex(steady_state_tex,steady_state_table);

%% Productivity discretization distributions

fig_grid = figure('Color','w','Visible','off', ...
    'Name','Productivity Discretization Diagnostics', ...
    'Position',[120 100 980 560]);
tl_grid = tiledlayout(fig_grid,1,2,'TileSpacing','loose','Padding','loose');
tl_grid.Units = 'normalized';
tl_grid.Position = [0.075 0.295 0.90 0.615];

ax = nexttile(tl_grid);
plot(ax,log(one_tau.theta_grid(:)),100*one_tau.Q_entrant(:), ...
    'k-','LineWidth',1.8);
hold(ax,'on');
plot(ax,log(one_rou.theta_grid(:)),100*one_rou.Q_entrant(:), ...
    '-','Color',[0.00 0.20 0.70],'LineWidth',1.8);
title(ax,'One-industry benchmark','FontWeight','normal','FontSize',13);
ylabel(ax,'Probability mass (\%)','Interpreter','latex');
xlabel(ax,'Log productivity, $z$','Interpreter','latex');
format_appendix_axis(ax);

ax = nexttile(tl_grid);
plot(ax,log(two_tau.theta_grid(:)),100*two_tau.Q_entrant(:), ...
    'k-','LineWidth',1.8);
hold(ax,'on');
plot(ax,log(two_rou.theta_grid(:)),100*two_rou.Q_entrant(:), ...
    '-','Color',[0.00 0.20 0.70],'LineWidth',1.8);
title(ax,'Two-type model','FontWeight','normal','FontSize',13);
ylabel(ax,'Probability mass (\%)','Interpreter','latex');
xlabel(ax,'Log productivity, $z$','Interpreter','latex');
format_appendix_axis(ax);

add_method_legend(fig_grid);
export_appendix_figure(fig_grid,appendix_dir, ...
    'Figure_A1_productivity_discretization');
close(fig_grid);

%% C. Aditional robustness

rho_low_irf = select_irf(rho_data.rho_robust_irfs,'rho_pe',0.50);
rho_high_irf = select_irf(rho_data.rho_robust_irfs,'rho_pe',0.95);
shock_low_irf = select_irf(shock_data.shock_size_robust_irfs,'shock_size',0.10);
shock_high_irf = select_irf(shock_data.shock_size_robust_irfs,'shock_size',0.40);

make_full_irf_figure(rho_low_irf,appendix_dir, ...
    'Figure_A2_full_irf_rho_050');
make_full_irf_figure(rho_high_irf,appendix_dir, ...
    'Figure_A3_full_irf_rho_095');
make_full_irf_figure(shock_low_irf,appendix_dir, ...
    'Figure_A4_full_irf_shock_010');
make_full_irf_figure(shock_high_irf,appendix_dir, ...
    'Figure_A5_full_irf_shock_040');

%% LOCAL FUNCTIONS

function path = latest_result_file(folder,pattern)
    files = dir(fullfile(folder,pattern));
    if isempty(files)
        error(['No file matching %s was found in %s. Run test.m once ' ...
            'before generating the appendix.'],pattern,folder);
    end
    [~,idx] = max([files.datenum]);
    path = fullfile(files(idx).folder,files(idx).name);
end

function residual = stationarity_residual(S)
    require_fields(S,{'Q_incumbent','Q_entrant'},'steady-state file');
    q = S.Q_entrant(:);
    q = q/sum(q);
    residual = max(abs(S.Q_incumbent'*q-q));
end

function require_fields(S,names,source_name)
    for i = 1:numel(names)
        if ~isfield(S,names{i})
            error('%s is missing required field %s.',source_name,names{i});
        end
    end
end

function irf = select_irf(irfs,field_name,target)
    values = arrayfun(@(x) x.(field_name),irfs);
    [distance,idx] = min(abs(values-target));
    assert(distance < 1e-10, ...
        'Requested %s = %.4f is unavailable.',field_name,target);
    irf = irfs(idx);
end

function make_full_irf_figure(irf,outdir,base_name)
    required = {'years','rho_pe','shock_size','Y_index','Y_low_index', ...
        'Y_high_index','E_index','E_low_index','E_high_index', ...
        'exit_rate_agg_pct','exit_rate_low_pct','exit_rate_high_pct', ...
        'mass_active_index','mass_active_low_index','mass_active_high_index'};
    require_fields(irf,required,'robustness IRF');

    n_plot = min(50,numel(irf.years)-1);
    plot_idx = 2:(n_plot+1);
    t_raw = 0:(n_plot-1);
    t_fine = linspace(0,n_plot-1,8*n_plot);

    log_shock = log(1+irf.shock_size);
    pe_dev = 100*(exp(log_shock*(irf.rho_pe.^(0:n_plot-1)))-1);

    exit_agg = irf.exit_rate_agg_pct(plot_idx)-irf.exit_rate_agg_pct(1);
    exit_low = irf.exit_rate_low_pct(plot_idx)-irf.exit_rate_low_pct(1);
    exit_high = irf.exit_rate_high_pct(plot_idx)-irf.exit_rate_high_pct(1);
    mass_agg = irf.mass_active_index(plot_idx)-100;
    mass_low = irf.mass_active_low_index(plot_idx)-100;
    mass_high = irf.mass_active_high_index(plot_idx)-100;
    y_agg = irf.Y_index(plot_idx)-100;
    y_low = irf.Y_low_index(plot_idx)-100;
    y_high = irf.Y_high_index(plot_idx)-100;
    e_agg = irf.E_index(plot_idx)-100;
    e_low = irf.E_low_index(plot_idx)-100;
    e_high = irf.E_high_index(plot_idx)-100;

    fig = figure('Color','w','Visible','off', ...
        'Name',base_name,'Position',[160 60 980 780]);
    tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
    tl.Units = 'normalized';
    tl.Position = [0.065 0.205 0.91 0.745];

    black = [0 0 0];
    blue = [0.00 0.20 0.70];
    red = [0.70 0.05 0.10];

    plot_aggregate_irf(nexttile(tl),t_raw,t_fine,pe_dev,black, ...
        '$P_e$','\% dev. from ss',false,n_plot);
    plot_aggregate_irf(nexttile(tl),t_raw,t_fine,exit_agg,black, ...
        '$\delta$','p.p. from ss',false,n_plot);
    plot_type_irf(nexttile(tl),t_raw,t_fine,exit_low,exit_high,blue,red, ...
        '$\delta_s$','p.p. from ss',false,n_plot);
    plot_aggregate_irf(nexttile(tl),t_raw,t_fine,mass_agg,black, ...
        '$M$','\% dev. from ss',false,n_plot);
    plot_type_irf(nexttile(tl),t_raw,t_fine,mass_low,mass_high,blue,red, ...
        '$M_s$','\% dev. from ss',false,n_plot);
    plot_aggregate_irf(nexttile(tl),t_raw,t_fine,y_agg,black, ...
        '$Y$','\% dev. from ss',false,n_plot);
    plot_type_irf(nexttile(tl),t_raw,t_fine,y_low,y_high,blue,red, ...
        '$Y_s$','\% dev. from ss',true,n_plot);
    plot_aggregate_irf(nexttile(tl),t_raw,t_fine,e_agg,black, ...
        '$E$','\% dev. from ss',true,n_plot);
    plot_type_irf(nexttile(tl),t_raw,t_fine,e_low,e_high,blue,red, ...
        '$E_s$','\% dev. from ss',true,n_plot);

    add_firm_type_legend(fig);
    export_appendix_figure(fig,outdir,base_name);
    close(fig);
end

function plot_aggregate_irf(ax,t_raw,t_fine,series,color,panel_title, ...
        y_label,show_x,n_plot)
    curve = interp1(t_raw,series(:),t_fine,'pchip');
    plot(ax,t_fine,curve,'Color',color,'LineWidth',1.8);
    hold(ax,'on');
    yline(ax,0,'--','Color',[0.45 0.45 0.45],'LineWidth',1.2);
    title(ax,panel_title,'Interpreter','latex','FontWeight','normal');
    ylabel(ax,y_label,'Interpreter','latex');
    if show_x
        xlabel(ax,'Years');
    end
    format_irf_axis(ax,n_plot);
end

function plot_type_irf(ax,t_raw,t_fine,low_series,high_series, ...
        blue,red,panel_title,y_label,show_x,n_plot)
    low_curve = interp1(t_raw,low_series(:),t_fine,'pchip');
    high_curve = interp1(t_raw,high_series(:),t_fine,'pchip');
    plot(ax,t_fine,low_curve,'Color',blue,'LineWidth',1.8);
    hold(ax,'on');
    plot(ax,t_fine,high_curve,'Color',red,'LineWidth',1.8);
    yline(ax,0,'--','Color',[0.45 0.45 0.45],'LineWidth',1.2);
    title(ax,panel_title,'Interpreter','latex','FontWeight','normal');
    ylabel(ax,y_label,'Interpreter','latex');
    if show_x
        xlabel(ax,'Years');
    end
    format_irf_axis(ax,n_plot);
end

function format_irf_axis(ax,n_plot)
    xlim(ax,[0 n_plot-1]);
    grid(ax,'off');
    box(ax,'off');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 13;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
end

function format_appendix_axis(ax)
    grid(ax,'off');
    box(ax,'off');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 13;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
end

function add_method_legend(fig)
    box_pos = [0.255 0.045 0.490 0.085];
    annotation(fig,'rectangle',box_pos,'Color',[0 0 0],'LineWidth',1.0);
    annotation(fig,'line',[0.300 0.355],[0.087 0.087], ...
        'Color',[0 0 0],'LineWidth',1.8);
    annotation(fig,'textbox',[0.365 0.067 0.135 0.040], ...
        'String','Tauchen','EdgeColor','none','FontName','Times New Roman', ...
        'FontSize',11,'VerticalAlignment','middle');
    annotation(fig,'line',[0.515 0.570],[0.087 0.087], ...
        'Color',[0.00 0.20 0.70],'LineWidth',1.8);
    annotation(fig,'textbox',[0.580 0.067 0.145 0.040], ...
        'String','Rouwenhorst','EdgeColor','none','FontName','Times New Roman', ...
        'FontSize',11,'VerticalAlignment','middle');
end

function add_firm_type_legend(fig)
    box_pos = [0.130 0.062 0.740 0.060];
    annotation(fig,'rectangle',box_pos,'Color',[0 0 0],'LineWidth',1.0);
    annotation(fig,'line',[0.160 0.205],[0.092 0.092], ...
        'Color',[0 0 0],'LineWidth',1.8);
    annotation(fig,'textbox',[0.215 0.074 0.130 0.034], ...
        'String','Aggregate','EdgeColor','none','FontName','Times New Roman', ...
        'FontSize',11,'VerticalAlignment','middle');
    annotation(fig,'line',[0.385 0.430],[0.092 0.092], ...
        'Color',[0.00 0.20 0.70],'LineWidth',1.8);
    annotation(fig,'textbox',[0.440 0.074 0.175 0.034], ...
        'String','Low-energy firms','EdgeColor','none', ...
        'FontName','Times New Roman','FontSize',11, ...
        'VerticalAlignment','middle');
    annotation(fig,'line',[0.635 0.680],[0.092 0.092], ...
        'Color',[0.70 0.05 0.10],'LineWidth',1.8);
    annotation(fig,'textbox',[0.690 0.074 0.180 0.034], ...
        'String','High-energy firms','EdgeColor','none', ...
        'FontName','Times New Roman','FontSize',11, ...
        'VerticalAlignment','middle');
end

function plot_cutoff_value_panel(ax,theta,baseline,impact,color,panel_title)
    plot(ax,theta,baseline,'k--','LineWidth',1.5);
    hold(ax,'on');
    plot(ax,theta,impact,'Color',color,'LineWidth',1.8);
    yline(ax,0,'--','Color',[0.45 0.45 0.45],'LineWidth',1.0);
    title(ax,panel_title,'Interpreter','latex','FontWeight','normal');
    ylabel(ax,'Continuation value');
    format_appendix_axis(ax);
end

function plot_exit_probability_panel(ax,theta,baseline,impact,color,panel_title)
    plot(ax,theta,baseline,'k--','LineWidth',1.5);
    hold(ax,'on');
    plot(ax,theta,impact,'Color',color,'LineWidth',1.8);
    title(ax,panel_title,'Interpreter','latex','FontWeight','normal');
    ylabel(ax,'Exit probability');
    xlabel(ax,'Log productivity, $z$','Interpreter','latex');
    ylim(ax,[-0.02 1.02]);
    format_appendix_axis(ax);
end

function add_cutoff_legend(fig)
    box_pos = [0.130 0.055 0.740 0.080];
    annotation(fig,'rectangle',box_pos,'Color',[0 0 0],'LineWidth',1.0);
    annotation(fig,'line',[0.165 0.210],[0.095 0.095], ...
        'Color',[0 0 0],'LineStyle','--','LineWidth',1.5);
    annotation(fig,'textbox',[0.220 0.077 0.145 0.034], ...
        'String','Steady state','EdgeColor','none','FontName','Times New Roman', ...
        'FontSize',11,'VerticalAlignment','middle');
    annotation(fig,'line',[0.390 0.435],[0.095 0.095], ...
        'Color',[0.00 0.20 0.70],'LineWidth',1.8);
    annotation(fig,'textbox',[0.445 0.077 0.175 0.034], ...
        'String','Low-energy impact','EdgeColor','none', ...
        'FontName','Times New Roman','FontSize',11, ...
        'VerticalAlignment','middle');
    annotation(fig,'line',[0.645 0.690],[0.095 0.095], ...
        'Color',[0.70 0.05 0.10],'LineWidth',1.8);
    annotation(fig,'textbox',[0.700 0.077 0.175 0.034], ...
        'String','High-energy impact','EdgeColor','none', ...
        'FontName','Times New Roman','FontSize',11, ...
        'VerticalAlignment','middle');
end

function export_appendix_figure(fig,outdir,base_name)
    png_path = fullfile(outdir,[base_name '.png']);
    pdf_path = fullfile(outdir,[base_name '.pdf']);
    exportgraphics(fig,png_path,'Resolution',300);
    try
        exportgraphics(fig,pdf_path,'ContentType','vector');
    catch
        exportgraphics(fig,pdf_path,'ContentType','image','Resolution',300);
    end
end

function write_discretization_latex(path,T)
    fid = fopen(path,'w');
    assert(fid >= 0,'Could not write %s.',path);
    cleaner = onCleanup(@() fclose(fid)); 
    fprintf(fid,'\\begin{table}[htbp]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\caption{Productivity discretization diagnostics}\n');
    fprintf(fid,'\\label{tab:appendix-discretization}\n');
    fprintf(fid,'\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llrrrrrr@{}}\n');
    fprintf(fid,'\\toprule\n');
    fprintf(fid,'Specification & Method & $N_\\theta$ & $\\rho_\\theta$ & $\\sigma_\\theta$ & $\\theta_{\\min}$ & $\\theta_{\\max}$ & Stationarity error \\\\\n');
    fprintf(fid,'\\midrule\n');
    for i = 1:height(T)
        fprintf(fid,'%s & %s & %d & %.4f & %.4f & %.4f & %.4f & %.2e \\\\\n', ...
            T.specification(i),T.method(i),T.n_theta(i),T.rho_theta(i), ...
            T.sigma_theta(i),T.theta_min(i),T.theta_max(i), ...
            T.stationarity_error(i));
    end
    fprintf(fid,'\\bottomrule\n');
    fprintf(fid,'\\end{tabular*}\n');
    fprintf(fid,'\\begin{minipage}{\\textwidth}\n');
    fprintf(fid,'\\footnotesize\\textit{Notes:} The entrant distribution is the stationary distribution of the incumbent transition matrix. The stationarity error is $\\max_i |(Q''q)_i-q_i|$. The one-industry persistence and innovation volatility are taken directly from \\texttt{ss\\_one\\_industry.m}; the remaining entries are stored in the baseline files.\n');
    fprintf(fid,'\\end{minipage}\n');
    fprintf(fid,'\\end{table}\n');
end

function write_steady_state_latex(path,T)
    fid = fopen(path,'w');
    assert(fid >= 0,'Could not write %s.',path);
    cleaner = onCleanup(@() fclose(fid)); 
    fprintf(fid,'\\begin{table}[htbp]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\caption{Steady-state moments across discretization methods}\n');
    fprintf(fid,'\\label{tab:appendix-steady-state-methods}\n');
    fprintf(fid,'\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llrrrrrr@{}}\n');
    fprintf(fid,'\\toprule\n');
    fprintf(fid,'Specification & Method & $Y$ & $E$ & $M$ & $\\delta$ & $\\delta_L$ & $\\delta_H$ \\\\\n');
    fprintf(fid,'\\midrule\n');
    for i = 1:height(T)
        fprintf(fid,'%s & %s & %.4f & %.4f & %.4f & %.3f & %s & %s \\\\\n', ...
            T.specification(i),T.method(i),T.output(i),T.energy_use(i), ...
            T.firm_mass(i),T.exit_aggregate_pct(i), ...
            latex_number(T.exit_low_pct(i),'%.3f'), ...
            latex_number(T.exit_high_pct(i),'%.3f'));
    end
    fprintf(fid,'\\bottomrule\n');
    fprintf(fid,'\\end{tabular*}\n');
    fprintf(fid,'\\begin{minipage}{\\textwidth}\n');
    fprintf(fid,'\\footnotesize\\textit{Notes:} Exit rates are percentages per model period. Low- and high-energy exit rates apply only to the two-type specification. The table is a numerical-method check; cross-specification differences are not treated as robustness estimates.\n');
    fprintf(fid,'\\end{minipage}\n');
    fprintf(fid,'\\end{table}\n');
end

function value = latex_number(x,format_spec)
    if isnan(x)
        value = '--';
    else
        value = sprintf(format_spec,x);
    end
end
