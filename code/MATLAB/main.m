clc
clear
close all

%% define user input variables

modes = ["fusion", "switch"]; active_mode = modes(1);
flight = "20210910_F01";
flight_type = "orbits";
set = 2; orbit = 1;

%% import data

disp('>>> Importing data ... ');

if ispc
    delimiter = "\";
elseif ismac || isunix
    delimiter = "/";
end

path = string(erase(pwd, delimiter + "code" + delimiter + "MATLAB")); path = path + delimiter + "data" + delimiter; 
warning('off');
data = readtable(path + flight + '.csv');
warning('on');

% import timestamps for analysis

timestamps = readtable(string(pwd) + delimiter + "input" + delimiter + "timestamps" + delimiter + flight + ".csv");
timestamps.START = datetime(eraseBetween(flight + " " + string(timestamps.START), 9, 12),'InputFormat','yyyyMMdd HH:mm:ss');
timestamps.END = datetime(eraseBetween(flight + " " + string(timestamps.END), 9, 12),'InputFormat','yyyyMMdd HH:mm:ss');

% overwrite orbit variable if best fit params have not been computed

is_params_computed = ~isempty(dir(pwd + delimiter + "output" + delimiter + flight + delimiter + "corr" + delimiter + "best-fit-params-ps-" + "set-" + sprintf('%02d', set) + ".mat"));
if ~is_params_computed
    warning('Best-fit parameters have not been computed for this set. Please allow current execution to finish and then re-run program.'); 
    orbit = nan;
else
    load(pwd + delimiter + "output" + delimiter + flight + delimiter + "corr" + delimiter + "best-fit-params-ps-" + "set-" + sprintf('%02d', set) + ".mat");
    load(pwd + delimiter + "output" + delimiter + flight + delimiter + "corr" + delimiter + "best-fit-params-pd-" + "set-" + sprintf('%02d', set) + ".mat");
    params_ps = table2cell(params_ps); params_pd = table2cell(params_pd);
end

% remove data points according to time period specified

if isnan(orbit)
    period = find(timestamps.SET == set & isnan(timestamps.ORBIT));
else
    period = find(timestamps.SET == set & timestamps.ORBIT == orbit);
end
idx_start = find(data.TIME == timestamps.START(period));
idx_end = find(data.TIME == timestamps.END(period));
data = data(idx_start:idx_end,:);

%% extract parameters for analysis

disp('>>> Extracting parameters ... ');

t = datenum(data.TIME);
lat = data.lat; 
lon = data.lon; 
alt = data.alt;
p_s = [data.ps_pitot, data.ps_858, data.ps_hgppt2, data.ps_aero, data.ps_aimms]; 
p_d = [data.pd_pitot, data.pd_858, data.pd_hgppt2, data.pd_aero, data.pd_aimms];

% remove outliers from static and dynamic pressure measurements

no_pressure_instr = size(p_s, 2);

for i=1:no_pressure_instr
    [~, TF1] = rmoutliers(p_s(:,i),'mean'); [~, TF2] = rmoutliers(p_d(i,i),'mean');
    TF = or(TF1, TF2);
    data(TF, :) = []; t(TF) = []; lat(TF) = []; lon(TF) = []; alt(TF) = []; p_s(TF,:) = []; p_d(TF,:) = [];
end

clear TF TF1 TF2

%% compute linear regression model parameters

if ~is_params_computed
    disp('>>> Computing linear regression model parameters ... ');
    
    params_ps = cell(no_pressure_instr, no_pressure_instr); params_pd = cell(no_pressure_instr, no_pressure_instr);
    
    for i=1:no_pressure_instr
        for j=1:no_pressure_instr
            x = p_s(:,i); y = p_s(:,j);
            temp = polyfit(x, y, 1);
            params_ps{i,j} = [temp(1) temp(2)];
            
            x = p_d(:,i); y = p_d(:,j);
            temp = polyfit(x, y, 1);
            params_pd{i,j} = [temp(1) temp(2)];
        end
    end
    
    params_ps = cell2table(params_ps); params_pd = cell2table(params_pd);
    params_ps.Properties.VariableNames = {'pitot', '858', 'hgppt2', 'aero', 'aimms'}; params_pd.Properties.VariableNames = params_ps.Properties.VariableNames; 
    params_ps.Properties.RowNames = params_ps.Properties.VariableNames; params_pd.Properties.RowNames = params_pd.Properties.VariableNames;
    
    save(pwd + delimiter + "output" + delimiter + flight + delimiter + "corr" + delimiter + "best-fit-params-ps-" + "set-" + sprintf('%02d', set), 'params_ps');
    save(pwd + delimiter + "output" + delimiter + flight + delimiter + "corr" + delimiter + "best-fit-params-pd-" + "set-" + sprintf('%02d', set), 'params_pd');
    
    disp('>>> Program execution complete!');
    
    clear i j x y temp
    
    return
end

%% compute regression errors

disp('>>> Computing regression errors ...');

no_data_pts = length(t);
ps_reg_errors = zeros(no_data_pts, no_pressure_instr); pd_reg_errors = zeros(no_data_pts, no_pressure_instr);

for i=1:no_data_pts
    for j=1:no_pressure_instr
        psx_measured = p_s(i, setdiff(1:no_pressure_instr, j));
        psy_measured = p_s(i, j);
        temp = params_ps(j, setdiff(1:no_pressure_instr, j)); temp = cell2mat(temp); temp = transpose(reshape(temp, 2, []));
        psy_expected = temp(:, 1) .* psx_measured' + temp(:, 2);
        ps_reg_errors(i, j) = sum(abs(psy_measured - psy_expected));
        
        pdx_measured = p_d(i, setdiff(1:no_pressure_instr, j));
        pdy_measured = p_d(i, j);
        temp = params_pd(j, setdiff(1:no_pressure_instr, j)); temp = cell2mat(temp); temp = transpose(reshape(temp, 2, []));
        pdy_expected = temp(:, 1) .* pdx_measured' + temp(:, 2);
        pd_reg_errors(i, j) = sum(abs(pdy_measured - pdy_expected));
    end
end

clear psx_measured pdx_measured psy_measured pdy_measured psy_expected pdy_expected temp

%% plot flight paths

disp('>>> Plotting flight paths ...');

keys = {'time';
    'ps_pitot_err';
    'ps_858_err';
    'ps_hgppt2_err';
    'ps_aero_err';
    'ps_aimms_err';
    'pd_pitot_err';
    'pd_858_err';
    'pd_hgppt2_err';
    'pd_aero_err';
    'pd_aimms_err'};
values = {t;
    ps_reg_errors(:,1);
    ps_reg_errors(:,2);
    ps_reg_errors(:,3);
    ps_reg_errors(:,4);
    ps_reg_errors(:,5);
    pd_reg_errors(:,1)
    pd_reg_errors(:,2)
    pd_reg_errors(:,3)
    pd_reg_errors(:,4)
    pd_reg_errors(:,5)};

M = containers.Map(keys, values);

colormap_labels = ["Time (datenum)", "Regression Errors (mbar)"];

screen_size = get(0,'screensize'); w = screen_size(3); h = screen_size(4);
fig = uifigure('Name', 'Flight Paths', 'HandleVisibility', 'on');
fig.Position = [w/2 0 w/2 h];

ax1 = uiaxes('Parent',fig,'Position',[40 h/2+80 w/3+40 h/3]); grid(ax1,'on'); view(ax1,45,45); h1 = colorbar(ax1); h1.Label.String = colormap_labels(1); h1.FontSize = 10;
ax2 = uiaxes('Parent',fig,'Position',[40 140 w/3+40 h/3]); grid(ax2,'on'); view(ax2,45,45); h2 = colorbar(ax2); h2.Label.String = colormap_labels(2); h2.FontSize = 10;

p1 = scatter3(ax1, lat, lon, alt, 10, values{1}, 'filled');
p1.CData = values{1};
xlabel(ax1,'Lat'); ylabel(ax1,'Lon'); zlabel(ax1,'Alt (m)');
p2 = scatter3(ax2, lat, lon, alt, 10, ps_reg_errors(:,1), 'filled');
p2.CData = values{2};
xlabel(ax2,'Lat'); ylabel(ax2,'Lon'); zlabel(ax2,'Alt (m)');
h2.Limits = [min(ps_reg_errors,[],'all') max(ps_reg_errors,[],'all')];

dd1 = uidropdown(fig,...
    'Position',[2*w/5 2*h/3+80 160 40],...
    'Items',keys,...
    'Value','time',...
    'ValueChangedFcn',@(dd1,event) path_var_selection(dd1,p1,M,h1,colormap_labels,ps_reg_errors,pd_reg_errors));

dd2 = uidropdown(fig,...
    'Position',[2*w/5 h/4+80 160 40],...
    'Items',keys,...
    'Value','ps_pitot_err',...
    'ValueChangedFcn',@(dd2,event) path_var_selection(dd2,p2,M,h2,colormap_labels,ps_reg_errors,pd_reg_errors));

%% plot static and dynamic pressures

disp('>>> Plotting static and dynamic pressures ...');

fig = figure;
fig.Position = [w/2 0 w/2 h];
subplot(2,1,1);
plot(data.TIME, p_s, 'LineWidth', 1); grid on
xlabel('Time (HH:mm)'); ylabel('Static Pressure (mbar)');
legend("pitot", "858", "hgppt2", "aero", "aimms");

subplot(2,1,2);
plot(data.TIME, p_d, 'LineWidth', 1); grid on
xlabel('Time (HH:mm)'); ylabel('Dynamic Pressure (mbar)');
legend("pitot", "858", "hgppt2", "aero", "aimms");

saveas(gcf, '.\output\' + flight + '\png\pressures-set-' + num2str(set,'%02.f') + "-orbit-" + num2str(orbit,'%02.f') + '.png');

%% plot correlation matrices

disp('>>> Plotting correlation matrices ...');

fig = figure;
fig.Position = [w/2 0 w/2 h];
corrplot(p_s, 'varNames', {'pitot', '858', 'hgppt', 'aero', 'aimms'});
saveas(gcf, '.\output\' + flight + '\png\static-pressures-corr-set-' + num2str(set,'%02.f') + "-orbit-" + num2str(orbit,'%02.f') + '.png');

fig = figure;
fig.Position = [w/2 0 w/2 h];
corrplot(p_d, 'varNames', {'pitot', '858', 'hgppt', 'aero', 'aimms'});
saveas(gcf, '.\output\' + flight + '\png\dynamic-pressures-corr-set-' + num2str(set,'%02.f') + "-orbit-" + num2str(orbit,'%02.f') + '.png');
        
disp('>>> Program execution complete!');

clear i j

%% function definitions

function path_var_selection(dd, p, M, h, colormap_labels, ps_reg_errors, pd_reg_errors)
    val = M(dd.Value);
    p.CData = val;
    if string(dd.Value) == "time"
        h.Label.String = colormap_labels(1);
        h.Limits = [min(M('time'),[],'all') max(M('time'),[],'all')];
        h.YTickLabel = extractBefore(h.YTickLabel,6);
    elseif strncmpi(string(dd.Value), "ps", 2)
        h.Label.String = colormap_labels(2);
        h.Limits = [min(ps_reg_errors,[],'all') max(ps_reg_errors,[],'all')];
    elseif strncmpi(string(dd.Value), "pd", 2)
        h.Label.String = colormap_labels(2);
        h.Limits = [min(pd_reg_errors,[],'all') max(pd_reg_errors,[],'all')];
    end
end
