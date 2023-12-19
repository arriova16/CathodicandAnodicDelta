%Loading Anodic and Cathodic mat folder

tld = 'B:\ProjectFolders\DARPA\Data\ProcessedData\Pinot\Cathodic_Anodic';


%% Loading mat files
mat_files = dir(fullfile(tld, '*.mat'));
data = struct();ii=1;

%going through mat files only
for b = 1:size(mat_files,1)
    fname_split = strsplit(mat_files(b).name,'_');
    % organize data struct by electrode 
    %have to only get numbers won't work with and in there
    %try to remove double response tables

    electrode_numbers = fname_split{2};
     if contains(fname_split{2}, 'and') % Pairs
        and_idx = strfind(fname_split{2}, 'and');
        ee = [str2double(fname_split{2}(1:and_idx-1)), str2double(fname_split{2}(and_idx+3:end))];
     end

    data(ii).Electrodes= ee;
     
     pulse_idx = string(fname_split{6}(1:2));
     if pulse_idx == "Ca"
         data(ii).Pulse = 'Cathodic';
     else 
         data(ii).Pulse = 'Anodic';
     end
     temp = load(fullfile(mat_files(b).folder, mat_files(b).name));
    
     data(ii).ResponseTable = temp.bigtable;

    ii=ii+1;

end
        
 %% Analysis 
 %need to get detection tables first %try to write function
 %then get the fit the curve using charles function


for i = 1:length(data)
    %need to create a new function in rewrite %include charles fitsigmoid
    [detection_table{i}, dprime_table{i}, coeff_table{i}] = AnalyzeHybridTable(data(i).ResponseTable);
    
    data(i).DetectionTable = detection_table{i};
    data(i).DprimeTable = dprime_table{i};
    data(i).CoeffTable = coeff_table{i};
       pulse_string = convertCharsToStrings(data(i).Pulse); 
   data(i).Pulse = pulse_string;
    
    
 end
    
%% mech threshold
%need to plot out both with and without icms
sigfun = @(c,x) (c(3) .* (1./(1 + exp(-c(1).*(x-c(2)))))) + c(4);
dprime_threshold = 1.35;

for d = 1:size(data,2)
    for j = 1:size(data(d).DetectionTable,1)
        for n = 1:size(data(d).DprimeTable,1)
            x = str2double(data(d).DetectionTable.Properties.RowNames);
             varnames = data(d).CoeffTable.Properties.VariableNames;
             y = data(i).DprimeTable{:,:}; 
             
             %changed because some points werent reaching threshold
              % xq = linspace(x(1), x(end));
             xq = linspace(0, 0.3);
             yq = sigfun(data(d).CoeffTable.dPrime{1}, xq);
             [~,b] = min(abs(yq-dprime_threshold));
            plot out thresholds to make sure its correct
            hold on
            plot(xq, yq, 'LineWidth', 2)

            plot([0 xq(b) xq(b)], [dprime_threshold, dprime_threshold -1], 'LineStyle','--')

            data(d).mechthreshold = xq(b);
        
     
            
        end
    end
    
end

%% delta anodic - cathodic

% each of the unique electrodes and if its cathodic and anodic subtract
%doing this outside of loop because it would only give me the unique
%electrodes of each row. not useful

electrodes_1 = vertcat(data(:).Electrodes);
electrodes = unique(electrodes_1, 'rows');


pulse_data = vertcat(data(:).Pulse);

cath_idx = strcmpi(pulse_data, 'Cathodic');
an_idx = strcmpi(pulse_data, 'Anodic');


for e = 1:size(electrodes, 1)
    
    
   % In data, find which rows contain those electrodes
   %doing find here messed it up and only went through rows then columns
   %one at a time instead of at the same time

     e_idx = ismember(vertcat(data(:).Electrodes),electrodes(e, :), "rows");

     %combining the cathodic trials and the same electrodes
     cath_e_idx = e_idx & cath_idx ;
     an_e_idx = e_idx & an_idx;
     delta_cath(e) = data(cath_e_idx).mechthreshold - data(an_e_idx).mechthreshold;
     delta_an(e) = data(an_e_idx).mechthreshold - data(cath_e_idx).mechthreshold;
    delta_an_sorted = sort(delta_an);
    delta_cath_sorted = sort(delta_cath);
end

%% PLotting 

hold on
ax = gca;
ax.FontSize = 15;
sz = 100;
scatter(delta_an, delta_cath,sz, 'd', 'filled','LineWidth',1.5)
ylabel('Cathodic-Anodic Threshold')
xlabel('Anodic-Cathodic Threshold')
xlim([-.2 .2])
ylim([-0.2 .2])

% change axis from 0- 0.2

% xaxis = delta threshold during anodic condition
% yaxis = delta threshold during cathodic condition
% the delta is the difference between icms and w/o icms in either cathodic
% or anodic