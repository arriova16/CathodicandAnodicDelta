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
    
 end
    
%% mech threshold

sigfun = @(c,x) (c(3) .* (1./(1 + exp(-c(1).*(x-c(2)))))) + c(4);
dprime_threshold = 1.35;

for d = 1:size(data,2)
    for j = 1:size(data(d).DetectionTable,1)
        for n = 1:size(data(d).DprimeTable,1)
            x = str2double(data(d).DetectionTable.Properties.RowNames);
            varnames = data(d).CoeffTable.Properties.VariableNames;
            y = data(i).DprimeTable{:,:};   
            % xq = linspace(x(1), x(end +.05);
            xq = linspace(x(1), x(end));
            yq = sigfun(data(d).CoeffTable.dPrime{1}, xq);
            [~,b] = min(abs(yq-dprime_threshold));
            % plot out thresholds to make sure its correct
            hold on
             plot(xq, yq, 'LineWidth', 2)

            plot([0 xq(b) xq(b)], [dprime_threshold, dprime_threshold -1], 'LineStyle','--')

            data(d).mechthreshold = xq(b);
        
     
            
        end
    end
    
end

%% anodic - cathodic
%trouble getting started with logic
% 42 and 31 Anodic not reaching threshold incorrect
%need to change constraints? ask charles

for i = 1:length(data)


end
