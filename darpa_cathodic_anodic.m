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
            % plot out thresholds to make sure its correct
            hold on
             % plot(xq, yq, 'LineWidth', 2)

            % plot([0 xq(b) xq(b)], [dprime_threshold, dprime_threshold -1], 'LineStyle','--')

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
    e_idx = find(vertcat(data(:).Electrodes) == electrodes(e, :));

    
    % cath_e_idx = e_idx & cath_idx == rows which match in electrode and are cathodic
    % an_e_idx = e_idx & an_idx == rows which match in electrode and are anodic
        % deltcath(e) = data(cath_e_idx).mechthreshold - data(an_e_idx).mechthreshold;

   %weher
    % if strcmpi(data(i).Pulse, 'Cathodic')
    % 
    % 
    % end
end


%    if data(i).Pulse == 'Cathodic'
%        cath_idx = find(data(i).Pulse == 'Cathodic');
%    else
%        an_idx = find(data(i).Pulse == 'Anodic');
%    end
% %wrong try again
% 
%     if data(i).Electrodes == electrodes 
%         deltcath(e) = data(cath_e_idx).mechthreshold - data(an_e_idx).mechthreshold;
%     else
%         deltan = data(i).mechthreshold - data(cath_idx).mechthreshold;
%    end
% 

    % if data(i).Electrodes == electrodes 
    %      deltcath = data(cath_idx).mechthreshold - data(an_idx).mechthreshold;
    % end


 

%condition to do this only if electrodes are matching but pulses are
%different
% data(i).mechthreshold - data(i).mechthreshold
