%Loading Anodic and Cathodic mat folder

tld = 'C:\Users\arrio\Box\BensmaiaLab\UserData\UserFolders\ToriArriola\DARPA_updated\ProcessedData\Pinot\Cathodic_Anodic';


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


%will probably have to go over!!!!



for d = 1:size(data,2)
    for j = 1:size(data(d).DetectionTable,1)
        for n = 1:size(data(d).DprimeTable,1)
            x = str2double(data(d).DetectionTable.Properties.RowNames);
             varnames = data(d).CoeffTable.Properties.VariableNames;
             xq = linspace(0, 0.3);
             %without icms
             yq = sigfun(data(d).CoeffTable.dPrime{1}, xq);
             %with icms
             aq = sigfun(data(d).CoeffTable.dPrime{2}, xq);
             [~,b] = min(abs(yq-dprime_threshold));
             [~,c] = min(abs(aq-dprime_threshold));
           % plot out thresholds to make sure its correct
            hold on
             plot(xq, yq, 'LineWidth', 2)

            plot([0 xq(b) xq(b)], [dprime_threshold, dprime_threshold -1], 'LineStyle','--')
            plot([0 xq(c) xq(c)], [dprime_threshold, dprime_threshold -1], 'LineStyle','--')
            
            % previously mechthreshold
            data(d).woicmsthreshold = xq(b);

            data(d).wicmsthreshold = xq(c);
     
            
        end
    end
    
end

%% delta anodic - cathodic


electrodes_1 = vertcat(data(:).Electrodes);
electrodes = unique(electrodes_1, 'rows');


pulse_data = vertcat(data(:).Pulse);

cath_idx = strcmpi(pulse_data, 'Cathodic');
an_idx = strcmpi(pulse_data, 'Anodic');


for e = 1:size(electrodes, 1)
       
   % In data, find which rows contain those electrodes
     e_idx = ismember(vertcat(data(:).Electrodes),electrodes(e, :), "rows");

    %  %combining the cathodic trials and the same electrodes
    %  cath_e_idx = e_idx & cath_idx ;
    %  an_e_idx = e_idx & an_idx;
    %  delta_cath(e) = data(cath_e_idx).mechthreshold - data(an_e_idx).mechthreshold;
    %  delta_an(e) = data(an_e_idx).mechthreshold - data(cath_e_idx).mechthreshold;
    % delta_an_sorted = sort(delta_an);
    % delta_cath_sorted = sort(delta_cath);
end

%% Analysis of w/icms and w/o icms thresholds within conditions

for t = 1:length(data)
    data(t).diff = data(t).wicmsthreshold - data(t).woicmsthreshold;
  
    %need to rewrote to not be hard coded
    cath_idx_1 = [1, 4, 6, 8, 10];
    an_idx_1 = [2, 3, 5, 7, 9];
  
    delta_cath = [data(cath_idx_1).diff];
    delta_an = [data(an_idx_1).diff];
end    
   
%% Plotting anodic and cathodic deltas

hold on
ax = gca;
ax.FontSize = 15;
sz = 100;
scatter(delta_an, delta_cath,sz, 'd', 'filled','LineWidth',1.5)
plot([0,0], [-0.04 0.04] , 'Color', [.6 .6 .6], 'LineStyle', '--')
plot( [-0.04 0.04],[0,0] , 'Color', [.6 .6 .6], 'LineStyle', '--')
ylabel('Delta Cathodic Threshold')
xlabel('Delta Anodic Threshold')
xlim([-.04 .04])
ylim([-0.04 .04])
axis square;
set(gcf, 'position', [5, 5, 400, 400])

%% plot of thresholds with and witout icms for cathodic condition

w_o_icms = vertcat(data(cath_idx).woicmsthreshold);
w_icms = vertcat(data(cath_idx).wicmsthreshold);

hold on
scatter(w_o_icms,w_icms,sz, 'd', 'filled','LineWidth',1.5)
title('Cathodic Thresholds')
xlabel('Without ICMS(mm)')
ylabel('With ICMS(mm)')
xlim([0 .1])
ylim([0 .1])

%% plotting psychometric curve for darpa
%electrodes 22 and 24, cathodic
box off
dprime_threshold = 1.35;

SetFont('Arial', 18)
sgt = sgtitle('Electrode 22 and 24');
sgt.FontSize = 25;
subplot(1,2,1); hold on 


    mechamps = str2double(data(6).DetectionTable.Properties.RowNames);
   %this can be reduced to one line
    plot(mechamps, data(6).DetectionTable{:,1}, 'o-', 'Color', rgb(66, 66, 66),'LineWidth', 4)
    plot(mechamps, data(6).DetectionTable{:,2}, 'o-', 'Color', rgb(198, 40, 40),'LineWidth', 4)
    text(.03, .37, 'Without ICMS', 'Color', rgb(66, 66, 66),'FontSize',18);
    text(.03, .4, 'With ICMS', 'Color', rgb(198, 40, 40), 'FontSize',18);
    
    ylabel('pDetect')
    xlabel('Stimulus Amplitude (mm)')
   

subplot(1,2,2); hold on


   plot(mechamps, data(6).DprimeTable{:,1}, 'o-', 'Color', rgb(66, 66, 66),'LineWidth', 4)
   plot(mechamps, data(6).DprimeTable{:,2}, 'o-', 'Color', rgb(198, 40, 40),'LineWidth', 4)
   yline(1.35,'-', 'Threshold', 'FontSize',18);
   text( .03, .9,'0.019', 'Color', rgb(66, 66, 66),'FontSize',18);
   text(.03, 1, '0.017', 'Color', rgb(198, 40, 40), 'FontSize',18);
   ylabel('d''')
   xlabel('Stimulus Amplitude (mm)')

%%
subplot(1,2,1); hold on

axis square
w_o_icms = vertcat(data(cath_idx).woicmsthreshold);
w_icms = vertcat(data(cath_idx).wicmsthreshold);

hold on

scatter(w_icms,w_o_icms,sz, 'filled','LineWidth',1.5)

title('Cathodic Thresholds')
ylabel('Without ICMS(mm)')
xlabel('With ICMS(mm)')
xlim([0 .09])
ylim([0 .09])
plot(xlim,ylim,'Color', [.8 .8 .8], 'LineStyle','--')
xticks(0:0.02:0.08)
yticks(0:0.02:0.08)

% x1 = 0:.1:0.01;
% y1= x1;
% plot(x1, y1)

% .01 .02 .03 .04 .05 .06 .07 .08 .09 .1

subplot(1,2,2); hold on

axis square

chuck=w_icms./w_o_icms;

both = mean(chuck);

Swarm(.6,chuck, 'DS', 'Box', 'ShowStats',true,'SwarmMarkerSize', 60)
set(gca, 'xtick', [])
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation','origin')
ylabel('T(ICMS) / T(Con)')
% text(.6, 1.03,'Group "0.6": Median (P(25), P(75)) = 0.71 (0.42, 0.86)', 'Color', rgb(66, 66, 66),'FontSize',12);


ylim([0 1.01])
xlim([.2 1])

