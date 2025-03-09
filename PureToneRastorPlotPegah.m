% prepare structure of data (filter and segment data) for Pure Tone Protocol
% Author:Maria Isabel Carreno-Munoz. If any problem contact:macamu1988@gmail.com

env = envInit();

matDi = dir('/media/gaba/MARIA/DATA/SYN/'); % ****  DATA DIRECTORY  **** % Get a list of all folders (sessions).

ID= 'SYN';% Project ID
ID2='All';%PV
Areas='AUD';%, Select the areas to analyze: {'Visual'},{'PPC'},{'AUD'},{'HIP'}
typeAUD='5KHz';
ChaStr='L5';
Areas='AUD'
St=5;
resultTable = table();

%load the '_Process2_PT.mat'

%SST
% group.wtNN=[{'JK0004-01'},{'JK0004-03'},{'JK0007-04'},{'JK0005-01'},{'JK0005-02'},{'JK0005-03'},{'JK0001-01'},{'JK0001-02'},{'JK0001-03'},{'JK0009-01'},{'JK0009-02'},{'JK0009-03'}];%,{'JK0004-02'},
% group.KONN=[{'JK0011-01'},{'JK0011-03'},{'JK0011-04'},{'JK0002-01'},{'JK0002-02'},{'JK0002-03'},{'JK0003-01'},{'JK0003-02'},{'JK0003-03'},{'JK0003-06'},{'JK0006-02'},{'JK0006-03'},{'JK0008-01'},{'JK0008-02'},{'JK0010-01'},{'JK0010-02'},{'JK0010-03'},{'JK0012-01'},{'JK0012-02'}];

% %PV
% group.wtNN=[{'JK0026-01'},{'JK0026-04'},{'JK0026-05'},{'JK0022-01'},{'JK0022-02'},{'JK0022-03'},{'JK0022-04'},{'JK0024-01'},{'JK0024-02'},{'JK0024-03'},{'JK0024-04'},{'JK0017-01'},{'JK0017-02'},{'JK0016-01'},{'JK0016-02'},{'JK0016-04'},{'JK0014-01'},{'JK0014-02'},{'JK0014-03'},{'JK0014-04'}];%{'JK0017-03'},{'JK0017-04'},{'JK0026-03'},
% group.KONN=[{'JK0025-01'},{'JK0025-02'},{'JK0025-03'},{'JK0023-01'},{'JK0023-02'},{'JK0023-03'},{'JK0023-04'},{'JK0015-01'},{'JK0015-03'},{'JK0015-04'},{'JK0019-01'},{'JK0019-02'},{'JK0019-03'},{'JK0019-04'},{'JK0018-01'},{'JK0018-02'},{'JK0018-03'},{'JK0018-04'}]; %,on the ball,{'JK0015-02'},


% group.w19NN=[{'AE0023-01'}]%,{'AE0019-05'},{'AE0019-06'},{'AE0019-07'},{'AE0019-08'},{'AE0019-09'},{'AE0019-10'},{'AE0019-11'},{'AE0019-30'},{'AE0019-31'},{'AE0019-32'},{'AE0019-33'},{'AE0019-34'},{'AE0019-35'}]%
 group.wtNN=[{'SYN077-19'}]%,{'AE0019-05'},{'AE0019-06'},{'AE0019-07'},{'AE0019-08'},{'AE0019-09'},{'AE0019-10'},{'AE0019-11'},{'AE0019-30'},{'AE0019-31'},{'AE0019-32'},{'AE0019-33'},{'AE0019-34'},{'AE0019-35'}]%


fn=fieldnames(group);
 FilteredMinTime=[]

    AllGroupFilteredMinTime = struct('Group', [], 'StimFreq', []);

for g=1:size(fn,1)% each group
    vfPeriodPzLfp=[];
    AllGroupFilteredMinTime(g).StimFreq = struct();
    AllGroupFilteredMinTime(g).Group = fn{g};


    for k=1:length(group.(fn{g}))
        if all(typeAUD(1:2)=='5K')
            nameprocesfile='_Process2_PT.mat';%
            typeDATAperiod=['DATAperiods5KHz'];
            typeAUDSS=typeDATAperiod(12:end);%'S1';
        end
        if all(typeAUD(1:2)=='10')
            nameprocesfile='_Process2_PT.mat';%
            typeDATAperiod=['DATAperiods10KHz'];
            typeAUDSS=typeDATAperiod(12:end);%
        end

        matD =  dir(strcat(matDi(1).folder,'/', group.(fn{g}){1,k},'/**/','*',nameprocesfile));
        load(strcat(matD.folder,'/',matD.name));
        ParamContainer=strcat(matD.folder,'/',matD.name,'/**/EXP1*.mat');
        filename=matD.name;
        directoryname = matD.folder;
        filename = strrep(filename, '_Process_PT.mat', '');
        Strfilename=strcat(directoryname,'/',filename,'_Process2_PT.mat');
        cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
        ChIsL=cellfun(cellfind(ChaStr),Data.(Areas));
        [rowCh,colCh]=find(double(ChIsL));

        

        fieldNames = fieldnames(Data.AUD{5, 2});
        for j=1:length(fieldNames);

            fieldName = string(fieldNames{j});  % Convert field name to string
            allERPMinpeak=[];
            allStdERP=[]
           
            namevar=fieldNames{j};
            for h=1:length(Data.(Areas){5,2}.(fieldNames{j}));
                dataCh3=[Data.(Areas){St,colCh}.(fieldNames{j}){1,h}'];
                dataAAA(:,h)=dataCh3;
            end
           
MatG=dataAAA';
            wantedWindow=[5,50];
            [vali,Idx1i]=min(abs(Data.timeline-wantedWindow(1)));% Detect index when time is 7
            [val2i,Idx2i]=min(abs(Data.timeline-wantedWindow(2)));% Detect index when time is 20
            selectdata=MatG(:,Idx1i:Idx2i);
            selecteddataTime=Data.timeline(Idx1i:Idx2i);
            [ERPMinpeak,MinIndx]=min(selectdata,[],2); %ERPpeak is the minimum value (ERP)
            MinTime=selecteddataTime(MinIndx);  % Calculate the timing of the minimum value for each of the 60 trials


SingleERPValues = MatG; % Copy data to a new variable (if needed)
numRepetitions = size(MatG, 1); % Number of repetitions (49)

for i = 1:numRepetitions
    figure; % Create a new figure for each repetition
    plot(Data.timeline, SingleERPValues(i, :), 'b', 'LineWidth', 1); % Plot the i-th repetition
    set(gca, 'XLim', [-500, 500]); % Set x-axis limits
    xline(0, 'LineStyle', '--', 'Color', 'k'); % Add a vertical line at 0
    ylabel('Voltage (\muV)'); % Label the y-axis
    xlabel('Time (ms)'); % Label the x-axis
saveas(gcf, fullfile(env.outputDir, ID, namevar, ...
    [group.(fn{g}){k}(1:9), '_ERPMinpeak_', num2str(i), '.png']), 'png');
end
            

            figure(1),clf;
            plot(1:length(ERPMinpeak), ERPMinpeak, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
            xticks(1:50); % Set x-axis ticks to match each point (1 to 50)
            xlabel('Trigger Repetitions'); % Label for x-axis
            ylabel('ERP Min Peak'); % Label for y-axis
            title(['ERPMinpeak Range', group.(fn{g}){k}(1:9),fieldName]);
            saveas(gcf,[env.outputDir,ID2,'/',namevar,'/',group.(fn{g}){k}(1:9),'ERPMinpeak Range','.png'],'png');
            %            userInput = input('Press Enter to continue...', 's');
            close


            stdERP=std(ERPMinpeak);%time at which we have the min value (ERP)
            % Define the threshold to filter out unwanted ERPMinpeak values
            threshold = -10; % in microvolts
            % Remove ERPMinpeak values that are more positive than the threshold
            validIdx = ERPMinpeak < threshold; % Logical index to keep only values less than -10 microvolts
            ERPMinpeak_filtered = ERPMinpeak(validIdx);
            % Calculate the standard deviation for the filtered ERPMinpeak values
            stdERPFiltered = std(ERPMinpeak_filtered);
            % Optional: also update MinIndx to keep only corresponding indices
            MinIndx_filtered = MinIndx(validIdx);
            % Quantify the failure rate: percentage of data points that were filtered out
            failureRate = (1 - (numel(ERPMinpeak_filtered) / numel(ERPMinpeak))) * 100;
            animalIdentity=group.(fn{g}){k}(1:9);
            animalIdentity = string(animalIdentity);  % Convert to string
            % Extract the current field name
            currentRow = table(animalIdentity,  fieldName, failureRate, 'VariableNames', {'AnimalIdentity', 'FieldName', 'FailureRate'});
            resultTable = [resultTable; currentRow];  % Add new row to the existing table
            % Write the table to the Excel file
            FilteredMinTime=selecteddataTime(MinIndx_filtered);


            figure(2), clf;
            subplot(311);
            stem(ERPMinpeak, 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Marker', 'o');
            title(['ERP Minpeaks/All Trials',group.(fn{g}){k}(1:9)]);
            % Calculate overall mean and standard deviation
            meanERP = mean(ERPMinpeak);
            meanERPFiltered = mean(ERPMinpeak_filtered);
            % Define mean and standard deviation for both original and filtered data
            meanValues = [meanERP, meanERPFiltered];
            stdValues = [stdERP, stdERPFiltered];
            % Discard data points outside ±2 SD for filtered ERP data
            validERPFiltered = ERPMinpeak_filtered(ERPMinpeak_filtered >= (meanERPFiltered - 2*stdERPFiltered) & ERPMinpeak_filtered <= (meanERPFiltered + 2*stdERPFiltered));

            subplot(312);
            stem(ERPMinpeak_filtered, 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Marker', 's');
            title(['Filtered Minpeaks All Trials',group.(fn{g}){k}(1:9),fieldName]);

            subplot(313);
            bar_handle = bar(meanValues, 'FaceColor', [0.2 0.6 0.5], 'LineWidth', 1);
            % Customize bar colors
            bar_handle.FaceColor = 'flat';
            bar_handle.CData(1, :) = [0.2 0.6 0.5]; % This sets the color of the first bar to a custom RGB color [0.2 0.6 0.5], a shade of teal.
            bar_handle.CData(2, :) = [0.8 0.2 0.2]; % This sets the color of the second bar to [0.8 0.2 0.2], a reddish shade.
            hold on;
            % Add error bars for both means
            errorbar(1:2, meanValues, stdValues, 'k', 'LineStyle', 'none', 'LineWidth', 1);
            % Customize the plot
            xticks(1:2);
            xticklabels({'Original Data', 'Filtered Data'});
            xlim([0.5 2.5]); % Adjust x-axis limits
            ylabel('Mean ERP Minpeak');
            title(['Mean and SDV Minpeak Original and Filtered Data',group.(fn{g}){k}(1:9)]);
                        saveas(gcf, [env.outputDir,ID2,'/',namevar,'/',group.(fn{g}){k}(1:9),'StdvTrials', '.png'], 'png');
            hold off;
            close;

            figure(4), clf; % Create a new figure and clear previous content
            % Annotate each point with its corresponding trial number, adding an offset
            xOffset = 1.5;  % Adjust as needed for x-direction spacing
            yOffset = 1.0;  % Adjust as needed for y-direction spacing
            scatter(MinTime, ERPMinpeak, 40, [0.5, 0.5, 0.5], 'o', 'DisplayName', 'All Trials'); % Light gray markers
            hold on;
            % Overlay filtered trials with filled blue markers
            scatter(FilteredMinTime, ERPMinpeak_filtered, 80, 'r', 'o', 'DisplayName', 'Filtered Trials'); % Blue filled markers
            % Highlight filtered out trials with 'X' markers
            threshold = -10; % Example threshold
            filteredOutIdx = ERPMinpeak >= threshold; % Find indices of trials not filtered out
            scatter(MinTime(filteredOutIdx), ERPMinpeak(filteredOutIdx), 80, 'k', 'X', 'DisplayName', 'Filtered Out Trials'); % Red X markers
            % Annotate filtered pointsmatD
            for j = 1:length(ERPMinpeak)
                text(MinTime(j), ERPMinpeak(j), sprintf('%d', j), ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
                    'FontSize', 8, 'Color', 'k');
            end
            % Customize the plot
            xlabel('Time of ERP (ms)'); % Adjust units if necessary
            ylabel('ERP Min(\muV)');
            grid on;
            % Adjust the axis limits
            ylim([min(ERPMinpeak) - 5, max(ERPMinpeak) + 5]);
            xlim([min(MinTime) - 5, max(MinTime) + 5]);
            hold off;
            % Optionally add a legend500
            legend('All Trials', 'Filtered Trials', 'Filtered Out Trials', 'Location', 'best');
            title(['ERP Minpeaks time',group.(fn{g}){k}(1:9),fieldName]);
            saveas(gcf, [env.outputDir,ID2,'/',namevar,'/',group.(fn{g}){k}(1:9),'MinTime', '.png'], 'png');
            %                         userInput = input('Press Enter to continue...', 's');

            close;

            figure(5), clf;

            subplot(221);
            MeanValue1=mean(dataAAA,2);
            plot(Data.timeline,MeanValue1','b','linew',1);
            set(gca, 'xlim', [-500, 500]);
            xline(0, 'LineStyle', '--', 'Color', 'k');
            title(['ERP',group.(fn{g}){k}(1:9),fieldName]);
            ylabel('Voltage (\muV)');
            
            subplot(222);
            DataMatrix=dataAAA';
            imagesc(Data.timeline', [], DataMatrix);
            colormap(jet);       % jet parula Set colormap
            colorbarHandle = colorbar; % Add colorbar
            caxis([-500, 500]); % Replace minValue and maxValue with your desired range
            xline(0, 'LineStyle', '--', 'Color', 'k');
            set(gca, 'xlim', [-500, 500]);
            xlabel('Data Point Index');
            ylabel('Trials');
            title(['All trial dynamics', group.(fn{g}){k}(1:9),fieldName]);
          
            subplot(223);
            FilteredData = dataAAA(:, validIdx');
            disp(size(FilteredData));  %
            MeanValue2=mean(FilteredData,2);
            plot(Data.timeline,MeanValue2','b','linew',1);
            set(gca, 'xlim', [-500, 500]);
            xline(0, 'LineStyle', '--', 'Color', 'k');
            title(['ERP from FilteredData',group.(fn{g}){k}(1:9),fieldName]);
            ylabel('Voltage (\muV)');

            subplot(224);
            imagesc(Data.timeline, [], FilteredData');
            colormap(jet);       % jet parula Set colormap
            colorbarHandle = colorbar; % Add colorbar            xline(0, 'LineStyle', '--', 'Color', 'k');
            %           set(gca,'clim', [mean(min(FilteredData)), mean(max(FilteredData))]);
            set(gca,'clim', [-500, 500]);
            set(gca, 'xlim', [-500, 500]);
            xline(0, 'LineStyle', '--', 'Color', 'k');
            xlabel('Time');
            ylabel('FilteredTrials');
            title(['Filtered trial dynamics',group.(fn{g}){k}(1:9),fieldName]);
            saveas(gcf,[env.outputDir,ID2,'/',namevar,'/',group.(fn{g}){k}(1:9),'FilteredRastorPlot.png'],'png');



                 % Assign FilteredMinTime to the corresponding frequency
        if ~isfield(AllGroupFilteredMinTime(g).StimFreq, fieldName)
            AllGroupFilteredMinTime(g).StimFreq.(fieldName) = [];
        end
            % Append data to the corresponding field
            AllGroupFilteredMinTime(g).StimFreq.(fieldName) = ...
            [AllGroupFilteredMinTime(g).StimFreq.(fieldName), FilteredMinTime];
            FilteredData=FilteredData';
            [amplitude_tags, amplitude_tagsA] = f_detect_artefact_ongoing_Maria(FilteredData,5,1);
            FilteredData(:,amplitude_tags)=[];
            %           FilteredDataStruct.FilteredData = FilteredData;
            Data.AUD{6, 1} ={ strcat('ThresholdStruct')};
            %           namevar=strcat('typeDATAperiod',fieldName);
            Data.AUD{6, colCh}.(namevar)= FilteredData;
            %             Strfilename=strcat(directoryname,'/',filename,'_Process2_PT.mat');
            Strfilename=strcat(directoryname,'/',filename,'_Process3_PT.mat');
            save(Strfilename,'Data');

        end
        clear FilteredMinTime
    end
    % clearvars -except env matDi g k Areas ID pathDATA Overwrite sampleTimeBef sampleTime sessionstoADD

    outputPath = fullfile(env.outputDir, 'FailureRate.xlsx');
    % Write the final result table to the Excel file
    writetable(resultTable, outputPath);
 



end


% Define the number of bins for the histogram
numBins = 70;

% Define the frequency fields for 5kHz and 10kHz
freqFields = {'DATAperiods5KHz', 'DATAperiods10KHz'};  % List the frequencies you want to analyze

% Define group colors and labels
groupColors = {[1, 0, 0], [0, 0, 1]};  % RGB values for red (Ctrl) and blue (cKO)
groupLabels = {'Ctrl', 'cKO'};  % Labels for legend

% Loop through the two frequency fields (5kHz and 10kHz)
for f = 1:length(freqFields)
    freqField = freqFields{f};  % Current frequency (e.g., '5kHz' or '10kHz')
    
    % Initialize arrays to store data separately for Ctrl and cKO groups
    ctrlData = [];
    ckoData = [];

    % Loop through all the animals to gather data for both groups
    for n = 1:length(AllGroupFilteredMinTime)
        % Check if the field for the specified frequency exists
        if isfield(AllGroupFilteredMinTime(n).StimFreq, freqField)
            % Append data to the appropriate group based on the group label
            if strcmp(AllGroupFilteredMinTime(n).Group, 'wtNN')
                ctrlData = [ctrlData, AllGroupFilteredMinTime(n).StimFreq.(freqField)];  % Ctrl group data
            elseif strcmp(AllGroupFilteredMinTime(n).Group, 'KONN')
                ckoData = [ckoData, AllGroupFilteredMinTime(n).StimFreq.(freqField)];  % cKO group data
            end
        end
    end
    
    % Create a new figure for each frequency
    figure(f), clf;
    hold on;
    
    % Plot histogram for Ctrl group for the current frequency
    if ~isempty(ctrlData)
        histogram(ctrlData, 'Normalization', 'pdf', ...
            'FaceColor', groupColors{1}, 'EdgeColor', 'k', 'FaceAlpha', 0.5, ...
            'NumBins', numBins, 'DisplayName', sprintf('%s Ctrl %s', groupLabels{1}, freqField));
    else
        warning('No data found for Ctrl group for frequency %s.', freqField);
    end

    % Plot histogram for cKO group for the current frequency
    if ~isempty(ckoData)
        histogram(ckoData, 'Normalization', 'pdf', ...
            'FaceColor', groupColors{2}, 'EdgeColor', 'k', 'FaceAlpha', 0.5, ...
            'NumBins', numBins, 'DisplayName', sprintf('%s cKO %s', groupLabels{2}, freqField));
    else
        warning('No data found for cKO group for frequency %s.', freqField);
    end

    % Add labels and title
    xlabel('Frequency (Hz)');
    ylabel('Probability Density');
    title(sprintf('Histogram of Filtered Time Data for %s', freqField));
    
    % Add legend
    legend('Location', 'best');
    
    hold off;
        saveas(gcf, [env.outputDir,ID2,'/',group.(fn{g}){k}(1:9), freqFields{f}, 'ERP TimeDistribution ', '.png'], 'png');
 close 
end


 % Assuming `freqField` is defined as 'DATAperiods5KHz' or 'DATAperiods10KHz'
for freqIdx = 1:2
    if freqIdx == 1
        freqField = 'DATAperiods5KHz';
        figureTitle = '5kHz';
    else
        freqField = 'DATAperiods10KHz';
        figureTitle = '10kHz';
    end

    figure(10 + freqIdx); % Separate figures for 5kHz and 10kHz
    hold on

    for m = 1:size(fn, 1) % Iterate over groups (Ctrl and cKO)
        % Retrieve data for the current group and frequency
        groupData = [];
        for n = 1:length(AllGroupFilteredMinTime)
            if strcmp(AllGroupFilteredMinTime(n).Group, groupLabels{m}) && ...
               isfield(AllGroupFilteredMinTime(n).StimFreq, freqField)
                groupData = [groupData, AllGroupFilteredMinTime(n).StimFreq.(freqField)];
            end
        end

        % Skip if no data for the current group and frequency
        if isempty(groupData)
            continue;
        end

        % Compute mean and standard deviation for Gaussian
        mu = mean(groupData);
        sigma = std(groupData);

        % Generate x values for Gaussian curve
        x = linspace(min(groupData), max(groupData), 100);

        % Gaussian formula
        gaussianCurve = (1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma).^2);

        % Scale Gaussian curve to match histogram y-axis
        binWidth = (max(groupData) - min(groupData)) / numBins;
        gaussianCurve = gaussianCurve / sum(gaussianCurve) * numel(groupData) * binWidth;

        % Plot Gaussian curve
        plotHandles(m) = plot(x, gaussianCurve, 'Color', groupColors{m}, 'LineWidth', 2, ...
            'DisplayName', sprintf('%s Gaussian', groupLabels{m}));
    end

    % Finalize plot
    legend('show');
    xlabel('Filtered Min Time');
    ylabel('Density');
    title(['Gaussian Fit of ERP Min Peaks - ', figureTitle]);
    hold off;

    % Save figure
    saveas(gcf, fullfile(env.outputDir, [figureTitle, '_Gaussian_Distribution.png']), 'png');
end

