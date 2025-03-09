% prepare structure of data (filter and segment data) for Oddball Protocol
% Author:Pegah Chehrazi. If any problem contact:pegah.chehrazi@gmail.com

env = envInit();
% It works for all, short and long triggers

%%%%%%%%%%%%%%%%%%%%%%% PATAMETERS TO MODIFY/CHECK %%%%%%%%%%%%%%%%%%%%%%%%%

matDi = dir('/E:/JK0/*'); % ****  DATA DIRECTORY  **** % Get a list of all folders (sessions).
% follow the same path organization /home/gaba/Desktop/GABA/DATA/JK*'



ID= 'JK0';% Project ID
ID2='all';

Areas='AUD';%, Select the areas to analyze: {'Visual'},{'PPC'},{'AUD'},{'HIP'}

typeDATAperiod=[{'DATAperiodsS1'},{'DATAperiodsS5'}]; % Define the sounds for which you want to analyze
group.wtNN=[{'JK0025-04'},{'JK0018-01'},{'JK0015-02'},{'JK0019-04'},{'JK0023-01'}];%5-10khz
group.KONN=[{'JK0016-06'},{'JK0017-02'},{'JK0014-03'},{'JK0024-03'},{'JK0026-03'},{'JK0022-04'}]; %on the ball
% group.wtNN=[{'AE0020-03'},{'AE0020-04'},{'AE0021-03'},{'AE0021-04'}];% 10-12khz





fn=fieldnames(group); %groups names
typeDG='OddBall1';%strcat Select the typeof Oddball protocol you want to analyze: 'OddBall1';'OddBall2'...
typeDGL=str2double(typeDG(end));
xlima=[-0.1 0.5];% X axis limits to plot
ylima=[-150 70];% Y axis limits to plot
bl=-0.1;% Baseline for fieldtrip BL correction to use recomended:-0.15 (use -0.015 when the baseline is crappy)
tpv=1;% threshold for p-value in the plots of traces (1=0.05 and 2=0.01)
SmoothWin=2;% Use window to smooth (max 60s)


arCh=Areas(1:2);
if all(arCh=='HI')
    ChaStr='Rad';
    ChaStr2='Pyr';
end
if all(arCh=='Vi') || all(arCh=='PP') || all(arCh=='M1')
    ChaStr='L5';
    ChaStr2='L4';
end
if all(arCh=='AU')
    ChaStr='L5';
end
if all(arCh=='PF')
    ChaStr='PL';
end
for g=1:size(fn,1)% each group
    for k=1:length(group.(fn{g})) %each animal
        animal=group.(fn{g}){k}
        if typeDGL==1
            nameprocesfile='_Process_OddBall.mat';%
        end
        
        matD =  dir(strcat(matDi(1).folder,'/',group.(fn{g}){1,k},'/**/',group.(fn{g}){1,k},'*',nameprocesfile));
        load(strcat(matD.folder,'/',matD.name));
        St=5;
        Row=1
        
        
        cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
        ChIsL=cellfun(cellfind(ChaStr),Data.(Areas));
        [rowCh,colCh]=find(double(ChIsL));
        
        
        fieldNames = fieldnames(Data.AUD{5, 2});
        ERPMinpeak_sound1 = [];
        ERPMinpeak_sound2 = [];
        ERPMinpeak_sound3 = [];
        ERPMinpeak_sound4 = [];
        ERPMinpeak_sound5 = [];
        ERPMinpeak_sound6 = [];
        ERPMinpeak_sound7 = [];
        ERPMinpeak_sound8 = [];
        ERPMinpeak_sound9 = [];
        ERPMinpeak_sound10 = [];
        numSounds=10
        for soundIndex=1:length(fieldNames);
            
            allERPMinpeak=[]
            allMinTime=[]
            allStdERP=[]
            
            for h=1:length(Data.(Areas){5,2}.(fieldNames{soundIndex}));
                dataCh3=[Data.(Areas){St,colCh}.(fieldNames{soundIndex}){1,h}];
                dataAAA(:,h)=dataCh3';
                namevar=[group.(fn{g}){1,k},'-',Data.(Areas){Row,colCh},'-',fieldNames{soundIndex}];
                
            end
            MeanValue=mean(dataAAA,2);
            
            
        end
        
        
        
        figure(1), clf
        subplot(211);
        plot(Data.timeline,MeanValue','b','linew',1);
        set(gca, 'xlim', [-500, 500]);
        %             xticks([-500:10:500]); % This sets ticks from -500 to 500 with a step of 100
        %             xticklabels([]);
        xline(0, 'LineStyle', '--', 'Color', 'k');
        title(['ERP from channel','-',namevar]);
        ylabel('Voltage (\muV)');
        
        % Set colormap and axis labels
        
        subplot(212);
        imagesc(Data.timeline, [], dataAAA');
        xline(0, 'LineStyle', '--', 'Color', 'k');
        set(gca,'clim', [mean(min(dataAAA)), mean(max(dataAAA))]);
        set(gca, 'xlim', [-500, 500]);
        xticks([-500:50:500]); % This sets ticks from -500 to 500 with a step of 100
        %             xticklabels([]);
        xline(0, 'LineStyle', '--', 'Color', 'k');
        xlabel('Time');
        ylabel('Trials');
        title(['All the trial dynamics']);
        saveas(gcf,[env.outputDir,namevar,'RastorPlot.png'],'png');
        close
        
        
        dataAAAA=dataAAA'
        wantedWindow=[60,120]; %in ms
        [vali,Idx1i]=min(abs(Data.timeline-wantedWindow(1)));% Detect index when time is 60
        [val2i,Idx2i]=min(abs(Data.timeline-wantedWindow(2)));% Detect index when time is 120
        selectedData=dataAAAA(:,Idx1i:Idx2i);
        selectedDataTime=Data.timeline(Idx1i:Idx2i);
        [ERPMaxpeak,MaxIndx]=max(selectedData,[],2); %max(selectedData, [], 2) calculates the minimum across rows of the original matrix.
        MaxTime = selectedDataTime(MaxIndx);
        
        stdERP=std(ERPMaxpeak);
        % Define the threshold to filter out unwanted ERPMinpeak values
        threshold = +50; % in microvolts
        % Remove ERPMinpeak values that are more positive than the threshold
        validIdx = ERPMaxpeak > threshold; % Logical index to keep only values less than -10 microvolts
        ERPMaxpeak_filtered = ERPMaxpeak(validIdx);
        % Calculate the standard deviation for the filtered ERPMinpeak values
        stdERPFiltered = std(ERPMaxpeak_filtered);
        % Optional: also update MinIndx to keep only corresponding indices
       MaxIndx_filtered = MaxIndx(validIdx);
        % Quantify the failure rate: percentage of data points that were filtered out
        failureRate = (1 - (numel(ERPMaxpeak_filtered) / numel(ERPMaxpeak))) * 100;
        
        % Display the failure rate text in subplot(312)
        figure(3), clf;
        % Subplot for ERP Minpeaks for All Trials
        subplot(311);
        stem(ERPMaxpeak, 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Marker', 'o');
        title('ERP Maxpeaks for All Trials');
        % Calculate overall mean and standard deviation
        meanERP = mean(ERPMaxpeak);
        meanERPFiltered=mean(ERPMaxpeak_filtered);
        stdERPFiltered=std(ERPMaxpeak_filtered);
        
        % Define mean and standard deviation for both original and filtered data
        meanValues = [meanERP, meanERPFiltered];
        stdValues=[stdERP,stdERPFiltered]
        % Create bar plot with two bars side-by-side
        subplot(312);
        bar_handle = bar(meanValues, 'FaceColor', [0.2 0.6 0.5], 'LineWidth', 1);
        % Customize bar colors
        bar_handle.FaceColor = 'flat';
        bar_handle.CData(1, :) = [0.2 0.6 0.5]; % Original data color
        bar_handle.CData(2, :) = [0.8 0.2 0.2]; % Filtered data color
        hold on;
        % Add error bars for both means
        errorbar(1:2, meanValues, stdValues, 'k', 'LineStyle', 'none', 'LineWidth', 1);
        % Customize the plot
        xticks(1:2);
        xticklabels({'Original Data', 'Filtered Data'});
        xlim([0.5 2.5]); % Adjust x-axis limits
        ylabel('Mean ERP Minpeak');
        title('Mean and Standard Deviation of ERP Minpeak for Original and Filtered Data');
        % Display the failure rate in the subplot
        text(1.5, max(meanValues + stdValues) * 1.1, ...
            sprintf('Failure Rate: %.2f%%', failureRate), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'Bottom',  'FontSize', 10, 'Color', 'k');
        hold off;
        % Subplot for filtered ERP Minpeaks
        subplot(313);
        stem(ERPMaxpeak_filtered, 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Marker', 'o');
        title('Filtered ERP Minpeaks for All Trials');
        hold off;
        % Save the figure with the failure rate text displayed
        saveas(gcf, [env.outputDir, namevar, 'StdvTrials', '.png'], 'png');
        % Close the figure if needed
        close;
        % Calculate the timing of the minimum value for each of the 60 trials
            FilteredMaxTime=selectedDataTime(MaxIndx_filtered);
        % Annotate each point with its corresponding trial number, adding an offset
        xOffset = 0.5;  % Adjust as needed for x-direction spacing
        yOffset = 1.0;  % Adjust as needed for y-direction spacing
   
            % Assuming MinTime and ERPMinpeak are both 1x60 vectors
            % MinTime: time at which the minimum value occurred for each trial
            % ERPMinpeak: the minimum value for each trial
            
            figure(3), clf; % Create a new figure and clear previous content
            % Plot all trials in light gray
            scatter(MaxTime, ERPMaxpeak, 60, [0.5, 0.5, 0.5], 'o', 'DisplayName', 'All Trials'); % Light gray markers
            hold on;
            % Overlay filtered trials with filled blue markers
            scatter(FilteredMaxTime, ERPMaxpeak_filtered, 20, 'r', 'o', 'DisplayName', 'Filtered Trials'); % Blue filled markers
            % Highlight filtered out trials with 'X' markers
            threshold = +50; % Example threshold
            filteredOutIdx = ERPMaxpeak <= threshold; % Find indices of trials not filtered out
            scatter(MaxTime(filteredOutIdx), ERPMaxpeak(filteredOutIdx), 60, 'r', 'X', 'DisplayName', 'Filtered Out Trials'); % Red X markers
            % Annotate filtered points
            for i = 1:length(FilteredMaxTime)
                text(FilteredMaxTime(i), ERPMaxpeak_filtered(i), sprintf('%d', i), ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
                    'FontSize', 8, 'Color', 'k');
            end
            % Customize the plot
            xlabel('Time of ERP Values (ms)'); % Adjust units if necessary
            ylabel('ERP MAX Value (\muV)');
            title('Timing of ERP Values with Filtered Points Marked');
            grid on;
            % Adjust the axis limits
            ylim([min(ERPMaxpeak) - 5, max(ERPMaxpeak) + 5]);
            xlim([min(MaxTime) - 5, max(MaxTime) + 5]);
            hold off;
            % Optionally add a legend
            legend('All Trials', 'Filtered Trials', 'Filtered Out Trials', 'Location', 'best');
            % Save the figure
            saveas(gcf, [env.outputDir, namevar, 'MaxTime', '.png'], 'png');
            close;
    end
            clear  dataAAA selectedDataTime selectedData stdERP ERPMinpeak
            close
        end
    
     % ERPMinpeak_sound1, ERPMinpeak_sound2, ERPMinpeak_sound3
        soundData = {ERPMinpeak_sound1, ERPMinpeak_sound2, ERPMinpeak_sound3, ERPMinpeak_sound4, ERPMinpeak_sound5, ERPMinpeak_sound6, ERPMinpeak_sound7, ERPMinpeak_sound8, ERPMinpeak_sound9, ERPMinpeak_sound10};
        soundLabels = {'Sound 1', 'Sound 2', 'Sound 3','Sound 4','Sound 5','Sound 6','Sound 7''Sound 8','Sound 9','Sound SD'};  % Labels for the different sounds
        % Create a new figure for the histograms
        figure(3), clf;
        hold on;
        % Define colors for each sound for better visualization
        colors = [
            1, 0, 0;     % Red
            0, 1, 0;     % Green
            0, 0, 1;     % Blue
            1, 0.5, 0;   % Orange
            0.5, 0, 0.5; % Purple
            0, 1, 1;     % Cyan
            1, 1, 0;     % Yellow
            0.5, 0.5, 0; % Olive
            0, 0.5, 1;   % Light Blue
            1, 0, 1;     % Magenta
            ];
        % Compute mean and standard deviation for each group
        meanERP = zeros(1, 10);
        stdERP = zeros(1, 10);
        % Plot individual data points using scatter for each group
        for i = 1:10
            soundDataGroup = soundData{i}; % Access sound data using cell indexing
            scatter(repmat(i, size(soundDataGroup)), soundDataGroup, 60, 'MarkerEdgeColor', colors(i, :), 'jitter', 'on', 'jitterAmount', 0.15);
            meanERP(i) = mean(soundDataGroup);
            stdERP(i) = std(soundDataGroup);
        end
        % Create a bar plot for the means
        bar_handle = bar(meanERP, 'FaceColor', 'flat'); % Create bar plot
        % Set different colors for each bar
        for i = 1:10
            bar_handle.CData(i, :) = colors(i, :); % Assign colors to each bar
        end
        % Add error bars to represent the standard deviation
        errorbar(1:10, meanERP, stdERP, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        % Customize the plot
        xticks(1:10);
        xticklabels({'Sound 1', 'Sound 2', 'Sound 3', 'Sound 4', 'Sound 5', ...
            'Sound 6', 'Sound 7', 'Sound 8', 'Sound 9', 'Sound 10'});
        xlabel('Sound');
        ylabel('Mean ERP Minpeak');
        title('Distribution, Mean and Standard Deviation of ERP Minpeak for Different Sounds');
        legend({'Sound 1', 'Sound 2', 'Sound 3', 'Sound 4', 'Sound 5', ...
            'Sound 6', 'Sound 7', 'Sound 8', 'Sound 9', 'Sound 10'}, ...
            'Location', 'northeast');
        hold off;
        % Save the figure if needed
        saveas(gcf, [env.outputDir, namevar, 'DistributionOfERPMinpeaks', '.png'], 'png');
     
        
        
        
        clear  ERPMinpeak_sound1 ERPMinpeak_sound2 ERPMinpeak_sound3 ERPMinpeak_sound4 ERPMinpeak_sound5 ERPMinpeak_sound6 ERPMinpeak_sound7 ERPMinpeak_sound8 ERPMinpeak_sound9 ERPMinpeak_soundSD
        
    
    clearvars dataAAA -except matDi env k ID pathDATA Overwrite sampleTimeBef sampleTime ParamContainer p Fs sampleTimeBef sampleTimeBefMs ...
        sampleTimeBefSam sampleTime filter Areas lengthPer sampleTimeMs sampleTimeSam
