%% Section 0 - Script Explanation and Participant Declaration

% This script should be ran from the sgrahamlab/labdata/NPTF2F folder
% Files Required (in Additional Scripts Folder):
%   - NPTF2F_RemoveErrors.m
%   - NPTF2F_SpeedData.m
%   - NPTF2F_SignalsData.m
%   - getAverageForce.m
%   - getTotalDistance.m
%   - sigfilt1.m
%   - spikeRemoval.m
%   - zeroX.m
cd("/sgrahamlab/labdata/NPTF2F");
addpath("Additional Scripts");

%Add the participant IDs and their corresponding dates below
sub_list = {'846'};
date_list = {'2024-03-12'};

% Toggle between the following Sequence Orders for each participant
seqOrder = ['INI' 'EPI'];
% seqOrder = ['EPI' 'INI'];

% Choose which tests to analyze - It is recommended to do every test so do
% not change this.
subns = str2double(sub_list(:));     
seqns = ['1' '2'];    %['1' '2'] 
parts = ['A' 'B'];    %['A' 'B']
runns = [1 2];        %[1 2]
testns= [1 2];        %[1 2]

%Change the videoDisplayRate to watch videos faster or slower when
%assigning TMT Link Scores
videoDisplayRate = 200; 
%1000 for 5s intervals, 200 for 1s intervals, 20 for 0.1s intervals

%% Section 1 - ProcessData - Read Data from the Output .txt Files
% This section calls the NPTF2F_ReadFromOutputTXTFiles() function for each
% Run because the data for each Run is stored in a different file.
for i = 1:length(sub_list)
    for seqi = 1:2      
        for runi = 1:2
            %FORMAT: NPTF2F_ReadFromOutputTXTFiles(subn, seqn, date, runn, prep) 
            NPTF2F_ReadFromOutputTXTFiles(str2double(sub_list{i}), num2str(seqi), date_list{i}, num2str(runi), false);
        end
    end
end
disp('-----------------ProcessData Complete-----------------');

%% SECTION 2 - VIDEOS AND LINK SCORES
% Watch the Videos and Fill in the Number of Total Links, Correct Links, and Incorrect Links
%-------------------------------------------------------------------------
% This section runs through each test and displays them in the figure that
% pops up. Press 'Enter' for a perfect trial, or input the Total Links,
% Correct Links, and Errors for each test. 
% Notes: # of errors reported should be individual groupings. I.e. 1-3-2-4
% is only 1 error (considered a swap of 2 and 3). Also, the number of
% correct links should be reported UP UNTIL THE FIRST ERROR! 
% Ex: [...]-21-23-24-25 is 23 total links (missed 22), 20 correct links (before missing 22) and 1 error 

for b = 1:length(seqns)
    for c = 1:length(parts)
        for d = 1:length(runns)
            for e = 1:length(testns)
                subn = num2str(subns(1));
                seqn = seqns(b);
                part = parts(c);
                runn = runns(d);
                testn= testns(e);
                
                % analysis parameters
                seqID = [num2str(subn,'%03.f') '_Seq' num2str(seqn)]; 
                runID = [seqID '_Run' num2str(runn)];
                parent_dir = [pwd '/'];
                save_dir = [parent_dir 'Participant Data/' seqID '/TMT_Output/'];
                stim_dir = [parent_dir 'gaudino-tmt-stim/'];
                

                % load data
                load([save_dir runID '.mat']);

                
                % get coordinate data
                if part == 'A'
                    data = trial_data.a{testn};
                    stim = stim_order.a{testn};
                    times = [1, 7999];
                elseif part == 'B'
                    data = trial_data.b{testn};
                    stim = stim_order.b{testn};
                    times = [1, 11997];
                end
                
                % Create and display Trial ID (Ex: B7)
                trialID= [subn ' - ' part num2str((4*(str2double(seqn)-1)+2*(runn-1)+testn))];
                disp(trialID)

                % Downsample to create 5ms grid and interp1 -> approximate realtime
                cursortimei = double(data(1,3) : 5 : data(end,3))';
                cursorposi = interp1(data(:,3), data(:,4:5), cursortimei);
                
                    % for loop is not needed here, but the times variable can be modified to
                    % play certain sections of each test, if desired.
                for time = 1:length(times)-1                    
                    %To play a portion of the video:
                    startTime = times(time);
                    endTime   = times(time + 1);


                    % animate response using ink
                    imagesc(imread(strcat(stim_dir, stim)));
                    hold on
                    ax = gca;
                    an = animatedline(ax,'Color','b','LineWidth',3);
                    set(gca, 'XAxisLocation', 'top', 'box', 'off', 'XTick', [])
                    set(gca, 'YAxisLocation', 'left', 'box', 'off', 'YTick', [])
                    set(gca, 'position', [0 0 1 1], 'units', 'normalized')                        

                    
                    %this goes from 1 to 8000 (A) or 12000 (B), which means that it is sampled every 200ms
                    %To get each interval of 5s, go to the next thousand
                    for idx = startTime:endTime-1 %1:1000 for 0s to 5s, 1:2000 for 0s to 10s, 1000:2000 for 5s to 10s
                        t = cursortimei(idx,1);
                        x = cursorposi(idx,1);
                        y = cursorposi(idx,2);
                        
                        %The tablet data stores the starting location
                        %at the centre of the screen (512, 384). Skip
                        %plotting the points at the beginning of the
                        %test if they have not yet touched the tablet
                        if not(and(x == 512, y == 384))
                            addpoints(an,x,y);
                        end 

                        %Add additional labels, if desired: 
                        %title(ax,sprintf('Time Elapsed:%6.1fs',t/1000));
                        %xlabel(ax,stim);

                        %Change the videoDisplayRate (Section 0) to display at different
                        %intervals. 
                        if mod(idx, videoDisplayRate) == 1 
                            drawnow
                        end
                            
                    end  
                    
                    %Prompt the TMT Score for each test
                    disp("Input the number of Total, Correct, and Incorrect Links below");
                    currentInput = input("Press Enter for a perfect trial or enter the Total number of Links: ");
                    
                    % Pressing Enter will input a perfect trial
                    if isempty(currentInput) 
                        if part == 'A'
                            links.Total.a{testn} = 24;
                            links.Correct.a{testn} = 24;
                            links.Incorrect.a{testn} = 0;
                        else %part == 'B'
                            links.Total.b{testn} = 24;
                            links.Correct.b{testn} = 24;
                            links.Incorrect.b{testn} = 0; 
                        end                      
                    else
                        if part == 'A'
                            links.Total.a{testn} = currentInput;
                            links.Correct.a{testn} = input("CORRECT (Up until 1st error): ");
                            links.Incorrect.a{testn} = input("INCORRECT (Number of major errors): ");
                        else
                            links.Total.b{testn} = currentInput;
                            links.Correct.b{testn} = input("CORRECT (Up until 1st error): ");
                            links.Incorrect.b{testn} = input("INCORRECT (Number of major errors): ");
                        end
                    end
                    
                    save([save_dir runID '.mat'], "links", "-append");
                    clc
                    
                end
            end
        end
    end
end
disp('-----------------LinkScores Complete-----------------');

%% Section 3 - RemoveErrors
% This section calls the NPTF2F_RemoveErrors() function for each
% TMT trial with errors recorded from Section 2. Select the starting and
% ending positions of the erroneous section to be removed. Right click the
% ending position if it is meant to remove all data for the remainder of
% the test (i.e., the error is continued throughout the entire test and
% should be removed up until the end).

for b = 1:length(seqns)
    for c = 1:length(runns)
         numPerfectTests = 0;
         for d = 1:length(parts)
            for e = 1:length(testns)
                subn = num2str(subns(1));
                seqn = seqns(b);
                runn = runns(c);
                part = parts(d);
                testn= testns(e);
                
                % analysis parameters
                seqID = [subn '_Seq' seqn];
                runID = [seqID '_Run' num2str(runn)];
                parent_dir = [pwd '/']; 
                save_dir = [parent_dir 'Participant Data/' seqID '/TMT_Output/'];
                stim_dir = [parent_dir 'gaudino-tmt-stim/'];
                

                % load data
                load([save_dir runID '.mat']);
                
                % Checks if errors were recorded
                if part == 'A' && links.Incorrect.a{testn} > 0
                    NPTF2F_RemoveErrors(subn, strcat(part, num2str((4*(str2double(seqn)-1)+2*(runn-1)+testn))));                  
                elseif part == 'B' && links.Incorrect.b{testn} >0
                    NPTF2F_RemoveErrors(subn, strcat(part, num2str((4*(str2double(seqn)-1)+2*(runn-1)+testn))));
                else 
                    numPerfectTests = numPerfectTests+1;
                end

                % Records number of perfect tests
                if numPerfectTests == 4
                    disp(['** No Errors in Sequence ' seqn ' | Run ' num2str(runn) ' **']);
                    clear numPerfectTests subn seqn part runn testn;
                    %Save everything except the variables listed below
                    save([save_dir runID ' - Errors Removed.mat'],'-regexp', '^(?!(b|c|d|e|date_list|parent_dir|parts|runID|runi|runns|save_dir|seqID|seqOrder|seqi|seqns|stim_dir|sub_list|testns|videoDisplayRate)$).');
                end

            end
        end
    end
end
disp('-----------------RemoveErrors Complete-----------------');

%% Section 4 - NPTF2F_SpeedData - Calculate Statistics from Trial Data
% This section calls the NPTF2F_SpeedData() function for each RUN

for seqi = 1:2
    for runi = 1:2
        %FORMAT: process_data(subn, seqn, date, runn, errorData) 
        NPTF2F_SpeedData(str2double (sub_list{1}), num2str(seqi), num2str(runi), false); %Full Trial Analysis
        NPTF2F_SpeedData(str2double (sub_list{1}), num2str(seqi), num2str(runi), true);  %Errors Removed Analysis

    end
end
disp('-----------------SpeedData Complete-----------------')

%% Section 5 - NPTF2F_SignalsData - Process Signals
% This section takes a long time to run. ~25 mins/Run * 8 Runs
% This section calls the NPTF2F_SignalsData() function for each Run 
% Make sure to check that the number of signals files failed = 0 every time

signalsFailed = 0;
for seqi = 1:2
    for runi = 1:2
        try
            NPTF2F_SignalsData(num2str(sub_list{1}), date_list{1}, num2str(seqi), num2str(runi), 'SRI', false);
        catch
            disp(['FULL TRIAL Signals failed for Seq: ' num2str(seqi) ' and Run ' num2str(runi)]);
            signalsFailed = signalsFailed + 1;
        end

        try 
            NPTF2F_SignalsData(num2str(sub_list{1}), date_list{1}, num2str(seqi), num2str(runi), 'SRI', true);
        catch
            disp(['ERRORS REMOVED Signals failed for Seq: ' num2str(seqi) ' and Run ' num2str(runi)]);
            signalsFailed = signalsFailed + 1;
        end
    end
end

disp(['Number of Signals Files failed: ' num2str(signalsFailed)]);
disp('-----------------SignalsData Complete-----------------')

%% Section 6 - Gather Data
% This section will gather the data into a usable format including: 
% - variablesFullTrial (all data from tests)
% - variablesErrorsRemoved (data calculations without erorrs)
% -  AnalyzedVariables.mat file (save easy-access variables)

%Setup the order of how the data will be displayed
testOrder = ["A1" "A2" "A3" "A4" "B1" "B2" "B3" "B4" "A5" "A6" "A7" "A8" "B5" "B6" "B7" "B8"];
testIndices = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];
testIndex = containers.Map(testOrder, testIndices);

%Header for the Variables Tables
variablesFullTrial (1,:)     = {"ParticipantID" "Sequence" "Configuration" "TMT Test" "Start Time" "Duration" "Total Links" "Correct Links" "Errors" "Mean Link Speed (px/ms)" "SPL"  "Mean Linking Duration" "Mean Non-Linking Duration" "Distance" "Force"};
variablesErrorsRemoved (1,:) = {"ParticipantID" "Sequence" "Configuration" "TMT Test" "Start Time" "Duration" "Total Links" "Correct Links" "Errors" "Mean Link Speed (px/ms)" "SPCL" "Mean Linking Duration" "Mean Non-Linking Duration" "Distance" "Force"};


% Obtain relevant force and distance data from NPTF2F_SignalsData.m output
forceFullTrial = NPTF2F_TrialForce(subns, false);
distanceFullTrial = NPTF2F_TrialDistance(subns, false);
forceErrorsRemoved = NPTF2F_TrialForce(subns, true);
distanceErrorsRemoved = NPTF2F_TrialDistance(subns, true);

%Perform data gathering for every trial
for seqn = 1:length(seqns) 
    seqID = [num2str(subns) '_Seq' seqns(seqn)]; %Participant ID
    seqType = seqOrder(3*str2double(seqns(seqn))-2:3*str2double(seqns(seqn)));
for runn = 1:length(runns)
for TMTType = 1:2
    switch TMTType
        case 1
            type = 'A';
        case 2
            type = 'B';
    end 
for testInRun =1:2    
    testNum = 4*(str2num(seqns(seqn))-1) + 2*(runns(runn)-1) + testInRun;
    TMTTest = [type num2str(testNum)];
    
    % Load and store TMT_Ouput data
    outputFolder = [pwd '/Participant Data/' seqID '/TMT_Output'];
    outputFile   = [outputFolder '/' seqID '_Run' num2str(runns(runn)) '.mat'];

    load(outputFile);
    
    %Save Full Trial Data
    if  type == 'A'
        variablesFullTrial(testIndex(TMTTest),:) = [ num2str(subns) seqType stim_order.a(testInRun) TMTTest start_end.a(testInRun, 1) completion_time.a(testInRun) links.Total.a(testInRun) links.Correct.a(testInRun) links.Incorrect.a(testInRun) mean_link_speed.a(testInRun) spl.a(testInRun) mean_link_duration.a(testInRun) mean_nonlink_duration.a(testInRun) distanceFullTrial(testNum + 1) forceFullTrial(testNum + 1)];

    elseif type == 'B'
        variablesFullTrial(testIndex(TMTTest),:) = [ num2str(subns) seqType stim_order.b(testInRun) TMTTest start_end.b(testInRun, 1) completion_time.b(testInRun) links.Total.b(testInRun) links.Correct.b(testInRun) links.Incorrect.b(testInRun) mean_link_speed.b(testInRun) spl.b(testInRun) mean_link_duration.b(testInRun) mean_nonlink_duration.b(testInRun) distanceFullTrial(testNum + 9) forceFullTrial(testNum + 9)];
    end
    
    
    %Now Perform the same thing with the Error Removed Data
    try
        outputFile   = [outputFolder '/' seqID '_Run' num2str(runns(runn)) ' - Errors Removed.mat'];
        load(outputFile);
    catch %Do nothing because the original file is already loaded
        disp('***** ERROR DATA FILE NOT FOUND *****')
    end
    
    %Save ErrorsRemoved Data
    if  type == 'A'
        variablesErrorsRemoved(testIndex(TMTTest),:) = [ num2str(subns) seqType stim_order.a(testInRun) TMTTest start_end.a(testInRun, 1) completion_time.a(testInRun) links.Total.a(testInRun) links.Correct.a(testInRun) links.Incorrect.a(testInRun) mean_link_speed.a(testInRun) spcl.a(testInRun) mean_link_duration.a(testInRun) mean_nonlink_duration.a(testInRun) distanceErrorsRemoved(testNum + 1) forceErrorsRemoved(testNum + 1)];
    
    elseif type == 'B'
        variablesErrorsRemoved(testIndex(TMTTest),:) = [ num2str(subns) seqType stim_order.b(testInRun) TMTTest start_end.b(testInRun, 1) completion_time.b(testInRun) links.Total.b(testInRun) links.Correct.b(testInRun) links.Incorrect.b(testInRun) mean_link_speed.b(testInRun) spcl.b(testInRun) mean_link_duration.b(testInRun) mean_nonlink_duration.b(testInRun) distanceErrorsRemoved(testNum + 9) forceErrorsRemoved(testNum + 9)];
    end
end
end 
end
end

savedFile = [pwd '/Participant Data/' num2str(subns) ' - AnalyzedVariables.mat'];
clear TMTTest TMTType subns seqns seqn runn runns pid testInRun type file folder deletePictures outputFile outputFolder parts seqType subn sub_list testIndices testIndex testns testNum testOrder videoDisplayRate;
save(savedFile);
disp("-----------------GatherData Complete-----------------")


%% FUNCTION 1 - Process data from tablet log and split into trials
function NPTF2F_ReadFromOutputTXTFiles(subn, seqn, date, runn, prep) 
disp(['Subject ' num2str(subn) ' | Sequence ' seqn ' | Run ' runn])  
  
% tmt design parameters
ctln = 0; % number of control trials per repetition
an = 1; % number of tmt a trials 
bn = 1; % number of tmt b trials
stim_count = ctln + an + bn; % total count of stimuli per repetition
tmt_rep = 2; % number of tmt repetition per part
% tmt_rep = 4; % EEG
duration.a = 40000; % ms
duration.b = 60000;
% analysis parameters
seqID = [num2str(subn,'%03.f') '_Seq' num2str(seqn)]; % NPTF2F 
% sid = ['s' subn]; % EEG
rid = [seqID '_Run' runn];
parent_dir = [pwd '/']; % path to /sgraham_data/NPTF2F
data_dir = [parent_dir 'Participant Data/' seqID '/'];
mkdir( [data_dir 'TMT_Output']);
save_dir = [data_dir 'TMT_Output/'];

%%% preprocessing %%%

disp('Reading tmt data......')    
% read tablet tmt log data
dataFileName = ['TMTRun' runn 'Log_S' num2str(subn) '_' date '_' seqn '.txt']; % NPTF2F data
% dataFileName = ['TMTRun1Log_S' subn '_' date '_' runn '.txt']; % EEG
disp([data_dir dataFileName])
fid = fopen([data_dir dataFileName],'rt');
fileData = textscan(fid, '%d %s %d %d %d %d %d %d %d', 'Headerlines',0);
fclose(fid);

% select and log raw tablet data
stim.index = fileData{:,1};
stim.name = fileData{:,2};
time.runOnset = fileData{:,3}; % total run time of the session
time.stimOnset = fileData{:,4}; % time relative to stim onset
coordinates.x = fileData{:,5};
coordinates.y = fileData{:,6};
forceFlag1= fileData{:,7}; 
force = fileData{:,8}; 
forceFlag2 = fileData{:,9}; 

% This will overwrite the current data, so do not rerun on an existing
% participant
save([save_dir rid '.mat'], 'stim','time','coordinates')

% split data based on each tmt trial
coorData = double([stim.index time.runOnset time.stimOnset coordinates.x coordinates.y forceFlag1 force forceFlag2]);
for trialn = 1:tmt_rep
    a = coorData(stim.index==stim_count*trialn-1,:);
    stim_order.a{trialn} = unique(string(stim.name(find(stim.index==stim_count*trialn-1),:))); % log stim order
    start_end.a(trialn,:) = [a(1,2)-a(1,3),duration.a-a(end,3)+a(end,2)]; % log start and end time for sync
    trial_data.a{trialn} = a; % log trial data

    b = coorData(stim.index==stim_count*trialn,:);
    stim_order.b{trialn} = unique(string(stim.name(find(stim.index==stim_count*trialn),:)));
    start_end.b(trialn,:) = [b(1,2)-b(1,3),duration.b-b(end,3)+b(end,2)];
    trial_data.b{trialn} = b;
end
save([save_dir rid '.mat'], 'stim_order','start_end','trial_data', '-append');
end

%% FUNCTION 2 - AverageForce
function averageForce = NPTF2F_TrialForce(subn, ErrorsRemoved)

    for seqn = 1:2
        for runn = 1:2
        
        % Load force data from TMT_Signals data
        seqID = [num2str(subn) '_Seq' num2str(seqn)]; %Participant ID
        signalsFolder = [pwd '/Participant Data/' seqID '/TMT_Output'];
                    
        
        %Load FullTrial signals file
        signalsFile = [signalsFolder '/' seqID '_Run' num2str(runn) '.mat'];
        load(signalsFile);  

        
        % Load ErrorsRemoved Data, if prompted. If it does not exist, it
        % will still load Full Trial file
        if ErrorsRemoved
            try
                signalsFile = [signalsFolder '/' seqID '_Run' num2str(runn) ' - Errors Removed.mat'];
                load(signalsFile);  
            catch
                disp(['FORCE    ' num2str(seqn) '|' num2str(runn) '             - No Errors were found for SEQUENCE ' num2str(seqn) ' and RUN ' num2str(runn)]);
            
            end
        end

        %Obtain the force readings and remove the base force that is always
        %present
        forces = zeros(10004,1);
        baseForce = mode(col_f2);
        forces(tmt_df.f1 == 1) = tmt_df.f2(tmt_df.f1 == 1) - baseForce;
        forces( forces < 2 ) = 0; %Removes any small fractional or negative values

        %Calculate the starting and ending times of each test based on the
        %NaN indices
        NaNIndices = find(isnan(tmt_df.f1));
        startIndices = [NaNIndices(1) NaNIndices(3) NaNIndices(5) NaNIndices(7)] + 1;
        endIndices   = [NaNIndices(2) NaNIndices(4) NaNIndices(6) NaNIndices(8)] - 1;

        % Save the values in the force array. 
        % Left side: save into positions matching the format above. 
        % Right side: specific start and end times
        averageForce(4*(seqn-1) + 2*runn    )  = nanmean(forces(startIndices(1) : startIndices(1) + ceil((cell2mat(completion_time.a(1))/20))-2));
        averageForce(4*(seqn-1) + 2*runn + 1)  = nanmean(forces(startIndices(2) : startIndices(2) + ceil((cell2mat(completion_time.a(2))/20))-2));
        averageForce(4*(seqn-1) + 2*runn + 8)  = nanmean(forces(startIndices(3) : startIndices(3) + ceil((cell2mat(completion_time.b(1))/20))-2));
        averageForce(4*(seqn-1) + 2*runn + 9)  = nanmean(forces(startIndices(4) : startIndices(4) + ceil((cell2mat(completion_time.b(2))/20))-2));
        end
    end
end

%% FUNCTION 3 - Distance 
function distance = NPTF2F_TrialDistance(subn, ErrorsRemoved)

    for seqn = 1:2
        for runn = 1:2
        
        % Load force data from TMT_Signals data
        seqID = [num2str(subn) '_Seq' num2str(seqn)]; %Participant ID
        signalsFolder = [pwd '/Participant Data/' seqID '/TMT_Output'];     
        signalsFile = [signalsFolder '/' seqID '_Run' num2str(runn) '.mat'];
        load(signalsFile); 

%       Try to load Errors Removed Data, if prompted when calling the
%       function
        if ErrorsRemoved
            try
                signalsFile = [signalsFolder '/' seqID '_Run' num2str(runn) ' - Errors Removed.mat'];
                load(signalsFile);  
               
            catch        
                disp(['DISTANCE ' num2str(seqn) '|' num2str(runn) '             - No Errors were found for SEQUENCE ' num2str(seqn) ' and RUN ' num2str(runn)]);
            end
        end

%       Find NaN start and end indices from the coordinates 
        NaNIndices = find(isnan(tmt_df.x));
        startIndices = [NaNIndices(1) NaNIndices(3) NaNIndices(5) NaNIndices(7)] + 1;
        endIndices   = [NaNIndices(2) NaNIndices(4) NaNIndices(6) NaNIndices(8)] - 1;

%       Save the values in the distance array.
%       Left side: save into positions matching the format above. 
%       Right side: test start and ends from NaN rows
        distance(4*(seqn-1) + 2*runn    )  = nansum(tmt_df.rel_dist(startIndices(1):endIndices(1)));
        distance(4*(seqn-1) + 2*runn + 1)  = nansum(tmt_df.rel_dist(startIndices(2):endIndices(2)));
        distance(4*(seqn-1) + 2*runn + 8)  = nansum(tmt_df.rel_dist(startIndices(3):endIndices(3)));
        distance(4*(seqn-1) + 2*runn + 9)  = nansum(tmt_df.rel_dist(startIndices(4):endIndices(4)));

        end
    end
%     disp('Distance acquired')
end
