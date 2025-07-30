% Process Key Signals from TMT_data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used by Sean Rose (Winter 2024), taken from previous work in tmt_analysis 
% with the following note:
    % fohara, Winter 2023
    % OVERVIEW: Process a .csv dataframe from raw tmt_data files for each
    % subject. Speed and other preliminary signals are calculated using
    % Zhongmin's previous scripts. New features such as proximity to hitboxes
    % are added to processing.
%Note: This script takes ~22  minutes to run for only one .txt file!
function NPTF2F_SignalsData(subn, date, seqn, runn, site, errors)
% Run from NPTF2F folder
% cd('/sgraham_data/NPTF2F')
% clear
% clc
if errors
    disp(['Signals - Subject ' num2str(subn) ' | Sequence ' seqn ' | Run ' runn ' - Errors Removed'])  
else
    disp(['Signals - Subject ' num2str(subn) ' | Sequence ' seqn ' | Run ' runn ' - Full Trial'])  
end
%% Main

% subj_list = ['NPTF2F_899V1'];
seqID = strcat(subn, '_Seq', seqn);
filename = ['TMTRun', runn, 'Log_S' subn, '_', date, '_', seqn, '.txt'];
% filename_list = ['TMTRun1Log_S899_2023-08-03_1.txt'];
% 
% site = 'SRI'; % 'SRI' or 'SMH', difference in number of columns

    % Trial info
    duration.a = 40000; % ms
    duration.b = 60000;
    ctln = 8; % number of control trials per repetition
    an = 1; % number of tmt a trials 
    bn = 1; % number of tmt b trials
    stim_count = ctln + an + bn; % total count of stimuli per repetition
    tmt_rep = 2; % number of tmt repetition per part

    % Setup info
    Fs = 1000;   % Sampling frequency = 1000Hz to match downsampled EEG 
    w = 191; % sigma filter window width, default 151

    % Processing variables
    parent_dir = [pwd '/']; % path to /tablet-tmt-analysis
    data_dir = [parent_dir 'Participant Data/' seqID '/'];
    save_dir = [data_dir 'TMT_Output/'];
%     disp([data_dir filename])
        

    %% Open
    % Load file
    if strcmp(site, 'SRI')
        % Load file more less columns
        fid = fopen([data_dir filename],'rt');
%         file_rawdata = fopen(fid, '%d %s %d %d %d %d %d %d %d', 'Headerlines',0);
        file_rawdata = textscan(fid, '%d %s %d %d %d %d %d %d %d', 'Headerlines',0);
        fclose(fid);
        
        stim.index = file_rawdata{:,1};
        stim.name = file_rawdata{:,2};
        time.runOnset = file_rawdata{:,3}; % total run time of the session
        time.stimOnset = file_rawdata{:,4}; % time relative to stim onset
        coordinates.x = file_rawdata{:,5};
        coordinates.y = file_rawdata{:,6};
        force.one = file_rawdata{:,7};
        force.two = file_rawdata{:,8};
        force.three = file_rawdata{:,9};
    elseif strcmp(site, 'SMH')
        % Load file with less columns
        fid = fopen([data_dir filename],'rt');
        file_rawdata = textscan(fid, '%s %d %d %d %d', 'Headerlines',0);
        fclose(fid);
        
        % Swich columns from the initial format
        stim.name = file_rawdata{:,1};
        time.runOnset = file_rawdata{:,2}; % total run time of the session
        time.stimOnset = file_rawdata{:,3}; % time relative to stim onset
        coordinates.x = file_rawdata{:,4};
        coordinates.y = file_rawdata{:,5};
        force.one = zeros(length(file_rawdata{:,1}), 1);
        force.two = zeros(length(file_rawdata{:,1}), 1);
        force.three = zeros(length(file_rawdata{:,1}), 1);
        
        % Calculate stim index from filename
        stim.index = zeros(length(file_rawdata{:,1}), 1);
        stim.index(strcmp('Gaudino-TMTB_n-small.png', file_rawdata{:,1}))   = 9;
        stim.index(strcmp('Gaudino-TMTA_nl-small.png', file_rawdata{:,1}))  = 10;
        stim.index(strcmp('Gaudino-TMTA-small.png', file_rawdata{:,1}))     = 19;
        stim.index(strcmp('Gaudino-TMTB-small.png', file_rawdata{:,1}))     = 20;
        stim.index(strcmp('Gaudino-TMTA_rot-small.png', file_rawdata{:,1})) = 9;
        stim.index(strcmp('Gaudino-TMTB_rot-small.png', file_rawdata{:,1})) = 10;
        stim.index(strcmp('Gaudino-TMTB_n_rot-small.png', file_rawdata{:,1})) = 19;
        stim.index(strcmp('Gaudino-TMTA_nl_rot-small.png', file_rawdata{:,1})) = 20;
    end

    % Reformat
    coorData = double([stim.index time.runOnset time.stimOnset coordinates.x coordinates.y, force.one force.two force.three]);
    for trialn = 1:tmt_rep
        a = coorData(stim.index== (2*trialn - 1), :); %stim_count*trialn-1,:);
        stim_order.a{trialn} = unique(string(stim.name(find(stim.index==stim_count*trialn-1),:))); % log stim order
        start_end.a(trialn,:) = [a(1,2)-a(1,3),duration.a-a(end,3)+a(end,2)]; % log start and end time for sync
        trial_data.a{trialn} = a; % log trial data    

        b = coorData(stim.index== (2*trialn), :); %stim_count*trialn,:);
        stim_order.b{trialn} = unique(string(stim.name(find(stim.index==stim_count*trialn),:)));
        start_end.b(trialn,:) = [b(1,2)-b(1,3),duration.b-b(end,3)+b(end,2)];
        trial_data.b{trialn} = b;
    end

    if errors
        try 
            load(strcat(data_dir, "TMT_Output/", seqID, '_Run', runn, " - Errors Removed.mat"));        
%             disp(strcat(data_dir, "TMT_Output/", seqID, '_', runn, " - Errors Removed.mat"));
%             disp('Error Data Used');
    
        catch 
            disp('**** ERROR DATA NOT USED ****');
        end
    else
        load(strcat(data_dir, "TMT_Output/", seqID, '_Run', runn, '.mat'));
    end

%     %% Speed - Performed in SpeedData
%     % Referenced from Zlin code
% %     disp('Calculating Cursor Speed...')
%     for parti = 1:2
%         part_ls = ['a','b'];
%         for trialn = 1:tmt_rep
% 
% %             %UNCOMMENT THIS SECTION AND ADD THE SPECIFIC TRIAL THAT BOARD
% %             %WIDTH NEEDS TO BE INCREASED FOR. 
% %             if parti == 2 && trialn == 1 
% %                 w = 235; %VARY THIS VALUE
% %                 disp(['Board Width was changed for parti = ' num2str(parti) ' and trialn = ' num2str(trialn)]);
% %                 disp('!!!!!!!!!!Make sure to chance the settings back! Lines 120 - 126!!!!!!!!!!!!!!');
% %             else 
% %                 w = 191; %DEFAULT - KEEP
% %             end 
% 
% 
% 
% 
%             %%%% calculate %%%%
%             td = trial_data.(part_ls(parti)){trialn};
%             % oversample
%             cursortimei = double(td(1,3) : 1 : td(end,3))';
%             % interpolate using pchip
%             cursorposi = interp1(td(:,3), td(:,4:5), cursortimei,'pchip');
%             % first derivative as speed, second as acceleration
%             dl = hypot(diff(cursorposi(:,1)),diff(cursorposi(:,2)));
%             dldt = dl./diff(cursortimei);
%             dldt2 = diff(dldt)./diff(cursortimei(2:end));
%             % log data
%             speed.(part_ls(parti)){trialn} = dldt;
%             sampling_delay.(part_ls(parti)){trialn} = [cursortimei(1), cursortimei(end)];
%             trial_time.(part_ls(parti)){trialn} = cursortimei;
%             acceleration.(part_ls(parti)){trialn} = dldt2;
% 
%             %%%% filter %%%%
%             % 1-D sigma filter
%             w = 191;
%             windowWidthTooSmall = true;
%             while windowWidthTooSmall
%                 [avgx,~] = sigfilt1(dldt,2,w,floor(w/2),true);
%                 % adjust parameters to get the minimum window width to not get NaN
%                 if sum(isnan(avgx)) > 0
% %                     error('Sigma filter error: NaN detected, increase the window width!')
%                     w = w + 2;
%                 else 
%                     windowWidthTooSmall = false;
%                 end
%             end
% % 
% %             if w > 191
% %                 message = ['Window Width increased to ' num2str(w) ' for Sequence ' num2str(seqn) ' | Run ' num2str(runn) part_ls(parti) num2str(trialn)];
% %                 disp(message);%                     disp(parti);
% %                     disp(trialn);
% %                 OutputMessages(messagesFileName, message)
% %             end
%             % spikeRemoval
%             [x,~] = spikeRemoval(avgx,'debug','');
%             % low pass filter at 5 Hz
%             lp = lowpass(x,5,Fs);
%             lp(lp<0)=0; % speed cannot be negative. skip for acceleration
%             % log filtered data
%             speedf.(part_ls(parti)){trialn} = lp;
% 
%             %%%% completion time %%%%
%             k = find(lp,1,'last'); % finds the last 1 index corresponding to nonzero elements in filtered speed
%             if k == length(lp) % if last data is nonzero, assume ceiling completion time
%                 completion_time.(part_ls(parti)){trialn} = duration.(part_ls(parti));
%             else % else 
%                 completion_time.(part_ls(parti)){trialn} = cursortimei(k+1);  % k+1 is the start time of the last zero speed segment      
%             end
%         end
%     end


    %% Resample Data into Output Table

    % The raw data comes at 1000 Hz, but is logged per change in coordinate. I
    % will resample to 50 Hz to reduce the number of proximity calculations.
    % Therefore, will convert to increment time by 20 ms instead of 1 ms. Will
    % also interpolate trials to last the full duration.
    
%     disp('Resampling data to 50 Hz...')

    % Calculate new time columns
    sample_period = 20; %ms; change here for sensitivity analysis
    time_a = 0:sample_period:duration.a;
    time_b = 0:sample_period:duration.b;

    % Interpolate raw signals to sampling frequency
    col_time = [time_a'; time_a'; time_b'; time_b'];

    col_x   = [ interp1(trial_data.a{1,1}(:,3), trial_data.a{1,1}(:,4), time_a)'; 
                interp1(trial_data.a{1,2}(:,3), trial_data.a{1,2}(:,4), time_a)'; 
                interp1(trial_data.b{1,1}(:,3), trial_data.b{1,1}(:,4), time_b)'; 
                interp1(trial_data.b{1,2}(:,3), trial_data.b{1,2}(:,4), time_b)'];

    col_y   = [ interp1(trial_data.a{1,1}(:,3), trial_data.a{1,1}(:,5), time_a)'; 
                interp1(trial_data.a{1,2}(:,3), trial_data.a{1,2}(:,5), time_a)'; 
                interp1(trial_data.b{1,1}(:,3), trial_data.b{1,1}(:,5), time_b)'; 
                interp1(trial_data.b{1,2}(:,3), trial_data.b{1,2}(:,5), time_b)'];

    col_f1  = [ interp1(trial_data.a{1,1}(:,3), trial_data.a{1,1}(:,6), time_a)'; 
                interp1(trial_data.a{1,2}(:,3), trial_data.a{1,2}(:,6), time_a)'; 
                interp1(trial_data.b{1,1}(:,3), trial_data.b{1,1}(:,6), time_b)'; 
                interp1(trial_data.b{1,2}(:,3), trial_data.b{1,2}(:,6), time_b)'];

    col_f2  = [ interp1(trial_data.a{1,1}(:,3), trial_data.a{1,1}(:,7), time_a)'; 
                interp1(trial_data.a{1,2}(:,3), trial_data.a{1,2}(:,7), time_a)'; 
                interp1(trial_data.b{1,1}(:,3), trial_data.b{1,1}(:,7), time_b)'; 
                interp1(trial_data.b{1,2}(:,3), trial_data.b{1,2}(:,7), time_b)'];

    col_f3  = [ interp1(trial_data.a{1,1}(:,3), trial_data.a{1,1}(:,8), time_a)'; 
                interp1(trial_data.a{1,2}(:,3), trial_data.a{1,2}(:,8), time_a)'; 
                interp1(trial_data.b{1,1}(:,3), trial_data.b{1,1}(:,8), time_b)'; 
                interp1(trial_data.b{1,2}(:,3), trial_data.b{1,2}(:,8), time_b)'];

    % Interpolate speed and acceleration signals to sampling frequency
    col_speed  = [ interp1(trial_time.a{1,1}(1:length(speed.a{1,1})), speed.a{1,1}, time_a)'; 
                   interp1(trial_time.a{1,2}(1:length(speed.a{1,2})), speed.a{1,2}, time_a)'; 
                   interp1(trial_time.b{1,1}(1:length(speed.b{1,1})), speed.b{1,1}, time_b)'; 
                   interp1(trial_time.b{1,2}(1:length(speed.b{1,2})), speed.b{1,2}, time_b)'];

    col_speedf = [ interp1(trial_time.a{1,1}(1:length(speedf.a{1,1})), speedf.a{1,1}, time_a)'; 
                   interp1(trial_time.a{1,2}(1:length(speedf.a{1,2})), speedf.a{1,2}, time_a)'; 
                   interp1(trial_time.b{1,1}(1:length(speedf.b{1,1})), speedf.b{1,1}, time_b)'; 
                   interp1(trial_time.b{1,2}(1:length(speedf.b{1,2})), speedf.b{1,2}, time_b)'];

    col_accel  = [ interp1(trial_time.a{1,1}(1:length(acceleration.a{1,1})), acceleration.a{1,1}, time_a)'; 
                   interp1(trial_time.a{1,2}(1:length(acceleration.a{1,2})), acceleration.a{1,2}, time_a)'; 
                   interp1(trial_time.b{1,1}(1:length(acceleration.b{1,1})), acceleration.b{1,1}, time_b)'; 
                   interp1(trial_time.b{1,2}(1:length(acceleration.b{1,2})), acceleration.b{1,2}, time_b)'];

    % Include trial labels
    col_sid   = repmat(seqID, (2*length(time_a)+2*length(time_b)), 1);
    col_subj  = repmat(str2double(subn), (2*length(time_a)+2*length(time_b)), 1);
    col_run   = repmat(str2double(runn), (2*length(time_a)+2*length(time_b)), 1);
    col_part  = [repmat('A', 2*length(time_a), 1); repmat('B', 2*length(time_b), 1)];
    col_trial = [repmat(1, length(time_a), 1); 
                 repmat(2, length(time_a), 1); 
                 repmat(1, length(time_b), 1); 
                 repmat(2, length(time_b), 1)];
    col_site   = repmat(site, (2*length(time_a)+2*length(time_b)), 1);
    
             

    %% Create preliminary table
    tmt_df = table;
    tmt_df.sid = col_sid;
    tmt_df.subj = col_subj;
    tmt_df.run = col_run;
    tmt_df.part = col_part;
    tmt_df.trial = col_trial;
    tmt_df.site = col_site;
    tmt_df.time = col_time;
    tmt_df.x = col_x;
    tmt_df.y = col_y;
    tmt_df.f1 = col_f1;
    tmt_df.f2 = col_f2;
    tmt_df.f3 = col_f3;
    tmt_df.speed = col_speed;
    tmt_df.speedf = col_speedf;
    tmt_df.accel = col_accel;

    %Set the position and force readings as NaNs at the start and end of
    %each trial
    tmt_df(    1, 8:14) = {NaN NaN NaN NaN NaN NaN NaN};
    tmt_df( 2001, 8:14) = {NaN NaN NaN NaN NaN NaN NaN};
    tmt_df( 2002, 8:14) = {NaN NaN NaN NaN NaN NaN NaN};
    tmt_df( 4002, 8:14) = {NaN NaN NaN NaN NaN NaN NaN};
    tmt_df( 4003, 8:14) = {NaN NaN NaN NaN NaN NaN NaN};
    tmt_df( 7003, 8:14) = {NaN NaN NaN NaN NaN NaN NaN};
    tmt_df( 7004, 8:14) = {NaN NaN NaN NaN NaN NaN NaN};
    tmt_df(10004, 8:14) = {NaN NaN NaN NaN NaN NaN NaN};


    %Remove any NaNs within each trial, but not between trials
    A1 = tmt_df(2:2000,:);
    A2 = tmt_df(2003:4001,:);
    B1 = tmt_df(4004:7002,:);
    B2 = tmt_df(7005:10003,:);
    testSeparations = [tmt_df(1,:); tmt_df(2001,:); tmt_df(2002,:); tmt_df(4002,:); tmt_df(4003,:); tmt_df(7003,:); tmt_df(7004,:); tmt_df(10004,:)];

    tmt_df = [testSeparations(1,:); rmmissing(A1); testSeparations(2:3,:); rmmissing(A2); testSeparations(4:5,:); rmmissing(B1); testSeparations(6:7,:); rmmissing(B2); testSeparations(8,:)];
    clear A1 A2 B1 B2 testSeparations;


    %% Board Region
    % The original board is split into a 3x3 grid labelled 1-9 moving left to
    % right and top to bottom. This is labelled board_region. Create by +1 when 
    % moving left to right and +3 moving down
    
%     disp('Calculating Board Region...')

    board_x_lim = 1024;
    board_y_lim = 768;
    tmt_df.board_region = ones(length(tmt_df.x), 1);
    tmt_df.board_region = tmt_df.board_region + (tmt_df.x>=board_x_lim/3) + (tmt_df.x>=2*board_x_lim/3);
    tmt_df.board_region = tmt_df.board_region + 3*(tmt_df.y>=board_y_lim/3) + (tmt_df.y>=2*board_y_lim/3);

    %% Calculate Relative Distance
    % In order for a desired cumulative distance feature to be easily
    % calculated, the distance travelled between time instances will be
    % simply calculated. The euclidean distance will be calculated and
    % saved into rel_dist. Because interpolation causes the first and last
    % coordinates to be NaN, there is no bleed over between trials.
    
%     disp('Calculating Relative Distance...')
    
    tmt_df.rel_dist = zeros(length(tmt_df.x), 1);
    for i=2:length(tmt_df.rel_dist)
        tmt_df.rel_dist(i) = sqrt((tmt_df.y(i)-tmt_df.y(i-1))^2 + (tmt_df.x(i)-tmt_df.x(i-1))^2);
    end
    
    %% Label Proximity
    % For each sample, the n-nearest borders pixels will be found and the
    % distance and angle to that border will be recorded. For completeness, n
    % will be the maximum number of targets. The labelled maps, their format,
    % and how they were made can be found in
    % tablet-tmt-analysis/gaudino-tmt-stim

    % Based on png names in TMT_data:
    %       run 1 part A trial 1 (9)  is TMTA-small 
    %       run 1 part A trial 2 (19) is TMTB_n-small 
    %       run 1 part B trial 1 (10) is TMTB-small 
    %       run 1 part B trial 2 (20) is TMTA_nl-small 
    %       run 2 part A trial 1 (9)  is TMTA_rot-small 
    %       run 2 part A trial 2 (19) is TMTB_n_rot-small 
    %       run 2 part B trial 1 (10) is TMTB_rot-small 
    %       run 2 part B trial 2 (20) is TMTA_nl_rot-small 

    % Instantiate all new columns as -1. There are a maximum of 25 targets.
    % These will be set as nearest 1 through 25 (n1, n2, n3, ... n25). For each
    % n, the label, distance, and angle will be saved.
%     disp('Instantiating Proximity Columns...')
    num_labels=25;
    for n_val=1:num_labels
        eval(strcat('tmt_df.n', string(n_val), '_label=nan(length(tmt_df.x), 1);'))
        eval(strcat('tmt_df.n', string(n_val), '_dist=nan(length(tmt_df.x), 1);'))
        eval(strcat('tmt_df.n', string(n_val), '_angle=nan(length(tmt_df.x), 1);'))
    end
    tmt_df.over_label = zeros(length(tmt_df.x), 1);
    tmt_df.center_prox = nan(length(tmt_df.x), 1);

    %% For each time instant: load the proper map, find closest point, save distance 
    % and angle, remove label, repeat until empty. Also consider special case
    % when distance is 0 in which a new map is loaded and the proximity to the
    % center is saved.

%     disp('Calculating Proximity Statistics...')
    numNaNs = 0;
    for i=1:length(tmt_df.x)

        % Load the proper trial board
        trialid = strcat(string(tmt_df.run(i)), tmt_df.part(i), string(tmt_df.trial(i)));
        switch trialid
            case "1A1"
                board_fn='gaudino-tmt-stim/Gaudino-TMTA/gaudino_tmta_ref.mat';
            case "1A2"
                board_fn='gaudino-tmt-stim/Gaudino-TMTB-n/gaudino_tmtb_n_ref.mat';
            case "1B1"
                board_fn='gaudino-tmt-stim/Gaudino-TMTB/gaudino_tmtb_ref.mat';
            case "1B2"
                board_fn='gaudino-tmt-stim/Gaudino-TMTA-nl/gaudino_tmta_nl_ref.mat';
            case "2A1"
                board_fn='gaudino-tmt-stim/Gaudino-TMTA-rot/gaudino_tmta_rot_ref.mat';
            case "2A2"
                board_fn='gaudino-tmt-stim/Gaudino-TMTB-n-rot/gaudino_tmtb_n_rot_ref.mat';
            case "2B1"
                board_fn='gaudino-tmt-stim/Gaudino-TMTB-rot/gaudino_tmtb_rot_ref.mat';
            case "2B2"
                board_fn='gaudino-tmt-stim/Gaudino-TMTA-nl-rot/gaudino_tmta_nl_rot_ref.mat'; 
        end
        % Get cursor coordinate and board
        load(board_fn)
        cursor_x = round(tmt_df.x(i))+1;
        cursor_y = round(tmt_df.y(i))+1;

        % Correct for nan coordinates
        if isnan(cursor_x) || isnan(cursor_y)
            cursor_x = 1;
            cursor_y = 1;
%             disp(strcat('Warning: df index ', string(i), ' has NaN coordinates. Setting to zero'))
            numNaNs = numNaNs + 1;
        end

        for n=1:num_labels
            % Get near-nonzero map
            [dist_map, idx] = bwdist(map_filled);
            idx_curr = idx(cursor_y, cursor_x);
            [row,col] = ind2sub (size(dist_map), idx_curr);

            % Save features
            eval(strcat("tmt_df.n", string(n), "_label(i) = map_filled(row, col);"));
            eval(strcat("tmt_df.n", string(n), "_dist(i) = dist_map(cursor_y, cursor_x);"));
            eval(strcat("tmt_df.n", string(n), "_angle(i) = atan2(row-cursor_y, col-cursor_x);")); 
            % NOTE: angle relative to CCS, not vision

            % Consider zero distance case
            if dist_map(cursor_y, cursor_x)==0
                tmt_df.over_label(i)=1;

                [center_dist_map, center_idx] = bwdist(map_center);
                tmt_df.center_prox(i) = center_dist_map(cursor_y, cursor_x);
            end

            % Remove the label from the map
            map_filled(map_filled(:,:)==map_filled(row, col)) = 0;
        end
    end
    
    if numNaNs ~= 8
        disp(['There is an incorrect number of NaN rows. There should be 8 (start and end of each of the 4 trials) but instead there are ' num2str(numNaNs)]);
    else
        clear numNaNs
    end
    %% Save the Dataframe

    mkdir(strcat('Participant Data/', seqID, '/TMT_Output/CSVs'))
    %Make output directory in tmt_signals
    runID = strcat(seqID, '_Run', runn);
    if errors
       OutputFileName = strcat('Participant Data/', seqID, '/TMT_Output/CSVs/', runID, ' - Errors Removed.csv');
    else
       OutputFileName = strcat('Participant Data/', seqID, '/TMT_Output/CSVs/', runID, '.csv');
    end
        
    writetable(tmt_df, OutputFileName)
%     disp(['Creating output file NPTF2F/' OutputFileName]) 
    



%% Section ##
% Saves the variables in the workspace as a separate file, which can be
% accessed later and used in the GatherDate.m file

if errors
%     SavedWorkspace = [save_dir 'NPTF2F_' seqID(end-4:end) ' - Run ' filename(7) ' - Errors Removed'];
%     SavedVariables = strcat('Participant Data/', seqID, '/TMT_Output/', seqID, '_', runn, ' - Errors Removed - Useful');
    SavedVariables = strcat('Participant Data/', seqID, '/TMT_Output/', runID, ' - Errors Removed');
else 
%     SavedWorkspace = [save_dir 'NPTF2F_' seqID(end-4:end) ' - Run ' filename(7)];
%     SavedVariables = strcat('Participant Data/', subj, '/TMT_Output/', subj, '_', runn, ' - Useful');
    SavedVariables = strcat('Participant Data/', seqID, '/TMT_Output/', runID);
end

% save(SavedWorkspace);

save(SavedVariables, "completion_time", "col_f2", "tmt_df", '-append');


% disp('Complete!')



% function OutputMessages(fileName, message)
%     try
%         load(fileName)
%         messageSaved = [messageSaved message];
%     catch
%         messageSaved = [message];
%     end
%     try
%         save(fileName, "messageSaved", '-append');
%     catch
%         save(fileName, "messageSaved");
%     end
%     disp(['Writing Increased Board Size Message to' fileName]);
% end 
