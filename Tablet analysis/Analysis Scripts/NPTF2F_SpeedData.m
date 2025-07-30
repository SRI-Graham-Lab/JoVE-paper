% Process data from tablet log and split into trials
% example: process_data(64, '1', '2021-10-28'), if prep is true, run
% preprocessing only, else run the rest only


% Before running:
% Ensure that NPTF2F_TMTScores.csv is updated. 
% Format: (Unknown as of right now - check convert_score.m)
% Note: This script takes ~ 4 min/participant (~ 1 min/run)


function NPTF2F_SpeedData(subn, seqn, runn, errorData) 
if errorData
    disp(['Speed - Subject ' num2str(subn) ' | Sequence ' seqn ' | Run ' runn ' - Errors Removed'])  
else
    disp(['Speed - Subject ' num2str(subn) ' | Sequence ' seqn ' | Run ' runn ' - Full Trial'])  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialize variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
runID = [seqID '_Run' runn];
parent_dir = [pwd '/']; % path to /sgrahamlab/labdata/NPTF2F
data_dir = [parent_dir 'Participant Data/' seqID '/'];
mkdir( [data_dir 'TMT_Output']);
save_dir = [data_dir 'TMT_Output/'];
messagesFileName = [save_dir runID ' - Messages'];

load([save_dir runID '.mat']);

if errorData
    try 
%         disp(['Trying to load data with errors removed from: ' strcat(save_dir, runID, ' - Errors Removed.mat')] )
        load(strcat(save_dir, runID, ' - Errors Removed.mat'));
%         disp('Data with errors removed loaded successfully!')
        fileName = strcat(save_dir, runID, ' - Errors Removed.mat');    
%         save([save_dir runID ' - Errors Removed.mat'], 'stim_order','start_end','trial_data', '-append')

    catch
        disp('************* ERROR - DATA WITH ERRORS REMOVED NOT FOUND - ERROR *************')
        errorData = false;
    end
else
%     save([save_dir runID '.mat'], 'stim_order','start_end','trial_data', '-append')
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Speed %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     disp('Processing speed data......')
    % filter parameters
    Fs = 1000;   % Sampling frequency = 1000Hz to match downsampled EEG 
    w = 191; % sigma filter window width, default 151
    % loop to obtain filtered speed and other metrics    
    for parti = 1:2
        part_ls = ['a','b'];
        for trialn = 1:tmt_rep
            
%             %UNCOMMENT THIS SECTION AND ADD THE SPECIFIC TRIAL THAT BOARD
%             %WIDTH NEEDS TO BE INCREASED FOR. 
%             if parti == 2 && trialn == 1 && seqn == '2' && runn == '2'
%                 w = 235; %VARY THIS VALUE
%             else 
%                 w = 191; %DEFAULT - KEEP
%             end 

            %%%% calculate %%%%
            TDWithNans = trial_data.(part_ls(parti)){trialn};
            TDNoNans = rmmissing(TDWithNans);

%           Remove the last 2 lines of TDNoNans (data not accurate if
%           errors are removed)
            TDNoNans(end-1:end, :) = [];

            % oversample
            cursortimei = double(TDNoNans(1,3) : 1 : TDNoNans(end,3))';
            % interpolate using pchip
            cursorposi = interp1(TDNoNans(:,3), TDNoNans(:,4:5), cursortimei,'pchip');
            % first derivative as speed, second as acceleration
            dl = hypot(diff(cursorposi(:,1)),diff(cursorposi(:,2)));
            dldt = dl./diff(cursortimei);
            dldt2 = diff(dldt)./diff(cursortimei(2:end));
            % log data
            speed.(part_ls(parti)){trialn} = dldt;
            sampling_delay.(part_ls(parti)){trialn} = [cursortimei(1), cursortimei(end)];
            trial_time.(part_ls(parti)){trialn} = cursortimei;
            acceleration.(part_ls(parti)){trialn} = dldt2;

            %%%% filter %%%%
            % 1-D sigma filter
            w = 191;
            windowWidthTooSmall = true;
            while windowWidthTooSmall
                [avgx,~] = sigfilt1(dldt,2,w,floor(w/2),true);
                % adjust parameters to get the minimum window width to not get NaN
                if sum(isnan(avgx)) > 0
%                     disp(parti);
%                     disp(trialn);
%                     error('Sigma filter error: NaN detected, increase the window width!')
                    w = w + 2;
                else 
                    windowWidthTooSmall = false;
                end
            end

            if w > 191
                message = ['Window Width increased to ' num2str(w) ' for Sequence ' num2str(seqn) ' | Run ' num2str(runn) part_ls(parti) num2str(trialn)];
                disp(message);
                OutputMessages(messagesFileName, message)
            end

            % spikeRemoval
            [x,~] = spikeRemoval(avgx,'debug','');
            % low pass filter at 5 Hz
            lp = lowpass(x,5,Fs);
            lp(lp<0)=0; % speed cannot be negative. skip for acceleration
            % log filtered data
            speedf.(part_ls(parti)){trialn} = lp;

            %%%% completion time %%%%
            k = find(lp,1,'last'); % finds the last 1 index corresponding to nonzero elements in filtered speed
%             if k == length(lp) % if last data is nonzero, assume ceiling completion time
%                 completion_time.(part_ls(parti)){trialn} = duration.(part_ls(parti));
%             else % else 
%                 completion_time.(part_ls(parti)){trialn} = cursortimei(k+1);  % k+1 is the start time of the last zero speed segment      
%             end
            completion_time.(part_ls(parti)){trialn} = cursortimei(k+1);

        end
    end
    % save speed related metrics
    if errorData
        save([save_dir runID ' - Errors Removed.mat'], 'completion_time', 'speed','sampling_delay','trial_time','acceleration','speedf', '-append')
    else
        save([save_dir runID '.mat'], 'completion_time', 'speed','sampling_delay','trial_time','acceleration','speedf', '-append')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SPL %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     disp('Processing seconds per link data......')


    %Calculate SPL and store links data (total, correct, errors)
    for parti = 1:2 %sequence
        part_ls = ['a','b'];    
        for trialn = 1:tmt_rep
            testNumber = 4*(str2double(seqn)-1) + 2*(str2double(runn) -1) + trialn;
%             trialid = str2double([num2str(subn) seqn runn num2str(parti) num2str(testNumber)]);
            % split trial data
%             total_link.(part_ls(parti)){trialn} = links(row_idx,2);
%             correct_link.(part_ls(parti)){trialn} = links(row_idx,3);
%             error_link.(part_ls(parti)){trialn} = links(row_idx,4);
            spl.(part_ls(parti)){trialn}  = (completion_time.(part_ls(parti)){trialn} * .001) / links.Total.(part_ls(parti)){trialn};
            spcl.(part_ls(parti)){trialn} = (completion_time.(part_ls(parti)){trialn} * .001) / links.Correct.(part_ls(parti)){trialn};
            % tabulate trial data
%             links(row_idx,5) = completion_time.(part_ls(parti)){trialn};
%             links(row_idx,6) = spl.(part_ls(parti)){trialn};

        end
    end
    
    if errorData
        save([save_dir runID ' - Errors Removed.mat'], "links", 'spl', 'spcl', '-append')
    else
        save([save_dir runID '.mat'], "links",'spl','-append')
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Linking and non-linking periods %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     disp('Processing link timing data......')

    % initialize variables for an incrementing first threshold and a fixed
    % second threshold
    step.a = 0.002;
    step.b = 0.005;
    mean_speed_thre.a = 0.030;
    mean_speed_thre.b = 0.067;

    % incremental speed thresholding, speed unit: px/ms
    part_ls = ['a','b']; 
    for parti = 1:2   
        for trialn = 1:tmt_rep        
            %%% initialize variables %%%
            thre = 0; % 1st threshold
            % time indices for neuroimaging sync/analysis
            link_prd = {};
            nonlink_prd = {};
            % average speed per link
            mean_speed = {};
            % difference between true number of links and those identified by
            % threshold
            thre_diff = [];

            % input speed and time and link number 
            k = find(speedf.(part_ls(parti)){trialn},1,'last'); % drop post completion time
            spd = speedf.(part_ls(parti)){trialn}(1:k);
            t = trial_time.(part_ls(parti)){trialn}(1:k);        
            lk_num = links.Total.(part_ls(parti)){trialn};

            %%% apply 1st threshold %%%
            thre_spd = spd - thre;
            % threshold intercept timings
            link_idx_full = zeroX(t,thre_spd); % obtain time indices of speed zero crossings 
            link_idx_cut = link_idx_full(1:floor(length(link_idx_full)/2)*2); % retain max even number length
            link_start_end = [link_idx_cut(:,1:2:end)' link_idx_cut(:,2:2:end)']; % organize into link start and link end

            % calculate average speed within each link %
            lk_avg = [];
            for i = 1:size(link_start_end,1)
                lk_avg(i) = mean(spd(link_start_end(i,1):link_start_end(i,2)));
            end

            %%% apply link speed cutoff (2nd threshold) %%%
            link_spd = lk_avg(lk_avg>mean_speed_thre.(part_ls(parti)));

            %%% find the best match threshold %%%
            while thre < 1 && length(link_spd) ~= lk_num 

                % increase 1st threshold by pre-defined step size            
                thre = thre + step.(part_ls(parti));
                thre_idx = int32(thre/step.(part_ls(parti)));
                % recalculate speed and time
                thre_spd = spd - thre;

                % average speed per link
                link_idx_full = zeroX(t,thre_spd); % zero crossing time indices
                link_idx_cut = link_idx_full(1:floor(length(link_idx_full)/2)*2);
                link_start_end = [link_idx_cut(:,1:2:end)' link_idx_cut(:,2:2:end)']; % organize into link start and link end

                % calculate average speed for each link
                lk_avg = [];
                for i = 1:size(link_start_end,1)
                    lk_avg(i) = mean(spd(link_start_end(i,1):link_start_end(i,2)));
                end
                % apply 2nd threshold
                link_spd = lk_avg(lk_avg>mean_speed_thre.(part_ls(parti)));
                mean_speed{thre_idx} = link_spd;

                % reorganize start end time markers based on 2nd threshold
                link_prd{thre_idx} = link_start_end(find(lk_avg>mean_speed_thre.(part_ls(parti))),:);
                nonlink_prd{thre_idx} = [[0; link_prd{thre_idx}(1:end-1,2)] link_prd{thre_idx}(:,1)];

                % threshold error
                thre_diff(thre_idx) = length(link_spd) - lk_num;
            end

            %%% log data %%%
            % determine the threshold that produce closest link number match
            [threshold_err.(part_ls(parti)){trialn},I] = min(abs(thre_diff));
            threshold.(part_ls(parti)){trialn} = I * step.(part_ls(parti)); 

            % time periods (start, end) for link and non-link
            link_time_index.(part_ls(parti)){trialn} = link_prd{I};
            nonlink_time_index.(part_ls(parti)){trialn} = nonlink_prd{I};

            % per link metrics
            % average speed
            link_speed.(part_ls(parti)){trialn} = mean_speed{I};
            % linking durations
            link_duration.(part_ls(parti)){trialn} = link_prd{I}(:,2) - link_prd{I}(:,1);
            % non-linking durations
            nonlink_duration.(part_ls(parti)){trialn} = nonlink_prd{I}(:,2) - nonlink_prd{I}(:,1);

            % per trial metrics
            % average durations
            mean_link_speed.(part_ls(parti)){trialn} = mean(mean_speed{I});
            mean_link_duration.(part_ls(parti)){trialn} = mean(link_prd{I}(:,2) - link_prd{I}(:,1));
            mean_nonlink_duration.(part_ls(parti)){trialn} = mean(nonlink_prd{I}(:,2) - nonlink_prd{I}(:,1));                
        end
    end
    
    if errorData
        save([save_dir runID ' - Errors Removed.mat'], 'threshold_err','threshold','link_speed','link_duration','nonlink_duration','link_time_index','nonlink_time_index','-append')
        save([save_dir runID ' - Errors Removed.mat'], 'mean_link_speed','mean_link_duration','mean_nonlink_duration','-append')
    else
        save([save_dir runID '.mat'], 'threshold_err','threshold','link_speed','link_duration','nonlink_duration','link_time_index','nonlink_time_index','-append')
        save([save_dir runID '.mat'], 'mean_link_speed','mean_link_duration','mean_nonlink_duration','-append')
    end



end


function OutputMessages(fileName, message)
    try
        load(fileName)
        messageSaved = [messageSaved message];
    catch
        messageSaved = [message];
    end
    try
        save(fileName, "messageSaved", '-append');
    catch
        save(fileName, "messageSaved");
    end
    disp(['Writing Increased Board Size Message to' fileName]);
end 