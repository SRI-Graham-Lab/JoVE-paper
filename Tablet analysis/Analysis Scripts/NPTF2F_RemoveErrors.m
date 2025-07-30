%% NPTF2F_RemoveErrors.m

%%Section 0: Script Explanation
% Author: Sean Rose - email s25rose@uwaterloo.ca with any questions
%This file requires the tmt_output files to have already been created,
%which can be done by running process_data.m first.

%% Section 1: Creating the images from the 'videos', for each interval
%Manually input the following prior to running this section of the script
%The ouput will be in the same folder as this file (sgraham_data/NPTF2F)
%Leave these files in this folder until AFTER you run Section 2
% clear;

function NPTF2F_RemoveErrors(subn, testWithErrors)
    
    % subn = 899;      %
    % testsWithErrors = ["B1" "B2"]; %
    
    %1000 for 5s intervals, 200 for 1s intervals, 20 for 0.1s intervals
    videoDisplayRate = 200; 
    deletePictures = false;
    global firstClick
    global cursorposi
    global cursortimei
    global save_dir
    global runID
    global testInRun
    global type
    global trial_data
    global data
    
    firstPoint = [0 0];
    secondPoint = [0 0];
    %-------------------------------------------------------------------------
    
%     for testn = 1:length(testWithErrors)
        % Define the Test to be Corrected
        
%         test = testWithErrors(testn).char();
 
        type = testWithErrors(1);
        testNumber = str2double(testWithErrors(2));
        testInRun = mod(testNumber+1,2)+1;
        if mod(testNumber, 4) == 1 || mod(testNumber, 4) == 2 
            runn = 1;
        else 
            runn = 2;
        end
        seqn = ceil(testNumber/4);
        disp(testWithErrors)
        
        
        % Load Run Data
        seqID = [num2str(subn) '_Seq' num2str(seqn)];
        runID = [seqID '_Run' num2str(runn)];
        parent_dir = [pwd '/']; % path to \tablet-tmt-analysis
        save_dir = [parent_dir 'Participant Data/' seqID '/TMT_Output/'];
        stim_dir = [parent_dir 'gaudino-tmt-stim/'];
        
    
        try load([save_dir runID ' - Errors Removed.mat']);
            ErrorRemovedData.a{1} = trial_data.a{1};
            ErrorRemovedData.a{2} = trial_data.a{2};
            ErrorRemovedData.b{1} = trial_data.b{1};
            ErrorRemovedData.b{2} = trial_data.b{2};
            clear trial_data;
    
        catch %Do nothing
            ErrorRemovedData.a{1} = [];
            ErrorRemovedData.a{2} = [];
            ErrorRemovedData.b{1} = [];
            ErrorRemovedData.b{2} = [];
        end
        
        load([save_dir runID '.mat']);
        
        if ~isempty(ErrorRemovedData.a{1})
            trial_data.a{1} = ErrorRemovedData.a{1};
        end
    
        if ~isempty(ErrorRemovedData.a{2})
            trial_data.a{2} = ErrorRemovedData.a{2};
        end
    
        if ~isempty(ErrorRemovedData.b{1})
            trial_data.b{1} = ErrorRemovedData.b{1};
        end
    
        if ~isempty(ErrorRemovedData.b{2})
            trial_data.b{2} = ErrorRemovedData.b{2};
        end
    
        % get stim and time data
        if type == 'A'
            data = trial_data.a{testInRun};
            stim = stim_order.a{testInRun};
            times = [1, 7999];
                           
        elseif type == 'B'
            data = trial_data.b{testInRun};
            stim = stim_order.b{testInRun};
            times = [1, 11999];    
        end
    
    
        
        % Downsample to create 5ms grid and interp1 -> approximate realtime
        cursortimei = double(data(1,3) : 5 : data(end,3));
        positions = interp1(data(:,3), data(:,4:5), cursortimei);
        cursorposi = round(positions,4);
    %     colors(1:length(cursorposi)) = 0.5;
    
    
        firstClick = true;
        
        imshow(imread(strcat(stim_dir, stim)));
        hold on
        plot = scatter(cursorposi(:,1), cursorposi(:,2), 20, 'blue', 'filled');  
    
        %Shows figure
        set(plot, 'ButtonDownFcn', @PointSelected);
        drawnow;
        input("Press Enter to go to the next test");
%     end
end
                    

function PointSelected(plotObject, EventData)
    global firstPoint
    global secondPoint
    global firstClick 
    
    endPoint = false;
    coordinates = EventData.IntersectionPoint(1:2);

    disp(['Point Selected: ' num2str(coordinates)])
    
    if firstClick % Store the value of the first point selected
        firstPoint = round(coordinates, 4);

    else % Remove the data from the first and second points selected
        try
            if get(gcf, 'SelectionType') == 'alt'
        %Check if the end was meant to be selected
%         if isEmpty(input('Press Enter if you selected the END of the trial, press any other key (then enter) if it was not the end '))
            disp('End Selected');
            endPoint = true; 
            end
        catch
            disp('End NOT Selected');
        end
        secondPoint = round(coordinates, 4);
        RemovePoints(firstPoint, secondPoint, endPoint);
        firstPoint = [0 0];
        secondPoint = [0 0];
    end
    firstClick = not(firstClick);
end

function RemovePoints(p1, p2, endPoint)

    global cursorposi
    global data
    global cursortimei
    global save_dir
    global runID
    global testInRun
    global type
    global trial_data

%     Find indexes in cursorposi where the x and y coordinates match the
%     points selected to be removed
    t1 = intersect(find(cursorposi(:, 1) == p1(1)), find(cursorposi(:, 2) == p1(2)))
    if endPoint % Point 2 should be the end
        t2 = length(cursorposi)
    else % Point 2 is in the middle of the trial
        t2 = intersect(find(cursorposi(:, 1) == p2(1)), find(cursorposi(:, 2) == p2(2)))  
    end
    

    %Swap timepoints if t1 comes after t2
    if t1 > t2 
        temp = t1;
        t1 = t2;
        t2 = temp;
        clear temp;
    end

    % Remove all error data between t1(1) and t2(1)
    scatter(cursorposi(t1(1):t2(1),1), cursorposi(t1(1):t2(1),2), 'red');  
    cursorposi(t1:t2, :) = nan; %This can possibly be removed now
    
    [~, T1] = min(abs(data(:, 3) - cursortimei(t1)));
    [~, T2] = min(abs(data(:, 3) - cursortimei(t2)));
    
    if endPoint
        data(T1:end, 4:8) = nan;
    else
        data(T1:T2, 4:8) = nan;
    end

    if type == 'A'
        trial_data.a{testInRun} = data;
    elseif type == 'B'
        trial_data.b{testInRun} = data;
    end
    load([save_dir runID], "stim", "time", "coordinates", "stim_order", "start_end", "links");
    save([save_dir runID ' - Errors Removed'], "coordinates", "links", "start_end", "stim", "stim_order", "time", "trial_data");
    disp('Saved! Continue removing errors or press ENTER to go to the next test')

end