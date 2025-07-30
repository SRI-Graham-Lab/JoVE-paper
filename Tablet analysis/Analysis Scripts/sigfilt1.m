% 1-D sigma filter. input has to be 1-D vector, square window width w have to be an odd number >= 3, K must <= (w-1)/2, sm is the smoothing indicator (logical) offering the option to skip smoothing/averaging in non-outlier case (M>K). output is the filtered signal, twosig is the signal within specified sigma cutoff. Zhongmin Lin. May 8, 2020.  
function [output,sigcut] = sigfilt1 (input,nsig,w,K,sm) 

% validating inputs
narginchk(1,5); 
if nargin == 1 % default filtering parameters
    nsig = 2; 
    w = 3; 
    K = 1;
    sm = true;
end
if isvector(input) == false
    error('input must be a vector')
end
wrange = 3:2:length(input);
if ismember(w,wrange) == false
    error(strcat('window width w must be an odd number from 3 to',num2str(wrange(end))))
end
w2 = floor(w/2);
w3 = ceil(w/2);
if K <= 0 || floor(K) ~= ceil(K)
    error('K should be a positive integer')
elseif K >= w3
    error(strcat('K must <',num2str(w3)))
end
if islogical(sm) == false
    error('Smoothing indicator should be a logical (true or false)')
end

% establish intensity range
xij = mean(input);
sigma = std(input);
delta = nsig * sigma; % 2 sigma
lb = xij - delta;
ub = xij + delta;

% get pixels within 2 sigma range, replace the rest with NaN
for j = 1:length(input)
    if input(j) >= lb && input(j) <= ub
        sigcut(j) = input(j);
    else
        sigcut(j) = NaN;
    end
end
% NaN padding
% nanpad(1:wdwidth-1,:) = NaN;
nanpad = padarray(sigcut,[w2 w-1],NaN,'both');

% move the window through the input array
for i = w3:length(nanpad)-w2
    
    % define window
    wd = nanpad(1:w,i-w2:i+w2);
    
    % sum and average within range pixels within the window
    sigsum = sum(wd,'all','omitnan');
    sigavg = mean(wd,'all','omitnan'); % sigavg include central pixel
    M = sum(~isnan(wd),'all'); % number of pixels within sigma range

    % Determine whether to exclude central pixel
    if M > K % many pixels within sigma range, include, smooth 
        if sm == true % if want smoothing
            output(i) = sigavg; % smooth/average            
        else % if don't want smoothing
            output(i) = wd(w3,w3); % leave pixel as it is, no averaging
        end
    else % M <= K. few pixels within sigma range, exclude, despike
        % calculate nb avg that exclude central pixel
        if isnan(wd(w3,w3)) == true % if central is outlier, then remove it from average don't make a difference
            output(i) = sigavg; % nbavg = sigavg
        else % if central is not outlier, remove it from nb avg calculation
            output(i) = (sigsum-wd(w3,w3))/(M-1);
        end        
    end
end
% output the trimmed and aligned data, same size as the input
output = output(w:end-w2);
end

% % original script with preset nsig,w,K. May 6, 2020
% xij = mean(spd2filt);
% sigma = std(spd2filt);
% % establish intensity range
% delta = 2 * sigma; % 2 sigma
% lb = xij - delta;
% ub = xij + delta;
% % sum all pixels lie within the range in a 1x3 window
% n=0;
% m=1;
% % 2 sigma range
% for j = 1:length(spd2filt)
%     if spd2filt(j) >= lb && spd2filt(j) <= ub
%         twosig(j) = spd2filt(j);
%     else
%         twosig(j) = NaN;
%     end
% end
% % NaN padding
% sigspd = [NaN;NaN;twosig';NaN;NaN];
% for i = 2:length(sigspd)-1
%     %     % pixels within range
%     %     if sigspd(i) >= lb && sigspd(i) <= ub
%     %         sigspd(i) = sigspd(i);
%     %     else
%     %         sigspd(i) = NaN;
%     %     end
%     % define window
%     wd = sigspd(i-1:i+1);
%     
%     % sum and average within range pixels within the window
%     sigsum = sum(wd,'omitnan');
%     sigavg = mean(wd,'omitnan'); % sigavg include central pixel
%     M = sum(~isnan(wd));
%     if isnan(sigspd(i)) == true % nbavg exclude central pixel
%         nbavg = sigavg;
%     else
%         nbavg = sigsum-wd(2)/(M-1);
%     end
%     K = 1;
%     if M > K
%         avgx(i) = sigavg;
%     else
%         avgx(i) = nbavg;
%     end
% end