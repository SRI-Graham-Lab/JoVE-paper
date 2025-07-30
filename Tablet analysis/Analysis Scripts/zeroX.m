function ZXi = zeroX(t,y)
% theta = [0:7:360*4,1440]; % Angle Vector
% y = sind(theta);          % Signal
% y = data;

UpZCi = @(v) find(v(1:end-1) <= 0 & v(2:end) > 0);	% Returns Up Zero-Crossing Indices
DownZCi = @(v) find(v(1:end-1) >= 0 & v(2:end) < 0);    % Returns Down Zero-Crossing Indices
ZeroX = @(x0,y0,x1,y1) x0 - (y0.*(x0 - x1))./(y0 - y1); % Interpolated x value for Zero-Crossing

ZXi = sort([UpZCi(y),DownZCi(y)]);
ZX = ZeroX(t(ZXi),y(ZXi),t(ZXi+1),y(ZXi+1));

% === Checking for zero at the ignored value ===
if y(end)==0
    ZX(end+1) = t(end);
end
% ==============================================

% check if ZX: index of zero crossings
%
% figure(1)
% plot(t, y, '-b')
% hold on;
% plot(ZX,zeros(1,length(ZX)),'ro')
% grid on;
% legend('Signal', 'Interpolated Zero-Crossing')

end