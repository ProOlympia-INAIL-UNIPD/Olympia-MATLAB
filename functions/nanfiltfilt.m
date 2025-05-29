function out = nanfiltfilt(b, a, in)
toomit = isnan(in(:, 1));
temp = diff(toomit);
start = find(temp == -1)+1;
stop = find(temp == 1);

out = nan(size(in));
if isempty(start) && isempty(stop)
    out = filtfilt(b,a, in);
elseif isempty(start) && ~isempty(stop)
    out(1:stop(1), :) = filtfilt(b, a, in(1:stop(1), :));
elseif ~isempty(start) && isempty(stop)
    out(start(1):end, :) = filtfilt(b, a, in(start(1):end, :));
else
    for i = 1:length(start)
        out(start(i):stop(i), :) = filtfilt(b, a, in(start(i):stop(i), :));
    end
end

end
%{
% inputmask=isnan(in);
% in(inputmask)=0;
% out=filtfilt(b,a,in);
% out(inputmask)=nan;
%return
toomit = isnan(in(:, 1));
temp = diff(toomit);
start = find(temp == -1)+1;
stop = find(temp == 1);

out = nan(size(in));
if isempty(start) && isempty(stop)
    out = filtfilt(b,a, in);
elseif isempty(start) && ~isempty(stop)
    out(1:stop(1), :) = filtfilt(b, a, in(1:stop(1), :));
elseif ~isempty(start) && isempty(stop)
    out(start(1):end, :) = filtfilt(b, a, in(start(1):end, :));
else
    for i = 1:length(start)
        out(start(i):stop(i), :) = filtfilt(b, a, in(start(i):stop(i), :));
    end
end

end
%}