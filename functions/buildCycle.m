function [cycle, stance, swing] = buildCycle(in, type, context)
arguments
    in Event
    type string {mustBeMember(type, ["point", "analog"])} = "analog";
    context string {mustBeMember(context, ["Left", "Right", ["Left", "Right"]])} = ["Left", "Right"];
end
cycle = struct();
stance = struct();
swing = struct();
temp = in.exportEvents(type,false);
for i = 1:length(context)
    [cycle.(context(i)), stance.(context(i)), swing.(context(i))] = eventscycle(temp.(context(i)).Foot_Strike, temp.(context(i)).Foot_Off);
end


function [pt, p1, p2] = eventscycle(eventpre, eventpost)
%EVENTSCYCLE Cycle finding
%   [PT, P1, P2] = EVENTSCYCLE(EVENTPRE, EVENTPOST) finds the
%   pattern cycle P1-P2-P1 from the first array of events
%   EVENTPRE and the second one EVENTPOST.
%   It outputs the cycles PT and the sub-cycles P1, P2.
%   EVENTPRE and EVENTPOST must be column arrays.
%   Typically, the input arrays are footstrike and footoff, and the
%   outputs are cycle, contact, and swing periods.
%
%   See also c3dheader, c3dparameters.
% 
%   Author(s): F. Patanè, 2000

% Initialization
pt = [];
p1 = [];
p2 = [];

% Check for empty input
if isempty(eventpre) || isempty(eventpost)
    return;
end

% Ensure inputs are column vectors
eventpre = eventpre(:);
eventpost = eventpost(:);

% Combine and label events
et = [eventpre; eventpost]; % Total events
et = [et, [ones(size(eventpre)); zeros(size(eventpost))]];

% Sort events by time
[~, b] = sort(et(:, 1));
et = et(b, :); % Sorted events

% Find initial contact indices
ic = find(diff(diff(et(:, 2))) == 2); % ic: initial contact

% Extract time events
et = et(:, 1);

first = (ic:2:length(et));
last = (ic+1:2:length(et));
stop = min(length(first), length(last));


% Define output cycles
pt = [et(ic), et(ic + 2)];
p1 = [et(first(1:stop)), et(last(1:stop))];
p2 = [et(ic + 1), et(ic + 2)];
end
end


% ===========old version (less robust)
% pt=[eventpre eventpre];
% pt(1:end-1,2)=pt(2:end,2);
% pt(end,:)=[];
% p1=pt;
% for i=1:length(eventpre)-1
%     temp=find(eventpost>pt(i,1));
%     p1(i,2)=eventpost(temp(1));
% end
% p2=[p1(:,2),pt(:,2)];