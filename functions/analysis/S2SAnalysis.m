function s2s=S2SAnalysis(trial)
%S2SANALYSIS Computes Step-2-Step transition model.
% The function computes the mechanical work of the COM considering the
% leading and trailing leg model [Donelan, et al. (2002) J Biomech
% 35:117-124]. Work values are stored as Trial metadata, while
% instantaneous powers are stored as Trial Scalars.

%% validation
arguments
    trial (1,1) Trial
end

switch trial.TaskType
    case "walking"
        disp("Step-2-Step transition analysis...");
    otherwise
        return
end

forceplat=trial.ForcePlatform;
if numel(forceplat)==0
    warning("Trial Contains no Force Platforms!, Step-2-Step analysis can't be performed!");
    return
end


if getEventCount(trial.Events)==0 %needed because metadata isn't updated (btk bug?)
    warning("Trial Contains no events!, Step-2-Step analysis can't be performed!");
    return
else
    events=trial.Events;
end

COM = trial.Points(trial.findByLabel("Points","COM"));
if size(COM, 2)==0
    warning("Trial contains no COM! Step-2-Step analysis can't be performed!");
    return
end

%% calculate mechanical power when foot on plate
k = 1 + strcmp(COM.Units,"mm")*999; %if COM in mm, transform in m
m = trial.Metadata.ANTROPOMETRY.mass;

for kc=trial.ConfigFile.KinematicModel.KinematicChain
    %GRF(length(forceplat)).(kc.("group"+trial.XMLatt)) = struct("value", []);
    try
        s=kc.("group"+trial.XMLatt);
        ev_analog =events.exportEvents("analog");
        for fp = 1:length(forceplat)
            evnames = fieldnames(ev_analog.(kc.("group"+trial.XMLatt)));
            onoff = zeros(trial.NSamples, 1);
            for j = 1:length(ev_analog.(kc.("group"+trial.XMLatt)).(evnames{2}))
                onoff(ev_analog.(kc.("group"+trial.XMLatt)).(evnames{1})(j): ev_analog.(kc.("group"+trial.XMLatt)).(evnames{2})(j)) = 1;
            end

            temp = [trial.ConfigFile.MarkerSet.Marker.("label"+trial.XMLatt)];
            platefoot = nan(trial.NFrames, 3, 2);
            for j = 1:2
                platefoot(:, : ,j) = trial.Points(trial.findByLabel("Points", temp(strcmpi([trial.ConfigFile.MarkerSet.Marker.("eventDetection"+trial.XMLatt)], strrep(strrep(evnames(j),"_","")," ","")) &...
                    strcmpi([trial.ConfigFile.MarkerSet.Marker.("group"+trial.XMLatt)], s) == 1))).Coordinates;
            end

            footin = footInPlate(trial, forceplat(fp), platefoot);

            onoff = and(onoff, footin);
            up = find(diff(onoff) == 1);
            down = find(diff(onoff) == -1);

            % Delete contacts lasting less than 0.5 s
            for ii = 1:min(length(up), length(down))
                if down(ii) - up(ii) < 0.5*forceplat(fp).SampleRate
                    onoff(up + 1: down) = 0;
                end
            end

            GRF.(s).value(:,:,fp) = forceplat(fp).GRF .* repmat(onoff, 1, 3);
        end
        GRF.(s) = resample(sum(GRF.(s).value, 3, "omitmissing"), trial.NFrames, trial.NSamples);
        P.(s) = dot(GRF.(s), COM.Velocity/k, 2)/m;
        if exist("GRFon", "var")
            GRFon = GRFon + GRF.(s);
        else
            GRFon = GRF.(s);
        end
    catch
    end
end

%% CONSIDER APPENDING THE MECHANICAL POWER to the c3d Scalars

% figure;
% plot(P.(context{1}))
% hold on
% xline(s2s_events.lead.Left, 'r-')
% xline(s2s_events.trail.Right(1), 'g--')
% plot(s2s_events.trail.Right(1), 0, 'go','LineWidth',2)
% xline(int_instant, 'k--')
% plot(int_instant, 0, 'kx','LineWidth',2)
% yline(0, 'k', 'LineWidth', 1.5)
% xline(s2s_events.trail.Right(2), 'Color', "#D55E00")
% plot(s2s_events.trail.Right(2), 0, 'MarkerFaceColor', [0.84, 0.37, 0], 'o', 'LineWidth',2)
% plot(P.(context{2}), 'm-.')

%%
GRFon = abs(GRFon(:,3)) > 0;
GRFoff = find(diff(GRFon) == -1) + 1;
GRFon = find(diff(GRFon) == 1) + 1;

events.selectFootContacts();
ev_point = events.exportEvents("point");

% delete all the events before GRFon and after GRFoff
context = fieldnames(ev_point);
context = context(contains(context,["left","right"],'IgnoreCase',true));

for s=1:length(context)
    ev = fieldnames(ev_point.(context{s}));
    for e = 1:length(ev)
        ev_point.(context{s}).(ev{e}) = ev_point.(context{s}).(ev{e})(and(ev_point.(context{s}).(ev{e}) >= GRFon, ev_point.(context{s}).(ev{e}) <= GRFoff));
    end
    if ev_point.(context{s}).Foot_Off(1) < ev_point.(context{s}).Foot_Strike(1)
        ev_point.(context{s}).Foot_Off(1) = [];
    end
end

% check how many lead-trail couples
% Left
for i = 1:numel(ev_point.(context{1}).Foot_Strike)
    for j = 1:min(numel(ev_point.(context{2}).Foot_Strike), numel(ev_point.(context{2}).Foot_Off))
        try
            if and(ev_point.(context{1}).Foot_Strike(i) < ev_point.(context{2}).Foot_Off(j), ...
                    ev_point.(context{1}).Foot_Strike(i) < ev_point.(context{2}).Foot_Strike(j+1))
                lead = context{1};
                trail = context{2};
                %
                s2s_events.lead.(lead) = ev_point.(context{1}).Foot_Strike(i);
                s2s_events.trail.(trail) = [ev_point.(context{2}).Foot_Off(j), ev_point.(context{2}).Foot_Strike(j+1)];
            end
        catch
        end
    end
end

% Right
for i = 1:numel(ev_point.(context{2}).Foot_Strike)
    for j = 1:min(numel(ev_point.(context{1}).Foot_Strike), numel(ev_point.(context{1}).Foot_Off))
        try
            if and(ev_point.(context{2}).Foot_Strike(i) < ev_point.(context{1}).Foot_Off(j), ...
                    ev_point.(context{2}).Foot_Strike(i) < ev_point.(context{1}).Foot_Strike(j+1))
                lead = context{2};
                trail = context{1};
                %
                s2s_events.lead.(lead) = ev_point.(context{2}).Foot_Strike(i);
                s2s_events.trail.(trail) = [ev_point.(context{1}).Foot_Off(j), ev_point.(context{1}).Foot_Strike(j+1)];
            end
        catch
        end
    end
end

% cut the instantaneous mechanical power to the leading and trailing legs
context = fieldnames(P);
for c = 1:numel(context)
    try
    %if numel(s2s_events.(c).(lt)) > 1
        %P.("lead").(context{c}) = P.(context{c})(s2s_events.);
        %P.("trail").(context{~matches(context, context{c})}) = P.(context{~matches(context, context{c})});
        idx = find(diff(sign(P.(context{c}))) ~= 0);
        [~, temp] = min(abs(idx - s2s_events.("trail").(context{~matches(context, context{c})})(1)));
        int_instant = idx(temp);
        
        W.unit = "J/kg";
        W.DS.lead.context(c) = string(context{c});
        W.DS.lead.value(c) = trapz(P.(context{c})(s2s_events.("lead").(context{c}):int_instant,:))/trial.Metadata.POINT.RATE;

        W.SS.lead.context(c) = string(context{c});
        W.SS.lead.value(c) = trapz(P.(context{c})(int_instant:s2s_events.("trail").(context{~matches(context, context{c})})(2),:))/trial.Metadata.POINT.RATE;

        W.DS.trail.context(c) = string(context{~matches(context, context{c})});
        W.DS.trail.value(c) = trapz(P.(context{~matches(context, context{c})})(s2s_events.("lead").(context{c}):s2s_events.("trail").(context{~matches(context, context{c})})(1)))/trial.Metadata.POINT.RATE;

        W.step.net(c) = trapz(P.(context{c})(s2s_events.("lead").(context{c}):s2s_events.("trail").(context{~matches(context, context{c})})(2)))/trial.Metadata.POINT.RATE + ... %Plead
                     trapz(P.(context{~matches(context, context{c})})(s2s_events.("lead").(context{c}):s2s_events.("trail").(context{~matches(context, context{c})})(2)))/trial.Metadata.POINT.RATE; %Ptrail

        W.step.abs(c) = trapz(abs(P.(context{c})(s2s_events.("lead").(context{c}):s2s_events.("trail").(context{~matches(context, context{c})})(2))))/trial.Metadata.POINT.RATE + ... %Plead
                     trapz(abs(P.(context{~matches(context, context{c})})(s2s_events.("lead").(context{c}):s2s_events.("trail").(context{~matches(context, context{c})})(2))))/trial.Metadata.POINT.RATE; %Ptrail
    catch
    end
end

trial.Metadata.S2S=W;
%trial.setC3DMetaData;
end
