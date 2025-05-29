function [fp] = splitForcePlates(obj)
%SPLITFORCEPLATES Splits the forceplates in groups based on existing events
%   This function takes multiple forceplates as input and uses event data
%   to split the data into different groups (one for each event context
%   containing Foot_Strike and Foot_Off events)
ev=obj(1).Parent.Events.exportEvents("analog",false);
ev=rmfield(ev,"units");
context=string(fieldnames(ev)');
k=1;
for c=context
    try
    FC=ev.(c).Foot_Strike;
    FO=ev.(c).Foot_Off;
    mask=true(1,length(obj(1).GRF(:,1)));
    for i=1:length(FC)
        mask(FC(i):FO(i))=false;
    end
    forcep=obj;
    for f=1:numel(obj)
    forcep(f).Channels.Force(mask,:)=0;
    forcep(f).Channels.Moment(mask,:)=0;
    forcep(f).getGRF;
    end
    fp(k)=forcep.combineFP;
    fp(k).Label="FP_"+c;
    k=k+1;
    catch
    end
end
end