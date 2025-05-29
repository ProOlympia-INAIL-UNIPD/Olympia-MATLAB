function out=makeQuickReportTable(Info,stepside,kinDB,theta,stiffness,speed)
%Tfile=writeSessionTableV2(Info,kinDB,lumped,stepside,Tfile,existingdatamode,speed)
%creates/updates a Session Excel Spreadsheet containing the info, kinetics and/or lumped results
%for each step existing in the input trial
% the function splits the existing template table into three blocks: INFO, AMTI,
% LUMPED which can be filled/overwrited depending on user inputs
% Inputs (Required):
% Info: structure with fields: Athlete, Session, Trial
%       Athlete contains ID, mass, Amputation Side
%       Session contains ID, Date
%       Trial contains ID of the run and ID of the trial
% kinDB: struct containing forces, impulses and other kinetic measurements
% lumped: struct containing 1-D leg lumped model results
% stepside: cell array containing the order of step laterality
% Inputs (Optional):
% Tfile: directory in to which save the .xlsx file (defaults to
% current directory), or specific destination file
% existingdatamode: defines how the function behaves when user tries to
%                   write data of a step which was already in the table 
%                   -Ignore: skip existing data entries, write missing
%                            entries (if any new)
%                   -Overwrite: overwrite the FULL row in which data exist
%                   -Append: move to the first empty row and write again
%                            the data
% speed: the average speed during the trial
arguments
Info struct   =[];
stepside cell {mustBeMember(stepside,{'L','R','Left','Right'})}={'L'};
kinDB struct  =[];
theta struct  =[];
stiffness struct  =[];
speed double =nan;
end

%% initialize missing input to empty
if isempty(kinDB)
   kinDB=struct('Left',[],'Right',[]);
end

if isempty(theta)
   theta=struct('Left',[],'Right',[]);
end

if isempty(stiffness)
   stiffness=struct('Left',[],'Right',[]);
end

%% Initialize Table 
Tvicon=table(nan,'VariableNames',{'speed'});
Tamti=table(nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,'VariableNames', ...
    {'meanGRFv','minGRFh','maxGRFh','maxGRFv',...
    'stpropIh','brakIh','propIh','Iv','I80','I90','I','contactTime'});
Ttheta=table(nan,nan,nan,'VariableNames',{'ThetaFS','ThetaFO','LGT'});
[Tathl, Tprost, Ttrial]=infoTable(Info);

out.Parameters=table();    
%% loop through steps

for i=1:length(stepside)
            side=stepside{i};
            ii=count(strjoin(stepside(1:i)),side);
            Ttrial.StepID=i;
            Ttrial.StepSide={side};

            cycles=nan(101,8); %initialize cycles (grf, cop, GT)
            if not(isempty(kinDB.(side))) %if AMTI is available create the table
            Tamti.meanGRFv=kinDB.(side).parameters.vertical.mean.cycle(ii);
            Tamti.maxGRFv=kinDB.(side).parameters.vertical.max.cycle(ii);
            Tamti.minGRFh=kinDB.(side).parameters.horizontal.min.cycle(ii);
            Tamti.maxGRFh=kinDB.(side).parameters.horizontal.max.cycle(ii);
            Tamti.stpropIh=kinDB.(side).impulse.startprop.cycle(ii);
            Tamti.brakIh=kinDB.(side).impulse.braking.cycle(ii);
            Tamti.propIh=kinDB.(side).impulse.propulsive.cycle(ii);
            Tamti.Iv=kinDB.(side).impulse.vertical.cycle(ii);
            Tamti.I80=nan;
            Tamti.I90=nan;
            Tamti.I=nan;
            Tamti.contactTime=kinDB.(side).ContactTime(ii);
            cycles(:,1)=kinDB.(side).GRF.cycle(:,2,ii);
            cycles(:,2)=kinDB.(side).GRF.cycle(:,3,ii);
            cycles(:,3)=vecnorm(cycles(:,1:2),2,2);
            cycles(:,4:5)=kinDB.(side).COP.cycle(:,2:-1:1,ii);
            end
            
            try
             %if lumped is available update table
            Ttheta.ThetaFS=theta.(side).GT.cycle(1,ii);
            Ttheta.ThetaFO=theta.(side).GT.cycle(end,ii);
            Ttheta.LGT=stiffness.(side).L0.all(ii);
            cycles(:,6:7)=stiffness.(side).PP.cycle(:,[2 3],ii);
            catch
                warning('%s Kinematics unavailable or corrupted!',side)
            end
            
                b=cycles(:,4:5)-cycles(:,6:7);
                b(:,3)=0;
                F=cycles(:,1:2);
                F(:,3)=0;
                Mhip=cross(b,F);
                cycles(:,8)=Mhip(:,3);

        Tvicon.speed=speed; %add speed of the trial
        cyclesvarnames={'GRFh','GRFv','GRFmod','COPh','COPl','GTh','GTv','Mhip'};
        Tcycles=array2table(cycles,'VariableNames',cyclesvarnames);
        sname=[char(Ttrial.RunID{:}) '_Cycles_S' num2str(i),char(side)];
        out.Parameters(i,:)=[cell2table({Info.Trial.ID},'VariableNames',"TrialID") Tathl Tprost Ttrial Tvicon Tamti Ttheta];
        out.(sname)=[array2table((1:101)','VariableNames',{'% Stance'}) Tcycles];
end %loop end

end

function [Tathl, Tprost, Ttrial]=infoTable(Info)
%initialize tables
Tathl=table({'none'},0,{'none'},0,'VariableNames',{'Name','mass','AffSide','UnaffGTheight'});
Tprost=table({'none'},0,0,{'none'},0,0,'VariableNames',{'Model','Cat','AffGTheight','SocketConfigID','TAP','SAP'});
Ttrial=table(datetime('today','Format','yyyy-MM-dd'),{''},0,{''},'VariableNames',{'Date','RunID','StepID','StepSide'});


%% fill table: common
% athlete
try
Tathl.Name={Info.Athlete.ID};
Tathl.mass=double(Info.Athlete.Mass);
Tathl.AffSide={Info.Athlete.AmputationSide};
catch
    warning('No Athlete Info Available')
end

%med device, if any
if not(isempty(Info.Athlete.AmputationSide))
try
Tprost.Model=Info.MedDev.Foot.Model;
Tprost.Cat=Info.MedDev.Foot.CAT;
Tprost.SocketConfigID=Info.MedDev.Socket.ID;
Tprost.TAP=Info.MedDev.Foot.TAP;
Tprost.SAP=Info.MedDev.Socket.SAP;
catch
    warning('No prostesis data detected!')
end
end

%trial
Ttrial.Date={Info.Session.Date};
Ttrial.RunID={Info.Trial.RunID};

end