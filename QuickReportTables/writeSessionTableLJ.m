function writeSessionTableLJ(GeneralInfo,kinDB,lumped,jumpdata,jumpside,sessiondir,existingdatamode)
arguments
GeneralInfo struct
kinDB struct  =[];
lumped struct  =[];
jumpdata struct  =[];
jumpside char {mustBeMember(jumpside,{'Left','Right',''})}='';
sessiondir char =cd;
existingdatamode char {mustBeMember(existingdatamode,{'','Ignore','Overwrite','Append'})} ='';
end

% initialize missing input to empty
if isempty(lumped)
   theta=struct('L',[],'R',[]);
   stiffness=struct('L',[],'R',[]);
else
   theta=lumped.Theta;
   stiffness=lumped.Stiffness;
end
if isempty(kinDB)
   kinDB=struct('L',[],'R',[]);
end


% identify jump side
if isempty(jumpside)
switch GeneralInfo.amputation_side
    case 'left'
    jumpside='L';
    case 'right'
    jumpside='R';
    otherwise
    error('Take-off Leg required for Bi-lateral and Able-Bodied athletes!');
end
end

Tdir=fullfile(cd,'Module14_Tables');
%sessiondir=uigetdir('Select Session directory');

hdrlines=4;
GeneralInfo.trial_info{2}='LJ';
%% check layout file
layoutT='Layout quick reporting LJ.xlsx';
if ~exist(layoutT,'file')
   warning(['Template not found in ' Tdir,' Please select a new file']);
   [file, path]=uigetfile('.xlsx','Select Template');
   copyfile(fullfile(path,file),fullfile(Tdir,layoutT))
end
%% check session file
sessionID=GeneralInfo.session_info{1};
Tfile=[sessionID '_Quick Report Table LJ.xlsx'];
Tfile=fullfile(sessiondir,Tfile);
if ~exist(Tfile,'file')
   copyfile(fullfile(Tdir,layoutT),Tfile);
   disp(['Session file created as: ' Tfile]);
end

Told=readtable(Tfile,'NumHeaderLines',hdrlines);
htable=height(Told);



%% Initialize tables
cycles=nan(101,6);
cycles(:,1)=0:100;
cyclevarnames={'Sample','VerticalGRF (Nm/BW)','HorizontalGRF (Nm/BW)','Trunk Lean (deg)','GTvert (mm)','GThor (mm)'};

Tcycles=array2table(cycles(:,:),'VariableNames',cyclevarnames);
Tgrf=Tcycles(:,2:3);
Tgt=Tcycles(:,4:6);

%% fill tables
[Tathl, Tprost, Ttrial]=infoTable(GeneralInfo);

try
exrows=find(contains(Told.Var12,Ttrial.JumpID),1);
catch
exrows=[];
end

if not(isempty(exrows))
       str=[string(['Data already exists in table for Trial: ',Ttrial.JumpID]);
            "-Info (Athlete, Session, and Trial ID)";
            "-video+AMTI (Jump length)";
            "-Vicon (Angles)"
            "-AMTI (GRF and Impulses";
            "-Vicon+AMTI (GT kinematics)"];

       existingMD=true;
       existingVideoAMTI=any(not(isnan(Told{exrows,15:17})));
       existingVicon=any(not(isnan(Told{exrows,18:33})));
       existingAMTI=any(not(isnan(Told{exrows,34:43})));
       existingViconAMTI=any(not(isnan(Told{exrows,44:51})));
       str=str([true existingMD existingVideoAMTI existingVicon existingAMTI existingViconAMTI]);
       if isequal(existingdatamode,'')
       existingdatamode=questdlg(str,'Action required','Ignore','Overwrite','Append','Ignore');
       end
       switch existingdatamode
           case 'Ignore'
           row=exrows;
           writeMD=false;
           writeVideoAMTI=not(existingVideoAMTI);
           writeVicon=not(existingVicon);
           writeAMTI=not(existingAMTI);
           writeViconAMTI=not(existingViconAMTI);
           case 'Overwrite'
           row=exrows;
           writeMD=true;
           writeVideoAMTI=true;
           writeVicon=true;
           writeAMTI=true;
           writeViconAMTI=true;
           case 'Append'
           row=htable+1;
           writeMD=true;
           writeVideoAMTI=true;
           writeVicon=true;
           writeAMTI=true;
           writeViconAMTI=true;
       end
else
    writeMD=true;
    writeVideoAMTI=true;
    writeVicon=true;
    writeAMTI=true;
    writeViconAMTI=true;
    row=htable+1;
end
row=num2str(hdrlines+row);
%% fill&write athlete table
if writeMD
range=['A' row];
writetable(Tathl,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
range(1)='E';
writetable(Tprost,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
range(1)='K';
writetable(Ttrial,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
end
%% fill&write videoAMTI table
TvideoAMTI=table(0,0,0,'VariableNames',{'FP2land','TakeOff2FP','JumpLength'});
%FP2land: distanza (misurata in campo) tra bordo anteriore pedana e atterraggio
%Toff2FP: distanza tra COP al takeoff e bordo anteriore pedana
%JumpLength: lunghezza effettiva salto

if not(isempty(jumpdata)) && writeVideoAMTI
    TvideoAMTI.FP2land=nan;
    TvideoAMTI.TakeOff2FP=jumpdata.deltaD;
    TvideoAMTI.JumpLength=nan;
range=['O' row];
writetable(TvideoAMTI,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);   
end

%% fill&write Vicon table
Tvicon=table(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,'VariableNames',{'GTh2ndlast','GThTD','GThTO','GTv2ndlast','GTvTD','GTvTO','thTrTD','thTrTO','thThTD','thShTD','hGT2ndlast','hGTTD','hGTTO','hLGTmax','hRGTmax','secGT'});
%'GTh2ndlast: media velocità del GT orizzontale nel penultimo passo
% GThTD: velocità orizzontale GT al touchdown
% GThTO: velocità orizzontale GT al TakeOff
% GTv2ndlast: come sopra ma verticale
% GTvTD','
% GTvTO','
% thTrTD: angolo Tronco al Touchdown
% thTrTO: angolo Tronco al TakeOff
% thThTD: angolo Tronco al Touchdown
% thShTD: angolo Tronco al Touchdown
% hGT2ndlast: altezza GT media nel penultimo passo
% hGTTD: altezza GT al Touchdown
% hGTTO
% hLGTmax
% hRGTmax

if not(isempty(jumpdata)) && writeVicon
Tvicon.GTh2ndlast=jumpdata.Vh_mean_2ndLast;% media velocità del GT orizzontale nel penultimo passo
Tvicon.GThTD=jumpdata.Vh_TD;% velocità orizzontale GT al touchdown
Tvicon.GThTO=jumpdata.Vh_TO;% velocità orizzontale GT al TakeOff
Tvicon.GTv2ndlast=jumpdata.Vv_mean_2ndLast;% come sopra ma verticale
Tvicon.GTvTD=jumpdata.Vv_TD;%
Tvicon.GTvTO=jumpdata.Vv_TO;%
Tvicon.thTrTD=jumpdata.trunk_angle_TD;% angolo Tronco al Touchdown
Tvicon.thTrTO=jumpdata.trunk_angle_TO;% angolo Tronco al TakeOff
Tvicon.thThTD=jumpdata.thigh_angle_TD;% angolo Coscia al Touchdown
Tvicon.thShTD=jumpdata.shank_angle_TD;% angolo Gamba al Touchdown
Tvicon.hGT2ndlast=jumpdata.H_GT_mean_2ndLast;% altezza GT media nel penultimo passo
Tvicon.hGTTD=jumpdata.H_GT_TD;% altezza GT al Touchdown
Tvicon.hGTTO=jumpdata.H_GT_TO;%
Tvicon.hLGTmax=jumpdata.H_LGT_max;%
Tvicon.hRGTmax=jumpdata.H_RGT_max;%
Tvicon.secGT=jumpdata.angleGTsec;




Tgt{:,1}=jumpdata.trunk_resampled; %trunk
Tgt{:,2}=jumpdata.GTy_resampled; %GTh
Tgt{:,3}=jumpdata.GTz_resampled; %GTz

range=['R' row];
writetable(Tvicon,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
%add COM(GT)-tip
range=['AX' row];
writetable(array2table([jumpdata.GT_tip_angle_TD, jumpdata.GT_tip_angle_TO]),Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
writecycle(Tfile,Ttrial,Tgt,'D1',existingdatamode)

end

%% fill&write AMTI
Tamti=table(0,0,0,0,0,0,0,0,0,0,'VariableNames',{'meanGRFbrak','meanGRFprop',...
    'meanGRFvert','maxGRFbrak','maxGRFprop','maxGRFvert','Istprop','Ibrak','Iprop','Ivert'});

if not(isempty(kinDB.(jumpside))) && writeAMTI
    Tamti.meanGRFbrak=nan;
    Tamti.meanGRFprop=nan;
    Tamti.meanGRFvert=kinDB.(jumpside).parameters.vertical.mean.cycle(end);

    Tamti.maxGRFbrak=kinDB.(jumpside).parameters.horizontal.min.cycle(end);
    Tamti.maxGRFprop=kinDB.(jumpside).parameters.horizontal.max.cycle(end);
    Tamti.maxGRFvert=kinDB.(jumpside).parameters.vertical.max.cycle(end);

    Tamti.Istprop=kinDB.(jumpside).impulse.braking.cycle(end);
    Tamti.Ibrak=kinDB.(jumpside).impulse.braking.cycle(end);
    Tamti.Iprop=kinDB.(jumpside).impulse.propulsive.cycle(end);
    Tamti.Ivert=kinDB.(jumpside).impulse.vertical.cycle(end);

    Tgrf{:,1}=kinDB.(jumpside).GRFz.cycle(:,end);
    Tgrf{:,2}=kinDB.(jumpside).GRFy.cycle(:,end);
    
    range=['AH' row];
    writetable(Tamti,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
    writecycle(Tfile,Ttrial,Tgrf,'B1',existingdatamode)
end

%% fill&write ViconAMTI
Tlumped=table(0,0,0,0,0,0,'VariableNames',{'lGT_TD','lGT_min','lGT_TO','theta_TD','theta_TO','LhGT'});

if not(isempty(theta.(jumpside)) && isempty(stiffness.(jumpside))) && writeViconAMTI
    Tlumped.lGT_TD=stiffness.(jumpside).GT.L0.all(end);
    Tlumped.lGT_min=stiffness.(jumpside).GT.Lmin.all(end);
    Tlumped.lGT_TO=stiffness.(jumpside).GT.Lfin.all(end);
    Tlumped.theta_TD=theta.(jumpside).GT.cycle(1,end);
    Tlumped.theta_TO=theta.(jumpside).GT.cycle(end,end);
    Tlumped.LhGT=stiffness.(jumpside).GT.deltaAPprox.all(end);
    range=['AR' row];
    writetable(Tlumped,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
end


end

function writecycle(Tfile,Ttrial,Tcycles,range,existingdatamode)

sname=['LastStepCycle ',Ttrial.JumpID];
if ~any(contains(sheetnames(Tfile),sname))
   writetable(array2table((1:101)','VariableNames',{'Sample'}),Tfile,'Sheet',sname);
   writetable(Tcycles,Tfile,'Range',range,'Sheet',sname);    
else
switch existingdatamode
    case 'Overwrite' 
        writetable(array2table((1:101)','VariableNames',{'Sample'}),Tfile,'Sheet',sname);
        writetable(Tcycles,Tfile,'Range',range,'Sheet',sname);
    case 'Append'
        j=1;
        sname=['LastStepCycle ',Ttrial.JumpID,' (',num2str(j),')']; 
        while any(sheetnames(Tfile)==sname)
        j=j+1;
        sname=['LastStepCycle ',Ttrial.JumpID,' (',num2str(j),')'];   
        end
        writetable(array2table((1:101)','VariableNames',{'Sample'}),Tfile,'Sheet',sname);
        writetable(Tcycles,Tfile,'Range',range,'Sheet',sname);      
    case 'Ignore'
end
end

end

function [Tathl, Tprost, Ttrial]=infoTable(GeneralInfo)
% table maker
Tathl=table({'none'},0,{'none'},0,'VariableNames',{'Name','mass','AffSide','UnaffGTheight'});
Tprost=table({'none'},0,0,{'none'},0,0,'VariableNames',{'Model','Cat','AffGTheight','SocketConfigID','TAP','SAP'});
Ttrial=table(datetime('today','Format','yyyy-MM-dd'),{'none'},{'none'},{'none'},'VariableNames',{'Date','JumpID','NrunUP','LrunUP'});

%% fill table: common
% athlete
Tathl.Name={GeneralInfo.athlete_name};
Tathl.mass=str2double(GeneralInfo.session_info{2});
Tathl.AffSide={GeneralInfo.amputation_side};

% prostesis
catsep=strfind(GeneralInfo.foot{1},'cat');
Tprost.Model=GeneralInfo.foot{1}(1:catsep-1);
Tprost.Cat=str2double(GeneralInfo.foot{1}(catsep+3:end));
Tprost.SocketConfigID=GeneralInfo.socket{1};
Tprost.TAP=str2double(GeneralInfo.foot{2});
Tprost.SAP=str2double(GeneralInfo.foot{3});

% trial
Ttrial.Date=datetime(GeneralInfo.session_info{3},'InputFormat','yyyyMMdd','Format','dd/MM/yyyy');
Ttrial.JumpID=strjoin(GeneralInfo.trial_info(1:3),'_');
Ttrial.NrunUP=nan;
Ttrial.LrunUP=nan;

end