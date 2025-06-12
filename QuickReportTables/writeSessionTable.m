function Tfile=writeSessionTable(GeneralInfo,eventsForce,kinDB,lumped,sessiondir,existingdatamode,speed)
%Tfile=writeSessionTable(GeneralInfo,eventsForce,kinDB,lumped,sessiondir,existingdatamode,speed)
%creates/updates a Session Excel Spreadsheet containing the info, kinetics and/or lumped results
%for each step existing in the input trial
% the function splits the existing template table into three blocks: INFO, AMTI,
% LUMPED which can be filled/overwrited depending on user inputs
% Inputs (Required):
% GeneralInfo: structure built from GeneralInfo.txt file of the trial
% eventsForce: structure of Foot Strike and Foot Off events
% kinDB: struct containing forces, impulses and other kinetic measurements
% lumped: struct containing 1-D leg lumped model results
% Inputs (Optional):
% sessiondir: directory in to which save the .xlsx file (defaults to
% current directory)
% existingdatamode: defines how the function behaves when user tries to
%                   write data of a step which was already in the table 
%                   -Ignore: skip existing data entries, write missing
%                            entries (if any new)
%                   -Overwrite: overwrite the FULL row in which data exist
%                   -Append: move to the first empty row and write again
%                            the data
% speed: write down vicon speed for the trial
arguments
GeneralInfo struct
eventsForce struct
kinDB struct  =[];
lumped struct  =[];
sessiondir char =cd;
existingdatamode char {mustBeMember(existingdatamode,{'','Ignore','Overwrite','Append'})} ='';
speed double =nan;
end
%% set initial variables
Tdir=fullfile(cd,'Module14_Tables');
%sessiondir=uigetdir('Select Session directory');
hdrlines=4; % number of header rows in the layout table
def=existingdatamode; %needed to assign an initial value to def (default behavior for existingdatamode)
%% initialize missing input to empty
if isempty(kinDB)
   kinDB=struct('L',[],'R',[]);
   theta=lumped.Theta;
   stiffness=lumped.Stiffness;
elseif isempty(lumped)
   theta=struct('L',[],'R',[]);
   stiffness=struct('L',[],'R',[]);
elseif isempty(kinDB)&&isempty(lumped)
   theta=struct('L',[],'R',[]);
   kinDB=struct('L',[],'R',[]);
   stiffness=struct('L',[],'R',[]);
end

%% check layout file
layoutT='Layout quick reporting SSR.xlsx';
if ~exist(layoutT,'file')
   warning(['Template not found in ' Tdir,' Please select a new file']);
   [file, path]=uigetfile('.xlsx','Select Template');
   copyfile(fullfile(path,file),fullfile(Tdir,layoutT))
end
%% check session file
sessionID=GeneralInfo.session_info{1};
Tfile=[sessionID '_Quick Report Table.xlsx'];
Tfile=fullfile(sessiondir,Tfile);
if ~exist(Tfile,'file')
   copyfile(fullfile(Tdir,layoutT),Tfile);
   disp(['Session file created as: ' Tfile]);
end

Told=readtable(Tfile,'NumHeaderLines',hdrlines);
htable=height(Told);

%% fix inputs and divide steps
affside=GeneralInfo.amputation_side;
FSleft=eventsForce.Left_Foot_Strike;
FSright=eventsForce.Right_Foot_Strike;
nsteps=length(FSleft)+length(FSright);
if FSleft(1)<FSright(1)
    firststep='Left';
    secondstep='Right';
else
    firststep='Right';
    secondstep='Left';
end
if contains(affside,firststep,'IgnoreCase',true)
    AffSteps=strcat(firststep,' ',num2str((1:2:nsteps)'));
    UnaffSteps=strcat(secondstep,' ',num2str((2:2:nsteps)'));
else
    AffSteps=strcat(secondstep,' ',num2str((2:2:nsteps)'));
    UnaffSteps=strcat(firststep,' ',num2str((1:2:nsteps)'));
end

if contains(affside,'left')
    aside='L';
    uside='R';
else
    aside='R';
    uside='L';
end

%AffSteps(end,:)=[];

%% Initialize tables
Tvicon=table(nan,'VariableNames',{'speed'});
Tamti=table(nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,'VariableNames',{'meanGRFvAS','meanGRFvUS',...
    'minGRFhAS','maxGRFhAS','maxGRFvAS','minGRFhUS','maxGRFhUS','maxGRFvUS',...
    'stpropIhAS','brakIhAS','propIhAS','IvAS','stpropIhUS','brakIhUS','propIhUS','IvUS'});
Ttheta=table(nan,nan,nan,nan,nan,nan,'VariableNames',{'ThetaFS_AS','ThetaFO_AS','LGT_AS','ThetaFS_US','ThetaFO_US','LGT_US'});


for i=1:ceil(nsteps/2)

    [Tathl, Tprost, Ttrial]=infoTable(GeneralInfo);
    if i<=size(AffSteps,1)
        Ttrial.AffStepID=AffSteps(i,:);
    end

    if i<=size(UnaffSteps,1)
        Ttrial.UnaffStepID=UnaffSteps(i,:);
    end


    % cerco le righe che contengono eventuali dati duplicati
    try
        exrun=contains(Told.Var12,Ttrial.RunID,'IgnoreCase',true);
        exstep=contains(Told.Var13,Ttrial.AffStepID,'IgnoreCase',true) | contains(Told.Var14,Ttrial.UnaffStepID,'IgnoreCase',true);
        exrow=exrun & exstep;
    catch
        exrow=false;
    end

    if any(exrow) % se ce almeno una riga di duplicati entra nel if-else
        str=[string(['Data already exists in table for run ',strcat(GeneralInfo.trial_info{2:3}), ', steps : ' Told.Var13{i} ' and ' Told.Var14{i}]);
            "-Info";
            "-AMTI";
            "-Theta"];
        exrow=find(exrow,1);
        existingMD=true;
        existingAMTI=any(not(isnan(Told{exrow,16:31}))); %colonne di amti per la riga con i dati esistenti
        existingTH=any(not(isnan(Told{exrow,32:37})));   %colonne di theta per la riga con i dati esistenti
        str=str([true existingMD existingAMTI existingTH]);

        if isequal(existingdatamode,'')
            def='';
            existingdatamode=questdlg(str,'Action required','Ignore','Overwrite','Append','Ignore');
        end
        switch existingdatamode
            case 'Ignore'
                row=exrow;
                writeMD=false;
                writeAMTI=not(existingAMTI);
                writeTH=not(existingTH);
            case 'Overwrite'
                row=exrow;
                writeMD=true;
                writeAMTI=true;
                writeTH=true;
            case 'Append'
                row=htable+i;
                writeMD=true;
                writeAMTI=true;
                writeTH=true;
                
        end
    else
        row=htable+i;
        writeMD=true;
        writeAMTI=true;
        writeTH=true;       
    end
    %% fill table: affected
    if i<=size(AffSteps,1)
        if isempty(kinDB.(aside))
            Tamti{:,[1,3:5,9:12]}=nan;
        else
            Tamti.meanGRFvAS=kinDB.(aside).parameters.vertical.mean.cycle(i);
            Tamti.maxGRFvAS=kinDB.(aside).parameters.vertical.max.cycle(i);
            Tamti.minGRFhAS=kinDB.(aside).parameters.horizontal.min.cycle(i);
            Tamti.maxGRFhAS=kinDB.(aside).parameters.horizontal.max.cycle(i);
            Tamti.stpropIhAS=kinDB.(aside).impulse.startprop.cycle(i);
            Tamti.brakIhAS=kinDB.(aside).impulse.braking.cycle(i);
            Tamti.propIhAS=kinDB.(aside).impulse.propulsive.cycle(i);
            Tamti.IvAS=kinDB.(aside).impulse.vertical.cycle(i);

            grf(:,1)=kinDB.(aside).GRFz.cycle(:,i);
            grf(:,2)=kinDB.(aside).GRFy.cycle(:,i);
        end

        if isempty(theta.(aside))
            Ttheta{:,1:3}=nan;
        else
            Ttheta.ThetaFS_AS=theta.(aside).GT.cycle(1,i);
            Ttheta.ThetaFO_AS=theta.(aside).GT.cycle(end,i);
            Ttheta.LGT_AS=stiffness.(aside).GT.L0.all(i);
        end
    else
        grf(:,1)=nan(101,1);
        grf(:,2)=nan(101,1);
        Tamti{:,[1,3:5,9:12]}=nan;
        Ttheta{:,1:3}=nan;
    end

    %% fill table: unaffected
    if i<=size(UnaffSteps,1)
        if isempty(kinDB.(uside))
            Tamti{:,[2,6:8,13:16]}=nan;
        else
            Tamti.meanGRFvUS=kinDB.(uside).parameters.vertical.mean.cycle(i);
            Tamti.maxGRFvUS=kinDB.(uside).parameters.vertical.max.cycle(i);
            Tamti.minGRFhUS=kinDB.(uside).parameters.horizontal.min.cycle(i);
            Tamti.maxGRFhUS=kinDB.(uside).parameters.horizontal.max.cycle(i);
            Tamti.stpropIhUS=kinDB.(uside).impulse.startprop.cycle(i);
            Tamti.brakIhUS=kinDB.(uside).impulse.braking.cycle(i);
            Tamti.propIhUS=kinDB.(uside).impulse.propulsive.cycle(i);
            Tamti.IvUS=kinDB.(uside).impulse.vertical.cycle(i);

            grf(:,3)=kinDB.(uside).GRFz.cycle(:,i);
            grf(:,4)=kinDB.(uside).GRFy.cycle(:,i);
        end
        if isempty(theta.(uside))
            Ttheta{:,4:6}=nan;
        else
            Ttheta.ThetaFS_US=theta.(uside).GT.cycle(1,i);
            Ttheta.ThetaFO_US=theta.(uside).GT.cycle(end,i);
            Ttheta.LGT_US=stiffness.(uside).GT.L0.all(i);
        end
    else
        grf(:,3)=nan(101,1);
        grf(:,4)=nan(101,1);
        Tamti{:,[2,6:8,13:16]}=nan;
        Ttheta{:,4:6}=nan;
    end

    %% general
    Tvicon.speed=speed;
    range=['A' num2str(hdrlines+row)];
    if writeMD
        writetable(Tathl,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
        range(1)='E';
        writetable(Tprost,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
        range(1)='K';
        writetable(Ttrial,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
        range(1)='O';
        writetable(Tvicon,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
    end

    if not(isempty(kinDB.(aside)) && isempty(kinDB.(uside))) && writeAMTI
        range(1)='P';
        writetable(Tamti,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
    end

    if not(isempty(theta.(aside)) && isempty(theta.(uside))) && writeTH
        range=['AF' range(2:end)];
        writetable(Ttheta,Tfile,'Range',range,'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
    end

    if (not(isempty(kinDB.(uside))) || not(isempty(kinDB.(aside))))
        grfvarnames={'Sample','Affected GRFv','Affected GRFh','Unaffected GRFv','Unaffected GRFh'};
        Tgrf=array2table([(1:101)',grf(:,:)],'VariableNames',grfvarnames);
        grf=[];
        sname=['GRF ',Ttrial.RunID, ',StepPair_' num2str(i)];
        if ~any(contains(sheetnames(Tfile),sname))
            writetable(Tgrf,Tfile,'Sheet',sname);
        else
            switch existingdatamode
                case 'Overwrite'
                    writetable(Tgrf,Tfile,'Sheet',sname);
                case 'Append'
                    j=1;
                    sname=['GRF ',Ttrial.RunID, ',Step_', AffSteps(i,:), '&' UnaffSteps(i,:),' (',num2str(j),')'];
                    while any(sheetnames(Tfile)==sname)
                        j=j+1;
                        sname=['GRF ',Ttrial.RunID, ',Step_',  AffSteps(i,:), '&' UnaffSteps(i,:),' (',num2str(j),')'];
                    end
                    writetable(Tgrf,Tfile,'Sheet',sname);
                case 'Ignore'

            end
        end
    end
existingdatamode=def;
end

end

function [Tathl, Tprost, Ttrial]=infoTable(GeneralInfo)
% table maker
Tathl=table({'none'},0,{'none'},0,'VariableNames',{'Name','mass','AffSide','UnaffGTheight'});
Tprost=table({'none'},0,0,{'none'},0,0,'VariableNames',{'Model','Cat','AffGTheight','SocketConfigID','TAP','SAP'});
Ttrial=table(datetime('today','Format','yyyy-MM-dd'),{'none'},{'none'},{'none'},'VariableNames',{'Date','RunID','AffStepID','UnaffStepID'});

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
Ttrial.RunID=strrep(strjoin(GeneralInfo.trial_info(1:3),''),' ' ,'');
end