function TableReport(file,OpenAfterProcessing,mass,ampside)

if not(exist('OpenAfterProcessing','var')) % if the user wants to see the excel
   OpenAfterProcessing=false;
end

if not(exist('file','var')) % then select one or more files
    [filelist, path]=uigetfile('*.c3d','MultiSelect','on');
    filelist=string(filelist);
elseif isfolder(file) %then start from the folder
    [filelist, path]=uigetfile(fullfile(file,'*.c3d'),'MultiSelect','on');
    filelist=string(filelist);
else % divide path and file
    [path,file,ext]=fileparts(file);
    file=[file,ext];
    filelist=string(file);
end

for file=filelist %loop through files
    
    file=char(file);
    H=btkReadAcquisition(char(fullfile(path,file)));
    md=btkGetMetaData(H);

    if not(exist('mass','var')) %if mass is not set
        try % to read from c3d
            mass=md.children.PROCESSING.children.mass.info.values;
        catch %enter manually
            mass= inputdlg('enter Subject mass (kg)');
            mass=str2double(mass);
        end
        if any(mass==[0,nan])
        error('mass must be a finite real number!');
        end
    end
if not(exist('ampside','var'))
try
Lamp=md.children.PROCESSING.children.LAmputation.info.values;
Ramp=md.children.PROCESSING.children.RAmputation.info.values;
if Lamp && Ramp
    ampside='Both';
elseif Lamp
    ampside='Left';
elseif Ramp
    ampside='Right';
else
    ampside='none';
end
catch
    ampside=string(listdlg("ListString",["Left", "Right", "Both", "none"],"PromptString","Amputation Side not Set, please select:","ListSize",[180 100],"SelectionMode","single"));
end
end

btkCloseAcquisition(H);

%% FILL INFO
exp='(?<subject>\w+)_(?<session>\w+)_T(?<tID>\d+)_(?<trial>\w+)'; %expression to match trial name
info=regexp(file,exp,'names'); %information retrieved from string

Info.Athlete.Mass=mass;
Info.Athlete.ID=info.subject;
Info.Athlete.AmputationSide=ampside;
try
Info.Session.Date=getTrialDate(fullfile(path,[file(1:end-4),'.Trial.enf']));
catch
    warning('enf file missing or date not available!')
    Info.Session.Date='n/a';
end
Info.Trial.ID=file(1:end-4);
Info.Trial.RunID=info.trial;

%% PROCESS TABLE
exmode='Ignore';
if contains(file,'JUMP')
    Info.Athlete.AmpLevel=regexprep(info.subject,'\d+','');
    Info.Trial.TakeOffLeg=Info.Athlete.AmputationSide;
    Info.Trial.JumpID=Info.Trial.RunID;
    [Tfile]=tableCoreProcessorLJ(fullfile(path,file),Info,fullfile(path,'..'),exmode);
else
    [Tfile,tab]=tableCoreProcessor(fullfile(path,file),Info,fullfile(path,'..'),exmode);
end
end

if OpenAfterProcessing
   winopen(Tfile)
end
end
%% SUBFUNCTIONS
function date=getTrialDate(enffile)
fid=fopen(enffile);
data=fread(fid,'*char')';
fclose(fid);
dc=regexp(data,'CREATIONDATEANDTIME=(?<y>\d+),(?<m>\d+),(?<d>\d+),\w+','names');
date=[dc.y,'/',dc.m,'/',dc.d];
end
