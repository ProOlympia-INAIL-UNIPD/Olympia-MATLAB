function [Tfile]=writeSessionTableXLSX(Itable,Tfile,existingdatamode,RJselector)
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
Itable struct  =[];
Tfile char =cd;
existingdatamode char {mustBeMember(existingdatamode,{'','Ignore','Overwrite','Append'})} ='';
RJselector char {mustBeMember(RJselector,{'TSSR','JUMP'})} ='TSSR';
end
%% trial type selection
if RJselector(1)=='T'
hdrlines=3; % number of header rows in the layout table
layoutT='Layout quick reporting SSR_v2.xlsx'; %name of the template
ds=16;
else
hdrlines=4; % number of header rows in the layout table
layoutT='Layout quick reporting LJ.xlsx'; %name of the template
ds=19;
end
%needed to assign an initial value to def (default behavior for existingdatamode)
%% initialize missing input to empty
%% check layout file
Tdir=mfilename('fullpath'); %the template should be in the same folder of this function
Tdir=fileparts(Tdir);

if ~exist(layoutT,'file') % if template is missing, allows users to browse for it
   warning(['Template not found in ' Tdir,' Please select a new file']);
   [file, path]=uigetfile('.xlsx','Select Template');
   copyfile(fullfile(path,file),fullfile(Tdir,layoutT))
end
%% check session file
filename='Quick Report Tablev2.xlsx'; %base file name
    try % to add the session date as a suffix
        subsess=Itable.Parameters.TrialID{1};
        subsess=regexp(subsess,'(\w+_?\w+_)T\d+_\w+','tokens');
        filename=[char(subsess{:}), filename];
    catch    
    end

if isfolder(Tfile) %if the input file is a directory, create the excel
Tfile=fullfile(Tfile,filename);
else
end
if ~exist(Tfile,'file') % if the file doesn't exist, create a fresh one from template
   copyfile(fullfile(Tdir,layoutT),Tfile);
   sprintf(join("Session file created as: ", Tfile));
end
warning off
Told=readtable(Tfile,'NumHeaderLines',hdrlines); %reads the existing table (which is empty if newly created)
warning on
if all(size(Told)==[0 0]) %empty excel
   copyfile(fullfile(Tdir,layoutT),Tfile);
   sprintf(join("Session file created as: ", Tfile));    
   Told=readtable(Tfile,'NumHeaderLines',hdrlines);
end
htable=height(Told);
Tnew=Itable.Parameters;
wid=width(Tnew);
Told=Told(:,1:wid);
Told.Properties.VariableNames=Itable.Parameters.Properties.VariableNames;

%% Initialize Table 
    
if htable>0
    oldrows=erase(strcat(string(Told.TrialID),string(Told.RunID),string(Told.StepID),string(Told.StepSide)),' ');

    for i=height(Tnew):-1:1
        wsteps=erase(strcat(Tnew.TrialID(i),string(Tnew.RunID(i)),string(Tnew.StepID(i)),Tnew.StepSide(i)),' ');
        exsteps=find(contains(oldrows,wsteps))';
        
    for j=exsteps % if any data already exist, check more in detail
        str=['Data already exists in table for Trial ',Told.RunID{i}, ', step' string(Told.StepID(i)) ' (' Told.StepSide{i} ')'];
        
        oldDATA=isfinite(Told{j,ds:wid});  %dati già presenti
        newDATA=isfinite(Tnew{i,ds:wid});         %dati che possono essere aggiunti
        writeData=newDATA & xor(oldDATA,newDATA);%i dati devono essere scrivibili (&) ma non già presenti (xor)
        writeLine=[false(1,ds-1) writeData];
        Told(j,writeLine)=Tnew(i,writeLine);%aggiungo i dati non presenti alla tabella
        if any(oldDATA & newDATA) %gestisco i dati se c'e sovrapposizione
           writeData= oldDATA & newDATA;
           writeLine=[false(1,ds-1) writeData];

            if isequal(existingdatamode,'') % if you didn't set a default action, ask the user what to do
               existingdatamode=questdlg(str,'Action required','Ignore','Overwrite','Append','Ignore');
            end
            switch existingdatamode
                case 'Ignore'
                    Tnew(i,:)=[]; %posso cancellare la riga
                case 'Overwrite'
                    Told(j,writeLine)=Tnew(i,writeLine); %sovrascrivo
                case 'Append'              
            end
        else
        Tnew(i,:)=[]; %posso cancellare la riga
        end 
    end  
    end
    Told=[Told;Tnew];
else
    Told=Tnew;
end

    writetable(Told,Tfile,'Range',['A',num2str(hdrlines+1)],'AutoFitWidth',false,'WriteVariableNames',false,'PreserveFormat',true);
    
    %% cycle sheets
    Itable=rmfield(Itable,'Parameters');
    fns=fieldnames(Itable);
    exsheet=sheetnames(Tfile);
    for i=1:length(fns)
        newdata=Itable.(fns{i});
        if any(contains(exsheet,fns{i}))
           olddata=readtable(Tfile,'Sheet',fns{i},'VariableNamingRule','preserve');
           oldvar=olddata.Properties.VariableNames;
           newvar=newdata.Properties.VariableNames;
           exvar=find(contains(oldvar,newvar));

           for j=exvar
               old=olddata.(oldvar{j});
               new=newdata.(oldvar{j});
               if (all(abs((old-new)./old)<1e-3)||all(isnan(new)))%se i dati sono uguali o comunque non disponibili
               elseif all(isnan(old))
                    olddata.(oldvar{j})=newdata.(oldvar{j});
               elseif all(isfinite(old+new)) % entrambi i dati sono presenti ma diversi
                    if isequal(existingdatamode,'') % if you didn't set a default action, ask the user what to do
                       str=[fns{i} 'Sheet already contains data for some of the desired variables'];
                       existingdatamode=questdlg(str,'Action required','Ignore','Overwrite','Append','Ignore');
                    end
                    switch existingdatamode
                        case 'Ignore'
                        case 'Overwrite'
                            olddata.(oldvar{j})=newdata.(oldvar{j});
                        case 'Append'
                            k=1;
                            fns{i}=[fns{i} '(',int2str(k) ')'];
                            while any(contains(exsheet,fns{i})) && k<10
                            k=k+1;
                            fns{i}(end-1)=int2str(k);
                            end
                            olddata.(oldvar{j})=newdata.(oldvar{j});
                    end %switch
               end %if equal
           end %exvars
        else
            olddata=newdata;
        end
    writetable(olddata,Tfile,'Sheet',fns{i});
    end

    return


