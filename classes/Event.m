classdef Event
    %Event is a class containing properties and methods for Events in a
    %biomechanical acquisition.
    %   Events are stored in a c3d file usually divided in General, Left
    %   and Right context. Events can be represented in seconds, point, or
    %   analog samples. Time can be relative to the start of the
    %   acquisition or to the first frame of the trim acquisition.

    properties
        Parent Trial        %Reference to the Trial
        General struct      %Context for General Events
        Left struct         %Context for Left Events
        Right struct        %Context for Right Events
        Units string {mustBeMember(Units,["seconds","analog","point"])}="seconds";
        TimeOfFirstFrame=0; %Time between the start of the trial and the first indexed frame
        FirstFrame=1;
        PointRate=1;    %Sample Rate for points
        AnalogRate=1;   %Sample Rate for analog signals
        Absolute logical =false; %indicates wheter event are refered to the TimeOfFirstFrame or to the Absolute time of the acquisition
    end

    methods

        function obj = Event(c3dfile)
            %obj = Event(c3dfile) reads event from a c3dfile
            %   Detailed explanation goes here
        obj.General=struct();
        obj.Left=struct();
        obj.Right=struct();
        %obj.c3dfile=c3dfile;
        obj.Units="seconds";
        closeflag=true;
        if isnumeric(c3dfile)
            H=c3dfile;
            closeflag=false;
        elseif isfile(c3dfile)
            H=btkReadAcquisition(c3dfile);
        else
            [c3dfile,p]=uigetfile('*.c3d');
            H=btkReadAcquisition(fullfile(p,c3dfile));
        end
        obj.PointRate=btkGetPointFrequency(H);
        obj.AnalogRate=btkGetAnalogFrequency(H);
        obj.TimeOfFirstFrame=(btkGetFirstFrame(H)-1)/btkGetPointFrequency(H); %t=(frame-firstframe)/point_freq
        obj.FirstFrame=btkGetFirstFrame(H);
        [ev,info]=btkGetEvents(H);
        sub=info.subjects;
        fns=string(fieldnames(ev))';
        for cev=fns
            label=char(strrep(cev,strrep(sub.(cev)," ","_"),""));
            if label(1)=='_'
                label(1)=[];
            end
            lowbar=strfind(label,'_');
            if isstruct(obj.(label(1:lowbar(1)-1)))
                if isfield(obj.(label(1:lowbar(1)-1)),(label(lowbar(1)+1:end)))
                   e=obj.(label(1:lowbar(1)-1)).(label(lowbar(1)+1:end));
                   obj.(label(1:lowbar(1)-1)).(label(lowbar(1)+1:end))=sort([e,ev.(cev)]);
                   continue
                end
            end
            obj.(label(1:lowbar(1)-1)).(label(lowbar(1)+1:end))=ev.(cev);
        end
        obj.Absolute=true;
        if closeflag
        btkCloseAcquisition(H);
        end
        end
        
        %% EVENT I/O
        function obj = appendEvent(obj,context,label,value,mode)
            %  APPENDEVENT add event to the list
            % append Events (in seconds) to the current structure with info specified by
            % context, label, value, units, mode, absolute.
            arguments
                obj Event
                context char {mustBeMember(context,{'General','Left','Right'})}
                label char
                value {mustBeNumeric(value)}
                mode char {mustBeMember(mode,{'append','replace','merge'})}='append';
            end
            %warning('Using this function could result in unexpected behavior, please use appendEventStruct instead!')

            switch mode
                case 'append'
                    if isfield(obj.(context),label)
                    obj.(context).(label)=sort([obj.(context).(label) value]);
                    else
                    obj.(context).(label)=sort(value);
                    end
                case 'replace'
                    obj.(context).(label)=sort(value);
                case 'merge'
                    if isfield(obj.(context),label)
                    obj.(context).(label)=sort([obj.(context).(label) value]);
                    else
                    obj.(context).(label)=sort(value);
                    end
                    obj=obj.mergeDuplicates;
            end

        end
        function obj = appendEventStruct(obj,evstruct,mode,absolute)
            % APPENDEVENTSTRUCT add set of events from a structure
            % reads event from the input structure and parses them into the
            % event object. The input structure must contain a field named
            % units containing the units of the events.
            arguments
                obj Event
                evstruct struct
                mode char {mustBeMember(mode,{'append','overwrite'})}='append';
                absolute logical =false;
            end
            
                units=evstruct.units;
                evstruct=rmfield(evstruct,'units');
                if not(absolute) &&obj.Absolute
                    t0=obj.TimeOfFirstFrame;
                elseif absolute && not(obj.Absolute)
                    t0=-obj.TimeOfFirstFrame;
                else
                    t0=0;
                end
                
                switch units
                    case "seconds"
                        fun=@(x) sort(x+t0);
                    case "analog"
                        fs=obj.AnalogRate;
                        fun=@(x) sort((x-1)/fs+t0);
                    case "point"
                        fs=obj.PointRate;
                        fun=@(x) sort((x-1)/fs+t0);
                    otherwise
                        error('invalid units, must be; seconds, analog, or point');
                end

           
            for context=string(fieldnames(evstruct))'
                mustBeMember(context,["General","Left","Right"]);
                for label=string(fieldnames(evstruct.(context)))'
                    switch mode
                        case 'append'
                                if isfield(obj.(context),label)
                                    obj.(context).(label)=sort([obj.(context).(label) fun([evstruct.(context).(label)])]);
                                else
                                    obj.(context).(label)=fun([evstruct.(context).(label)]);
                                end  
                        case 'overwrite'
                                    obj.(context).(label)=fun([evstruct.(context).(label)]);     
                    end
                end
            end

                obj=obj.mergeDuplicates;


            obj.updateC3D(obj.Parent.C3DHandle,mode)
            obj.Parent.setC3DMetaData();
        end

        function eventStruct=exportEvents(obj,outputunits,absolute)
        % EXPORTEVENTS exports the events to a structure
        % each events is separated into a field for each context. Output units
        % could be seconds, point or analog and relative to the start of
        % the trial or of the trim region.
            arguments
                obj
                outputunits="seconds";
                absolute=false;
            end
            if absolute && obj.Absolute
                t0=0;
            else
                t0=obj.TimeOfFirstFrame;
            end
            eventStruct.units=outputunits;
        switch outputunits
            case "seconds"
            fun=@(t) t-t0;
            case "analog"
            fs=obj.AnalogRate;
            fun=@(t) round((t-t0)*fs+1);
            case "point"
            fs=obj.PointRate;
            fun=@(t) round((t-t0)*fs+1);
            case "perc"
                fun = @(context) stridePerc(obj, context);
            otherwise
            error('Units not valid!')
        end

        for context=["General","Left","Right"]
            if outputunits == "perc"
                eventStruct.(context) = fun(context);
            else
                eventStruct.(context)=structfun(fun,obj.(context),"UniformOutput",false);
            end
        end
        end

        function out = stridePerc(obj, context)
            % Custom processing for Foot_Strike and Foot_Off in the Left context
            out = struct();

            if isfield(obj.(context), 'Foot_Strike') && isfield(obj.(context), 'Foot_Off')
                % array alternating 0 and 100 for Foot_Strike
                n = length(obj.(context).Foot_Strike);
                out.Foot_Strike = repmat([0, 100], 1, ceil(n / 2));
                out.Foot_Strike = out.Foot_Strike(1:n); % Trim to correct length

                % compute Foot_Off
                footStrike = obj.(context).Foot_Strike;
                footOff = obj.(context).Foot_Off;
                out.Foot_Off = zeros(1, length(footOff)); % Preallocate output

                for i = 1:length(footOff)
                    if i < length(footStrike)
                        out.Foot_Off(i) = footOff(i) / (footStrike(i + 1) - footStrike(i));
                    else
                        out.Foot_Off(i) = NaN; % Handle last element if no next Foot_Strike
                    end
                end
            else
                warning('Required fields Foot_Strike and Foot_Off are missing in %s.', context);
            end
        end
        %% EVENT MODIFIER
        function obj = mergeDuplicates(obj)
            % check the object for duplicate events (i.e., elements having
            % the same label and separated by time less than the resolution
            % of the system)
            tol=1/obj.AnalogRate;
            for context=["Left","Right","General"]
            try
            obj.(context)=structfun(@(x) x([true; diff(x)>tol]),obj.(context),"UniformOutput",false);
            catch
            end
            end
        end
        
        function obj=detectfromForcePlates(obj,mode,threshold,timebase,startdetectionoffset)
        %DETECTFROMFORCEPLATES use force platform data to detect Foot
        %Strike and FootOff
        %to search FootStrike and FootOff events. COP data is used to
        %divide left and right contacts.
            arguments
                obj Event
                mode char {mustBeMember(mode,{'append','overwrite'})}='append';
                threshold=20;
                timebase string {mustBeMember(timebase,["analog","point"])}="analog";
                startdetectionoffset=0;
            end
            fp=obj.Parent.ForcePlatform;
            fp=fp.cleanSignal;
            fp=fp.combineFP;

            [FC, FO]=fp.getEvents("samples",threshold,startdetectionoffset);
            COP=fp.COP;
            for i=length(FC):-1:1
            COPm(i,:)=mean(COP(FC(i):FO(i),:));
            end
            
            R=pca(COPm); %this should align the COP to the running/walking direction
            if isempty(R)
                warning('Not Enough Foot Contacts for PCA analysis')
            else
                COPm=COPm*R;% first column is now Saggital (max variation), second ML, third Vertical
            end
            avML=mean(COPm(:,2));
            COPm=COPm(:,2)-avML;
            lr=COPm<0;
            %FC=ceil(FC/(obj.AnalogRate/obj.PointRate));
            %FO=floor(FO/(obj.AnalogRate/obj.PointRate));
            events.Left.Foot_Strike=FC(lr)';
            events.Left.Foot_Off=FO(lr)';
            events.Right.Foot_Strike=FC(~lr)';
            events.Right.Foot_Off=FO(~lr)';
            events.units=timebase;
            obj=obj.appendEventStruct(events,mode,false);
        end

        function obj=setUnitsandOffset(obj,newunits,absolute)
        %DONT USE THIS METHOD, EXPORT EVENTS IF NEED FOR DIFFERENT UNITS
        %obj=setUnitsandOffset(obj,newunits,absolute) changes the internal
        %units and the offset in which events are stored
            arguments
                obj Event
                newunits string {mustBeMember(newunits,["seconds","point","analog"])}
                absolute logical =obj.Absolute
            end

            error("This method is not available!")
        
                if (newunits==obj.Units)&&(absolute==obj.Absolute)
                    return
                end
                if not(absolute) &&obj.Absolute
                    t0=-obj.TimeOfFirstFrame;
                elseif absolute && not(obj.Absolute)
                    t0=obj.TimeOfFirstFrame;
                else
                    t0=0;
                end
                t0=t0-1/obj.PointRate;
                
                switch newunits
                    case "seconds"
                        mult=1;
                        t0=t0+1/obj.PointRate;
                    case "analog"
                        mult=1/obj.AnalogRate;
                    case "point"
                        mult=1/obj.PointRate;
                    otherwise
                        error('invalid units, must be; seconds, analog, or point');
                end

                 switch obj.Units
                     case "seconds"
                         oldmult=1;
                         t0=t0+1/obj.PointRate;
                     case "analog"
                         oldmult=1/obj.AnalogRate;
                     case "point"
                         oldmult=1/obj.PointRate;
                 end

            obj.Units=newunits;
            obj.Absolute=absolute;

            for i=["General","Left","Right"]
                try
                    if isequal(newunits,"seconds")
                obj.(i)=structfun(@(x) (x*oldmult+t0)/mult, obj.(i),'UniformOutput',false);
                    else
                obj.(i)=structfun(@(x) round((x*oldmult+t0)/mult)+1, obj.(i),'UniformOutput',false);
                    end
                catch
                end
            end
        
        end

        function nev=getEventCount(obj)
            nev=0;
            for ctx=["General","Left","Right"]
                fns=string(fieldnames(obj.(ctx)))';
                for f=fns
                    nev=nev+numel(obj.(ctx).(f));
                end
            end
        end

        %% C3D I/O

        function updateC3D(obj,H,mode)
        %UPDATEC3D updates the data in the btk handle specified by H with the desired mode
        %"overwrite" overwrites existing data and adds new data
        %"append" ignores existing data and puts all the new points
        %"newfile" resets all the marker data in the acquisition
            arguments
                obj Event
                H double
                mode string {mustBeMember(mode,["overwrite","newfile","append"])}="overwrite"
            end
            % if mode=="newfile"
            %    btkClearEvents(H);
            % end
            EV=btkGetEvents(H);

            for context=["General","Left","Right"]
                for ev=string(fieldnames(obj.(context)))'
        switch mode
            case "overwrite"
                
                for i=1:length(obj.(context).(ev))
                    try
                        if any(EV.(context).(ev)==obj.(context).(ev)(i))
                        else
                           btkAppendEvent(H,char(strrep(ev,'_',' ')),obj.(context).(ev)(i),char(context));
                           %btkAppendEvent(H, char(strrep(ev,'_',' ')), obj.(context).(ev)(i), char(context), 'MAFIO', 'bel micetto', 1)
                        end
                    catch
                           btkAppendEvent(H,char(strrep(ev,'_',' ')),obj.(context).(ev)(i),char(context));
                           %btkAppendEvent(H, char(strrep(ev,'_',' ')), obj.(context).(ev)(i), char(context), 'MAFIO', 'bel micetto', 1)

                    end
                        
                end
            case "append"
                for i=1:length(obj.(context).(ev))
                btkAppendEvent(H,char(strrep(ev,'_',' ')),obj.(context).(ev)(i),char(context));
                end
        end %switch
        
                end 
            end
            obj.Parent.readC3DMetaData;
        end %updateC3D
        
        function writetoC3D(obj,c3dfile,mode)
            % writes the events to the specified c3dfile
            h=btkReadAcquisition(c3dfile);
            obj.updateC3D(h,mode);
            btkWriteAcquisition(h,c3dfile);
            btkCloseAcquisition(h);
        end
    end %methods
end %class