classdef Analog
    %ANALOG Contains Analog signals
    %   This class contains Analog data and information

    properties
        Label (1,1) string      %Label of the Analog
        Description (1,1) string%Description of the Analog
        Data (:,1) double       %Data contained in the Analog
        Units (1,1) string      %Units assigned to the Analog channel
        SampleRate (1,1) double %Sample Rate of the Analog Channel
        Parent (1,1) Trial      %Trial containing the Analog Channel
    end

    methods
        function obj = Analog(c3dfile)
            %ANALOG Construct an instance of this class
            %   Use a C3D file/handle to retrieve analog channels. Analog
            %   channels relative to forceplate data are excluded and
            %   stored in forcePlatformType2 object.
            % See Also FORCEPLATFORMTYPE2
            if not(exist('c3dfile','var'))
                return
            end
            if isnumeric(c3dfile)
               [an,info]=btkGetAnalogs(c3dfile);
               fpdata=btkGetMetaData(c3dfile,'FORCE_PLATFORM');
               fpanalogs=fpdata.children.CHANNEL.info.values;
               
            elseif isfile(c3dfile)
               H=btkReadAcquisition(c3dfile);
               [an,info]=btkGetAnalogs(H);
               fpdata=btkGetMetaData(H,'FORCE_PLATFORM');
               fpanalogs=fpdata.children.CHANNEL.info.values;
               btkCloseAcquisition(H);
            else
               return
            end

            fns=string(fieldnames(an))';
            fns(fpanalogs(:))=[];

            if isempty(fns)
               obj=obj.empty;
               return
            end

            for i=length(fns):-1:1
                obj(i).Label=info.label.(fns(i));
                obj(i).Description=info.description.(fns(i));
                obj(i).Data=an.(fns(i));
                obj(i).Units=info.units.(fns(i));
                obj(i).SampleRate=info.frequency;
            end
               
        end

        function obj = filtfilt(obj,b,a)
            %FILTFILT Overload the filtfilt function
            %   Detailed explanation goes here
            for i=1:length(obj)
                obj(i).Data=filtfilt(obj(i).Data,b,a);
            end
        end
        function obj = resample(obj,targetSR)
            %RESAMPLE Downsample by an integer ratio
            %   Detailed explanation goes here
            for i=1:length(obj)
                obj(i).Data=resample(obj(i).Data,obj(i).SampleRate,targetSR);
            end
        end
        function obj = mean(obj)
            %MEAN Overload of the mean function
            %   Detailed explanation goes here
            for i=1:length(obj)
                obj(i).Data=mean(obj(i).Data);
            end
        end        
        function updateC3D(obj,H,mode)
            % UPDATEC3D Update the C3D handle

        switch mode
            case 'newfile'
              fpdata=btkGetMetaData(H,'FORCE_PLATFORM');
               fpanalogs=fpdata.children.CHANNEL.info.values;
               nan=btkGetAnalogNumber(H);
                for i=1:nan
                    if any(fpanalogs==i)
                    else
                    btkRemoveAnalog(H,i);
                    end
                end
                for i=1:length(obj)
                btkAppendAnalog(H,char(obj(i).Label),obj(i).Data,char(obj(i).Description));
                end
                if obj(1).NSamples==1
                   btkSetAnalogSampleNumberPerFrame(H,1);%the average file could be written
                end

            case 'overwrite'
                for i=1:length(obj)
                btkAppendAnalog(H,char(obj(i).Label),obj(i).Data,char(obj(i).Description));
                end
            case 'append'
                for i=1:length(obj)
                btkAppendAnalog(H,char(obj(i).Label),obj(i).Data,char(obj(i).Description));
                end
         end
    end
function [c,cmedio,cstd]=timeNorm(obj,npoints,base)
        % TIMENORM Normalize the Analog data based on time events
        %[c,cmedio,cstd]=timeNorm(obj,npoints,base) returns the angle, angular
        %velocity, power and moment in the joint, normalized cycles (npoints) for the stance or the stride
        %See Also EVENT
            arguments
                obj Analog
                npoints=101;
                base string {mustBeMember(base,["stance","stride"])}="stride"
            end
            ev=exportEvents(obj.Parent.Events,'point',false);

            group=["Left", "Right"];
            
            for g=group
                S=[]; E=[];
                S=[S ev.(g).Foot_Strike];
                if base=="stride"
                    E=[E ev.(g).Foot_Strike(2:end)];
                    S(end)=[];
                else
                    E=[E ev.(g).Foot_Off];
                end
            
            
            for i=(length(S)):-1:1
                try
                c.(g+"Data")(:,:,i)=time2cycle([],obj.Data(S(i):E(i),:),npoints);
                catch
                end
            end
            end 
            cmedio=structfun(@(x) mean(x,3,'omitnan'),c,"UniformOutput",false);
            cstd=structfun(@(x) std(x,1,3,'omitnan'),c,"UniformOutput",false);
        end
    end
end