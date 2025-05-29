classdef Scalar
    %SCALAR Contains scalar data of an acquisition

    properties
        Label (1,1) string
        Description (1,1) string
        Data (:,3) double
        Units (1,1) string
        SampleRate (1,1) double
        Parent (1,1) Trial
    end

    methods
        function obj = Scalar(varargin)
            %SCALAR Construct an instance of this class from properties or
            %a c3d file
            if nargin>1
               pp=string(properties(obj))';
               for i=1:nargin
                   obj.(pp(i))=varargin{i};
               end
                return
            elseif nargin==1
                c3dfile=varargin{1};
            else
                return
            end
            if isnumeric(c3dfile)
               [scalar,info]=btkGetScalars(c3dfile);               
            elseif isfile(c3dfile)
               H=btkReadAcquisition(c3dfile);
               [scalar,info]=btkGetScalars(H);
               btkCloseAcquisition(H);
            else
               return
            end
            
            fns=string(fieldnames(scalar))';

            if isempty(fns)
               obj=obj.empty;
               return
            end

            for i=length(fns):-1:1
                obj(i).Label=(fns(i));
                obj(i).Description="";%info.description.(fns(i));
                obj(i).Data=scalar.(fns(i));
                obj(i).Units=info.units.ALLSCALARS;
                obj(i).SampleRate=info.frequency;
            end
               
        end

        function obj = filtfilt(obj,b,a)
            %FILTFILT Overload filtfilt for scalar
            for i=1:length(obj)
                obj(i).Data=filtfilt(obj(i).Data,b,a);
            end
        end
        function obj = resample(obj,targetSR)
            %RESAMPLE overload resample for this class
            for i=1:length(obj)
                obj(i).Data=resample(obj(i).Data,obj(i).SampleRate,targetSR);
            end
        end
        function obj = mean(obj)
            %MEAN overloads the mean function
            for i=1:length(obj)
                obj(i).Data=mean(obj(i).Data);
            end
        end        
        function updateC3D(obj,H,mode)
            % UPDATEC3D updates the btk handle associated with the scalar
            % data

        switch mode
            case 'overwrite'
                for i=1:length(obj)
                btkAppendPoint(H,'scalar',char(obj(i).Label),obj(i).Data,obj(i).Data(:,1)*0,char(obj(i).Description));
                end
            case 'append'
                for i=1:length(obj)
                btkAppendPoint(H,'scalar',char(obj(i).Label),obj(i).Data,obj(i).Data(:,1)*0,char(obj(i).Description));
                end
         end
        end

        function [c,cmedio,cstd]=timeNorm(obj,npoints,base)
        %TIMENORM returns the scalar data normalized in stance or stride
        arguments
                obj Scalar
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