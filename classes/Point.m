classdef Point<handle
    %POINT represents a point in space defined by its three carthesian
    %coordinates. Additional properties define the appearance and the
    %relationship of the point with other objects in a C3D file analysis

    properties
        Label (1,1) string%Label of the point
        XData (:,1) double%X coordinates
        YData (:,1) double%Y coordinates
        ZData (:,1) double%Z coordinates
        Color (:,3) double%Color in the plot
        Units (1,1) string {mustBeMember(Units,["m","mm"])}="mm"; %Length units
        SampleRate (1,1) double %Sample rate of the points
        Segment (1,1) string %Label of the Segment to which the point belongs
        Type (1,1) string %Type of the point (technical, virtual,...)
        Cluster (1,1) logical=1;%Flag to indicate if the point has to be used in cluster optimization
        Group (1,1) string      %Group of the Segment
        Parent Trial            %Trial containing the point
    end
    properties (Dependent)
        Coordinates
        Velocity
        Acceleration
    end

    properties (Access=private)
        NFrames    
    end

    methods
        function obj = Point(varargin)
            % POINT create an istance of point from a C3D file or
            % properties

            % if the input is a c3dfile
            if nargin==1
                if isnumeric(varargin{1})
                    try
                        [mkr,info]=btkGetMarkers(varargin{1});
                        sr=info.frequency;
                        units=info.units.ALLMARKERS;
                    catch ME
                        disp(ME.message)
                    end
                elseif isfile(varargin{1})
                    try
                        H=btkReadAcquisition(varargin{1});
                        [mkr,info]=btkGetMarkers(H);
                        sr=info.frequency;
                        units=info.units.ALLMARKERS;
                        btkCloseAcquisition(H);
                    catch ME
                        disp(ME.message)
                    end
                end

                mkr=structfun(@zero2nan,mkr,'UniformOutput',false); %change btk 0 to nan
                fns=fieldnames(mkr);
                if nargin<2
                    color=repmat([0 0 0],length(fns),1);
                end
                for i=length(fns):-1:1
                    obj(i).Label = fns{i};
                    obj(i).XData = mkr.(fns{i})(:,1);
                    obj(i).YData = mkr.(fns{i})(:,2);
                    obj(i).ZData = mkr.(fns{i})(:,3);
                    obj(i).Color = color(i,:);
                    obj(i).Units = char(units);
                    obj(i).SampleRate=sr;
                    obj(i).NFrames=size(obj(i).XData,1);
                end
                return
            end
            % if the input is point properties
            P=properties(obj);
            if length(P)*2<length(varargin)
                error('Too many input arguments!')
            end
            
            for j=1:2:length(varargin)
                for i=1:length(P)
                    if strcmp((P{i}),varargin{j})
                        obj.(P{i})=varargin{j+1};
                        break;
                    else
                        continue
                    end
                end
            end

        end
        %% POINT MODIFIER
        function obj=clearUnlabeled(obj,key)
        %CLEARUNLABELED removes unlabeled points specified by key
            arguments
                obj Point
                key char ='C_';
            end
            i_rmv=contains([obj.Label],key);
            obj(i_rmv)=[];
        end
        
        function obj=changeCoordinates(obj,T)
        %CHANGECOORDINATES applies the coordinate transformation specified
        %by T to all the points
        if any(not(size(T)==[4,4]))||not(isnumeric(T))||not(any(isfinite(T),'all'))
           error("T must be a finite double matrix of size 4x4x1!");
        end
            for i=1:length([obj.Label])
                xyz=[obj(i).Coordinates];
                xyz(:,4)=1;
                XYZ=(T*xyz')';
                obj(i).XData=XYZ(:,1);
                obj(i).YData=XYZ(:,2);
                obj(i).ZData=XYZ(:,3);
            end
        
        end

        function xyz=point2local(obj,T)
        %POINT2LOCAL transform the coordinates of a point
        arguments
            obj Point
            T (4,4,:) double
        end
        xyz=points2local(obj.Coordinates,T);         
        end

        function [obj,XYZ]=setGlobalCoords(obj,xyz,label,T)
        %SETGLOBALCOORDS transforms the coordinates of a point specified by
        %label to the system described by T
        arguments
            obj Point
            xyz (1,3) double
            label char
            T (4,4,:) double
        end
            XYZ=pagemtimes(T,[xyz 1]')';
            obj=obj.appendPoint(obj,XYZ,label);
        end
        
        function obj=setUnits(obj,newunits)
            % SETUNITS change the length units of the Point
            mustBeMember(newunits,["mm","m"]);
            for i=1:length(obj)
            if isequal(obj(i).Units,newunits)
            elseif isequal(obj(i).Units,"mm") && isequal(newunits,"m")
                obj(i).XData=obj(i).XData/1000;
                obj(i).YData=obj(i).YData/1000;
                obj(i).ZData=obj(i).ZData/1000;
            elseif isequal(obj(i).Units,"m") && isequal(newunits,"mm")
                obj(i).XData=obj(i).XData*1000;
                obj(i).YData=obj(i).YData*1000;
                obj(i).ZData=obj(i).ZData*1000;
            end
            obj(i).Units=newunits;
            end
            obj(i).Parent.Metadata.POINT.UNITS=newunits;
        end
        
        function obj=fillgaps(obj,method,NameValue)
            % FILLGAPS fill gaps in Points using spline
            arguments
                obj Point
                method char
                NameValue.MaxGap=20/100; % 20% max gap is allowed
                NameValue.AllowTailPrediction=false; %tails are removed
            end
            switch method
                case 'spline'
                    prop=namedargs2cell(NameValue);
                    obj.XData=time2cycle([],obj.XData,obj.NFrames,1,prop{:});
                    obj.YData=time2cycle([],obj.YData,obj.NFrames,1,prop{:});
                    obj.ZData=time2cycle([],obj.ZData,obj.NFrames,1,prop{:});
            end
        end
        
        function p=appendPoint(obj,label,XYZ,unit,segment,type,cluster,group)
            %APPENDPOINT adds a Point specified by label and XYZ coordinates to the object
            % obj=APPENDPOINT(obj,label,XYZ,segment,type,cluster,group) adds a
            % point to the end of the Point array. The new point need at least
            % a label and a set of XYZ coordinates to be created. Additional
            % specification provide context to the point.
            arguments
                obj Point
                label char
                XYZ (:,3) double
                unit string {mustBeMember(unit,["m","mm"])}="m";
                segment string="";
                type string ="";
                cluster logical=false;
                group string="";
            end
            % check if point already exist, if exist overwrites
            if ~any(matches([obj.Units],unit))
                if strcmp(unit,"m")
                    %convert from m to mm
                    XYZ = XYZ*1000;
                else
                    %convert from mm to m
                    XYZ = XYZ/1000;
                end
            end
            i=matches([obj.Label],label);
            if any(i)
                warning("A point in the acquisition has already the label %s!",label);
                i=find(i);
                p=obj(i);
                return

                sr=obj(i).SampleRate;
                tr=obj(i).Parent;
            else
                i=length(obj)+1;
                sr=obj(1).SampleRate;
                tr=obj(1).Parent;
            end
            % append point
            obj(i).Label=label;
            obj(i).XData=XYZ(:,1);
            obj(i).YData=XYZ(:,2);
            obj(i).ZData=XYZ(:,3);
            obj(i).Segment=segment;
            obj(i).Type=type;
            obj(i).Units=unit;
            obj(i).NFrames=size(XYZ,1);
            obj(i).SampleRate=sr;
            obj(i).Cluster=cluster;
            obj(i).Group=group;
            obj(i).Parent=tr;

            plist=[tr.ConfigFile.MarkerSet.Marker.("label"+tr.XMLatt)];
            if sum(matches(plist,label))<1 %and in case add it
                tr.ConfigFile.MarkerSet.Marker(end+1).("label"+tr.XMLatt)=label;
                tr.ConfigFile.MarkerSet.Marker(end).("segment"+tr.XMLatt)=segment;
                tr.ConfigFile.MarkerSet.Marker(end).("group"+tr.XMLatt)=group;
                tr.ConfigFile.MarkerSet.Marker(end).("type"+tr.XMLatt)=type;
                tr.ConfigFile.MarkerSet.Marker(end).("cluster"+tr.XMLatt)=double(cluster);
            end



            p=obj(i);
        end
        %% BASIC OPERATIONS
        function obj=sort(obj)
            % SORT sort points in alphabetic order
            [~,ord]=sort([obj.Label]);
            obj=obj(ord);
        end
        
        function obj=filtfilt(obj,b,a)
            %FILTFILT overloads filtfilt for class point
            try
            for i=1:length(obj)
            obj(i).XData=nanfiltfilt(b,a,obj(i).XData);
            obj(i).YData=nanfiltfilt(b,a,obj(i).YData);
            obj(i).ZData=nanfiltfilt(b,a,obj(i).ZData);
            end
            catch ME
                warning("%s: %s",obj(i).Label,ME.message);
            end

        end

        function obj=trim(obj,range)
            %TRIM cuts the point to the range specified by range
            obj.XData = obj.XData(:,range);
            obj.YData = obj.YData(:,range);
            obj.ZData = obj.ZData(:,range);
        end

        function obj=resample(obj,new_samplerate)
            %RESAMPLE sets the new sample rate using built-in resample
			for i=1:length(obj)
            p=obj(i).SampleRate;
            obj(i).XData = resample(obj(i).XData,p,new_samplerate);
            obj(i).YData = resample(obj(i).XData,p,new_samplerate);
            obj(i).ZData = resample(obj(i).XData,p,new_samplerate);
            obj(i).SampleRate=new_samplerate;
            obj(i).NFrames=obj(i).NFrames*new_samplerate/p;
			end
        end
		
		function obj=spline(obj,new_samplerate)
            %SPLINE use spline to increase the number of points to a
            %higher number
			for i=1:length(obj)
            p=obj(i).SampleRate;
			np=obj(i).NFrames;
			newp=np*new_samplerate/p;
			tb_old=linspace(1,np,np);
			tb_new=linspace(1,np,newp);
            obj(i).XData = spline(tb_old,obj(i).XData,tb_new);
            obj(i).YData = spline(tb_old,obj(i).YData,tb_new);
            obj(i).ZData = spline(tb_old,obj(i).ZData,tb_new);
            obj(i).SampleRate=new_samplerate;
            obj(i).NFrames=newp;
			end
        end

        function obj=mean(obj)
            %MEAN overloads the mean function
            for i=1:length(obj)
            obj(i).XData=mean(obj(i).XData,'omitnan');
            obj(i).YData=mean(obj(i).YData,'omitnan');
            obj(i).ZData=mean(obj(i).ZData,'omitnan');
            obj(i).SampleRate=nan;
            obj(i).NFrames=1;
            end
        end
        
        function [P, V, A]=PointStruct(obj)
        %POINTSTRUCT returns a structure of point position, velocity and acceleration with a field for each
        %element in label
            for i=1:length({obj.Label})
                label=matlab.lang.makeValidName(obj(i).Label,'Prefix','C_');
                P.(obj(i).Label)=obj(i).Coordinates;
                V.(obj(i).Label)=obj(i).Velocity;
                A.(obj(i).Label)=obj(i).Acceleration;
            end
        end

        function xyz=get.Coordinates(obj)
            xyz=[obj.XData,obj.YData,obj.ZData];
        end

        function V=get.Velocity(obj)
            V=[diff(obj.Coordinates,1,1); nan(1,3)]*obj.SampleRate;
        end
        
        function A=get.Acceleration(obj)
            A=[diff(obj.Coordinates,2,1); nan(2,3)]*(obj.SampleRate^2);
        end
        %% GRAPHICS

        function ln=plot(obj,frame,NameValue)
        %PLOT overloads the plot function for Point, represents the point in
        %the specified frame. If frame is a
        % vector, it shows the trace for the specified frames
            arguments
                obj Point
                frame=1;
                NameValue.?matlab.graphics.primitive.Line;
            end
            prop=namedargs2cell(NameValue);    
            x=[obj.XData];
            y=[obj.YData];
            z=[obj.ZData];
            x=x(frame,:);
            y=y(frame,:);
            z=z(frame,:);
            ln=plot3(x(:),y(:),z(:),'LineStyle','none','Marker','o',prop{:});
        end
        function txt=displayLabels(obj,ln,TextSpec)
        %DISPLAYLABELS shows the label of each point for an existing plot
        %and returns a text handle. If the point contains more than 1 frame
        %the text is displayed at the first frame
            arguments
                obj Point
                ln matlab.graphics.chart.primitive.Line
                TextSpec.?matlab.graphics.primitive.Text;
            end
            
            textspec=namedargs2cell(TextSpec);
            ax=ln.Parent;
            for i=length([obj.Label]):-1:1
            txt(i)=text(ax,ln.XData(i),ln.YData(i),ln.ZData(i),obj(i).Label,textspec{:});
            end
        end

        function updatehg(obj,frame,ln,txt)
        %UPDATEHG updates the graphics handles to a new frame. If frame is a
        % vector, it shows the trace for the specified frames
        arguments
            obj Point
            frame double
            ln matlab.graphics.chart.primitive.Line
            txt matlab.graphics.primitive.Text=text();
        end
            x=[obj.XData];
            y=[obj.YData];
            z=[obj.ZData];
            ln.XData=[x(frame,:)];
            ln.YData=[y(frame,:)];
            ln.ZData=[z(frame,:)];
        if nargin==4
            nt=length(txt);
            for i=1:nt
            txt(i).Position(1)=ln.XData(i);
            txt(i).Position(2)=ln.YData(i);
            txt(i).Position(3)=ln.ZData(i);
            end
        end
        end

        %% C3D I/0

        function updateC3D(obj,H,mode)
        %UPDATEC3D updates the data in the btk handle specified by H with the desired mode
        %"overwrite" overwrites existing data and adds new data
        %"append" ignores existing data and puts all the new points
        %"newfile" resets all the marker data in the acquisition
            arguments
                obj Point
                H double
                mode string {mustBeMember(mode,["overwrite","newfile","append"])}="overwrite"
            end
        
        btkSetPointsUnit(H,'marker',char(obj(1).Units));
        switch mode
            case "overwrite"
            for i=string(obj.Label)    
                try
                btkSetPoint(H,char(i),obj.PointStruct.(i));
                catch
                btkAppendPoint(H,'marker',char(i),obj.PointStruct.(i));
                end
            end
            case "append"
            for i=string([obj.Label])
                try
                btkAppendPoint(H,'marker',char(i),obj.PointStruct.(i));
                catch ME
                warning([ME.message, ', point won''t be created']);
                end
            end
        end %switch
        obj(1).Parent.readC3DMetaData;
        end %updateC3D

        function list(obj)
        % lists the point and their status
                fprintf("Index\t     Label\t        Segment\tGroup\tLabeling (%%)\n");
                i=1;
            for p=obj
                try
                cperc=100*sum(isfinite(p.XData))/p.NFrames;
                catch
                    cperc=nan;
                end
                fprintf("%5i\t%10s\t%15s\t%4s\t%10.0f\n",i,p.Label,p.Segment,p.Group,cperc)
                i=i+1;
            end
        end
    end %methods
end
