classdef Segment< matlab.mixin.Copyable
    %Segment is an object representing a body described by the position of
    %discrete points and a structure defining a coordinate system
    properties
        Label (1,1) string  %Unique label given to the segment
        Group (1,1) string  %Group of the segment
        Points (1,:) Point  %Points belonging to the segment
        Parent (1,1) Trial  %Source trial for the data
        Child  (1,1) string %Actually unused
        SampleRate (1,1) double=nan %Sample Rate of the data
        Type (1,1) string
        EndPoints (1,2) Point   %Points to which the segment is connected to other segments
        COM (1,:) Point         %Center of Mass of the segment
        Mass (1,1) double       %Mass of the segment
        Icm (3,3) double        %Baricentric moment of inertia in local coordinates
        InertialUnits (1,2) string %Units for Mass and Icm
        AngleUnits (1,1) {mustBeMember(AngleUnits,["rad","deg"])}="deg";
        %AngleAlias=AngleSequence; %Convention for names
        AngleSequence (1,3) char ='XYZ';
        CoordsysDef (1,1) struct %Structure containing the info to obtain the Transformation matrix from points
        RotationOffset=eye(3);
        SignConvention (1,3) double =[1 1 1]; %Convention for sign
        AngleAlias="";
    end
    properties (Access=private)
        ln matlab.graphics.chart.primitive.Line
        
    end
    properties (Dependent)
        TransformMat
        Orientation (3,3,:) double%=nan(3,3,1)
        Origin (:,3) double%=nan(1,3)%diventa un label di uno dei points(in caso da mettere dopo?)
        EulerAngles (:,3) double
        
        OriginSpeed (:,3) double
        OriginAcceleration (:,3) double
        AngularVelocity (:,3) double
        AngularAcceleration (:,3) double
        
    end
    methods
        function obj = Segment(varargin)
            %SEGMENT Construct an instance of this class

            P=properties(obj);
            if length(P)<length(varargin)
            error('Too many input arguments!')
            end

            for i=1:length(varargin)
            obj.(P{i})=varargin{i};
            end

        end
        
        %% CORE FUNCTIONS
        function T=get.TransformMat(obj)
            try
            T=buildCS(obj.Points.PointStruct,obj.CoordsysDef);
            catch
            T=nan(4,4,length(obj.Points(1).XData));
            end
        end

        function R=get.Orientation(obj)
            R=obj.TransformMat(1:3,1:3,:);
        end

        function O=get.Origin(obj)
            O=permute(obj.TransformMat(1:3,4,:),[3 1 2]);
        end

        function angle=get.EulerAngles(obj)
            angle=rotm2eul(pagemtimes(obj.RotationOffset(1:3,1:3),obj.Orientation),obj.AngleSequence);
            %[~,ord]=sort(double(char(obj.AngleSequence)));
            %ord=1:3;
            obj.Label
            if isequal(obj.AngleUnits,"deg")
            angle=angle*180/pi;
            end
            angle=angle.*obj.SignConvention;
            %angle=angle(:,ord);
        end
        % function alias=get.AngleAlias(obj)
        %     xyz=["Roll","Yaw","Pitch"];
        %     seq=obj.AngleSequence-'X';
        %     alias=xyz(seq+1);
        %     alias=strjoin(alias,", ");
        % end
        function omega=get.AngularVelocity(obj)
        %getAngularVelocity returns the angular velocity in the fixed or
        %moving frame and in the units specified by the user
        frame='fixed';
        %
        
        Rdot=diff(obj.Orientation,1,3);
        Rdot(:,:,end+1)=nan; %adds a fake last sample equal to nan
        if frame(1)=='f'
           S=pagemtimes(Rdot,'none',obj.Orientation,'transpose');
        else
           S=pagemtimes(obj.Orientation,'transpose',Rdot,'none');
        end
            omega(:,:)=permute([ S(3,2,:) S(1,3,:) S(2,1,:)],[3 2 1])*obj.SampleRate;
        if isequal(obj.AngleUnits,"deg")
            omega=omega*180/pi;
        end
            if all(isfinite(obj.Parent.KinematicFilter))
                b=obj.Parent.KinematicFilter(1,:);
                a=obj.Parent.KinematicFilter(2,:);
                omega=nanfiltfilt(b,a,omega);
            end
        end

        function alpha=get.AngularAcceleration(obj)
            omega=obj.AngularVelocity;
            alpha=diff(omega)*obj.SampleRate;
            alpha(end+1,:)=nan;
            if all(isfinite(obj.Parent.KinematicFilter))
                b=obj.Parent.KinematicFilter(1,:);
                a=obj.Parent.KinematicFilter(2,:);
                alpha=nanfiltfilt(b,a,alpha);
            end
        end

        function speed=get.OriginSpeed(obj)
        %getOriginSpeed returns the speed of the origin    
            speed=[diff(obj.Origin);nan(1,3)]*obj.SampleRate;
        end

        function acc=get.OriginAcceleration(obj)
        %getOriginSpeed returns the acceleration of the origin     
            acc=[diff(obj.Origin,2);nan(2,3)]*(obj.SampleRate^2);
        end
        
        function [angle]=strideNorm(obj,FS,npoints)
        for i=(length(FS)-1):-1:1
            try
        angle(:,:,i)=time2cycle([],obj.EulerAngles(FS(i):FS(i+1),:),npoints);
            catch
            end
        end
        end
        %% GRAPHICS
        function ln=show(obj,ax,frame,NameValue)
        %show Creates the Segment in the axes specified by ax for the specified
        % frame or set of frames
            arguments
                obj
                ax=axes();
                frame=1;
                NameValue.ScaleFactor=100;
                NameValue.?matlab.graphics.primitive.Line;
            end
            ii=length(obj);
            
            for i=length(obj):-1:1
            cT=obj(i).TransformMat;
            if isempty(cT)
            else
            T(:,:,:,ii)=cT;
            ii=ii-1;
            end
            end
            prop=namedargs2cell(NameValue);
            ln=plotframe(ax,T(:,:,frame,:),prop{:});  
        end

        function updatehg(obj,frame,ln)
        %updatehg update the graphics to the new frame
        arguments 
            obj
            frame=1
            ln=obj.show(gca,frame,'ScaleFactor',100);         
        end
          
        ii=length(obj);
        for i=length(obj):-1:1
        try
        T(:,:,:,ii)=obj(i).TransformMat;
        ii=ii-1;
        catch
        end
        end
        R=T(1:3,1:3,:,:);
        O=permute(T(1:3,4,:,:),[3 1 4 2]);
        for i=1:3
        V=vec2linedata(O(frame,:,:),R(i,:,frame,:),scale);
        ln(i).XData=V(:,1);
        ln(i).YData=V(:,2);
        ln(i).ZData=V(:,3);
        end
        end
        
        %% C3D I/O
        function updateC3D(obj,H,mode,jointsuffix)
        %updateC3D updates the data in the btk handle specified by H with the desired mode
        %"overwrite" overwrites existing data and adds new data
        %"append" ignores existing data and puts all the new angle
        %"newfile" resets all the angle data in the acquisition
            arguments
                obj Segment
                H double
                mode string {mustBeMember(mode,["overwrite","newfile","append"])}="overwrite"
                jointsuffix string ="Joint"
            end
        % Segment Angle
        for o=obj

        switch mode
            case "overwrite" 
                try
                btkSetPoint(H,char(o.Label),o.EulerAngles);
                catch
                btkAppendPoint(H,'angle',char(o.Label),o.EulerAngles,o.EulerAngles(:,1)*0,char(o.AngleAlias));
                end

            case "append"
                try
                btkAppendPoint(H,'angle',char(o.Label),o.EulerAngles,o.EulerAngles(:,1)*0,char(o.AngleAlias));
                catch ME
                warning([ME.message, ', point won''t be created']);
                end
            case "newfile"
                ang=btkGetAngles(H);
                 newlabel=[o.Label];
                for i=1:length(newlabel)
                    try
                    btkAppendPoint(H,'angle',char(o.Label),o.EulerAngles,o.EulerAngles(:,1)*0,char(o.AngleAlias));
                    catch ME
                    end
                end
        end
        end
        end
        
        %% TIME NORMALIZATION
        function [c,cmedio,cstd]=timeNorm(obj,npoints,base)
        %TIMENORM returns the angle, angular
        %velocity, power and moment in the joint, normalized cycles (npoints) for the stance or the stride
            arguments
                obj Segment
                npoints=101;
                base string {mustBeMember(base,["stance","stride"])}="stride"
            end
            ev=exportEvents(obj.Parent.Events,'point',false);
            group=obj.Group;
            if matches(group,["Left","Right"])
            else
            group=["Left", "Right"];
            %group = "Left";
            end
            S=[]; E=[];
            for g=group
                S=[S ev.(g).Foot_Strike];
                if base=="stride"
                    E=[E ev.(g).Foot_Strike(2:end)];
                    S(end)=[];
                else
                    E=[E ev.(g).Foot_Off];
                end
            
            % if length(group)<2
            %    name="angle";
            % else
            %     name="angle"+g;
            % end
                for i=(length(S)):-1:1
                    try
                    c.angle(:,:,i)=time2cycle([],obj.EulerAngles(S(i):E(i),:),npoints);
                    catch
                    end
                end
            end 
            cmedio=structfun(@(x) mean(x,3,'omitnan'),c,"UniformOutput",false);
            cstd=structfun(@(x) std(x,1,3,'omitnan'),c,"UniformOutput",false);
        end
    end
end