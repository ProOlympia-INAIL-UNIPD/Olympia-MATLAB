classdef Joint
    %JOINT contains the data of a Joint in a biomechanical
    %acquisition/model

    properties
        Label               %Label of the Joint
        Group (1,1) string  %Group of the Joint
        Parent              %Parent Segment for the joint
        Child               %Child Segment for the joint
        Type                %Type of Joint (Free, Hinge, ...)
        AngleSequence       %Sequence of angles for Euler/Cardan
        AngleUnits (1,1) {mustBeMember(AngleUnits,["rad","deg"])}="deg"; %Unit of measure for the angles
        JointCenter (1,1) Point % 3D Point corresponding to the Joint Center
        Force (:,3) double  % Joint Force (frames,XYZ) in the Parent Reference System 
        ForceUnits="N";     % Units of Force
        Moment (:,3) double % Joint Moment (frames,XYZ) in the Parent Reference System  
        MomentUnits (1,1) {mustBeMember(MomentUnits,["Nm","Nmm"])}="Nmm"; %Units of Torque
        Acquisition Trial   % Trial containing the data for the Joint
        
    end
    properties (Dependent)
        Rotation (3,3,:) double         %Rotation matrix of the Joint R=RP*RC'
        CardanAngles (:,3) double       %Cardan Angles (frames,AngleSequence)
        AngularVelocity (:,3) double    %Angular Velocity of the Child Segment in the Parent Reference System
        AngularAcceleration (:,3) double%Angular Acceleration of the Child Segment in the Parent Reference System
        Power (:,3) double              %Joint Power
        
    end
    properties (Access=private)
        LabelSuffix (1,1) string="Joint"
        PowerUnits (1,1)="W";
    end
    methods
        function obj = Joint(varargin)
        % JOINT Construct an istance of Joint
            P=properties(obj);
            if length(P)<length(varargin)
            error('Too many input arguments!')
            end

            for i=1:length(varargin)
            obj.(P{i})=varargin{i};
            end
            
        end
        %% DEPENDENT PROPERTIES
        function R=get.Rotation(obj)
            prox=obj.Parent.Orientation;
            dist=obj.Child.Orientation;
            R=pagemtimes(prox,'transpose',dist,'none');
        end

        function angle=get.CardanAngles(obj)
            angle=rotm2eul(obj.Rotation,obj.AngleSequence);
            if isequal(obj.AngleUnits,"deg")
            angle=angle*180/pi;
            end
            %[~,ord]=sort(double(char(obj.AngleSequence)));
            %angle=angle(:,ord);
        end

        function omega=get.AngularVelocity(obj)
        %getAngularVelocity returns the angular velocity in the fixed or
        %moving frame and in the units specified by the user
        frame='fixed';
        Rdot=diff(obj.Rotation,1,3);
        Rdot(:,:,end+1)=nan; %adds a fake last sample equal to nan
        if frame(1)=='f'
           S=pagemtimes(Rdot,'none',obj.Rotation,'transpose');
        else
           S=pagemtimes(obj.Rotation,'transpose',Rdot,'none');
        end
            omega(:,:)=[ S(3,2,:) S(1,3,:) S(2,1,:)]*obj.Parent.SampleRate;
            omega=omega';
        if isequal(obj.AngleUnits,"deg")
            omega=omega*180/pi;
        end
        end

        function alpha=get.AngularAcceleration(obj)

            omega=obj.AngularVelocity;
            alpha=diff(omega)*obj.Parent.SampleRate;
            alpha(end+1,:)=nan;
        end

        function power=get.Power(obj)
            try                
                au=obj.AngleUnits;
                obj.AngleUnits="rad";
                power=obj.AngularVelocity.*obj.Moment;
                if obj.MomentUnits=="Nmm"
                    power=power/1000;
                end
                %obj.AngleUnits=au;
            catch
            power=[];
            end
        end

        function updateC3D(obj,H,mode)
        %UPDATEC3D Update the btk Handle 
        %UPDATEC3D updates the data in the btk handle specified by H with the desired mode
        %"overwrite" overwrites existing data and adds new data
        %"append" ignores existing data and puts all the new jointangle
        %"newfile" resets all the angle data in the acquisition
            arguments
                obj Joint
                H double
                mode string {mustBeMember(mode,["overwrite","newfile","append"])}="overwrite"
            end
        for quantity=["Angle","Moment","Power"]
            try
            if isequal(quantity,"Angle")
                descr=obj.AngleSequence;
                [~,ord]=sort(double(char(obj.AngleSequence)));
                data=obj.CardanAngles(:,ord);
            elseif isequal(quantity,"Moment")
			    descr=sprintf("expressed in %s CS",obj.Parent.Label);
                data=obj.Moment;
            elseif isequal(quantity,"Power")
			    descr="";
                data=obj.Power;
            end
            if isempty(data)
               continue %nothing to write
            end
			label=char(obj.Label+obj.LabelSuffix+quantity);
            unit=char(obj.(quantity+"Units"));
            cunits=btkGetPointsUnit(H,char(quantity));
            if isequal(cunits,unit)||not(isequal(mode,"append"))
            else
                fhandle=sprintf("btkGet%ss",quantity);
                m=feval(fhandle,H);
                if isempty(fieldnames(m))
                    btkSetPointsUnit(H,char(quantity),cunits)
                else
                error("Existing  Units in C3D (%s) mismatch current Joint units (%s)!",cunits,unit);
                end
            end

        switch mode
            case "overwrite"  
                try
                btkSetPoint(H,label,data);
                catch
                btkAppendPoint(H,char(lower(quantity)),label,data,data*0,char(descr));
                end
        case "append"
                btkAppendPoint(H,char(lower(quantity)),label,data,data(:,1)*0,char(descr));
        end
            catch ME
                warning("Error when writing %sJoint%s: %s !",obj.Label,quantity,ME.message);
            end
        end %quantity
        obj.Acquisition.readC3DMetaData;
        end %updateC3D
    
    %% TIME NORMALIZATION
    function [c,cmedio,cstd]=timeNorm(obj,npoints,base)
    %[c,cmedio,cstd]=timeNorm(obj,npoints,base) returns the angle, angular
    %velocity, power and moment in the joint, normalized cycles (npoints) for the stance or the stride
        arguments
            obj Joint
            npoints=101;
            base string {mustBeMember(base,["stance","stride"])}="stride"
        end
        ev=exportEvents(obj.Acquisition.Events,'point',false);
        group=obj.Group;
        % if matches(group,["Left","Right"])
        % else
        % %group=["Left", "Right"];
        % group = "Left";
        % end
        S=[]; E=[];
        for g=group
            S=[S ev.(g).Foot_Strike];
            if base=="stride"
                E=[E ev.(g).Foot_Strike(2:end)];
                S(end)=[];
            else
                E=[E ev.(g).Foot_Off];
            end
        end 
        
        for i=(length(S)):-1:1
            try
            c.angle(:,:,i)=time2cycle([],obj.CardanAngles(S(i):E(i),:),npoints);
            c.omega(:,:,i)=time2cycle([],obj.AngularVelocity(S(i):E(i),:),npoints);
            c.moment(:,:,i)=time2cycle([],obj.Moment(S(i):E(i),:),npoints);
            c.power(:,:,i)=time2cycle([],obj.Power(S(i):E(i),:),npoints);
            catch
            end
        end
        cmedio=structfun(@(x) mean(x,3,'omitnan'),c,"UniformOutput",false);
        cstd=structfun(@(x) std(x,1,3,'omitnan'),c,"UniformOutput",false);
    end

    end %methods
end %Joint