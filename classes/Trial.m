 classdef Trial<handle
    %Trial creates a container/manipulator for a C3D trial
    %   This class manages the contents of a biomechanical acquisition
    %   captured with marker based systems and allows the manipulation of
    %   its elements following standard procedures in motion analysis

    properties
        Metadata %struct containing C3D Metadata
        Points Point %points in the acquisisition 
        ForcePlatform forcePlatformType2 %force platforms in the acquisition
        Analogs Analog %Analog channels contained in the acquisition
        Events Event %Events stored in the acquisition
        Segments %Segments defined in the biomechanical model
        Joints %Joints defined in the biomechanical model
        Scalars Scalar %Scalars contained in the acquisition
        c3dfile char %Filename of the acquisition
        ConfigFile struct %Configuration file
        configfile (1,1) string %filename of the configuration
        XMLatt (1,1) string ="_att"; %attribute to be used in the configuration
        useJSON (1,1) logical =false;% selector to use JSON configuration file
        NFrames %Number of frames in the C3D
        NSamples%Number of analog samples in the C3D
     end
     properties (Transient)
        C3DHandle%Memory address of the loaded C3D file
     end

     methods
        function obj = Trial(c3dfile,config)
            %Trial Construct an instance of Trial
            %   this methods builds an instance of the class Trial. The
            %   first input must be a c3dfile. The second optional input is
            %   a configuration structure containing info to parse the raw
            %   information (e.g., markers) to structured groups such as
            %   Segments and Joints
            if nargin==0
                return
            end
            obj.c3dfile=c3dfile;
            H=btkReadAcquisition(obj.c3dfile);  
            obj.C3DHandle=H;
            obj=obj.readC3DMetaData();
            obj.Points=Point(H);
            [obj.Points.Parent]=deal(obj);
            obj.ForcePlatform=forcePlatformType2(H);
            [obj.ForcePlatform.Parent]=deal(obj);
            obj.Analogs=Analog(H);
            [obj.Analogs.Parent]=deal(obj);
            obj.Scalars=Scalar(H);
            [obj.Scalars.Parent]=deal(obj);
            obj.Events=Event(H);
            obj.Events.Parent=obj;
            obj.NFrames=btkGetPointFrameNumber(H);
            obj.NSamples=btkGetAnalogFrameNumber(H);
            
            
            if nargin==2
                if isstruct(config)
                elseif isfile(config)
                    obj.configfile=config;
                    try
                    config=readstruct(config,'FileType','xml','AttributeSuffix',obj.XMLatt);
                    catch
                    config=readstruct(config,'FileType','json');
                    obj.XMLatt="";
                    obj.useJSON=true;
                    end  
                else
                    error('2-nd argument must be a valid configuration file or structure!');
                end
            obj.applyConfiguration(config);
            end
        end
        
        function obj=applyConfiguration(obj,config)
        % applyConfiguration applies a configuraiion file to a Trial
        %   obj=applyConfiguration(obj,config) applies a configuration
        %   structure to the trial setting point properties and building the
        %   skeleton, returning the updated trials
            for i=1:length(obj)
                obj(i).ConfigFile=config;
                for m=obj(i).ConfigFile.MarkerSet.Marker
                p=matches([obj(i).Points.Label],m.("label"+obj(i).XMLatt));
                if sum(p)==1
                obj(i).Points(p).Type=m.("type"+obj(i).XMLatt);
                obj(i).Points(p).Segment=m.("segment"+obj(i).XMLatt);
                obj(i).Points(p).Cluster=m.("cluster"+obj(i).XMLatt);
                obj(i).Points(p).Group=m.("group"+obj(i).XMLatt);
                end
                end   
                obj(i).buildSkeleton;
            end
        end
        %% INERTIAL PROPERTIES
        %set of function operating on the subject parameters and on the
        %inertial properties of segments

        function obj=setSubjectAntropometry(obj,mass,gender)
            %setSubjectAntropometry Set the subject antropometry
            %   obj=setSubjectAntropometry(obj,mass,mf) sets the subject
            %   antropometric characteristics (mass, gender) to specified
            %   user inputs. If antropometric metadata is already present in
            %   the PROCESSING group (e.g., if already defined in Vicon
            %   Nexus), then this data is used.
            arguments
                obj Trial
                mass=mean(vecnorm(sum(cat(3,obj.ForcePlatform.GRF),3),2,2))/9.81;
                gender char {mustBeMember(gender,["Male","Female"])}='Female';
            end

            if isfield(obj.Metadata,'PROCESSING')
               if isfield(obj.Metadata.PROCESSING,'mass')
                   mass=obj.Metadata.PROCESSING.mass;
               end
            end
            obj.Metadata.ANTROPOMETRY.gender=gender;
            obj.Metadata.ANTROPOMETRY.mass=mass;
            obj.Metadata.ANTROPOMETRY.units='kg';
            obj.setC3DMetaData;
        end
        function obj=scaleInertialProp(obj)
            %SCALEINERTIALPROP Scale inertial properties based on De Leva
            %   obj=SCALEINERTIALPROP(obj) use antropometric metadata
            %   contained in obj and scales the De Leva-Tsairkovsky table
            %   accordingly. Segments COM, mass, and I are added to the
            %   biomechanical model in obj.
            %   
            % See Also SETSUBJECTANTROPOMETRY

        obj.Metadata.SEGMENTS_MASS.units="kg";
        obj.Metadata.SEGMENTS_ICM.units="kgm^2";

        if not(isfield(obj.Metadata.ANTROPOMETRY,'gender'))
            obj.Metadata.ANTROPOMETRY.gender='Male';        
        elseif isempty(obj.Metadata.ANTROPOMETRY.gender)
            obj.Metadata.ANTROPOMETRY.gender='Male';
        end
        
        try
            mainpath=mfilename('fullpath');
            mainpath=fileparts(mainpath);
            dltable=fullfile(mainpath,'..','..', 'Models','Olympia','XML','DeLeva.inertialparam');
            if isfile(dltable)
            else
                [f,p]=uigetfile('*.inertialparam','Select De Leva Table');
                dltable=fullfile(p,f);
            end
        DLtable=readstruct(dltable,'FileType','xml','AttributeSuffix',obj.XMLatt);
        DLtable=DLtable.(string(obj.Metadata.ANTROPOMETRY.gender));
        catch
        error('Unable to open Antropometric Table')
        end
        
        for group=string(fieldnames(obj.Segments))'
            for seg=string(fieldnames(obj.Segments.(group)))'
            
                segment=obj.Segments.(group).(seg); %current segment
                alias=obj.ConfigFile.KinematicModel.Segment([obj.ConfigFile.KinematicModel.Segment.("label"+obj.XMLatt)]==seg).("inertialalias"+obj.XMLatt);
                if alias==""||ismissing(alias)
                   continue
                end
                try

                % select the row for the current segment
                tmpDL=DLtable.(alias);

                % extract endpoints, calculate CM, and add to output
                try
                P=segment.EndPoints(1).Coordinates;
                D=segment.EndPoints(2).Coordinates;
                catch
                error(['Required EndPoints not found in trial for ' segment.Label]);    
                end
                % COM definition and append to trial
                com=P+tmpDL.CM.("value"+obj.XMLatt)/100*(D-P);
                try
                comlabel=obj.ConfigFile.KinematicModel.Segment([obj.ConfigFile.KinematicModel.Segment.("label"+obj.XMLatt)]==seg).("comlabel"+obj.XMLatt);
                catch
                comlabel=segment.Label+"COM";
                end
                if comlabel=="" || ismissing(comlabel)
                    comlabel=segment.Label+"COM";
                end
                obj.Points(end+1)=obj.Points.appendPoint(comlabel,com,segment.Label,'V',0,group);
                % verify that the COM is in the markerset
                plist=[obj.ConfigFile.MarkerSet.Marker.("label"+obj.XMLatt)];
                if sum(matches(plist,comlabel))<1 %and in case add it
                obj.ConfigFile.MarkerSet.Marker(end+1).("label"+obj.XMLatt)=comlabel;
                obj.ConfigFile.MarkerSet.Marker(end).("segment"+obj.XMLatt)=segment.Label;
                obj.ConfigFile.MarkerSet.Marker(end).("cluster"+obj.XMLatt)=0;
                end
                obj.Segments.(group).(seg).Points(end+1)=obj.Points(matches([obj.Points.Label],comlabel));
                obj.Segments.(group).(seg).COM=obj.Points(matches([obj.Points.Label],comlabel));
                % Segment Mass
                mi=tmpDL.Mass.("value"+obj.XMLatt)/100*obj.Metadata.ANTROPOMETRY.mass;
                obj.Segments.(group).(seg).Mass=mi;
                % Segment ICM
                stl2xyz=[1 3 2];
                %stl2xyz=[1 2 3];
                lsegm=vecnorm(D-P,2,2);
                
                Ii(stl2xyz)=([tmpDL.gyration_radii.("rSagg"+obj.XMLatt) tmpDL.gyration_radii.("rTrans"+obj.XMLatt) tmpDL.gyration_radii.("rLong"+obj.XMLatt)]/100*mean(lsegm)).^2*mi;
                if obj.Metadata.POINT.UNITS=="mm"
                   Ii=Ii*1e-6;
                elseif obj.Metadata.POINT.UNITS=="m"
                else
                   error("Length Units must be either m or mm");
                end
                obj.Segments.(group).(seg).Icm=diag(Ii);

                obj.Segments.(group).(seg).InertialUnits=["kg","kgm^2"];
                obj.Metadata.SEGMENTS_MASS.(seg)=mi;
                obj.Metadata.SEGMENTS_ICM.(seg)=diag(Ii);
                catch ME
                    warning("Error when setting inertial properties for %s: %s",segment.Label,ME.message);
                end
            end
        end

        obj.setC3DMetaData;
        end
        
        function obj=moveInertialProperties(obj,src)
        % MOVEINERTIALPROPERTIES move segment inertial properties to a
        % different trial
        %   obj=moveInertialProperties(obj,src) assign the Inertial
        %   properties from a trial (e.g., static) to the current trial
        % 
        % See Also SCALEINERTIALPROP
        i=1;
        obj.Metadata.ANTROPOMETRY=src.Metadata.ANTROPOMETRY;
        obj.Metadata.SEGMENTS_MASS.units=src.Metadata.SEGMENTS_MASS.units;
        obj.Metadata.SEGMENTS_ICM.units=src.Metadata.SEGMENTS_ICM.units;
        for group=string(fieldnames(src.Segments))'
            for seg=string(fieldnames(src.Segments.(group)))'
            try
            obj.Segments.(group).(seg).COM=obj.Points(matches([obj.Points.Label],src.Segments.(group).(seg).COM.Label));
            obj.Segments.(group).(seg).Mass=src.Metadata.SEGMENTS_MASS.(seg);
            obj.Metadata.SEGMENTS_MASS.(seg)=src.Metadata.SEGMENTS_MASS.(seg);
            obj.Segments.(group).(seg).Icm=src.Metadata.SEGMENTS_ICM.(seg);
            obj.Metadata.SEGMENTS_ICM.(seg)=src.Metadata.SEGMENTS_ICM.(seg);
            m(i)=obj.Segments.(group).(seg).Mass;    
            com(:,:,i)=obj.Segments.(group).(seg).COM.Coordinates*m(i);
            i=i+1;
            catch
            end %try
            end %seg
        end %side
        if exist('com','var')
        COM=sum(com,3)/sum(m);
        obj.Points(end+1)=obj.Points.appendPoint("COM",COM);
        obj.setC3DMetaData;
        end
        end
        
        function obj=readC3DInertialProperties(obj)
                 for group=string(fieldnames(obj.Segments))'
                    for seg=string(fieldnames(obj.Segments.(group)))'
                        segment=obj.Segments.(group).(seg); %current segment
                        try
                        segment.Mass=obj.Metadata.SEGMENTS_MASS.(seg);
                        segment.Icm=obj.Metadata.SEGMENTS_ICM.(seg);
                        comlabel=obj.ConfigFile.KinematicModel.Segment(matches([obj.ConfigFile.KinematicModel.Segment.("label"+obj.XMLatt)],seg)).("comlabel"+obj.XMLatt);
                        segment.COM=obj.Points(matches([obj.Points.Label],comlabel));
                        catch
                        warning('No Inertial data stored in C3D for %s',seg);
                        end
                    end
                 end
        end
        %% POINT RECONSTRUCTION AND OPTIMIZATION
        function obj=appendVirtualMarker(obj,method,label,sourcemarkers,segment,type,wand)
        %APPENDVIRTUALMARKER Add virtual/calibrated markers to the
        %acquisition
        %appendVirtualMarker(obj,method,label,sourcemarkers,segment,type) 
        % adds a virtual marker to the acquisition, calculated according to several methods
        % Inputs:
        %   obj             - Trial to operate
        %   method          - Method of reconstruction
        %   label           - Label to assign to the new marker
        %   sourcemarkers   - Markers to be used as a base for the virtual
        %                     marker
        %   segment         - Segment to which the marker will be assigned
        %   type            - Type of marker
        %   wand            - Characteristics of the wand used for calibrated markers
        %
        % See Also SVDRECONSTRUCT
            sourcemarkers=strtrim(split(sourcemarkers,','));
            try
            switch method
                case 'centroid'
                    for i=length(sourcemarkers):-1:1
                        pt(:,:,i)=obj.Points.PointStruct.(sourcemarkers(i));
                    end
                    pt=mean(pt,3);
                    obj.Points(end+1)=obj.Points.appendPoint(label,pt,segment,type);
                case 'Bell'
                    pelvis=obj.ConfigFile.KinematicModel.Segment;
                    pelvis=pelvis(matches([pelvis.("label"+obj.XMLatt)],segment));
                    [P]=buildCS(obj.Points.PointStruct,pelvis.CoordSys);
                    asis=pelvis.CoordSys.origin;
                    PW=mean(vecnorm(obj.Points.PointStruct.(asis(1))-obj.Points.PointStruct.(asis(2)),2,2),'omitnan');
                    if contains(label,"L")
                        pm=-1;
                    elseif contains(label,"R")
                        pm=1;
                    end
                    HJC = squeeze(pagemtimes(P,[[-0.19, -0.30, pm*0.36]*PW, 1]'))';            
                    obj.Points(end+1)=obj.Points.appendPoint(label,HJC(:,1:3),segment,type);
                case 'wand' %obsolete
                    path=fileparts(obj.c3dfile);
                    wtr=Trial(char(fullfile(path,sourcemarkers)));
                    wtr.ConfigFile=obj.ConfigFile;
                    obj=obj.wandReconstruct(wtr,label,wand);
                    obj.Points(matches([obj.Points.Label],label)).Segment=segment;
                    obj.Points(matches([obj.Points.Label],label)).Type=type;
                otherwise 
                    error('Unrecognized Method: %s', method);
            end
            catch ME
                warning(ME.message);
            end
        end

        function obj=reconstructVirtualMarkers(obj)
        % RECONSTRUCTVIRTUALMARKERS Reconstruct virtual/calibrated markers
        % defined in the configuration file
        % obj=reconstructVirtualMarkers(obj) applies the
        % APPENDVIRTUALMARKER function to all the virtual and calibrated
        % markers defined in the configuration file
        % See Also APPENDVIRTUALMARKER
            try
            mkr=reshape(obj.ConfigFile.MarkerSet.VirtualMarkerDef,1,[]);
            for cm=mkr
                m=obj.ConfigFile.MarkerSet.Marker(matches([obj.ConfigFile.MarkerSet.Marker.("label"+obj.XMLatt)],cm.("label"+obj.XMLatt)));    
                if length(m)~=1
                   warning("Virtual marker with label: %s is not part of the defined markerset and won't be created!",cm.("label"+obj.XMLatt));
                else
                    try
                    obj=obj.appendVirtualMarker(cm.("method"+obj.XMLatt),cm.("label"+obj.XMLatt),cm.("source"+obj.XMLatt),m.("segment"+obj.XMLatt),m.("type"+obj.XMLatt));
                    catch ME
                        warning("Error occurred during reconstruction of %s",cm.("label"+obj.XMLatt));
                        warning(ME);
                    
                    
                    end
                end
            end
            catch
                warning('Virtual Markers are absent or not properly defined!')
            end
        end
        function obj = wandReconstruct(obj,wand_trial,wand,label)
            % WANDRECONSTRUCT Reconstruct a wand-calibrated marker
            % obj=wandReconstruct(obj,wand_trial,point,wand) Reconstruct a
            % point from a wand trial, if the wand has multiple frames, a
            % mean function is applied.
            % Inputs:
            % obj           - Currently open trial to which add the marker
            % wand_trial    - Trial containing wand coordinates used to
            %                 calibrate the marker
            % wand          - Struct containing labels of the wand tip and
            %                 tail points, and the distance of the tail marker 
            %                 to the phisical tip of the wand 
            % label         - Label to be assigned to the marker
            arguments 
                obj Trial
                wand_trial Trial
                wand struct
                label char
            end
            if ~isequal(obj.Metadata.POINT.UNITS,wand_trial.Metadata.POINT.UNITS)
                error("Both Wand and Current Trial must have the same length units!");
            end
            if wand.TailMarker==wand.TipMarker
            newpt=wand_trial.Points.PointStruct.(wand.TailMarker);
            else
            newpt=wand_trial.Points.PointStruct.(wand.TailMarker)+wand.Tail2RealTip*normalize(wand_trial.Points.PointStruct.(wand.TipMarker)-wand_trial.Points.PointStruct.(wand.TailMarker),2,'norm');
            end
            if size(newpt,1)>1
            newpt=mean(newpt,'omitnan');
            wand_trial=mean(wand_trial);
            end
            pts=obj.ConfigFile.MarkerSet.Marker; %points in the trial
            cl=logical([obj.ConfigFile.MarkerSet.Marker.("cluster"+obj.XMLatt)]); % points to be used as cluster
            seg=pts(matches([pts.("label"+obj.XMLatt)],label)).("segment"+obj.XMLatt);  %segment of the point to be reconstructed
            group=pts(matches([pts.("label"+obj.XMLatt)],label)).("group"+obj.XMLatt);
            pts=pts(cl);
            plist=[pts(matches([pts.("segment"+obj.XMLatt)],seg)).("label"+obj.XMLatt)];  %points belonging to the segment
            plist=plist(matches(plist,[wand_trial.Points.Label])&matches(plist,[obj.Points.Label])); %points available in both trials
            p_src=cat(1,wand_trial.Points(matches([wand_trial.Points.Label],plist)).sort.Coordinates);
            p_tgt=permute(cat(3,obj.Points(matches([obj.Points.Label],plist)).sort.Coordinates),[3 2 1]);

            [~,~,T]=svdSolidification(p_src,p_tgt);
            newpt=pagemtimes(T,[newpt 1]');
            newpt=permute(newpt,[3 1 2]);
            obj.Points(end+1)=obj.Points.appendPoint(label,newpt(:,1:3),seg,'C',false,group);
        end

        
       
    function obj=segmentSolidification(obj,stat,existingpoints)
        %SEGMENTSOLIDIFICATION Apply SVD to segments to reconstruct and
        %optimize marker positions
        %
        %   obj=segmentSolidification(obj,stat) use SVD to reconstruct
        %   anatomical segments in a dynamic trial starting from a static
        %   trial
        % Inputs:
        % obj   - Trial which will be manipulated
        % stat  - Trial used as source for markers
        % existingpoints    - 'replace'|'ignore' select what to do with
        %                     points which already exist
        % See Also SVDRECONSTRUCT
        arguments
            obj Trial
            stat Trial
            existingpoints string {mustBeMember(existingpoints,["replace","ignore"])}="replace";
        end
        
        if ~isequal(obj.Metadata.POINT.UNITS,stat.Metadata.POINT.UNITS)
           error("Both Source and Current Trial must have the same length units!");
        end
        if stat.Metadata.POINT.FRAMES~=1
            error('Static Trial must have a single frame!')
        end
        segments=unique([obj.Points.Segment]);
        segments(segments=="")=[];
        for s=segments
            try
            p_dyn=obj.Points(matches([obj.Points.Segment],s));%&[obj.Points.Cluster]);% points in the dynamic trial in the segment and used in cluster definition
            p_dyn=p_dyn.sort; %sorted in alphabetical order
            p_stat=stat.Points(matches([stat.Points.Segment],s)); %points in the static trial belonging to the segment
            p_stat=p_stat.sort;%sorted in alphabetical order                
            p_statdyn=p_stat(matches([p_stat.Label],[p_dyn.Label])); % points matching in both trials
            
            p_statsvd=p_statdyn([p_statdyn.Cluster]);% (base for SVD)
            p_dynsvd=p_dyn(matches([p_dyn.Label],[p_statsvd.Label]));% (base for SVD)
            
            p_dyn=p_dyn(matches([p_dyn.Label],[p_statdyn.Label])); % points ONLY in dynamic trial (will be optimized)
            p_stat=p_stat(not(matches([p_stat.Label],[p_dyn.Label]))); % points ONLY in static trial (will be reconstructed)
                         
            % if length(p_dyn)==3
            %    warning('%s: At Least 4 points common to the two trials are needed to use svd optimization!',s);
            %    existingpoints="ignore";
            % else
            if length(p_dyn)<3
               error("%s: Unable to reconstruct local points in segments with less than 3 points!",s)
            end
            cl_dyn=permute(cat(3,p_dynsvd.Coordinates),[3 2 1]); %moving cluster page-framed (can have multiple frames)
            cl_statdyn=cat(1,p_statsvd.Coordinates); %source cluster (must have a single frame)
            cl_stat=cat(1,p_stat.Coordinates); %cluster of points to be added in the trial

            [~,~,T]=svdSolidification(cl_statdyn,cl_dyn); %compute T matrix moving cl_statdyn to cl_dyn

            % add missing points from static trial
            for i=1:length(p_stat)
               p=p_stat(i);
               XYZ=permute((pagemtimes(T,[cl_stat(i,:) 1]')),[3 1 2]);
               obj.Points(end+1)=obj.Points.appendPoint(p.Label,XYZ(:,1:3),p.Segment,p.Type,false,p.Group);
            end

            %if specified, replace collected points with svd version
            if isequal(existingpoints,"replace") %this run only if existing points need to be optimized
                cl_statdyn=cat(1,p_statdyn.Coordinates);
                for i=1:length(p_dyn)
                   p=p_dyn(i);
                   XYZ=permute((pagemtimes(T,[cl_statdyn(i,:) 1]')),[3 1 2]);
                   p.XData=XYZ(:,1);
                   p.YData=XYZ(:,2);
                   p.ZData=XYZ(:,3);
                   %obj.Points(end+1)=obj.Points.appendPoint(p.Label,XYZ(:,1:3),p.Segment,p.Type,false,p.Group);
               end
            end
            catch ME
                 warning("%s : %s",s, ME.message);
            end %try catch
        end %segment loop
        end %fun


        %% KINEMATIC MODEL AND ANGLES

        function obj=buildSkeleton(obj)
            %BUILDSKELETON Creates Segment and Joints from the
            %configuration file
            %   obj=BUILDSKELETON(obj)creates the segments and joints from the information inserted
            %   in the configuration file under KinematicModel. Each segment
            %   and joint is parsed into a structure in different fields
            %   based on their group attribute, and then into different
            %   fields based on their label attribute.
            
            %%-SEGMENTS CREATION
            for cseg=obj.ConfigFile.KinematicModel.Segment
                     pt=split(cseg.("points"+obj.XMLatt),',')';
                     pt=obj.Points(matches([obj.Points.Label],pt));
                try % to build the segment

                    %first info from config
                    seg=Segment();
                    seg.Label=cseg.("label"+obj.XMLatt);
                    seg.Points=pt;
                    seg.Parent = obj;
                    seg.SampleRate=obj.Points.SampleRate;

                    try
                    % insert Coordsys definition
                    seg.CoordsysDef=cseg.CoordSys;
                    catch
                    warning('Unable to generate Coordinate System for: %s',cseg.("label"+obj.XMLatt))
                    end
                    % finally place the segment in context
                    try
                        ep=split(cseg.("endpoints"+obj.XMLatt),',');
                        seg.EndPoints=obj.Points(matches([obj.Points.Label],ep));
                    catch
                    end
                    % Group selector
                    if ismissing(cseg.("group"+obj.XMLatt))
                        group="Undefined";
                    elseif cseg.("group"+obj.XMLatt)==""
                        group="Undefined";
                    else
                        group=cseg.("group"+obj.XMLatt);
                    end
                    seg.Label=cseg.("label"+obj.XMLatt);
                    seg.Group=group;
                    obj.Segments.(group).(seg.Label)=seg;
                catch ME
                    warning(ME.message)
                end
            end %segment loop
        
            %%--JOINTS
            for joint=obj.ConfigFile.KinematicModel.Joint %loop trough each joint
                try
                    % group parser block
                    if ismissing(joint.("group"+obj.XMLatt))
                        group="Undefined";
                    elseif joint.("group"+obj.XMLatt)==""
                        group="Undefined";
                    else
                        group=joint.("group"+obj.XMLatt);
                    end
  
                % selection of the proximal segment
                prox=joint.("prox"+obj.XMLatt);
                cseg=obj.ConfigFile.KinematicModel.Segment(matches([obj.ConfigFile.KinematicModel.Segment.("label"+obj.XMLatt)],prox));
                    if ismissing(cseg.("group"+obj.XMLatt))
                        tmp="Undefined";
                    else
                        tmp=cseg.("group"+obj.XMLatt);
                    end
                proxSegm=obj.Segments.(tmp).(cseg.("label"+obj.XMLatt));
                %selection of distal segment
                dist=joint.("dist"+obj.XMLatt);
                cseg=obj.ConfigFile.KinematicModel.Segment(matches([obj.ConfigFile.KinematicModel.Segment.("label"+obj.XMLatt)],dist));
                    if ismissing(cseg.("group"+obj.XMLatt))
                        tmp="Undefined";
                    else
                        tmp=cseg.("group"+obj.XMLatt);
                    end
                distSegm=obj.Segments.(tmp).(cseg.("label"+obj.XMLatt));
                try
                JC=obj.Points(matches([obj.Points.Label],joint.("JC"+obj.XMLatt)));
                if isempty(JC)
                   JC=Point;
                end
                catch
                   JC=Point;
                end
                units="N"+proxSegm.Points(1).Units;
                obj.Joints.(group).(joint.("label"+obj.XMLatt))=Joint(joint.("label"+obj.XMLatt),group,proxSegm,distSegm,joint.("type"+obj.XMLatt),joint.("angleseq"+obj.XMLatt),"deg",JC,units);
                obj.Joints.(group).(joint.("label"+obj.XMLatt)).Acquisition=obj;
                catch ME
                    warning(ME.message);
                end
                
            end %joint loop
        end

       
        %% CROSS-OBJECT TRIAL ANALYSIS
        function [obj, grf] = GRFAnalysis(obj,normalize)
         %GRFANALYSIS Computes impulses and GRF parameters
         %  [obj, grf] = GRFAnalysis(obj,normalize) computes impulses and GRF
         %  parameters from the forcePlatform data stored in the acquisition,
         %  assigning it into the Metadata.
         %  Inputs:
         %  obj         - Trial
         %  normalize   - Select what to use for mass normalization
         %  Outputs:
         %  obj - Updated Trial
         %  grf - structure containing the calculated GRF data

            arguments
            obj (1,1) Trial
            normalize string {mustBeMember(normalize,["no","bodymass","bodyweight"])}="bodyweight";
            end
        
        forceplat=obj.ForcePlatform;
        if numel(forceplat)==0
           warning("Trial Contains no Force Platforms!, GRF analysis can't be performed!");
           return
        end   

        
        if btkGetEventNumber(obj.C3DHandle)==0 %needed because metadata isn't updated (btk bug?)
           warning("Trial Contains no events!, GRF analysis can't be performed!");
           return
        else
           events=obj.Events;
        end

 

        try %import mass from trial and adjust units
        mass=obj.Metadata.ANTROPOMETRY.mass;
        if obj.Metadata.ANTROPOMETRY.units=="t"
           mass=mass*1000;
        end
        mustBePositive(mass);
        switch normalize
            case "no"
                mass=1;
                unit="N%s";
            case "bodymass"
                %mass=mass;
                unit="N%s/kg";
            case "bodyweight"
                mass=mass*9.81;
                unit="N%s/N";
        end
        catch
            mass=1;
            unit="N%s";
        end

   
        GRF=sum(cat(3,forceplat.GRF),3)/mass;
        R=forceplat.align2ISB();

        GRF=GRF*R'; %rotate the GRF to the new directions
        AP=1; % now the first direction has the most variability (i.e., is the running direction)
        %ML=2;
        V=2;  % in a typical flat laboratory, the vertical should change least (i.e., is vertical)
        
        
        f_force=forceplat.SampleRate;
        events=events.selectFootContacts;
        events=exportEvents(events,'analog',false);
        events=rmfield(events,'units');
        
        for ctx=string(fieldnames(events))' %run the analysis for each context of events
            if isempty(fieldnames(events.(ctx)))
                continue
            end
            try
            FS=events.(ctx).Foot_Strike;
            FO=events.(ctx).Foot_Off;
            for j=length(FS):-1:1
            fc=FS(j);
            fo=FO(j);
            
            %%-MAX, MIN, MEAN
            % anteroposterior
            grf.GRF_PARAM.Units=sprintf(unit,"");
            grf.GRF_PARAM.(ctx+"AnteroPosteriorMax")(j)=mean(GRF(fc:fo,AP));
            grf.GRF_PARAM.(ctx+"AnteroPosteriorMin")(j)=min(GRF(fc:fo,AP));
            grf.GRF_PARAM.(ctx+"AnteroPosteriorMaxMean")(j)=max(GRF(fc:fo,AP));
            % vertical
            grf.GRF_PARAM.(ctx+"VerticalMax")(j)=max(GRF(fc:fo,V));
            grf.GRF_PARAM.(ctx+"VerticalMean")(j)=mean(GRF(fc:fo,V));
            %%-IMPULSES
            grf.GRF_IMPULSES.Units=sprintf(unit,"s");
            grf.GRF_IMPULSES.(ctx+"Vertical")(j)=trapz(GRF(fc:fo,V)-1)/f_force; %(Fy-BW)/BW= Fy/BW-1
            %%-Braking/Propulsive impulses X for each step of the same limb
            GRF_AP_e=GRF(fc:fo,AP);
            zerocross=[ne(diff(sign(GRF_AP_e)),0);false]; %find zero crossing
            Ic=cumtrapz(GRF_AP_e)/f_force;
            Ic=[Ic(zerocross)', Ic(end)];
            In=Ic-[0,Ic(1:end-1)];
            
            grf.GRF_IMPULSES.(ctx+"AnteroPosteriorStartPropulsion")(j)=sum(In)-(sum(In(In<0))+max(In));
            grf.GRF_IMPULSES.(ctx+"AnteroPosteriorBraking")(j)=sum(In(In<0));
            grf.GRF_IMPULSES.(ctx+"AnteroPosteriorPropulsive")(j)=max(In);
            grf.GRF_IMPULSES.(ctx+"AnteroPosteriorNet")(j)=sum(In);    

           
            
            end
            catch ME
                warning(ME.message)
            end
        end
            if exist('grf','var')
                if isfield(grf,"GRF_PARAM")
            obj.Metadata.GRF_PARAM=grf.GRF_PARAM;
                end
                if isfield(grf,"GRF_IMPULSES")
            obj.Metadata.GRF_IMPULSES=grf.GRF_IMPULSES;
                end
            obj.setC3DMetaData;
            else
               warning("GRF analysis not successfull, check events and forceplate data!")
               grf=struct();
            end
        end

        function SPATIOTEMP=spatioTemporalAnalysis(obj)
         %SPATIOTEMPORALANALYSIS Computes spatial-temporal parameters
            arguments
            obj (1,1) Trial
            end

        forceplat=obj.ForcePlatform;
        if numel(forceplat)==0
           warning("Trial Contains no Force Platforms!, GRF analysis can't be performed!");
           return
        end   

        
        if getEventCount(obj.Events)==0 %needed because metadata isn't updated (btk bug?)
           warning("Trial Contains no events!, GRF analysis can't be performed!");
           return
        else
           events=obj.Events;
        end

        fp=forceplat.combineFP;
   
        %GRF=sum(cat(3,forceplat.GRF),3)/mass;
        COP=fp.COP;
        R=forceplat.align2ISB();

        %GRF=GRF*R'; %rotate the GRF to the new directions
        COP=COP*R';
        AP=1; % now the first direction has the most variability (i.e., is the running direction)
        %ML=2;
        V=2;  % in a typical flat laboratory, the vertical should change least (i.e., is vertical)
        events=events.selectFootContacts;
    ev_force=events.exportEvents("analog");
    ev_point=events.exportEvents("point");
    ev_time=events.exportEvents("seconds");
side=["Left", "Right"];
ss=["L","R"];
for s=1:2
    lr=side(s);                  %considered side
    cl=side(not(side==side(s))); %controlateral
    %stance time
    SPATIOTEMP.(lr+"ContactTime")=ev_time.(lr).Foot_Off-ev_time.(lr).Foot_Strike;
    
    %flight time
    FO_lr=ev_time.(lr).Foot_Off;
    FS_cl=ev_time.(cl).Foot_Strike;
    FS_cl(FS_cl<FO_lr(1))=[];
    ns=min(numel(FO_lr),numel(FS_cl));
    SPATIOTEMP.(lr+"FlightTime")=FS_cl(1:ns)-FO_lr(1:ns);
    % stance length 
    P=obj.Points.PointStruct.(ss(s)+"GT")*R';
    P=P(:,AP);

    SPATIOTEMP.(lr+"ContactLengthGT")=P(ev_point.(lr).Foot_Off)-P(ev_point.(lr).Foot_Strike);
    

    % stride length
    FS_lr=ev_force.(lr).Foot_Strike;
    FS_cl=ev_force.(cl).Foot_Strike;
    FS_cl(FS_cl<FS_lr(1))=[];
    ns=min(numel(FS_lr),numel(FS_cl));
    COP_lr=COP(FS_lr(1:ns),AP);
    COP_cl=COP(FS_cl(1:ns),AP);

    SPATIOTEMP.(lr+"StrideLength")=COP_cl-COP_lr;
end
obj.Metadata.SPATIOTEMP=SPATIOTEMP;
obj.setC3DMetaData;

end


        function obj=inverseDynamics(obj,useID,g)
        % INVERSEDYNAMICS Calculate Joint force and moment
        %   obj=inverseDynamics(obj,useID,g)
        % resolves each kinematic chain described in the configuration and
        % returns the updated trial
            arguments
                obj Trial
                useID logical=true
                g=[0 0 -9.806]; %default is Z up
            end

            u=obj.Metadata.POINT.UNITS;
            
            obj.setUnits("m"); %inverse dynamics work in m, rad, kg
            for kc=obj.ConfigFile.KinematicModel.KinematicChain
                try
            s=kc.("group"+obj.XMLatt);
            joints=strrep(split(kc.("joints"+obj.XMLatt),',')," ","");
            endpoint=kc.("endbody"+obj.XMLatt);
                % select the endpoint type
                if endpoint=="ForcePlatform"
                    force=obj.ForcePlatform;
                    if length(force)>1
                    force=force.combineFP;
                    end
                    force=force.resample(obj.Points(1).SampleRate);
                    F=force.GRF;
                    M=force.GRM*0;
                    COP=force.COP;
                    ev=obj.Events.exportEvents('point',false);
                    FC=ev.(kc.("group"+obj.XMLatt)).Foot_Strike;
                    FO=ev.(kc.("group"+obj.XMLatt)).Foot_Off;
                    ctc=false(size(F));
                    for i=1:length(FC)
                        ctc(FC(i):FO(i),:)=true;
                    end
                    Fd=F.*ctc;
                    Md=M.*ctc;
                    DP=COP.*ctc;
                elseif isequal(endpoint,'Free')||ismissing(endpoint)
    
                    DP=obj.Joints.(s).(joints(end)).Child.EndPoints(end).Coordinates;
                    Fd=0*DP;
                    Md=0*DP;
                elseif isfield(obj.Joints.(s),endpoint)
                    jj=obj.Joints.(s).(endpoint);
                    Fd=pagemtimes(jj.Force,jj.Parent.Orientation);
                    Md=pagemtimes(jj.Moment,jj.Parent.Orientation);
                    DP=jj.JointCenter.Coordinates;  
                else 
                    error('endbody Attribute must be specified either as ForcePlatform, Free, or with a valid Joint');
                end
            
            for j=length(joints):-1:1
                joint=obj.Joints.(s).(joints(j));
                proxSegm=joint.Parent;
                JC=joint.JointCenter.Coordinates;
                [Fp,Mp]=DistDynCalc(joint,Fd,Md,g,DP,useID);
                %prepare for next iteration
                Fd=-Fp; % Newton third law
                Md=-Mp;
                DP=JC;
                %assign to joint
                Rprox=proxSegm.Orientation;
                mp=Mp;
                fp=Fp;
                for i=1:size(Rprox,3)
                    mp(i,:)=Mp(i,:)*Rprox(:,:,i);
                    fp(i,:)=Fp(i,:)*Rprox(:,:,i);
                end
                if u=="mm"
                obj.Joints.(s).(joint.Label).Moment=mp*1000;
                obj.Joints.(s).(joint.Label).MomentUnits="Nmm";
                elseif u=="m"
                obj.Joints.(s).(joint.Label).Moment=mp;
                obj.Joints.(s).(joint.Label).MomentUnits="Nm";
                else
                error("Invalid units, must  be m or mm!")
                end
                obj.Joints.(s).(joint.Label).Force=fp;
            end 
                catch ME
                    warning (ME.message)
                end
            end

        obj.setUnits(u);
        end
        %% ARITMETICS
        function obj=mean(obj)
        % MEAN Overloads the mean function
        %obj=MEAN(obj) applies the mean to the data in the acquisition and
        %returns a single frame object
                if obj.NFrames>1
                obj.Points=mean(obj.Points);
                obj.NFrames=1;
                end
                if obj.NSamples>1
                obj.ForcePlatform=mean(obj.ForcePlatform);
                obj.Analogs=mean(obj.Analogs);
                obj.Scalars=mean(obj.Scalars);
                obj.NSamples=1;
                end  
                obj.Metadata.POINT.FRAMES=1;
        end
        
        function obj=setUnits(obj,lengthunits)
        % SETUNITS Change the trial lenght units
        % obj=setUnits(obj,lengthunits) sets the units in the data to
        % either mm or m
            obj.ForcePlatform=obj.ForcePlatform.setUnits(lengthunits);
            obj.Points=obj.Points.setUnits(lengthunits);
            obj.Metadata.POINT.UNITS=char(lengthunits);
        end

        function obj=changeCoordinates(obj,T)
        % CHANGECOORDINATES Change the global coordinates
            obj.ForcePlatform=obj.ForcePlatform.changeCoordinates(T);
            obj.Points=obj.Points.changeCoordinates(T);
        end

        %% FILE I/O
        
        function targetfile=saveConfiguration(obj,targetfile)
        % SAVECONFIGURATION Save the configuration to an XML/JSON file
        % targetfile=SAVECONFIGURATION(obj,targetfile) saves the
        % configuration structure to an XML/JSON file
            if nargin==1
                try
            [p,f,e]=fileparts(obj.configfile);
            [f,p]=uiputfile(fullfile(p,[f e]),'Select Configuration file destination:');
                catch
                [f,p]=uiputfile('.configuration','Select Configuration file destination:'); 
                end
            if f==0 %no file selected
               return
            end
            targetfile=fullfile(p,f);
            end
            if obj.useJSON
            writestruct(obj.ConfigFile,targetfile,'FileType','json');
            else
            writestruct(obj.ConfigFile,targetfile,'FileType','xml','AttributeSuffix',obj.XMLatt);
            end
        end
        
        function targetfile=writeC3D(obj,targetfile,mode)
        % WRITEC3D Write the Trial to a C3D file
            arguments
                obj Trial
                targetfile="";
                mode="newfile";
            end
        if isempty(obj.C3DHandle)
        H=btkNewAcquisition(numel(obj.Points),obj.NFrames);
        btkSetFrequency(H,obj.Metadata.POINT.RATE);
        btkSetFirstFrame(H,obj.Metadata.TRIAL.ACTUAL_START_FIELD(1),0);
        obj.C3DHandle=H;
        else
        H=obj.C3DHandle;
        end
        if mode=="newfile"
        btkClearPoints(H); %remove point, angles and power
        btkClearAnalogs(H);%remove analogs
        try
        btkRemoveMetaData(H,'FORCE_PLATFORM'); %remove forceplates
        catch
        end
        btkClearEvents(H);%remove events
        btkSetFrameNumber(H,obj.NFrames);
        btkSetAnalogSampleNumberPerFrame(H,obj.NSamples/obj.NFrames);
        mode="append";


        end
        if targetfile==""
            [p,f,e]=fileparts(obj.c3dfile);
            [f,p]=uiputfile(fullfile(p,[f e]),'Select C3D file destination:');
            targetfile=fullfile(p,f);
            if f==0 %no file selected
               return
            end
        end

        fns=string(properties(obj))';
        mask=["Metadata","ConfigFile","XMLatt","configfile","c3dfile","useJSON","NFrames","NSamples","C3DHandle"];
        fns(matches(fns,mask))=[];
        for f=fns
                try
                    if matches(f,["Segments","Joints"]) %special behavior because data is structured
                       for s=string(fieldnames(obj.(f)))'
                           for i=string(fieldnames(obj.(f).(s)))'
                               obj.(f).(s).(i).updateC3D(H,mode);
                           end
                       end
                    else
                    obj.(f).updateC3D(H,mode);
                    end
                    
                catch ME
                    warning(ME.message);
                end
        end
        setC3DMetaData(obj); %also metadata has a special method
        btkWriteAcquisition(H,targetfile);
        end
       
       %% PLOT 
       function [line, cline]=stridePlot(obj,ax,var)
       ev=obj.Events.setUnitsandOffset('point',false);
       np=ax.NextPlot;
       ax.NextPlot='add';
       for side=["Left","Right","General"]
           try
           FS=ev.(side).Foot_Strike;
       for fs=1:length(FS)-1
           c=time2cycle([],var(FS(fs):FS(fs+1),:),101);
           line=plot(ax,0:100,c);
           cline=xline(ev.(side).Foot_Off,'FootOff');
       end
           catch
           end
       end
       ax.NextPlot=np;
       end
     
        function obj=readC3DMetaData(obj)
        % READC3DMETADATA Reads the C3D Metadata
            try
                groups=btkGetMetaData(obj.C3DHandle);
            catch ME
                warning(ME)
                return
            end
            groups=string(fieldnames(groups.children))';

           for g=groups
               labels=btkGetMetaData(obj.C3DHandle,char(g));
               labels=string(fieldnames(labels.children))';
               for l=labels
                   v=getfield(getfield(btkGetMetaData(obj.C3DHandle,char(g),char(l)),'info'),'values');
               MD.(g).(l)=v;
               end
               obj.Metadata=MD;
           end
        end

        function setC3DMetaData(obj)
        % SETC3DMETADATA Updateds the C3D Handle with current Metadata
            groups=string(fieldnames(obj.Metadata))';
            H=obj.C3DHandle;
           for g=groups
               labels=string(fieldnames(obj.Metadata.(g)))';
               for l=labels
                   value=obj.Metadata.(g).(l);
                   format=class(value);
                   format(1)=upper(format(1));
                   if format(1)=='U' %int8
                       format='Byte';
                   elseif format(1)=='D'
                       format='Real';
                   elseif format(1)=='C' || format(1)=='S'
                       format='Char';
                       value={char(value)};
                   end

                   INFO = btkMetaDataInfo(format, value);
                   try
                   btkAppendMetaData(H,char(g),char(l),INFO);
                   
                   catch %in the case the trial is not opened
                   end

               end
           end
           obj.readC3DMetaData;
        end   
     end

     methods (Access=private)
         function delete(obj)
             obj.C3DHandle=[];
         end
     end
 end

