classdef forcePlatformType2
    %   forcePlatformType2 is a Force Platform object following Type2
    %   specs.
    %   Type2 Force Platform are described exhaustively in https://www.c3d.org/docs/C3D_User_Guide.pdf
    %   The elements characterizing a force platforms are its position with
    %   respect to the global reference frame, defined by the four vertexs
    %   (Corners). Each platform has six physical channels each defining a
    %   component of the force (ch1 to ch3) and of the moment vector (ch4
    %   to ch6). Forces and Moments are referred to the Origin of the force
    %   platform and their components are oriented according to a Reference
    %   System built on the four corners. The Physical Origin of the force platform (where an applied force would give 0 moment) may differ from the
    %   centroid of the corners. This could be the case when the platform
    %   surface is covered and therefore the working plane (defined by the
    %   Corners) differs from the platform surface.
    properties
        Label       (1,1) string% name of the forceplate
        Corners     (4,3) double% position of the four corners (corner,XYZ)
        Origin      (1,3) double% Origin with respect to working plate (XYZ)
        Channels    (1,1) struct% physical channels of Force(sample,xyz) and Moment(sample,xyz) in local coordinates
        GRF         (:,3) double% Ground Reaction Forces in absolute coordinates (sample,XYZ)
        GRM         (:,3) double% Ground Reaction Moments in absolute coordinates (sample,XYZ)
        COP         (:,3) double% COP position in absolute coordinates (sample,XYZ)
        Units       (1,3) string% Forces, Moment, and COP units
        SampleRate  (1,1) double% sample rate of the acquisition
        Color       (1,3) double% display color of the forceplate in the 3D plots
        NSamples    (1,1) double% number of analog samples
        Parent      (1,1) Trial % reference to the Trial containing force platform data
    end

    methods
        function obj=forcePlatformType2(varargin)
            %constructs an instance of this class by passing a C3D file or btk handle
            closeflag=true;
            if nargin==0
               return
            end
            
            c3dfile=varargin{1};
            
            if isnumeric(c3dfile) %btk Handle
                H=c3dfile;
                closeflag=false;
            elseif isfile(c3dfile)%file
                H=btkReadAcquisition(c3dfile);
            % 
            % else %nargin=0
            %     [f,p]=uigetfile('*.c3d');
            %     if f==0 %no file
            %        return
            %     end
            %     c3dfile=fullfile(p,f);
            end
                   
            [FP,fpinfo]=btkGetForcePlatforms(H);
            %platform check
            if isempty(FP)
                warning('No ForcePlatform data available!')
                obj(:)=[];
                return
            elseif any(not([FP.type]==2))
                error('one or more Force Platforms differ from Type2 specifications!')
            end

            %color attribution
            color=colororder;
            while size(color,1)<length(FP)
               color=[color; color];
            end
            color=color(1:length(FP),:);

            for i=length(FP):-1:1
                fpunits=fpinfo(i).units;
                unname=fieldnames(fpunits);
                for u=1:length(unname)
                    units{u}=fpunits.(unname{u});
                end
                units=unique(units);
                units{3}=btkGetPointsUnit(H,'marker');      
                   
    
                
                oldnames=fieldnames(FP(i).channels);
                for j=1:3
                ch.Force(:,j)=FP(i).channels.(oldnames{j});
                ch.Moment(:,j)=FP(i).channels.(oldnames{j+3});
                end

                obj(i).Label =sprintf("FP%i",i);
                obj(i).Channels=ch;
                obj(i).Corners=transpose(FP(i).corners);
                obj(i).Origin=FP(i).origin';
                obj(i).Color = color(i,:);
                obj(i).Units = string(units);
                obj(i).SampleRate=btkGetAnalogFrequency(H);
                obj(i).NSamples=btkGetAnalogFrameNumber(H);
                obj(i)=obj(i).getGRF;
            end
                
            if closeflag
            btkCloseAcquisition(H);
            end
            
        end

        function varargout=extractSOR(obj)
            %   obtains the Force Platform SOR and returns either 
            %   Rotation (3,3) Matrix and Origin (XYZ) or Transformation (4,4) Matrix
            for i=numel(obj):-1:1
            corners=obj(i).Corners;
            FPR(:,1,i)=normalize(corners(1,:)-corners(2,:),2,'norm');
            FPR(:,2,i)=normalize(corners(1,:)-corners(4,:),2,'norm');
            FPR(:,3,i)=normalize(cross(FPR(:,1,i),FPR(:,2,i)),1,'norm');
            OFP(:,i)=mean(corners(:,:));
            T(:,:,i)=[FPR(:,:,i) OFP(:,i);0 0 0 1];
            end

            if nargout==1
               varargout{1}=T;
            elseif nargout==2
                varargout{1}=FPR;
                varargout{2}=OFP';
            end
        end

        function obj = getGRF(obj)
        %   computes the GRF, GRM, COP properties using Corners, Origin and Channels 
            for i=length(obj):-1:1
                [FPR,OFP]=extractSOR(obj(i));
                floc=obj(i).Channels.Force; %local forces
                mloc=obj(i).Channels.Moment;%local moments (with respect to physical origin)
                mloc=mloc+(skew(obj(i).Origin)*floc')';%local moments (with respect to woring plane origin)
                obj(i).GRF=floc*FPR'; %Ground Reaction Forces
                mgnd=mloc*FPR'; %Ground Reaction Moment (in FP working origin)
                obj(i).GRM=mgnd+(skew(OFP)*obj(i).GRF')'; %GRM in global origin  (0,0,0);
                cop(:,1)=-mloc(:,2)./floc(:,3); %COP in FP working plane
                cop(:,2)=mloc(:,1)./floc(:,3);
                cop(:,3)=0;                     %COP lies in working plane!
                obj(i).COP=OFP+cop*FPR'; %COP in absolute system
            end
          
        end

        function obj = cleanSignal(obj,NameValue)
            % cleans the force platform data reducing noise and cross-talk
            % NameValue properties control the parameters. 
            % -MaxRadius: sets the maximum allowed half-duration of a contact (in samples),
            % -ActiveThreshold: sets the minimum force to consider the
            % forceplatform as active, 
            % -MaxNumContacts: sets the maximum number of contacts which will 
            % be considered for each force platform
            % -FilterOrder: sets the order for the Butterworth low-pass
            % filter (will be doubled in filtfilt)
            % -FilterCutOff: sets the cut-off frequency in Hz for the
            % filter.
            arguments
                obj
                NameValue.MaxRadius=200;      %expected half-duration of the contact (in samples) 
                NameValue.FilterCutOff=100;   %cut-off frequency for low pass filter
                NameValue.FilterOrder=2;      %low pass filter order (will be doubled in filtfilt)
                NameValue.ActiveThreshold=100;%minimum threshold to accept FP as used during the acquisition
                NameValue.MaxNumContacts=100;   %maximum number of hits accepted for each FP
                NameValue.Reflect=true;
            end
            threshold=20;
            maxradius=NameValue.MaxRadius;
             %build filter
            for i=length(obj):-1:1
                [b,a]=butter(NameValue.FilterOrder,NameValue.FilterCutOff/(obj(i).SampleRate/2));
                Fnow=-obj(i).Channels.Force; %work on the normal force (i.e., vertical)
                isactive=not(Fnow(:,3)==0);  %missing values are treated as 0, non 0 elements are signals from the FP)
                %isactive=abs(Fnow(:,3))>threshold;
                if numel(obj)>1 %cross talk removal is meaningful only with multiple forceplates
                warning off %suppress warnings from findpeaks in case no peak is detected
                [~, loc]=findpeaks(Fnow(:,3),"NPeaks",NameValue.MaxNumContacts,"MinPeakHeight",NameValue.ActiveThreshold,"MinPeakDistance",NameValue.MaxRadius);
                warning on
                maxzone=[];

                for j=1:length(loc)
                    maxzone=[maxzone max(1,loc(j)-maxradius):min(loc(j)+maxradius,length(isactive))];
                end

                ismaxzone=false(length(isactive),1); %convert maxzone to logical 1) set all elements to false
                ismaxzone(maxzone)=1; %convert maxzone to logical 2) set maxzone elements to true
                else
                    ismaxzone=isactive;
                end
                isvalid=isactive & ismaxzone;        % active and valid signal
                changeState=[false; diff(isvalid); false];
                fc=find(changeState==1)-1;
                fo=find(changeState==-1)-1;
                dur=fo-fc;
                dur=find(dur<=5);
                for d=dur'
                    isvalid(fc(d):fo(d))=false;
                end
                
                changeState=[false; diff(isvalid); false];
                fc=find(changeState==1)-1;
                fo=find(changeState==-1)-1;

                F=double(obj(i).Channels.Force);
                M=double(obj(i).Channels.Moment);
                if NameValue.Reflect
                   changeState=[false; diff(isvalid); false];
                   fc=find(changeState==1)-1;
                   fo=find(changeState==-1)-1;

                   for j=1:length(fc)
                       w=floor((fo(j)-fc(j))/2);
                       pre=fc(j):-1:max(1,fc(j)-w);
                       post=fo(j):min(length(F(:,1)),fo(j)+w);
                       F(pre,:)=-F(fc(j):fc(j)+length(pre)-1,:);
                       F(post,:)=-F(fo(j):-1:fo(j)-length(post)+1,:);
                       M(pre,:)=-M(fc(j):fc(j)+length(pre)-1,:);
                       M(post,:)=-M(fo(j):-1:fo(j)-length(post)+1,:);
                   end

                end
                obj(i).Channels.Force=filtfilt(b,a,F).*isvalid;
                obj(i).Channels.Moment=filtfilt(b,a,M).*isvalid;
            end
            obj=obj.getGRF;
        end

        function obj = combineFP(obj,newcorners,neworigin,index)
            %combines multiple force platform into a single one.
            % this function takes a forceplatform array and combines them
            % together all its elements, otherwise specified by index.
            % the corners and the origin can be either selected or the
            % program will auto calculate them to generate a force platform
            % covering the selected force platforms.
            arguments
                obj
                newcorners=[];
                neworigin=[];
                index=1:numel(obj);
            end
            index=sort(index,'ascend');
            if isscalar(obj)
               warning('The ForcePlatform number is already 1!')
               return
            elseif numel(index)<2
                warning('At least 2 ForcePlatform must be selected to be combined!');
                return
            end
            if isempty(newcorners)
            T=obj(1).extractSOR;
            cloc=zeros(4,3,index(end));
            for i=index               
                for j=1:4
                cloc(j,:,i)=points2local(obj(i).Corners(j,:),T);
                end
            end
            if any(diff(cloc(:,3,:)),'all')
               error('Unable to combine ForcePlates with different working planes!"');
            end
            C=[ max(cloc(:,1,:),[],'all'),max(cloc(:,2,:),[],'all'),0;
                max(cloc(:,1,:),[],'all'),min(cloc(:,2,:),[],'all'),0;
                min(cloc(:,1,:),[],'all'),min(cloc(:,2,:),[],'all'),0;
                min(cloc(:,1,:),[],'all'),max(cloc(:,2,:),[],'all'),0];
            C(:,4)=1;
            
            C=(T*C')';
            newcorners=C(:,1:3);
            end
            if isempty(neworigin)
                neworigin=obj(1).Origin;
            end
            old_obj=obj;
            [RFP,OFP]=old_obj.extractSOR;
            obj=obj(1);
            obj.Corners=newcorners;
            obj.Label="FP"+num2str(index(1));

            obj.Origin=neworigin;
            [RFPn,~]=obj.extractSOR;
            T=obj.extractSOR;
            OFP(:,end+1)=1;
            OFP=(T\OFP')';
            OFP(:,4)=[];
            floc=zeros(obj.NSamples,3,index(end));
            mloc=zeros(obj.NSamples,3,index(end));

         for i=index
             R=RFPn'*RFP(:,:,i);
             floc(:,:,i)=old_obj(i).Channels.Force*R';
             mloc(:,:,i)=old_obj(i).Channels.Moment*R';
             mloc(:,:,i)=mloc(:,:,i)+(skew(OFP(i,:))*floc(:,:,i)')';       
         end
         ch.Force=sum(floc,3);
         ch.Moment=sum(mloc,3);
         obj.Channels=ch;
         
         obj=obj.getGRF;
         old_obj(index)=[];
         obj=[old_obj,obj];
        end

        function [FC, FO]=getEvents(obj,units,treshold,startdetectionoffset)
            % detects contacts on the force plates using user-specified thresholds 
            % the function returns the foot contact and foot off events in
            % samples, seconds, or milliseconds. Contacts are defined when
            % the normal force module is over the defined threshold.
            % startdetectionoffset could be used to remove events before a
            % certain sample.
            arguments
                obj
                units {mustBeMember(units,{'samples','seconds','milliseconds'})}='samples';
                treshold=10;
                startdetectionoffset=0;
            end

            R=[obj.Channels];
            R=[R.Force(:,3)];
            R=sqrt(sum(R.^2,2));
            iscontact=any(R>treshold,2);
            FC=find(diff(iscontact)==1)+1;
            FO=find(diff(iscontact)==-1);
            FC(FC<startdetectionoffset)=[];
            FO(FO<startdetectionoffset)=[];
            if FC(1)>FO(1)
               FC=[1; FC];
            end
            if FO(end)<FC(end)
               FO=[FO; size(R,1)];
            end
            switch units
                case 'samples'
                case 'seconds'
                    FC=(FC-1)/obj.SampleRate;
                    FO=(FO-1)/obj.SampleRate;
                case 'milliseconds'
                    FC=1000*(FC-1)/obj.SampleRate;
                    FO=1000*(FO-1)/obj.SampleRate;
            end
           
        end

        function obj=resample(obj,targetSampleRate)
                 % downsamples Force plate signal to a desired rate.
                 srratio=obj(1).SampleRate/targetSampleRate;
                 if mod(obj(1).SampleRate,targetSampleRate)>0 || srratio<1
                     error('Only integer Downsampling ratio is allowed!');
                 end
                 for i=1:numel(obj)
                     obj(i).Channels.Force=obj(i).Channels.Force(1:srratio:end,:);
                     obj(i).Channels.Moment=obj(i).Channels.Moment(1:srratio:end,:);
                     obj(i).SampleRate=targetSampleRate;
                     obj(i).NSamples=obj(i).NSamples/srratio;
                 end
                 obj=obj.getGRF;

        end

        function obj=mean(obj)
            % average force platform data
            for i=length(obj):-1:1
            newch=structfun(@(x) mean(x,1),obj(i).Channels,'UniformOutput',false);
            
            obj(i).Channels=newch;
            obj(i)=obj(i).getGRF;
            obj(i).SampleRate=nan;
            obj(i).NSamples=1;
            end
        end

        function obj=changeCoordinates(obj,T)
        % apply transformation matrix to force platform position
            for i=1:length(obj)
                c=obj(i).Corners;
                c(:,4)=1;
                c=(T*c')';
                obj(i).Corners=c(:,1:3);
                o=obj(i).Origin;
                o(:,4)=1;
                o=(T*o')';
                obj(i).Origin=o(:,1:3);
            end
            obj=obj.getGRF; 
        end

        function obj=setUnits(obj,newunits)
            % sets the length units for the force platform to m or mm
            mustBeMember(newunits,["mm","m"]);
            for i=1:numel(obj)
                if obj(i).Units(3)==newunits
                   sF=1;
                elseif obj(i).Units(3)=="m"
                    sF=1000;
                elseif obj(i).Units(3)=="mm"
                    sF=1e-3;
                end
                obj(i).Channels.Moment=obj(i).Channels.Moment*sF;
                obj(i).Corners=obj(i).Corners*sF;
                obj(i).Origin=obj(i).Origin*sF;
                obj(i).Units(2)=strcat(obj(i).Units(1),newunits);
                obj(i).Units(3)=newunits;
            end
            obj=obj.getGRF;    
        end

        %% UTILITIES
        function R=align2ISB(obj)
                % compute the matrix bringing the coordinates to ISB
                % convention.
                % Use the platform data to align the current acquisition to ISB conventions: X
                % antero-posterior, Y vertical, Z medio-lateral. The
                % function use the pca function to determine the running
                % direction (maximum COP variability), lateral direction
                % (intermediate COP variability) and vertical direction (no
                % COP variability).
                
        if isscalar(obj)
            COP=obj.COP;
            % align the COP to the principal directions
            R=pca(COP);
            R=R(:,[1 3 2])';
        else %ask the user to manually select the directions
            notok=true;
                while notok
                xyz=inputdlg(["AnteroPosterior:","MedioLateral:","Vertical:"],"Enter trial direction",[1 30],["+Y","+X","+Z"]);
                xyz=char(xyz);
                s=xyz(:,1);
                s=s=='-';
                xyz=double(xyz(:,end)-'X')+1;
                R=eye(3);
                R=R(:,xyz);
                R(:,s)=-R(:,s);
                    if det(R)~=1
                        switch questdlg("Invalid set of axes, Try again?")
                        case "Yes"
                            continue
                        case "No"   
                            return
                        otherwise
                            return
                        end
                    end
                notok=false;
                end
        end
        end
        %% GRAPHICS
        function [p, ln]=show(obj,ax)
            % show a 3D-view of the Force plate in the current or user specified axes
            if nargin==1
                ax=gca();
            end
            ax.NextPlot='add';
            ax.View=[30 -30];
            ax.DataAspectRatio=[1 1 1];
            T=obj.extractSOR;
            for i=length(obj):-1:1
            p(i)=fill3(ax,obj(i).Corners(:,1),obj(i).Corners(:,2),obj(i).Corners(:,3),obj(i).Color);
            end
            ln=plotframe(ax,T,'ScaleFactor',min(range(obj(i).Corners(:,1:2))/2),'LineWidth',1.2);

        end
        function [f,p,l]=showGRF(obj,ax,ff,scale)
                 % show force plate surface and force vector in the user-specified
                 % axes for the frame(s) defined by ff. Applies a scale
                 % factor to the force with respect to world coordinates
                 arguments
                     obj
                     ax=gca();
                     ff=1;
                     scale=50;
                 end
                 [p,l]=obj.show(ax);
                 ax=p.Parent;
                 ax.NextPlot='add';
                 COP=cat(3,obj.COP);
                 GRF=cat(3,obj.GRF);
                 f=vecplot(ax,COP(ff,:,:),GRF(ff,:,:),ScaleFactor=scale);
                 for i=1:length(f)
                 f(i).Color=obj(i).Color;
                 end
                 ax.NextPlot='add';
        end
        function updatehg(obj,f,scale,fv)
            %updates force vector in an open plot
        arguments 
            obj
            f=1
            scale=50;
            fv=obj.showGRF(f);         
        end
        for i=1:length(fv)
        V=vec2linedata(obj(i).COP(f,:),obj(i).GRF(f,:),scale);
        fv(i).XData=V(:,1);
        fv(i).YData=V(:,2);
        fv(i).ZData=V(:,3);
        end
        end
        function [f, m]=plotChannels(obj,fig)
         %plot force plate channels in a new or existing figure
                arguments
                    obj
                    fig=figure()
                 end
                 
                 figure(fig);
                 t=(0:size(obj(1).Channels.Force(:,1))-1)/obj.SampleRate;
                 for i=1:length(obj)
                 subplot(2,1,1)
                 f(i,:)=plot(t,obj(i).Channels.Force);
                 hold on
                 lgF(:,i)=strcat(obj.Label(i),{'.Fx','.Fy','.Fz'});
                 subplot(2,1,2)
                 m(i,:)=plot(t,obj(i).Channels.Moment);
                 hold on
                 lgM(:,i)=strcat(obj.Label(i),{'.Mx','.My','.Mz'});
                 end
                 subplot(2,1,1)
                 
                 title('Local Forces')
                 legend(lgF{:},'NumColumns',length(obj.Channels),'FontSize',8);
                 xlabel('time (s)');
                 ylabel(sprintf('Force (%s)',obj.Units{1}));

                 subplot(2,1,2)
                 title('Local Moments')
                 legend(lgM{:},'NumColumns',length(obj.Channels),'FontSize',8);
                 xlabel('time (s)');
                 ylabel(sprintf('Moment (%s)',obj.Units{2}));
        end 
        function [grf, grm, cop]=plotGRF(obj,fig)
        %plot GRFs in a new or existing figure
                arguments
                    obj
                    fig=figure()
                 end
                 
                 t=(0:obj(1).NSamples-1)/obj(1).SampleRate;
                 figure(fig);
                 for i=1:numel(obj)
                 subplot(3,1,1)
                 grf(i,:)=plot(t,obj(i).GRF);
                 hold on
                 lgF(:,i)=strcat(obj(i).Label,{'.GRFx','.GRFy','.GRFz'});
                 subplot(3,1,2);
                 grm(i,:)=plot(t,obj(i).GRM);
                 hold on
                 lgM(:,i)=strcat(obj(i).Label,{'.GRMx','.GRMy','.GRMz'});
                 subplot(3,1,3);
                 cop(i,:)=plot(t,obj(i).COP);
                 hold on
                 lgC(:,i)=strcat(obj(i).Label,{'.COPx','.COPy','.COPz'});
                 end
                 subplot(3,1,1)
                 title('Ground Reaction Forces')
                 legend(lgF{:},'NumColumns',length(obj),'FontSize',8);
                 xlabel('time (s)');
                 ylabel(sprintf('Force (%s)',obj(1).Units{1}));

                 subplot(3,1,2)
                 title('Ground Reaction Moments')
                 legend(lgM{:},'NumColumns',length(obj),'FontSize',8);
                 xlabel('time (s)');
                 ylabel(sprintf('Moment (%s)',obj(1).Units{2}));
                 
                 subplot(3,1,3)
                 title('Center of Pressure')
                 legend(lgC{:},'NumColumns',length(obj),'FontSize',8);
                 xlabel('time (s)');
                 ylabel(sprintf('COP (%s)',obj(1).Units{3}));
        end 

        %% C3D file manipulation
        function updateC3D(obj,H,mode)
            % updates the btk handle associated with the c3dfile 
            % the user can control the mode and decide if overwrite the
            % acquisition with the new data or append the forceplatform to
            % the existing acquisition

        switch mode
            case 'overwrite'
                fpdata=btkGetMetaData(H,'FORCE_PLATFORM');
                fpanalogs=fpdata.children.CHANNEL.info.values;
                fpanalogs_array=sort(fpanalogs(:));
                nfp=numel(obj);
                objanalogs=nfp*6;
                for i=1:nfp
                btkSetAnalogValues(H,fpanalogs(1,i),double(obj(i).Channels.Force(:,1))); %Fx
                btkSetAnalogValues(H,fpanalogs(2,i),double(obj(i).Channels.Force(:,2))); %Fy
                btkSetAnalogValues(H,fpanalogs(3,i),double(obj(i).Channels.Force(:,3))); %Fz
                btkSetAnalogValues(H,fpanalogs(4,i),double(obj(i).Channels.Moment(:,1)));%Mx
                btkSetAnalogValues(H,fpanalogs(5,i),double(obj(i).Channels.Moment(:,2)));%My
                btkSetAnalogValues(H,fpanalogs(6,i),double(obj(i).Channels.Moment(:,3)));%Mz
                btkSetAnalogUnit(H,fpanalogs(1,i),char(obj(i).Units(1)));
                btkSetAnalogUnit(H,fpanalogs(2,i),char(obj(i).Units(1)));
                btkSetAnalogUnit(H,fpanalogs(3,i),char(obj(i).Units(1)));
                btkSetAnalogUnit(H,fpanalogs(4,i),char(obj(i).Units(2)));
                btkSetAnalogUnit(H,fpanalogs(5,i),char(obj(i).Units(2)));
                btkSetAnalogUnit(H,fpanalogs(6,i),char(obj(i).Units(2)));
                end
                while objanalogs<numel(fpanalogs_array)
                    btkRemoveAnalog(H,fpanalogs_array(end));
                    fpanalogs_array(end)=[];
                end
            case 'append'
                %overwrite=true;
                na=btkGetAnalogNumber(H);
                for i=1:numel(obj)
                    warning off
                    btkAppendForcePlatformType2(H,double(obj(i).Channels.Force),double(obj(i).Channels.Moment),double(obj(i).Corners),double(obj(i).Origin),1);      
                    warning on
                    for x=1:6
                        if x<4
                    btkSetAnalogUnit(H,na+x,char(obj(i).Units(1)));
                        else
                    btkSetAnalogUnit(H,na+x,char(obj(i).Units(2)));
                        end
                    end
                    na=btkGetAnalogNumber(H);
                end
        end
    obj(1).Parent.readC3DMetaData;
    end

function [c,cmedio,cstd]=timeNorm(obj,npoints,base)
        %returns GRF and COP data time-normalized for each stride, as well
        %as the average cycle and its standard deviation
            arguments
                obj forcePlatformType2
                npoints=101;
                base string {mustBeMember(base,["stance","stride"])}="stance"
            end
            ev=exportEvents(obj.Parent.Events,'point',false);
            group = ["Left","Right"];
            
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
                c.(g+"GRF")(:,:,i)=time2cycle([],obj.GRF(S(i):E(i),:),npoints);
                c.(g+"COP")(:,:,i)=time2cycle([],obj.COP(S(i):E(i),:),npoints);
                catch
                end
            end
            end 
            cmedio=structfun(@(x) mean(x,3,'omitnan'),c,"UniformOutput",false);
            cstd=structfun(@(x) std(x,1,3,'omitnan'),c,"UniformOutput",false);
        end
    end
end


        % function c3dfile=write2C3D(obj,c3dfile,mode)
        %     %saves the acquisition to file
        %     arguments
        %         obj
        %         c3dfile char=obj.c3dfile;
        %         mode char {mustBeMember(mode,{'newfile','append','overwrite'})} ='newfile';
        %     end
        % 
        %     h=btkReadAcquisition(obj.c3dfile);
        %     updateC3D(obj,h,mode)
        %     %[path,file,ext]=fileparts(c3dfile);
        % 
        %     % if overwrite==false
        %     %     k=1;
        %     %     c3dfile=fullfile(path,[file,ext]);
        %     %     while isfile(c3dfile)
        %     %         c3dfile=fullfile(path,[file,'(',num2str(k),')',ext]);
        %     %         k=k+1;
        %     %     end
        %     % else
        %     % end
        %     btkWriteAcquisition(h,c3dfile);
        %     btkCloseAcquisition(h);
        % 
        % end