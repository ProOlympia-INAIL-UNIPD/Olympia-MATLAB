classdef mocapView < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Axes
        MarkerData
        ForcePlatformData
        LocalSOR
        CurrentFrame=1        
        Markers
        ForcePlatforms
        Forces
        movFrame
        DrawList={'Markers','Labels','Forces','SOR'}
        SORscaler
        Fscaler
    end
    properties (Access = private)
        MarkerSR=nan
        ForceSR=nan
        SRRatio
        normGRF
        MarkerAvailable=false;
        ForceAvailable=false;
        SORAvailable=false;
    end
    methods
        function obj = mocapView(varargin)
            if (0<nargin)<2
                obj.MarkerData=varargin{1};
                obj.MarkerSR=obj.MarkerData.SampleRate;
                obj.MarkerAvailable=true;
            end    
            if nargin<=2
                obj.MarkerData=varargin{1};
                obj.MarkerSR=obj.MarkerData.SampleRate;
                obj.MarkerAvailable=true;
                obj.ForcePlatformData=varargin{2};         
                obj.ForceSR=obj.ForcePlatformData.SampleRate;
                obj.ForceAvailable=true;
                allfp=[obj.ForcePlatformData.GRF];
                fmax=max([allfp(:,3)],[],'all');
                hmax=max(obj.MarkerData.ZData,[],'all');
                normGRF.F=obj.ForcePlatformData.GRF*hmax/fmax;
                normGRF.P=obj.ForcePlatformData.COP;
                
                obj.normGRF=normGRF;
                %obj.normGRF.F=obj.normGRF.F*hmax/fmax;
            end
            if nargin<=3  
                obj.MarkerData=varargin{1};
                obj.MarkerSR=obj.MarkerData.SampleRate;
                obj.MarkerAvailable=true;
                obj.ForcePlatformData=varargin{2};         
                obj.ForceSR=obj.ForcePlatformData.SampleRate;
                obj.ForceAvailable=true;
                allfp=[obj.ForcePlatformData.GRF];
                fmax=max([allfp(:,3)],[],'all');
                hmax=max(obj.MarkerData.ZData,[],'all');
                obj.Fscaler=2000*hmax/fmax;
                normGRF.F=obj.ForcePlatformData.GRF*hmax/fmax;
                normGRF.P=obj.ForcePlatformData.COP;
                obj.normGRF=normGRF;         
                obj.LocalSOR=varargin{3};  
                obj.SORAvailable=true;
                %hmax=max(obj.MarkerData.ZData,[],'all');
 
                obj.SORscaler=hmax/10;
            end
        
        obj.SRRatio=obj.ForceSR/obj.MarkerSR;
        end
        
        function drawscene(obj,ax,frame)
            arguments
                obj
                ax=gca()
                frame=obj.CurrentFrame;
            end
            ax.NextPlot='replaceall';
            obj.Axes=ax;

            figure(obj.Axes.Parent);
            
               fk=frame;
               if length(fk)==1
                   if fk==1
                      ff=1;
                   else
                   ff=floor(fk*obj.SRRatio);
                   end
               else
                   idx=1:obj.SRRatio;
                   ff=zeros(1,obj.SRRatio*length(fk));
                   for i=1:length(fk)
                   ff(idx)=fk(i)+idx-1;
                   idx=idx+obj.SRRatio;
                   end
               end

            obj.Axes.NextPlot='replaceall';
            delete(obj.Axes.Children);
            if any(contains(obj.DrawList,'Forces'))
             
            [obj.Forces, obj.ForcePlatforms]=obj.ForcePlatformData.showGRF(obj.Axes,ff,obj.Fscaler);
            plotframe(obj.Axes,eye(4),'ScaleFactor',obj.SORscaler*1.5,'LineWidth',1.5)
            obj.Axes.NextPlot='add';
            obj.Axes.DataAspectRatio=[1 1 1];    
            obj.Axes.XLim=[min(obj.MarkerData.XData,[],"all") max(obj.MarkerData.XData,[],"all")]*2;
            obj.Axes.YLim=[min(obj.MarkerData.YData,[],"all") max(obj.MarkerData.YData,[],"all")]*1.1;
            obj.Axes.ZLim=[min(obj.MarkerData.ZData,[],"all") max(obj.MarkerData.ZData,[],"all")]*1.5;
            obj.Axes.Box='on';
            end
            if any(contains(obj.DrawList,'Markers'))
            obj.Markers.Line=plot3(obj.Axes,obj.MarkerData.XData(fk,:),obj.MarkerData.YData(fk,:),obj.MarkerData.ZData(fk,:),'o');
            [obj.Markers.Line.MarkerFaceColor]=deal(obj.MarkerData.Color(1,:));
            [obj.Markers.Line.Color]=deal(obj.MarkerData.Color(1,:));
            [obj.Markers.Line.MarkerSize]=deal(2);
            end
            if any(contains(obj.DrawList,'Labels'))
            obj.Markers.Text=text(obj.Axes,obj.MarkerData.XData(fk(1),:),obj.MarkerData.YData(fk(1),:),obj.MarkerData.ZData(fk(1),:),obj.MarkerData.Label,...
                "FontSize",6,'VerticalAlignment','middle',HorizontalAlignment='right',Color=obj.MarkerData.Color(1,:));
            end
            if any(contains(obj.DrawList,'SOR'))
            obj.movFrame=obj.LocalSOR.show(obj.Axes,fk,obj.SORscaler);
            end
            obj.Axes.NextPlot='replaceall';
            end

        function updatescene(obj,frame)
               if nargin==2
                   obj.CurrentFrame=frame;
               end
               fk=obj.CurrentFrame;
               if length(fk)==1
                   if fk==1
                      ff=1;
                   else
                   ff=floor(fk*obj.SRRatio);
                   end
               else
                   idx=1:obj.SRRatio;
                   ff=zeros(1,obj.SRRatio*length(fk));
                   for i=1:length(fk)
                   ff(idx)=fk(i)+idx-1;
                   idx=idx+obj.SRRatio;
                   end
               end
            
            if any(contains(obj.DrawList,'Markers'))   
            x=obj.MarkerData.XData(fk,:);
            y=obj.MarkerData.YData(fk,:);
            z=obj.MarkerData.ZData(fk,:);
            obj.Markers.Line.XData=x(:);
            obj.Markers.Line.YData=y(:);
            obj.Markers.Line.ZData=z(:);
            end

            %for i=1:length(obj.Markers.Text)
            if any(contains(obj.DrawList,'Labels'))
            P=[obj.MarkerData.XData(fk(1),:)',obj.MarkerData.YData(fk(1),:)',obj.MarkerData.ZData(fk(1),:)'];
            for i=1:size(P,1)
            obj.Markers.Text(i).Position=P(i,:);
            end
            end

            if any(contains(obj.DrawList,'Forces'))
            for i=1:length(obj.Forces)
            obj.ForcePlatformData.updatehg(ff,obj.Fscaler,obj.Forces)
            end
            end
            if any(contains(obj.DrawList,'SOR'))
            for i=1:length(obj.LocalSOR)
            obj.LocalSOR.updatehg(fk,obj.SORscaler,obj.movFrame);
            end
            end
        end

        function play(obj,replayspeed)
                figure(obj.Axes.Parent);
                if nargin==1
                    replayspeed=0.5;
                end
                sps=replayspeed*obj.MarkerSR;
                dt=1/sps;
                l=length(obj.MarkerData.XData(:,1));
                obj.CurrentFrame=obj.CurrentFrame(1);
                frame=obj.CurrentFrame;
                for i=frame:l
                t=tic;
                obj.updatescene;
                obj.CurrentFrame=i;
                while toc(t)<dt
                pause(dt/10);
                end
                end
                obj.CurrentFrame=frame;
end
    end

end