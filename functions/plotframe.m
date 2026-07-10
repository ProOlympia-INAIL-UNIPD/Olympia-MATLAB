function [l]=plotframe(ax,T,NameValue)
arguments
    ax 
    T (4,4,:,:) double
    NameValue.ScaleFactor double =1;
    NameValue.?matlab.graphics.chart.primitive.Line
end
T=reshape(T,4,4,[],1);
%plot pose matrix T in 3D space input must be the 4x4 pose matrix
%Multiplier sets the magnitude of the arrow
%LineWidth sets the width of the arrow
%MarkerSize sets the size of the origin marker
scale=NameValue.ScaleFactor;
prop=namedargs2cell(rmfield(NameValue,'ScaleFactor'));
nrf=size(T,3);
X=nan(nrf*3-1,3);
Y=nan(nrf*3-1,3);
Z=nan(nrf*3-1,3);
idx=1:2;
for i=1:nrf
o=T(1:3,4,i)';
r=T(1:3,1:3,i)';
X(idx,:)=vec2linedata(o,r(1,:,:),scale);
Y(idx,:)=vec2linedata(o,r(2,:,:),scale);
Z(idx,:)=vec2linedata(o,r(3,:,:),scale);
idx=idx+3;
end
ax.NextPlot='add';
l(3)=plot3(ax,Z(:,1),Z(:,2),Z(:,3),prop{:});
l(3).Color=[0 0 1];
l(1)=plot3(ax,X(:,1),X(:,2),X(:,3),prop{:});
l(1).Color=[1 0 0];
l(2)=plot3(ax,Y(:,1),Y(:,2),Y(:,3),prop{:});
l(2).Color=[0 1 0];
ax.NextPlot='replaceall';


%% old version
%{

    p=scatter3(ax,Oframe(1),Oframe(2),Oframe(3),NameValue.MarkerSize,NameValue.MarkerFaceColor,'filled');
    hold on
    qx=quiver3(ax,Oframe(1),Oframe(2),Oframe(3),Rframe(1,1),Rframe(2,1),Rframe(3,1),'Color',[1 0 0],'LineWidth',NameValue.LineWidth,'AutoScaleFactor',NameValue.AutoScaleFactor);
    qy=quiver3(ax,Oframe(1),Oframe(2),Oframe(3),Rframe(1,2),Rframe(2,2),Rframe(3,2),'Color',[0 1 0],'LineWidth',NameValue.LineWidth,'AutoScaleFactor',NameValue.AutoScaleFactor);
    qz=quiver3(ax,Oframe(1),Oframe(2),Oframe(3),Rframe(1,3),Rframe(2,3),Rframe(3,3),'Color',[0 0 1],'LineWidth',NameValue.LineWidth,'AutoScaleFactor',NameValue.AutoScaleFactor);
    hold off
    xlabel('X[mm]');
    ylabel('Y[mm]');
    zlabel('Z[mm]'); 
    axis equal

end
%}