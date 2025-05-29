function showROS(GRFcycle, COPcycle, method, color, legendtag)
arguments
    GRFcycle (:,3,:) double
    COPcycle (:,3,:) double
    method char {mustBeMember(method,{'default','deltaROS'})}='default'
    color (1,3) double = [0,0,0];
    legendtag string = "";
end
s=2;

COPx = squeeze(COPcycle(:,1,:));
COPy = squeeze(COPcycle(:,2,:));

if isequal(method, 'deltaROS')
    COPx = COPx - COPx(1,1);
    COPy = COPy - COPy(1,1);
end

GRFx = squeeze(GRFcycle(:,1,:));
GRFy = squeeze(GRFcycle(:,2,:));

pl = plot(COPx, COPy, 'color', color);
if size(GRFcycle, 3) > 1
    i=2;
    while i<=size(GRFcycle, 3)
        pl(i).HandleVisibility='off';
        i=i+1;
    end
end
hold on

quiver(COPx(1:s:end,:),COPy(1:s:end,:), GRFx(1:s:end,:),GRFy(1:s:end,:), 'color', color, 'HandleVisibility', 'off')

axis equal
grid on;

ylabel("y-axis (mm)");
xlabel("x-axis (mm)");

FS=plot(COPx(1,:), COPy(1,:), '*k', 'MarkerFaceColor', color);
FO=plot(COPx(end,:) , COPy(end,:), 'ok', 'MarkerFaceColor', color);

legend(legendtag,'FS','FO')

end
