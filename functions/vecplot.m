function [l]=vecplot(ax,O,V,NameValue)
arguments
    ax 
    O (:,3,:) double
    V (:,3,:) double
    NameValue.ScaleFactor double =100;
    NameValue.?matlab.graphics.chart.primitive.Line
end

scale=NameValue.ScaleFactor;
prop=namedargs2cell(rmfield(NameValue,'ScaleFactor'));

[nvec]=size(O,3);
for i=nvec:-1:1
v=vec2linedata(O(:,:,i),V(:,:,i),scale);
l(i)=plot3(ax,v(:,1),v(:,2),v(:,3),prop{:});
end