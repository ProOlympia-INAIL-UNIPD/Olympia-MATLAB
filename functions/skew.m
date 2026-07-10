function [S]=skew(v)
arguments
    v (:,3) double
end
for i=1:length(v(:,1))
S(:,:,i)=[0      -v(i,3) v(i,2)
          v(i,3)  0     -v(i,1)
         -v(i,2)  v(i,1) 0];
end
end

