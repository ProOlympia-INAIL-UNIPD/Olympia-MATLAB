function v=vec2linedata(O,V,scale)
arguments
    O (:,3,:,:) double
    V (:,3,:,:) double
    scale=1
end
O=pagetranspose(O);
O=O(:,:)';
V=pagetranspose(V);
V=V(:,:)';
if any(size(O)~=size(V))
    error('Input must be the same size!')
end

npt=size(O,1);
v=nan(npt*3-1,3);
V=O+normalize(V,2,'norm')*scale;
idx=1:2;
for i=1:npt
v(idx,:)=[O(i,:);V(i,:)];
idx=idx+3;
end