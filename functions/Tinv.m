function Ti = Tinv (T)
Ti = nan(3,3,size(T,3));
for i = 1:size(T,3)
    Ti(1:3,1:3,i) = T(1:3,1:3,i)';
    Ti(1:3, 4, i) = -Ti(1:3,1:3,i)*T(1:3, 4, i);
end
Ti(4,4,:) = 1;
