function [gnd_T_loc, gnd_R_loc] = buildCS(points, CSstruct)
% --------------------- ORIGIN DEF -------------------- %
nPor = size(CSstruct.origin, 2);
nF = size(points.(CSstruct.origin(1)), 1);

origin = nan(nF, 3, nPor);
for i = 1:nPor
    origin(:, :, i) = points.(CSstruct.origin(i));
end
origin = mean(origin, 3, 'omitmissing');

points.origin = origin; %locally store the origin as some axes/temp vectors might be defined as a function of it

% -------------------- PRIMARY AXIS ------------------- %
nAxis = length(CSstruct.axis1.axis);
aux = nan(nF, 3, nAxis);

for iA = 1:nAxis
    nPtail = length(CSstruct.axis1.axis(iA).tail);
    nPtip = length(CSstruct.axis1.axis(iA).tip);
    
    tail = nan(nF, 3, nPtail);
    for i = 1:nPtail
        tail(:, :, i) = points.(CSstruct.axis1.axis(iA).tail(i));
    end
    tail = mean(tail, 3, 'omitmissing');

    tip = nan(nF, 3, nPtip);
    for i = 1:nPtip
        tip(:, :, i) = points.(CSstruct.axis1.axis(iA).tip(i));
    end
    tip = mean(tip, 3, 'omitmissing');

    aux(:, :, iA) = tip - tail;
end
switch nAxis
    case 1
        axis1 = normalize(aux, 2, 'norm');
    case 2 % nPtail and nPtip have same number of points but their length is grater than 1
        axis1 = normalize(cross(aux(:, :, 1), aux(:, :, 2)), 2, 'norm');
end

% ------------------- SECONDARY AXIS ------------------ %
nAxis = length(CSstruct.temp_axis.axis);
aux = nan(nF, 3, nAxis);

for iA = 1:nAxis
    nPtail = length(CSstruct.temp_axis.axis(iA).tail);
    nPtip = length(CSstruct.temp_axis.axis(iA).tip);
 
    tail = nan(nF, 3, nPtail);
    for i = 1:nPtail
        tail(:, :, i) = points.(CSstruct.temp_axis.axis(iA).tail(i));
    end
    tail = mean(tail, 3, 'omitmissing');

    tip = nan(nF, 3, nPtip);
    for i = 1:nPtip
        tip(:, :, i) = points.(CSstruct.temp_axis.axis(iA).tip(i));
    end
    tip = mean(tip, 3, 'omitmissing');

    aux(:, :, iA) = tip - tail;
end
switch nAxis
    case 1
        temp_axis = normalize(aux, 2, 'norm');
    case 2
        temp_axis = normalize(cross(aux(:, :, 1), aux(:, :, 2)), 2, 'norm');
end

axis2 = normalize(cross(axis1, temp_axis), 2, 'norm');

% --------------------- THIRD AXIS -------------------- %
axis3 = normalize(cross(axis2, axis1), 2, 'norm');

if ismember(CSstruct.seq, {'xyz', 'yzx', 'zxy'})
    axis3 = -axis3;
end

% ------------- NAME AXIS AS SEQUENCE ----------------- %
% Rotation matrix: CS attitude
% Version 1: haarder to interpret, but nicer to look at
axes = permute(cat(3, axis1, axis2, axis3), [2 3 1]);
seq = convertseq(CSstruct.seq{:});
gnd_R_loc = cat(2, axes(:, seq==1, :), axes(:, seq==2, :), axes(:, seq==3, :));

% Version 2: easier to read, but less elegant 
%eval([SoRstruct.seq{:}(1), '= axis1;'])
%eval([SoRstruct.seq{:}(2), '= axis2;'])
%eval([SoRstruct.seq{:}(3), '= axis3;'])

%gnd_R_loc = permute(cat(3, x, y, z), [2 3 1]);

% Homogeneus transformation matrix: CS pose
gnd_T_loc = gnd_R_loc;
gnd_T_loc(4, 4, :) = 1;
gnd_T_loc(1:3, 4, :) = permute(origin, [2 3 1]);

    function seqArray = convertseq(seq)
        seq=upper(seq); %switch lowercase to UPPERCASE in input
        seq=unique(seq, 'stable');
        if length(seq) ~= 3
            error("Sequences must be of three elements.");
        end
        mustBeMember(seq, ['X' 'Y' 'Z']);
        
        seqArray=double(seq)-87; %convert XYZ to 123 (XZY to 132 etc...)
    end
end