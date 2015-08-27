%%
scrsz = get(0,'screensize');
% frameRate = 10;

%%
q = 0;
if ~isempty(FTD)
    for j = 1:length(FTD(:,1))
        nosetail(FTD(j,1):FTD(j,2)) = [];
        frames(FTD(j,1):FTD(j,2)) = [];
        wormlength(FTD(j,1):FTD(j,2)) = [];
        q = length(FTD(j,1):FTD(j,2));
    end
end
% if N
%     FTF = FTF-q;
% end
for j = 1:length(FTF(:,1))
    for i = FTF(j,1):FTF(j,2)
        nosetail{i}(:,:) = flipud(nosetail{i}(:,:));
    end
end

NumFrames = length(nosetail);   %store the number of frames in the dataset
FRAMES = 1:NumFrames;           %set the range to be analyzed

TotSeg = 20;                    %how many segments to divide worm by
SegNum = TotSeg - 2;            %number of segments in the analysis
SegLngth = 100/TotSeg;          %the point length of each segment along worm spline.
SEGMENTS = 1:SegNum;            %set the range of segments for the analysis

%% ComputeAngles
AllHeadAngles = zeros(SEGMENTS(end),NumFrames);

for i = FRAMES
    for k = SEGMENTS
        v1_x = nosetail{i}(SegLngth*(k-1)+1,1);
        v1_y = nosetail{i}(SegLngth*(k-1)+1,2);
        v2_x = nosetail{i}(SegLngth*(k-1)+SegLngth,1);
        v2_y = nosetail{i}(SegLngth*(k-1)+SegLngth,2);
        v = [v2_x-v1_x, v2_y - v1_y]; 

        u1_x = nosetail{i}(SegLngth*(k+1)+1,1);
        u1_y = nosetail{i}(SegLngth*(k+1)+1,2);
        u2_x = nosetail{i}(SegLngth*(k+1)+SegLngth,1);
        u2_y = nosetail{i}(SegLngth*(k+1)+SegLngth,2);
        u = [u2_x-u1_x, u2_y - u1_y]; 

        cos_alpha = dot(v,u) / sqrt(dot(v,v) * dot(u,u)); 
        alpha = acos(cos_alpha);

        w = cross([u ,0], [v, 0]);
        alpha = alpha*(sign(w(3)));

        if ~isreal(alpha)
            alpha = real(alpha)
        end
        
        AllHeadAngles(k,i) = alpha;
    end;
end;

for j = SEGMENTS
    for i = 2:FRAMES(end)-2
        if isfinite(AllHeadAngles(j,i))
        else
            x1 = AllHeadAngles(j,i-2:i-1);
            x2 = AllHeadAngles(j,i+1:i+2);
            AllHeadAngles(j,i) = (mean(x1(isfinite(x1))) + mean(x2(isfinite(x2))))/2;
        end
    end
end

clear('v1_x','v1_y','v2_x','v2_y','v','u1_x','u1_y','u2_x','u2_y','u','cos_alpha','alpha','w','FTF');

%% Real Time axis,and Sliding Window Pointer arrays
RT = (0:frames(end));
FRAMESRT = 1:length(RT);
RT_L = ismember(RT,frames); %Positions along 10 hour 10Hz time where there is data.
RT_L = logical(RT_L);

clear('i','j','k','x','x1','x2');

%% Construct a nosetail and AllHeadAngles structures that corresponds to real time
nosetailRT = cell(100,2);
AllHeadAnglesRT = zeros(SegNum,length(RT));
for i = 1:length(RT)
    if RT_L(i) == 1
        x1 = find(frames==i-1);
        nosetailRT{i} = nosetail{x1};
        AllHeadAnglesRT(:,i) = AllHeadAngles(:,x1);
    else
        nosetailRT{i}(1:100,1) = NaN;
        nosetailRT{i}(1:100,2) = NaN;
        AllHeadAnglesRT(:,i) = NaN;
    end
end
%% Smooth the data
AHA_S05 = {};

for j = SEGMENTS
    AHA_S05{1}(j,FRAMES) = smooth(AllHeadAngles(j,FRAMES),3);
end
for i = 2:5
    for j = SEGMENTS
        AHA_S05{i}(j,FRAMES) = smooth(AHA_S05{i-1}(j,FRAMES),3);
    end
end

%% Worm Lengths
% WormLengths = wormlength;
for i = 1:NumFrames
    x = [];
    x = [0;cumsum(sqrt(diff(nosetail{i}(:,1)).^2 + diff(nosetail{i}(:,2)).^2))];
    WormLengths(i) = x(end);
end

WormLengthRT = zeros(1,length(RT));
for i = 1:length(RT)
    if RT_L(i) == 1
        x1 = find(frames==i-1);
        WormLengthRT(i) = WormLengths(x1);
    else
        WormLengthRT(i) = NaN;
    end
end

%% calculate the derivative of the angles
dAngles = [];
for i = SEGMENTS
    dAngles(i,:)= abs(diff(AHA_S05{5}(i,FRAMES)));
end

dAnglesRT = zeros(SegNum,length(RT_L))*NaN;
for i = 1:length(RT)-1
    if RT_L(i) == 1
        x1 = find(frames==i-1);
        if size(dAngles,2) < x1
            dAnglesRT(:,i) = NaN;
        else
            dAnglesRT(:,i) = dAngles(:,x1);
        end
    else
        dAnglesRT(:,i) = NaN;
    end
end
%% Center of Mass
%Calculate the center of mass
CoM = zeros(length(RT_L),2)*NaN;
for i = 1:length(RT_L)
    if RT_L(i) == 1
        x1 = find(frames==i-1);
        MeanX = mean(nosetail{x1}(:,1));
        MeanY = mean(nosetail{x1}(:,2));
        CoM(i,:) = [MeanX MeanY];
    end
end
%%
save([Worm ' AngleData'],'Worm','frames','frameRate','RT_L','nosetailRT','AllHeadAnglesRT','AHA_S05',...
    'WormLengthRT','CoM','dAngles','dAnglesRT');




