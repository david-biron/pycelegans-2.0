function [xs,ys] = ReconstructMidline(Angs)

data_Ang = Angs;
% data_Ang = [pi/2, pi/2, pi/2, pi/2, pi/2];
xs = 0:(length(data_Ang)+2);
ys = zeros(size(xs));

i = 1;
if Angs(i)>0
    theta = data_Ang(i)/2;
    R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
else
    theta = -data_Ang(i)/2;
    R = [cos(theta) sin(theta) ; -sin(theta) cos(theta)];
end;
xs = xs-xs(i+1);
ys = ys-ys(i+1);
u = [xs(i+2:end); ys(i+2:end)];
rotu = R * u;
xs(i+2:end) = rotu(1,:);
ys(i+2:end) = rotu(2,:);

theta2 = atan(ys(i+2)/xs(i+2));
u = [xs; ys];
if theta2 > 0
    rotu = [cos(theta2) sin(theta2) ; -sin(theta2) cos(theta2)] * u;
else
    %     theta2 = -theta2;
    rotu = [cos(theta2) sin(theta2) ; -sin(theta2) cos(theta2)] * u;
end
xs = rotu(1,:);
ys = rotu(2,:);
%%%%%%%%%%%%%%%

xs = xs-xs(i+2);
ys = ys-ys(i+2);
u = [xs(i+3:end); ys(i+3:end)];
rotu = R * u;
xs(i+3:end) = rotu(1,:);
ys(i+3:end) = rotu(2,:);


%%%%%%%%%%%%%%%
for i = 2:length(data_Ang)
    
    theta_a = atan2(ys(i+2), xs(i+2));
    
    
    theta2 = atan2(ys(i+2), xs(i+2));
    u = [xs; ys];
    if theta2 > 0
        rotu = [cos(theta2) sin(theta2) ; -sin(theta2) cos(theta2)] * u;
    else
        %         theta2 = -theta2;
        rotu = [cos(theta2) sin(theta2) ; -sin(theta2) cos(theta2)] * u;
    end
    xs = rotu(1,:);
    ys = rotu(2,:);
    
    
    if Angs(i)>0
        theta = data_Ang(i)-theta_a;
        R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
    else
        theta = -data_Ang(i)+theta_a;
        R = [cos(theta) sin(theta) ; -sin(theta) cos(theta)];
    end;
    
    
    xs = xs-xs(i+2);
    ys = ys-ys(i+2);
    u = [xs(i+3:end); ys(i+3:end)];
    rotu = R * u;
    xs(i+3:end) = rotu(1,:);
    ys(i+3:end) = rotu(2,:);
    
    
end;

%     subplot(2,1,2);
% figure;
% hold on;
% plot(xs, ys,'b');
% plot(xs(1),ys(1),'m.','markersize',10)
% axis equal;
return;