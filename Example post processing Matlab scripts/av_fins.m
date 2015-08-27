function [mns,sds,Ns] = av_fins(x)

mns = [];sds = [];Ns = [];
for j = 1:size(x,2)
    indx = isfinite(x(:,j));
    mns = [mns, mean(x(indx,j))];
    sds = [sds, std(x(indx,j))];
    Ns = [Ns, sum(indx)];
end;
return;