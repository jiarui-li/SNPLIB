function unrelated_ind = FindUnrelatedGroup(obj, threshold)
%FINDUNRELATEDGROUP Using a SNP list with HWE, maf filtered and 
%   Detailed explanation goes here
if nargin<2
    threshold = 0.044;
end
king = obj.CalcKING();
R = king>threshold;
r = sum(R);
I = find(r>0);
R = R(I,I);
r = r(I);
[~,idx] = max(r);
ind = [];
while true
    R(idx,:) = 0;
    R(:,idx) = 0;
    ind = [ind,idx];
    if isempty(ind)
        break;
    end
    r = sum(R);
    if all(r==0)
        break;
    end
    [~,idx] = max(r);    
end
unrelated_ind = setdiff(1:obj.nSamples, I(ind));
end

