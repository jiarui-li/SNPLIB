function [index1,index2] = IntersectSNP(obj1,obj2)
% INTERSECTSNP Summary of this function goes here
%   Detailed explanation goes here
if ~(isa(obj1,'SNPLIB') && isa(obj2,'SNPLIB'))
    error('Please provide a SNPLIB class object as reference!');
end
chrs = unique(obj1.CHR);
alleles = union(union(cellstr(obj1.ALT),cellstr(obj2.ALT)),union(cellstr(obj1.REF),cellstr(obj2.REF)));
AMAP = containers.Map(alleles,1:length(alleles));
index1 = [];
index2 = [];
for i=1:length(chrs)
    ind1 = find(obj1.CHR==chrs(i));
    ind2 = find(obj2.CHR==chrs(i));
    POS1 = uint64(obj1.POS(ind1));
    POS2 = uint64(obj2.POS(ind2));
    set_1 = POS1 + ...
        uint64(bitshift(cell2mat(AMAP.values(obj1.ALT(ind1))),32,'uint64')) + ...
        uint64(bitshift(cell2mat(AMAP.values(obj1.REF(ind1))),48,'uint64'));
    set_1r = POS1 + ...
        uint64(bitshift(cell2mat(AMAP.values(obj1.REF(ind1))),32,'uint64')) + ...
        uint64(bitshift(cell2mat(AMAP.values(obj1.ALT(ind1))),48,'uint64'));
    set_2 = POS2 + ...
        uint64(bitshift(cell2mat(AMAP.values(obj2.ALT(ind2))),32,'uint64')) + ...
        uint64(bitshift(cell2mat(AMAP.values(obj2.REF(ind2))),48,'uint64'));
    set_2r = POS2 + ...
        uint64(bitshift(cell2mat(AMAP.values(obj2.REF(ind2))),32,'uint64')) + ...
        uint64(bitshift(cell2mat(AMAP.values(obj2.ALT(ind2))),48,'uint64'));
    ind_1 = find(ismember(set_1,set_2));
    ind_2 = find(ismember(set_1r,set_2));
    % if there exists CHR:POS:A:G and CHR:POS:G:A, keep the one originally
    % matching another dataset.
    [~,~,i2] = intersect(set_1(ind_1),set_1r(ind_2));
    ind_2 = setdiff(ind_2,ind_2(i2));
    index1 = [index1;ind1(union(ind_1,ind_2))];
    ind_1 = find(ismember(set_2,set_1));
    ind_2 = find(ismember(set_2r,set_1));
    [~,~,i2] = intersect(set_2(ind_1),set_2r(ind_2));
    ind_2 = setdiff(ind_2,ind_2(i2));
    index2 = [index2;ind2(union(ind_1,ind_2))];
end
end

