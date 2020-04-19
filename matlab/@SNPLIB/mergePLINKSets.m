function mergePLINKSets(obj,bfile_list)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
bimSpec = '%s%s%f%f%s%s';
famSpec = '%s%s%s%s%f%s';
obj.CHR = [];
obj.RSID = [];
obj.POS = [];
obj.ALT = [];
obj.REF = [];
bfile = bfile_list(1);
filename = strcat(bfile,'.fam');
fileID = fopen(filename,'r');
dataArray = textscan(fileID, famSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
obj.FID = dataArray{1};
obj.IID = dataArray{2};
obj.PID = dataArray{3};
obj.MID = dataArray{4};
obj.Sex = dataArray{5};
L = ceil(obj.nSamples/4);
for i = 1:length(bfile_list)
    bfile = bfile_list(i);
    filename = strcat(bfile,'.bim');
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, bimSpec, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    chr = dataArray{1};
    CHR = str2double(chr);
    CHR(strcmp(chr,'X')) = 23;
    CHR(strcmp(chr,'Y')) = 24;
    CHR(strcmp(chr,'XY')) = 25;
    CHR(strcmp(chr,'MT')) = 26;
    obj.CHR = [obj.CHR;CHR];
    obj.RSID = [obj.RSID;dataArray{2}];
    obj.POS = [obj.POS;dataArray{4}];
    obj.ALT = [obj.ALT;upper(dataArray{5})];
    obj.REF = [obj.REF;upper(dataArray{6})];    
    fid = fopen(strcat(bfile, '.bed'),'rb');
    bin = fread(fid,inf,'uint8=>uint8');
    fclose(fid);
    nSNPs = length(CHR);    
    obj.GENO = [obj.GENO,reshape(bin(4:end),[L nSNPs])];
end
end

