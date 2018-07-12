clear;
fileID = fopen('Original.txt');
nfileID = fopen('FinalSolution.txt');
if all(fgetl(fileID)==-1)
    disp('Empty File')
end

pt =[];
while(~feof(fileID))
    line = fgetl(fileID);
    tmp = strsplit(line);
    tmp = [str2double(tmp{1}) str2double(tmp{2}) str2double(tmp{3})];
    pt = [pt;tmp];
end

nwpt = [];
while (~feof(nfileID))
    line = fgetl(nfileID);
    tmp = strsplit(line);
    tmp = [str2double(tmp{1}) str2double(tmp{2}) str2double(tmp{3})];
    nwpt = [nwpt;tmp];
end
scatter3(pt(:,1),pt(:,2),pt(:,3));
hold on;
scatter3(nwpt(:,1),nwpt(:,2),nwpt(:,3));