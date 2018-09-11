function outputM = callResponse(filePath, startFrame, endFrame)
fileFolder=fullfile(filePath);
dirOutput=dir(fullfile(fileFolder,'*.xlsx'));
fileNames={dirOutput.name}';
fprintf('Hi! I detect %d file(s) under the folder %s/.\n', length(fileNames), filePath)

% Calculation
t1 = startFrame;
t2 = endFrame;

for i = 1:length(fileNames)
    c = char(strcat(filePath,'/',fileNames(i)));
    [F(:,:,i), time, realTime] = callResponse(c, t1, t2);
end

meanF = mean(F,3);
stdF = std(F,0,3);
se = stdF/sqrt(3);

% Plot
clf

for i = 1:size(meanF,2)
    errorbar(time, meanF(:,i), se(:,i))
    hold on
end

xlabel('time/s');
ylabel('\DeltaF/F_0')

M = youandme(meanF,se);
outputM = [realTime time M];
csvwrite(strcat(filePath,'/','ROIs_meanF_se.csv'), outputM);
fprintf('\n=============================SPLIT LINE=============================\n\n')
fprintf('Congratulations! Please enjoy the script! \nIf you have any questions about it,\nfeel free to contact with the author. \ne-mail: zihua.liu666@gmail.com\n');
end


%% calResponse
function [deltaF_F0, time, realTime] = callResponse(Mfilepath, startFrame, endFrame)
tic;
% Loading data
M = xlsread(Mfilepath);
realTime = M(startFrame:end,1);
M = [M(:,2) M(:,5:end)];


% Divide M into time, ROIs and background
background = M(:,end);
ROI = M(:,2:end-1);

% ROI - backround
ROI = bsxfun(@minus,ROI,background);

% Baseline vertification
F0 = mean(ROI(startFrame:endFrame,:));

% Real F start
ROI = ROI(startFrame:end,:);
deltaF_F0 = bsxfun(@rdivide,bsxfun(@minus,ROI,F0),F0);
time = M(startFrame:end,1);

toc
end

%%youandme
function M = youandme(m1, m2)

[a, b] = size(m1);
M = zeros(a, 2*b);

j = 1;
k = 1;
for i = 1:2*b
    if mod(i,2) == 1
        M(:,i) = m1(:,j);
        j = j + 1;
    else
        M(:,i) = m2(:,k);
        k = k + 1;
    end
end
end