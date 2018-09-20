function SNP_heatmap(mfilepath,namelistpath,chr,classlist,labelnum)
% plot and save the SNP heatmap
%
% For example, SNP_heatmap('/Users/chr4A.phase.fre.vcf.plot','/Users/classlistwheat.txt','chr4A',[6 20 5 29 21],20)
%
% This script is very friendly to those who are not familiar with MATLAB.
%
% mfilepath:     The absoulute path of your input matrix file. In the matrix file,
%                the first column contains the SNP position info and the rest coulnms are SNP info.
%
% namelistpath:  The absoulute path of sample label file. Notice! The namelist file must be a txt file.
%                The format should be, for example, 
%                                     W1
%                                     W2
%                                     W3
%                                     W31
%                                     W38
%                                     ...
%                                     D5
%                                     
%                                     
% chr:           Chromosome name, for example, 'chr4A'.
%
% classlist:     It should be a list. For example, [6 20 5 29 21]. In this case, five populations are recorded.
%
% labelnum:      It decides the label number of Y-axis

%Color map
lblue = [216 240 245]./255;
dblue = [12 20 120]./255;
rred = [204 41 32]./255;
yyellow = [251 219 130]./255;
wwhite = [1 1 1];
mymap = [dblue;lblue;yyellow;rred;wwhite];

%Load matrix file and split info
m = importdata(mfilepath,'\t');
pos = m(:,1);
matrix = m(:,2:end);

%Load name list
namelist = importdata(namelistpath);

a = size(m,1);

%Class list and initialize all matrix
N = [];
namelistnew = [];
all_four = 4 * ones(a,2);
classlist = [0,classlist];
for i = 1:length(classlist)-1
    eval(['n',num2str(i),'=matrix(:,2*(sum(classlist(1:i)))+1:2*sum(classlist(1:i+1)));']);
    eval(['nl',num2str(i),'=namelist(sum(classlist(1:i))+1:sum(classlist(1:i+1)))']);
    N = [N eval(['n',num2str(i)]) all_four];
    namelistnew = [namelistnew; eval(['nl',num2str(i)]); ' '];
end
N = N(:,1:end-2);

%Plot
h = figure;
hh = imagesc(N);
colormap(mymap);

set(gca,'XTick',1.5:2:size(N,1));
set(gca,'XTicklabel',namelistnew)

%Y axis ticks labeling
rr = 1:round(a/(labelnum-1)):a;
set(gca,'YTick',rr);
yi = pos(rr)./1000000; %Mb
 
set(gca,'YTicklabel',yi);
set(gca,'XTickLabelRotation',30);
set(gcf,'color','white','paperpositionmode','auto');

ax = ancestor(hh,'axes');
yrule = ax.YAxis;
xrule = ax.XAxis;
yrule.FontSize = 12;
xrule.FontSize = 8;

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)+11,pos(4)+5]);
set(h,'Position',[0 0 7000 3000]);
ylabel('Mb')
title(sprintf('%s',chr),'FontSize',18);

% OUTPUT FILE PATH
mfilepathI = strfind(mfilepath,'/');
mfilepathII = mfilepath(mfilepathI(1):mfilepathI(end));
ofilepath = strcat(mfilepathII,nn);

print(h,ofilepath,'-dpdf','-r0');
