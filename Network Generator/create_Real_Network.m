try 
d.unload
catch ERR
end 
fclose all;clear class;
close all;clear all;clc;
addpath(genpath(pwd))

%% Load Network:
[inpname,dispname] = chooseNetwork([]);
d=epanet(inpname);
nn=d.getNodeCount;

%% Compare demand patterns:
% demands1Res = demandGenerator('Residential'); % put demand vector here
% demands2Com = demandGenerator('Commercial'); % put demand vector here
% pat = d.getPattern;
% pat1Res = pat(1,:);
% pat2Com = pat(2,:);
% 
% [min(pat1Res) mean(pat1Res) max(pat1Res)]
% [min(demands1Res) mean(demands1Res) max(demands1Res)]
% x=6:6:288;
% figure
% plot(x,pat1Res(1:288/6))
% hold all
% plot(demands1Res(1:288))
% ylim([0.1,1.5])
% [min(pat2Com) mean(pat2Com) max(pat2Com)]
% [min(demands2Com) mean(demands2Com) max(demands2Com)]
% x=6:6:288;
% figure
% plot(x,pat2Com(1:288/6))
% hold all
% plot(demands2Com(1:288))
% ylim([0.1,1.5])

%% Assign patterns to network:
%%% Get existing patterns:
%%%     d.getNodeDemandPatternIndex{categoryIndex}(nodeIndex)
demInd1 = double(d.getNodeDemandPatternIndex{1}(:)); 
demInd2 = double(d.getNodeDemandPatternIndex{2}(:));
baseDemInd1 = double(d.getNodeBaseDemands{1}(:));
baseDemInd2 = double(d.getNodeBaseDemands{2}(:));

%%% Delete two(2) old patterns:
d.deletePattern(1)
d.deletePattern(1) % pattern number 2 become 1 after deletion of 1

%%% Create new patterns:
resDemPatNum = 50;
for i = 1:resDemPatNum
    demands = demandGenerator('Residential'); % put residential demand vector here
    resDemPatInd(i)=d.addPattern(['P-Res-',num2str(i)],demands);
    disp(['Creating pattern Residential ',num2str(i)])
end
comDemPatNum = 50;
for i = 1:comDemPatNum
    demands = demandGenerator('Commercial'); % put commercial demand vector here
    comDemPatInd(i)=d.addPattern(['P-Com-',num2str(i)],demands);
    disp(['Creating pattern Commercial ',num2str(i)])
end

%%% Residential = 2
%%% Commercial = 1

% Demand Category 1:
for i=1:length(demInd1)
    disp(['Indexing pattern ',num2str(i),' category 1'])
    if baseDemInd1(i)==0
        continue 
    elseif demInd1(i)==2 % Residential
        ind = round(rand*resDemPatNum+0.5);
        demInd1(i)=resDemPatInd(ind);
    elseif demInd1(i)==1 % Commercial
        ind = round(rand*comDemPatNum+0.5);
        demInd1(i)=comDemPatInd(ind);
    else
        error('unknown demand pattern')
    end
end

% Demand Category 2:
for i=1:length(demInd2)
    disp(['Indexing pattern ',num2str(i),' category 1'])
    if baseDemInd2(i)==0
        continue 
    elseif demInd2(i)==2 % Residential
        ind = round(rand*resDemPatNum+0.5);
        demInd2(i)=resDemPatInd(ind);
    elseif demInd2(i)==1 % Commercial
        ind = round(rand*comDemPatNum+0.5);
        demInd2(i)=comDemPatInd(ind);
    else
        error('unknown demand pattern')
    end
end

%%% Assign new patterns:
for i=1:nn
    disp(['Assigning pattern ',num2str(i),' out of ',num2str(nn)])
%%%     d.setNodeDemandPatternIndex(nodeIndex, categoryIndex, patternIndex)     
    d.setNodeDemandPatternIndex(i, 1, demInd1(i))
    d.setNodeDemandPatternIndex(i, 2, demInd2(i))
end

%% Randomize base demands:
uncB = 0.1; % base demand uncertainty

baseDemInd1 = double(d.getNodeBaseDemands{1}(:));
bl1 = baseDemInd1-uncB*baseDemInd1;
bu1 = baseDemInd1+uncB*baseDemInd1;
baseDemNew1 = bl1+rand(size(baseDemInd1)).*(bu1-bl1);

baseDemInd2 = double(d.getNodeBaseDemands{2}(:));
bl2 = baseDemInd2-uncB*baseDemInd2;
bu2 = baseDemInd2+uncB*baseDemInd2;
baseDemNew2 = bl2+rand(size(baseDemInd2)).*(bu2-bl2);

for i=1:nn
    d.setNodeBaseDemands(i, 1, baseDemNew1(i))
    d.setNodeBaseDemands(i, 2, baseDemNew2(i))
end

%% Randomize pipe lenghts, roughness and diameters, tank diameter:
uncR=0.1; % roughness uncertainty
uncL=0.02; % length uncertainty
uncD=0.02; % diameter uncertainty
uncT=0.05; % tank diameter uncertainty

R=double(d.getLinkRoughnessCoeff);
Rl=R-uncR*R;
Ru=R+uncR*R;
R=Rl+rand(size(R)).*(Ru-Rl);
d.setLinkRoughnessCoeff(R);

L=double(d.getLinkLength);
Ll=L-uncL*L;
Lu=L+uncL*L;
L=Ll+rand(size(L)).*(Lu-Ll);
d.setLinkLength(L);

D=double(d.getLinkDiameter);
Dl=D-2*uncD*D; % In the real network the diameters are always smaller than the model
Du=D;
D=Dl+rand(size(D)).*(Du-Dl);
d.setLinkDiameter(D);

T=double(d.getNodeTankDiameter);
Tl=T-uncT*T;
Tu=T+uncT*T;
T=Tl+rand(size(T)).*(Tu-Tl);
d.setNodeTankDiameter(T)

%% Structural uncertainty:
closedLinkInd = d.getLinkIndex({'p37','p251'});
for i = closedLinkInd
    d.setLinkInitialStatus(i,0)
end

%% Correct times:
d.setTimeReportingStep(300)
d.setTimeHydraulicStep(300)
d.setTimePatternStep(300)

%% Save new input file:
newInpname = ['networks\',dispname,'_Real.inp'];
d.saveInputFile(newInpname);
disp('NETWORK READY!')
d.unload;


