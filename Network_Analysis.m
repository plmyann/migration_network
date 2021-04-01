clear; clc;

dirWeight = '/MigrationNetwork/Weight/';
dirnodes = '/MigrationNetwork/NodeID/';
dirImBet = '/MigrationNetwork/ImprovedNetwork/';

pacell = readmatrix('/MigrationNetwork/PA/pa_cell.csv');
pacell = unique(pacell);
orders = readmatrix('/MigrationNetwork/orders.csv', 'OutputType', 'string');
unprotect = readmatrix('/MigrationNetwork/unprotkeynode.csv', 'OutputType', 'string');
curord = char(orders(1));
disp(curord);

outfolder = [dirImBet, curord];
if ~exist(outfolder, 'dir')
    mkdir(outfolder)
end
 
%% Weight
allweight = readmatrix([dirWeight, 'w', curord, '.csv']);
%% NodeID
allnodes = readmatrix([dirnodes, curord, '.csv']);
%% Unprotected node for current order
curunp = unprotect(unprotect(:,3)==curord,:);
curunp = str2double(curunp(:,2));

parpool('local');

%% Betweenness by adding each node
IPart = 1; % Change IPart from 1 to 10
IStart = (IPart-1)*4 + 1;
if IPart == 10
    IEnd = size(curunp,1);
else
    IEnd = (IPart)*4;
end

parfor i = IStart:IEnd
    disp(int2str(i));
    
    newcell = [pacell; curunp(1:i)];
    
    % Select nodes
    weight = allweight(ismember(allweight(:,4), newcell) & ismember(allweight(:,5), newcell), :);

    weight(:,6) = 1./weight(:,6);
    curGraph = graph(weight(:,4), weight(:,5), weight(:,6));
    bet = centrality(curGraph, 'betweenness', 'Cost', curGraph.Edges.Weight);
    bet(:,2) = 1:length(bet);
    
    %% Normalize betweenness
    curcell = newcell(ismember(newcell, allnodes(:,1)));
    bet = bet(ismember(bet(:,2), curcell), :);
    temp1 = curcell(~ismember(curcell, bet(:,2)));
    if ~isempty(temp1)
        temp1(:,2) = 0;
        bet = [[temp1(:,2), temp1(:,1)]; bet];
    end
    npa = size(curcell,1);
    bet(:,1) = 2*bet(:,1)./((npa-2)*(npa-1));
    bet = sortrows(bet, 1, 'descend');

    writematrix(bet, [outfolder, '/Bet', int2str(i), '.csv']);
end
