function [ANOVAtbl, PHtbl] = UnbalancedANOVA2(Data, Grp1, Grp2)

% Data ... single vector contains data values
% Grp1&2 ... group ids for each factor

%% ANOVA
[~, ANOVAtbl, stats] = anovan(Data, {Grp1, Grp2}, 'model',2,'varnames',{'Area','Stimulation'});

%% post-hoc test
c = multcompare(stats, 'Ctype', 'hsd', 'Dimension', [1 2], 'Display', 'on');
c

%% extracting p values 
nGrps = max(c(:,2));
PHtbl = zeros(nGrps, nGrps);
for i = 1:size(c,1)
    PHtbl(c(i,1), c(i,2)) = c(i, end);
end

%%





