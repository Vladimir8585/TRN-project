






%% For 30 sec trace
%%


Nk5w = 0;
Nk5 = 0;
Nk5nw = 0;
Nk5n = 0;


K5w = [];
K5 = [];
K5nw = [];
K5n = [];


for a = 1:N3_rows-2
    if Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...  % On signal which lasts for 30 sec
            sum(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000))) ~= 0 && ...
            sum(Scorr((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) ~= 0; % after On signal 80 seconds state change
        Nk5w = Nk5w+1; 
        
        K5w(Nk5w,1) = find(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000)),1,'first'); % time from On signal till state change
        K5w(Nk5w,2) = 30000 - find(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000)),1,'last');
        
        K5w(Nk5w,3)= (a-1)/2 +1; % order of stimulation
        
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...  %
            sum(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000))) ~= 0 && ...
            sum(Scorr((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) == 0; % Same - only collecting numbers of traces with no state change
        Nk5 = Nk5+1;  
        
        K5(Nk5,1) = 61000; % define seconds for trace with non changed state
        K5(Nk5,2) = 30000 - find(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000)),1,'last');
        
        K5(Nk5,3)= (a-1)/2 +1; % order of the trace
       
    end
end


for a = 1:N3_rows-2
    if Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...  % On signal which lasts for 30 sec
            sum(ScorR((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000))) ~= 0 && ...
            sum(ScorR((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(ScorR((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) ~= 0; % after On signal 80 seconds state change
        Nk5nw = Nk5nw+1; 
        
        K5nw(Nk5nw,1) = find(ScorR((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000)),1,'first'); % time from On signal till state change
        K5nw(Nk5nw,2) = 30000 - find(ScorR((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000)),1,'last');
        
        K5nw(Nk5nw,3)= (a-1)/2 +1; % order of stimulation
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...  %
            sum(ScorR((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000))) ~= 0 && ...
            sum(ScorR((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(ScorR((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) == 0; % Same - only collecting numbers of traces with no state change
        Nk5n = Nk5n+1;  
        
        K5n(Nk5n,1) = 61000; % define seconds for trace with non changed state
        K5n(Nk5n,2) = 30000 - find(ScorR((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000)),1,'last');
        
        K5n(Nk5n,3)= (a-1)/2 +1; % order of the trace
    end
end


%
if size(K5,1) == 0;
    K5(:,1) = 0;
end
if size(K5w,1) == 0;
    K5w(:,1) = 0;
end
if size(K5n,1) == 0;
    K5n(:,1) = 0;
end
if size(K5nw,1) == 0;
    K5nw(:,1) = 0;
end

%

Pik = [];

Pik(1,1) = numel(K5(:,1));
Pik(1,2) = numel(K5w(:,1));
Pik(1,3) = numel(K5n(:,1));
Pik(1,4) = numel(K5nw(:,1));

printmat(Pik, 'Scorred', 'N',  'K5 K5w K5n K5nw')






