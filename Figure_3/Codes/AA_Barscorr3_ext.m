






%% For 30 sec trace
%%


Nk5w = 0;
Nk5 = 0;


K55w = [];
K5w = [];

K55 = [];
K5 = [];



for a = 1:N3_rows-2
    if Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...  % On signal which lasts for 30 sec
            sum(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000))) ~= 0 && ...
            sum(Scorr((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) ~= 0; % after On signal 80 seconds state change
        Nk5w = Nk5w+1; 
        
        K5w(Nk5w,1) = find(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000)),1,'first'); % time from On signal till state change
        K55w(Nk5w,1) = find(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000)),1,'last');
        
        K5w(Nk5w,2)= (a-1)/2 +1; % order of stimulation
        K55w(Nk5w,2)= (a-1)/2 +1; % order of stimulation
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...  %
            sum(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000))) ~= 0 && ...
            sum(Scorr((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+60000))) == 0; % Same - only collecting numbers of traces with no state change
        Nk5 = Nk5+1;  
        
        K5(Nk5,1) = 61000; % define seconds for trace with non changed state
        K55(Nk5w,1) = find(Scorr((Scorr_sig3(a,1)-30000):(Scorr_sig3(a,1)-10000)),1,'last');
        
        K5(Nk5,2)= (a-1)/2 +1; % order of the trace
        K55(Nk5,2)= (a-1)/2 +1;
    end
end


%% For 60 seconds
%%


Nk5w = 0;
Nk5 = 0;


K55w = [];
K5w = [];

K55 = [];
K5 = [];



for a = 1:N3_rows-2
    if Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...  % On signal which lasts for 30 sec
            sum(Scorr((Scorr_sig3(a,1)-40000):(Scorr_sig3(a,1)-10000))) ~= 0 && ...
            sum(Scorr((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+80000))) ~= 0; % after On signal 80 seconds state change
        Nk5w = Nk5w+1; 
        
        K5w(Nk5w,1) = find(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+80000)),1,'first'); % time from On signal till state change
        K55w(Nk5w,1) = find(Scorr((Scorr_sig3(a,1)-40000):(Scorr_sig3(a,1)-10000)),1,'last');
        
        K5w(Nk5w,2)= (a-1)/2 +1; % order of stimulation
        K55w(Nk5w,2)= (a-1)/2 +1; % order of stimulation
        
    elseif Scorr_sig3(a,2) == 1 && (Scorr_sig3(a+1,3))>=29999 && (Scorr_sig3(a+1,3))< 35000 &&...  %
            sum(Scorr((Scorr_sig3(a,1)-40000):(Scorr_sig3(a,1)-10000))) ~= 0 && ...
            sum(Scorr((Scorr_sig3(a,1)-10000):(Scorr_sig3(a,1)))) == 0 && ...
            sum(Scorr((Scorr_sig3(a,1)):(Scorr_sig3(a,1)+80000))) == 0; % Same - only collecting numbers of traces with no state change
        Nk5 = Nk5+1;  
        
        K5(Nk5,1) = 61000; % define seconds for trace with non changed state
        K55(Nk5w,1) = find(Scorr((Scorr_sig3(a,1)-40000):(Scorr_sig3(a,1)-10000)),1,'last');
        
        K5(Nk5,2)= (a-1)/2 +1; % order of the trace
        K55(Nk5,2)= (a-1)/2 +1;
    end
end


