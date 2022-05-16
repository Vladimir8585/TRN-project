
%%

% finding starting and ending points

Sl_start = [];
Sl_end = [];
Nz = 1;

for z = 2: tx
    if ((Scorr(z-1)) == 1 ) && ((Scorr(z)) == 0)
        Sl_start(Nz) = z;
        Nz = Nz +1;
    elseif ((Scorr(z-1)) == 0 ) && ((Scorr(z)) == 1) || ((Scorr(z-1)) == 0 ) && ((Scorr(z)) == 2)
        Sl_end(Nz) = z;
    end
end
%

% trimming unfinished points


if   (Sl_end(1))<(Sl_start(1)) 
          Sl_end(1) = [];
end

if     (Sl_end(end))<(Sl_start(end)) 
           Sl_start(end) = [];
end

size(Sl_end)
size(Sl_start)
kr = [Sl_end',Sl_start']




% Making matrix of the starting time/ending time/difference(in sec)

P_dif = (Sl_end - Sl_start);
P_time = (P_dif/1000)';

Sc_mat = [Sl_start' Sl_end' P_dif'];

Sc_mat


%% Rem scorring and counting

Rem_d = [];
Nf = 1;

for z = 2: tx
    if ((Scorr(z-1)) == 2 ) && ((Scorr(z)) == 0)
        Rem_d(Nf) = z;
        Nf = Nf +1;     
    end
end

Rem_N  = numel(Rem_d);
SWS_N = numel(P_time);
Re_ratio = Rem_N/SWS_N;

%
REM_start = [];
REM_end = [];
Nv = 1;

for z = 2: tx
    if ((Scorr(z-1)) == 1 ) && ((Scorr(z)) == 2)
        Nv = Nv +1;
        REM_start(Nv) = z;
    elseif ((Scorr(z-1)) == 2 ) && ((Scorr(z)) == 1) || ((Scorr(z-1)) == 2 ) && ((Scorr(z)) == 0)
        REM_end(Nv) = z;
    end
end

if   (REM_end(1))<(REM_start(1)) 
          REM_end(1) = [];
end
if     (REM_end(end))<(REM_start(end)) 
           REM_start(end) = [];
end
if     (REM_end(1))==(REM_start(1)) 
           REM_start(1) = []; 
           REM_end(1) = [];     
end


R_dif = REM_end - REM_start;
R_time = (R_dif/1000)';
Re_mat = [REM_start' REM_end' R_dif'];
RM_ratio = numel(R_time)/SWS_N;
%R_time(1) = [];

%%
REM_Und = [];
REM_Off = [];
REM_10 = [];
Fk = 1;
Zk = 1;
Fz = 1;


a = find(Scorr_sig3(:,2) == 1)

z = Scorr_sig3(a,1)

for ped = 1: numel(REM_start)
    for a = 2:numel(z)
         if  REM_start(ped) >= z(a-1) && REM_start(ped) < (z(a-1)+30000)
                REM_Und(Fk) = REM_start(ped);
                Fk = Fk +1;
         elseif REM_start(ped) >= (z(a-1)+30000) && REM_start(ped) < z(a)
                REM_Off(Zk) = REM_start(ped);
                Zk = Zk +1;
         end
    end
end

UndR = numel(REM_Und);
UndOf = numel(REM_Off);

%%

P_time;
R_time;
Ferd = [Rem_N Re_ratio RM_ratio UndR UndOf];








