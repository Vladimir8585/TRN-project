function [AL, AR, FL, FR] = VV_LoadFiles(ParentPath)

% load AL3
cd([ParentPath, '\AL_S']);
Finfo = dir('*_AL.mat');
tmp = load(Finfo(1).name);
AL = tmp.Vecto_AL3;

% load AR3
cd([ParentPath, '\AR_S']);
Finfo = dir('*_AR.mat');
tmp = load(Finfo(1).name);
AR = tmp.Vecto_AR3;

% load FL3
cd([ParentPath, '\FL_S']);
Finfo = dir('*_FL.mat');
tmp = load(Finfo(1).name);
FL = tmp.Vecto_FL3;

% load FR3
cd([ParentPath, '\FR_S']);
Finfo = dir('*_FR.mat');
tmp = load(Finfo(1).name);
FR = tmp.Vecto_FR3;

end

