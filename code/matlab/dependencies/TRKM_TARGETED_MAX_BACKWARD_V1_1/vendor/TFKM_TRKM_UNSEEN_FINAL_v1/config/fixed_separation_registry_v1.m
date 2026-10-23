function [delta,rho] = fixed_separation_registry_v1(pn,r,c,arm,difficulty)
% Calibrated anchors. Intermediate structures use deterministic scaling.
arm = string(arm); difficulty = string(difficulty);
levels = ["difficult","moderate","easy","perfect"];
k = find(levels==difficulty,1); assert(~isempty(k),'Unknown difficulty.');
% Base anchors for (r,c)=(2,3).
if pn <= 0.15
    if arm=="TRKM_favorable", d=[2 3 4 10]; rho0=.5; else, d=[3 4 5 10]; rho0=2; end
elseif pn <= 0.4
    if arm=="TRKM_favorable", d=[1.8 2.7 5.4 18]; rho0=.35; else, d=[4 6 12 40]; rho0=6; end
elseif pn <= 0.75
    if arm=="TRKM_favorable", d=[3 6 15 40]; rho0=.25; else, d=[5 10 24 64]; rho0=8; end
elseif pn <= 1.5
    if arm=="TRKM_favorable", d=[5 12 28 80]; rho0=.18; else, d=[8 18 40 96]; rho0=16; end
else
    if arm=="TRKM_favorable", d=[8 20 40 128]; rho0=.10; else, d=[12 28 64 128]; rho0=40; end
end
% Mild structure scaling, preserving the calibrated ordering.
scale = sqrt((r+c)/(2+3));
delta = d(k)*scale;
rho = rho0;
end
