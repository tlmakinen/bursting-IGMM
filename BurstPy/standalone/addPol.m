function [arrivalPol]=addPol(ton,tPol)
% Calculate the arrive time of Polymerase within the time window ton
arrivalPol=[];
tprev=0;
while ton-tprev>=tPol
    arrivalPol=[arrivalPol tprev+tPol];
    tprev=tprev+tPol;
end
if binornd(1,(ton-tprev)/tPol)
    arrivalPol=[arrivalPol tprev+(ton-tprev)*rand];
end