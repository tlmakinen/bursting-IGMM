function [ts_I]=numericaltwostate(kon,koff,tPol,ms,DT,T,ncell,isround)
% From the two state model kon, koff, with interphase duration T.
    % Input:
        % kon
        % koff
        % tPol: polymerase loading time (default: 6s)
        % ms: number of binding sites along the gene segment
        % DT: oupput sampling interval
        % T: interphase duration
        % ncell: number of cells
% Generate the time series:
        % ts_onoff: time series: ON-OFF window
        % ts_I: time series: spot intensity
    if ~exist('isround','var')
        isround=false;
    end
    dt=6;
    %% Generate the on-off window:
    tick(1)=0;
    val(1)=0;
    poltick=[];
    cnt=1;
    while tick(cnt)<=(T*ncell)
        % Generate the on-off window
        cnt=cnt+1;
        tick(cnt)=tick(cnt-1)+exprnd(1/kon);
        val(cnt)=1;
        cnt=cnt+1;
        toff=exprnd(1/koff);
        tick(cnt)=tick(cnt-1)+toff;
        val(cnt)=0;
        % Generate the polymerase arrival time
        poltick=[poltick tick(cnt-1)+addPol(toff,tPol)];
    end
    
    % Polymerase arrival time signal
    signal2=interp1(tick,val,0:DT:max(tick),'previous');
    % Intensity signal:
    L=[];
    ke=25;              % Elongation rate
    dlen=dt*ke;         % Polymerase size
    for i=1:numel(ms)/dlen
        L(i)=sum(ms((dlen*(i-1)+1):(dlen*i)))/dlen;
    end
    if i<numel(ms)/dlen
        L(i+1)=sum(ms((dlen*(i)+1):end))/dlen;
    end
    L=[L 0 0];
    signal3=conv(signal2,L);
    %% Break down the long time series into smaller ones
    ts_I={};
    for i=1:ncell
        ts_I{i}=signal3((i-1)*T/dt+1:DT/dt:i*T/dt);
    end