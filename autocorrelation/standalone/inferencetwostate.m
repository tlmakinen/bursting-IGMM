function []=inferencetwostate(T,ncell,pon_range,isrounded)
    sumk_range=[0.002 0.01 0.025:0.1 0.11:0.01:0.15];
    tPol=6; % RNAP loading time
    ke=25;  % Elongation rate
    DT=6;
    itrmax=5; % Iteration time
    ts_onoff={};
    ts_I={};
    % Create storage
    meanpon=[];
    errpon=[];
    % Load the binding site number along the gene sequence
    load('therightL');  % Load ms the MS2 configuration
        L=[];
        ke=25;              % Elongation rate
        sizePol=tPol*ke;    % Polymerase size
        for i=1:numel(ms)/sizePol
            L(i)=sum(ms((sizePol*(i-1)+1):(sizePol*i)))/sizePol;
        end
        if i<numel(ms)/sizePol
            L(i+1)=sum(ms((sizePol*(i)+1):end))/sizePol;
        end
        L=[L 0 0];
    %% Simulation
    konhat={};
    koffhat={};
    meanpon={};
    varpon={};
    errpon={};
    cnt1=0;
    for sumk=sumk_range
        cnt1=cnt1+1;
        cnt2=0;
        for pon=pon_range
            cnt2=cnt2+1;
            % Calculate kon and koff
            kon=sumk*pon;
            koff=sumk*(1-pon);
            % Begin the inference
            meanpon_=zeros(1,itrmax);
            sumkhat_=zeros(1,itrmax);
            parfor itr=1:itrmax
                itr
                % Create the samples
                [ts_I]=numericaltwostate(kon,koff,tPol,ms,DT,T,ncell,isrounded);
                meanpon_(itr)=pon;
                % Calculate autocorrelation function
                cnt=0;
                auto=[];
                for i=1:numel(ts_I)
                    take=ts_I{i};
                    cnt=cnt+1;
                    % Take the traces
                    if (sum(take>0)>=1)
                        % Subtract the traces to the mean value
                        take=take-mean(take);
                        autoco=[];
                        for j=1:numel(take)
                            autoco(j)=mean(take(1:(numel(take)-j+1)).*take(j:numel(take)));
                        end
                    else
                        autoco=NaN*ones(1,length(take));
                    end
                    % Normalization of autocorrelation function
                    auto(cnt,:)=autoco;
                end
                theauto=nanmean(auto,1);
                thevar=nanvar(auto,1);
                p_true=[log((1-sumk)*tPol) 1];
                % Perform the inference:
                fun1 = @(p) errorcorrectedtwostate4(p,theauto',pon,DT,0,numel(theauto),L,thevar',tPol);
                %problem1=opti('fun',fun1,'x0',[log(1-(kon+koff)*tPol) 1],'bounds',[-4 0.5],[-0.001 2],'options',opts);
                %[p_min,fval] = solve(problem1);
                %opts=optimset('TolFun',1e-1,'MaxFunEvals',400);
                [p_min,fvalhat]=fminsearchbnd(fun1,p_true,[-4 0.5],[-0.001 2]);
                fval=errorcorrectedtwostate4(p_true,theauto',pon,DT,0,numel(theauto),L,thevar',tPol);
                if fval<fvalhat
                    p_min=p_true;
                end
                if (p_min(1)<=-4+1e-3)|(p_min(1)>=-0.001)
                    p_min(1)=NaN;
                end
                sumkhat_(itr)=(1-exp(p_min(1)))/tPol;
                [(1-exp(p_min(1)))/tPol sumk]
            end
            meanpon{cnt1,cnt2}=meanpon_;
            sumkhat{cnt1,cnt2}=sumkhat_;
        end
    end
    %% Save the results
    save(['Results/twostate_' num2str(pon_range) '_T' num2str(T) '_ncell' num2str(ncell) '.mat']);

    