function [newU,order]=expm_4(V,D,invV,s,imax)
    % Calculate the matrix exponential of U*(-imax)..U*imax
    tmp1=s-imax:s+imax;
    [tmp2,order]=sort(abs(tmp1));
    order(order)=1:numel(order);
    newU=arrayfun(@(x) V*diag(exp(D*abs(x)))*invV,tmp2,'UniformOutput',false);
    
    
    
    
