function test_on_jaya
    clc;
    clear all;
    n = 25;%% population size;
    m = 30;%% feature size;
    maxIter = 20000;
    iter = 0;
    
    %% initialization
    lb = -30; ub = 30;
%     for i = 1:m
%         xp(:,i) = rand(n,1).*(ub-lb)+lb;
%     end
    xp = rand(n,m).*(ub-lb)+lb;
    [f]=myobj(xp);

    while iter < maxIter
        iter = iter+1;
        
        %% find best and worst
        [bestval,bestid] = min(f);
        [worstval,worstid] = max(f);
        
        xbest = xp(bestid,:);
        xworst = xp(worstid,:);
        
        %% 
        xnew = xp + rand(n,m).*(repmat(xbest,n,1)-abs(xp)) - rand(n,m).*(repmat(xworst,n,1)-abs(xp));
        
        idp = xnew(:,:) > ub;
        xnew(idp) = ub;
        idn = xnew(:,:) < lb;
        xnew(idn) = lb;
        
        [fnew]=myobj(xnew);
        ids = fnew < f;
        xp(ids,:) = xnew(ids,:);
        f(ids,:) = fnew(ids,:);
        
        disp('%%%%%%%% Final population %%%%%%%%%');    
        disp(xp);
        
    end
end

function [f]=myobj(x)
    [r,c]=size(x);
    
    for i=1:r
        x1 = x(i,1:(end-1));
        x2 = x(i,2:(end));
        y = sum(100.*(x1.^2-x2).^2+(1-x1).^2);
        z(i) = y;
    end
    f = z';
end