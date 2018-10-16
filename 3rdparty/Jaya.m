%% Jaya algorithm
%% Rosenbrock function
function Jaya()
    clc;
    clear all;
    RUNS=30;
    runs=0;
    while(runs<RUNS)
        pop=25; % population size
        var=30; % no. of design variables
        maxFes=500000;
        maxGen=floor(maxFes/pop);
        mini=-30*ones(1,var);
        maxi=30*ones(1,var);
        [row,var]=size(mini);
        x=zeros(pop,var);
        for i=1:var
            x(:,i)=mini(i)+(maxi(i)-mini(i))*rand(pop,1);
        end
        ch=1;
        gen=0;
        f=myobj(x);
        while(gen<maxGen)
            xnew=updatepopulation(x,f);
            xnew=trimr(mini,maxi,xnew);
            fnew=myobj(xnew);
            for i=1:pop
                if(fnew(i)<f(i))
                    x(i,:)=xnew(i,:);
                    f(i)=fnew(i);
                end
            end
            disp('%%%%%%%% Final population %%%%%%%%%');    
            disp([x,f]);
            fnew=[];xnew=[];
            gen=gen+1;
            fopt(gen)=min(f);
        end
        runs=runs+1;
        [val,ind]=min(fopt);
        Fes(runs)=pop*ind;
        best(runs)=val;
    end
    bbest=min(best);
    mbest=mean(best);
    wbest=max(best);
    stdbest=std(best);
    mFes=mean(Fes);
    stdFes=std(Fes);
    fprintf('\n best=%f',bbest);
    fprintf('\n mean=%f',mbest);
    fprintf('\n worst=%f',wbest);
    fprintf('\n std. dev.=%f',stdbest);
    fprintf('\n mean function evaluations=%f',mFes);
end

function[z]=trimr(mini,maxi,x)
    [row,col]=size(x);
    for i=1:col
        x(x(:,i)<mini(i),i)=mini(i);
        x(x(:,i)>maxi(i),i)=maxi(i);
    end
    z=x;
end

function [xnew]=updatepopulation(x,f)
    [row,col]=size(x);
    [t,tindex]=min(f);
    Best=x(tindex,:);
    [w,windex]=max(f);
    worst=x(windex,:);
    xnew=zeros(row,col);
    for i=1:row
        for j=1:col
            r=rand(1,2);
            xnew(i,j)=x(i,j)+r(1)*(Best(j)-abs(x(i,j)))-r(2)*(worst(j)-abs(x(i,j)));
        end
    end
end

function [f]=myobj(x)
    [r,c]=size(x);
    for i=1:r
        y=0;
        for j=1:c-1
            y=y+(100*(x(i,j)^2-x(i,j+1))^2+(1-x(i,j))^2);
        end
        z(i)=y;
    end
    f=z';
end