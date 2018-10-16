function test_on_tlbo
    % enable long format
    format long;
    
    %% initialize
    lb = [0 0 0 0 0 0 0 0 0 0 0 0 0];
    ub = [1 1 1 1 1 1 1 1 1 100 100 100 1];
    
    nsize = 100;
    msize = 13;
    maxiter = 1200;
    
    %% initialize particles
    for i = 1:msize
        pn(:,i) = rand(nsize,1).*(ub(i)-lb(i))+lb(i);
    end
    
    %% maintain diversity of the sampels
    [pn] = remove_duplicate(pn,ub,lb);
    %% compute scores
    fs = evaluatefunc(pn);
    
    [newfs, indices] = sort(fs,'ascend');
    newpn = pn(indices,:);
    pn = newpn;
    fs = newfs;
    
    validfs = fs(fs<inf);
    disp(strcat('avg_score: ', num2str(mean(validfs))));
    
    for iter=1:maxiter
        meanpn = mean(pn);
        [bestscore, bestid] = max(fs);
        teacher = pn(bestid,:);
        
        for j = 1:size(pn,1)
            TF = floor(1+rand(1));
            newpn1 = pn(j,:) + rand(1,msize).*(teacher-meanpn.*TF);
            
            %% constaint
            idp = newpn1 > ub;
            idn = newpn1 < lb;
            newpn1(idp) = ub(idp);
            newpn1(idn) = lb(idn);
            %% new f
            fnew = func(newpn1);
            
            if fnew < fs(j)
                fs(j) = fnew;
                pn(j,:) = newpn1;
            end
            
            %% interlearning
            ids = 1:size(pn,1);
            ids(j) = [];
            pairid = randi(size(ids));
            pair = pn(ids(pairid),:);
            
            if func(pair) < fs(j)
                newpn2 = pn(j,:) + rand(1,msize).*(pair-pn(j,:));
            else
                newpn2 = pn(j,:) + rand(1,msize).*(-pair+pn(j,:));
            end
            fnew2 = func(newpn2);
            %% constaint
            idp = newpn2 > ub;
            idn = newpn2 < lb;
            newpn2(idp) = ub(idp);
            newpn2(idn) = lb(idn);
            if fnew2 < fs(j)
                fs(j) = fnew2;
                pn(j,:) = newpn2;
            end
        end
        
        fs = evaluatefunc(pn);
        
        [newfs, indices] = sort(fs,'ascend');
        newpn = pn(indices,:);
        pn = newpn;
        fs = newfs;
        
        if rand(1) < 1
            [pn] = remove_duplicate(pn,ub,lb);
        end
        fs = evaluatefunc(pn);
        [newfs, indices] = sort(fs,'ascend');
        newpn = pn(indices,:);
        pn = newpn;
        fs = newfs;
        disp(pn);
    end
    
end


function fs = evaluatefunc(pn)
    fs = zeros(size(pn,1),1);
    for i = 1:size(pn,1)
        fs(i,1) = func(pn(i,:));
    end
end

function f = func(x)
    format long;
    
    s1 = 5*sum(x(1:4));
    s2 = -5*sum(x(1:4).^2);
    s3 = -sum(x(5:13));
    
    f1 = s1 + s2 + s3;
    
    p1=x(1);
    p2=x(2);
    p3=x(3);
    p4=x(4);
    p5=x(5);
    p6=x(6);
    p7=x(7);
    p8=x(8);
    p9=x(9);
    p10=x(10);
    p11=x(11);
    p12=x(12);
    p13=x(13);
    
    t1=2*p1+2*p2+p10+p11-10;
    t2=2*p1+2*p3+p10+p12-10;
    t3=2*p2+2*p3+p11+p12-10;
    t4=-8*p1+p10;
    t5=-8*p2+p11;
    t6=-8*p3+p12;

    t7=-2*p4-p5+p10;
    t8=-2*p6-p7+p11;
    t9=-2*p8-p9+p12;
    nc=9;
    g1(1)=t1;
    g1(2)=t2;
    g1(3)=t3;
    g1(4)=t4;
    g1(5)=t5;
    g1(6)=t6;
    g1(7)=t7;
    g1(8)=t8;
    g1(9)=t9;
    fun=0;
    cov=0;  
    for io=1:nc
        if g1(io)>0
            fun=fun+g1(io)^2;
            cov=cov+1;
        end
    end
    f = (f1) + (1e20*fun) + (1e15*cov);
end

function [pn] = remove_duplicate(pn,ub,lb)
   for i=1:size(pn,1)
       m1 = sort(pn(i,:));
       for j = i+1:size(pn,1)
           m2 = sort(pn(j,:));
           if isequal(m1,m2)
               switchbit = randi(size(pn,2));
               pn(j,switchbit) = lb(switchbit)+(ub(switchbit)-lb(switchbit)).*rand();
           end
       end
   end
end