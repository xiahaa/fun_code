%%%%%%%%%Implement%%%%%%%%%
function [ini_fun, result_fun, result_fun_new,opti_fun, opti_fun_new] = impliment 
    format long;
    ini_fun = @implimentInitialize ;
    result_fun = @implimentresult;
    result_fun_new = @implimentresult_new;
    opti_fun = @implimentopti;
    opti_fun_new = @implimentopti_new;
return;

function [upper_limit, lower_limit, Students, select] = implimentInitialize(select)
    global lower_limit upper_limit ll ul
    Granularity = 1;
    lower_limit = ll;
    upper_limit = ul;
    ll=[0 0 0 0 0 0 0 0 0 0 0 0 0];
    ul=[1 1 1 1 1 1 1 1 1 100 100 100 1];
    %lower_limit = ll;
    upper_limit = ul;
    for popindex = 1 : select.classsize
        for k = 1 : select.var_num
            mark(k) =(ll(k))+ ((ul(k) - ll(k)) * rand);
        end
        Students(popindex).mark = mark;
    end
    select.OrderDependent = true;
return;

function [Students] = implimentresult(select,Students)
    global lower_limit upper_limit
    classsize = select.classsize;
    for popindex = 1 : classsize
        for k = 1 : select.var_num
            x(k) = Students(popindex).mark(k);
        end
        Students(popindex).result = objective(x);
    end
return

function [Studentss] = implimentresult_new(select, Students)
    global lower_limit upper_limit
    classsize = select.classsize;
    for popindex = 1 : size(Students,1)
        for k = 1 : select.var_num
            x(k) = Students(popindex,k);
        end
        Studentss(popindex) = objective(x);
    end
return

function [Students] = implimentopti(select,Students)
    global lower_limit upper_limit ll ul
    for i = 1 : select.classsize
        for k = 1 : select.var_num
            Students(i).mark(k) = max(Students(i).mark(k), ll(k));
            Students(i).mark(k) = min(Students(i).mark(k), upper_limit(k));
        end
    end
return;

function [Students] = implimentopti_new(select,Students)
    global lower_limit upper_limit ll ul
    for i = 1 : size(Students,1)
        for k = 1 : select.var_num
            Students(i,k)= max(Students(i,k), ll(k));
            Students(i,k) = min(Students(i,k), upper_limit(k));
        end
    end
return;

%%%%%%%%%%Initialize%%%%%%%%%%
function [Students, select, upper_limit,lower_limit, ini_fun, min_result, avg_result,result_fun, opti_fun, result_fun_new,opti_fun_new] = Initialize(note1, obj_fun,RandSeed)
    format long;
    select.classsize =100;
    select.var_num = 13;
    select.itration =1200;
    if ~exist('RandSeed', 'var')
        rand_gen = round(sum(100*clock));
    end
    rand('state', rand_gen);
    [ini_fun, result_fun, result_fun_new, opti_fun,opti_fun_new,] = obj_fun();
    [upper_limit, lower_limit, Students, select] = ini_fun(select);
    Students = remove_duplicate(Students,upper_limit, lower_limit);
    Students = result_fun(select, Students);
    Students = sortstudents(Students);
    average_result = result_avg(Students);
    min_result = [Students(1).result];
    avg_result = [average_result];
return;

%%%%%%%%% Objective %%%%%%%%%%
function yy=objective(x)
    format long;
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
    ZZ=(5*(p1+p2+p3+p4)-5*(p1^2+p2^2+p3^2+p4^2)- ...
        (p5+p6+p7+p8+p9+p10+p11+p12+p13));
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
    yy=(ZZ)+(1e20*fun)+(1e15*cov);

%%%%%%%%% Output %%%%%%%%%%
function out_put(note1, select, Students, within_bound, min_result)
    format long;
    if note1
        duplicate_no = 0;
        for i = 1 : select.classsize
            Mark_1 = sort(Students(i).mark);
            for k = i+1 : select.classsize
                Mark_2 = sort(Students(k).mark);
                if isequal(Mark_1, Mark_2)
                    duplicate_no = duplicate_no + 1;
                end
            end
        end
        Mark = sort(Students(1).mark);
    end
return;

%%%%%%%% remove_duplicate%%%%%%
function [Students] = remove_duplicate(Students,upper_limit, lower_limit)
    format long;
    global ll ul
    for i = 1 : length(Students)
        Mark_1 = sort(Students(i).mark);
        for k = i+1 : length(Students)
            Mark_2 = sort(Students(k).mark);
            if isequal(Mark_1, Mark_2)
                m_new = floor(1+(length(Students(k).mark)-1)*(rand));
                if length(upper_limit)==1
                    Students(k).mark(m_new) = (lower_limit + (upper_limit - lower_limit) * rand);
                else
                    Students(k).mark(m_new) = (ll(m_new) + (upper_limit(m_new) - ll(m_new)) * rand);
                end
            end
        end
    end
return;

%%%%%%% result_avg %%%%%%%%%%
function [result_av, within_bound] = result_avg(Students)
    format long;
    Result = [];
    within_bound = 0;
    for i = 1 : length(Students)
        if Students(i).result < inf
            Result = [Result Students(i).result];
            within_bound = within_bound + 1;
        end
    end
    result_av = mean(Result);
return;

%%%%%%%% run_tlbo %%%%%%%%%
function run_tlbo()
    clc;
    run=1;
    format long;
    for i=1:run
        TLBO(@impliment);
    end

%%%%%% sortstudents %%%%%%%%%%
function [Students, indices] = sortstudents(Students)
    classsize = length(Students);
    Result = zeros(1, classsize);
    indices = zeros(1, classsize);
    for i = 1 : classsize
        Result(i) = Students(i).result;
    end
    [Result, indices] = sort(Result, 2, 'ascend');
    Marks = zeros(classsize, length(Students(1).mark));
    for i = 1 : classsize
        Marks(i, :) = Students(indices(i)).mark;
    end
    for i = 1 : classsize
        Students(i).mark = Marks(i,:);
        Students(i).result = Result(i);
    end

%%%%%%%% TLBO %%%%%%%%%%%
function TLBO(obj_fun, note1, note2)
    format long;
    global ll
    if ~exist('note1', 'var')
        note1 = true;
    end
    if ~exist('note2', 'var')
        note2 = true;
    end
    [Students, select, upper_limit, lower_limit, ini_fun, min_result, avg_result, result_fun, ...
        opti_fun,result_fun_new, opti_fun_new] = Initialize(note1,obj_fun);
    elite=0;
    for COMP = 1 : select.itration
        for i = 1 : elite
            markelite(i,:) = Students(i).mark;
            resultelite(i) = Students(i).result;
        end
        for i=1:length(Students)
            cs(i,:)=Students(i).mark;
            cs_result(i)=Students(i).result;
        end
        cs;
        cs_result;
        for i = 1 : length(Students)
            mean_result=mean(cs);

            TF=round(1+rand*(1));
            [r1 r2]=sort(cs_result);
            best=cs(r2(1),:);
            for k = 1 : select.var_num
                cs_new(i,k)=cs(i,k)+((best(1,k)-TF*mean_result(k))*rand);
            end
            cs_new(i,:) = opti_fun_new(select,cs_new(i,:));
            cs_new_result(i) = result_fun_new(select,cs_new(i,:));
            if cs_new_result(i)<Students(i).result
                Students(i).mark =cs_new(i,:);
                cs(i,:)=cs_new(i,:);
                Students(i).result=cs_new_result(i);
            end
            hh=ceil(length(Students)*rand);
            while hh==i
                hh=ceil(length(Students)*rand);
            end
            if Students(i).result<Students(hh).result
                for k = 1 : select.var_num
                    cs_new(i,k)= Students(i).mark(k) + ((Students(i).mark(k) - Students(hh).mark(k))*rand);
                end
            else
                for k = 1 : select.var_num
                    cs_new(i,k)= Students(i).mark(k) + ((Students(hh).mark(k) - Students(i).mark(k))*rand);
                end
            end
            cs_new(i,:) = opti_fun_new(select,cs_new(i,:));
            cs_new_result(i) = result_fun_new(select,cs_new(i,:));
            if cs_new_result(i)<Students(i).result
                Students(i).mark =cs_new(i,:);
                cs(i,:)=cs_new(i,:);
                Students(i).result=cs_new_result(i);
            end
        end
        n = length(Students);
        Students = opti_fun(select, Students);
        Students = result_fun(select, Students);
        Students = sortstudents(Students);
        for i = 1 : elite
            Students(n-(i-1)).mark = markelite(i,:);
            Students(n-(i-1)).result = resultelite(i);
        end
        if rand<1
            Students = remove_duplicate(Students, upper_limit, lower_limit);
        end
        Students = sortstudents(Students);

        [average_result, within_bound] = result_avg(Students);
        min_result = [min_result Students(1).result];
        avg_result = [avg_result average_result];
        Mark = (Students(1).mark);
        if note1
            disp([num2str(min_result(end))]);
            disp([num2str(Mark)]);
        end
    end
    fprintf('\n %e',min_result(end));
    fprintf('\n %6.10f',Mark);
    out_put(note1, select, Students, within_bound,min_result);