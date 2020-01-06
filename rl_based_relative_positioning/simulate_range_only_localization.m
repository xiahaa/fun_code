function simulate_range_only_localization
close all;
%% simulation parameters
numRob = 3;%the number of robots involved in this simulation
stdMotion = 0.01;% covariance for motion data and range data
stdRange = 0.1;
dim = 2;% 2D space
sampleTime = 0.01;% 1s
freq = 100;%100 hz
totalTime = 10;
t = 1/freq;
%% firstly, generate simulated data
lastSampleTime = -1.1;
initTime = 0;
currTime = initTime;
psi = zeros(numRob,1);%% heading
vel = zeros(numRob,1);%% velocity

p1 = [0 0 0];
p2 = [3 0 -pi/2];
p3 = [1 2 pi/2];

pos{1} = p1;
pos{2} = p2;
pos{3} = p3;

k = 1;

%% range computation
ranges{1} = [0 norm(p2(1:2)-p1(1:2)) norm(p3(1:2)-p1(1:2))];
ranges{2} = [norm(p1(1:2)-p2(1:2)) 0 norm(p3(1:2)-p2(1:2))];
ranges{3} = [norm(p1(1:2)-p3(1:2)) norm(p2(1:2)-p3(1:2)) 0];

while currTime < totalTime
    currTime = currTime + t;
    if (currTime - lastSampleTime) >= sampleTime
        lastSampleTime = currTime;
        %% random sample vel and heading
        psinew = -pi + 2*pi.*rand(3,1);
        delpsi = psinew - psi;
        id = abs(delpsi) < pi/4;
        psi(id) = psinew(id);
        psi(~id) = psi(~id) + pi/4.*sign(delpsi(~id));
        vel = 5.*rand(3,1);
    end
    
    %% sample new pos
    for j = 1:numRob
        p = [pos{j}(end,1) + t.*vel(j) * cos(psi(j)) ...
             pos{j}(end,2) + t.*vel(j) * sin(psi(j)) ...
             psi(j)];
        if j == 1
            p = [0 0 0];
        end
        pos{j}(k+1,:) = p;
    end
    
    %% sample new range
    for j = 1:numRob
        ranges_i = zeros(1,numRob);
        pj = pos{j}(k+1,1:2);
        for i = 1:numRob
            pii = pos{i}(k+1,1:2);
            d = norm(pii-pj);
            ranges_i(1,i) = d;
        end
        ranges{j}(k+1,:) = ranges_i;
    end
    k = k+1;
end

%% motion data
stdvar = stdMotion;
motions{1} = pos{1}(2:end,1:2) - pos{1}(1:end-1,1:2) + stdvar.*randn(size(pos{1},1)-1,2);
motions{2} = pos{2}(2:end,1:2) - pos{2}(1:end-1,1:2) + stdvar.*randn(size(pos{1},1)-1,2);
motions{3} = pos{3}(2:end,1:2) - pos{3}(1:end-1,1:2) + stdvar.*randn(size(pos{1},1)-1,2);

rangesn = ranges;
stdvar = stdRange;
rangesn{1}(:,2) = rangesn{1}(:,2) + stdvar.*randn(size(rangesn{1},1),1);
rangesn{1}(:,3) = rangesn{1}(:,3) + stdvar.*randn(size(rangesn{1},1),1);
rangesn{2}(:,3) = rangesn{2}(:,3) + stdvar.*randn(size(rangesn{1},1),1);

rangesn{2}(:,1) = rangesn{1}(:,2);
rangesn{3}(:,1) = rangesn{1}(:,3);
rangesn{3}(:,2) = rangesn{2}(:,3);

posr{1} = p1;
posr{2} = p2;
posr{3} = p3; 

for k = 1:1:size(motions{1},1)
    for j = 1:1:3
        posr{j}(k+1,1:2) = sum(motions{j}(1:k,:),1) + posr{j}(1,1:2);
    end
end

figure
plot(pos{1}(:,1),pos{1}(:,2),'r-');hold on;
plot(posr{1}(:,1),posr{1}(:,2),'r-o');hold on;
plot(pos{2}(:,1),pos{2}(:,2),'g-');hold on;
plot(posr{2}(:,1),posr{2}(:,2),'g-o');hold on;
plot(pos{3}(:,1),pos{3}(:,2),'b-');hold on;
plot(posr{3}(:,1),posr{3}(:,2),'b-o');hold on;

figure
subplot(3,1,1);
plot(ranges{1}(:,2),'r-');hold on;
plot(rangesn{2}(:,1),'g-');hold on;
subplot(3,1,2);
plot(ranges{1}(:,3),'r-');hold on;
plot(rangesn{3}(:,1),'g-');hold on;
subplot(3,1,3);
plot(ranges{2}(:,3),'r-');hold on;
plot(rangesn{3}(:,2),'g-');hold on;

%% relative positioning development
posrNTU{1} = p1;
posrNTU{2} = [0 0 0];
posrNTU{3} = [0 0 0];

kk = 1;
k = 1;
tralaterationNum = 100;

doDLT = 0;
doFineOpt = 0;
while k < size(rangesn{1},1)
    %% collect data
    ide = k + tralaterationNum;
    if ide > size(rangesn{1},1)
        ide = size(rangesn{1},1);
        tralaterationNum = ide - k;
    end
    
    if doDLT == 0
        doDLT = 1;
        doFineOpt = 1;
        %% compute 2
        d0 = rangesn{1}(k,2);
        cummOdom2 = zeros(tralaterationNum,2);
        for j = 1:1:tralaterationNum
            cummOdom2(j,:) = sum(motions{2}(k:(k+j-1),:),1);
        end
    
        b = rangesn{1}(k+1:ide,2).^2 - repmat(d0,ide-k, 1).^2 - cummOdom2(:,1).^2 - cummOdom2(:,2).^2;
        A = 2.*[cummOdom2(:,1) cummOdom2(:,2)];
        xx2 = inv(A'*A)*A'*b;
    
        %% compute 3
        d0 = rangesn{1}(k,3);
        cummOdom3 = zeros(tralaterationNum,2);
        for j = 1:1:tralaterationNum
            cummOdom3(j,:) = sum(motions{3}(k:(k+j-1),:),1);
        end
    
        b = rangesn{1}(k+1:ide,3).^2 - repmat(d0,ide-k, 1).^2 - cummOdom3(:,1).^2 - cummOdom3(:,2).^2;
        A = 2.*[cummOdom3(:,1) cummOdom3(:,2)];
        xx3 = inv(A'*A)*A'*b;
    else
        doFineOpt = 0;
        xx2 = posrNTU{2}(end,1:2) + motions{2}(k,:);xx2 = xx2';
        xx3 = posrNTU{3}(end,1:2) + motions{3}(k,:);xx3 = xx3';
        cummOdom2 = zeros(tralaterationNum,2);
        for j = 1:1:tralaterationNum
            cummOdom2(j,:) = sum(motions{2}(k:(k+j-1),:),1);
        end
        cummOdom3 = zeros(tralaterationNum,2);
        for j = 1:1:tralaterationNum
            cummOdom3(j,:) = sum(motions{3}(k:(k+j-1),:),1);
        end
    end
        
    %% nonlinear optimization
    sols = [xx2;xx3];
%     cummOdom = zeros(tralaterationNum,2);
%     for j = 1:1:tralaterationNum
%         cummOdom(j,:) = sum(motions{2}(k:(k+j-1),:),1);
%     end
    motionsopt{1} = cummOdom2;
%     cummOdom = zeros(tralaterationNum,2);
%     for j = 1:1:tralaterationNum
%         cummOdom(j,:) = sum(motions{3}(k:(k+j-1),:),1);
%     end
    motionsopt{2} = cummOdom3;
    
    rangesopt = [rangesn{1}(k:ide,2);rangesn{1}(k:ide,3);rangesn{2}(k:ide,3)];
    if doFineOpt == 1
        sols = nonlinOpt(sols, motionsopt, rangesopt, tralaterationNum, 1e-10, 1e-10, 1000);
    else
        sols = nonlinOpt(sols, motionsopt, rangesopt, tralaterationNum, 1e-2, 1e-2, 100);
    end
    xx2 = sols(1:2);
    xx3 = sols(3:4);
    posrNTU{1}(kk+1,:) = [0 0 0];
    posrNTU{2}(kk+1,:) = [xx2' 0];
    posrNTU{3}(kk+1,:) = [xx3' 0];
%     for ii = 1:1:tralaterationNum
%         posrNTU{1}(kk+1,:) = [0 0 0];
%         p = xx2' + cummOdom2(ii,1:2);
%         posrNTU{2}(kk+1,:) = [p 0];
%         p = xx3' + cummOdom3(ii,1:2);
%         posrNTU{3}(kk+1,:) = [p 0];
%         kk = kk + 1;
%     end
    
    k = k + tralaterationNum;
    kk = kk + 1;
    break;
end

figure
plot(pos{1}(:,1),pos{1}(:,2),'r-');hold on;
% plot(posrNTU{1}(1,1),posrNTU{1}(1,2),'r+');
hold on;plot(posrNTU{1}(2,1),posrNTU{1}(2,2),'r*');hold on;

plot(pos{2}(:,1),pos{2}(:,2),'g-');hold on;
% plot(posrNTU{2}(1,1),posrNTU{2}(1,2),'g+');
hold on;plot(posrNTU{2}(2,1),posrNTU{2}(2,2),'g*');hold on;

plot(pos{3}(:,1),pos{3}(:,2),'b-');hold on;
% plot(posrNTU{3}(1,1),posrNTU{3}(1,2),'b+');
hold on;plot(posrNTU{3}(2,1),posrNTU{3}(2,2),'b*');hold on;

%% now, we have the initialized position, we do EKF
kk = size(posrNTU{2},1);
k = 1;
covP = eye(2,2);
covR = stdRange*stdRange;
covQ = stdMotion * stdMotion * eye(2);
    while k < size(rangesn{1},1)
        pr = posrNTU{2}(kk+k-1,1:2);
        %% prediction
        pp = pr + motions{2}(k,:);pp = pp';
        covP = covP + covQ;
        %% update
        z = rangesn{1}(k,2);
        C = [1/norm(pp)*pp(1) 1/norm(pp)*pp(2)];
        K = covP*C'*inv(C*covP*C'+covR);
        pc = pp + K * (z - norm(pp));
        posrNTU{2}(kk+k,1:2) = pc';
        covP = (eye(2) - K*C)*covP;
        k = k + 1;
    end

kk = size(posrNTU{3},1);
k = 1;
covP = 0.1*eye(2,2);
covR = stdRange*stdRange;
covQ = stdMotion * stdMotion * eye(2);
while k < size(rangesn{1},1)
    pr = posrNTU{3}(kk+k-1,1:2);
    %% prediction
    pp = pr + motions{3}(k,:);pp = pp';
    covP = covP + covQ;
    %% update
    z = rangesn{1}(k,3);
    C = [1/norm(pp)*pp(1) 1/norm(pp)*pp(2)];
    K = covP*C'*inv(C*covP*C'+covR);
    pc = pp + K * (z - norm(pp));
    posrNTU{3}(kk+k,1:2) = pc';
    covP = (eye(2) - K*C)*covP;
    k = k + 1;
end

figure
plot(pos{2}(:,1),pos{2}(:,2),'g-','LineWidth',2);hold on;
hold on;plot(posrNTU{2}(2:end,1),posrNTU{2}(2:end,2),'g-o');hold on;
plot(pos{3}(:,1),pos{3}(:,2),'b-','LineWidth',2);hold on;
hold on;plot(posrNTU{3}(2:end,1),posrNTU{3}(2:end,2),'b-o');hold on;
legend('ground truth of drone 1','estimation of drone 1', 'ground truth of drone 2','estimation of drone 2');
grid on;
title('Simulation of range only relative localization');
end

function [sols] = nonlinOpt(sols, motions, ranges, motionCnt, epsilon1, epsilon2, cnt)
    k = 0; v = 2; x = sols;
    J = calcJacobian(x, motions, motionCnt, ranges);
    f = evaluationf(x, motions, ranges, motionCnt);
    A = J'*J; g = J'*f;
    
    if (max(g)<=epsilon1) return; end
    
    miu = 1e-6 * max(diag(A));
    
    iter = 1;
    
    while iter < cnt
        cost = f'*f;
        disp(strcat('iter:',num2str(iter),'; cost:',num2str(cost)));
        iter = iter + 1;
        hlm = (A+miu*eye(size(A,1)))\(-g);
        if norm(hlm) <= epsilon2*(norm(x)+epsilon2) break; end
        xnew = x + hlm;
        fnew = evaluationf(xnew, motions, ranges, motionCnt);
        rho = (f'*f - fnew'*fnew) * 2 / (hlm'*(miu.*hlm-g));
        if rho > 0
            x = xnew;
            J = calcJacobian(x, motions, motionCnt, ranges);
            f = evaluationf(x, motions, ranges, motionCnt);
            A = J'*J; g = J'*f;
            if (max(g)<=epsilon1) break; end
            miu = miu * max(1/3,1-(2*rho-1)^3); v = 2;
        else
            miu = miu * v; v = v * 2;
        end
    end
    
    sols = x;
end

function  J = calcJacobian(sols, motions, motionCnt, ranges)
    n1 = numel(sols)/2;
    J = zeros(numel(ranges),numel(sols));
    
    %% compute all frame positions
    for i = 1:1:n1
        p0 = sols(i*2-1:i*2)';
        pm = repmat(p0,motionCnt,1) + motions{i}(:,1:2);
        pos{i} = [p0;pm];
    end
    
    %% i to 0
    k = 1;
    for i = 1:1:n1
        p = pos{i}(:,1:2);
        rp = vecnorm(p,2,2);
        for j = 1:1:numel(rp)
            J(k,i*2-1) = -1/rp(j)*p(j,1);
            J(k,i*2) = -1/rp(j)*p(j,2);
            k = k + 1;
        end
    end
    %% i to j
    for i = 1:1:n1
        p1 = pos{i}(:,1:2);
        for j = i+1:1:n1
            p2 = pos{j}(:,1:2);
            pd = p1 - p2;
            rpd = vecnorm(pd,2,2);
            for l = 1:1:numel(rpd)
                J(k,i*2-1) = -1/rpd(l)*pd(l,1);
                J(k,i*2) = -1/rpd(l)*pd(l,2);
                
                J(k,j*2-1) = 1/rpd(l)*pd(l,1);
                J(k,j*2) = 1/rpd(l)*pd(l,2);
                
                k = k + 1;
            end
        end
    end
end

function feval = evaluationf(sols, motions, ranges, motionCnt)
    n1 = numel(sols)/2;
    feval = [];
    
    %% compute all frame positions
    for i = 1:1:n1
        p0 = sols(i*2-1:i*2)';
        pm = repmat(p0,motionCnt,1) + motions{i}(:,1:2);
        pos{i} = [p0;pm];
    end
        
    k = 1;
    %% i to 0
    for i = 1:1:n1
        p = pos{i}(:,1:2);
        rp = vecnorm(p,2, 2);
        ff = ranges(k:k+numel(rp)-1) - rp;
        feval = [feval;ff];
        k = k + numel(rp);
    end
    %% i to j
    for i = 1:1:n1
        p1 = pos{i}(:,1:2);
        for j = i+1:1:n1
            p2 = pos{j}(:,1:2);
            pd = p1 - p2;
            rpd = vecnorm(pd,2, 2);
            ff = ranges(k:k+numel(rpd)-1) - rpd;
            feval = [feval;ff];
            k = k + numel(rpd);
        end
    end
end


