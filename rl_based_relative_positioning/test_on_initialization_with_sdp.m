function test_on_initialization_with_sdp
    clear all;
    close all;
    clc;
    
    dim = 2;
    %%
    nx = -10;
    ny = -10;
    px = 10;
    py = 10;
    num1 = 1000;
    bx = linspace(nx,px,num1)';
    by = linspace(ny,py,num1)';
    %% workspace
    boders = [[bx ones(numel(bx),1).*ny];[bx ones(numel(bx),1).*py];[ones(numel(by),1).*nx by];[ones(numel(by),1).*px by]];
    figure(1);
    plot(boders(1:num1,1),boders(1:num1,2),'r-','LineWidth',3);hold on; grid on;
    plot(boders(num1+1:num1*2,1),boders(num1+1:num1*2,2),'r-','LineWidth',3);
    plot(boders(num1*2+1:num1*3,1),boders(num1*2+1:num1*3,2),'r-','LineWidth',3);
    plot(boders(num1*3+1:num1*4,1),boders(num1*3+1:num1*4,2),'r-','LineWidth',3);

    %% anchor
    num_of_anchors = 5;
    pos_anchors = zeros(num_of_anchors,dim);
    for i = 1:num_of_anchors
        sx = rand()*(px-nx)+nx;
        sy = rand()*(py-ny)+ny;
        pos_anchors(i,:) = [sx sy];
    end
    figure(1);
    plot(pos_anchors(:,1),pos_anchors(:,2),'bs','MarkerSize',10);
    
    %% node
    num_of_nodes = 1;
    pos_nodes = zeros(num_of_nodes,dim);
    for i = 1:num_of_nodes
        while(1)
            sx = rand()*(px-nx)+nx;
            sy = rand()*(py-ny)+ny;

            ex = sx - pos_anchors(:,1);
            ey = sy - pos_anchors(:,2);
            e_err = sqrt(ex.^2+ey.^2);
            redo = ~isempty(find(e_err < 1));
            if redo == true
                continue;
            else
                break;
            end
        end
        pos_nodes(i,:) = [sx sy];
    end
    figure(1);
    plot(pos_nodes(:,1),pos_nodes(:,2),'rd','MarkerSize',10);
    
    %% compute ranges
    ranges = zeros(num_of_nodes, num_of_anchors);
    range_noise_var = 1.0;
    for i = 1:num_of_nodes
        ranges(i,:) = sqrt((pos_anchors(:,1)-pos_nodes(i,1)).^2+(pos_anchors(:,2)-pos_nodes(i,2)).^2)';
        ranges(i,:) = ranges(i,:) + randn(1,numel(ranges(i,:))).*sqrt(range_noise_var);
    end
    
    for i = 1:num_of_nodes
        for j = 1:num_of_anchors
            lx = [pos_nodes(i,1), pos_anchors(j,1)];
            ly = [pos_nodes(i,2), pos_anchors(j,2)];
            plot(lx,ly,'g-','LineWidth',2);
            text(mean(lx), mean(ly), num2str(ranges(i,j)), 'Interpreter', 'latex');
        end
    end
    
    color = jet(128);
    
    %% without initialization solve with steepest descent
    pos_nodes_est1 = pos_nodes;
    for i = 1:num_of_nodes
        r = ranges(i,:)';
        tic
        x1 = solveBySteepestDescentMethod(pos_anchors, r); 
        toc
        pos_nodes_est1(i,:) = x1';
    end
    figure(1);
    h1 = plot(pos_nodes_est1(:,1),pos_nodes_est1(:,2),'Marker','*','MarkerSize',10, 'Color', color(1,:));
    
    pos_nodes_est2 = pos_nodes;
    for i = 1:num_of_nodes
        r = ranges(i,:)';
        tic
        x1 = solveByGaussNewton(pos_anchors, r); 
        toc
        pos_nodes_est2(i,:) = x1';
    end
    figure(1);
    h2 = plot(pos_nodes_est2(:,1),pos_nodes_est2(:,2),'Marker','*','MarkerSize',10, 'Color', color(20,:));
    
    pos_nodes_est3 = pos_nodes;
    for i = 1:num_of_nodes
        r = ranges(i,:)';
        tic
        x1 = solveByLevenbergMarquardt(pos_anchors, r); 
        toc
        pos_nodes_est3(i,:) = x1';
    end
    figure(1);
    h3 = plot(pos_nodes_est3(:,1),pos_nodes_est3(:,2),'Marker','o','MarkerSize',10, 'Color', color(40,:));
    
    xb = [-5,5];
    yb = [-5,5];
    
    pos_nodes_est4 = pos_nodes;
    for i = 1:num_of_nodes
        r = ranges(i,:)';
        tic
        xd0 = dlt(pos_anchors, r); 
%         x1 = sdp(pos_anchors, r, xb, yb, x0); 
        x1 = solveByLevenbergMarquardt(pos_anchors, r, xd0); 
        toc
        pos_nodes_est4(i,:) = x1';
    end
    figure(1);
    h4 = plot(pos_nodes_est4(:,1),pos_nodes_est4(:,2),'Marker','+','MarkerSize',20, 'Color', color(100,:));
    
    pos_nodes_est5 = pos_nodes;
    for i = 1:num_of_nodes
        r = ranges(i,:)';
%         x0 = dlt(pos_anchors, r); 
        xs0 = sdp(pos_anchors, r, xb, yb, [0 0]'); 
%         x1 = solveByLevenbergMarquardt(pos_anchors, r, xs0'); 
        pos_nodes_est5(i,:) = xs0';
    end
    figure(1);
    h5 = plot(pos_nodes_est5(:,1),pos_nodes_est5(:,2),'Marker','o','MarkerSize',20, 'Color', color(120,:));
    
    legend([h1,h2,h3,h4,h5],'SteepestDescent','GaussNewtom', 'LevenbergMarquardt','DLT','SDP','Interpreter','latex');
   
    bx = linspace(-20,20,1000)';
    by = linspace(-20,20,1000)';
    
    F = zeros(numel(bx),numel(by));
    for i = 1:numel(by)
        ys = by(i);
        for j = 1:size(pos_anchors,1)
            f = (ranges(1,j) - sqrt((pos_anchors(j,1)-bx(:)).^2+(pos_anchors(j,2)-repmat(ys,size(bx,1),size(bx,2))).^2));
            f2(:,j) = f.^2;
        end
        F(i,:) = sum(f2,2);
    end
    
    figure(2);
    mesh(bx,by,F);
    hold on;
    grid on;
    
    [v,id] = min(F(:));
    [my, mx] = ind2sub(size(F),id);
    plot3(bx(mx),by(my),F(my,mx),'ms','MarkerSize',10);
    
    cost = ranges(1,:)' - sqrt((pos_anchors(:,1)-pos_nodes(1)).^2+(pos_anchors(:,2)-pos_nodes(2)).^2);
    cost = cost'*cost;
    plot3(pos_nodes(:,1),pos_nodes(:,2),cost,'rd','MarkerSize',10);
%     view(3);
    
end

function grad = calcGradient(x, pa)
    numerator1 = pa(:,1) - x(1);
    numerator2 = pa(:,2) - x(2);
    denominator = sqrt(numerator1.^2+numerator2.^2);
    grad = [(numerator1./denominator)';(numerator2./denominator)'];
end

function jaco = calcJacobian(x, pa)
    grad = calcGradient(x, pa);
    jaco = grad';
end

function f = calcf(x,pa,r)
    f = r - sqrt((pa(:,1)-x(1)).^2+(pa(:,2)-x(2)).^2);
end

function a = lineSearch(x0, h, pa, r)
    a = 1;
    beta = 0.8;
    f0 = calcf(x0,pa,r);
    F0 = f0'*f0;
    k = 1;
    while 1
        xnew = x0 + a.*h;
        f = calcf(xnew,pa,r);
        F = f'*f;
        if F < F0 * 0.95
            break;
        else
            a = a * beta;
            k = k + 1;
        end
        if k > 1e5
            a = 0;
            break;
        end
    end
end

function x1 = solveBySteepestDescentMethod(varargin)
    if nargin() == 3
        x0 = varargin{3};
    else
        x0 = [0;0];
    end
    pa = varargin{1}; 
    r = varargin{2};
    tol = 1e6;
    maxIter = 1e6;
    k = 1;
    oldtol = 1e-6;
    while abs(tol - oldtol) > 1e-10 && k < maxIter
        f = calcf(x0,pa,r);
        jaco = calcJacobian(x0, pa);
        Fdot = (f'*jaco)';
        h = -Fdot;
        a = lineSearch(x0, h, pa, r);
        x0 = x0 + h * a;
        k = k + 1;
        oldtol = tol;
        tol = f'*f;
    end
    x1 = x0;
end

function x1 = solveByGaussNewton(varargin)
    if nargin() == 3
        x0 = varargin{3};
    else
        x0 = [0;0];
    end
    pa = varargin{1}; 
    r = varargin{2};
    tol = 1e6;
    maxIter = 1e6;
    k = 1;
    oldtol = 1e-6;
    while abs(tol - oldtol) > 1e-10 && k < maxIter
        f = calcf(x0,pa,r);
        jaco = calcJacobian(x0, pa);
        
        H = jaco'*jaco;
        b = -jaco'*f;
        h = H\b;
        
        a = lineSearch(x0, h, pa, r);
        x0 = x0 + h * a;
        k = k + 1;
        oldtol = tol;
        tol = f'*f;
    end
    x1 = x0;
end

function x1 = solveByLevenbergMarquardt(varargin)
    if nargin() == 3
        x0 = varargin{3};
    else
        x0 = [0;0];
    end    
    pa = varargin{1}; 
    r = varargin{2};
    k = 0; v = 2;
    J = calcJacobian(x0, pa);
    f = calcf(x0,pa,r);
    A = J'*J; g = J'*f;
    epsilon1 = 1e-10;
    epsilon2 = 1e-10;
    if (abs(max(g))<=epsilon1) 
        x1 = x0; 
        return; 
    end
    
    miu = 1e-6 * max(diag(A));
    iter = 1;
    maxIter = 1e6;
    while iter < maxIter
        cost = f'*f;
        disp(strcat('iter:',num2str(iter),'; cost:',num2str(cost)));
        iter = iter + 1;
        hlm = (A+miu*eye(size(A,1)))\(-g);
        if norm(hlm) <= epsilon2*(norm(x0)+epsilon2) break; end
        xnew = x0 + hlm;
        fnew = calcf(xnew,pa,r);
        rho = (f'*f - fnew'*fnew) * 2 / (hlm'*(miu.*hlm-g));
        if rho > 0
            x0 = xnew;
            J = calcJacobian(x0, pa);
            f = calcf(x0,pa,r);
            A = J'*J; g = J'*f;
            if (abs(max(g))<=epsilon1) 
                break; 
            end
            miu = miu * max(1/3,1-(2*rho-1)^3); v = 2;
        else
            miu = miu * v; v = v * 2;
        end
    end
    
    x1 = x0;
end

function x0 = dlt(pa, r)
    n1 = size(pa,1);
    n2 = round(n1*(n1-1)*0.5);
    A = zeros(n2,2);
    b = zeros(n2,1);
    k = 1;
    for i = 1:n1
        for j = i+1:n1
            di = r(i);
            dj = r(j);
            xi = pa(i,1);yi = pa(i,2);
            xj = pa(j,1);yj = pa(j,2);
            b(k) = di^2-dj^2-(xi^2-xj^2)-(yi^2-yj^2);
            A(k,:) = [2*xj-2*xi, 2*yj-2*yi];
            k = k + 1;
        end
    end
    x0 = A\b;
end

function x0 = sdp(pa, r, xb, yb, xinit)
    %% use sdr to solve this problem
    % calc A
    n1 = size(pa,1);
    A = zeros(5,5);
    dim = 2;
    
    import mosek.fusion.*
    M = Model('cqo1');
%     e = M.variable('e', n1, Domain.unbounded());
    c = M.variable('c', n1, Domain.greaterThan(0.0));
    x = M.variable('x', dim, Domain.unbounded());
    t = M.variable('t', 1, Domain.greaterThan(0.0));
    % create the aliases
    a1 = Expr.sub(Var.repeat(x.index(1),n1),pa(:,1));
    b1 = Expr.sub(Var.repeat(x.index(2),n1),pa(:,2));
    e = Expr.sub(r,c);
%     d = Expr.mul(x.index(2),5);
    M.constraint('c1',Expr.hstack(c,a1,b1),Domain.inQCone());
    M.constraint('c2',Expr.hstack(t,Expr.constTerm(1, 1.0/2.0),Expr.transpose(e)),Domain.inRotatedQCone());
    v1 = ones(n1,1);
    M.objective("obj", ObjectiveSense.Minimize, t);
    M.solve();
    
    disp(toString(M.getPrimalSolutionStatus));
    disp(num2str(t.level()));
    x0 = x.level();
    x0 = x0';
  
    %% randomization
%     ndim = size(X,1);
%     mu = zeros(ndim,1);
%     Sigma = X;
%     [V,D]=eig(Sigma);
%     A = V*sqrt(max(D,0));
%     
%     K = 200;
%     val = 1e6;xhat = zeros(ndim,1);
%     for k = 1:K
%         n = size(mu,1);
%         z = randn(n,1);
%         xr = mu + A*z;
%         cost = xr'*A*xr;
%         if val > cost
%             val = cost;
%             xhat = xr;
%         end
%     end
%     x0 = [xhat(3)/xhat(2);xhat(5)/xhat(4)];
end

