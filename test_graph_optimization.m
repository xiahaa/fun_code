function test_graph_optimization
    %% read data and display
    [rows, cols, entries, rep, field, symm] = mminfo('jagmesh1.mtx');
    [A, rows, cols, entries] = mmread('jagmesh1.mtx');
    
    %% draw
    x = rand(rows,2)*10-5;

    figure(1);
    scatter(x(:,1),x(:,2),'s');
    hold on;grid on;
    links = {};
    for i = 1:rows
        a = [];
        for j = 1:cols
            if A(i,j) > 0
%                 plot([x(i,1),x(j,1)],[x(i,2),x(j,2)],'r-');
                if i ~= j
                    a = [a;j];
                end
            end
        end
        links{i} = a;
    end
    
    converged = false;
    step = 0.2;
    energy = inf;
    d = 10 / rows;
    progress = 0;
    t = 0.9;
    iter = 1;
    while converged ~= true
        disp(iter);
        iter = iter + 1;
        %% mem
        xold = x;
        energyold = energy;
        %% initialize
        energy = 0;
        f = [0 0];
        %% calc new energy
        for i = 1:rows
            rindices = 1:cols;
            rindices(i) = [];
            xr1 = (x(rindices,1) - repmat(x(i,1),cols-1,1));
            yr1 = (x(rindices,2) - repmat(x(i,2),cols-1,1));
            rnorm = sqrt(xr1.^2+yr1.^2);
            fr = rnorm-repmat(d,cols-1,1);
            vr1 = [xr1,yr1]./rnorm;
            
            %% another force
            aindices = links{i};
            xa1 = (x(aindices,1) - repmat(x(i,1),numel(aindices),1));
            ya1 = (x(aindices,2) - repmat(x(i,2),numel(aindices),1));
            anorm = sqrt(xa1.^2+ya1.^2);
            fa = anorm-repmat(d,numel(aindices),1);
            va1 = [xa1,ya1]./anorm;
            
            f = sum(fr .* vr1) + sum(fa.*va1);
            
            x(i,:) = x(i,:) + step.*f./norm(f);
            
            energy = energy + norm(f)^2;
        end
        
        %% update step
        if energy < energyold
            progress = progress + 1;
            if progress >= 5
                progress = 0;
                step = step/t;
            end
        else
            progress = 0;
            step = step * t;
        end
        
        if (sum(norm(x-xold)) < 1e-3)
            converged = true;
        end
    end
    figure(2);
    scatter(x(:,1),x(:,2),'s');
    hold on;grid on;
%     for i = 1:rows
%         for j = 1:cols
%             if A(i,j) > 0
%                 plot([x(i,1),x(j,1)],[x(i,2),x(j,2)],'r-');
%             end
%         end
%     end
end