function test_on_pso

    f = @(x)(x^3-3*x-1);

    S = 100;
    dim = 1;

    lb = -5;
    ub = 5;

    particles = rand(S,dim).*(ub-lb)+lb;
    pi = particles;
    fpi = arrayfun(f,pi);
    [bestg,bestid] = min(fpi);
    g = pi(bestid,:);

    vmax = ub-lb;
    vi = rand(S,dim).*(2*vmax)-vmax;

    %% tuning parameters
    w = 0.1;phip = 0.3;phig = 0.5;
    
    maxiter = 1e3;
    lastbest = -1e6
    iter = 1;
    tol = 1e-6;
    while iter < maxiter
        for i = 1:size(particles,1)
            for j = 1:dim
                rp = rand(1);
                rg = rand(1);
                vi(i,j) = w*vi(i,j)+phip*rp*(pi(i,j)-particles(i,j))+phig*rg*(g-particles(i,j));
            end
            particles(i,:) = particles(i,:) + vi(i,:);
            fx = f(particles(i,:));
            if fx < fpi(i)
                fpi(i) = fx;
                pi(i,:) = particles(i,:);
                if fx < bestg
                    bestg = fx;
                    g = particles(i,:);
                end
            end
        end
        
        iter = iter + 1;
    end
    
    g

end


