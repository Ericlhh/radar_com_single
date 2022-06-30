function n = theoretic(m,s,d)
    syms t;
    syms u;
    f = s*(m+t^2)+(d-s)*int((u-t)^2*u^(m-1)*exp(-u^2/2)/(2^(m/2-1)*gamma(m/2)),u,t,inf);
    g = diff(f,t);
    t1 = solve(g);
    n = s*(m+t1^2)+(d-s)*int((u-t1)^2*u^(m-1)*exp(-u^2/2)/(2^(m/2-1)*gamma(m/2)),u,t1,inf);