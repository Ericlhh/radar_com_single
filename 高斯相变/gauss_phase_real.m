N = 100;
gauss_phase_res = zeros(100,100);
for col = 1:100
    %æÿ’Û…˙≥…
    for p =1:100
        suc = 0;
        for loop = 1:50
            x1 = zeros(N,1);
            q = randperm(N,p);
            x1(q) = randn(p,1);
            fai = randn(col,N);
            b = fai*x1;
            cvx_begin quiet
                variable x(N)
                minimize( norm( x, 1 ) )
                subject to
                    fai * x == b
            cvx_end
            %disp((norm(x-x1,1)))
            if (norm(x-x1,1))<10e-5
                 suc = suc+1;
             end
        end
        gauss_phase_res(col,p)=suc/50;
    end
end
save gauss_phase_real;
