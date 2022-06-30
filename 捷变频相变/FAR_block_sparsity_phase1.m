close all;
clear all;
clc;
M = 4;
N = 128;
%block_sparsity = 1;
tol = 1e-5;
trial = 50;
epi = 0.02;
result = zeros(N,25);
for col = 4:4:128
    for block_sparsity = 10:18
        success_count = 0;
        for loop = 1:trial
            FAR_model = zeros(N,M*N);
            %Cn = randperm(M)-1
            for n = 0 : N-1
                Cn = floor(rand()*M);
                for q = 0 : N-1
                    for p = 0:M-1
                        FAR_model(n+1,q*M+p+1) = exp(1i*2*pi*p/M*Cn+1i*2*pi*q/N*n*(1+Cn*epi));
                    end
                end
            end
            col_choose = randperm(N,col);
            FAR_model = FAR_model(col_choose,:);
            sparse_signal = zeros(M,N);
            block = randperm(N,block_sparsity);
            sparse_signal(:,block) = exp(1i*2*pi*rand(M,block_sparsity));
            y = FAR_model * sparse_signal(:);
            cvx_begin
            variable x(M,N) complex
            norm21 = 0;
            for i = 1:N
                norm21 = norm21 + norm(x(:,i));
            end
            minimize(norm21)
            subject to
            FAR_model * x(:) == y
            cvx_end
            if norm(x(:)-sparse_signal(:))<tol
                success_count = success_count+1;
            end
        end
        result(col,block_sparsity) = success_count/trial;
    end
end
save('FARblockepsilon2.mat');
