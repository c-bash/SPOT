numiters = 10000;
speed1 = zeros(numiters,1);
speed2 = zeros(numiters,1);

for inter = 1:numiters
    A = random_cov(135);
    tic
        CF_ANS = chol(A,'lower');
    speed1(inter) = toc;
    
    tic
        CF_TEST = my_cf(A);
    speed2(inter) = toc;
    
end

disp(mean(speed1))
disp(mean(speed2))
disp(mean(speed2)/mean(speed1))

function L = my_cf(A)
    L = zeros(size(A));
    
    iters = size(A,1);
    for rank=1:iters
        lambda = sqrt(A(rank,rank));
        l = A(rank+1:end,rank)/lambda;
        
        L(rank,rank) = lambda;
        L(rank+1:end,rank) = l;
        
        if rank < iters
            A(rank+1:end,rank+1:end) = A(rank+1:end,rank+1:end) - l*l.';
        end
    end
    
end

function L = my_cf_bidiag(A,n) 
    % n is the size of the Lii blocks
    L = zeros(size(A));
    
    L(1:n,1:n) = my_cf(A(1:n,1:n));
    
    iters = size(A,1)/n;
    for row=2:iters % i.e., going by row
        ii = (row-2)*n;
        ip1 = (row-1)*n;
%         disp(ip1+n)
        temp = forwardsubstitution(L(ii+1:ii+n,ii+1:ii+n),A(ii+1:ii+n,ip1+1:ip1+n));
        L(ip1+1:ip1+n,ii+1:ii+n) = temp.';
        
        L(ip1+1:ip1+n,ip1+1:ip1+n) = my_cf(A(ip1+1:ip1+n,ip1+1:ip1+n) - L(ip1+1:ip1+n,ii+1:ii+n)*temp);
    end
    
end


function A = random_cov(n)

    Q = randn(n,n);

    eigen_mean = 13; 
    % can be made anything, even zero 
    % used to shift the mode of the distribution

    A = Q' * diag(abs(eigen_mean+randn(n,1))) * Q;

end 