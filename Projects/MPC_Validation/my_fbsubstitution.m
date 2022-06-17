
FSTEST1a = forwardsubstitution(A,b);
FSTEST1b = forwardsubstitution2(A,b);

%%
numiters = 10000;
speed1 = zeros(numiters,1);
speed2 = zeros(numiters,1);
speed2b = zeros(numiters,1);
speed3 = zeros(numiters,1);
speed4 = zeros(numiters,1);

for inter = 1:numiters
    A = random_cov(9);
    A = tril(A);
    b = randn(9,6);
    
    tic
        FSTEST2 = A\b;
    speed1(inter) = toc;
    
    tic
        FSTEST1 = forwardsubstitution(A,b);
    speed2(inter) = toc;
    
    tic
        FSTEST1b = forwardsubstitution2(A,b);
    speed2b(inter) = toc;
    
    A = random_cov(9);
    A = triu(A);
    b = randn(9,6);

    tic
        BSTEST2 = A\b;
    speed3(inter) = toc;

    tic
        BSTEST1 = backwardsubstitution(A,b);
    speed4(inter) = toc;


end

disp([mean(speed1),mean(speed2),mean(speed2b)])
disp(mean(speed2)/mean(speed1))

disp([mean(speed3),mean(speed4)])
disp(mean(speed4)/mean(speed3))

function x = forwardsubstitution(A,b)
    % A must be square, lower triangular
    x = zeros(size(b));
    for col = 1:size(b,2)
        for row = 1:size(b,1)
            temp = 0;
            if row > 1
                for aind = 1:row-1
                    temp = temp - A(row,aind)*x(aind,col);
                end
            end
            x(row,col) = (b(row,col) + temp) / A(row,row);
        end
    end
end

function x = forwardsubstitution2(A,b)
    % A must be square, lower triangular
    x = zeros(size(b));
    for col = 1:size(b,2)
        for row = 1:size(b,1)
            temp = -A(row,:) * x(:,col);
            x(row,col) = (b(row,col) + temp) / A(row,row);
        end
    end
end

function x = backwardsubstitution(A,b)
    % A must be square, upper triangular
    x = zeros(size(b));
    n = size(b,1);
    for col = 1:size(b,2)
        for row = size(b,1):-1:1
            temp = 0;
            if row < n
                for aind = n:-1:row+1
                    temp = temp - A(row,aind)*x(aind,col);
                end
            end
            x(row,col) = (b(row,col) + temp) / A(row,row);
        end
    end
end

function A = random_cov(n)

    Q = randn(n,n);

    eigen_mean = 5; 
    % can be made anything, even zero 
    % used to shift the mode of the distribution

    A = Q' * diag(abs(eigen_mean+randn(n,1))) * Q;

end 