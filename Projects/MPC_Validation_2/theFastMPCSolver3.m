function [zout,nuout] = theFastMPCSolver3(z,nu,zt,H,C,b,P,h,epsilon1,ALPHA,BETA,maxiter)
    
    % Create a range of kappa values to decrease by. Kappa is the weight of
    % the log barrier in the objective function. A better approximation of
    % the indicator function is computed when kappa is small. Decrease
    % kappa by 10 each time the loop runs, from WangBoyd2010.
%     maxiter = 10;
%     cost = [];
%     count = 1;
    feas_flag = 0;
    N = length(z)/9;
    numcon = length(h)/N;
    
    % Ensure a feasible start
    if any((h - P*z) <= 0)
%         disp('Started with infeasible solution.');
        infeas0 = find((h - P*z) <= 0);
%         disp(infeas0)
        

            % Find a new feasible z-val for the violated constraint(s)
            for loopind = 1:length(infeas0)
                ind = infeas0(loopind);
                
                if numcon == 8 % no ec constraints yet
                    ind_startrange = 7 + 8*idivide(ind,int16(8+1));
                    ind_endrange = ind_startrange + 1;
                else
                    ind_startrange = 7 + 10*idivide(ind,int16(10+1));
                    ind_endrange = ind_startrange + 3;
                end
%                 disp(ind_startrange)
%                 disp(ind_endrange)
                    
                z_infeas_ind = find(P(ind,:)); % optimization variable(s) to vary
    %             disp('z(infeas)=')
    %             disp(z_infeas_ind)

                newind = z_infeas_ind - 9;
                
                search_range = 1;
                search_stepnum = 1001;
                
                while (1)
                    if newind < 0
%                         disp("Sampling...")

                        % make a grid of search points
                        try_x = linspace(z(z_infeas_ind(1))-search_range,z(z_infeas_ind(1))+search_range,search_stepnum);
                        try_y = linspace(z(z_infeas_ind(2))-search_range,z(z_infeas_ind(2))+search_range,search_stepnum);
                        [xx,yy] = meshgrid(try_x,try_y);

                        % calcuate the constraint at each grid point
%                         temp = zeros(search_stepnum,search_stepnum);
                        temp2 = zeros(search_stepnum*search_stepnum,3);
                        for i1=1:search_stepnum
                            for i2 = 1:search_stepnum
                                temp = h(ind_startrange:ind_endrange) - P(ind_startrange:ind_endrange,z_infeas_ind)*[xx(i1,i2);yy(i1,i2)];

                                % take note of grid points that satisfy the constraint
                                if ~any(temp <= 0)
                                    temp2((i1-1)*search_stepnum+i2,1) = min(temp);
                                    temp2((i1-1)*search_stepnum+i2,2) = i1;
                                    temp2((i1-1)*search_stepnum+i2,3) = i2;
                                end
                            end
                        end

                        % get rid of data that does not satisfy the constraints
                        temp2( ~any(temp2,2), : ) = [];  %rows
                        
                        if isempty(temp2)
%                             if search_stepnum > search_stepnum * 2 * 2
%                                 break
%                             end
%                             search_range = 3;
%                             search_stepnum = search_stepnum * 2;
%                             continue
                            break
                        end

                        % choose one of the grid points that satisfy the constraint
                        % currently finds the first minimum value... could be further improved
                        ind_of_min = find(temp2(:,1) == min(temp2(:,1)));
                        ind_chosen = ind_of_min(1);

                        % update the x0 with the new grid point
                        [i1i2_chosen] = temp2(ind_chosen,[2,3]);
                        i1_chosen = i1i2_chosen(1);
                        i2_chosen = i1i2_chosen(2);
                        z(z_infeas_ind) = [xx(i1_chosen,i2_chosen);yy(i1_chosen,i2_chosen)];                     
                        
                    else
                        z(z_infeas_ind) = z(newind); % re-tracing previous sol'n
                    end

                    if (h(ind)- P(ind,:)*z) <= 0
                        if newind < 0
                            continue
                        else
                            newind = newind - 9;
                        end
                    else
%                         disp("Resolved.")
                        break
                        
                    end
                end
            end
%         % Find a new feasible z-val for the violated constraint(s)
%         for ind = infeas0.'
%             z_infeas_ind = find(P(ind,:)); % optimization variable(s) to vary
% %             disp('z(infeas)=')
% %             disp(z_infeas_ind)
%             
%             newind = z_infeas_ind - 9;
%             while (1)
%                 z(z_infeas_ind) = z(newind); % re-tracing previous sol'n
% 
%                 if (h(ind)- P(ind,:)*z) <= 0
%                     newind = newind - 9;
%                     
%                 else
%                     break
%                 end
%             end
%         end
        
    end
    
    for kappa = 0.01*10 % [10,1,0.1,0.01]
%         disp('--------------------------------')
%         disp(strcat('Kappa...',string(kappa)))
%         fi = P*z-h;
%         disp(fi(1:9).')

        for iiter = 0:maxiter

            % Compute the dual and primal residuals. Then compute the norm of 
            % the stacked variable r=(r_d,r_p), ||r||.
            d = 1./(h - P*z);
            
            r_d = 2*H*(z-zt) + kappa*P.'*d + C.'*nu ;
            r_p = C*z - b;
            r = sqrt(r_d.'*r_d + r_p.'*r_p );
            
            if (r < epsilon1) 
                break 
            end

            % Compute the matrices needed to calculate delta_z and delta_nu.
            % These matrices are based off of solving Equation (8) of
            % WangBoyd2010. Note that in WangBoyd2010, step 3) to determine
            % delta_z is incorrect since the derivation from Equation (8) calls
            % for using delta_nu instead of nu.
            d2 = diag(d)^2;
            Phi = 2*H + kappa*P.'*d2*P;
            
            [delta_z,delta_nu] = primaldualnewtonstep(Phi,r_p,r_d,C,N);
            
%             TS_dz(:,count) = delta_z;
%             count = count+1;
            
            % Feasibility search
            t = 1;
            feas_count = 0;
            while (1)
                feasz = z + t*delta_z;
%                 fi = (h - P*feasz);
                if any((h - P*feasz) <= 0)
                    t = BETA*t;
                else
                    break
                end
                
                feas_count = feas_count + 1;
                if feas_count > 100
                    feas_flag = 1;
                    infeas = find((h - P*feasz) <= 0);
%                     disp('No feasibility.')
%                     disp(infeas)
                    break
                end
            end
            
            if feas_flag == 1
                break
            end
            
            ztmp = z + t*delta_z;
            nutmp = nu + t*delta_nu;

            % Backtracking line search. From Boyd's Convex Optimization
            % textbook, Section 10.3.2, Algorithm 10.2.
            while (1)
                %while ||r(x + t∆xnt, ν + t∆νnt)|| > (1 − αt)||r(x, ν)||, t := βt.
                d = 1./(h - P*ztmp);
                
                r_d_tmp = 2*H*(ztmp-zt) + kappa*P.'*d + C.'*nutmp ;
                r_p_tmp = C*ztmp - b;
                r_tmp = sqrt(r_d_tmp.'*r_d_tmp + r_p_tmp.'*r_p_tmp );
                
                if (r_tmp <= (1-ALPHA*t)*r)
                    break
                end
                
                t = t*BETA;
                
                ztmp = z + t*delta_z;
                nutmp = nu + t*delta_nu;
                
            end

            z = ztmp;
            nu = nutmp;
        
        end
        
        if feas_flag == 1
            break
        end

    end

%     cost(1) = z.'*H*z - 2*zt.'*H*z + zt.' * H * zt; 
%     cost(2) = -kappa*sum(log(h - P*z));
%     cost(3) = kappa; 
%     temp1 = 2*H*(z-zt);
%     temp2 = kappa*P.'*d;
%     cost(4) = temp1(4);
%     cost(5) = temp2(4);

    zout = z;
    nuout = nu;
        
end

function [delta_z,delta_nu] = primaldualnewtonstep(Phi,rp,rd,C,N)
    Ad = -C(7:12,4:9);
    Bd = -C(1:6,1:3);
    n = 6;
    m = 3;
    
    % make PhiRinvBt, PhiQinv, PhiQinvAt
    invPhiQ = zeros(N*n,n);
    invPhiQAt = zeros(N*n,n);
    invPhiR = zeros(N*m,m);
    invPhiRBt = zeros(N*m,n);

    for ind = 1:N
    %         L = chol(PHI(1+(ind-1)*9:3+(ind-1)*9,1+(ind-1)*9:3+(ind-1)*9),'lower');
    %         temp1 = forwardsubstitution(L,Bd.');
    %         invPhiRBt(1+(ind-1)*3:ind*3,:) = backwardsubstitution(L.',temp1);
%         L = chol(Phi(1+(ind-1)*9:3+(ind-1)*9,1+(ind-1)*9:3+(ind-1)*9),'lower');
        L = my_cf(Phi(1+(ind-1)*9:3+(ind-1)*9,1+(ind-1)*9:3+(ind-1)*9));
        temp1 = forwardsubstitution(L,eye(m,m));
        temp2 = backwardsubstitution(L.',temp1);
        invPhiR(1+(ind-1)*3:ind*3,:) = temp2;
        invPhiRBt(1+(ind-1)*3:ind*3,:) = temp2*Bd.';

%         L = chol(Phi(4+(ind-1)*9:ind*9,4+(ind-1)*9:ind*9),'lower');
        L = my_cf(Phi(4+(ind-1)*9:ind*9,4+(ind-1)*9:ind*9));
        temp1 = forwardsubstitution(L,eye(n,n));
        temp2 = backwardsubstitution(L.',temp1);
        invPhiQ(1+(ind-1)*6:ind*6,:) = temp2;
        invPhiQAt(1+(ind-1)*6:ind*6,:) = temp2*Ad.';
    end
    
    % make Y
    Y = zeros(6*N,6*N);
    Y(1:6,1:6) = Bd*invPhiRBt(1:3,:) + invPhiQ(1:6,1:6);
    for ind=2:N
        Y(1+(ind-1)*6:ind*6,1+(ind-1)*6:ind*6) = Ad*invPhiQAt(1+(ind-2)*6:(ind-1)*6,:) +...
            Bd*invPhiRBt(1+(ind-2)*3:(ind-1)*3,:) + invPhiQ(1+(ind-1)*6:(ind)*6,:);
        Y(1+(ind-2)*6:(ind-1)*6,1+(ind-1)*6:ind*6) = -invPhiQAt(1+(ind-2)*6:(ind-1)*6,:);
        Y(1+(ind-1)*6:ind*6,1+(ind-2)*6:(ind-1)*6) = -Ad*invPhiQ(1+(ind-2)*6:(ind-1)*6,:);
    end
    
    % make CPhiinv
    CPhiinv = zeros(n*N,(n+m)*N);
    CPhiinv(1:n,1:m) = -Bd*invPhiR(1:m,:);
    CPhiinv(1:n,m+1:n+m) = invPhiQ(1:n,:);
    for ind=2:N
        rowst = (ind-1)*n+1;
        rowen = ind*n;
        CPhiinv(rowst:rowen,m+1+(ind-2)*(n+m):(ind-1)*(n+m)) = -Ad*invPhiQ(1+(ind-2)*n:(ind-1)*n,:);
        CPhiinv(rowst:rowen,1+(ind-1)*(n+m):m+(ind-1)*(n+m)) = -Bd*invPhiR(1+(ind-1)*m:(ind)*m,:);
        CPhiinv(rowst:rowen,1+m+(ind-1)*(n+m):(ind)*(n+m)) = invPhiQ(1+(ind-1)*n:(ind)*n,:);
    end
    
    % make beta
    beta = -rp + CPhiinv * rd;
    
    % cholesky factorization of Y
    LY = my_cf_bidiag(Y,n);
    
    % delta nu
    temp = forwardsubstitution(LY,-beta);
    delta_nu = backwardsubstitution(LY.',temp);
    
    % make Phiinv from invPhiR,invPhiQ
    Phiinv = zeros((n+m)*N,(n+m)*N);
    for ind = 0:N-1
        nm = n+m;
        Phiinv(1+ind*nm:m+ind*nm,1+ind*nm:m+ind*nm) = invPhiR(1+ind*m:(ind+1)*m,:);
        Phiinv(1+m+ind*nm:(ind+1)*nm,1+m+ind*nm:(ind+1)*nm) = invPhiQ(1+ind*n:(ind+1)*n,:);
    end
    
    % delta_z
    mrdmCtdnu = -rd - C.'*delta_nu;
    delta_z = Phiinv * mrdmCtdnu;
    
end

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
    
%     L(1:n,1:n) = chol(A(1:n,1:n),'lower');
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