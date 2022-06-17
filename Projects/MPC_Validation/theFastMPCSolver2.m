function [zout,nuout] = theFastMPCSolver2(z,nu,zt,H,C,b,P,h,epsilon1,ALPHA,BETA,maxiter)
    
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
    
    if any((h - P*z) <= 0)
%         disp('Started with infeasible solution.');
        infeas0 = find((h - P*z) <= 0);
%         disp(infeas0)
        
        % Find a new feasible z-val for the violated constraint(s)
        for ind = infeas0.'
            z_infeas_ind = find(P(ind,:)); % optimization variable(s) to vary
%             disp('z(infeas)=')
%             disp(z_infeas_ind)
            
            newind = z_infeas_ind - 9;
            while (1)
                z(z_infeas_ind) = z(newind); % re-tracing previous sol'n

                if (h(ind)- P(ind,:)*z) <= 0
                    newind = newind - 9;
                    
                else
                    break
                end
            end
        end
        
    end
    
    for kappa = 0.01*5 % [10,1,0.1,0.01]
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
            Y = C*(Phi\C.');
            beta = -r_p + C*(Phi\r_d);
    %         disp(r)
    %         disp(sum(abs(r_d)))
    %         disp(sum(abs(r_p)))
    
            delta_nu = -Y\beta;
            delta_z = Phi\(-r_d - C.'*delta_nu); % NOTE: error in WangBoyd2010 on this
            
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
                    disp('No feasibility.')
                    disp(infeas)
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