function [xk , history] = ADMMGroupLasso(A , b , BlockIndex , LenB , lambda, rho, opts)
%% ADMM to solve the group Lasso
[RA , CA]  = size(A);

xk       = rand(CA,1);
zk       = rand(CA,1);
uk       = rand(CA,1);

history  = struct();

% group information
numgroup         = CA-LenB*length(BlockIndex)+length(BlockIndex); 
BlockIndicator   = zeros(1,numgroup);
BlockIndicator(BlockIndex)  = 1;


Atb      = A'*b;
% pre-factor
[L U]    = factor(A, rho);

if ~opts.QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective','loss','regularization');
end

for iter = 1 : opts.itermax
    % x-update
    qk = Atb + rho*(zk - uk);    % temporary value
    if( RA >= CA )    % if skinny
       xk = U \ (L \ qk);
    else            % if fat
       xk = qk/rho - (A'*(U \ ( L \ (A*qk) )))/rho^2;
    end
    % z-update
    zold  = zk;
    xhat  = xk + uk;
%     for itergroup = 1 : CA/LenB
%         groupidx      = (itergroup-1)*LenB+1 : itergroup*LenB;
%         zk(groupidx)  = shrinkage(xhat(groupidx) , lambda/rho);
%     end
    currentindex = 0;
    for itergroup = 1 : numgroup
        if BlockIndicator(itergroup)
            groupidx  = currentindex+1 : currentindex+LenB;
        else
            groupidx  = currentindex+1;
        end
        currentindex  = groupidx(end);
        zk(groupidx)  = shrinkage(xhat(groupidx) , lambda/rho);
    end
    uk = uk + xk - zk;    
    % diagnostics, reporting, termination checks
    [history.objval(iter) , history.objloss(iter) , history.objreg(iter)]  = objective(A, b, lambda, LenB , xk);
    history.r_norm(iter)  = norm(xk - zk);
    history.s_norm(iter)  = norm(-rho*(zk - zold));

    history.eps_pri(iter) = sqrt(CA)*opts.abstol + opts.reltol*max(norm(xk), norm(-zk));
    history.eps_dual(iter)= sqrt(CA)*opts.abstol + opts.reltol*norm(rho*uk);


    if ~opts.QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\t%10.2f\t%10.2f\n', iter, ...
            history.r_norm(iter), history.eps_pri(iter), ...
            history.s_norm(iter), history.eps_dual(iter), history.objval(iter), history.objloss(iter) , history.objreg(iter));
    end

    if (history.r_norm(iter) < history.eps_pri(iter) && ...
       history.s_norm(iter) < history.eps_dual(iter))
         break;
    end
end

        
                       
function z = shrinkage(x, kappa)
    z = max(1 - kappa/norm(x),0)*x;
end

function [L U] = factor(A, rho)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
    end
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end

end
    


    
    
     
    










