function  [X_trdoff] =  tradeoff_comrad_per_ant(rou,H,Y,power,X_arbi)
    [N,K] = size(H);
    [~,L] = size(Y);
    p1 = sqrt(rou);
    p2 = sqrt(1-rou);
    A = sqrt((L*power)/N)*[p1*H,p2*eye(N)].';
    A = A';
    B = [p1*sqrt(power)*Y;p2*X_arbi];
    B = B';


    % Create the problem structure.
    M = obliquecomplexfactory(L,N);
    problem.M = M;
     
    % Define the problem cost function and its gradient.
     
    problem.cost = @cost;
    function f=cost(X)
        f = norm(X*A-B,'fro')^2;
    end
    problem.grad = @(X) problem.M.egrad2rgrad(X,egrad(X));
    function g = egrad(X)
       g = 2*(X*A-B)*A';
    end
%     figure;
%     checkgradient(problem);   
    % Execute the optimization
    X_temp = conjugategradient(problem);
    X_trdoff = sqrt((L*power)/N)*X_temp';
end