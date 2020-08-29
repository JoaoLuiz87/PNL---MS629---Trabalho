function [L, gradL] = Laum(x, lambda, data, m, n, rho, C)

%Vetor com a avaliação das restrições
h = data(:,1).*(data(:,2:end)*x(1:n) + x(end)) + x((n+1):(n+m)).^2 - x((n+m+1):(n+m+m)).^2 - 1;

%AVALIAÇÃO DO LAGRANGEANO AUMENTADO
L = dot(x(1:n),x(1:n))/2 + C*dot(x((n+1):(n+m)),x((n+1):(n+m))) + dot(lambda,h) + rho*dot(h,h);

if(nargout>1)
    %AVALIAÇÃO DO GRADIENTE DO LAGRANGEANO AUMENTADO
    grad_wh = data(:,2:end)'.*data(:,1)';
    grad_csih = diag(2*C*x(n+1:n+m));
    grad_th = diag(-2*x(n+m+1:n+m+m));
    grad_bh = data(:,1)';

    grad_h = [grad_wh; grad_csih; grad_th; grad_bh];

    gradL = zeros(n + m + m + 1, 1);

    gradL(1:n) = x(1:n);
    gradL(n+1:n+m) = 2*x(n+1:n+m);

    gradL = gradL + grad_h*lambda + rho*(grad_h*h);
end

end

