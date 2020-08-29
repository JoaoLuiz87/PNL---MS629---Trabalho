function [x] = solveSVM(data, k_max)

[m, n] = size(data);
n = n - 1;
options = optimoptions('fminunc','SpecifyObjectiveGradient',true, 'Display', 'off');

%PONTO INICIAL(origem)
%w = zeros(n, 1) -> x(1:n)
%csi = zeros(m, 1) -> x(n+1:n+m)
%t = zeros(m, 1) -> x(n+m+1:n+m+m)
%b = 0; -> x(n+m+m+1)
%x = (w,csi,t,b)

x = zeros(n + m + m + 1, 1);
lambda = zeros(m, 1);

rho_k = 2;
for k = 1:k_max

    x = fminunc(@(x)Laum(x, lambda, data ,m, n, rho_k,1000), x, options);
     
    h = data(:,1).*(data(:,2:end)*x(1:n) + x(end)) + x((n+1):(n+m)).^2 - x((n+m+1):(n+m+m)).^2 - 1;
    lambda = lambda + rho_k*h; 
    rho_k = 2*rho_k;
end

end

