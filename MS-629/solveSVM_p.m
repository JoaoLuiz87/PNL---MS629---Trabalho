function [x] = solveSVM_p(data, k_max)

[m, n] = size(data); 
n = n - 1;
options = optimoptions('fminunc','SpecifyObjectiveGradient',true, 'Display','off');

%PONTO INICIAL(origem)
%w = zeros(n, 1) -> x(1:n)
%csi = zeros(m, 1) -> x(n+1:n+m)
%t = zeros(m, 1) -> x(n+m+1:n+m+m)
%b = 0; -> x(n+m+m+1)
%x = (w,csi,t,b) VETOR c COMPOSTO POR w, csi, t, b. 

x = zeros(n + m + m + 1, 1);

rho_k = 2;
for k = 1:k_max

    x = fminunc(@(x)penalty(x, data ,m, n, rho_k,1000), x, options);
    rho_k = 2*rho_k;
end

end

