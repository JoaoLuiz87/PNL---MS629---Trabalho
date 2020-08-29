clear 
load('dados.mat')

data = data(:,1:3);
x = solveSVM_p(data,300);

figure,
hold on
for i=1:569
    if(data(i,1)==-1)
        plot(data(i,2),data(i,3),'og')
    else
        plot(data(i,2),data(i,3),'xr')
    end
    
    aux = min(data(:,2)):0.01:max(data(:,3));
    plot(aux, -(x(1)*aux + x(end))/x(2))
end