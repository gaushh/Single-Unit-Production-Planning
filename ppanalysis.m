function [R_utilized, B_utilized, product_used, prodprofit,product_produced, process_used, procprofit,amt_produced]  = ppanalysis(x, l, m, h, Cl, Cm, Ch, A4, A5,saleprice,product,NbyP) 
global bp_val


np = length(l);
x(x(1:np)>h) = h(x(1:np)>h);

nonzero_process = find(x(1:np));
process_used = {};
amt_produced = [];
for i = 1:length(nonzero_process)
    process_used{end+1} =strcat('P',num2str(nonzero_process(i)));
    amt_produced = [amt_produced;x(nonzero_process(i))];
end

PC = Cl.*x(np+1:2*np) + Cm.*x(2*np+1:3*np) + Ch.*x(3*np+1:4*np);

procprofit = saleprice.*x(1:np)-PC;

product_used = {};prodprofit = zeros(1,max(product));
product_produced = zeros(1,max(product));
for i=1:length(product)
    product_produced(product(i)) = x(i)+product_produced(product(i));  % each product from all process
    prodprofit(product(i)) = procprofit(i) + prodprofit(product(i));
end
b = find(product_produced);
for i = 1:length(b)
    product_used{end+1} =strcat('T',num2str(b(i)));
end

prodprofit = prodprofit(prodprofit>0);
procprofit = procprofit(procprofit>0);

if bp_val==1
    A5 = [A5 zeros(size(A5,1),NbyP)];
    A4 = [A4 zeros(size(A4,1),NbyP)];
end
    
R_utilized = A5*x;
B_utilized = sum(A4*x);