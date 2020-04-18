function [data] = import_MNIST(file, digit)
%Import the csv data 'mnist_test.dat'
% make it into an list of n 28 by 28 matrices

M = readtable(file);
M = table2array(M);
M=M(M(:,1)==digit,:);
M=M(1:100,2:785);
temp = M;
data = reshape(temp(1,:),[28,28]);
for i = 2:100
    data = cat(28,data,reshape(temp(i,:),[28,28]));
end
end

