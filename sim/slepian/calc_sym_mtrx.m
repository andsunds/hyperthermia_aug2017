function [ A ] = calc_sym_mtrx( fun_of_2_args, vector )
%Given a symmetric function of two variables in the same value inteval and
%a vector containing the values of the inteval, this function calculates
%the matrix 
%   A(i,j) = fun_of_2_args(vector(i), vector(j));
%It is implicitly understood that 
%   fun_of_2_args(x,y) = fun_of_2_args(x,y)



L=length(vector); %We get the discretization size from x.

tri=@(i) i.*(i+1)/2; %the i'th triangle number
%An indexed lower triangle of a matrix together with row and col numbers.
%   m:  _  
%   1  |1|_  :j
%   2  |2 3|_ 
%   3  |4 5 6|__
%   4  |7_8_9_10|
%    n: 1 2 3 4
%Given an index j, the row (m) and col (n) can be calculated as below. To
%calculate m, all we have to do is invert the triangle number formula
J = 1:tri(L); %all indeces
m = ceil(-.5 + sqrt(.25 +2*J));
n = J -tri(m-1);

%For each j, we need the values x(m(j)) and y(n(j)). To optimize the
%parallellization, new (and longer) vectors X and Y are created precisely
%so that X(j) = x(m(j)) and Y(j) = y(n(j)).
v1=vector(m); 
v2=vector(n); %x(n) is sufficient since x and y are discretized identically

a=zeros(1,tri(L)); %init
for j=J
    a(j)=fun_of_2_args(v1(j), v2(j));
end

%Transfering the parallellized (1D) results to the matrix A
A=zeros(L);
for j=J
    A(m(j),n(j)) = a(j);
end
%A is symmetric, which is why we only calculate half of A and then flip it
%over on itself; if just A+A', we get double the value on the diagonal
%which is why -diag(diag(A)) is needed. 
A=A+A.'-diag(diag(A));
% note: important to use .'. If rot_order =/= 0, A is not hermetian, only
% symmetric.

end

