function [J_x J_A] = MVAR_JacCSD(A,x_ext,p)
%% Complex step differentiation
%%% Written by: Amir Omidvarnia
%%% See also: http://www.mathworks.com.au/matlabcentral/fileexchange/18176-complex-step-jacobian/content/jacobiancsd.m

% x_ext: (Mp x 1) Input vector of the MVAR model, x_ext(k-1) = [x(k-1); x(k-2); ...; x(k-p)]
% A: (Mp x Mp) All MVAR matrices in a vector column, A = [A1 A2 ... Ap; I 0 0 .. 0; 0 I 0 ... 0; ...; 0 0 ... I 0]
% F: (Mp x 1) Output of the model
% x_ext(k) = F(x_ext(k-1)) = A * x_ext(k-1)
%
% J_x: (Mp x Mp) Jacobian matrix through Complex Step Differentiation (CSD) --> df/dx
% J_A: (Mp x M*Mp) Jacobian matrix through Complex Step Differentiation (CSD) --> df/dA
%
% f(x0 + ih) = f(x0) + ih*f'(x0) + (ih)^2*f''(x0)/2! + ...
% --> f'(x0) = imag( f(x0+ih)/h ) (Approximately)

M = length(x_ext)/p; % Number of states (Mp)
[Mp Mp] = size(A); % MVAR coefficients matrix (Mp x Mp)
h = .00000001;
J_x = zeros(Mp,Mp); % df/dx
J_A = zeros(Mp,M*Mp);

%% df/dx
for i = 1 : Mp
    tmp_x = x_ext;
    tmp_x(i) = tmp_x(i) + 1i*h; % x0 + ih
    tmp_fx = A * tmp_x + randn(Mp,1); % f(x0+ih)
    J_x(:,i) = imag(tmp_fx/h); % f'(x0) = imag( f(x0+ih)/h ) (Approximately) ----> df/dx
end

%% df/dA ----> DfDa = -J_A (why?)
a = reshape(A(1:M,:)',M*Mp,1);
for i = 1 : M*Mp
    tmp_a = a;
    tmp_a(i) = tmp_a(i) + 1i*h; % a0 + ih
    tmp_A = reshape(tmp_a, M*p, M)';
    for r = 2 : p
        tmp_A((r-1)*M+1:r*M,(r-2)*M+1:(r-1)*M) = eye(M);
    end
    tmp_fa = tmp_A * x_ext + randn(Mp,1); % f(a0+ih)
    J_A(:,i) = imag(tmp_fa/h);
end