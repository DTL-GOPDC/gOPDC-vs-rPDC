function [A,Sigma_e] = DEKF3(inp_model)
% A: Estimated time-varying parameters, A = [A1 A2 ... Ar]
% inp_model.data: (L x CH) data matrix
% inp_model.order: Model order
% inp_model.param{i}.vector: Data vector of the i'th time-varying parameter
% inp_model.param{i}.name: Name of the i'th time-varying parameter
% UC: Update Coefficient or Forgetting Factor Coefficient
%% Written by: Amir Omidvarnia
y = inp_model.data;
p = inp_model.order;
L = size(y,1); % Number of samples
%% 
y = y';
M = size(y,1);                    % Number of states (here, M = N)
LEN = size(y,2);                  % Number of the multivariate observations
%% Initial parameters for Dual Extended Kalman Filter
%%%%% (EKF 1)
xh = zeros(M*p,LEN);            % (EKF 1) Initial a-posteriori states (Mp x 1)
Px = .1*eye(M*p);               % (EKF 1) Initial a-posteriori state covariance matrix
% R = zeros(M,M,LEN);             % (EKF 1,2) Measurement error covariance matrix
% R(:,:,p+1) = 100*eye(M);          % (EKF 1,2) Initial observation noise covariance matrix
R = eye(M);
Q = 10*eye(M);                  % (EKF 1,2) Initial process noise covariance matrix
B = zeros(M*p,M);               % (EKF 1) Relationship between states and the process noise ( x(k) = F[x(k-1)] + B*v(k) )
B(1:M,:) = eye(M);              % (EKF 1) B = [I 0 ... 0]'
C = B';                         % (EKF 1,2) Measurement matrix (identity matrix, C = B')
% Kx = zeros(M*p,M,LEN);
% ah  = zeros(M*M*p,LEN);       % (EKF 2) Initial a-posteriori parameters estimaes
Ah = zeros(M*p,M*p,LEN);        % (EKF 2) Initial a-posteriori parameters estimates (matrix form of 'ah' plus identity matrices)
Ah(1:M,1:M*p,p) = .5*randn(M,M*p);     % (EKF 2) Initial a-posteriori parameters estimates (matrix form of 'ah' plus identity matrices)
for r = 2 : p
    Ah((r-1)*M+1:r*M,(r-2)*M+1:(r-1)*M, p) = eye(M);
end
%%%%% EKF 2
Pa = eye(M*M*p);                % (EKF 2) Initial a-posteriori parameters covariance matrix
% Ka = zeros(M*M*p,M,LEN);
for r = 1 : p
    xh((r-1)*M+1:r*M,p+1) = y(:,p-r+1);
end
Sigma_e = zeros(M,M,LEN);
%% DEKF starts ....
for i = p+1 : LEN
        
    [J_x J_A] = MVAR_JacCSD(Ah(:,:,i-1),xh(:,i-1),p); % xh(k) = F(A(k-1) * xh(k-1)) = Ah(k-1) * xh(k-1)
    Ah_ = Ah(:,:,i-1);                                 % Ah_(k) = Ah(k-1)
    %% EKF 1 ---> States estimation    
    %---------- Time Update (EKF1) ----------
    Rv = B * Q * B';                          % According to Haykin's book
    xh_ = Ah_ * xh(:,i-1);                    % xh_(k) = A_h(k-1) * xh(k-1)
    Px_ = J_x * Px * J_x' + Rv;               % Px_(k) = A_h(k-1) * Px(k-1) * A_h(k-1)' + B * Q * B'
    
    %---------- Measurement Update (EKF1) ----------
    Rn = R; %R(:,:,i-1);                                            % According to Haykin's book
    Kx = Px_ * C' * inv(C * Px_ * C' + Rn);            % Kx(k)  = Px_(k) * C' * inv(C * Px_(k) * C' + R)
    Px = (eye(M*p) - Kx * C) * Px_;                    % Px(k)  = (I - Kx(k) * C) * Px_(k)
    e = y(:,i) - C * xh_;          % inov(k) = y(k) - C * Ah_(k) * xh(k-1)
    xh(:,i) = xh_ + Kx * e;           % xh(k)  = xh_(k) + Kx(k) * (y(k) - C * xh_(k)) 
    
    
    %% EKF 2 ---> Parameters estimation
    %---------- Time Update (EKF2) ----------
    ah_ = reshape(Ah_(1:M,:)',M*M*p,1);                % ah_ = vec(Ah(k-1))
    Rr = .02*Pa;                                      % Rr = lambda * Pa(k-1)
    Pa_ = Pa + Rr;                                     % Pa_(k) = Pa(k-1) + Rr
    %---------- Measurement Update (EKF2) ----------
    %%%%%% Compute DfDa and H (H(k) = C*DfDa(k-1))
    H = C * (-J_A);                                    % J_A = -DfDa; ---> Why?
    %%%%%%%%%%%%%%%
    Re = (Rn + Q);                                        % According to Haykin's book
    Ka = Pa_ * H' * inv(H * Pa_ * H' + Re);            % Ka(k) = Pa_(k) * H(k) * inv(H(k) * Pa_(k) * H(k)' +R + Q)
    Pa = (eye(M*M*p) - Ka * H) * Pa_;                  % Pa(k) = (I - Ka(k) * H(k)) * Pa_(k)
    ah = ah_ + Ka * (y(:,i)- C * Ah_ * xh(:,i-1));                    % ah(k) = ah_(k) + Ka(k) * (y(k) - yh_(k))
      
    % Re-arrange vector ah(k) into the matrix Ah(k)
    Ah(1:M,1:M*p,i) = reshape(ah,M*p,M)';
    for r = 2 : p
        Ah((r-1)*M+1:r*M,(r-2)*M+1:(r-1)*M, i) = eye(M);
    end
%     Update the estimation error covariance matrix for dDTF, Partial
%     Coherence and adaptive AIC
        Sigma_e(:,:,i) = C * Px_ * C' + R; % Ref: Advanced digital Signal Processing and noise redution (4th ed.), Saeed Vaseghi, Eq. 7-18
%       i        
end
A = Ah(1:M,:,:); % Estimated time-varying parameters, A = [A1 A2 ... Ar]
