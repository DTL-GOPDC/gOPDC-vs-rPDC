function [rPDC,GPDC,OPDC,PDC,GOPDC,S] = PDC_dDTF_imag(A,C,p_opt,Fs,Fmax,Nf,Covariance)

%% Written by: Amir Omidvarnia
%%% Please see also the functions provided in BioSig Toolbox
%%% (available at: http://biosig.sourceforge.net)

% A = [A1 A2 ... Ap]
% Nf: Number of frequency points
% Fmax: maximum frequency limit (should be less than Fs/2)
syms l m w;
Matrix=@(l,m,w)[cos(l*w)*cos(m*w) cos(l*w)*sin(m*w);sin(l*w)*cos(m*w) sin(l*w)*sin(m*w)];
[CH,L,T] = size(A);
R=inv(Covariance);
PDC   = zeros(CH,CH,Nf,T); % Orthogonalized Partial Directed Coherence
OPDC   = zeros(CH,CH,Nf,T); % Orthogonalized Partial Directed Coherence
GPDC   = zeros(CH,CH,Nf,T); % Orthogonalized Partial Directed Coherence
GOPDC   = zeros(CH,CH,Nf,T); % Orthogonalized Partial Directed Coherence

S   = zeros(CH,CH,Nf,T); % Partial Directed Coherence
% f = (0:Nf-1)*((Fs/2)/Nf); % Frequency span ---> original code
f = (0:Nf-1)*(Fmax/Nf); % Frequency span
z = 1i*2*pi/Fs;

for t = 1 : T % Number of time points
    t
    clear A2 
    A2 = [eye(CH) -squeeze(A(:,:,t))];
    
    C_tmp = squeeze(C(:,:,t));
    Cd = diag(diag(C_tmp)); % Cd is useful for calculation of DC
    invCd = abs(inv(Cd));% invCd is useful for calculation of GPDC
    if(sum(sum(isinf(invCd)))>0)
        invCd = zeros(CH);
    end

    for n = 1 : Nf
        A_f = zeros(CH);
        for k = 1 : p_opt+1
            A_f = A_f + (A2(:,k*CH+(1-CH:0))*exp(z*(k-1)*f(n))); % sum(Ar(f)), r=0:p --> A(f) (if B==I)
        end
        
    
        H_f = inv(A_f);
        
        %% Power Spectral Density matrix
        S(:,:,n,t) = (H_f * C_tmp * H_f');
        
        for ii = 1:CH
            a_ii = squeeze(A_f(:,ii)); % ii'th column of A(f) --> ak
            a_denum(ii) = sqrt(a_ii'*a_ii); % hermitian(ak) * ak (k=1:CH)
            a2_denum(ii) = sqrt(a_ii'*invCd*a_ii); % for the GPDC - uses diagonal covariance information
        end

        PDC(:,:,n,t)  = A_f./a_denum(ones(1,CH),:); % Aij/(sqrt(aj'*aj))
        OPDC(:,:,n,t)  = (real(A_f).*(imag(A_f))./(a_denum(ones(1,CH),:)).^2); % Aij/(sqrt(aj'*aj))
        
        GPDC(:,:,n,t)  = sqrt(invCd)*A_f./a2_denum(ones(1,CH),:); % Aij/(sqrt(aj'*aj))
        GOPDC(:,:,n,t)  = (invCd)*(real(A_f).*(imag(A_f))./(a2_denum(ones(1,CH),:)).^2); % Aij/(sqrt(aj'*aj))
        

      
        for in_1=1:size(A_f,1)
            for in_2=1:size(A_f,2)
                V=zeros(2);
                for c1=1:p_opt
                    for c2=1:p_opt

                        V=V+(R(in_2,in_2)*C_tmp(in_1,in_1)*Matrix(c1,c2,2*pi*f(n)));
                    end
                end
                X=[real(A_f(in_1,in_2));imag(A_f(in_1,in_2))];
                
                inv_V=inv(V);
                if(sum(sum(isinf(inv_V)))>0)
                  inv_V = zeros(size(V,1));
                end
                lambda(in_1,in_2)=X'*inv_V*X;
            end 
        end
           rPDC(:,:,n,t)=lambda;  
    end

end


