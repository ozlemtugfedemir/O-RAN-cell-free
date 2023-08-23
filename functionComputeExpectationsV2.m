function [signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
    functionComputeExpectationsV2(Hhat,H,D,C,nbrOfRealizations,N,K,L,p)

%Store the N x N identity matrix
eyeN = eye(N);

%Obtain the diagonal matrix with UE transmit powers as its diagonal entries
PowMat = diag(p);

%Scale C by power coefficients
Cp = zeros(size(C));
for k=1:K
    Cp(:,:,:,k) = p(k)*C(:,:,:,k);
end

%Prepare to store simulation results

signal_LP_MMSE = zeros(K,K,L);
signal2_LP_MMSE = zeros(K,K,L);
scaling_LP_MMSE = zeros(L,K);

%% Compute scaling factors for combining/precoding

%Go through all channel realizations
for n=1:nbrOfRealizations
    
    
    %Go through all APs
    for l = 1:L
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(l-1)*N:l*N,n,:),[N K]);
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        %Obtain the statistical matrices used for
        %computing partial combining/precoding schemes
        Cpserved = reshape(sum(Cp(:,:,l,servedUEs),4),[N N]);
        Pserved = PowMat(servedUEs,servedUEs);
        
        %Compute MR combining scaled by square root of transmit powers
        Vp_MR = Hhatallj(:,servedUEs)*sqrt(Pserved);
        Vp_MRAll = Hhatallj*sqrt(PowMat);

        %Compute LP-MMSE combining
        V_LP_MMSE = (((Vp_MR*Vp_MR')+Cpserved+eyeN)\Vp_MRAll)*sqrt(PowMat);
        
        %Go through all UEs served by the AP
        for k = 1:K
            
            %Normalize LP-MMSE precoding
            w = V_LP_MMSE(:,k);
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in the SE expressions and
            %update Monte-Carlo estimates 
            signal2_LP_MMSE(:,k,l) = signal2_LP_MMSE(:,k,l) + abs(Hallj'*w).^2/nbrOfRealizations;
            
            signal_LP_MMSE(:,k,l) = signal_LP_MMSE(:,k,l) + Hallj'*w/nbrOfRealizations;
            
            scaling_LP_MMSE(l,k) = scaling_LP_MMSE(l,k) + sum(abs(w).^2,1)/nbrOfRealizations;

            
        end
    end
   
end