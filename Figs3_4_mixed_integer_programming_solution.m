for seedd = 1:30
    seedd
    rng(seedd)

    %Number of APs
    L = 16;

    %Number of antennas per AP
    N = 4;

    %Number of UEs in the network
    K = 8;

    %Number of DUs
    W = 4;


    SEAll = 0.25:0.25:10;
    numberSetup = length(SEAll);

    Axx = zeros(K,L,numberSetup);
    Add = zeros(W,numberSetup);
    Azz = zeros(L,numberSetup);
    All = zeros(W,numberSetup);
    Arho = zeros(K*L,numberSetup);

    Bxx = zeros(K,L,numberSetup);
    Bdd = zeros(W,numberSetup);
    Bzz = zeros(L,numberSetup);
    Bll = zeros(W,numberSetup);
    Brho = zeros(K*L,numberSetup);


    TotalPower = zeros(2,numberSetup);


    %Number of channel realizations per setup
    nbrOfRealizations = 500;

    Nsmooth = 12;
    Nslot = 16;
    %Length of coherence block
    tau_c = Nsmooth*Nslot;


    %Length of pilot sequences
    tau_p = K;

    %Compute the prelog factor assuming only downlink data transmission
    preLogFactor = (tau_c-tau_p)/tau_c;

    %Angular standard deviation in the local scattering model (in radians)
    ASD_varphi = deg2rad(15);  %azimuth angle
    ASD_theta = deg2rad(15);   %elevation angle

    %Total uplink transmit power per UE (mW)
    p = 100;

    %Total downlink transmit power per AP (mW)
    rho_tot = 1000;

    %Generate one setup with UEs at random locations
    [gainOverNoisedB,R,pilotIndex,D,D_small,APpositions,UEpositions] = generateSetup(L,K,N,tau_p,1,0,ASD_varphi,ASD_theta);

    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cell-free APs
    [Hhat,H,B,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);

    % Full uplink power for the computation of precoding vectors using
    % virtual uplink-downlink duality
    p_full = p*ones(K,1);

    %Define the case when all APs serve all UEs
    D_all = ones(L,K);

    %Obtain the expectations 
    [signal_P_MMSE, signal2_P_MMSE, scaling_P_MMSE,...
        signal_P_RZF, signal2_P_RZF, scaling_P_RZF,...
        signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
        functionComputeExpectationsV2(Hhat,H,D_all,C,nbrOfRealizations,N,K,L,p_full);

    
    bk = zeros(L,K);
    Ck = zeros(L,L,K,K);

    %Go through all UEs
    for k = 1:K


        bk(:,k) = real(vec(signal_LP_MMSE(k,k,:)))./sqrt(scaling_LP_MMSE(:,k));

        %Go through all UEs
        for i = 1:K

            if i==k
                Ck(:,:,k,k) = bk(:,k)*bk(:,k)';
            else
                Ck(:,:,k,i) = diag(1./sqrt(scaling_LP_MMSE(:,i)))...
                    *(vec(signal_LP_MMSE(k,i,:))...
                    *vec(signal_LP_MMSE(k,i,:))')...
                    *diag(1./sqrt(scaling_LP_MMSE(:,i)));


            end

            for j = 1:L
                Ck(j,j,k,i) = signal2_LP_MMSE(k,i,j)/scaling_LP_MMSE(j,i);
            end
        end


    end



    %Take the real part (in the SINR expression,the imaginary terms cancel
    %each other)
    Ck = real(Ck);
    bb = zeros(L*K,K);

    CC = zeros(L*K,L*K,K);
    CC2 = zeros(L*K,L*K,K);
    for k = 1:K
        for i = 1:K
            CC((i-1)*L+1:i*L,(i-1)*L+1:i*L,k) = Ck(:,:,k,i);
        end
        CC2(:,:,k) = sqrtm(CC(:,:,k));
        bb((k-1)*L+1:k*L,k) = vec(bk(:,k));
    end

    %%
    for sss = 1:numberSetup
        sss
        gamma = 2^(SEAll(sss)/preLogFactor)-1;
        PAP0 =  6.8*N;
        DeltaTr = 4;
        PONU = 7.7;
        sigmaCool = 0.9;
        POLT = 20;
        Pdisp = 120;
        Pproc0 = 20.8;
        DeltaProc = 74;
        GOPSmax = 180;
        fs = 30.72*10^6;
        Ts = 71.4*10^(-6);
        NDFT = 2048;
        Nused = 1200;

        Nbits = 12;

        Cfilter = 40*N*fs/(10^9);
        CDFT = (8*N*NDFT*log2(NDFT)/Ts)/(10^9);

        SE0 = SEAll(sss);
        CprecodingAP = Nused/(Ts*tau_c*10^9)*...
            ((8*N*tau_p+8*N^2)*tau_p + (4*N^2+4*N)*tau_p + 8*(N^3-N)/3);
        CotherAP = ((Nbits/16)^(1.2))*(1.3*N) ...
            + ((Nbits/16)^(0.2))*(2.7*sqrt(N));
        ZLP = Cfilter + CDFT + CprecodingAP + CotherAP;

        CprecodingUE = Nused*(tau_c-tau_p)/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*(8*N^2);

        FLP = ((Nbits/16)^(1.2))*( 1.3*((SE0/6)^(1.5))*K )...
            + ((Nbits/16)^(1.2))*( 1.3*SE0/6*K )...
            + 8*SE0/6*K;

        XLP = CprecodingUE;

        Rmax = 10*10^9;
        Rfronthaul = 2*fs*Nbits*N;
        Wmax = floor(Rmax/Rfronthaul);

        PPl = PAP0 + PONU + DeltaProc/GOPSmax*ZLP/sigmaCool;


        %% Cell-free C-RAN

        cvx_begin quiet
        variable xx(K,L) binary
        variable dd(W,1) binary
        variable zz(L,1) binary
        variable ll(W,1) binary
        variable rho(K*L,1)
        minimize Pdisp + PPl*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10 ...
            + POLT*(1:W)*ll/sigmaCool + Pproc0*(1:W)*dd/sigmaCool ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool
        subject to
        sum(xx,2) >= ones(K,1);
        sum(zz) <= sum(sum(xx));
        for k = 1:K
            norm([CC2(:,:,k)*10*rho; 1]) <= sqrt((gamma+1)/gamma)*bb(:,k)'*10*rho;
        end
        sum(zz) <= Wmax*W;
        zz <= (sum(xx,1)).';
        K*zz >= (sum(xx,1)).';
        sum(zz)/Wmax <= (1:W)*ll;
        sum(zz)/Wmax >= (0:W-1)*ll;

        ZLP*sum(zz) + XLP*sum(sum(xx)) + FLP <= GOPSmax*(1:W)*dd;

        sum(dd) == 1;
        sum(ll) == 1;


        (1:W)*dd >= (1:W)*ll;
        zeros(K*L,1) <= rho;
        rho <= sqrt(rho_tot/100)*vec(xx.');

        for ell = 1:L
            norm(rho(ell:L:end,1)) <= sqrt(rho_tot/100)*zz(ell,1);
        end


        cvx_end


        TotalPower(1,sss) = Pdisp + PPl*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10 ...
            + POLT*(1:W)*ll/sigmaCool + Pproc0*(1:W)*dd/sigmaCool ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool;

        Axx(:,:,sss) = xx;
        Add(:,sss) = dd;
        Azz(:,sss) = zz;
        All(:,sss) = ll;
        Arho(:,sss) = rho;

        if isnan(TotalPower(1,sss))
            break
        end

        %% Small-cell C-RAN
        cvx_begin quiet
        variable xx(K,L) binary
        variable dd(W,1) binary
        variable zz(L,1) binary
        variable ll(W,1) binary
        variable rho(K*L,1)
        minimize Pdisp + PPl*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10 ...
            + POLT*(1:W)*ll/sigmaCool + Pproc0*(1:W)*dd/sigmaCool ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool
        subject to
        sum(xx,2) == ones(K,1);
        sum(zz) <= sum(sum(xx));
        for k = 1:K
            norm([CC2(:,:,k)*10*rho; 1]) <= sqrt((gamma+1)/gamma)*bb(:,k)'*10*rho;
        end
        sum(zz) <= Wmax*W;
        zz <= (sum(xx,1)).';
        K*zz >= (sum(xx,1)).';
        sum(zz)/Wmax <= (1:W)*ll;
        sum(zz)/Wmax >= (0:W-1)*ll;

        ZLP*sum(zz) + XLP*sum(sum(xx)) + FLP <= GOPSmax*(1:W)*dd;

        sum(dd) == 1;
        sum(ll) == 1;


        (1:W)*dd >= (1:W)*ll;
        zeros(K*L,1) <= rho;
        rho <= sqrt(rho_tot/100)*vec(xx.');

        for ell = 1:L
            norm(rho(ell:L:end,1)) <= sqrt(rho_tot/100)*zz(ell,1);
        end


        cvx_end


        TotalPower(2,sss) = Pdisp + PPl*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10 ...
            + POLT*(1:W)*ll/sigmaCool + Pproc0*(1:W)*dd/sigmaCool ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool;

        Bxx(:,:,sss) = xx;
        Bdd(:,sss) = dd;
        Bzz(:,sss) = zz;
        Bll(:,sss) = ll;
        Brho(:,sss) = rho;

    end
    save(strcat('MIPsim',num2str(seedd),'.mat'))
end