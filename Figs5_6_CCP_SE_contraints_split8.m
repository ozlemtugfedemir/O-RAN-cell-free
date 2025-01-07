for seedd = 1:30
    seedd
    rng(seedd)

    %Number of APs
    L = 36;

    %Number of antennas per AP
    N = 4;

    %Number of UEs in the network
    K = 16;

    %Number of DUs
    W = 40;

    SEAll = 1:4;
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
    tau_p = 8;

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

    


    [signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
        functionComputeExpectationsV2(Hhat,H,D,C,nbrOfRealizations,N,K,L,p_full);

    
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
            + Nused/(Ts*tau_c*10^9)*(8*N^2)*2;

        FLP = ((Nbits/16)^(1.2))*( 1.3*((SE0/6)^(1.5))*K )...
            + ((Nbits/16)^(1.2))*( 1.3*SE0/6*K )...
            + 8*SE0/6*K;

        XLP = CprecodingUE;

        Rmax = 10*10^9;
        Rfronthaul = 2*fs*Nbits*N;
        Wmax = floor(Rmax/Rfronthaul);

        PPz = PAP0 + PONU + DeltaProc/GOPSmax*ZLP/sigmaCool;


        PPz2 = PPz + Pproc0/GOPSmax*ZLP/sigmaCool ...
            +(POLT+Pproc0)/sigmaCool/Wmax;

        PPrho2 = (Pproc0 + DeltaProc)/GOPSmax*XLP/sigmaCool;
        %%  Cell-Free Low complex CCP

        Powerkapi = inf;

        for randomInitial = 1:5

            kapii2 = 1;
            rho0 = rand(L*K,1)*sqrt(rho_tot/100);
            zz0 = rand(L,1)*sqrt(rho_tot/100);


            kapii = 1;
            aExp = 3+(4-SE0);




            Gzz = aExp.*exp(-aExp*zz0);
            Grho = aExp.*exp(-aExp*rho0);

            Fzz = 1-exp(-aExp*zz0);
            Frho = 1-exp(-aExp*rho0);


            diff = 1;
            iterr = 0;
            ObjOld = PPz2*sum(Fzz) + DeltaTr*quad_form(rho0,eye(K*L))/10 ...
                + PPrho2*sum(Frho);

            while ((diff>0.00001)||(iterr<10))&&(iterr<50)
                [seedd sss diff]
                iterr = iterr+1;

                cvx_begin quiet
                variable rho(K*L,1)
                variable zz(L,1)

                minimize PPz2*(Gzz'*zz) + DeltaTr*quad_form(rho,eye(K*L))/10 ...
                    + PPrho2*(Grho'*rho)

                subject to

                for k = 1:K
                    norm([CC2(:,:,k)*10*rho; 1]) <= sqrt((gamma+1)/gamma)*bb(:,k)'*10*rho;
                end

                zeros(K*L,1) <= rho;

                for ell = 1:L
                    norm(rho(ell:L:end,1)) <= zz(ell,1);
                end

                for ell = 1:L
                    norm(rho(ell:L:end,1)) <= sqrt(rho_tot/100);
                end
                cvx_end

                if cvx_status(1)~='S'
                    kapii2 = 0;
                    break
                end

                rho0 = rho;
                zz0 = zz;

                Gzz = aExp.*exp(-aExp*zz0);
                Grho = aExp.*exp(-aExp*rho0);



                Fzz = 1-exp(-aExp*zz0);
                Frho = 1-exp(-aExp*rho0);


                ObjNew = PPz2*sum(Fzz) + DeltaTr*quad_form(rho0,eye(K*L))/10 ...
                    + PPrho2*sum(Frho);


                diff = abs(ObjNew-ObjOld)^2/abs(ObjOld)^2;

                ObjOld = ObjNew;




            end

            indexRho  = find(rho/max(rho)<=0.001);


            Powerkapi2 = inf;
            while kapii2>0
                rho(indexRho) = 0;
                zz = zeros(L,1);
                xx = sign(reshape(rho,[L,K])).';
                dd = zeros(W,1);
                ll = zeros(W,1);



                for ell = 1:L
                    if norm(rho(ell:L:end,1))>10^(-9)
                        zz(ell) = 1;
                    end


                end

                NumAPzero = L-sum(zz);

                llIndex = ceil(sum(zz)/Wmax);
                ll(llIndex) = 1;

                ddIndex = ceil((ZLP*sum(zz) + XLP*sum(sum(xx)) + FLP)/GOPSmax);
                ddIndex = max(ddIndex,llIndex);
                dd(ddIndex) = 1;

                cvx_begin quiet
                variable rho(K*L,1)
                minimize DeltaTr*quad_form(rho,eye(K*L))/10
                subject to

                for k = 1:K
                    norm([CC2(:,:,k)*10*rho; 1]) <= sqrt((gamma+1)/gamma)*bb(:,k)'*10*rho;
                end

                rho(indexRho) == 0;



                zeros(K*L,1) <= rho;

                for ell = 1:L
                    norm(rho(ell:L:end,1)) <= sqrt(rho_tot/100);
                end


                cvx_end
                if cvx_status(1)~='S'
                    break
                end
                if Powerkapi2 > Pdisp + PPz*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10 ...
                        + POLT*(1:length(ll))*ll/sigmaCool + Pproc0*(1:length(dd))*dd/sigmaCool ...
                        + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
                        + DeltaProc/GOPSmax*FLP/sigmaCool

                    Powerkapi2 = Pdisp + PPz*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10 ...
                        + POLT*(1:length(ll))*ll/sigmaCool + Pproc0*(1:length(dd))*dd/sigmaCool ...
                        + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
                        + DeltaProc/GOPSmax*FLP/sigmaCool;
                    Sxx = xx;
                    Sdd = dd;
                    Szz = zz;
                    Sll = ll;
                    Srho = rho;
                end

                [~,indexZZ] = sort(sum(reshape(rho,[L,K]).^2,2),'ascend');
                for zi = indexZZ(1:NumAPzero+1).'
                    for xi = zi:L:K*L
                        if any(indexRho == xi)
                        else
                            indexRho = [indexRho; xi];
                        end
                    end
                end

            end
            if Powerkapi > Powerkapi2

                TotalPower(1,sss) = Powerkapi2;

                Axx(:,:,sss) = Sxx;
                Add(:,sss) = Sdd;
                Azz(:,sss) = Szz;
                All(:,sss) = Sll;
                Arho(:,sss) = Srho;
                Powerkapi =   TotalPower(1,sss);
            end
        end



    end
    save(strcat('J_SE2_K',num2str(K),'_sim_8_',num2str(seedd),'.mat'))
end
