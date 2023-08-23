scen = [ 1:30 ];

nbrOfRuns = length(scen);
totPow = zeros(6,3);


load('J_sumSEpen5_8.mat')
totPow(1,1) = mean(TotalPower(1,:));
load('J_sumSEpen5_72.mat')
totPow(2,1) = mean(TotalPower(1,:));

load('J_sumSE_pen50_8.mat')
totPow(3,1) = mean(TotalPower(1,:));
load('J_sumSE_pen50_72.mat')
totPow(4,1) = mean(TotalPower(1,:));

load('JOnlySumSE_sim_8.mat')
totPow(5,1) = mean(TotalPower(1,:));
load('JOnlySumSE_sim_72.mat')
totPow(6,1) = mean(TotalPower(1,:));


load('J_sumSEpen5_8.mat')

for nnn = 1:nbrOfRuns


    xx = Axx(:,:,nnn);
    zz = Azz(:,nnn);
    activeAPs = find(zz>0);
    SE = Aratee(:,nnn);
    rho = Arho(:,nnn);
    for randomtrial = 1:5
        orderOfLCs = randperm(36);
        APsLC1 = orderOfLCs(1:3);
        APsLC2 = orderOfLCs(4:6);
        APsLC3 = orderOfLCs(7:9);
        APsLC4 = orderOfLCs(10:12);
        APsLC5 = orderOfLCs(13:15);
        APsLC6 = orderOfLCs(16:18);
        APsLC7 = orderOfLCs(19:21);
        APsLC8 = orderOfLCs(22:24);
        APsLC9 = orderOfLCs(25:27);
        APsLC10 = orderOfLCs(28:30);
        APsLC11 = orderOfLCs(31:33);
        APsLC12 = orderOfLCs(34:36);

        APsDU1 = activeAPs(ismember(activeAPs,APsLC1)>0);
        APsDU2 = activeAPs(ismember(activeAPs,APsLC2)>0);
        APsDU3 = activeAPs(ismember(activeAPs,APsLC3)>0);
        APsDU4 = activeAPs(ismember(activeAPs,APsLC4)>0);
        APsDU5 = activeAPs(ismember(activeAPs,APsLC5)>0);
        APsDU6 = activeAPs(ismember(activeAPs,APsLC6)>0);
        APsDU7 = activeAPs(ismember(activeAPs,APsLC7)>0);
        APsDU8 = activeAPs(ismember(activeAPs,APsLC8)>0);
        APsDU9 = activeAPs(ismember(activeAPs,APsLC9)>0);
        APsDU10 = activeAPs(ismember(activeAPs,APsLC10)>0);
        APsDU11 = activeAPs(ismember(activeAPs,APsLC11)>0);
        APsDU12 = activeAPs(ismember(activeAPs,APsLC12)>0);
        NumberOfLCs = 12-isempty(APsDU1)-isempty(APsDU2)-isempty(APsDU3)...
            -isempty(APsDU4)-isempty(APsDU5)-isempty(APsDU6) ...
            -isempty(APsDU7)-isempty(APsDU8)-isempty(APsDU9) ...
            -isempty(APsDU10)-isempty(APsDU11)-isempty(APsDU12);
        NumberOfDUs =  ceil((ZLP*sum(zz(APsDU1))+XLP*sum(sum(xx(:,APsDU1))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU2))+XLP*sum(sum(xx(:,APsDU2))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU3))+XLP*sum(sum(xx(:,APsDU3))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU4))+XLP*sum(sum(xx(:,APsDU4))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU5))+XLP*sum(sum(xx(:,APsDU5))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU6))+XLP*sum(sum(xx(:,APsDU6))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU7))+XLP*sum(sum(xx(:,APsDU7))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU8))+XLP*sum(sum(xx(:,APsDU8))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU9))+XLP*sum(sum(xx(:,APsDU9))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU10))+XLP*sum(sum(xx(:,APsDU10))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU11))+XLP*sum(sum(xx(:,APsDU11))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU12))+XLP*sum(sum(xx(:,APsDU12))))/GOPSmax);



        dd = Add(:,nnn);

        NumberOfDUs = max((1:W)*dd, NumberOfDUs);
        NumberOfDUs = max(NumberOfDUs,NumberOfLCs);



        CprecodingAP = Nused/(Ts*tau_c*10^9)*...
            ((8*N*tau_p+8*N^2)*tau_p + (4*N^2+4*N)*tau_p + 8*(N^3-N)/3);
        CotherAP = ((Nbits/16)^(1.2))*(1.3*N) ...
            + ((Nbits/16)^(0.2))*(2.7*sqrt(N));
        ZLP = Cfilter + CDFT + CprecodingAP + CotherAP;

        CprecodingUE = Nused*(tau_c-tau_p)/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*(8*N^2)*2;

        FLP = sum(((Nbits/16)^(1.2))*( 1.3*((SE/6).^(1.5)) )...
            + ((Nbits/16)^(1.2))*( 1.3*SE/6 )...
            + 8*SE/6);

        XLP = CprecodingUE;

        NumberOfDUs2 = max(NumberOfDUs,ceil((ZLP*36+XLP*K*36+FLP)/GOPSmax));
        NumberOfDUs2 = max(NumberOfDUs2,12);


        tempor1 = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
        tempor2 = PONU*sum(zz) + POLT*NumberOfLCs/sigmaCool;
        tempor3 = Pdisp;
        tempor4 = Pproc0*NumberOfDUs/sigmaCool;
        tempor4b = Pproc0*NumberOfDUs2/sigmaCool;
        tempor5 = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool;
        totPow(1,2) = totPow(1,2) + (tempor1 + tempor2 + tempor3 + tempor4 + tempor5)/(5*30);
        totPow(1,3) = totPow(1,3) + (tempor1 + tempor2 + tempor3 + tempor4b + tempor5)/(5*30);




    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('J_sumSE_pen50_8.mat')

for nnn = 1:nbrOfRuns


    xx = Axx(:,:,nnn);
    zz = Azz(:,nnn);
    activeAPs = find(zz>0);

    rho = Arho(:,nnn);
    SE = Aratee(:,nnn);
    for randomtrial = 1:5
        orderOfLCs = randperm(36);
        APsLC1 = orderOfLCs(1:3);
        APsLC2 = orderOfLCs(4:6);
        APsLC3 = orderOfLCs(7:9);
        APsLC4 = orderOfLCs(10:12);
        APsLC5 = orderOfLCs(13:15);
        APsLC6 = orderOfLCs(16:18);
        APsLC7 = orderOfLCs(19:21);
        APsLC8 = orderOfLCs(22:24);
        APsLC9 = orderOfLCs(25:27);
        APsLC10 = orderOfLCs(28:30);
        APsLC11 = orderOfLCs(31:33);
        APsLC12 = orderOfLCs(34:36);

        APsDU1 = activeAPs(ismember(activeAPs,APsLC1)>0);
        APsDU2 = activeAPs(ismember(activeAPs,APsLC2)>0);
        APsDU3 = activeAPs(ismember(activeAPs,APsLC3)>0);
        APsDU4 = activeAPs(ismember(activeAPs,APsLC4)>0);
        APsDU5 = activeAPs(ismember(activeAPs,APsLC5)>0);
        APsDU6 = activeAPs(ismember(activeAPs,APsLC6)>0);
        APsDU7 = activeAPs(ismember(activeAPs,APsLC7)>0);
        APsDU8 = activeAPs(ismember(activeAPs,APsLC8)>0);
        APsDU9 = activeAPs(ismember(activeAPs,APsLC9)>0);
        APsDU10 = activeAPs(ismember(activeAPs,APsLC10)>0);
        APsDU11 = activeAPs(ismember(activeAPs,APsLC11)>0);
        APsDU12 = activeAPs(ismember(activeAPs,APsLC12)>0);
        NumberOfLCs = 12-isempty(APsDU1)-isempty(APsDU2)-isempty(APsDU3)...
            -isempty(APsDU4)-isempty(APsDU5)-isempty(APsDU6) ...
            -isempty(APsDU7)-isempty(APsDU8)-isempty(APsDU9) ...
            -isempty(APsDU10)-isempty(APsDU11)-isempty(APsDU12);
        NumberOfDUs =  ceil((ZLP*sum(zz(APsDU1))+XLP*sum(sum(xx(:,APsDU1))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU2))+XLP*sum(sum(xx(:,APsDU2))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU3))+XLP*sum(sum(xx(:,APsDU3))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU4))+XLP*sum(sum(xx(:,APsDU4))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU5))+XLP*sum(sum(xx(:,APsDU5))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU6))+XLP*sum(sum(xx(:,APsDU6))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU7))+XLP*sum(sum(xx(:,APsDU7))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU8))+XLP*sum(sum(xx(:,APsDU8))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU9))+XLP*sum(sum(xx(:,APsDU9))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU10))+XLP*sum(sum(xx(:,APsDU10))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU11))+XLP*sum(sum(xx(:,APsDU11))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU12))+XLP*sum(sum(xx(:,APsDU12))))/GOPSmax);



        dd = Add(:,nnn);

        NumberOfDUs = max((1:W)*dd, NumberOfDUs);
        NumberOfDUs = max(NumberOfDUs,NumberOfLCs);



        CprecodingAP = Nused/(Ts*tau_c*10^9)*...
            ((8*N*tau_p+8*N^2)*tau_p + (4*N^2+4*N)*tau_p + 8*(N^3-N)/3);
        CotherAP = ((Nbits/16)^(1.2))*(1.3*N) ...
            + ((Nbits/16)^(0.2))*(2.7*sqrt(N));
        ZLP = Cfilter + CDFT + CprecodingAP + CotherAP;

        CprecodingUE = Nused*(tau_c-tau_p)/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*(8*N^2)*2;

        FLP = sum(((Nbits/16)^(1.2))*( 1.3*((SE/6).^(1.5)) )...
            + ((Nbits/16)^(1.2))*( 1.3*SE/6 )...
            + 8*SE/6);

        XLP = CprecodingUE;

        NumberOfDUs2 = max(NumberOfDUs,ceil((ZLP*36+XLP*K*36+FLP)/GOPSmax));
        NumberOfDUs2 = max(NumberOfDUs2,12);

        tempor1 = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
        tempor2 = PONU*sum(zz) + POLT*NumberOfLCs/sigmaCool;
        tempor3 = Pdisp;
        tempor4 = Pproc0*NumberOfDUs/sigmaCool;
        tempor4b = Pproc0*NumberOfDUs2/sigmaCool;
        tempor5 = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool;
        totPow(3,2) = totPow(3,2) + (tempor1 + tempor2 + tempor3 + tempor4 + tempor5)/(5*30);
        totPow(3,3) = totPow(3,3) + (tempor1 + tempor2 + tempor3 + tempor4b + tempor5)/(5*30);





    end
end


%%%%%%%%%%%%%%%%%%%%%%%


load('JOnlySumSE_sim_8.mat')

for nnn = 1:nbrOfRuns


    xx = Axx(:,:,nnn);
    zz = Azz(:,nnn);
    activeAPs = find(zz>0);
    SE = Aratee(:,nnn);
    rho = Arho(:,nnn);
    for randomtrial = 1:5
        orderOfLCs = randperm(36);
        APsLC1 = orderOfLCs(1:3);
        APsLC2 = orderOfLCs(4:6);
        APsLC3 = orderOfLCs(7:9);
        APsLC4 = orderOfLCs(10:12);
        APsLC5 = orderOfLCs(13:15);
        APsLC6 = orderOfLCs(16:18);
        APsLC7 = orderOfLCs(19:21);
        APsLC8 = orderOfLCs(22:24);
        APsLC9 = orderOfLCs(25:27);
        APsLC10 = orderOfLCs(28:30);
        APsLC11 = orderOfLCs(31:33);
        APsLC12 = orderOfLCs(34:36);

        APsDU1 = activeAPs(ismember(activeAPs,APsLC1)>0);
        APsDU2 = activeAPs(ismember(activeAPs,APsLC2)>0);
        APsDU3 = activeAPs(ismember(activeAPs,APsLC3)>0);
        APsDU4 = activeAPs(ismember(activeAPs,APsLC4)>0);
        APsDU5 = activeAPs(ismember(activeAPs,APsLC5)>0);
        APsDU6 = activeAPs(ismember(activeAPs,APsLC6)>0);
        APsDU7 = activeAPs(ismember(activeAPs,APsLC7)>0);
        APsDU8 = activeAPs(ismember(activeAPs,APsLC8)>0);
        APsDU9 = activeAPs(ismember(activeAPs,APsLC9)>0);
        APsDU10 = activeAPs(ismember(activeAPs,APsLC10)>0);
        APsDU11 = activeAPs(ismember(activeAPs,APsLC11)>0);
        APsDU12 = activeAPs(ismember(activeAPs,APsLC12)>0);
        NumberOfLCs = 12-isempty(APsDU1)-isempty(APsDU2)-isempty(APsDU3)...
            -isempty(APsDU4)-isempty(APsDU5)-isempty(APsDU6) ...
            -isempty(APsDU7)-isempty(APsDU8)-isempty(APsDU9) ...
            -isempty(APsDU10)-isempty(APsDU11)-isempty(APsDU12);
        NumberOfDUs =  ceil((ZLP*sum(zz(APsDU1))+XLP*sum(sum(xx(:,APsDU1))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU2))+XLP*sum(sum(xx(:,APsDU2))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU3))+XLP*sum(sum(xx(:,APsDU3))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU4))+XLP*sum(sum(xx(:,APsDU4))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU5))+XLP*sum(sum(xx(:,APsDU5))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU6))+XLP*sum(sum(xx(:,APsDU6))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU7))+XLP*sum(sum(xx(:,APsDU7))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU8))+XLP*sum(sum(xx(:,APsDU8))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU9))+XLP*sum(sum(xx(:,APsDU9))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU10))+XLP*sum(sum(xx(:,APsDU10))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU11))+XLP*sum(sum(xx(:,APsDU11))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU12))+XLP*sum(sum(xx(:,APsDU12))))/GOPSmax);



        dd = Add(:,nnn);

        NumberOfDUs = max((1:W)*dd, NumberOfDUs);
        NumberOfDUs = max(NumberOfDUs,NumberOfLCs);


        CprecodingAP = Nused/(Ts*tau_c*10^9)*...
            ((8*N*tau_p+8*N^2)*tau_p + (4*N^2+4*N)*tau_p + 8*(N^3-N)/3);
        CotherAP = ((Nbits/16)^(1.2))*(1.3*N) ...
            + ((Nbits/16)^(0.2))*(2.7*sqrt(N));
        ZLP = Cfilter + CDFT + CprecodingAP + CotherAP;

        CprecodingUE = Nused*(tau_c-tau_p)/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*(8*N^2)*2;

        FLP = sum(((Nbits/16)^(1.2))*( 1.3*((SE/6).^(1.5)) )...
            + ((Nbits/16)^(1.2))*( 1.3*SE/6 )...
            + 8*SE/6);

        XLP = CprecodingUE;

        NumberOfDUs2 = max(NumberOfDUs,ceil((ZLP*36+XLP*K*36+FLP)/GOPSmax));
        NumberOfDUs2 = max(NumberOfDUs2,12);

        tempor1 = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
        tempor2 = PONU*sum(zz) + POLT*NumberOfLCs/sigmaCool;
        tempor3 = Pdisp;
        tempor4 = Pproc0*NumberOfDUs/sigmaCool;
        tempor4b = Pproc0*NumberOfDUs2/sigmaCool;
        tempor5 = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool;
        totPow(5,2) = totPow(5,2) + (tempor1 + tempor2 + tempor3 + tempor4 + tempor5)/(5*30);
        totPow(5,3) = totPow(5,3) + (tempor1 + tempor2 + tempor3 + tempor4b + tempor5)/(5*30);





    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


load('J_sumSEpen5_72.mat')


for nnn = 1:nbrOfRuns

    xx = Axx(:,:,nnn);
    zz = Azz(:,nnn);
    activeAPs = find(zz>0);
    SE = Aratee(:,nnn);
    rho = Arho(:,nnn);
    for randomtrial = 1:5
        orderOfLCs = randperm(36);
        APsLC1 = orderOfLCs(1:6);
        APsLC2 = orderOfLCs(7:12);
        APsLC3 = orderOfLCs(13:18);
        APsLC4 = orderOfLCs(19:24);
        APsLC5 = orderOfLCs(25:30);
        APsLC6 = orderOfLCs(31:36);

        APsDU1 = activeAPs(ismember(activeAPs,APsLC1)>0);
        APsDU2 = activeAPs(ismember(activeAPs,APsLC2)>0);
        APsDU3 = activeAPs(ismember(activeAPs,APsLC3)>0);
        APsDU4 = activeAPs(ismember(activeAPs,APsLC4)>0);
        APsDU5 = activeAPs(ismember(activeAPs,APsLC5)>0);
        APsDU6 = activeAPs(ismember(activeAPs,APsLC6)>0);

        NumberOfLCs = 6-isempty(APsDU1)-isempty(APsDU2) ...
            -isempty(APsDU3)-isempty(APsDU4) ...
            -isempty(APsDU5)-isempty(APsDU6);
        NumberOfDUs =  ceil((ZLP*sum(zz(APsDU1))+XLP*sum(sum(xx(:,APsDU1))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU2))+XLP*sum(sum(xx(:,APsDU2))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU3))+XLP*sum(sum(xx(:,APsDU3))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU4))+XLP*sum(sum(xx(:,APsDU4))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU5))+XLP*sum(sum(xx(:,APsDU5))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU6))+XLP*sum(sum(xx(:,APsDU6))))/GOPSmax);


        dd = Add(:,nnn);

        NumberOfDUs = max((1:W)*dd,NumberOfDUs);
        NumberOfDUs = max(NumberOfDUs,NumberOfLCs);


        CprecodingAP = Nused/(Ts*tau_c*10^9)*...
            ((8*N*tau_p+8*N^2)*tau_p + (4*N^2+4*N)*tau_p + 8*(N^3-N)/3);
        CotherAP = ((Nbits/16)^(1.2))*(1.3*N) ...
            + ((Nbits/16)^(0.2))*(2.7*sqrt(N));
        SLP = Cfilter + CDFT;
        ZLP = CprecodingAP + CotherAP;

        CprecodingUE = Nused*(tau_c-tau_p)/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*(8*N^2)*2;

        FLP = sum(((Nbits/16)^(1.2))*( 1.3*((SE/6).^(1.5)) )...
            + ((Nbits/16)^(1.2))*( 1.3*SE/6 )...
            + 8*SE/6);

        XLP = CprecodingUE;
        NumberOfDUs2 = max(NumberOfDUs,ceil((ZLP*36+XLP*K*36+FLP)/GOPSmax));
        NumberOfDUs2 = max(NumberOfDUs2,6);

        tempor1 = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
        tempor2 = PONU*sum(zz) + POLT*NumberOfLCs/sigmaCool;
        tempor3 = Pdisp;
        tempor4 = Pproc0*NumberOfDUs/sigmaCool;
        tempor4b = Pproc0*NumberOfDUs2/sigmaCool;
        tempor4c = PRUproc0*sum(zz);

        tempor5 = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool;
        tempor5b = DeltaRUProc/GOPSmax*SLP*sum(zz);
        totPow(2,2) = totPow(2,2)+(tempor1 + tempor2 + tempor3 + tempor4 + tempor5 ...
            + tempor4c + tempor5b)/(5*30);
        totPow(2,3) = totPow(2,3)+(tempor1 + tempor2 + tempor3 + tempor4b + tempor5 ...
            + tempor4c + tempor5b)/(5*30);




    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('J_sumSE_pen50_72.mat')


for nnn = 1:nbrOfRuns

    xx = Axx(:,:,nnn);
    zz = Azz(:,nnn);
    activeAPs = find(zz>0);
    SE = Aratee(:,nnn);
    rho = Arho(:,nnn);
    for randomtrial = 1:5
        orderOfLCs = randperm(36);
        APsLC1 = orderOfLCs(1:6);
        APsLC2 = orderOfLCs(7:12);
        APsLC3 = orderOfLCs(13:18);
        APsLC4 = orderOfLCs(19:24);
        APsLC5 = orderOfLCs(25:30);
        APsLC6 = orderOfLCs(31:36);

        APsDU1 = activeAPs(ismember(activeAPs,APsLC1)>0);
        APsDU2 = activeAPs(ismember(activeAPs,APsLC2)>0);
        APsDU3 = activeAPs(ismember(activeAPs,APsLC3)>0);
        APsDU4 = activeAPs(ismember(activeAPs,APsLC4)>0);
        APsDU5 = activeAPs(ismember(activeAPs,APsLC5)>0);
        APsDU6 = activeAPs(ismember(activeAPs,APsLC6)>0);

        NumberOfLCs = 6-isempty(APsDU1)-isempty(APsDU2) ...
            -isempty(APsDU3)-isempty(APsDU4) ...
            -isempty(APsDU5)-isempty(APsDU6);
        NumberOfDUs =  ceil((ZLP*sum(zz(APsDU1))+XLP*sum(sum(xx(:,APsDU1))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU2))+XLP*sum(sum(xx(:,APsDU2))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU3))+XLP*sum(sum(xx(:,APsDU3))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU4))+XLP*sum(sum(xx(:,APsDU4))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU5))+XLP*sum(sum(xx(:,APsDU5))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU6))+XLP*sum(sum(xx(:,APsDU6))))/GOPSmax);


        dd = Add(:,nnn);

        NumberOfDUs = max((1:W)*dd,NumberOfDUs);
        NumberOfDUs = max(NumberOfDUs,NumberOfLCs);


        CprecodingAP = Nused/(Ts*tau_c*10^9)*...
            ((8*N*tau_p+8*N^2)*tau_p + (4*N^2+4*N)*tau_p + 8*(N^3-N)/3);
        CotherAP = ((Nbits/16)^(1.2))*(1.3*N) ...
            + ((Nbits/16)^(0.2))*(2.7*sqrt(N));
        SLP = Cfilter + CDFT;
        ZLP = CprecodingAP + CotherAP;

        CprecodingUE = Nused*(tau_c-tau_p)/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*(8*N^2)*2;

        FLP = sum(((Nbits/16)^(1.2))*( 1.3*((SE/6).^(1.5)) )...
            + ((Nbits/16)^(1.2))*( 1.3*SE/6 )...
            + 8*SE/6);

        XLP = CprecodingUE;
        NumberOfDUs2 = max(NumberOfDUs,ceil((ZLP*36+XLP*K*36+FLP)/GOPSmax));
        NumberOfDUs2 = max(NumberOfDUs2,6);

        tempor1 = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
        tempor2 = PONU*sum(zz) + POLT*NumberOfLCs/sigmaCool;
        tempor3 = Pdisp;
        tempor4 = Pproc0*NumberOfDUs/sigmaCool;
        tempor4b = Pproc0*NumberOfDUs2/sigmaCool;
        tempor4c = PRUproc0*sum(zz);

        tempor5 = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool;
        tempor5b = DeltaRUProc/GOPSmax*SLP*sum(zz);
        totPow(4,2) = totPow(4,2)+(tempor1 + tempor2 + tempor3 + tempor4 + tempor5 ...
            + tempor4c + tempor5b)/(5*30);
        totPow(4,3) = totPow(4,3)+(tempor1 + tempor2 + tempor3 + tempor4b + tempor5 ...
            + tempor4c + tempor5b)/(5*30);




    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('JOnlySumSE_sim_72.mat')


for nnn = 1:nbrOfRuns

    xx = Axx(:,:,nnn);
    zz = Azz(:,nnn);
    activeAPs = find(zz>0);
    SE = Aratee(:,nnn);
    rho = Arho(:,nnn);
    for randomtrial = 1:5
        orderOfLCs = randperm(36);
        APsLC1 = orderOfLCs(1:6);
        APsLC2 = orderOfLCs(7:12);
        APsLC3 = orderOfLCs(13:18);
        APsLC4 = orderOfLCs(19:24);
        APsLC5 = orderOfLCs(25:30);
        APsLC6 = orderOfLCs(31:36);

        APsDU1 = activeAPs(ismember(activeAPs,APsLC1)>0);
        APsDU2 = activeAPs(ismember(activeAPs,APsLC2)>0);
        APsDU3 = activeAPs(ismember(activeAPs,APsLC3)>0);
        APsDU4 = activeAPs(ismember(activeAPs,APsLC4)>0);
        APsDU5 = activeAPs(ismember(activeAPs,APsLC5)>0);
        APsDU6 = activeAPs(ismember(activeAPs,APsLC6)>0);

        NumberOfLCs = 6-isempty(APsDU1)-isempty(APsDU2) ...
            -isempty(APsDU3)-isempty(APsDU4) ...
            -isempty(APsDU5)-isempty(APsDU6);
        NumberOfDUs =  ceil((ZLP*sum(zz(APsDU1))+XLP*sum(sum(xx(:,APsDU1))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU2))+XLP*sum(sum(xx(:,APsDU2))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU3))+XLP*sum(sum(xx(:,APsDU3))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU4))+XLP*sum(sum(xx(:,APsDU4))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU5))+XLP*sum(sum(xx(:,APsDU5))))/GOPSmax) ...
            +ceil((ZLP*sum(zz(APsDU6))+XLP*sum(sum(xx(:,APsDU6))))/GOPSmax);



        dd = Add(:,nnn);

        NumberOfDUs = max((1:W)*dd,NumberOfDUs);
        NumberOfDUs = max(NumberOfDUs,NumberOfLCs);


        CprecodingAP = Nused/(Ts*tau_c*10^9)*...
            ((8*N*tau_p+8*N^2)*tau_p + (4*N^2+4*N)*tau_p + 8*(N^3-N)/3);
        CotherAP = ((Nbits/16)^(1.2))*(1.3*N) ...
            + ((Nbits/16)^(0.2))*(2.7*sqrt(N));
        SLP = Cfilter + CDFT;
        ZLP = CprecodingAP + CotherAP;

        CprecodingUE = Nused*(tau_c-tau_p)/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*8*N ...
            + Nused/(Ts*tau_c*10^9)*(8*N^2)*2;

        FLP = sum(((Nbits/16)^(1.2))*( 1.3*((SE/6).^(1.5)) )...
            + ((Nbits/16)^(1.2))*( 1.3*SE/6 )...
            + 8*SE/6);

        XLP = CprecodingUE;
        NumberOfDUs2 = max(NumberOfDUs,ceil((ZLP*36+XLP*K*36+FLP)/GOPSmax));
        NumberOfDUs2 = max(NumberOfDUs2,6);

        tempor1 = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
        tempor2 = PONU*sum(zz) + POLT*NumberOfLCs/sigmaCool;
        tempor3 = Pdisp;
        tempor4 = Pproc0*NumberOfDUs/sigmaCool;
        tempor4b = Pproc0*NumberOfDUs2/sigmaCool;
        tempor4c = PRUproc0*sum(zz);

        tempor5 = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
            + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
            + DeltaProc/GOPSmax*FLP/sigmaCool;
        tempor5b = DeltaRUProc/GOPSmax*SLP*sum(zz);
        totPow(6,2) = totPow(6,2)+(tempor1 + tempor2 + tempor3 + tempor4 + tempor5 ...
            + tempor4c + tempor5b)/(5*30);
        totPow(6,3) = totPow(6,3)+(tempor1 + tempor2 + tempor3 + tempor4b + tempor5 ...
            + tempor4c + tempor5b)/(5*30);




    end

end


set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

set(groot,'defaultLegendInterpreter','latex');





X = categorical({'Low-weight (FS-8)','Low-weight (FS-7.2)','Medium-weight (FS-8)','Medium-weight (FS-7.2)','Only-sumSE-max (FS-8)','Only-sumSE-max (FS-7.2)'});
X = reordercats(X,{'Low-weight (FS-8)','Low-weight (FS-7.2)','Medium-weight (FS-8)','Medium-weight (FS-7.2)','Only-sumSE-max (FS-8)','Only-sumSE-max (FS-7.2)'});
figure
bar(X,totPow);
legend({'End-to-end', 'Local coordination', 'Radio-only'},'Interpreter','Latex','Location','Best');
ylabel('Power consumption (Watts)','Interpreter','Latex');
ylim([0 4000])
set(gca,'fontsize',18);