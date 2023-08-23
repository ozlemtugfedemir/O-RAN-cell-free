scen = 1:30;
sinrlength = 40;
nbrOfRuns = length(scen);
totPow1 = zeros(nbrOfRuns,sinrlength);
totPow2 = zeros(nbrOfRuns,sinrlength);
totPow3 = zeros(nbrOfRuns,sinrlength);
totPow4 = zeros(nbrOfRuns,sinrlength);

totPow1b = zeros(nbrOfRuns,5);
totPow2b = zeros(nbrOfRuns,5);
totPow3b = zeros(nbrOfRuns,5);
totPow4b = zeros(nbrOfRuns,5);

totPow1c = zeros(nbrOfRuns,5);


for nnn = 1:nbrOfRuns
    seedd = scen(nnn);
    load(strcat('MIPsim',num2str(seedd),'.mat'));
    totPow1(nnn,:) = TotalPower(1,:);
    totPow2(nnn,:) = TotalPower(2,:);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for xxx = 1:sinrlength


        if isnan(totPow1(nnn,xxx)) || (totPow1(nnn,xxx)==0)
        else
            xx = Axx(:,:,xxx);
            zz = Azz(:,xxx);
            activeAPs = find(zz>0);

            rho = Arho(:,xxx);
            for randomtrial = 1:5
                orderOfLCs = randperm(16);
                APsLC1 = orderOfLCs(1:3);
                APsLC2 = orderOfLCs(4:6);
                APsLC3 = orderOfLCs(7:9);
                APsLC4 = orderOfLCs(10:12);
                APsLC5 = orderOfLCs(13:15);
                APsLC6 = orderOfLCs(16);
                APsDU1 = activeAPs(ismember(activeAPs,APsLC1)>0);
                APsDU2 = activeAPs(ismember(activeAPs,APsLC2)>0);
                APsDU3 = activeAPs(ismember(activeAPs,APsLC3)>0);
                APsDU4 = activeAPs(ismember(activeAPs,APsLC4)>0);
                APsDU5 = activeAPs(ismember(activeAPs,APsLC5)>0);
                APsDU6 = activeAPs(ismember(activeAPs,APsLC6)>0);
                NumberOfLCs = 6-isempty(APsDU1)-isempty(APsDU2)-isempty(APsDU3)...
                    -isempty(APsDU4)-isempty(APsDU5)-isempty(APsDU6);
                NumberOfDUs =  ceil((ZLP*sum(zz(APsDU1))+XLP*sum(sum(xx(:,APsDU1))))/GOPSmax) ...
                    +ceil((ZLP*sum(zz(APsDU2))+XLP*sum(sum(xx(:,APsDU2))))/GOPSmax) ...
                    +ceil((ZLP*sum(zz(APsDU3))+XLP*sum(sum(xx(:,APsDU3))))/GOPSmax) ...
                    +ceil((ZLP*sum(zz(APsDU4))+XLP*sum(sum(xx(:,APsDU4))))/GOPSmax) ...
                    +ceil((ZLP*sum(zz(APsDU5))+XLP*sum(sum(xx(:,APsDU5))))/GOPSmax) ...
                    +ceil((ZLP*sum(zz(APsDU6))+XLP*sum(sum(xx(:,APsDU6))))/GOPSmax);
               

                dd = Add(:,xxx);
                NumberOfDUs = max((1:W)*dd,NumberOfDUs);
                NumberOfDUs = max(NumberOfDUs,NumberOfLCs);

                SE0 = SEAll(xxx);
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
                NumberOfDUs2 = max(NumberOfDUs,ceil((ZLP*16+XLP*K*16+FLP)/GOPSmax));
                NumberOfDUs2 = max(NumberOfDUs2,6);

                tempor1 = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
                tempor2 = PONU*sum(zz) + POLT*NumberOfLCs/sigmaCool;
                tempor3 = Pdisp;
                tempor4 = Pproc0*NumberOfDUs/sigmaCool;
                tempor4b = Pproc0*NumberOfDUs2/sigmaCool;

                tempor5 = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
                    + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
                    + DeltaProc/GOPSmax*FLP/sigmaCool;
                totPow3(nnn,xxx) = totPow3(nnn,xxx) + (tempor1 + tempor2 + tempor3 + tempor4 + tempor5)/5;
                totPow4(nnn,xxx) = totPow4(nnn,xxx) + (tempor1 + tempor2 + tempor3 + tempor4b + tempor5)/5;

                if xxx == 5
                    totPow3b(nnn,1) = totPow3b(nnn,1) + tempor1/5;
                    totPow3b(nnn,2) = totPow3b(nnn,2) + tempor2/5;
                    totPow3b(nnn,3) = totPow3b(nnn,3) + tempor3/5;
                    totPow3b(nnn,4) = totPow3b(nnn,4) + tempor4/5;
                    totPow3b(nnn,5) = totPow3b(nnn,5) + tempor5/5;

                    totPow4b(nnn,1) = totPow4b(nnn,1) + tempor1/5;
                    totPow4b(nnn,2) = totPow4b(nnn,2) + tempor2/5;
                    totPow4b(nnn,3) = totPow4b(nnn,3) + tempor3/5;
                    totPow4b(nnn,4) = totPow4b(nnn,4) + tempor4b/5;
                    totPow4b(nnn,5) = totPow4b(nnn,5) + tempor5/5;
                end

            end
        end


    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xx = Axx(:,:,5);
    dd = Add(:,5);
    zz = Azz(:,5);
    ll = All(:,5);
    rho = Arho(:,5);

    SE0 = SEAll(5);
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

    totPow1b(nnn,1) = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
    totPow1b(nnn,2) = PONU*sum(zz) + POLT*(1:W)*ll/sigmaCool;
    totPow1b(nnn,3) = Pdisp;
    totPow1b(nnn,4) = Pproc0*(1:W)*dd/sigmaCool;
    totPow1b(nnn,5) = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
        + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
        + DeltaProc/GOPSmax*FLP/sigmaCool;

    %%%%%%%%%%%%
    xx = Axx(:,:,10);
    dd = Add(:,10);
    zz = Azz(:,10);
    ll = All(:,10);
    rho = Arho(:,10);

    SE0 = SEAll(10);
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


    totPow1c(nnn,1) = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
    totPow1c(nnn,2) = PONU*sum(zz) + POLT*(1:W)*ll/sigmaCool;
    totPow1c(nnn,3) = Pdisp;
    totPow1c(nnn,4) = Pproc0*(1:W)*dd/sigmaCool;
    totPow1c(nnn,5) = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
        + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
        + DeltaProc/GOPSmax*FLP/sigmaCool;

    %%%%%%%%%%%%
    xx = Bxx(:,:,5);
    dd = Bdd(:,5);
    zz = Bzz(:,5);
    ll = Bll(:,5);
    rho = Brho(:,5);

    SE0 = SEAll(5);
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


    totPow2b(nnn,1) = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
    totPow2b(nnn,2) = PONU*sum(zz) + POLT*(1:W)*ll/sigmaCool;
    totPow2b(nnn,3) = Pdisp;
    totPow2b(nnn,4) = Pproc0*(1:W)*dd/sigmaCool;
    totPow2b(nnn,5) = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
        + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
        + DeltaProc/GOPSmax*FLP/sigmaCool;
    %%%%%%%%%%%%%%%%%%%%

end
PowEnd1 = zeros(sinrlength,1);
PowEnd2 = zeros(sinrlength,1);
PowEnd3 = zeros(sinrlength,1);
PowEnd4 = zeros(sinrlength,1);

PowEnd1b = zeros(1,5);
PowEnd2b = zeros(1,5);
PowEnd3b = zeros(1,5);
PowEnd4b = zeros(1,5);



t1 = nbrOfRuns*ones(sinrlength,1);
t2 = nbrOfRuns*ones(sinrlength,1);
for x = 1:sinrlength

    for n = 1:nbrOfRuns
        if isnan(totPow1(n,x)) || (totPow1(n,x)==0)
            t1(x) = t1(x)-1;
        else
            PowEnd1(x) = PowEnd1(x) + totPow1(n,x);
            PowEnd3(x) = PowEnd3(x) + totPow3(n,x);
            PowEnd4(x) = PowEnd4(x) + totPow4(n,x);


            if x == 5
                PowEnd1b = PowEnd1b + totPow1b(n,:);
                PowEnd3b = PowEnd3b + totPow3b(n,:);
                PowEnd4b = PowEnd4b + totPow4b(n,:);

            end
        end

        if isnan(ttotPow2(n,x)) || (ttotPow2(n,x)==0)
            t2(x) = t2(x) -1;
        else
            PowEnd2(x) = PowEnd2(x) +totPow2(n,x);
           
            if x == 5
                PowEnd2b = PowEnd2b + totPow2b(n,:);

            end
        end
    end
end
%

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

set(groot,'defaultLegendInterpreter','latex');

figure


plot(SEAll(1:13),PowEnd1(1:13)./t1(1:13))
hold on
plot(SEAll(1:7),PowEnd2(1:7)./t2(1:7))
plot(SEAll(1:13),PowEnd3(1:13)./t1(1:13))
plot(SEAll(1:13),PowEnd4(1:13)./t1(1:13))
legend({'Cell-free (end-to-end)', 'Small-cell (end-to-end)', 'Cell-free (local coordination)', 'Cell-free (radio-only)'},'Interpreter','Latex','Location','Best');



Y = zeros(4,5);
Y(1,:) = PowEnd1b/t1(5);
Y(2,:) = PowEnd2b/t2(5);
Y(3,:) = PowEnd3b/t1(5);
Y(4,:) = PowEnd4b/t1(5);

X = categorical({'Cell-free (end-to-end)','Small-cell (end-to-end)','Cell-free (local coordination)','Cell-free (radio-only)'});
X = reordercats(X,{'Cell-free (end-to-end)','Small-cell (end-to-end)','Cell-free (local coordination)','Cell-free (radio-only)'});
figure
bar(X,Y);
legend({'RU', 'Fronthaul', 'Controller', 'GPP idle', 'GPP baseband'},'Interpreter','Latex','Location','Best');
ylabel('Power consumption (Watts)','Interpreter','Latex');
ylim([0 300])
set(gca,'fontsize',18);