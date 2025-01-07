scen = [ 1:30 ];
Userlength = 4;
nbrOfRuns = length(scen);
totPow1 = zeros(nbrOfRuns,Userlength,2);
totPow3 = zeros(nbrOfRuns,Userlength,2);
totPow4 = zeros(nbrOfRuns,Userlength,2);

totPow1b = zeros(nbrOfRuns,6,2);
totPow3b = zeros(nbrOfRuns,6,2);
totPow4b = zeros(nbrOfRuns,6,2);



for nnn = 1:nbrOfRuns
    seedd = scen(nnn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for xxx = 1:Userlength
        if xxx == 1
            load(strcat('J_SE2_K4_sim_8_',num2str(seedd),'.mat'));
            totPow1(nnn,xxx,1) = TotalPower(1,:);

        elseif xxx == 2
            load(strcat('J_SE2_K8_sim_8_',num2str(seedd),'.mat'));
            totPow1(nnn,xxx,1) = TotalPower(1,:);

        elseif xxx == 3
            load(strcat('J_SE2_K12_sim_8_',num2str(seedd),'.mat'));
            totPow1(nnn,xxx,1) = TotalPower(1,:);

        else
            load(strcat('J_SE2_K16_sim_8_',num2str(seedd),'.mat'));
            totPow1(nnn,xxx,1) = TotalPower(1,2);

        end
        if isnan(totPow1(nnn,xxx,1)) || (totPow1(nnn,xxx,1)==0)
        else
            if xxx == 4
                xx = Axx(:,:,2);
                zz = Azz(:,2);
                activeAPs = find(zz>0);
                dd = Add(:,2);
                ll = All(:,2);
                rho = Arho(:,2);
            else
                xx = Axx;
                zz = Azz;
                activeAPs = find(zz>0);
                dd = Add;
                ll = All;
                rho = Arho;
            end

                
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



                

                NumberOfDUs = max((1:W)*dd, NumberOfDUs);
                NumberOfDUs = max(NumberOfDUs,NumberOfLCs);


                SE0 = 2;
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
                %NumberOfDUs2 = max(NumberOfDUs,12*ceil((ZLP*3)/GOPSmax));
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
                totPow3(nnn,xxx,1) = totPow3(nnn,xxx,1) + (tempor1 + tempor2 + tempor3 + tempor4 + tempor5)/5;
                totPow4(nnn,xxx,1) = totPow4(nnn,xxx,1) + (tempor1 + tempor2 + tempor3 + tempor4b + tempor5)/5;

                if xxx == 2
                    totPow1b(nnn,1,1) = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
                    totPow1b(nnn,2,1) = PONU*sum(zz) + POLT*(1:W)*ll/sigmaCool;
                    totPow1b(nnn,3,1) = Pdisp;
                    totPow1b(nnn,4,1) = Pproc0*(1:W)*dd/sigmaCool;
                    totPow1b(nnn,5,1) = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
                    + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
                    + DeltaProc/GOPSmax*FLP/sigmaCool;

                
                    totPow3b(nnn,1,1) = totPow3b(nnn,1,1)+tempor1/5;
                    totPow3b(nnn,2,1) = totPow3b(nnn,2,1)+tempor2/5;
                    totPow3b(nnn,3,1) = totPow3b(nnn,3,1)+tempor3/5;
                    totPow3b(nnn,4,1) = totPow3b(nnn,4,1)+tempor4/5;
                    totPow3b(nnn,5,1) = totPow3b(nnn,5,1)+tempor5/5;

                    totPow4b(nnn,1,1) = totPow4b(nnn,1,1)+tempor1/5;
                    totPow4b(nnn,2,1) = totPow4b(nnn,2,1)+tempor2/5;
                    totPow4b(nnn,3,1) = totPow4b(nnn,3,1)+tempor3/5;
                    totPow4b(nnn,4,1) = totPow4b(nnn,4,1)+tempor4b/5;
                    totPow4b(nnn,5,1) = totPow4b(nnn,5,1)+tempor5/5;
                end
            end
        end
    end

    
  


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

for nnn = 1:nbrOfRuns
    seedd = scen(nnn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for xxx = 1:Userlength
         if xxx == 1
            load(strcat('J_SE2_K4_sim_72_',num2str(seedd),'.mat'));
            totPow1(nnn,xxx,2) = TotalPower(1,:);

        elseif xxx == 2
            load(strcat('J_SE2_K8_sim_72_',num2str(seedd),'.mat'));
            totPow1(nnn,xxx,2) = TotalPower(1,:);

        elseif xxx == 3
            load(strcat('J_SE2_K12_sim_72_',num2str(seedd),'.mat'));
            totPow1(nnn,xxx,2) = TotalPower(1,:);

        else
            load(strcat('J_SE2_K16_sim_72_',num2str(seedd),'.mat'));
            totPow1(nnn,xxx,2) = TotalPower(1,2);

        end
        




        if isnan(totPow1(nnn,xxx,2)) || (totPow1(nnn,xxx,2)==0)
        else
            if xxx == 4
                xx = Axx(:,:,2);
                zz = Azz(:,2);
                activeAPs = find(zz>0);
                dd = Add(:,2);
                ll = All(:,2);
                rho = Arho(:,2);
            else
                xx = Axx;
                zz = Azz;
                activeAPs = find(zz>0);
                dd = Add;
                ll = All;
                rho = Arho;
            end

           
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



                NumberOfDUs = max((1:W)*dd,NumberOfDUs);
                NumberOfDUs = max(NumberOfDUs,NumberOfLCs);


                SE0 = 2;
                CprecodingAP = Nused/(Ts*tau_c*10^9)*...
                    ((8*N*tau_p+8*N^2)*tau_p + (4*N^2+4*N)*tau_p + 8*(N^3-N)/3);
                CotherAP = ((Nbits/16)^(1.2))*(1.3*N) ...
                    + ((Nbits/16)^(0.2))*(2.7*sqrt(N));
                SLP = Cfilter + CDFT;
                ZLP = CprecodingAP + CotherAP;

                CprecodingUE = Nused*(tau_c-tau_p)/(Ts*tau_c*10^9)*8*N ...
                    + Nused/(Ts*tau_c*10^9)*8*N ...
                    + Nused/(Ts*tau_c*10^9)*(8*N^2)*2;

                FLP = ((Nbits/16)^(1.2))*( 1.3*((SE0/6)^(1.5))*K )...
                    + ((Nbits/16)^(1.2))*( 1.3*SE0/6*K )...
                    + 8*SE0/6*K;

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
                totPow3(nnn,xxx,2) = totPow3(nnn,xxx,2)+(tempor1 + tempor2 + tempor3 + tempor4 + tempor5 ...
                    + tempor4c + tempor5b)/5;
                totPow4(nnn,xxx,2) = totPow4(nnn,xxx,2)+(tempor1 + tempor2 + tempor3 + tempor4b + tempor5 ...
                    + tempor4c + tempor5b)/5;




                if xxx == 2

                    totPow1b(nnn,1,2) = PAP0*sum(zz) + DeltaTr*quad_form(rho,eye(K*L))/10;
                    totPow1b(nnn,2,2) = PONU*sum(zz) + POLT*(1:W)*ll/sigmaCool;
                    totPow1b(nnn,3,2) = Pdisp;
                    totPow1b(nnn,4,2) = Pproc0*(1:W)*dd/sigmaCool;
                    totPow1b(nnn,5,2) = DeltaProc/GOPSmax*ZLP/sigmaCool*sum(zz) ...
                        + DeltaProc/GOPSmax*XLP*sum(sum(xx))/sigmaCool ...
                        + DeltaProc/GOPSmax*FLP/sigmaCool;
                    totPow1b(nnn,6,2) = PRUproc0*sum(zz)+DeltaRUProc/GOPSmax*SLP*sum(zz);

                    totPow3b(nnn,1,2) = totPow3b(nnn,1,2)+tempor1/5;
                    totPow3b(nnn,2,2) = totPow3b(nnn,2,2)+tempor2/5;
                    totPow3b(nnn,3,2) = totPow3b(nnn,3,2)+tempor3/5;
                    totPow3b(nnn,4,2) = totPow3b(nnn,4,2)+tempor4/5;
                    totPow3b(nnn,5,2) = totPow3b(nnn,5,2)+tempor5/5;
                    totPow3b(nnn,6,2) = totPow3b(nnn,6,2)+tempor4c/5 + tempor5b/5;


                    totPow4b(nnn,1,2) = totPow4b(nnn,1,2)+tempor1/5;
                    totPow4b(nnn,2,2) = totPow4b(nnn,2,2)+tempor2/5;
                    totPow4b(nnn,3,2) = totPow4b(nnn,3,2)+tempor3/5;
                    totPow4b(nnn,4,2) = totPow4b(nnn,4,2)+tempor4b/5;
                    totPow4b(nnn,5,2) = totPow4b(nnn,5,2)+tempor5/5;
                    totPow4b(nnn,6,2) = totPow4b(nnn,6,2)+tempor4c/5 + tempor5b/5;
                end
            end
        end
    end
    

    
    

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555555
tugfe1 = zeros(Userlength,2);
tugfe3 = zeros(Userlength,2);
tugfe4 = zeros(Userlength,2);

tugfe1b = zeros(2,6);
tugfe3b = zeros(2,6);
tugfe4b = zeros(2,6);


t1 = nbrOfRuns*ones(Userlength,2);
for x = 1:Userlength

    for n = 1:nbrOfRuns
        if isnan(totPow1(n,x,1)) || (totPow1(n,x,1)==0)
            t1(x,1) = t1(x,1)-1;
        else
            tugfe1(x,1) = tugfe1(x,1) + totPow1(n,x,1);
            tugfe3(x,1) = tugfe3(x,1) + totPow3(n,x,1);
            tugfe4(x,1) = tugfe4(x,1) + totPow4(n,x,1);


            if x == 2
                tugfe1b(1,:) = tugfe1b(1,:) + totPow1b(n,:,1);
                tugfe3b(1,:) = tugfe3b(1,:) + totPow3b(n,:,1);
                tugfe4b(1,:) = tugfe4b(1,:) + totPow4b(n,:,1);


            end
        end

        if isnan(totPow1(n,x,2)) || (totPow1(n,x,2)==0)
            t1(x,2) = t1(x,2)-1;
        else
            tugfe1(x,2) = tugfe1(x,2) + totPow1(n,x,2);
            tugfe3(x,2) = tugfe3(x,2) + totPow3(n,x,2);
            tugfe4(x,2) = tugfe4(x,2) + totPow4(n,x,2);


            if x == 2
                tugfe1b(2,:) = tugfe1b(2,:) + totPow1b(n,:,2);
                tugfe3b(2,:) = tugfe3b(2,:) + totPow3b(n,:,2);
                tugfe4b(2,:) = tugfe4b(2,:) + totPow4b(n,:,2);


            end
        end


    end
end
%

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

set(groot,'defaultLegendInterpreter','latex');

figure


plot(4:4:16,tugfe1(1:4,1)./t1(:,1))
hold on
plot(4:4:16,tugfe3(1:4,1)./t1(:,1))
plot(4:4:16,tugfe4(1:4,1)./t1(:,1))

plot(4:4:16,tugfe1(1:4,2)./t1(:,2))
plot(4:4:16,tugfe3(1:4,2)./t1(:,2))
plot(4:4:16,tugfe4(1:4,2)./t1(:,2))
legend({'FS-8, end-to-end','FS-8, local coordination','FS-8, radio-only','FS-7.2, end-to-end','FS-7.2, local coordination','FS-7.2, radio-only'});



Y = zeros(6,6);
Y(1,:) = tugfe1b(1,:)/t1(2,1);
Y(2,:) = tugfe3b(1,:)/t1(2,1);
Y(3,:) = tugfe4b(1,:)/t1(2,1);

Y(4,:) = tugfe1b(2,:)/t1(2,2);
Y(5,:) = tugfe3b(2,:)/t1(2,2);
Y(6,:) = tugfe4b(2,:)/t1(2,2);

X = categorical({'FS-8, end-to-end','FS-8, local coordination','FS-8, radio-only','FS-7.2, end-to-end','FS-7.2, local coordination','FS-7.2, radio-only'});
X = reordercats(X,{'FS-8, end-to-end','FS-8, local coordination','FS-8, radio-only','FS-7.2, end-to-end','FS-7.2, local coordination','FS-7.2, radio-only'});
figure
bar(X,Y);
legend({'RU', 'Fronthaul', 'Controller', 'GPP idle', 'GPP baseband', 'RU baseband'},'Interpreter','Latex','Location','Best');
ylabel('Power consumption (Watts)','Interpreter','Latex');
ylim([0 300])
set(gca,'fontsize',18);
