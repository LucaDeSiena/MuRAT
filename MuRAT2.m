% PROGRAM MuRAT for Multi-resolution attenuation tomography analysis
%
% PROGRAM AIM: 3D attenuation tomography - the coda-normalization method.
%
% Author:  L. De Siena.
%
%           January 2017 
% Changes necessary:
% 1. Work with SAC files
% 2. Include ray-tracing
% 3. Include commands files in Header
% 4. Change for varying node spacing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./Utilities_Matlab')

clear
close all

%==========================================================================
% The following prompts (file commands.m in v1.0) are included TO BE 
% EDITED before running!
%==========================================================================

%COMMANDS for MuRAT

% Choose between P- (2) and S-direct (3) waves
PorS = 2; 

% Work with 1 vertical (1) or 2 horizontal (2) recordings, or with the
% three components(3)
compon = 1;

% Step of the inversion - set to 2 in case of double grid, 1 otherwise
nstep= 1;

% Name of the file containing the 1st velocity model - according to PorS
modv = load('modv.txt','%f');%edit
v0= mean(modv(:,4));
resol2 = abs(modv(2,3)-modv(1,3))/2;

% Name of the file containing the picking of origin time, P, and S waves
tempi = load('tempi.txt','%f');%unnecessary if with SAC files

% Length of the window used to measure P- or S- wave energy in
% seconds
fin = 1;%edit

% Central frequency
cf = 18;%edit

% Lapse time (from the origin time of the event) must be >> 2*tS
tC = 15;%edit

% Total coda window length
tW = 10;%edit

%UTM coordinates of the origin of the velocity model - MSH
originWE=538311;
originSN=5092338;
originz=3350;

%Rays are measured in meters. If the ray length is already in km set um =1;
um = 1000;

%In order to check the images, one can change the smoothing parameter
%(set smoot=1)
smoot=0;

%The minimum coda-to-noise energy ratio for the weighted inversion
soilnoise = 1.4;%edit

%Size anomaly for testing: twice(2)/four(4) times
sizea=4;

if nstep == 2
%     Choose to select just part of the data affected by anomalous
%     behavior of the energy ratios with travel times
%    tA1=0;
%    tA2=2;
%     Otherwise impose a second grid in the regions showing this behavior 
     modv2 = load('modv2.txt','%f');
     resol22 = abs(modv2(2,3)-modv2(1,3))/2;
     smoot2=0;
end

display('STEP 1 - Measures related to the rays and inversion matrix')

%==========================================================================
% In the first analysis the ray-path contained in the folder "rays" are
% analyzed. The files are ascii files containing the cartesian coordinates
% of the rays in the reference system of the velocity model, starting at
% the source and with positive depth. The velocity model is a x,y,z grid,
% possibly with different steps in the three directions.
%==========================================================================

%Function creating the parameters for the inversion matrix
[lunparz, blocchi, luntot,inblocchi,sb] =...
    segments2(modv,originz,um);
save raysvar.mat lunparz blocchi luntot inblocchi sb

load raysvar.mat
signal = zeros(length(luntot),1); % The direct wave energies
coda = zeros(length(luntot),1); % The coda wave energies
rapsp = zeros(length(luntot),1); % The spectral ratios
rapspcn = zeros(length(luntot),1); % The coda versus noise ratios
noise = zeros(length(luntot),1); % The noise energies
index1=0;

display('STEP 2 - Seismic attributes')
warning off all
index = 0;

%==========================================================================
% A loop to measure the energy ratios for each seismic trace located in
% the folder "traces" with the coda-normalization method.
% The seismograms may be the output of both single and three-component
% seismic stations. In case of more than one component the rays must be
% ordered in the folder starting with the WE component, then SN and finally
% the Z component for each ray.
% 
% The program accepts 2 columns ascii data, where the first column
% is the time and the second column is displacement, velocity, or
% accelleration.
%==========================================================================

fin1=1; %time window for noise
sn=5; %start time to consider noise - before P-arrival
Wn = ([cf-cf/3 cf+cf/3]/50); %frequency band
[z,p,k] = butter(6,Wn,'bandpass'); %butter filter
[sos,g] = zp2sos(z,p,k);      % Convert to SOS form
Hd = dfilt.df2tsos(sos,g);   % Create a dfilt object   

% Name of the file containing the seismograms - in folder traces
list=dir('./traces');  %get name of files in current folder traces
isfile=~[list.isdir]; %determine index of files vs folders
filenames={list(isfile).name}; %create cell array of file names
%creates a list of the waveforms and calculates their number
lista = filenames';
for i = 1:length(lista)
    lista{i,1}=cat(2,'./traces/',lista{i});
end
ll = length(lista);

% %For loop on the waveforms
% for i = 1:ll
%     index = index+1;
%     waveform = index;
%     display(waveform);% number of the waveform
%     si = lista{i,1};% start reading waveform
%     fileID = fopen(si);
%     [sis] = textscan(fileID,'%f %f');% use textscan to get the two columns
%     fclose(fileID);
%     tempis = sis{1};% first column is time
%     sisma = sis{2};% second is measurement
%     
%     srate = (tempis(5)-tempis(4))^-1;% sampling frequency
%     int = fin*srate;% number of sample for the chosen window
%     intn = fin1*srate;% number of sample for the noise
%     t00=tempis(1); %starting time of the waveform
%     
%     %direct wave
%     cursor1 = (tempi(i,PorS)-t00)*srate;%starting-sample direct window
%     cursor2 = cursor1 + int-1; %end-sample of the window
%     tsisma = sisma(cursor1:cursor2);% tapering direct wave
%     fsisma = filter(Hd,tsisma);% filtering direct wave
%     spamp = smooth(rms(fsisma),10);% direct energy
%     %choose if P or S and measure noise
%     if PorS == 2
%         cursorn1 = (tempi(i,2)-t00-sn)*srate; %starting sample for noise
%     elseif PorS == 3
%         cursorn1 = (t00)*srate+1; %starting sample for noise
%     end
%     cursorn2 = cursorn1 + intn-1; % end-sample of the noise-window
%     tsisman = sisma(cursorn1:cursorn2);%tapering noise
%     fsisman = filter(Hd,tsisman);%filtering noise
%     spampn = smooth(rms(fsisman),10); %noise energy
%     
%     %coda average energy
%     maxtime = tW;
%     cursorc1 = (tempi(i,1)-t00+tC-1)*srate;%coda start sample
%     cursorc2 = cursorc1 + maxtime*srate-1;%coda end sample
%     if cursorc2 > length(sisma)
%         lsis=length(sisma);
%         tsismac = sisma(cursorc1:lsis);%in case coda is too long
%     else
%         tsismac = sisma(cursorc1:cursorc2);%tapering
%     end
%     fsismac = filter(Hd,tsismac); %filter coda
%     spampc = smooth(rms(fsismac),10);%energy coda
%     
%     rmes = spamp/spampc;%spectral ratios signal-coda
%     rmcn = spampc/spampn;%spectral ratios coda-noise
%     
%     signal(index,1) = spamp; %direct-wave energies
%     coda(index,1) = spampc; %coda-wave energies
%     rapsp(index,1) = rmes; %spectral ratios
%     rapspcn(index,1)= rmcn; %coda-to-noise ratios
%     noise(index,1)= spampn; %noise energies
% end
% %save quantities for the data vector
% save energies.mat signal coda rapsp rapspcn noise
% 
%==========================================================================
% The average quality factor is obtained from the data.
%==========================================================================

display('STEP 3 - Average geometrical spreading factor and Q')
load energies.mat
icomp = 0;

% 1 component (tipically only vertical)
if compon ==  1
    lsig=ll;
    signal1 = signal;
    coda1 = coda;
    rapsp1 = rapsp;
    rapspcn1 = rapspcn;
% 2 components (tipically the two horizontals)
elseif compon ==  2
    lsig = ll/2;
    signal1 = zeros(lsig,1); % The average direct wave energies
    coda1 = zeros(lsig,1); % The average coda wave energies
    rapsp1 = zeros(lsig,1); % The average spectral ratios
    rapspcn1 = zeros(lsig,1); % The average coda versus noise ratios
    for i = 1:2:(ll-1)
        icomp = icomp+1;
        signal1(icomp,1) = (signal(i,1)+signal(i+1))/2;
        coda1(icomp,1) = (coda(i)+coda(i+1))/2;
        rapsp1(icomp,1) = (rapsp(i)+rapsp(i+1))/2;
        rapspcn1(icomp,1) = (rapspcn(i)+rapspcn(i+1))/2;
    end
% 3 components (WE, SN, Z)
elseif compon == 3
    lsig = ll/3;
    signal1 = zeros(lsig,1); % The average direct wave energies
    coda1 = zeros(lsig,1); % The average coda wave energies
    rapsp1 = zeros(lsig,1); % The average spectral ratios
    rapspcn1 = zeros(lsig,1); % The average coda versus noise ratios
    for i = 1:3:(ll-2)
        icomp = icomp+1;
        signal1(icomp,1) = ((signal(i)+signal(i+1))/2 + signal(i+2))/2;
        coda1(icomp,1) = ((coda(i)+coda(i+1))/2 + coda(i+2))/2;
        rapsp1(icomp,1) = ((rapsp(i)+rapsp(i+1))/2 + rapsp(i+2))/2;
        rapspcn1(icomp,1) = ((rapspcn(i)+rapspcn(i+1))/2 + rapspcn(i+2))/2;
    end
end

mcn = find(rapspcn1<soilnoise);% set the weigtht
% weighting
time0=tempi(1:compon:end-compon+1,PorS)-tempi(1:compon:end-compon+1,1);
W1=rapspcn1<soilnoise;
W2=rapsp1>1000;
W3=rapsp1<1.4;
W4=W3+W2+W1;
W=diag(1./(W4+1));% weights
d0=log(rapsp1)/pi/cf; %data of the inverse problem

dW=W*d0;
K2 = ones(lsig,1); %the K in the CN method
G1=K2; %matrix creation
G1(:,2)=-log(luntot')/pi/cf;
%time0=sum(lunparz.*sb)';
G1(:,3)=-time0;

G=W*G1;
%The 3 parameters are the constant, the geometrical spreading, and the
%average Q and they are contained in constQmean.
constQmean(1,1:3)=lsqlin(G,dW(:,1));% damped least square inversion
cova = (G'*G)^(-1)*G'*cov(dW)*G*(G'*G)^(-1); %covariance matrix
er = sqrt(diag(cova)); %error from the covariance matrix
constQmean(2,:)=er;

%save the three average measurements
save -ascii Qmean.txt constQmean


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS RELATED TO THE CN METHOD: The large cyan dots represent the left-hand side of the
% CN method linear equation. The red line is the fit of the three
% parameters. The black lines are the uncertainties.
% In the lower panel the coda-to noise ratios, to see the effe ct of noise
% on the data.

% Right hand side of the CN linear equation
x = constQmean(1,1)...
    -constQmean(1,2)*log(luntot')/pi/cf...
    -time0*constQmean(1,3);
xer1 = (constQmean(1,1)-er(1))...
    -(constQmean(1,2)-er(2))*log(luntot')/pi/cf...
    -time0*(constQmean(1,3)-er(3));
xer2 = (constQmean(1,1)+er(1))...
    -(constQmean(1,2)+er(2))*log(luntot')/pi/cf...
    -time0*(constQmean(1,3)+er(3));

%Plot of the left and right hand sides of the CN equation. The average
%inverse Q and geometrical spreading coefficient are schown above each
%panel
figure
subplot(2,1,1)
plot(time0,d0, 'o',...
    'MarkerSize',6,'MarkerEdgeColor',[0 0.6 1],'MarkerFaceColor',[0 0.6 1])
hold on
plot(time0,x,'r.')
plot(time0,xer1,'k.')
plot(time0,xer2,'k.')
hold off
xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Logarithm of the direct-to-coda energy ratio','FontSize',12,...
    'FontWeight','bold','Color','k')

title(['Average inverse quality factor = ', num2str(constQmean(1,3)),...
    ' +/- ', num2str(er(3)), ' and geometrical spreading = ',...
    num2str(constQmean(1,2)),' +/- ', num2str(constQmean(2,2))],...
    'FontSize',12,'FontWeight','bold','Color','k')
h_legend=legend('Direct-to-coda energy ratios','Inverse average Q',...
    'Inverse Q spreading-related uncertainties');
set(gca,'XTick',0:round(max(time0))) ;
set(gca,'YTick',min(log(rapsp1)+constQmean(1,2)*log(luntot'))/pi/cf:...
    (max(log(rapsp1)+constQmean(1,2)*log(luntot'))/pi/cf-...
    min(log(rapsp1)+constQmean(1,2)*log(luntot'))/pi/cf)/5:...
    max(log(rapsp1)+constQmean(1,2)*log(luntot'))/pi/cf) ;
set(gca,'FontSize',10,'FontWeight','bold');
set(h_legend,'FontSize',10,'FontWeight','bold');

% %Plot of coda-to-noise energy ratios. The residual geometrical
% spreading should be around zero, or coherent phases could be
% included in the coda-energy; it also provides the average coda-to-noise
% ratio
subplot(2,1,2)
plot(time0,(log(rapspcn1)),'o','MarkerSize',6,'MarkerEdgeColor',...
    [.2 .8 0],'MarkerFaceColor',[.2 .8 0])
gc=polyfit(time0,log(rapspcn1),1);
mc = mean(rapspcn1);
mno=length(mcn);
perno= mno/length(rapspcn1)*100;
hold on
plot(time0(mcn),log(rapspcn1(mcn)), 'o','MarkerSize',6,...
    'MarkerEdgeColor',[1 .5 0],'MarkerFaceColor',[1 .5 0])
hold off
xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Logarithm of the coda-to-noise energy ratio','FontSize',12,...
    'FontWeight','bold','Color','k')
title(['Residual coda geometrical spreading = ',num2str(gc(1)), ...
    ' while the percentage of measures < ',num2str(soilnoise), 'is ',...
    num2str(perno), '%'],'FontSize',12,'FontWeight','bold','Color','k');
h_legend=legend('Coda-to-noise energy ratios','Ratios below soil');
set(gca,'XTick',0:round(max(time0))) ;
set(gca,'YTick',min(log(rapspcn1)):...
    (max(log(rapspcn1))-min(log(rapspcn1)))/5:...
    max(log(rapspcn1)));
set(gca,'FontSize',10,'FontWeight','bold');
set(h_legend,'FontSize',10,'FontWeight','bold');

text(max(time0), min(log(rapspcn1)),...
[' Lapse time = ',num2str(tC),' s, coda window = ', num2str(tW), ' s'],...
'VerticalAlignment','bottom',...
'HorizontalAlignment','right',...
'FontSize',12)

% END PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF THE FIRST TOMOGRAPHIC INVERSIOn
display('STEP 4 - Tomographic results')

%data creation
d1 = d0 - constQmean(1,1) + constQmean(1,2)*log(luntot')/pi/cf...
    + time0*constQmean(1,3);

%INVERSION MATRIX: A

% if the ray crosses  a block not solved by the inversion, the
% block is characterized by the average quality factor, and the data
% vector is updated.
A = zeros(length(time0),length(inblocchi(:,1)));
for j = 1:length(sb(1,:))
    for i = 2:length(sb(:,1))
        bl = blocchi(i,j);
        slo = sb(i,j);
        l = lunparz(i,j);
        fbl = find(inblocchi(:,5)==bl);
        if isempty(fbl)==0
            A(j,fbl)=-l*slo;% inversion matrix element in sec
        end
        continue
    end
    continue
end
A1=A;
A=W*A1;

% NEW DATA VECTOR

dW1=W*d1;

% tikhonov inversion - by using the programs in HANSEN et al. 1994
[U,S,V]=svd(A);
figure
tik0_reg=l_curve(U,diag(S),dW1,'Tikh');

% piocard plot
figure
picard(U,diag(S),dW1);

%sets the smoothing parameter - only if user defined
if smoot ==1
    tik0_reg=input('Your personal smoothing parameter ');
end

%results
mtik0=tikhonov(U,diag(S),V,dW1,tik0_reg);
mQ1=[inblocchi(:,1:3) mtik0];

%save result in matlab format
%save mQ1.mat mQ1

%Simple resolution test: produces two files which can be compared.

display('STEP 5 - Testing')
%FIRST TEST: CHECKERBOARD
%INPUT: Checkerboard structure - node spacing of the anomalies
% doubled with respect to the grid.

%Depth
passox = find(modv(:,1)~=modv(1,1),1,'first')-1;
passoy = find(modv(:,2)~=modv(1,2),1,'first')-1;

sizea2=2*sizea;
sizeap=sizea*passoy;
sizeap1=(sizea+1)*passoy;

for k=1:sizea
    modv(k:sizea2:passoy-sizea+k,6)=.02;
    modv(sizea+k:sizea2:passoy-sizea+k,6)=0.001;
end
for k=1:sizea-1
    modv(k*passoy+1:(k+1)*passoy,6)=modv(1:passoy,6);
end
for k=1:sizea
    modv(sizeap+k:sizea2:sizeap1-sizea+k,6)=.001;
    modv(sizeap+sizea+k:sizea2:sizeap1-sizea+k,6)=0.02;
end

py4  = 2*sizeap;
for k=1:sizea-1
    modv((sizea+k)*passoy+1:(sizea+k+1)*passoy,6)=modv(sizeap+1:sizeap1,6);
end
z = (passox-mod(passox,py4))/py4;
for i = 1:(z-1)
    modv(i*py4+1:(i+1)*py4,6)=modv(1:py4,6);
end
if ~isequal(mod(passox,py4),0)
    modv(z*py4+1:z*py4+mod(passox,py4),6)= modv(1:mod(passox,py4),6);
end

%Along y
sizeapx=sizea*passox;

for k=1:sizea-1
    modv(k*passox+1:(k+1)*passox,6)=modv(1:passox,6);
end

for k = 1:sizeapx
    if modv(k,6)==.02
        modv(sizeapx+k,6)=.001;
    elseif modv(k,6)==.001;
        modv(sizeapx+k,6)=.02;
    end
end

%Along x
px4  = 2*sizea*passox;

z2= (length(modv(:,1))-mod(length(modv(:,1)),px4))/px4;
for i = 1:(z2-1)
    modv(i*px4+1:(i+1)*px4,6)=modv(1:px4,6);
end
if ~isequal(mod(passox,py4),0)
    modv(z*px4+1:z*px4+mod(length(modv(:,1)),px4),6)=...
        modv(1:mod(length(modv(:,1)),px4),6);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS OF THE FIRST TOMOGRAPHIC INVERSION + tests
% In the original UTM reference system we put the quality factors in
% the center of each block. The average quality factor characterizes the
% blocks not solved in the inversion.
Q3D=modv(:,1:3);
Q3D(:,4)=constQmean(1,3);
% create a file of input
in1=zeros(length(mtik0),1);
input01 = zeros(length(mtik0),1);

for i = 1:length(in1)
    in1(i)=find(modv(:,1)==inblocchi(i,1)  & modv(:,2)==inblocchi(i,2)...
        & modv(:,3)==inblocchi(i,3));
    Q3D(in1(i),4)=mtik0(i)+constQmean(1,3);
    input01(i)=modv(in1(i),6);
end

% result in the original reference system
Q3D(:,1)=originWE+modv(:,1)+resol2;
Q3D(:,2)=originSN+modv(:,2)+resol2;
Q3D(:,3)=modv(:,3)+resol2;

% input save
Qin=[Q3D(:,1:3) modv(:,6)];

% output file and save
re = A1*input01;
Qou(:,1:3)=Q3D(:,1:3);
Qou(:,4)=0;
mcheck=tikhonov(U,diag(S),V,W*(re-mean(input01)),tik0_reg);
Qou(in1,4)=mcheck;

% AS SECOND TEST WE ALSO OUTPUT THE DIAGONAL OF THE RESOLUTION MATRIX

% using the filter functions
sS=size(S);
fil_reg = fil_fac(diag(S),tik0_reg);
if sS(2)>sS(1)
    fil_reg(end+1:sS(2),1)=0;
end
R=V*diag(fil_reg)*V';
dR =diag(R);
Qou(in1,5)=dR;

% AS THIRD AND FINAL TEST, THE USER CAN SET SYNTHETIC ANOMALIES.
%First is with the entire set of results.
%Second is user-defined.
lmQ1=length(mQ1(:,4));
syQ0=constQmean(1,1)-constQmean(1,2)*log(luntot')/pi/cf-...
    time0.*constQmean(1,3)+A*mQ1(:,4);

K2 = ones(lsig,1); %the K in the CN method
G1=K2; %matrix creation
G1(:,2)=-log(luntot')/pi/cf;
%time0=sum(lunparz.*sb)';
G1(:,3)=-time0;

dsW=W*syQ0;
G=W*G1;

%The 3 parameters are the constant, the geometrical spreading, and the
%average Q and they are contained in synthQmean.
synthQmean(1,1:3)=lsqlin(G,dsW(:,1));% damped least square inversion
cova = (G'*G)^(-1)*G'*cov(dsW)*G*(G'*G)^(-1); %covariance matrix
er = sqrt(diag(cova)); %error from the covariance matrix
synthQmean(2,:)=er;

%data creation
syQ1 = W*(syQ0 - synthQmean(1,1) + synthQmean(1,2)*log(luntot')/pi/cf...
    + time0*synthQmean(1,3));

[syU,syS,syV]=svd(A);
figure
syk0_reg=l_curve(syU,diag(syS),syQ1,'Tikh');

% picard plot
figure
picard(U,diag(S),syQ1);

%sets the smoothing parameter - only if user defined
%results
diftik0=tik0_reg-syk0_reg;
display(diftik0)
sytik0=tikhonov(syU,diag(syS),syV,syQ1,tik0_reg);
Qou(:,6)=synthQmean(1,3);
Qou(in1,6)=synthQmean(1,3)+sytik0;
Q3D(:,5)=Qin(:,4);
Q3D(:,6:8)=Qou(:,4:6);

%Create a medium with average Q and same results at depth
Qplus=find(inblocchi(:,3)>-3500);
Qminus=find(inblocchi(:,3)<-3500);
syQ2=mQ1;
syQ2(Qplus,4)=mean(mQ1(:,4));
syQ2(Qminus,4)=0.02;

d0=0.1-0.3*log(luntot')/pi/cf-time0*0.02;
dW=W*d0;
K2 = ones(lsig,1); %the K in the CN method
G1=K2; %matrix creation
G1(:,2)=-log(luntot')/pi/cf;
%time0=sum(lunparz.*sb)';
G1(:,3)=-time0;

G=W*G1;
%The 3 parameters are the constant, the geometrical spreading, and the
%average Q and they are contained in constQmean.
constQmean(1,1:3)=lsqlin(G,dW(:,1));% damped least square inversion
cova = (G'*G)^(-1)*G'*cov(dW)*G*(G'*G)^(-1); %covariance matrix
er = sqrt(diag(cova)); %error from the covariance matrix
constQmean(2,:)=er;


%syQ2(:,4)=0.02;
syQ3=0.1-0.3*log(luntot')/pi/cf-time0*0.02;

dsW2=W*syQ3;

synthQ2mean(1,1:3)=lsqlin(G,dsW2(:,1));% damped least square inversion
synthQ2mean(2,:)=er;

syQ4=W*(syQ3 - synthQ2mean(1,1) + synthQ2mean(1,2)*log(luntot')/pi/cf...
    + time0*synthQ2mean(1,3));

syk0_reg2=l_curve(syU,diag(syS),syQ4,'Tikh');
sytik02=tikhonov(syU,diag(syS),syV,syQ4,tik0_reg);

Qou(:,7)=constQmean(1,3);
Qou(:,8)=synthQ2mean(1,3);
Qou(in1,7)=syQ2(:,4);
Qou(in1,8)=synthQ2mean(1,3)+sytik02;
Q3D(:,9:10)=Qou(:,7:8);

% save
save -ascii Q3D.txt Q3D

