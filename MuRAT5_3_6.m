% PROGRAM MuRAT for Multi-resolution attenuation tomography analysis
%
% PROGRAM FOR: 3D direct- and 2D peak delay and coda-wave (with and without
% kernels) attenuation tomography
%
% Author:  L. De Siena, September 2018 
% 
% Some reference papers:
%
% Peak-delays and 2D coda attenuation: Takahashi et al. 2007; Calvet et al.
% 2013b; De Siena et al. 2016
%
% Kernel-based 2D coda attenuation: Del Pezzo et al. 2016; De Siena et al.
% 2017 (GRL)
%
% Direct wave coda-normalized attenuation: Del Pezzo et al. 2006; De Siena
% et al. 2010 (JGR); De Siena et al. 2014a (MuRAT, JVGR) and b (JGR)
% 
% The CN method has changed in that the Qc measurements are used to
% calculate the denominator of the equation, assuming coda waves are
% comprised of either body or surface waves.
%
% BEFORE RUNNING DEFINE AN INPUT FILE FROM TEMPLATES
%
% Necessary edits:
% 
% 1. Add MLTWA - Simona
% 2. Set Qc kernel inversion - Luca, Panayiota
% 3. Set diffusive non-linear inversion- Panayiota
% 4. Moving the program to python - Luca
% 5. Writing documentation - Luca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS
addpath('/Users/lucadesiena/Documents/Utilities_Matlab')
clear
close all
clc

disp('Input Section')
tic
% Write the name of the input file in the command window
input('Name of .m input file ')

% The following prompts will measure peak-delays (pd) ad Qc depending on
% the input files.
%% RAY TRACING METHOD: RAY BENDING - for evst=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With option evst=1, we use two external files for event and station
% locations. The SAC file name must have a specified format for recognition
           
if evst==1 
%     tic
    disp('Inversion Section for evst=1')
    
%Predefine the inversion matrix for pd and Qc imaging
    
    Ac=zeros(lls,lxy);
    Apd=zeros(lls,lxy);
    D=zeros(lls,1);
    evestaz=zeros(lls,4);

% This loop creates both rays and elements of Apd and Ac - for peak delays
% and Qc imaging in case of files of events and stations outside SAC haeder

    indexray=0; % Increases every time a ray is created.
    
    %Sets figure for rays
    rays=figure('Name','Rays','NumberTitle','off','visible',visib);
     
    if pa==3
        lunparz = zeros(100,lls);
        blocchi = zeros(100,lls);
        sb = zeros(100,lls);
        luntot = zeros(lls,1);
        
        if createrays ==1
            ma1=zeros(100,5,lls);
        end
    end
    % Loop over events in even.txt file
    for ii=1:numev
        ray(1,1)=even(ii,1);
        ray(2,1)=even(ii,2);
        ray(3,1)=-even(ii,3);
        namee1=nameeven{ii};
        namee=namee1(1:12);
        
        % Loop over stations in staz.txt file
        for ir=1:numst
            namest1=namesta{ir};
            namest=namest1(1:3);
            c1=cat(2,namee,namest);
            
            % Pattern recognition from the name of the SAC file
            if find(strncmp(lista(:,1),c1,length(c1)))>0
                indexray=indexray+1;
                if isequal(mod(indexray,100),0)
                    disp(['Ray number is ', num2str(indexray)])
                end
                ray(1,2) = staz(ir,1);
                ray(2,2) = staz(ir,2);
                ray(3,2) = -staz(ir,3);
                xx=XY(:,1);
                yy=XY(:,2);    
                miix=min(ray(1,1),ray(1,2));
                maax=max(ray(1,1),ray(1,2));
                miiy=min(ray(2,1),ray(2,2));
                maay=max(ray(2,1),ray(2,2));
                fxy=find(xx>miix & xx<maax & yy>miiy & yy<maay);
                             
                Apd(indexray,fxy)=1;% Creates pd matrix
                sst = [ray(1,1) ray(2,1) ray(1,2) ray(2,2)];
                evestaz(indexray,1:4)=sst;% Creates file of selected coords
                
                % Epicantral distance
                D(indexray,1)=sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2)/1000;
                
                if pa==1
                    
                    % Same pd and Qc inversion matrix
                    Ac=Apd;
   
                elseif pa>1
                    
                    % Operations to calculate normalised kernels
                    D1=sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2);
                    deltaxy=0.2;
                    F1=1/2/pi/deltaxy^2/D1^2;
                    F2=1/deltaxy^2/D1^2;
                    F3=F1*(0.5*exp(-abs((xx-(sst(1)+sst(3))/2).^2*F2/2+...
                        (yy-(sst(2)+sst(4))/2).^2*F2/0.5))+...
                        exp(-abs((xx-sst(1)).^2*F2/2+...
                        (yy-sst(2)).^2*F2/2))+...
                        exp(-abs((xx-sst(3)).^2*F2/2+...
                        (yy-sst(4)).^2*F2/2)));
                    if find(F3)>0
                        F=F3/sum(F3);
                    else
                        F=F3;
                    end
                    no=F<0.0001;
                    F(no)=0;
                    
                    %Inversion matrix for Qc
                    Ac(indexray,:)=F;
    
                end
                
                % In the case a velocity model is available
                if pa==3
                    
                    % Function for ray-bending
                    rma=tracing(ray,gridD,pvel);
                    
                    % Set this to create reay-files
                    if createrays ==1
                        lrma=length(rma(:,1));
                        ma1(1:lrma,:,indexray) = rma;
                    end
                
                    %Creates figure with rays
                    hold on
                    subplot(2,2,1)
                    plot(rma(:,2),rma(:,3),'k');
                    xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
                    ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
                    hold off
                    hold on
                    subplot(2,2,2)
                    plot(rma(:,4),rma(:,3),'k');
                    xlabel('Depth UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
                    ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
                    hold off
                    hold on
                    subplot(2,2,3)
                    plot(rma(:,2),-rma(:,4),'k');
                    xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
                    ylabel('Depth UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
                    hold off
                    
    %===================================================================
    % We calculate segments in a grid, defined from the velocity model
    % The rays are in the cartesian coordinates and in the
    % reference system of the velocity
    % model, first point at the source, and with positive depth.
    % The velocity model is a x,y,z grid,
    % with the same steps in the three directions. The outputs are
    % necessary to invert for direct-wave attenuation.
    %===================================================================

                    
                    [lunpar, blocch, lunto, s, modv] =...
                        segments_single(modv,originz,um,rma);
                    lunparz(1:length(lunpar),indexray)=lunpar;
                    blocchi(1:length(blocch),indexray)=blocch;
                    luntot(indexray,1)=lunto;
                    sb(1:length(s),indexray)=s;
                    
                end
            end
        end
    end
    save(inputinv,'Apd', 'Ac', 'D');
        
    if pa<3
        
        %Create figure with 2D rays
        for nn=1:lls
            hold on
            plot([evestaz(nn,1) evestaz(nn,3)],...
                [evestaz(nn,2) evestaz(nn,4)],'k-')
        end
        hold on
        scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
        hold off
        
        grid on
        ax = gca;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        ax.GridAlpha = 1;
        ax.LineWidth = 1;
        
    elseif pa==3
        
        hold on
        
        %Create figure with 3D rays
        subplot(2,2,1)        
        scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        hold on
        subplot(2,2,2)        
        scatter(-even(:,3),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(-staz(:,3),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        hold on
        subplot(2,2,3)        
        scatter(even(:,1),even(:,3),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,3),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        
        % Only solve for the blocks crossed by at least 1 ray
        fover = modv(:,5)>=1; 
        inblocchi = modv(fover,:);
        inblocchi(:,5)=find(fover);

        %INVERSION MATRIX for direct waves: A
        A = zeros(length(Ac(:,1)),length(inblocchi(:,1)));
        for j = 1:length(sb(1,:))
            for i = 2:length(sb(:,1))
                bl = blocchi(i,j);
                slo = sb(i,j);
                l = lunparz(i,j);
                fbl = find(inblocchi(:,5)==bl);
                if isempty(fbl)==0
                    % check the inversion matrix element IN SECONDS
                    A(j,fbl)=-l*slo;
                end
            end
        end
    
        save(inputinv,'A', 'lunparz', 'blocchi', 'luntot', 'inblocchi',...
            'sb', '-append');
    
        if createrays == 1
            save(inputinv,'ma1','-append')
        end
    end
    
    % Saves the figure of RAYS
    FName = 'Rays';
    saveas(rays,fullfile(FPath, FLabel, FName), fformat);
    
    if pa>1
        
        % Kernel sensitivity figure
        Qcsen=figure('Name','Qc sensitivity, first source-station pair',...
            'NumberTitle','off','visible',visib);
        Qcss=Ac(1,:);
        Qcs=zeros(nxc,nyc);

        index=0;
        for i=1:length(x)
            for j=1:length(y)
                index=index+1;        
                Qcs(i,j)=Qcss(index);
            end
        end
        
        [X,Y]=meshgrid(x,y);
        contourf(X,Y,Qcs')
        
        colorbar
        grid on
        ax = gca;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        ax.GridAlpha = 1;
        ax.LineWidth = 1;
        
        FName = 'Qc_sensitivity';
        saveas(Qcsen,fullfile(FPath, FLabel, FName), fformat);
    end
end
% rtoc=toc;
%% Seismic attributes for peak delay, Qc and Qp,s imaging
clear
% tic
disp('Data Section')
load inputs.mat

%==========================================================================
% A loop to measure the Qc, peak-delay, and energy ratios for
% each seismic trace located in
% the folder "traces" with the coda-normalization method.
% The seismograms may be the output of both single and three-component
% seismic stations. In case of more than one component the rays must be
% ordered in the folder starting with the WE component, then SN and finally
% the Z component for each ray.
% 
% The program accepts SAC files and uses the information in the header
% The header must include peaking of the direct phase of interest, marker
% "a" for P and "t0" for S. If evst=2 you can get event and station info
% directly from the header. This will create the rays and related files
% automatically

if pa==3
    signal = zeros(lls,1); % The direct wave energies
    coda = zeros(lls,1); % The coda wave energies
    rapsp = zeros(lls,1); % The spectral ratios
    rapspcn = zeros(lls,1); % The coda versus noise ratios
end

tempi= zeros(lls,4); %Variables containing time origin (column 1), 
Qm = zeros(lls,1); %Qc with linearised or non-linear theory
RZZ = zeros(lls,1); %Correlation coefficient with respect to linear
ttheory=zeros(lls,1);
peakd=zeros(lls,1);
constQmean=zeros(2,2);
tlapse=(tCm+nW/2:nW:tCm+tWm-nW/2)';

index=0;
indexlonger=0;
    
if evst==1
    load(inputinv)
    
elseif evst==2
    disp('& Inversion Section for evst = 2')
    %Predefine the inversion matrix for peak-delay and coda-Q imaging
    even=zeros(lls,3);
    staz=zeros(lls,3);
    Ac=zeros(lls,lxy);
    Apd=zeros(lls,lxy);
    D=zeros(lls,1);
    evestaz=zeros(lls,1);
    rays=figure('Name','Rays','NumberTitle','off','visible',visib);
        
    if pa==3
        lunparz = zeros(100,lls);
        blocchi = zeros(100,lls);
        sb = zeros(100,lls);
        luntot = zeros(lls,1);
        if createrays ==1
            ma1=zeros(100,5,lls);
        end
    end
end

for i=1:lls
    index = index+1;
    if isequal(mod(index,100),0)
        disp(['Waveform number is ', num2str(index)])
    end
    
    %Read seismogram and get peaking/event/station info
    [tempis,sisma,SAChdr] = fget_sac(listasac{i}); % Converts SAC files
    if SAChdr.times.a == -12345
        tempi(i,2)=SAChdr.times.t0;
    else
        tempi(i,2)=SAChdr.times.a; %P-wave peaking
        tempi(i,3)=SAChdr.times.t0; %S-wave peaking
    end
    srate=1/SAChdr.times.delta; %sampling frequency
    
    if i==seisi && seesp==1 %to see the seisi spectrogram
        spect=figure('Name','Spectrogram','NumberTitle','off','visible',visib);
        view(2)
        title('Spectrogram of the first seism','FontSize',12,...
            'FontWeight','bold','Color','k');
        spectrogram(sisma,50,30,50,srate,'yaxis')
        FName = 'Spectrogram';
        saveas(spect,fullfile(FPath, FLabel, FName), fformat);
    end
    
    % Filter
    Wn = ([cf-cf/3 cf+cf/3]/srate*2); %frequency band
    [z,p,k] = butter(4,Wn,'bandpass'); %butter filter
    % [z,p,k] = cheby1(4,4,Wn,'bandpass'); %chebyshev filter
    [sos,g] = zp2sos(z,p,k); % Convert to SOS form
    Hd = dfilt.df2tsos(sos,g); % Create a dfilt object   

    if evst==2
        even1=[SAChdr.event.evlo SAChdr.event.evla SAChdr.event.evdp];
        staz1=[SAChdr.station.stlo SAChdr.station.stla SAChdr.station.stel];
        even(i,1:3)=even1; %#ok<SAGROW>
        staz(i,1:3)=staz1; %#ok<SAGROW>
        ray(1,1) = even1(1);
        ray(2,1) = even1(2);
        ray(3,1) = -even1(3);
        ray(1,2) = staz1(1);
        ray(2,2) = staz1(2);
        ray(3,2) = -staz1(3);
        xx=XY(:,1);
        yy=XY(:,2);    
        miix=min(even1(1),staz1(1));
        maax=max(even1(1),staz1(1));
        miiy=min(even1(2),staz1(2));
        maay=max(even1(2),staz1(2));
        fxy=find(xx>miix & xx<maax & yy>miiy & yy<maay);
        Apd(i,fxy)=1;
        sst = [even1(1) even1(2) staz1(1) staz1(2)];
        evestaz(i,1:4)=sst;
        D(i,1)=sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2);
        if pa==1
            
            % Same pd and Qc inversion matrix
            Ac=Apd;
            
        elseif pa>1
            
            % Operations to calculate normalised kernels
            D1=sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2);
            deltaxy=0.2;
            F1=1/2/pi/deltaxy^2/D1^2;
            F2=1/deltaxy^2/D1^2;
            F3=F1*(0.5*exp(-abs((xx-(sst(1)+sst(3))/2).^2*F2/2+...
                (yy-(sst(2)+sst(4))/2).^2*F2/0.5))+...
                exp(-abs((xx-sst(1)).^2*F2/2+...
                (yy-sst(2)).^2*F2/2))+...
                exp(-abs((xx-sst(3)).^2*F2/2+...
                (yy-sst(4)).^2*F2/2)));
            if find(F3)>0
                F=F3/sum(F3);
            else
                F=F3;
            end
            no=F<0.0001;
            F(no)=0;
            
            %Inversion matrix for Qc
            Ac(i,:)=F;
            
        end
        
        if pa==3
            rma=tracing(ray,gridD,pvel);
            if createrays ==1
                lrma=length(rma(:,1));
                ma1(1:lrma,:,i) = rma;
            end
                
            lrma=length(rma(:,1));
            
            %Creates figure with rays
            hold on
            subplot(2,2,1)
            plot(rma(:,2),rma(:,3),'k');
            xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
            ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
            hold off
            hold on
            subplot(2,2,2)
            plot(rma(:,4),rma(:,3),'k');
            xlabel('Depth UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
            ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
            hold off
            hold on
            subplot(2,2,3)
            plot(rma(:,2),-rma(:,4),'k');
            xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
            ylabel('Depth UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
            hold off
            
            %Function creating the parameters for the inversion matrix
            [lunpar, blocch, lunto, s, modv] =...
                segments_single(modv,originz,um,ma1);
            lunparz(1:length(lunpar),i)=lunpar;
            blocchi(1:length(blocch),i)=blocch;
            luntot(i)=lunto;
            sb(1:length(s),i)=s;
            
        end
    end
    
    pktime = tempi(i,PorS); %Picked time from SAC file
    
    if pa<3
        if SAChdr.times.o == -12345
            if SAChdr.times.b == -12345
                if evst==1
                    ttheory(i,1)=D(i)/vth;
                elseif evst==2
                    ttheory(i,1)=D(i)*degorutm/vth;
                end
                tempi(i,1)=pktime-ttheory(i,1);
            else
                tempi(i,1)=SAChdr.times.b;
            end   
        elseif SAChdr.times.o ~= -12345
            tempi(i,1)=SAChdr.times.o;
        end
    elseif pa==3
        ttheory(i,1)=sum(lunparz(:,i).*sb(:,i));
        tempi(i,1)=pktime-ttheory(i,1);
    end
       
    t00=tempis(1); %starting time of the waveform
    
    %starting-sample direct window
    cursor1 = floor((tempi(i,PorS)-t00)*srate);
    
    %Look for peak delay
    pdcursor=floor(cursor1+maxtpde*srate); %end sample to pick the pick delay
    tpdsisma = sisma(cursor1:pdcursor);% tapering direct wave
    fpdsisma = filter(Hd,tpdsisma);% filtering to find peak
    hpdsp = hilbert(fpdsisma); %hilbert
    pdsp = smooth(abs(hpdsp),nf/cf*srate);%ms of the filtered waveform
    [mspm,tspm]=max(pdsp);
    peakd(i,1)=tspm/srate;

    %Compute Qc from late lapse times (tCm)
    cursorc3 = floor((tempi(i,1)-t00+tCm-1)*srate);%coda start sample
    cursorc4 = floor(cursorc3 + tWm*srate-1);%coda end sample
    lsis=length(sisma);
    if cursorc4 > lsis
        tsismacm = sisma(cursorc3:lsis);%in case coda is too long
        indexlonger=indexlonger+1;
    else
        tsismacm = sisma(cursorc3:cursorc4);%tapering
    end
    L=length(tsismacm);
    tu=tukeywin(L,.05);
    tsismacm=tu.*tsismacm;    
    fsismacm = filter(Hd,tsismacm); %filter coda
    hspcm = hilbert(fsismacm); %hilbert
    spcm = smooth(abs(hspcm),nf/cf*srate);%ms of the filtered waveform
    lspm = length(spcm);
    tm = (tCm+1/srate:1/srate:tCm+lspm/srate)';
    
    %Only evaluate central time series
    edgeno=floor(0.05*length(tm));
    tm1=tm(edgeno:end-edgeno);
    spcm1=spcm(edgeno:end-edgeno);
    
    if nonlinear==0

        %THIS IS THE DATA VECTOR WITH LINEARISED THEORY - STANDARD
    
        EWz=spcm1.*tm1.^sped; 
        lspmz = log(EWz)/2/pi/cf; % source-station data
        Rz=corrcoef([tm1,lspmz]); %sets uncertainty
        polyz = polyfit(tm1,lspmz,1); %Qc^-1
    
        if polyz(1)<0
            Qm(i,1)=-polyz(1);
            RZZ(i,1)=abs(Rz(1,2));
        end    
    
    elseif nonlinear==1

    %THIS IS THE SYSTEM OF EQUATIONS FOR THE NON-LINEAR SOLUTION - USING
    %THE LINEAR INVERSE QC AS STARTING MODEL
        d_obs=zeros(ntW,1);
        ntm=0;
        for k=1:ntW
            lntm=length(ntm)+nW*srate;
            if lntm>lspm
                break
            end
            ntm = (k-1)*nW*srate + 1:k*nW*srate;
        
            nt= tm(floor(ntm));
            d_obs(k,1)=mean(spcm(floor(ntm)));    
        end
        d_obs1=d_obs(1:end-1)/d_obs(end);
    
        E=zeros(L1,1);    
        for n=1:L1
            d_pre=tlapse.^(-sped).*exp(-2*pi*cf.*tlapse*m1a(n));
            d_pre1=d_pre(1:end-1)/d_pre(end);
            E(n,1)=sum(abs(d_obs1-d_pre1));
        end
        [Emin, indexE] = min(E);
        Qm(i,1)=m1a(indexE);
        RZZ(i,1)=1/Emin;
    end
    
    if pa==3
        
        %Direct-wave intensity
        int = fin1*srate;% number of sample for the chosen window
        intn = fin1*srate;% number of sample for the noise
        cursor2 = floor(cursor1 + int-1); %end-sample of the window for direct
        tsisma = sisma(cursor1:cursor2);% tapering direct wave
        fsisma = filter(Hd,tsisma);% filtering direct wave
        hsp = hilbert(fsisma); %hilbert
        sp = smooth(abs(hsp),nf/cf*srate);%ms of the filtered waveform
        spamp = trapz(sp); %direct energy
       
        %Coda measure for coda-normalization
        cursorc1 = floor((tempi(i,1)-t00+tCm-1)*srate);%coda start sample
        cursorc2 = floor(cursorc1 + tWm*srate-1);%coda end sample
        if cursorc2 > length(sisma)
            lsis=length(sisma);
            tsismac = sisma(cursorc1:lsis);%in case coda is too long
        else
            tsismac = sisma(cursorc1:cursorc2);%tapering
        end
        fsismac = filter(Hd,tsismac); %filter coda
        hspc = hilbert(fsismac); %hilbert
        spc = smooth(abs(hspc),nf/cf*srate);%ms of the filtered waveform
        spampc = trapz(spc); %coda energy
    
        %Noise
        if PorS == 2
            cursorn1 = floor((tempi(i,2)-t00-sn)*srate); %starting sample for noise
        elseif PorS == 3
            cursorn1 = floor((t00)*srate+1); %starting sample for noise
        end
        cursorn2 = floor(cursorn1 + intn-1); % end-sample of the noise-window
        tsisman = sisma(cursorn1:cursorn2);%tapering noise
        fsisman = filter(Hd,tsisman);%filtering noise
        hspn = hilbert(fsisman); %hilbert
        spn = smooth(abs(hspn),nf/cf*srate);%ms of the filtered waveform
        spampn = trapz(spn); %noise energy
    
        rmes = spamp/spampc;%spectral ratios signal-coda
        rmcn = spampc/spampn;%spectral ratios coda-noise
    
        %STORE
        signal(index,1) = spamp; %direct-wave energies
        coda(index,1) = spampc; %coda-wave energies
        rapsp(index,1) = rmes; %spectral ratios
        rapspcn(index,1)= rmcn; %coda-to-noise ratios
    
    end
end

noQm=sum(Qm==0)/length(Qm)*100;
if nonlinear==0
    noRZZ=sum(RZZ<=RZZ2)/length(RZZ)*100;
elseif nonlinear == 1
    noRZZ=sum(RZZ>=RZZ2)/length(RZZ)*100;
end
if noQm>20
    Qless=['Attention: ',num2str(noQm),...
        '% of your Q are <=0 and ',num2str(noRZZ),...
        '% of your uncertainty coefficients are lower than treshold'];
    
    warning(Qless)
    
else
    if noRZZ>20
        Rless=['Attention: ',num2str(noRZZ),...
            '% of your coerrelation coefficients are lower than treshold'];
        warning(Rless)
    end
end
    
if evst==2
    figure(rays)
    if degorutm==111
        load coastlines
        hold on
        geoshow(coastlat,coastlon);
        xlim([8 18]);ylim([36 46]);
    end
    if pa<3
        
        for nn=1:lls
            hold on
            plot([evestaz(nn,1) evestaz(nn,3)],...
                [evestaz(nn,2) evestaz(nn,4)],'k-')
        end
        hold on
        scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1)
        hold off
        grid on
        ax = gca;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        ax.GridAlpha = 1;
        ax.LineWidth = 1;
        
        save(inputinv,'Apd', 'Ac', 'D');
        save('inputs.mat', 'even', 'staz','-append');
    
    elseif pa==3
    
        save('inputs.mat', 'even', 'staz','-append');
    
        hold on
        
        subplot(2,2,1)        
        scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        hold on
        subplot(2,2,2)        
        scatter(-even(:,3),even(:,2),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(-staz(:,3),staz(:,2),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
        hold on
        subplot(2,2,3)        
        scatter(even(:,1),even(:,3),sz,'c','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold on
        scatter(staz(:,1),staz(:,3),sz,'^','MarkerEdgeColor',...
            [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1)
        hold off
    
    %INVERSION MATRIX for direct waves: A - in case of event-station in SAC
    % Only store the blocks crossed by at least 1 ray
        fover = modv(:,5)>=1; 
        inblocchi = modv(fover,:);
        inblocchi(:,5)=find(fover);

        A = zeros(length(Ac(:,1)),length(inblocchi(:,1)));
        for j = 1:length(sb(1,:))
            for i = 2:length(sb(:,1))
                bl = blocchi(i,j);
                slo = sb(i,j);
                l = lunparz(i,j);
                fbl = find(inblocchi(:,5)==bl);
                if isempty(fbl)==0
                    A(j,fbl)=-l*slo;% inversion matrix element in sec
                end
            end    
        end
    
        save(inputinv,'Apd', 'Ac', 'D','A', 'lunparz', 'blocchi', 'luntot',...
            'inblocchi', 'sb');
    
        if createrays == 1
            save(inputinv,'ma1','-append')
        end
    
    end
    
    FName = 'Rays';
    saveas(rays,fullfile(FPath, FLabel, FName), fformat);
    
    if pa>1
        
        Qcsen=figure('Name','Qc sensitivity, first source-station pair',...
            'NumberTitle','off','visible',visib);
        Qcss=Ac(30,:);
        Qcs=zeros(nxc,nyc);

        index=0;
        for i=1:length(x)
            for j=1:length(y)
                index=index+1;        
                Qcs(i,j)=Qcss(index);
            end
        end
        
        [X,Y]=meshgrid(x,y);
        contourf(X,Y,Qcs')
        
        colorbar
        grid on
        ax = gca;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        ax.GridAlpha = 1;
        ax.LineWidth = 1;
        hold on
        geoshow(coastlat,coastlon);
        xlim([8 18]);ylim([36 46]);
        
        FName = 'Qc_sensitivity';
        saveas(Qcsen,fullfile(FPath, FLabel, FName), fformat);
    end
    
end  

% Setting up the inversion vector in case of 2- and 3-components data
icomp = 0;
ll=lls;

%Matrix A
if evst==2 && compon>1
    lA=length(A(:,1));
    Ac=Ac(1:compon:lA-compon,:);
    Apd=Apd(1:compon:lA-compon,:);
    A=A(1:compon:lA-compon,:);
end

% 2 components (tipically the two horizontals)
if compon ==  2
    lsig = ll/2;
    signal1 = zeros(lsig,1); % The average direct wave energies
    coda1 = zeros(lsig,1); % The average coda wave energies
    rapsp1 = zeros(lsig,1); % The average spectral ratios
    rapspcn1 = zeros(lsig,1); % The average coda versus noise ratios
    Qm1 = zeros(lsig,1); % The average coda versus noise ratios
    RZZ1 = zeros(lsig,1); % The average coda versus noise ratios
    peakd1 = zeros(lsig,1); % The average coda versus noise ratios
    for i = 1:2:(ll-1)
        icomp = icomp+1;
        signal1(icomp,1) = (signal(i,1)+signal(i+1))/2;
        coda1(icomp,1) = (coda(i)+coda(i+1))/2;
        rapsp1(icomp,1) = (rapsp(i)+rapsp(i+1))/2;
        rapspcn1(icomp,1) = (rapspcn(i)+rapspcn(i+1))/2;
        Qm1=(Qm(i,1)+Qm(i+1))/2;
        RZZ1=(Rzz(i,1)+Rzz(i+1))/2;
        peakd1=(peakd(i,1)+lpeakd(i+1))/2;
    end    
    signal=signal1;
    coda=coda1;
    rapsp=rapsp1;
    rapspcn=rapspcn1;
    Qm=Qm1;
    RZZ=RZZ1;
    peakd=peakd1;
    
% 3 components (WE, SN, Z)
elseif compon == 3
    lsig = ll/3;
    signal1 = zeros(lsig,1); % The average direct wave energies
    coda1 = zeros(lsig,1); % The average coda wave energies
    rapsp1 = zeros(lsig,1); % The average spectral ratios
    rapspcn1 = zeros(lsig,1); % The average coda versus noise ratios
    Qm1 = zeros(lsig,1); % 
    RZZ1 = zeros(lsig,1); % 
    peakd1 = zeros(lsig,1); %
    
    for i = 1:3:(ll-2)
        icomp = icomp+1;
        signal1(icomp,1) = ((signal(i)+signal(i+1))/2 + signal(i+2))/2;
        coda1(icomp,1) = ((coda(i)+coda(i+1))/2 + coda(i+2))/2;
        rapsp1(icomp,1) = ((rapsp(i)+rapsp(i+1))/2 + rapsp(i+2))/2;
        rapspcn1(icomp,1) = ((rapspcn(i)+rapspcn(i+1))/2 + rapspcn(i+2))/2;
        Qm1(icomp,1) = ((Qm1(i)+Qm1(i+1))/2 + Qm1(i+2))/2;
        RZZ1(icomp,1) = ((RZZ1(i)+RZZ1(i+1))/2 + RZZ1(i+2))/2;
        peakd1(icomp,1) = ((peakd(i)+peakd(i+1))/2 + peakd(i+2))/2;
    end
    signal=signal1;
    coda=coda1;
    rapsp=rapsp1;
    rapspcn=rapspcn1;
    Qm=Qm1;
    RZZ=RZZ1;
    peakd=peakd1;
end

if pa<3
    lu=D;
    if evst==1
        time0=D/vth;
    elseif evst==2
        time0=D*degorutm/vth;
    end
    
elseif pa==3
    %Calculate average geometrical spreading factor and Q
    %sets the constant in the CN method
    lu=luntot;
    
    constCNQc=(tCm+tWm/2)^sped.*exp(Qm*2*pi*cf*(tCm+tWm/2));
    mcn = find(rapspcn<tresholdnoise);% set the weigth

    % weighting
    time0=tempi(1:compon:end-compon+1,PorS)-tempi(1:compon:end-compon+1,1);
    W1=rapspcn<tresholdnoise;
    W2=rapsp>1000;
    W3=rapsp<tresholdnoise;
    W4=W3+W2+W1;
    W=diag(1./(W4+1));% weights
    d0=log(rapsp./constCNQc)/2/pi/cf; %data of the inverse problem

    dW=W*d0;
    G1=-log(lu)/pi/cf; %matrix creation
    G1(:,2)=-time0;

    G=W*G1;
    %The 3 parameters are the constant, the geometrical spreading, and the
    %average Q and they are contained in constQmean.
    constQmean(1,1:2)=lsqlin(G,dW(:,1));% damped least square inversion
    cova = (G'*G)^(-1)*G'*cov(dW)*G*(G'*G)^(-1); %covariance matrix
    er = sqrt(diag(cova)); %error from the covariance matrix
    constQmean(2,:)=er;

end

%Peak-delay, for reference see e.g. Calvet et al. 2013, Tectonophysics
l10l=log10(lu);
l10pd=log10(peakd);
fitrobust = fit(l10l,l10pd,'poly1','Robust','on');
pab=[fitrobust.p1 fitrobust.p2];
l10=linspace(min(l10l),max(l10l));
lpd=pab(1)*l10+pab(2);
l10pdt=polyval(pab,l10l);
lpdelta=l10pd-l10pdt;

% To remove outliers
I = abs(lpdelta) > 2*std(lpdelta) | peakd >= maxtpde | peakd < mintpde;
outlierspd = excludedata(l10l,lpdelta,'indices',I);

%Same for Qc
%Remove anomalous Qm remembering they have a log-normal distribution
if nonlinear==0
    retainQm=find(Qm>0 & RZZ>RZZ2);
    mQm=mean(Qm(retainQm));
    outliersc = Qm>0 & Qm > mQm+2*std(Qm(retainQm));
    retainQm=find(Qm>0 & RZZ>RZZ2 & outliersc==0);
    discardQm=find(Qm<=0 | RZZ<RZZ2 | outliersc==1);
    discardQm2= find(Qm>0 & (RZZ<RZZ2 | outliersc==1));
    
elseif nonlinear==1
    retainQm=find(Qm>0 & RZZ<RZZ2);
    mQm=mean(Qm(retainQm));
    outliersc = Qm>0 & Qm > mQm+2*std(Qm(retainQm));
    retainQm=find(Qm>0 & RZZ<RZZ2 & outliersc==0);
    discardQm=find(Qm<=0 | RZZ>RZZ2 | outliersc==1);
    discardQm2= find(Qm>0 & (RZZ>RZZ2 | outliersc==1));
end

Qmdis=Qm;
Qm(discardQm)=mQm;

%save quantities for the data vector
save(inputdata,'Qm','RZZ','lpdelta','outlierspd','time0','mQm')

if pa==3
    save(inputdata,'signal','coda','rapsp',...
        'rapspcn','constQmean','d0','W','-append');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS 

%plot to check that Qc is constant with travel time and peak delays
%increase with travel time

Qcpd=figure('Name','Qc and peak-delays',...
    'NumberTitle','off','visible',visib);
subplot(2,1,1)
plot(time0(retainQm),Qm(retainQm),'o',...
    'MarkerSize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
hold on
plot(time0(discardQm2),Qmdis(discardQm2),'o',...
    'MarkerSize',6,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0])
hold off
title('Dependence of inverse Qc on travel time');
xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Inverse Qc','FontSize',12,...
    'FontWeight','bold','Color','k')

subplot(2,1,2)
plot(fitrobust,'k--',l10l,l10pd,'ko',outlierspd,'r*')
title('Dependence of peak delays on travel times');
xlabel('Log. Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Log. Peak delay (s)','FontSize',12,...
    'FontWeight','bold','Color','k')
FName = 'Qc_Peak_Delay';
saveas(Qcpd, fullfile(FPath, FLabel, FName), fformat);

if pa==3
% The cyan dots represent the left-hand side of the
% CN method linear equation. The red line is the fit of the three
% parameters. The black lines are the uncertainties.
% In the lower panel the coda-to noise ratios, to see the effe ct of noise
% on the data.

    % Residual geometrical spreading in the coda
    gc=polyfit(time0,log(rapspcn),1);
    totgs=-(constQmean(1,1));
    totgs1=-(constQmean(1,1)-er(1));
    totgs2=-(constQmean(1,1)+er(1));
    % Right hand side of the CN linear equation - corrected by residual
    % geometrical spreading
    x = totgs*log(luntot)/pi/cf...
        -time0*constQmean(1,2);
    xer1 = totgs1*log(luntot)/pi/cf...
        -time0*(constQmean(1,2)-er(2));
    xer2 = totgs2*log(luntot)/pi/cf...
        -time0*(constQmean(1,2)+er(2));

%Plot of the left and right hand sides of the CN equation. The average
%inverse Q and geometrical spreading coefficient are schown above each
%panel
    CN=figure('Name','Coda-normalized intensities vs travel times',...
        'NumberTitle','off','visible',visib);
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

    title(['Average inverse quality factor = ', num2str(constQmean(1,2)),...
        ' +/- ', num2str(er(2)), ' and geometrical spreading = ',...
        num2str(constQmean(1,1)),' +/- ', num2str(constQmean(2,1))],...
        'FontSize',12,'FontWeight','bold','Color','k')
    h_legend=legend('Direct-to-coda energy ratios','Inverse average Q',...
        'Inverse Q spreading-related uncertainties');
    set(gca,'XTick',0:round(max(time0))) ;
    set(gca,'YTick',min(log(rapsp/constCNQc)+...
        constQmean(1,1)*log(luntot))/pi/cf:...
        (max(log(rapsp/constCNQc)+...
        constQmean(1,1)*log(luntot))/pi/cf-...
        min(log(rapsp/constCNQc)+constQmean(1,1)*log(luntot))/pi/cf)/5:...
        max(log(rapsp/constCNQc)+constQmean(1,1)*log(luntot))/pi/cf) ;
    set(gca,'FontSize',10,'FontWeight','bold');
    set(h_legend,'FontSize',10,'FontWeight','bold');

% %Plot of coda-to-noise energy ratios. The residual geometrical
% spreading should be around zero, or coherent phases could be
% included in the coda-energy; it also provides the average coda-to-noise
% ratio
    subplot(2,1,2)
    plot(time0,(log(rapspcn)),'o','MarkerSize',6,'MarkerEdgeColor',...
        [.2 .8 0],'MarkerFaceColor',[.2 .8 0])
    mc = mean(rapspcn);
    mno=length(mcn);
    perno= mno/length(rapspcn)*100;
    hold on
    plot(time0(mcn),log(rapspcn(mcn)), 'o','MarkerSize',6,...
    'MarkerEdgeColor',[1 .5 0],'MarkerFaceColor',[1 .5 0])
    hold off
    xlabel('Travel time (s)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('Logarithm of the coda-to-noise energy ratio','FontSize',12,...
        'FontWeight','bold','Color','k')
    title(['Residual coda geometrical spreading = ',num2str(gc(1)), ...
        ' while the percentage of measures < ',num2str(tresholdnoise), ' is ',...
        num2str(perno), '%'],'FontSize',12,'FontWeight','bold','Color','k');
    h_legend1=legend('Coda-to-noise energy ratios','Ratios below treshold');
    set(gca,'XTick',0:round(max(time0))) ;
    set(gca,'YTick',min(log(rapspcn)):...
        (max(log(rapspcn))-min(log(rapspcn)))/5:...
        max(log(rapspcn)));
    set(gca,'FontSize',10,'FontWeight','bold');
    set(h_legend1,'FontSize',10,'FontWeight','bold');

    text(max(time0), min(log(rapspcn)),...
    [' Lapse time = ',num2str(tCm),' s, coda window = ',...
    num2str(tWm), ' s'],'VerticalAlignment','bottom',...
    'HorizontalAlignment','right','FontSize',12)
    FName = 'CN_plot';
    saveas(CN,fullfile(FPath, FLabel, FName), fformat);
end
% END PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rtoc=toc;
%%  2D peak-delay and Qc and 3D CN TOMOGRAPHIC INVERSIONS
clear
load inputs.mat
load(inputinv)
load(inputdata)
pd=XY(:,1);
pd(:,2)=XY(:,2);
pd(:,3)=-1000;
Qc=pd;

%Peak delay mapping
lpdelta_o=lpdelta(outlierspd==0);
Apd=Apd(outlierspd==0,:);
Apd(Apd~=0)=1;
lApd=size(Apd);
Apd1=Apd;
pdelay=zeros(lApd(2),1);
for i=1:lApd(1)
    Apd1(i,1:end)=Apd1(i,1:end)*lpdelta_o(i);
end
for j=1:lApd(2)
    a1=Apd1(:,j);
    sipd=find(a1);
    pdelay(j,1)=mean(a1(sipd));
end
pdelay(isnan(pdelay))=0;
pd(:,4)=pdelay;

% save peak-delay
save -ascii peakdelay.txt pd


%Qc mapping
%Remove anomalous Qm
retainQm=find(Qm~=mQm);
stat=Qm(retainQm);

if pa==1

    p = [0.1 0.9];
    y = quantile(stat,p);
    zk = y;
    mm = [mQm median(stat)];
    sk = [skewness(stat) kurtosis(stat)];

    lAc=size(Ac);
    Ac1=Ac;

    Qcf=zeros(lAc(2),1);
    for i=1:lAc(1)
        Ac1(i,1:end)=Ac1(i,1:end)*Qm(i);
    end
    for j=1:lAc(2)
        a1=Ac1(:,j);
        sic=find(a1);
        Qcf(j,1)=mean(a1(sic));
    end
    Qcf(isnan(Qcf))=0;
    Qcf(Qcf==0)=mQm;
    Qc(:,4)=Qcf;

elseif pa>1
    Ac1=Ac(retainQm,:);
    RZZ1=RZZ(retainQm,1);
    W1=RZZ1<0.3;
    W2=RZZ1<0.5;
    W3=RZZ1<0.7;
    W4=W3+W2+W1;
    Wc=diag(1./(W4+1));% weights
    
    dcW=Wc*stat;
    Gc=Wc*Ac1;
    [lA,llA]=size(Gc);
    [Uc,Sc,Vc]=svd(Gc);

    index1=0;
    
    % smooting parameter is user defined
    LcQc=figure('Name','L-curve Qc','NumberTitle','off');
    l_curve(Uc,diag(Sc),stat,'Tikh')
    tik0_regC=input('Your personal smoothing parameter for coda ');
    FName = 'Lc_Qc';
    saveas(LcQc,fullfile(FPath, FLabel, FName), fformat);
    close(LcQc)
    
    % picard plot
    PpQc=figure('Name','Picard-plot Qc',...
        'NumberTitle','off','visible','off');
    picard(Uc,diag(Sc),stat);
    FName = 'Picard_Qc';
    saveas(PpQc,fullfile(FPath, FLabel, FName), fformat);
    
    mtik0C=tikhonov(Uc,diag(Sc),Vc,stat,tik0_regC);
    Qc(:,4)=mtik0C;

    %Testing - Creating 2D checkerboard matrix
    nxc1=nxc/sizea;
    nyc1=nyc/sizea;
    I = imresize(xor(mod(1:nyc1, 2).', mod(1:nxc1, 2)),...
        [sizea*nyc1 sizea*nxc1], 'nearest');
    Qc(:,5)=I(1:end);
    Qc(Qc(:,5)==1,5)=latt;
    Qc(Qc(:,5)==0,5)=hatt;
    
    Qc5= Qc(:,5);
    re = Gc*Qc5;
    mcheckc=tikhonov(Uc,diag(Sc),Vc,re,tik0_regC);
    Qc(:,6)=mcheckc;
    Qc(:,7)=Gc(1,:);
    Qc(:,8)=Gc(2,:);

end

% save Qc
save -ascii Qc.txt Qc

if pa==3
% Direct wave attenuation

    %Data creation, removing the pre-calculated parameters
    d1 = d0  + constQmean(1,1)*log(luntot)/pi/cf...
    + time0*constQmean(1,2);


    % if the ray crosses a block not solved by the inversion, the
    % block is characterized by the average quality factor, and the data
    % vector is updated.
    
    A1=A;
    A=W*A1;

    % NEW DATA VECTOR

    dW1=W*d1;

    % tikhonov inversion - by using the programs in HANSEN et al. 1994
    [U,S,V]=svd(A);

    %sets the smoothing parameter - always user defined
    LcCN=figure('Name','L-curve coda-normalization','NumberTitle','off');
    l_curve(U,diag(S),dW1,'Tikh')
    tik0_reg=input('Your personal smoothing parameter ');
    FName = 'Lc_CN';
    saveas(LcCN,fullfile(FPath, FLabel, FName), fformat);
    close(LcCN)
    
    % picard plot
    PpCN=figure('Name','Picard coda-norm.','NumberTitle','off',...
        'visible','off');
    picard(U,diag(S),dW1);
    FName = 'Picard_CN';
    saveas(PpCN,fullfile(FPath, FLabel, FName), fformat);
    
    %results
    mtik0=tikhonov(U,diag(S),V,dW1,tik0_reg);
    mQ1=[inblocchi(:,1:3) mtik0];

    % Simple resolution tests: produces columns in output.
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
        modv(k:sizea2:passoy-sizea+k,6)=hatt;
        modv(sizea+k:sizea2:passoy-sizea+k,6)=latt;
    end
    for k=1:sizea-1
        modv(k*passoy+1:(k+1)*passoy,6)=modv(1:passoy,6);
    end
    for k=1:sizea
        modv(sizeap+k:sizea2:sizeap1-sizea+k,6)=latt;
        modv(sizeap+sizea+k:sizea2:sizeap1-sizea+k,6)=hatt;
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
        if modv(k,6)==hatt
            modv(sizeapx+k,6)=latt;
        elseif modv(k,6)==latt
            modv(sizeapx+k,6)=hatt;
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
    Q3D(:,4)=constQmean(1,2);
    % create a file of input
    in1=zeros(length(mtik0),1);
    input01 = zeros(length(mtik0),1);

    for i = 1:length(in1)
        in1(i)=find(modv(:,1)==inblocchi(i,1)  & modv(:,2)==inblocchi(i,2)...
            & modv(:,3)==inblocchi(i,3));
        Q3D(in1(i),4)=mtik0(i)+constQmean(1,2);
        input01(i)=modv(in1(i),6);
    end

    % result in the original reference system
    Q3D(:,1)=modv(:,1)+resol2;
    Q3D(:,2)=modv(:,2)+resol2;
    Q3D(:,3)=modv(:,3)+resol2;

    % input save
    Qin=[Q3D(:,1:3) modv(:,6)];

    A1=W*A;
    % output file and save
    re = A1*input01;
    Qou(:,1:3)=Q3D(:,1:3);
    Qou(:,4)=0;
    mcheck=tikhonov(U,diag(S),V,W*(re-mean(input01)),tik0_reg);
    Qou(in1,4)=mcheck;

    % WE ALSO OUTPUT THE DIAGONAL OF THE RESOLUTION MATRIX
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
    %First is with the entire set of actual results - this is improper but
    %used in some studies
    lmQ1=length(mQ1(:,4));
    syQ0=constQmean(1,1)*log(luntot)/pi/cf+time0.*constQmean(1,2)+...
        A*mQ1(:,4);

    G1=-log(luntot)/pi/cf; %matrix creation
    G1(:,2)=-time0;

    dsW=W*syQ0;
    G=W*G1;

    %The 2 parameters are the geometrical spreading and the
    %average Q and they are contained in synthQmean.
    synthQmean(1,1:2)=lsqlin(G,dsW(:,1));% damped least square inversion
    cova = (G'*G)^(-1)*G'*cov(dsW)*G*(G'*G)^(-1); %covariance matrix
    er = sqrt(diag(cova)); %error from the covariance matrix
    synthQmean(2,:)=er;

    %data creation
    syQ1 = W*(syQ0 + synthQmean(1,1)*log(luntot)/pi/cf...
        + time0*synthQmean(1,2));

    [syU,syS,syV]=svd(A);
    
    
    % picard plot
    Ppsy=figure('Name','Picard synth','NumberTitle','off','visible','off');                
    picard(U,diag(S),syQ1);

    
    %results
    sytik0=tikhonov(syU,diag(syS),syV,syQ1,tik0_reg);
    Qou(:,6)=synthQmean(1,2);
    Qou(in1,6)=synthQmean(1,2)+sytik0;
    Q3D(:,5)=Qin(:,4);
    Q3D(:,6:8)=Qou(:,4:6);

    %Second is user-defined - better
    %As example, create a 2-layer medium with average Q latt and hatt
    Qplus=find(inblocchi(:,3)>dtreshold);
    Qminus=find(inblocchi(:,3)<dtreshold);
    A2=A;
    syQ2=mQ1;
    syQ2(Qplus,4)=latt;
    syQ2(Qminus,4)=hatt;
    A3=zeros(lls,1);
    for k=1:lls
        A3(k,1)=A2(k,:)*syQ2(:,4);
    end
        

    tsp=sped; %theoretical geom. spr.
    d0=tsp*log(luntot)/pi/cf-A3;
    dW=W*d0;
    G1=-log(luntot)/pi/cf;
    G1(:,2)=-time0;

    G=W*G1;
    syconstQmean(1,1:2)=lsqlin(G,dW(:,1));% damped least square inversion
    cova = (G'*G)^(-1)*G'*cov(dW)*G*(G'*G)^(-1); %covariance matrix
    er = sqrt(diag(cova)); %error from the covariance matrix
    syconstQmean(2,:)=er;

    syQ3=tsp*log(luntot)/pi/cf-A3+time0*constQmean(1,2);

    dsW2=W*syQ3;
    
    sytik02=tikhonov(syU,diag(syS),syV,dsW2,tik0_reg);

    inup=Qou(:,3)>dtreshold;
    indown=Qou(:,3)<dtreshold;
    Qou(inup,7)=hatt;
    Qou(indown,8)=latt;
    Qou(in1,7)=syQ2(:,4);
    Qou(in1,8)=syconstQmean(1,2)+sytik02;
    Q3D(:,9:10)=Qou(:,7:8);

    % save Q
    save -ascii Q3D.txt Q3D

end

%% Creating maps
load inputs.mat
load(inputinv)
load(inputdata)
load coastlines

Qc=load('Qc.txt');
pd=load('peakdelay.txt');

pdel=zeros(length(x),length(y));
QQc=zeros(length(x),length(y));
QQchi=zeros(length(x),length(y));
QQcho=zeros(length(x),length(y));
[X,Y]=meshgrid(x,y);
index=0;
for i=1:length(x)
    for j=1:length(y)
        index=index+1;        
        pdel(i,j)=pd(index,4);
        QQc(i,j)=Qc(index,4);
        if pa>1
            QQchi(i,j)=Qc(index,5);
            QQcho(i,j)=Qc(index,6);
        end
    end
end

pdmap=figure('Name','Peak-delay map','NumberTitle','off','visible',visib);
contourf(X,Y,pdel');
axis equal
view(2)
hcb=colormap(autumn);
colorbar
xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
title('Logarithmic peak-delay variations',...
    'FontSize',12,'FontWeight','bold','Color','k');
hold on
scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
hold on
scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
if degorutm==111
    hold on
    geoshow(coastlat,coastlon);
    xlim([8 18]);ylim([36 46]);
end
hold off
FName = 'Peak_delay_map';
saveas(pdmap,fullfile(FPath, FLabel, FName), fformat);

Qcmap=figure('Name','Qc map','NumberTitle','off','visible',visib);
contourf(X,Y,QQc');
axis equal
view(2)
colormap(copper)
colorbar
xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
title('Inverse Qc','FontSize',12,'FontWeight','bold','Color','k');
hold on
scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
hold on
scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
if degorutm==111
    hold on
    geoshow(coastlat,coastlon);
    xlim([8 18]);ylim([36 46]);
end
hold off
FName = 'Qc_map';
saveas(Qcmap,fullfile(FPath, FLabel, FName), fformat);

%Parameter analysis
if pa ==1
    pdd = pd(:,4)~=0 & Qc(:,4)~=mQm;

    pdef=pd(pdd,:);
    Qcef=Qc(pdd,:);
    pd(:,5)=0;
% pdef=pd;
% Qcef=Qc;

elseif pa>=2
    
    Qccheck=figure('Name','Qc checkerboard test','NumberTitle',...
        'off','visible',visib);
    subplot(1,2,1)
    contourf(X,Y,QQchi');
    axis equal
    view(2)
    colormap(gray)
    colorbar
    xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    title('Qc checherboard input',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    if degorutm==111
        hold on
        geoshow(coastlat,coastlon);
        xlim([8 18]);ylim([36 46]);
    end
    hold off
    subplot(1,2,2)
    contourf(X,Y,QQcho');
    axis equal
    view(2)
    colormap(gray)
    colorbar
    xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    title('Qc checherboard output',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1.5)
    if degorutm==111
        hold on
        geoshow(coastlat,coastlon);
        xlim([8 18]);ylim([36 46]);
    end
    hold off
    
    FName = 'Qc_checkerboard';
    saveas(Qccheck,fullfile(FPath, FLabel, FName), fformat);
    
    %Parameter analysis
    pdd =  pd(:,4)~=0;

    pdef=pd(pdd,:);
    Qcef=Qc(pdd,:);
    
end

pdef(:,4)=pdef(:,4)-mean(pdef(:,4));
Qcef(:,4)=Qcef(:,4)-mean(Qcef(:,4));
miQcm=min(Qcef(:,4));
maQcm=max(Qcef(:,4));
mipdm=min(pdef(:,4));
mapdm=max(pdef(:,4));
trepd=0.05*mapdm;
treQc=0.05*maQcm;

Qps=Qcef(:,4);
pdps=pdef(:,4);

param_plot=figure('Name','Parameter space separation','NumberTitle','off',...
    'visible',visib);
hax=axes;

par=pdef(:,1:2);
par(:,3)=pdef(:,3);

c=Qps<-treQc & pdps<-trepd;
par(c,4)=1;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[0 0.8 0])
hold on
line([0 0],[mipdm-trepd mapdm+trepd],'Color',[0 0 0],...
    'LineWidth',3)
hold on
line([miQcm-treQc maQcm+treQc],[0 0],'Color',[0 0 0],...
    'LineWidth',3)
hold on
c=Qps<-treQc & pdps>trepd;
par(c,4)=2;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[0 0.6 1])
hold on
c=Qps>treQc & pdps<-trepd;
par(c,4)=3;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[1 0.6 0])
hold on
c=Qps>treQc & pdps>trepd;
par(c,4)=4;
scatter(Qps(c),pdps(c),65,'filled','MarkerFaceColor',[1 0 0])
hold on
c=(Qps>-treQc & Qps<treQc) | (pdps>-trepd & pdps<trepd);
par(c,4)=0;
scatter(Qps(c),pdps(c),85,'filled','MarkerFaceColor',[0.7 0.7 0.7],...
    'MarkerEdgeColor',[1 1 1],'LineWidth',2)
hold off 
xlim([miQcm-treQc maQcm+treQc])
ylim([mipdm-trepd mapdm+trepd])
xlabel('Qc','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Log. peak delay','FontSize',12,'FontWeight','bold','Color','k')
title('Parameter space plot',...
    'FontSize',12,'FontWeight','bold','Color','k');
FName = 'Parameter_space_variations';
saveas(param_plot,fullfile(FPath, FLabel, FName), fformat);

para=pd(:,1:3);
for k=1:length(par(:,1))
    px=par(k,1);
    py=par(k,2);
    pp=par(k,4);
    pf= pd(:,1)==px & pd(:,2)==py;
    para(pf,4)=pp;
end

index=0;
param=zeros(size(QQc));
for i=1:length(x)
    for j=1:length(y)
        index=index+1;        
        param(i,j)=para(index,4)-5;
    end
end

mparam=figure('Name','Parameter separation map',...
    'NumberTitle','off','visible',visib);
surf(X,Y,param');
view(2)
un_X = unique(param);

fu=find(un_X==-1);
if isempty(fu)
    cmap = [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0];
    HTick={'Average','Ls La','Hs La','Ls Ha'};
else
    cmap = [0.7 0.7 0.7;  0 0.8 0; 0 0.6 1; 1 0.6 0; 1 0 0];
    HTick={'Average','Ls La','Hs La','Ls Ha','Hs Ha'};
end
colormap(cmap)

colorbar('Ticks',un_X,'TickLabels',HTick);
    
hold on
scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1.5)
hold on
scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
    [1 1 1], 'MarkerFaceColor',[.5 .5 .5], 'LineWidth',1.5)
if degorutm==111
    hold on
    geoshow(coastlat,coastlon);
    xlim([8 18]);ylim([36 46]);
end
hold off
xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')

title('Parameter separation','FontSize',12,'FontWeight','bold','Color','k');
FName = 'Parameter_map';
saveas(mparam,fullfile(FPath, FLabel, FName), fformat);

%Image velocity model
if pa==3
    xD=unique(modv(:,1));
    yD=unique(modv(:,2));
    zD=sort(unique(modv(:,3)),'descend');
    
    ix=length(xD);%numer of x layers,given the step of the grid
    iy=length(yD);%numer of y layers,given the step of the grid
    iz=length(zD);%numer of depths layers,given the step of the grid
    
    V=zeros(iy,ix,iz);
    Q=zeros(iy,ix,iz);
    Qi=zeros(iy,ix,iz);
    Qo=zeros(iy,ix,iz);
    index=0;
    for i=1:ix
        for j=1:iy
            for k=1:iz
                index=index+1;
                V(j,i,k)=modv(index,4);
                Q(j,i,k)=Q3D(index,4);
                Qi(j,i,k)=Q3D(index,5);
                Qo(j,i,k)=Q3D(index,6);
            end
        end
    end
    
    [X,Y,Z]=meshgrid(xD,yD,zD);
    
    V_Q=figure('Name','Velocity Model','NumberTitle','off','visible',visib);
    
    %Velocity model
    subplot(1,2,1)
    slice(X,Y,Z,V,WEi,SNi,zi,'spline');
    colormap(cool)
    colorbar
    axis equal
    view(2)
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold off
    xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('Depth (m)','FontSize',12,'FontWeight','bold','Color','k')
    title('Velocity model',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    %Attenuation model
    subplot(1,2,2)
    slice(X,Y,Z,Q,WEi,SNi,zi,'spline');
    colormap(autumn)
    colorbar
    axis equal
    view(2)
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold off
    xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('Depth (m)','FontSize',12,'FontWeight','bold','Color','k')
    title('Attenuation model',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    FName = 'V_Q_models';
    savefig(V_Q,fullfile(FPath, FLabel, FName));
    
    Qcheck =...
        figure('Name','Checkerboard 3D','NumberTitle','off','visible',visib);
    
    %3D checkerboard input
    subplot(1,2,1)
    slice(X,Y,Z,Qi,WEi,SNi,zi,'spline');
    colormap(gray)
    colorbar
    axis equal
    view(2)
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold off
    xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('Depth (m)','FontSize',12,'FontWeight','bold','Color','k')
    title('Q 3D checkerboard input',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    %3D checkerboard output
    subplot(1,2,2)
    slice(X,Y,Z,Qo,WEi,SNi,zi,'spline');
    colormap(gray)
    colorbar
    axis equal
    view(2)
    hold on
    scatter(even(:,1),even(:,2),sz,'c','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold on
    scatter(staz(:,1),staz(:,2),sz,'^','MarkerEdgeColor',...
        [1 1 1], 'MarkerFaceColor',[0.5 .5 .5], 'LineWidth',1.5)
    hold off
    xlabel('WE UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('SN UTM (WGS84)','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('Depth (m)','FontSize',12,'FontWeight','bold','Color','k')
    title('3D checkerboard output',...
        'FontSize',12,'FontWeight','bold','Color','k');
    
    FName = '3DQ_check';
    savefig(Qcheck,fullfile(FPath, FLabel, FName));
end

