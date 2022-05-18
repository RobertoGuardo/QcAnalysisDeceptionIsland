% PROGRAM AIM: Qcoda in several frequency band with weighting
%
% Author:  L. De Siena.
%
%           Jan 2017
%
% Adapted for Decepction Island Volcano by R. Guardo 
%
%           Jan-Aug 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  STEP 0 - Import and Setup

addpath('OutImage');
addpath('Rays_it0_20283');
addpath('Rays_it1_14972');
addpath('Rays_it2_13105');
addpath('Rays_it3_7895');
addpath('Rays_it4_7196');
addpath('Traces');
addpath('Utilities_Matlab/');
 clc
 clear

disp('START');
sT = datetime;
warning off all

% % % Name of the file containing the path of the trace-files -  in traces
% % % list=dir('./Traces/');  
% % % isfile=~[list.isdir]; %determine index of files vs folders
% % % filenames={list(isfile).name}; %create cell array of file names
% % % lista = filenames';
% % % for i = 1:length(lista)
% % %     lista{i,1}=cat(2,'./Traces/',lista{i});
% % % end
% % % ll = length(lista);


%ForS = 2;  % Choose between P-(2) and S-direct(3) waves

%%%%% DCP Setup  

% tC = 3;  % starting lapse-time
% tW = 6;   % length of the window

nW = 1.5;   % noise windows (ex "fin1")
cf = 3;     % central frequency = 6, change at 2n loop (Prudencio et al., 2015)
cW = 4;     % coda windows
mL = 10;    % total length of the signal
sn = 0;     % start time to consider noise - before P-arrival
keepCf = cf; 
srate = 100; % sampling rate (time(n)-time(n-1))^-1;

res = 2;    % select resolution 1 = 1km, 2 = 2km

% CODA Choose between: 
% 1: coda between 4 and 10 s
% 2: coda between 6 and 10 s
% 3: coda between 4 and 8 s
% any other value: 4s coda from the tP 
coda = 4; 

%UTM coordinates of the origin of the velocity model
originWE = 618284;	% DCP
originSN = 3015685;	% DCP
originz = 0;		% DCP

% Iteration, choose from 0 to 4
% 0 = Full dataset 20283 events
% 1 = 14972 rays
% 2 = 13105 rays
% 3 = 7895 rays
% 4 = 7196 rays
it = 2;

%% STEP 1 - Trace analysis
disp('STEP 1 - Trace analysis');

switch it
    case 1
        lista = textread('traces_it1_14972.txt','%q','delimiter','\n');
    case 2
        lista = textread('traces_it2_13105.txt','%q','delimiter','\n');
    case 3
        lista = textread('traces_it3_7895.txt','%q','delimiter','\n');
    case 4
        lista = textread('traces_it4_7196.txt','%q','delimiter','\n');
    otherwise
        lista = textread('traces_it0_20283.txt','%q','delimiter','\n');
end

ll = length(lista);

Qm3 = zeros(ll,6);
noise = zeros(ll,6); % The coda versus noise ratios
RZZ = zeros(ll,6); % The coda versus noise ratios
        
for i = 1:ll
    
    disp(['Waveform: ',num2str(i)]);  % number of the waveform
    
    si = lista{i,1};    % start reading the vertical seismograms
    [sis] = textread(si,'%f');  % use textread to get the two columns files
    
    tempis = sis(1:2:end-1);    % first column is time
    sisma = sis(2:2:end);       % second is measurement
    
    tempis = tempis(1:srate*mL);    % cut the trace according to the tW
    sisma = sisma(1:srate*mL);
    
    cf = keepCf;        %Necessary to reset the cf before the for loop
    nf=[4 4 4 8 8 8];
    for n = 1:6         % I'll keep n=1, 2, 4 and 8
        
        cf = cf + 3;                        % Increase frequency
        Wn = ([cf-cf/3 cf+cf/3]/50);        % Frequency band
        
        [z,p,k] = butter(6,Wn,'bandpass');  % Butter filter
        [sos,g] = zp2sos(z,p,k);            % Convert to SOS form
        Hd = dfilt.df2tsos(sos,g);          % Create a dfilt object   
    
        % NOISE %---------------------------------------------------------%
        intn = nW*srate;                    % Number of sample for the noise
		cursorN1 = 1;                       % Starting sample for noise
        cursorN2 = cursorN1 + intn-1; 		% End-sample of the noise-window
        tsisman = sisma(cursorN1:cursorN2);	% Tapering noise
        
        % ENVELOPE %------------------------------------------------------%
        fsisman = filter(Hd,tsisman);		% filtering noise
        hspn = hilbert(fsisman);			% hilbert transform
        spmn = abs(hspn.^2);				% hilbert
        spsn = smooth(spmn,nf(n)/cf*srate);   
        spampn = sqrt(spsn); 				% measure energy noise
        
        % ENTIRE - We now compute the entire envelope energy %-------------%
		fsismac = filter(Hd,sisma);         % filter coda
        hsp = hilbert(fsismac);             % hilbert
        spm = abs(hsp.^2);                  % ms of the filtered waveform
        sps = smooth(spm,nf(n)/cf*srate);
        
        % CODA %----------------------------------------------------------% 
        
        switch coda
            case 1 % Coda from 4s to the end of the trace
                sps2 = sps(4*srate:length(sisma));
                
            case 2 % Coda from 6s to the end of the trace
                sps2 = sps(6*srate:length(sisma)); 

            case 3 % Coda from 4s to 8s
                sps2 = sps(4*srate:4*srate*2);
 
            otherwise % Coda from tS to (tS + 4s)
                [maxsps,imaxsps] = max(sps);
                if imaxsps > length(sisma)-cW*srate
                    continue
                end
                cLim = imaxsps + cW*srate;          % coda limit
                
                sps2 = sps(imaxsps:cLim);
                               
        end
             
		% STORE %---------------------------------------------------------%
        lspm = length(sps2);
        tm = (imaxsps:lspm+imaxsps-1)'/srate;
        lspm3 = log(sps2.*tm.^1.5)/2/pi/cf;     % linearizzazione
                
        % Rajout Marie
        Rz=corrcoef([tm,lspm3]); 
        RZZ(i,n)=abs(Rz(1,2));
        polysph3 = polyfit(tm,lspm3,1);
        Qm3(i,n)= -polysph3(1);      % Coda Quality factor 
        noise(i,n) = mean(spm)/mean(spampn);  
    end   
end

save dataQc.mat Qm3 noise RZZ

disp('End STEP 1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 2 - Statistics
disp('STEP 2 - Statistics');

load dataQc.mat

mm=zeros(6,2);%Statistic of Qc
sk=zeros(6,2);%Statistic of Qc
z=zeros(6,2);%Statistic of Qc
Qm3(isnan(Qm3))=0;
RZZ(isnan(RZZ))=0;
noise(isnan(noise))=0;
soglia=0.3;% Minimum RZZ
index1=0;

for i = 1:6
    % Using the coefficients Marie
    Q1 = Qm3(:,i);
    R1 = RZZ(:,i);
    n1 = noise(:,i); %Noise is always higher than 10 before eruption, skipping
    no = R1<soglia;
    Q1(no) = 0;
    
    stat = Q1;
    %{   
    if i==1 || i==2 || i==4 || i==6
        index1=index1+1;
        figure(1)
        subplot(2,2,index1)
        hist(stat,30);
        figure(2)
        subplot(2,2,index1)
        probplot('normal',stat)
    end
    
        
    p = [0.1 0.9];
    y = quantile(stat,p);
    z(i,1:2) = y;
    mm(i,1:2) = [mean(stat) median(stat)];
    sk(i,1:2) = [skewness(stat) kurtosis(stat)]; 
    %}
            
end

disp('End STEP 2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 3 - INVERSION - Matrix creation
disp('STEP 3 - INVERSION - Matrix creation');

%Name of the file containing the rays -  in rays
switch it
    case 1
        list=dir('./Rays_it1_14972');  %get info of files/folders in current directory
    case 2
        list=dir('./Rays_it2_13105');
    case 3
        list=dir('./Rays_it3_7895');
    case 4
        list=dir('./Rays_it4_7196');
    otherwise
        list=dir('./Rays_it0_20283');
end

isfile=~[list.isdir];           %determine index of files vs folders
filenames={list(isfile).name};  %create cell array of file names
listar = filenames';

switch it
    case 1
        for i = 1:length(listar)
            listar{i,1}=cat(2,'./Rays_it1_14972/',listar{i});
        end
    case 2
        for i = 1:length(listar)
            listar{i,1}=cat(2,'./Rays_it2_13105/',listar{i});
        end
    case 3
        for i = 1:length(listar)
            listar{i,1}=cat(2,'./Rays_it3_7895/',listar{i});
        end
    case 4
        for i = 1:length(listar)
            listar{i,1}=cat(2,'./Rays_it4_7196/',listar{i});
        end
    otherwise
        for i = 1:length(listar)
            listar{i,1}=cat(2,'./Rays_it0_20283/',listar{i});
        end
end
        
nrays = length(listar);

if res == 1
    %set the origin
    X = -9000; % -9000 with stepg = 1000 
    Y = -9000;
    
    stepg = 1000; %step of the grid in meters
    
    %set the number of nodes from the origin
    nx = 21; % 21 with stepg = 1000
    ny = 21;
    disp(['The grid step is: ', num2str(stepg), ' m']);
else
    X = -11000; % -11000 with stepg = 2000 
    Y = -11000;
    
    stepg = 2000; %step of the grid in meters
    
    %set the number of nodes from the origin
    nx = 12; % 12 with stepg = 2000 
    ny = 12;
    disp(['The grid step is: ', num2str(stepg), ' m']);
end

origin = [X,Y];
index = 0;

%build the nodes - they are your parameter model locations
XY=zeros(nx*ny,2);
for i=1:nx
    for j=1:ny
        index=index+1;
        XY(index,1:2)=[X+stepg*i Y+stepg*j]; % are the grid lines
    end
end
lxy=length(XY(:,1));

%Build the inversion matrix
A=zeros(nrays,lxy);
bs=zeros(nx*ny,2);
bst=zeros(nx*ny,2);
lung=zeros(nrays,1);

%source-station coordinates
for i=1:nrays
    
    fileID = fopen(listar{i,1});
    C = textscan(fileID,'%f %f %f %f %f' );
    fclose(fileID);

    
    sisWE=C{1,2}; % X
    sisSN=C{1,3}; % Y
        
    sst = [sisWE(1) sisSN(1) sisWE(end) sisSN(end)]; % prendere prima e ultima riga della seconda e terza colonna
    D=sqrt((sst(3)-sst(1))^2+(sst(4)-sst(2))^2);
    F1=1/2/pi/0.04/D^2;
    F2=1/.04/D^2;
    F3=F1*(0.5*exp(-abs((XY(:,1)-(sst(1)+sst(3))/2).^2*F2/2+...
        (XY(:,2)-(sst(2)+sst(4))/2).^2*F2/0.5))+...
        exp(-abs((XY(:,1)-sst(1)).^2*F2/2+(XY(:,2)-sst(2)).^2*F2/2))+...
        exp(-abs((XY(:,1)-sst(3)).^2*F2/2+(XY(:,2)-sst(4)).^2*F2/2)));
    F=F3/sum(F3);
    A(i,:)=F;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save A2D.mat A

disp('End STEP 3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 4 - Tikhonov Inversion
disp('STEP 4 - Tikhonov Inversion');

load A2D.mat A

lA=length(A(1,:));
llA= length(A(:,1));
clear Qc
A1=A;
Qc(:,1)=originWE+XY(:,1);
Qc(:,2)=originSN+XY(:,2);
Qc(:,3)=-4000;

smallper = A < 0.0001;
A(smallper) = 0;

[U,S,V]=svd(A);


for k = 1:6
    
%     figure
     tik0_reg(:,k)=l_curve(U,diag(S),Qm3(:,k),'Tikh');

%   tik0_reg, 4th it., 6Hz, 1km = 1.1829 (0.118 - 11.8)
%   tik0_reg, 4th it., 15Hz, 1km = 1.1863 (0.118 - 11.8)
%   tik0_reg, 4th it., 18Hz, 1km = 1.1501 (0.115 - 11.5)
%   tik0_reg, 4th it., 6Hz, 2km = 2.1686 (0.216 - 21.6)
%   tik0_reg, 4th it., 15Hz, 2km = 1.5586 (0.155 - 15.5)
%   tik0_reg, 4th it., 18Hz, 2km = 1.5966 (0.159 - 15.9)
    
%     tik0_reg = 11.5;
    
    %picard plot
%     figure
%     picard(U,diag(S),Qm3(:,k));
    
    mtik0=tikhonov(U,diag(S),V,Qm3(:,k),tik0_reg(:,k));
    Qc(:,3+k)=mtik0;
    
end
    
save A Qm3 Qc;

disp('End STEP 4');

%% FIRST TEST: CHECKERBOARD
disp('FIRST TEST: CHECKERBOARD');

%INPUT: Checkerboard structure - node spacing of the anomalies
% doubled with respect to the grid.

%Size anomaly for testing: twice (2)/four (4) times the original.
sizea=2;

%Depth
passox = find(Qc(:,1)~=Qc(1,1),1,'first')-1;

sizea2=2*sizea;
sizeap=sizea*passox;
sizeap1=(sizea+1)*passox;

for k=1:sizea
    Qc(k:sizea2:passox-sizea+k,10)=.02;
    Qc(sizea+k:sizea2:passox-sizea+k,10)=0.001;
end
for k=1:sizea-1
    Qc(k*passox+1:(k+1)*passox,10)=Qc(1:passox,10);
end
for k=1:sizea
    Qc(sizeap+k:sizea2:sizeap1-sizea+k,10)=.001;
    Qc(sizeap+sizea+k:sizea2:sizeap1-sizea+k,10)=0.02;
end

py4  = 2*sizeap;
for k=1:sizea-1
    Qc((sizea+k)*passox+1:(sizea+k+1)*passox,10)=Qc(sizeap+1:sizeap1,10);
end
z = (passox-mod(passox,py4))/py4;
for i = 1:(z-1)
    Qc(i*py4+1:(i+1)*py4,10)=Qc(1:py4,10);
end
if ~isequal(mod(passox,py4),0)
    Qc(z*py4+1:z*py4+mod(passox,py4),10)= Qc(1:mod(passox,py4),10);
end

%Along y
sizeapx=sizea*passox;

for k=1:sizea-1
    Qc(k*passox+1:(k+1)*passox,10)=Qc(1:passox,10);
end

for k = 1:sizeapx
    if Qc(k,10)==.02
        Qc(sizeapx+k,10)=.001;
    elseif Qc(k,10)==.001
        Qc(sizeapx+k,10)=.02;
    end
end

%Along x
px4  = 2*sizea*passox;

z2= (length(Qc(:,1))-mod(length(Qc(:,1)),px4))/px4;
for i = 1:(z2-1)
    Qc(i*px4+1:(i+1)*px4,10)=Qc(1:px4,10);
end
if ~isequal(mod(passox,py4),0)
    Qc(z*px4+1:z*px4+mod(length(Qc(:,1)),px4),10)=...
        Qc(1:mod(length(Qc(:,1)),px4),10);
end

QcH10= Qc(:,10);
re = A*QcH10;
for k = 1:6
    mcheck=tikhonov(U,diag(S),V,re,tik0_reg(:,k));
    Qc(:,10+k)=mcheck;
end
Qc(:,17)=A(1,:);
save Qc.mat
save Qc.txt Qc -ascii
% save Qc.mat QcTry Qc
% save -ascii QcV.txt QcTry
% save -ascii QcH.txt Qc

disp('END CHECKERBOARD TEST');

%% Save Images
% disp('Save Images');
% 
% folder = getfield(list,'folder');
% cur = strfind(folder,'\');
% name = [folder(cur(end)+5:end),'_', num2str(res), 'km'];
% 
% fName = ['OutImage/Grid',name];
% 
% fig = get(groot,'CurrentFigure');
% openFigs = fig.Number;
% 
% for outImg = 1:openFigs
% 	figs = ['figure (',num2str(outImg),')'];
%     saveas(figure(outImg), fullfile(fName, num2str(outImg)), 'jpg');
%     close(figure(outImg));
% end
% 
% disp('Saved!');
disp('END');
eT = datetime;
d = eT-sT

