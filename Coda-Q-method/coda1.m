close all; clear all;

myFolder = './Seizmogrami_STON'; % folder with seismograms
filePattern = fullfile(myFolder, 'STON*.txt');
theFiles = dir(filePattern);
delimiter = '\t';

for i = 1 : length(theFiles)
        
    filename = 'Popis_potresa_STON.txt'; % list of earthquakes
    P = importdata(filename,delimiter,1);
    YYYY1 = P.data(i, 2); % hypocenter time
    MM1 = P.data(i, 3); 
    DD1 = P.data(i, 4); 
    HH1 = P.data(i, 5);
    mm1 = P.data(i, 6); 
    ss1 = P.data(i, 7); 
    h = P.data(i,10); % hypocenter depth
    e = P.data(i,20); % epicentral distance
    r = sqrt(h^2 + e^2); % hypocentral distance

    baseFileName = theFiles(i).name;
    fullFileName = fullfile(myFolder, baseFileName);
    P = importdata(fullFileName, delimiter,1);    

    YYYY = P.data(1,1); % beginning of the record
    MM = P.data(2,1);
    DD = P.data(3,1);
    HH = P.data(4,1);
    mm = P.data(5,1);
    dt = P.data(7, 1);

    N_S = P.data(8:end, 2); % components of interest
    E_W = P.data(8:end, 3);

    N_S = detrend(N_S); % removing linear trend
    E_W = detrend(E_W);

    N_S = detrend(N_S, 'constant'); % removing mean
    E_W = detrend(E_W, 'constant');

    n = length(N_S); % time axis for event
    ss = ones(n,1); 
    ss(1) = P.data(6,1);
    t = ones(n,1);
    t(1) = datenum(YYYY,MM,DD,HH,mm,ss(1));
    for j = 1:(n-1)
        ss(j+1)= ss(j) + dt;
        t(j+1) = datenum(YYYY,MM,DD,HH,mm,ss(j));
        if t(j+1) == datenum(YYYY1,MM1,DD1,HH1,mm1,ss1)
            N = j+1;
        end
    end 
    hip = linspace(-N*dt,n*dt,n); % hypocentral time

    ts = r / 3.5; % nastup S-faze
    tp = r / (3.5 * sqrt(3)); % P-phase pick

    fs = [1.5 3 6 9 12 18]; % center frequency
    fN = 1/(2*dt); % Nyquist frequency
    fd = [1 2 4 6 8 12]; % lower frequency band
    fg = [2 4 8 12 16 24]; % upper frequency band

    for j = 1:length(fs)
    
        if (4*fs/3) < fN
        
            [b,a] = butter(4,[fd(j) fg(j)]/fN); % Butterworth filter
            N_S1 = filter(b, a, N_S);
            E_W1 = filter(b, a, E_W);

            N_S1 = abs(hilbert(N_S1)); % Hilbert transform envelope
            E_W1 = abs(hilbert(E_W1));

            N_S2 = smooth(N_S1, 5/dt, 'moving'); % smoothing
            E_W2 = smooth(E_W1, 5/dt, 'moving');

            nemir1 = N + round(tp-5)/dt; % noise on each component
            nemir2 = N + round(tp)/dt; 
            noiseNS = mean([N_S2(nemir1:nemir2)]);
            noiseEW = mean([E_W2(nemir1:nemir2)]); 


            for k = 1:2
            
                if k == 1
                    qB = N + (2*ts)/dt; % beginning of the window 1 for coda
                else
                    qB = N + (2*ts + 10)/dt; % beginning of the window 2 for coda
                end
                qE = qB + 25/dt; % ending of the windows for coda
                
                S1 = mean([N_S2((qE-5/dt):qE)]); % mean noise
                S2 = mean([E_W2((qE-5/dt):qE)]);
                
                if S1 >= 2*noiseNS % removing noise
                    N_S2 = N_S2 - noiseNS;
                end
                
                if S2 >= 2*noiseEW
                    E_W2 = E_W2 - noiseEW;
                end

                lnAtN_S = log(N_S2.*hip'); % logarithm of the envelope*time
                lnAtE_W = log(E_W2.*hip');
            
                y1 = [lnAtN_S(qB:qE)]; % values
                y2 = [lnAtE_W(qB:qE)];
                window = [hip(qB:qE)]; % window for coda

                complete = [y1; y2]; % vector values
                wwindow = [window'; window']; 
                p = polyfit(wwindow,complete,1); % linear regression

                Q(k,j) = -pi * fs(j)/p(1); % Q-factor
                check = isreal(Q(k,j));
                
                if check == false || Q(k,j) < 0 % removing negative and imaginary values
                    Q(k,j) = NaN;
                end
            end
        end
    end
    
    if i == 1 % saving coda Q-factors
        Qc = Q;
    else
        Qc = [Qc, Q];
    end
end

fileID = fopen('faktori.txt','w'); % writing values to text
fprintf(fileID,'%s\t%s\r\n','1. prozor','2. prozor');
fprintf(fileID,'%.8f\t%.8f\r\n',Qc);
fclose(fileID);