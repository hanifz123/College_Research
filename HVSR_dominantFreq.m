clear all
clc
close all
data = importdata('D:\MAYA TA\DATA\ascii\36');
[p,q] = size(data); 
    tmax = 30; %10 12 13 14 16 17 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
%     tmax = 15; %06 07 08 09 11
%     tmax = 28; %18
%     tmax = 29; %15
t = tmax*60; 
    tt = 0.004;  %06 07 08 09 10 11 12 13 14 16 17 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
%     tt = 0.948; %18
%     tt = 37.1; %15
t = t+tt;
    fs = 250; 
Ts = 1/fs;
tseries = linspace(0,t,fs*t);

figure(1), hold on
z = data(:,1); subplot(3,2,1), plot(tseries,z), xlabel('sekon'), ylabel('Vertikal (Z)'), grid on
title('Time Series Mikrotremor')
n = data(:,2); subplot(3,2,3), plot(tseries,n), xlabel('sekon'), ylabel('Northing (N)'), grid on
e = data(:,3); subplot(3,2,5), plot(tseries,e), xlabel('sekon'), ylabel('Easting (E)'), grid on

% figure(2)
% pspectrum(z,fs,'FrequencyLimits',[0 1])
% 
% figure(3)
% [pzz,f]=pspectrum(z,fs,'FrequencyLimits',[0 1]);plot(f,10*log10(pzz)),grid on
% % 
wdw = 200; index = wdw*fs;
nlength = 1:index:p-index;

loop = 0;
for i = 1:index:p-index
loop = loop + 1;
figure(2);
zz = z(i:i+(index-1));    
[pzz,f]=pspectrum(zz,fs,'FrequencyLimits',[0 1]);
[pzz1,f1]=pspectrum(zz,fs,'FrequencyLimits',[1 10]);
pzz = [pzz;pzz1]; f = [f;f1];  yzp =pzz; 
% yzp = smoothSpectra(pzz); 
% minf = 1; maxf = 12;

nn = n(i:i+(index-1));    
[pnn,f] = pspectrum(nn,fs,'FrequencyLimits',[0 1]);
[pnn1,f1]=pspectrum(nn,fs,'FrequencyLimits',[1 10]);
pnn = [pnn;pnn1]; f = [f;f1]; ynp = pnn;
% ynp = smoothSpectra(pnn);

ee = e(i:i+(index-1));    
[pee,f] = pspectrum(ee,fs,'FrequencyLimits',[0 1]);
[pee1,f1]=pspectrum(ee,fs,'FrequencyLimits',[1 10]);
pee = [pee;pee1]; f = [f;f1]; yep = pee;
% yep = smoothSpectra(pee);
% f = f*((fs/2)/3.1416);

hvp = sqrt((ynp.^2+yep.^2)./2)./yzp; %hvp = 1/(fs*length(hvp))*hvp.^2;
hvp = smoothSpectra(hvp); hvp = abs(hvp');
% fdiv = sum(hvp(minf:maxf).*f(minf:maxf));
% divmean = sum((hvp(minf:maxf)));
fdiv = sum(hvp.*f); %freqlimit
divmean = sum(hvp);%freqlimit
fdiv = fdiv/divmean;
fmean(loop) = fdiv;

subplot(3,6,1), semilogx(f,10*log10(yzp)), hold on, grid on
title('PSD Vertikal (z)'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10])
xlabel('frequency (Hz)'),ylabel('PSD Vertikal (z) G^2/Hz'),xlim([0.2 10])
subplot(3,6,7), semilogx(f,10*log10(ynp)), hold on, grid on
title('PSD Northing (n)'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10])
xlabel('frequency (Hz)'),ylabel('PSD Northing (n) G^2/Hz'),xlim([0.2 10])
subplot(3,6,13), semilogx(f,10*log10(yep)), hold on, grid on
title('PSD Easting (e)'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10])
xlabel('frequency (Hz)'),ylabel('PSD Easting (e) G^2/Hz'),xlim([0.2 10])
subplot(3,6,[2 3 8 9 14 15]), semilogx(f,10*log10(hvp)), hold on, grid on
title('Kurva H/V Squared Average'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10]); %set(gca, 'YTick', [0:1:4]);
xlabel('frequency (Hz)'),ylabel('H/V G^2/Hz'),xlim([0.2 10])
% clear pnn;clear pee; clear pzz; clear hvp
end

f0 = mean(fmean,'all')
fstd = std(fmean)
subplot(3,6,[2 3 8 9 14 15]); f01 = xline(f0,'-',strcat('f0=  ',num2str(f0)),'LineWidth', 1.5);f02 = xline(f0+fstd, '-',strcat('±',num2str(fstd)),'LineWidth', 1);f03 = xline(f0-fstd, '-','LineWidth', 1);


d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.005,0.016,1,120);
Hd = design(d,'ifir');

tukey0 = tukeywin(length(data),0.05);
[tukey,ft]=pspectrum(tukey0,fs,'FrequencyLimits',[0 1]);
[tukey1,ft1]=pspectrum(tukey0,fs,'FrequencyLimits',[1 10]);
tukey = [tukey;tukey1]; f = [ft;ft1];

figure(1), hold on
zf = filter(Hd,z);subplot(3,2,2), plot(tseries,zf), xlabel('sekon'), ylabel('Vertikal (Z)'), grid on
title('Filtered Time Series'); 
zf = round(zf); 
% zf = filter(tukey,zf);
nf = filter(Hd,n);subplot(3,2,4), plot(tseries,nf), xlabel('sekon'), ylabel('Northing (N)'), grid on;
nf = round(nf);
% nf = filter(tukey,nf);
ef = filter(Hd,e);subplot(3,2,6), plot(tseries,ef), xlabel('sekon'), ylabel('Easting (E)'), grid on
ef = round(ef);
% ef = filter(tukey,ef);
hvsro = [zf nf ef];

loop = 0;
for i = 1:index:p-index
loop = loop+1;
figure(2);
zz = zf(i:i+(index-1));  
[pzz,f] = pspectrum(zz,fs,'FrequencyLimits',[0 1]); 
[pzz1,f1]=pspectrum(zz,fs,'FrequencyLimits',[1 10]);
pzz = [pzz;pzz1]; f = [f;f1];
pzz = pzz.*tukey;
yzpp = pzz;
% yzpp = smoothSpectra(pzz); 
% minf = find(f==0); maxf = find(f >=9.99 & f <=10.02);

nn = nf(i:i+(index-1));    
[pnn,f] = pspectrum(nn,fs,'FrequencyLimits',[0 1]); 
[pnn1,f1]=pspectrum(nn,fs,'FrequencyLimits',[1 10]);
pnn = [pnn;pnn1]; f = [f;f1];
pnn = pnn.*tukey;
ynpp = pnn;
% ynpp = smoothSpectra(pnn);

ee = ef(i:i+(index-1));    
[pee,f] = pspectrum(ee,fs,'FrequencyLimits',[0 1]); 
[pee1,f1]=pspectrum(ee,fs,'FrequencyLimits',[1 10]);
pee = [pee;pee1]; f = [f;f1];
pee = pee.*tukey;
yepp = pee;
% yepp = smoothSpectra(pee);
% f = f*((fs/2)/3.1416);

hvf = sqrt((ynpp.^2+yepp.^2)./2)./yzpp; %hvpf = 1/(fs*length(hvf))*hvf.^2;
hvf = smoothSpectra(hvf); hvf = abs(hvf');
% fdivf = sum(hvf(minf:maxf).*f(minf:maxf));
% divmeanf = sum((hvf(minf:maxf)));
fdivf = sum(hvf.*f); %freqlimit
divmeanf = sum(hvf);%freqlimit
fdivf = fdivf/divmeanf;
fmeanf(loop) = fdivf;

subplot(3,6,4), semilogx(f,10*log10(yzpp)), hold on, grid on
title('PSD Vertikal (z)'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10])
xlabel('frequency (Hz)'),ylabel('PSD Vertikal (z) G^2/Hz'),xlim([0.2 10])
subplot(3,6,10), semilogx(f,10*log10(ynpp)), hold on, grid on
title('PSD Northing (n)'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10])
xlabel('frequency (Hz)'),ylabel('PSD Northing (n) G^2/Hz'),xlim([0.2 10])
subplot(3,6,16), semilogx(f,10*log10(yepp)), hold on, grid on
title('PSD Easting (e)'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10])
xlabel('frequency (Hz)'),ylabel('PSD Easting (e) G^2/Hz'),xlim([0.2 10])
subplot(3,6,[5 6 11 12 17 18]), semilogx(f,10*log10(hvf)), hold on, grid on
title('Kurva H/V Squared Average'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10]); %set(gca, 'YTick', [0:1:4]);
xlabel('frequency (Hz)'),ylabel('H/V G^2/Hz'),xlim([0.2 10])

figure(3),semilogx(f,10*log10(hvf)), hold on, grid on
title('Kurva H/V Squared Average'); set(gca, 'XTick', [0.6 0.8 1 2 4 6 8 10]); %set(gca, 'YTick', [0:1:4]);
xlabel('frequency (Hz)'),ylabel('H/V G^2/Hz'),xlim([0.2 10])
end
% 
f0fil = mean(fmeanf,'all')
fstdfil = std(fmeanf)
figure(2);subplot(3,6,[5 6 11 12 17 18]); f01f = xline(f0fil,'-',strcat('f0=  ',num2str(f0fil)),'LineWidth', 1.5);f02f = xline(f0fil+fstdfil, '-',strcat('±',num2str(fstdfil)),'LineWidth', 1);f03f = xline(f0fil-fstdfil, '-','LineWidth', 1);

% 
% % fvtool(Hd); hold on