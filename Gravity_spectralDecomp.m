clear all 
clc 
close all
format long
%%
%=========================== INPUT DATA DISINI ============================
AA = importdata('D:\8_TUGAS AKHIR\2_SIDANG PROGRESS\GGMPLUS\PETREL\Vertikal.txt');
[pp,qq] = size(AA);


%%
ab   =  AA(:,1); 
or   =  AA(:,2);
g   =  AA(:,3);

if ab(2) == ab(1)
    x = or;
elseif or(2) == or(1)
    x = ab;
end

n=5;
nx=length(x); N = nx; dt = abs(x(2)-x(1));
g=(g);

Fsf=1/dt;   freqf = Fsf*(1:ceil((N-1)/2))/(N); lmd = 1./freqf;  %memastikan nilai k dengan 1/dt
Fs=1;   freq = 0:Fs/nx:Fs/2; freq = freq(2:length(freq)); %normalized frequency
lambda=(round(abs(logspace(0.01,2,11)))); %z=(1./freq(lambda)/10)';
save cutoff.mat freqf freq lambda
% pause
Yf      = fft(g,N);
amp_Y   = (1/(Fsf*N))*10*log10(abs(Yf(1:floor(N/2))));
amp_Y   = 2*(amp_Y.^2);
% pause
% z=log(abs(Yf(1:floor(N/2))))./(-2*freqf'); z = -z';
z=lmd./10;

% pause
figure(1)
subplot(3,2,[5 6]), semilogx(freq,amp_Y,'*','DisplayName','spektral awal'), grid on, hold on, title('domain frekuensi') 
xlabel('normalized frequency (pi rad/sample)'),ylabel('magnitude (dB)')
subplot(3,2,[1 2 3 4]), plot(x,g,'o-','DisplayName','profil awal'), grid on, hold on, title('hasil filter')

figure(2)
subplot(2,2,[1 2]), plot(x,g), grid on, hold on, title('model awal')
subplot(2,2,3), plot(x,g,'o'),grid on, hold on
subplot(2,2,4), plot(x,g,'o'),grid on, hold on

figure(3)
semilogx(freq,amp_Y,'o-','DisplayName','Profil awal'), grid on, hold on, title('domain frekuensi')
gdout = [0 0 0 0];
%%
for i = 2:length(lambda)
    j = lambda(i)+1;
    if j >= length(freq)
        break
    end
    d = fdesign.lowpass('Fp,Fst,Ap,Ast',freq(lambda(i)),freq(lambda(i)+1),1,30);
    Hd = design(d,'ifir');
    fvtool(Hd); hold on
    zd = z(lambda(i))*ones(1,length(x));;
        for ii = 3:3:qq
        ab   = AA(:,ii-2); 
        or   = AA(:,ii-1);
        g   =  AA(:,ii);
        gL  = ones(1,2*length(x))*g(1); 
        gR  = ones(1,2*length(x))*g(end);
        g = g';
        g = [gL g gR];
        gg= filter(Hd,g);
        delay = floor(mean(grpdelay(Hd)));
        g=g(length(gL)+1+delay:end-length(gR)+delay);
        gg=gg(length(gL)+1+delay:end-length(gR)+delay);
        ggg=0;

        Yfg      = fft(gg,N);
        amp_Yfg  = (1/(Fsf*N))*10*log10(abs(Yfg(1:floor(N/2))));
        amp_Yfg  = 2*(amp_Yfg.^2);

        figure(2)
    %     subplot(2,2,3), plot(x,gg), grid on, hold on, title('lowpass filtered'),xlabel('offset (m)'), ylabel('mGal')
        subplot(2,2,4), plot(x,ggg), grid on, hold on, title('residual'),xlabel('offset (m)'), ylabel('mGal')

        figure(3)
        semilogx(freq,amp_Yfg,'-','DisplayName',['k =' num2str(freq(lambda(i)))]), legend show
        stem(freq(lambda(i)),amp_Y(lambda(i)),'ro','DisplayName',['k = ' num2str(freq(lambda(i)))]),xlabel('normalized frequency (pi rad/sample)'), ylabel('magnitude (dB)')
    %     gd = gg.*(((((gg+gg)./gg))-1))-z(i); %plotting gg, so so
    %     gd = max(gg).*(((((gg+max(gg))./max(gg)))-1))-z(i); %plotting max(gg), so so
    %     gd = (lmd(i).*(((gg+max(gg))./max(gg))-1))-z(i); %plotting lambda, bad
        gd = ((gg+max(gg))./max(gg))-1;
    %     zlist = -z(i)*ones(length(ab),1);
    %     flist = freqf(lambda(i))*ones(length(ab),1);
    %     gdo = [zlist flist ab or gd'];
    %     gdout = [gdout; gdo];
        figure(1), subplot(3,2,[1 2 3 4]), plot(x,gd,'DisplayName',['k = ' num2str(freqf(lambda(i)))]), title('hasil filter'), grid on, hold on, legend show
        xlabel('offset (m)'), ylabel('pseudodepth'), axis([min(x) max(x) -5000 300])
        subplot(3,2,[5 6]), stem(freq(lambda(i)),amp_Y(lambda(i)),'ro','DisplayName',['k = ' num2str(freq(lambda(i)))]), legend show
        figure(2)
        subplot(2,2,3), plot(x,gd), grid on, hold on, title('lowpass filtered'),xlabel('offset (m)'), ylabel('mGal')
        pause
        gdo = [ab or gd' zd'];
        gdout = [gdout;gdo];
    %     pause 
    %     break
        end
    gdout = gdout(2:end,:)
    outputfile=['C:\Users\ASUS\Documents\1_TUGASKULIAH\8_TUGAS AKHIR\2_SIDANG PROGRESS\OUTPUT\600km_wvy',num2str(freq(lambda(i))),'.txt']
    fid=fopen(outputfile,'wt');
    [m,n]=size(gdout);
        for vv=1:m
            fprintf(fid,'%8.3f %8.3f %10.3f %6.3f \n',gdout(vv,1),gdout(vv,2),gdout(vv,3), gdout(vv,4)); %ganti nama z nya jadi nama apa
        end
    gdout = [0 0 0 0];
%     pause
end