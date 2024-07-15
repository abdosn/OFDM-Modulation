clear all;


%% OFDM SYSTEM INPUT PARAMs

nLevels = 11;

TxFileName = 'voicetest.wav';

RxFileName = 'recovered.wav';

nfft = input('Enter number of subcarriers : ');

Ncp = input('Enter number of cyclic prefix : ');

mod_order = input('1-->BPSK\n2-->QPSK\n3-->8PSK\n4-->16QAM\n5-->32QAM\n6-->64QAM\nEnter modulation order : ');


ch = input('1-->Noiseless\n2-->AWGN\n3-->Raylegih\nEnter channel type : ');

if ch == 3
   ch_taps = input('Enter no of channel taps : '); 
end

SNR = input('Enter SNR : ');


for mod_order = [1 2 3 4 5 6]
if mod_order == 1
    mod_method = 'BPSK';
end
if mod_order == 2
    mod_method = 'QPSK';
end
if mod_order == 3
    mod_method = '8PSK';
end
if mod_order == 4
    mod_method = '16QAM';
end
if mod_order == 5
    mod_method = '32QAM';
end
if mod_order == 6
    mod_method = '64QAM';
end





%sample rate ???
%freq spacing ???


%{
Rb = 100 ; % bit/sec
Tb = 1/Rb; % bit time
Time = [0:Tb:(Nbits-1)*Tb];
Time = Time.';
%}


%%  Transmitter

%generate data

[ser_data_bin , quan_Levels ,fs] = getvoicedata(TxFileName , nLevels);

%convert char type to int 

x = (ser_data_bin.')-48;

%symbol padding
 sym_rem = mod( mod_order - mod(length(x) , mod_order) , mod_order);
 padding = zeros(sym_rem,1);
 bin_padded = [x;padding];

bin_sym = reshape(bin_padded,mod_order,length(bin_padded)/mod_order)';
dec_sym = bin2dec(num2str(bin_sym));

    
    % ----- Modulation ----- %

%BPSK
if mod_order == 1
    mod_ind = 2^(mod_order-1); % =1
    n = 0:pi/mod_ind:2*pi-pi/mod_ind; % n [0 : pi : pi]
    in_phase = cos(n); 
    quadrature = sin(n); 
    symbol_book = (in_phase + quadrature*1i); % symbol_book [ cos0 +jsin0   cos180 + j sin180 ]
end

% QPSK & 8PSK
if mod_order == 2 || mod_order == 3
    
    mod_ind = 2^(mod_order-1); % 2 qpsk 4 8psk
    n = 0:pi/mod_ind:2*pi-pi/mod_ind; % n [0: pi/2 :3pi/2]    n [0 : pi/4 : 7pi/4 ]
    in_phase = cos(n+pi/4);  
    quadrature = sin(n+pi/4);
    symbol_book = (in_phase + quadrature*1i);
    
end

%16QAM, 64QAM
if mod_order == 4 || mod_order == 6
 mod_ind = sqrt(2^mod_order);
 %n = 0:pi/mod_ind:2*pi-pi/mod_ind;
 in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
 quadrature = repmat(linspace(-1,1,mod_ind)',1,mod_ind);
 symbol_book = (in_phase(:) + quadrature(:)*1i).';
end
%32QAM - Generates 6x6 constellation and removes corners
if mod_order == 5
 mod_ind = 6;
 %n = 0:pi/mod_ind:2*pi-pi/mod_ind;
 in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
 quadrature = repmat(linspace(-1,1,mod_ind)',1,mod_ind);
 symbol_book = (in_phase(:) + quadrature(:)*1i);
 symbol_book = symbol_book([2:5 7:30 32:35]).'; %corners are removed
end

% Modulate data
x_modulated = symbol_book(dec_sym+1);
 
fft_rem = mod(nfft-mod(length(x_modulated),nfft),nfft); 
X_padded = [x_modulated,zeros(1,fft_rem)]; 
X_blocks = reshape(X_padded,nfft,length(X_padded)/nfft);


% ifft of each subcarrier
IFFT_SIG = ifft(X_blocks);

%add cyclic prefix
cp =  IFFT_SIG( end-Ncp+1 : end , : );

OFDM_SIG = [cp;IFFT_SIG];


%Parallel to serial
OFDM_SIG = OFDM_SIG(:);

data_pwr = mean(abs(OFDM_SIG.^2));
%%
%PAPR = 10*log10((max((abs(OFDM_SIG)).^2)) / (mean((abs(OFDM_SIG)).^2))) % Change with number of bits --> ???

PAPR = 20*log10( peak2rms(OFDM_SIG) )


%% Channel

if ch == 1 % signal without noise
    Rec = OFDM_SIG;
    
end

if ch == 2 % AWGN
    sig_w_noise = awgn( OFDM_SIG , SNR , 'measured' );
    Rec = sig_w_noise;
end
    
if ch == 3 % Rayleigh
    
    noise_pwr = data_pwr/10^(SNR/10);
    noise = normrnd( 0 , sqrt(noise_pwr/2) , size(OFDM_SIG) ) + normrnd( 0 , sqrt(noise_pwr/2) , size(OFDM_SIG)) * 1i ;
    %sig_w_noise = awgn( OFDM_SIG , SNR , 'measured' );
    sig_w_noise = OFDM_SIG + noise;
    snr_meas = 10*log10(mean(abs(OFDM_SIG.^2))/mean(abs(noise.^2)));
   
    g = exp(-(0:ch_taps-1)); 
    g = g/norm(g); 
    sig_w_noise_fading = conv(sig_w_noise,g,'same');
    
    Rec = sig_w_noise_fading;
end

%% Recevier

%serial to parallel
Rec_s = reshape(Rec , nfft + Ncp ,length(Rec)/(nfft + Ncp));
%Remove cyclic prefix
Rec = Rec_s(Ncp+1 : end , :);

%fft
Rec_sym = fft(Rec);

% Estimate channel with least square method
G = Rec_sym(:,1) ./X_blocks(:,1) ;

Rec_sym2 = Rec_sym./repmat(G,1,size(Rec_sym,2));

Rec_sym = Rec_sym2(:);
Rec_sym = Rec_sym(1:end - fft_rem);



A=[real(symbol_book) imag(symbol_book)];

if (size(A,2)>2) %QPSK and 8PSK
    A=[real(symbol_book)' imag(symbol_book)'];
end


rec_syms = knnsearch(A,[real(Rec_sym) imag(Rec_sym)])-1;
rec_syms_cons = dec2bin(rec_syms);
rec_im_bin = reshape(rec_syms_cons',numel(rec_syms_cons),1); 
rec_im_bin = rec_im_bin(1:end-sym_rem);
y = rec_im_bin - 48;

%% 
[nerrbits , BER] = biterr(x,y)

%%

recvoicedata( RxFileName , (rec_im_bin.') , quan_Levels , nLevels , fs);


figure
subplot(1,2,1);
plot(x_modulated,'x','linewidth',2,'markersize',10);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In phase')
ylabel('Qudrature')
if ch == 3
 title(sprintf('\\bfTransmit Constellation\n\\rm%s Modulation\nMultipath Channel Taps: %d',mod_method,ch_taps));
else
 title(sprintf('\\bfTransmit Constellation\n\\rm%s Modulation\nAWGN',mod_method));
end
grid on
% Recovered constellation
subplot(1,2,2);
plot(Rec_sym(1:500:end),'x','markersize',3);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In phase')
ylabel('Qudrature')
if ch == 3
 title(sprintf('\\bfReceived Constellation\n\\rmMeasured SNR: %.2d dB',SNR));
else
 title(sprintf('\\bfReceived Constellation\n\\rmMeasured SNR: %.2d dB',SNR));
end
grid on
end