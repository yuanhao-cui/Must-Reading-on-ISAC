function [QPSK_symbols] = QPSK_mapper(bitseq)

%MIMO-OFDM Wireless Communications with MATLAB¢ç   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%2010 John Wiley & Sons (Asia) Pte Ltd

QPSK_table = [1+j -1+j 1-j -1-j]/sqrt(2);
for ii=1:length(bitseq)/2
    temp=2*bitseq(ii*2-1)+bitseq(ii*2);
    QPSK_symbols(ii)=QPSK_table(temp+1);
end