function beating_freqs = check_paceable(time,n_paces,pacingfreq,filtered_Ca,fps)
L = length(time(time<=n_paces/pacingfreq));
zm_force = filtered_Ca(time<=n_paces/pacingfreq) - mean(filtered_Ca(time<=n_paces/pacingfreq));
NFFT = 2^nextpow2(L);
fc = fft(zm_force, NFFT)/L;
fc = 2*abs(fc(1:NFFT/2+1));
f = fps/2 * linspace(0, 1, NFFT/2+1);
[~, max_i] = max(fc);
beating_freqs = f(max_i);
end