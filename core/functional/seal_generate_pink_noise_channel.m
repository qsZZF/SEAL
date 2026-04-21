%% pink noise
function noise_channel = seal_generate_pink_noise_channel(n_samples)
    if n_samples == 0, noise_channel = []; return; end
    n_fft = 2^nextpow2(n_samples); 
    
    f = linspace(0, 1, n_fft/2 + 1); 
    inv_f_spectrum = zeros(size(f));
    inv_f_spectrum(f > 0) = 1./sqrt(f(f > 0)); 
    inv_f_spectrum(1) = 0; 

    random_phases = exp(1i * 2 * pi * rand(size(f)));
    half_spectrum = inv_f_spectrum .* random_phases;
    
    full_spectrum = zeros(1, n_fft);
    full_spectrum(1:n_fft/2+1) = half_spectrum;
    if n_fft/2+2 <= n_fft % Ensure indices are valid for fliplr
        full_spectrum(n_fft/2+2:end) = conj(fliplr(half_spectrum(2:end-1)));
    end

    noise_channel_temp = real(ifft(full_spectrum, n_fft));
    noise_channel = noise_channel_temp(1:n_samples); 
    noise_channel = noise_channel(:)'; 
end