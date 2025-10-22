<<<<<<< HEAD
<<<<<<< HEAD
function t = otsu(S,nbins)
if nargin<2
    nbins = numel(S);
end

counts =hist(S,nbins);

p = counts/sum(counts);
omega = cumsum(p);
mu = cumsum(p .* (1:nbins));
mu_t = mu(end);
sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

maxval = max(sigma_b_squared);
idx = mean(find(sigma_b_squared == maxval));
% Normalize the threshold to the range [0, 1].
t = (idx - 1) / (nbins - 1);

=======
function t = otsu(S,nbins)
if nargin<2
    nbins = numel(S);
end

counts =hist(S,nbins);

p = counts/sum(counts);
omega = cumsum(p);
mu = cumsum(p .* (1:nbins));
mu_t = mu(end);
sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

maxval = max(sigma_b_squared);
idx = mean(find(sigma_b_squared == maxval));
% Normalize the threshold to the range [0, 1].
t = (idx - 1) / (nbins - 1);

>>>>>>> origin
=======
function t = otsu(S,nbins)
if nargin<2
    nbins = numel(S);
end

counts =hist(S,nbins);

p = counts/sum(counts);
omega = cumsum(p);
mu = cumsum(p .* (1:nbins));
mu_t = mu(end);
sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

maxval = max(sigma_b_squared);
idx = mean(find(sigma_b_squared == maxval));
% Normalize the threshold to the range [0, 1].
t = (idx - 1) / (nbins - 1);

>>>>>>> 3afa99beaec32453db712d4621fd05217f53fcfd
end