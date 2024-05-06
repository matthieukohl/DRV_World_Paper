function asym = Skew(w)
% Calculate the asym for w(x,y,t)
w_av = mean(mean(w,1),2);
w_prime = w-w_av;



asym = squeeze(mean(mean(w_prime.^3,1),2)./(mean(mean(w_prime.^2,1),2)).^(3/2));

end