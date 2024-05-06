function asym = Lambda(w)
% Calculate the asym for w(x,y,t)
w_av = mean(mean(w,1),2);
w_prime = w-w_av;

w_up = max(w,0);
w_up_av = mean(mean(w,1),2);
w_up_prime = w_up - w_up_av;


asym = squeeze(mean(mean(w_prime.*w_up_prime,1),2)./(mean(mean(w_prime.^2))));


end

