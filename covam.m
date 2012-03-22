function r= covam(theta,spr,K)
% r=covam(theta,spr,K) returns a covariance matrix of user with
% parameteres as:
% 
% theta: user's azimuth angle in degrees
% spr: angle spread in degrees
% K: number of antennas




theta = theta*pi/180;
spr=spr*pi/180;
for p=1:K
    for q=1:K
        r(p,q)=exp(j*pi*(q-p)*sin(theta))*exp(-(pi*(q-p)*spr*cos(theta))^2/2);
    end
end



