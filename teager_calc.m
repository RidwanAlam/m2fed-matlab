function [tmean,tstd,tmax] = teager_calc(z)
    % calculates mean, max, and std of Teager energy
    % needs 'rTeagerCompute.m'
    tz = rTeagerCompute(z,2);
    tmean = mean(tz);
    tstd = std(tz);
    tmax = max(tz);
end