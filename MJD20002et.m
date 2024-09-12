function [et] = MJD20002et(mjd2000)

    % MJD=0 @ 1858-Nov-17 00:00 TDB
    % MJD2000=0 @ 2000-Jan-01 12:00 TDB
    
    mjd=mjd2000+51544; % (+0.5)

    jd=mjd+2400000.5;

    et_date=datetime(jd,'convertfrom','juliandate','Format','yyyy-MM-dd HH:mm:ss.SSS');

    et_str=char(et_date);

    et=cspice_str2et(et_str);
end

