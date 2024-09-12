function [mjd2000] = et2MJD2000(et)

    % MJD=0 @ 1858-Nov-17 00:00 TDB
    % MJD2000=0 @ 2000-Jan-01 12:00 TDB
    
    et_str=cspice_et2utc(et,'ISOC',3);

    et_date=datetime(et_str,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS','TimeZone','UTC');

    jd=juliandate(et_date).';

    mjd=jd-2400000.5;

    mjd2000=mjd-51544; % (-0.5)
end