function [mjd2000] = et2MJD2000(et)
    
    et_str=cspice_et2utc(et,'ISOC',3);

    et_date=datetime(et_str,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS','TimeZone','UTC');

    jd=juliandate(et_date).';

    mjd2000=jd-2451544.5;
end