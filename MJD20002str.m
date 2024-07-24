function [et_str] = MJD20002str(mjd2000)

%     mjd=mjd2000+51544;
    jd=mjd2000+2451544.5;

    et_date=datetime(jd,'convertfrom','juliandate','Format','yyyy-MM-dd HH:mm:ss.SSS');

    et_str=char(et_date);

end