% 2017-08-13 20:06:54.321339192 +0200

A=[csarea(calib.zs0,calib.zb.A,1); csarea(calib.zs0,calib.zb.median,1)];
d=diff(A);
s = rms(d)
s./midrange(mean(A))

