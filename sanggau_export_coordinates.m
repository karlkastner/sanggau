% Thu  2 Nov 11:44:06 CET 2017

m = sanggau_metadata();
s = struct();
s.X = m.centre(1);
s.Y = m.centre(2);
s.zone     = '49N';
s.Geometry = 'Point';

filename = 'sanggau-coordinates';

Shp.export_gpx(s,filename) 
Shp.write(s,filename);

tar -cf sanggau-discharge-2017-11-02.tar sanggau-coordinates/ -C /home/pia/large/phd/src/discharge/mat/ sanggau-stage-discharge-2017-07-21.csv

