import epygram
epygram.init_env()

r=epygram.formats.resource('analysis.surf-surfex.tl1798-c22.fa',openmode='r')
r.listfields()
lon=6.133242475917926
lat=59.40401122652407
var='SFX.TEMP_OC1'
f=r.readfield(var)
p=f.extract_point(lon,lat)
print(p.data)

