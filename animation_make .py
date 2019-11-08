#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from novas import compat as novas
from novas.compat  import eph_manager
from novas  import constants
import math
#轨迹
from mpl_toolkits.basemap import Basemap
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
#儒略日转成北京时间
def jdutc2bt(jd_utc):
    t=novas.cal_date(jd_utc+8/24 )
    a = int(t[3])
    b = int((t[3] - a) * 60)
    c = (t[3] - a - b / 60.0)*3600.0
    return (t,a,b,c) 
k=0
coor=[]
lons=[]
lats=[]
with open("eclipse.txt",'r') as f:
    strs=f.readlines()
    for str in strs:
        temp1=str.strip('\n')
        temp2=temp1.split(' ')
        #print(float(temp2[0]))
        if k%3==0:
            for i in temp2:
                lon=float(i)
                lons.append(lon)
        if k%3==1:
            for i in temp2:
                lat=float(i)
                lats.append(lat)
        if k%3==2:
            for i in temp2:
                jd_utc=float(i)
            coor.append([lons,lats,jd_utc])
            lons=[]
            lats=[]
        k=k+1
print(len(coor))
duration=int(len(coor)/10)
fig=plt.figure()
def make_frame(t):
    plt.cla()
    lons=coor[int(t/0.1)][0]
    lats=coor[int(t/0.1)][1]
    jd=coor[int(t/0.1)][2]
    time=jdutc2bt(jd-8/24)
    my_map=Basemap(projection='ortho',lat_0=0,lon_0=100)
    my_map.drawcoastlines()
    my_map.drawcountries()
    my_map.drawmapboundary(fill_color = 'aqua')
    #my_map.fillcontinents(color = 'coral', lake_color = 'aqua')
    x,y = my_map(lons, lats)
    my_map.scatter(x, y, marker='.',color='r')
    date=datetime.utcnow()
    date=date.replace(int(time[0][0]),int(time[0][1]),int(time[0][2]),int(time[1]),int(time[2]),int(time[3]))
    my_map.nightshade(date)
    time=jdutc2bt(jd)
    date=date.replace(int(time[0][0]),int(time[0][1]),int(time[0][2]),int(time[1]),int(time[2]),int(time[3]))
    plt.title('solar eclipse %s'% date.strftime("%d %b %Y %H:%M:%S"))
    return mplfig_to_npimage(fig)

animation = VideoClip(make_frame, duration=duration)
animation.ipython_display(fps=20, loop=True, autoplay=True)


# In[ ]:




