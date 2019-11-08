#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
#import datetime
#时间初始值设置
year=2019
month=12
day=26
hour=2
#儒略日转成北京时间
def jdutc2bt(jd_utc):
    t=novas.cal_date(jd_utc+8/24 )
    a = int(t[3])
    b = int((t[3] - a) * 60)
    c = (t[3] - a - b / 60.0)*3600.0
    return (t,a,b,c) 
#叉乘计算,返回一个单位矢量
def chacheng(O_E,O_S):
    O_U=[]
    s=0
    O_U.append(O_E[1]*O_S[2]-O_E[2]*O_S[1])
    O_U.append(O_E[2]*O_S[0]-O_E[0]*O_S[2])
    O_U.append(O_E[0]*O_S[1]-O_E[1]*O_S[0])
    for i in range(3):
        s+=O_U[i]*O_U[i]
    s=math.sqrt(s)
    for i in range(3):
        O_U[i]=O_U[i]/s
    return O_U
#地理位置
    #佘山
latitude=31.084451#纬度
longitude=121.165869#经度
height=12.2#海拔
    #老家
#latitude=34.173863#纬度
#longitude=115.927278#经度
#height=42#海拔

tempetature=25.0#温度
pressure=1010.0#气压
#距离以AU为单位
#re = (constants.ERAD + 65000.0 )/ constants.AU#地球半径，考虑地球大气层的影响
re = (constants.ERAD) / constants.AU#地球半径，不考虑地球大气层的影响
rs = 696000000.0 / constants.AU#太阳半径
rm = 1738000.0 / constants.AU#月亮半径
ta = constants.AU / constants.C#从太阳到地球的光行时
dtheta=0.0
leap_second=37#闰秒
tt_tai=32.184#tt和国际原子时tai之差
#EOP参数
ut1_utc=0.06723
x_pole = -0.002
y_pole = +0.529
#tt-ut1
delta_t=tt_tai+leap_second-ut1_utc
#构造不同的儒略日时间变量
jd_utc=novas.julian_date(year,month,day,hour)#utc
jd_ut1=jd_utc+ut1_utc/86400.0#ut1
jd_tt=jd_utc+(leap_second+tt_tai)/86400.0#tt
jd=(jd_tt,0.0)
jd0=(jd_tt-ta/86400,0.0)
#打开de历表
jd_s,jd_e,num=eph_manager.ephem_open()
 #太阳和地球构造
sun=novas.make_object(0,10,'sun',None)
moon=novas.make_object(0,11,'moon',None)
earth=novas.make_object(0,3,'earth',None)
 #位置构造
location1=novas.make_on_surface(latitude,longitude,height,25.0,1013)
location=novas.make_observer_on_surface(latitude,longitude,height,25.0,1013)
#矢量和夹角初始化
O1_S1=[0.0,0.0,0.0]#月球半影锥点到太阳质心矢量坐标
O2_S1=[0.0,0.0,0.0]#月球全影锥点到太阳质心矢量坐标
O1_E=[0.0,0.0,0.0]#月球半影锥点到地心距离矢量坐标
O1_T=[0.0,0.0,0.0]#月球半影锥点到地球某一点距离矢量坐标
O2_E=[0.0,0.0,0.0]#月球全影锥点到地心距离矢量坐标
O2_T=[0.0,0.0,0.0]#月球全影锥点到地球某一点距离矢量坐标
O1_M=[0.0,0.0,0.0]#月球半影锥点到月心距离矢量坐标
O2_M=[0.0,0.0,0.0]#月球全影锥点到月心距离矢量坐标
#求食分
S1_O1=[0.0,0.0,0.0]#太阳质心到月球全影锥点矢量坐标
S1_E=[0.0,0.0,0.0]#太阳质心到地心矢量坐标
T_S1=[0.0,0.0,0.0]#地球上一点到太阳质心距离矢量坐标
T_M=[0.0,0.0,0.0]#地球上一点锥点到月心距离矢量坐标
#求地球上满足偏食条件的点
O1_X=[0.0,0.0,0.0]#垂直于太阳质心，地心，月心所在平面的单位矢量坐标
O1_Y=[0.0,0.0,0.0]#垂直于太阳质心，O1_X，月心所在平面的单位矢量坐标
O1_O=[0.0,0.0,0.0]#月球半影锥点到锥轴上一点（O点，刚进入偏食阶段地球上的点垂直于锥轴的点）距离矢量坐标
O_T=[0.0,0.0,0.0]#O点到刚进入偏食阶段地球上的点距离矢量坐标
#判断一点是否处于白天
E_O4=[0.0,0.0,0.0]
#经纬度转为ITRS下坐标
E_T=[re*math.cos(latitude*constants.DEG2RAD)*math.cos(longitude*constants.DEG2RAD),	 re*math.cos(latitude*constants.DEG2RAD)*math.sin(longitude*constants.DEG2RAD),	 re*math.sin(latitude*constants.DEG2RAD)     ]
#ITRS->GCRS
E_TT=novas.ter2cel(jd_ut1, 0.0, delta_t, x_pole, y_pole,E_T, 1, 0,0)#
thetaE_O1=0.0#地球半径在半影锥点O1角距
thetaM_O1=0.0#月球半径在半影锥点O1角距
thetaE_O2=0.0#地球半径在全影锥点O2角距
thetaM_O2=0.0#月球半径在全影锥点O2角距
thetaEM_O1=10.0#地心和月心在月球半影锥点为原点的夹角
thetaTM_O1=10.0#地球某一点和月心在月球半影锥点为原点的夹角
thetaEM_O2=1.0#地心和月球在月心全影锥点为原点的夹角
thetaTM_O2=1.0#地球某一点和月球在月心全影锥点为原点的夹角
O2_Elength=1#月球全影锥点可能在地球内部，需要判断，先赋初值足够大
O2E_re=1#月球全影锥点到地心距离和地球半径之差
#theta=constants.TWOPI
flag1=0#日偏食开始结束标志，用来结束循环
flag2=0#日全食发生标记
flag3=0#日环食标记
flag4=0#食甚标记
flag_a=0#日环食时月球全影锥点进出地球标记
#某一点
flag11=0#日偏食开始结束标志，用来结束循环
flag22=0#日全食发生标记
flag33=0#日环食标记
flag44=0#食甚标记
jdd=jd_utc
#经纬度转为itrs下坐标
def geo2itrs(longitude,latitude):
    E_T=[re*math.cos(latitude*constants.DEG2RAD)*math.cos(longitude*constants.DEG2RAD),	 re*math.cos(latitude*constants.DEG2RAD)*math.sin(longitude*constants.DEG2RAD),	 re*math.sin(latitude*constants.DEG2RAD)]
    E_TT=novas.ter2cel(jd_ut1, 0.0, delta_t, x_pole, y_pole,E_T, 1, 0,0)
    return E_TT
#判断地球上某一地点是否是白天
def day_time(longitude,latitude,E_O4):
    ETT=geo2itrs(longitude,latitude)
    vpd=0
    ETTlength=0
    E_O4length=0
    for i in range(3):
        vpd+=ETT[i]*E_O4[i]
        ETTlength+=ETT[i]*ETT[i]
        E_O4length+=E_O4[i]*E_O4[i]
    ETTlength=math.sqrt(ETTlength)
    E_O4length=math.sqrt(E_O4length)
    theta=math.acos(vpd/(ETTlength*E_O4length))
    if theta>(math.pi/2-math.asin(rs-re)):
        return True
    else:
        return False
coor=[]#存储时间点数据
f=open("eclipse.txt",'w')
while True:
    #partial
#准备工作
    dthetaEM_O1=thetaEM_O1-thetaE_O1-thetaM_O1
    dthetaTM_O1=thetaTM_O1-thetaM_O1
    #total
    dthetaEM_O2=thetaEM_O2-thetaE_O2-thetaM_O2
    dthetaTM_O2=thetaTM_O2-thetaM_O2
    #annular
    dthetaEM_O2A=thetaEM_O2+thetaE_O2+thetaM_O2-math.pi
    dthetaTM_O2A=thetaTM_O2+thetaM_O2-math.pi
    O2E_re=math.fabs(O2_Elength-re)#在地球内部设其值为1
    thetaEM=thetaEM_O1#上一步的地月夹角，食甚时夹角最小
    thetaTM=thetaTM_O1#上一步的地月夹角，食甚时夹角最小
    jd_utc+=0.1/86400#每次增加0.1秒
    jd_ut1=jd_utc+ut1_utc/86400.0#ut1
    jd_tt=jd_utc+(leap_second+tt_tai)/86400.0#tt
    jd=(jd_tt,0.0)#tt代替tdb，差别不大
    jd0=(jd_tt-ta/86400,0.0)#太阳发出光时的时间
    pos_earth0=novas.ephemeris(jd0,earth)#太阳发出光时的地球icrs坐标
    pos_earth=novas.ephemeris(jd,earth)#太阳到达地球时的地球icrs坐标
    pos_moon0=novas.ephemeris(jd0,moon)#太阳发出光时的月球icrs坐标
    pos_moon=novas.ephemeris(jd,moon)#太阳到达地球时的月球icrs坐标
    vpd_O1EM=0.0#O1E和O1M矢量积
    vpd_O1TM=0.0#O1T和O1M矢量积
    vpd_O2EM=0.0#O2E和O2M矢量积
    vpd_O2TM=0.0#O2T和O2M矢量积
    O1_Elength=0.0#月球半影锥点到地心距离
    O1_Tlength=0.0#月球半影锥点到地球某一点距离
    O1_Mlength=0.0#月球半影锥点到月心距离
    O2_Elength=0.0#月球全影锥点到地心距离
    O2_Tlength=0.0#月球全影锥点到地球某一点距离
    O2_Mlength=0.0#月球全影锥点到月心距离
    #ITRS->GCRS
    E_TT=novas.ter2cel(jd_ut1, 0.0, delta_t, x_pole, y_pole,E_T, 1, 0,0)#
    for i in range(3):
        O1_S1[i] = rm / (rm+rs) * pos_moon0[0][i]-pos_moon[0][i]
        O2_S1[i] = rm/  (rm-rs) * pos_moon0[0][i]-pos_moon[0][i] 
        O1_M[i] = O1_S1[i] + pos_moon[0][i]
        O2_M[i] = O2_S1[i] + pos_moon[0][i]
        O1_E[i] = O1_S1[i] + pos_earth[0][i]
        O1_T[i] = O1_E[i]+E_TT[i]#月球半影锥点到地球上某一点的矢量
        O2_E[i] = O2_S1[i] + pos_earth[0][i]
        O2_T[i] = O2_E[i]+E_TT[i]#月球全影锥点到地球上某一点的矢量
        E_O4[i] = re/(rs-re)*pos_earth0[0][i]#用于判断某一点是否处于白天
        O1_Elength+=O1_E[i]*O1_E[i]
        O1_Tlength+=O1_T[i]*O1_T[i]#月球半影锥点到地球上某一点的矢量长度
        O1_Mlength+=O1_M[i]*O1_M[i]
        O2_Elength+=O2_E[i]*O2_E[i]
        O2_Tlength+=O2_T[i]*O2_T[i]#月球全影锥点到地球上某一点的矢量长度
        O2_Mlength+=O2_M[i]*O2_M[i]
        vpd_O1EM+=O1_E[i]*O1_M[i]
        vpd_O1TM+=O1_T[i]*O1_M[i]
        vpd_O2EM+=O2_E[i]*O2_M[i]
        vpd_O2TM+=O2_T[i]*O2_M[i]
    #矢量的长度
    O1_Elength=math.sqrt(O1_Elength)#O1_E长度
    O1_Tlength=math.sqrt(O1_Tlength)#O1_T长度
    O1_Mlength=math.sqrt(O1_Mlength)#O1_M长度
    O2_Elength=math.sqrt(O2_Elength)#O2_E长度
    O2_Tlength=math.sqrt(O2_Tlength)#O2_T长度
    O2_Mlength=math.sqrt(O2_Mlength)#O2_M长度
    #地球和月球在O1,O2点为原点的半径的角距
    thetaE_O1=math.asin(re/O1_Elength)
    thetaM_O1=math.asin(rm/O1_Mlength)
    if O2_Elength<re:#全影锥点在地球内部
        thetaE_O2=math.pi/2#在地球内部设其值为1
    else:
        thetaE_O2=math.asin(re/O2_Elength)
    thetaM_O2=math.asin(rm/O2_Mlength)
    #矢量O1E,O1M的夹角和矢量O2E,O2M的夹角
    thetaEM_O1=math.acos(vpd_O1EM/(O1_Elength*O1_Mlength))
    thetaTM_O1=math.acos(vpd_O1TM/(O1_Tlength*O1_Mlength))
    thetaEM_O2=math.acos(vpd_O2EM/(O2_Elength*O2_Mlength))
    thetaTM_O2=math.acos(vpd_O2TM/(O2_Tlength*O2_Mlength))
#开始计算 
    #日偏食
    if dthetaEM_O1*(thetaEM_O1-thetaE_O1-thetaM_O1)<0: #（1）flag1==flag4==0：日偏食开始之前
        if flag1==0:                                                             #（2）flag1==flag4==1：食甚之后，日偏食结束之前
            t1=jdutc2bt(jd_utc)#开始
            flag1+=1
        elif flag1==1:
            t2=jdutc2bt(jd_utc)#结束
            break#日偏食结束，日食结束
     #日环食
    if dthetaEM_O2A*(thetaEM_O2+thetaE_O2+thetaM_O2-math.pi)*O2E_re*(O2_Elength-re)<0 and flag_a==0:#两个判断条件不能同时满足 ，防止后者引发前者的发生
        if flag3==0:
            t5=jdutc2bt(jd_utc)#开始
            flag3+=1
        elif flag3>0:
            t6=jdutc2bt(jd_utc)#结束
            flag3+=1
        if O2E_re*(O2_Elength-re)<0:
            flag_a+=1
    #日全食
    if dthetaEM_O2*(thetaEM_O2-thetaE_O2-thetaM_O2)<0 or (O2E_re*(O2_Elength-re)<0):
        if flag2==0:
            t3=jdutc2bt(jd_utc)#开始
            flag2+=1
        elif flag2>0:
            t4=jdutc2bt(jd_utc)#结束
            flag2+=1
    #食甚
    if thetaEM<thetaEM_O1 and flag4==0:#食甚按照半影影锥锥点来算，因为全影影锥和半影影锥的大小变化、同时达到最值的时间大小不一定相同
        t7=jdutc2bt(jd_utc)#食甚
        flag4+=1
        #计算食分
        S1_O1length=0.0
        S1_Elength=0.0
        T_S1length=0.0
        T_Mlength=0.0
        E_Tlength=0.0
        vpd_S1O1E=0.0
        vpd_TS1M=0.0
        for i in range(3):
            S1_O1[i] = pos_moon0[0][i]*(rs/(rs+rm))
            S1_E[i] = S1_O1[i]+O1_E[i]
            S1_Elength+=S1_E[i]*S1_E[i]
            vpd_S1O1E+=S1_O1[i]*S1_E[i]
        S1_O1length=O1_Mlength*rs/rm
        S1_Elength=math.sqrt(S1_Elength)
        if flag2!=0 or flag3!=0:#发生日环食或者日全食
            ratio=(vpd_S1O1E-math.sqrt(vpd_S1O1E**2-S1_O1length**2*(S1_Elength**2-re**2)))/(S1_O1length**2)
            thetaM_T=math.asin(rm/(ratio*S1_O1length-O1_Mlength*(rs+rm)/rm))
            thetaS_T=math.asin(rs/(ratio*S1_O1length))
            mag=(thetaM_T+thetaS_T)/(2*thetaS_T)#食分
        else:#只发生了偏食
            ratio=vpd_S1O1E/(S1_O1length*S1_O1length)
            for i in range(3):#这里的E_T指从地心垂直于锥轴于点T的矢量
                E_Tlength+=(ratio*S1_O1[i]-S1_E[i])*(ratio*S1_O1[i]-S1_E[i])
            E_Tlength=math.sqrt(E_Tlength)
            for i in range(3):
                T_S1[i]=-(re*ratio/E_Tlength*S1_O1[i]+(1-re/E_Tlength)*S1_E[i])
                T_M[i]=T_S1[i]+pos_moon0[0][i]
                vpd_TS1M+=T_S1[i]*T_M[i]
                T_S1length+=T_S1[i]*T_S1[i]
                T_Mlength+=T_M[i]*T_M[i]
            T_S1length=math.sqrt(T_S1length)
            T_Mlength=math.sqrt(T_Mlength)
            thetaM_T=math.asin(rm/T_Mlength)
            thetaS1_T=math.asin(rs/T_S1length)
            thetaT_S1M=math.acos(vpd_TS1M/(T_S1length*T_Mlength))
            mag=(thetaM_T+thetaS1_T-thetaT_S1M)/(2*thetaS1_T)#食分
    if flag1!=0 and jd_utc-jdd>1/60/24:
        jdd=jd_utc
        #print("进入")
        lons=[]#经度
        lats=[]#纬度
        E_Tlength=0
         #计算轨迹点
        O1_Y=chacheng(O1_E,O1_M)
        O1_X=chacheng(O1_Y,O1_M)
        rat=-O1_Mlength*math.sin(thetaM_O1)/math.cos(thetaM_O1)
        O1_Tlength=0.0#月球半影锥点到地球某一点距离
        vpd_O1TE=0.0#O1T和O1M矢量积
        for i in range(3):
            O1_T[i] = O1_M[i]+rat*O1_X[i]
            O1_Tlength+=O1_T[i]*O1_T[i]
            vpd_O1TE+=O1_T[i]*O1_E[i]
        O1_Tlength=math.sqrt(O1_Tlength)
        ratio =(vpd_O1TE-math.sqrt(vpd_O1TE**2-O1_Tlength**2*(O1_Elength**2-re**2)))/(O1_Tlength**2)
        ratio1=(vpd_O1TE+math.sqrt(vpd_O1TE**2-O1_Tlength**2*(O1_Elength**2-re**2)))/(O1_Tlength**2)
        dratio=ratio1-ratio
        for k in range(1000):
            for i in range(3):
                O1_O[i]=(ratio+k*dratio/1000)*O1_M[i]
            r=O1_Mlength*ratio*math.sin(thetaM_O1)/math.cos(thetaM_O1)#OT大小
            for n in range(360):#分成360份，一份一度
                #print(n)
                E_Tlength=0
                for i in range(3):
                    O_T[i]=r*(O1_X[i]*math.cos(n*math.pi/180)+O1_Y[i]*math.sin(n*math.pi/180))
                    O1_T[i]= O1_O[i]+O_T[i]
                    E_T[i]=O1_T[i]-O1_E[i]
                    E_Tlength+=E_T[i]*E_T[i]
                E_Tlength=math.sqrt(E_Tlength)
                #print(math.fabs(E_Tlength-re))
                if math.fabs(E_Tlength-re)<1e-8:#大概1.5km
                    #print("哈哈哈")
                    E_TT=novas.cel2ter(jd_ut1, 0.0, delta_t, x_pole, y_pole,E_T, 1, 1,0)
                    coor1=novas.vector2radec(E_TT)
                    if day_time(coor1[0]*15.0,coor1[1],E_O4):
                        lons.append(coor1[0]*15.0)
                        lats.append(coor1[1])
        if len(lons)>0:
            #print(len(lons))
            #t=jdutc2bt(jd_utc-8/24)
            coor.append([lons,lats,jd_utc])
            str1=' '.join(str(i) for i in lons)+"\n"
            f.write(str1)
            str1=' '.join(str(i) for i in lats)+"\n"
            f.write(str1)
            str1=str(jd_utc)+"\n"
            f.write(str1)
f.close()
print(len(coor))
duration=int(len(coor)/20)
fig=plt.figure()
def make_frame(t):
    plt.cla()
    lons=coor[int(t/0.1)][0]
    lats=coor[int(t/0.1)][1]
    jd=coor[int(t/0.1)][2]
    time=jdutc2bt(jd-8/24)
    my_map=Basemap(projection='ortho',lat_0=45,lon_0=130)
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




