import sys
sys.path.insert(0, r'C:\Scripts')
from math import *
from py_post import *
import os
from heapq import nsmallest
from heapq import nlargest
import time 
#
# first version with support for axisymmetric models...
# first python 3.5 script, same as v41 (python 2.7)
#------------------------------------------------------------------------------
# function definitions

def pgrid(nid,x,y,z,tmp_nas):

    tmp='GRID*   %16i                %16.7f%16.7f \n*       %16.7f \n' %(nid,x,y,z)
    tmp_nas.write(tmp)
def ptria(eid,pid,nid1,nid2,nid3,tmp_nas):

    tmp='CTRIA3  %8i%8i%8i%8i%8i\n' %(eid,pid,nid1,nid2,nid3)
    tmp_nas.write(tmp)
    
    
def pbeam(eid,pid,gid1,gid2,tmp_nas):
  
    tmp='CBAR    %8i%8i%8i%8i     1.0     0.0     0.0 \n' %(eid,pid,gid1,gid2)
    tmp_nas.write(tmp)

def Sprod(g1,g2,g3,o1,o2,o3):
    qq=(g1*o1+g2*o2+g3*o3)/(sqrt(g1*g1+g2*g2+g3*g3)*sqrt(o1*o1+o2*o2+o3*o3))
    if qq > 0.999999999999999 :
       angle=0.0
    else :
       angle=acos((g1*o1+g2*o2+g3*o3)/(sqrt(g1*g1+g2*g2+g3*g3)*sqrt(o1*o1+o2*o2+o3*o3)))*180.0/3.141562
    return angle

def SprodT(g1,g2,g3,o1,o2,o3):

    add=0.0
    angle=acos((g1*o1+g2*o2+g3*o3)/(sqrt(g1*g1+g2*g2+g3*g3)*sqrt(o1*o1+o2*o2+o3*o3)))*180.0/3.141562
    if g2 > 0.001 :
        angle=360.0-angle
    
    return angle
    
def SprodS(g1,g2,g3,o1,o2,o3,n1,n2,n3,negflag):

    angle=acos((g1*o1+g2*o2+g3*o3)/(sqrt(g1*g1+g2*g2+g3*g3)*sqrt(o1*o1+o2*o2+o3*o3)))*180.0/3.141562

    xAx=g1
    xAy=g2
    xAz=g3
        
    o1nA=xAy*float(o3)-xAz*float(o2)
    o2nA=xAz*float(o1)-xAx*float(o3)
    o3nA=xAx*float(o2)-xAy*float(o1)
        
    vlen=sqrt(o1nA*o1nA+o2nA*o2nA+o3nA*o3nA)
    o1n=o1nA/vlen
    o2n=o2nA/vlen
    o3n=o3nA/vlen
    
    sign1=(n1*o1n+n2*o2n+n3*o3n)/(sqrt(n1*n1+n2*n2+n3*n3)*sqrt(o1n*o1n+o2n*o2n+o3n*o3n))
    
    sangle=sign1*angle
    
    if negflag==1 and sign1 < 0.0:
        sangle=sangle+360.0
        

    return sangle
    
def get_stress(g1,g2,t,tens):
    
    
    for k in range(0,len(tens)):
        if int(g1)==tens[k].id:
            tens1=tens[k]
        if int(g2)==tens[k].id:
            tens2=tens[k]
    int_tens=[(tens1.t11+t*(tens2.t11-tens1.t11)),(tens1.t12+t*(tens2.t12-tens1.t12)),(tens1.t13+t*(tens2.t13-tens1.t13)),(tens1.t22+t*(tens2.t22-tens1.t22)),(tens1.t23+t*(tens2.t23-tens1.t23)),(tens1.t33+t*(tens2.t33-tens1.t33))]
#    print int_tens
    return int_tens
def get_plst(g1,g2,t,scal):
    
    
    for k in range(0,len(scal)):
        if int(g1)==scal[k].id:
            scal1=scal[k]
        if int(g2)==scal[k].id:
            scal2=scal[k]
    int_scal=scal1.value+t*(scal2.value-scal1.value)
#    print int_tens
    return int_scal    
def get_stressT(g1,g2,t,tens):
    
    
    for k in range(0,len(tens)):
        if int(g1)==tens[k].id:
            tens1=tens[k]
        if int(g2)==tens[k].id:
            tens2=tens[k]
    axial_stress=(tens1.t33+t*(tens2.t33-tens1.t33))
#    print int_tens
    return axial_stress
    
def Xprod(g1,g2,g3,o1,o2,o3):

    xAx=g1
    xAy=g2
    xAz=g3
        
    o1nA=xAy*float(o3)-xAz*float(o2)
    o2nA=xAz*float(o1)-xAx*float(o3)
    o3nA=xAx*float(o2)-xAy*float(o1)
        
    vlen=sqrt(o1nA*o1nA+o2nA*o2nA+o3nA*o3nA)
    o1n=o1nA/vlen
    o2n=o2nA/vlen
    o3n=o3nA/vlen
    
    return o1n,o2n,o3n,vlen        
    
def cross(absBB,orig,nvec,edge):
#k1=[n1,n2,(n2[1]-n1[1]),(n2[2]-n1[2]),(n2[3]-n1[3])]
    nx=nvec[0]  # unit vector in local normal direction
    ny=nvec[1]
    nz=nvec[2]
    xnoll=orig[0]
    ynoll=orig[1]
    znoll=orig[2]
    vx=float(edge[2])
    vy=float(edge[3])
    vz=float(edge[4])
    xnollp=float(edge[0][1])
    ynollp=float(edge[0][2])
    znollp=float(edge[0][3])
    scal=(nx*vx+ny*vy+nz*vz)
    elen=sqrt(vx**2+vy**2+vz**2)
    locX=xnollp-xnoll  # local vector from Origo of plane to node 1
    locY=ynollp-ynoll
    locZ=znollp-znoll
    if scal==0.0:
        t='parallel'
        return t
    else:
        dist=abs(locX*nx+locY*ny+locZ*nz)
        t=-1.0*(nx*(locX)+ny*(locY)+nz*(locZ))/scal
        if t>=0.0 and t<=1.0:
            xcr=xnollp+t*vx
            ycr=ynollp+t*vy
            zcr=znollp+t*vz
            if absBB == 0 :
                dist=t-0.5
            return t,xcr,ycr,zcr,dist
        else:
#            t='nocross'
            if absBB == 1 :
               t=dist
            else :
               t=t-0.5
            return 'n',t,elen

def getForceMom(n1,n2,n3,A,S1,S2,S3,AT):

    centx=(n1[1]+n2[1]+n3[1])/3
    centy=(n1[2]+n2[2]+n3[2])/3
    centz=(n1[3]+n2[3]+n3[3])/3
    
    a11=AT[0]
    a12=AT[1]
    a13=AT[2]
    a21=AT[3]
    a22=AT[4]
    a23=AT[5]
    a31=AT[6]
    a32=AT[7]
    a33=AT[8]
        
    sx1=a11*a11*S1[0]+a12*a12*S1[3]+a13*a13*S1[5]+2*a11*a13*S1[2]+2*a12*a13*S1[4]+2*a11*a12*S1[1]
    sy1=a11*a21*S1[0]+a12*a22*S1[3]+a13*a23*S1[5]+(a11*a22+a12*a21)*S1[1]+(a12*a23+a13*a22)*S1[4]+(a11*a23+a13*a21)*S1[2]
    sz1=a11*a31*S1[0]+a12*a32*S1[3]+a13*a33*S1[5]+(a11*a32+a12*a31)*S1[1]+(a12*a33+a13*a32)*S1[4]+(a11*a33+a13*a31)*S1[2]
    
    sx2=a11*a11*S2[0]+a12*a12*S2[3]+a13*a13*S2[5]+2*a11*a13*S2[2]+2*a12*a13*S2[4]+2*a11*a12*S2[1]
    sy2=a11*a21*S2[0]+a12*a22*S2[3]+a13*a23*S2[5]+(a11*a22+a12*a21)*S2[1]+(a12*a23+a13*a22)*S2[4]+(a11*a23+a13*a21)*S2[2]
    sz2=a11*a31*S2[0]+a12*a32*S2[3]+a13*a33*S2[5]+(a11*a32+a12*a31)*S2[1]+(a12*a33+a13*a32)*S2[4]+(a11*a33+a13*a31)*S2[2]

    sx3=a11*a11*S3[0]+a12*a12*S3[3]+a13*a13*S3[5]+2*a11*a13*S3[2]+2*a12*a13*S3[4]+2*a11*a12*S3[1]
    sy3=a11*a21*S3[0]+a12*a22*S3[3]+a13*a23*S3[5]+(a11*a22+a12*a21)*S3[1]+(a12*a23+a13*a22)*S3[4]+(a11*a23+a13*a21)*S3[2]
    sz3=a11*a31*S3[0]+a12*a32*S3[3]+a13*a33*S3[5]+(a11*a32+a12*a31)*S3[1]+(a12*a33+a13*a32)*S3[4]+(a11*a33+a13*a31)*S3[2]
    
    Xs=(sx1+sx2+sx3)/3
    Ys=(sy1+sy2+sy3)/3
    Zs=(sz1+sz2+sz3)/3
    
    s1=[sx1,sy1,sz1]
    s2=[sx2,sy2,sz2]
    s3=[sx3,sy3,sz3]
    
    Fx=Xs*A
    Fy=Ys*A
    Fz=Zs*A
    Ax=centx*A
    Ay=centy*A
    Az=centz*A
    
    

    return [centx,centy,centz,Xs,Ys,Zs,Fx,Fy,Fz,A,Ax,Ay,Az,[n1[1],n1[2],n1[3]],[n2[1],n2[2],n2[3]],[n3[1],n3[2],n3[3]],s1,s2,s3]

def getPlStA(A1,pl1,pl2,pl3):
    
    APl=A1*(pl1+pl2+pl3)/3
    
    return APl
    
def getForceMomAxi(y1,y2,yc,S1,S2,Sc):

    swap=0
    if y1 > yc :
       A1=pi*(y1**2-yc**2)
       s1_1=(S1[0]+Sc[0])/2.0
       s2_1=(S1[3]+Sc[3])/2.0
       s3_1=(S1[5]+Sc[5])/2.0
       A2=pi*(yc**2-y2**2)
       s1_2=(S2[0]+Sc[0])/2.0
       s2_2=(S2[3]+Sc[3])/2.0
       s3_2=(S2[5]+Sc[5])/2.0
       swap=1
    else:
       A1=pi*(yc**2-y1**2)
       s1_1=(S1[0]+Sc[0])/2.0
       s2_1=(S1[3]+Sc[3])/2.0
       s3_1=(S1[5]+Sc[5])/2.0
       A2=pi*(y2**2-yc**2)
       s1_2=(S2[0]+Sc[0])/2.0
       s2_2=(S2[3]+Sc[3])/2.0
       s3_2=(S2[5]+Sc[5])/2.0
    
    Fx=A1*s1_1+A2*s1_2
    Fy=A1*s2_1+A2*s2_2
    Fz=A1*s3_1+A2*s3_2
    
    Iyy=abs(pi*(y1**4+y2**4)/4)
    Izz=abs(pi*(y1**4+y2**4)/4)
    Ixx=Iyy+Izz
    A=A1+A2

    return [Fx,Fy,Fz,A,Ixx,Iyy,Izz,yc,s1_1,s2_1,s3_1,s1_2,s2_2,s3_2,swap]

    
def refinePatch(p1,p2,p3,s1,s2,s3,A):
    
    p12=[(p1[0]+p2[0])/2,(p1[1]+p2[1])/2,(p1[2]+p2[2])/2]
    p23=[(p3[0]+p2[0])/2,(p3[1]+p2[1])/2,(p3[2]+p2[2])/2]
    p31=[(p1[0]+p3[0])/2,(p1[1]+p3[1])/2,(p1[2]+p3[2])/2]
    s12=[(s1[0]+s2[0])/2,(s1[1]+s2[1])/2,(s1[2]+s2[2])/2]
    s23=[(s3[0]+s2[0])/2,(s3[1]+s2[1])/2,(s3[2]+s2[2])/2]
    s31=[(s1[0]+s3[0])/2,(s1[1]+s3[1])/2,(s1[2]+s3[2])/2]
    
    centx=[(p1[0]+p12[0]+p31[0])/3,(p12[0]+p2[0]+p23[0])/3,(p23[0]+p3[0]+p31[0])/3,(p12[0]+p23[0]+p31[0])/3]
    centy=[(p1[1]+p12[1]+p31[1])/3,(p12[1]+p2[1]+p23[1])/3,(p23[1]+p3[1]+p31[1])/3,(p12[1]+p23[1]+p31[1])/3]
    centz=[(p1[2]+p12[2]+p31[2])/3,(p12[2]+p2[2]+p23[2])/3,(p23[2]+p3[2]+p31[2])/3,(p12[2]+p23[2]+p31[2])/3]

    centSx=[(s1[0]+s12[0]+s31[0])/3,(s12[0]+s2[0]+s23[0])/3,(s23[0]+s3[0]+s31[0])/3,(s12[0]+s23[0]+s31[0])/3]
    centSy=[(s1[1]+s12[1]+s31[1])/3,(s12[1]+s2[1]+s23[1])/3,(s23[1]+s3[1]+s31[1])/3,(s12[1]+s23[1]+s31[1])/3]
    centSz=[(s1[2]+s12[2]+s31[2])/3,(s12[2]+s2[2]+s23[2])/3,(s23[2]+s3[2]+s31[2])/3,(s12[2]+s23[2]+s31[2])/3]

    centFx=[centSx[0]*A/4,centSx[1]*A/4,centSx[2]*A/4,centSx[3]*A/4]
    centFy=[centSy[0]*A/4,centSy[1]*A/4,centSy[2]*A/4,centSy[3]*A/4]
    centFz=[centSz[0]*A/4,centSz[1]*A/4,centSz[2]*A/4,centSz[3]*A/4]
    
    Ax=[centx[0]*A/4,centx[1]*A/4,centx[2]*A/4,centx[3]*A/4]
    Ay=[centy[0]*A/4,centy[1]*A/4,centy[2]*A/4,centy[3]*A/4]
    Az=[centz[0]*A/4,centz[1]*A/4,centz[2]*A/4,centz[3]*A/4]
    
    XpRet=[p1[0],p12[0],p31[0],p12[0],p2[0],p23[0],p23[0],p3[0],p31[0],p12[0],p23[0],p31[0]]
    YpRet=[p1[1],p12[1],p31[1],p12[1],p2[1],p23[1],p23[1],p3[1],p31[1],p12[1],p23[1],p31[1]]
    ZpRet=[p1[2],p12[2],p31[2],p12[2],p2[2],p23[2],p23[2],p3[2],p31[2],p12[2],p23[2],p31[2]]
    XsRet=[s1[0],s12[0],s31[0],s12[0],s2[0],s23[0],s23[0],s3[0],s31[0],s12[0],s23[0],s31[0]]
    YsRet=[s1[1],s12[1],s31[1],s12[1],s2[1],s23[1],s23[1],s3[1],s31[1],s12[1],s23[1],s31[1]]
    ZsRet=[s1[2],s12[2],s31[2],s12[2],s2[2],s23[2],s23[2],s3[2],s31[2],s12[2],s23[2],s31[2]]
    
    patch1=[centx[0],centy[0],centz[0],centSx[0],centSy[0],centSz[0],centFx[0],centFy[0],centFz[0],A/4,Ax[0],Ay[0],Az[0],[p1[0],p1[1],p1[2]],[p12[0],p12[1],p12[2]],[p31[0],p31[1],p31[2]],[s1[0],s1[1],s1[2]],[s12[0],s12[1],s12[2]],[s31[0],s31[1],s31[2]]]
    patch2=[centx[1],centy[1],centz[1],centSx[1],centSy[1],centSz[1],centFx[1],centFy[1],centFz[1],A/4,Ax[1],Ay[1],Az[1],[p12[0],p12[1],p12[2]],[p2[0],p2[1],p2[2]],[p23[0],p23[1],p23[2]],[s12[0],s12[1],s12[2]],[s2[0],s2[1],s2[2]],[s23[0],s23[1],s23[2]]]
    patch3=[centx[2],centy[2],centz[2],centSx[2],centSy[2],centSz[2],centFx[2],centFy[2],centFz[2],A/4,Ax[2],Ay[2],Az[2],[p23[0],p23[1],p23[2]],[p3[0],p3[1],p3[2]],[p31[0],p31[1],p31[2]],[s23[0],s23[1],s23[2]],[s3[0],s3[1],s3[2]],[s31[0],s31[1],s31[2]]]
    patch4=[centx[3],centy[3],centz[3],centSx[3],centSy[3],centSz[3],centFx[3],centFy[3],centFz[3],A/4,Ax[3],Ay[3],Az[3],[p12[0],p12[1],p12[2]],[p23[0],p23[1],p23[2]],[p31[0],p31[1],p31[2]],[s12[0],s12[1],s12[2]],[s23[0],s23[1],s23[2]],[s31[0],s31[1],s31[2]]]
    
    return [patch1,patch2,patch3,patch4]
    
def meanP(po,ps):
    
    xm=0
    ym=0
    zm=0
    ps11=0
    ps12=0
    ps13=0
    ps22=0
    ps23=0
    ps33=0
    
    for yy in range(0,(len(po))):
        xm=xm+po[yy][1]
        ym=ym+po[yy][2]
        zm=zm+po[yy][3]
        ps11=ps11+ps[yy][0]
        ps12=ps12+ps[yy][1]
        ps13=ps13+ps[yy][2]
        ps22=ps22+ps[yy][3]
        ps23=ps23+ps[yy][4]
        ps33=ps33+ps[yy][5]
        
    div=len(po)
    
    return xm/div,ym/div,zm/div,ps11/div,ps12/div,ps13/div,ps22/div,ps23/div,ps33/div
def meanP_pl(ps):
    
    ps11=0
    
    for yy in range(0,(len(ps))):
        ps11=ps11+ps[yy]
        
    div=len(ps)
    
    return ps11/div

def getLocalVec(c1,c2,c3,o1,o2,o3,AT):

    l1=(c1-o1)*AT[0]+(c2-o2)*AT[1]+(c3-o3)*AT[2]
    l2=(c1-o1)*AT[3]+(c2-o2)*AT[4]+(c3-o3)*AT[5]
    l3=(c1-o1)*AT[6]+(c2-o2)*AT[7]+(c3-o3)*AT[8]
    
    return l1,l2,l3

def nodord4(poi,a1,a2,n1,n2,N):

    sign1=(n1[0]*N[0]+n1[1]*N[1]+n1[2]*N[2])
    sign2=(n2[0]*N[0]+n2[1]*N[1]+n2[2]*N[2])
    alpha=[0,sign1*a1,sign2*a2]
    asort=sorted(alpha)
    sortind=[]
    for s in asort:
        indn=1+alpha.index(s)
        sortind.append(indn)
    return sortind
    
def nodord5(poi,a1,a2,a3,n1,n2,n3,N):

    sign1=(n1[0]*N[0]+n1[1]*N[1]+n1[2]*N[2])
    sign2=(n2[0]*N[0]+n2[1]*N[1]+n2[2]*N[2])
    sign3=(n3[0]*N[0]+n3[1]*N[1]+n3[2]*N[2])
    alpha=[0,sign1*a1,sign2*a2,sign3*a3]
    asort=sorted(alpha)
    sortind=[]
    for s in asort:
        indn=1+alpha.index(s)
        sortind.append(indn)
    return sortind
    
def nodord6(poi,a1,a2,a3,a4,n1,n2,n3,n4,N):

    sign1=(n1[0]*N[0]+n1[1]*N[1]+n1[2]*N[2])
    sign2=(n2[0]*N[0]+n2[1]*N[1]+n2[2]*N[2])
    sign3=(n3[0]*N[0]+n3[1]*N[1]+n3[2]*N[2])
    sign4=(n4[0]*N[0]+n4[1]*N[1]+n4[2]*N[2])
    alpha=[0,sign1*a1,sign2*a2,sign3*a3,sign4*a4]
    asort=sorted(alpha)
    sortind=[]
    for s in asort:
        indn=1+alpha.index(s)
        sortind.append(indn)
    return sortind

def pTenForce(nn,gid,sTVec,plVecStr,tmp_ten) :
    kk=0
    while kk < nn :
        tmp='FORCE,10,%8i,,%8.3f,%s\n' %((gid-nn+kk),sTVec[kk],plVecStr)
        tmp_ten.write(tmp)
        kk=kk+1
    return

def pPload(id,FM,gid,tmp_nas) :
    tmp_nas.write('PLOAD2,%i,%8.2e,%8i\n' %(id,-FM[3],gid))
    tmp_nas.write('PLOAD2,%i,%8.2e,%8i\n' %(id*2,-FM[4],gid))
    tmp_nas.write('PLOAD2,%i,%8.2e,%8i\n' %(id*3,-FM[5],gid))
    
    return
    
def pPloadB(id,FM,eid,tmp_nas) :
    if FM[14] == 1 :
       tmp_nas.write('PLOAD1,%i,%8i,FX,FR,0.0,%8.2e,1.0,%8.2e\n' %(id,eid,FM[8],FM[8]))
       tmp_nas.write('PLOAD1,%i,%8i,FY,FR,0.0,%8.2e,1.0,%8.2e\n' %(id*2,eid,FM[9],FM[9]))
       tmp_nas.write('PLOAD1,%i,%8i,FZ,FR,0.0,%8.2e,1.0,%8.2e\n' %(id*3,eid,FM[10],FM[10]))
       tmp_nas.write('PLOAD1,%i,%8i,FX,FR,0.0,%8.2e,1.0,%8.2e\n' %(id,eid+1,FM[11],FM[11]))
       tmp_nas.write('PLOAD1,%i,%8i,FY,FR,0.0,%8.2e,1.0,%8.2e\n' %(id*2,eid+1,FM[12],FM[12]))
       tmp_nas.write('PLOAD1,%i,%8i,FZ,FR,0.0,%8.2e,1.0,%8.2e\n' %(id*3,eid+1,FM[13],FM[13]))
    else :
       tmp_nas.write('PLOAD1,%i,%8i,FX,FR,0.0,%8.2e,1.0,%8.2e\n' %(id,eid+1,FM[8],FM[8]))
       tmp_nas.write('PLOAD1,%i,%8i,FY,FR,0.0,%8.2e,1.0,%8.2e\n' %(id*2,eid+1,FM[9],FM[9]))
       tmp_nas.write('PLOAD1,%i,%8i,FZ,FR,0.0,%8.2e,1.0,%8.2e\n' %(id*3,eid+1,FM[10],FM[10]))
       tmp_nas.write('PLOAD1,%i,%8i,FX,FR,0.0,%8.2e,1.0,%8.2e\n' %(id,eid,FM[11],FM[11]))
       tmp_nas.write('PLOAD1,%i,%8i,FY,FR,0.0,%8.2e,1.0,%8.2e\n' %(id*2,eid,FM[12],FM[12]))
       tmp_nas.write('PLOAD1,%i,%8i,FZ,FR,0.0,%8.2e,1.0,%8.2e\n' %(id*3,eid,FM[13],FM[13]))
    return
    
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def search_elem(var_list):
    datF,t16F,n1g1,n2g1,n3g1,n4g1,nodF,cbpr,si,txt,extp,sym,stressT,elset1x, \
    refL1,refLc,nas,tenCB,boxTol,immRem,calcPl,DiamTol,nskip=var_list
    stim_1=time.clock()
    if int(si)==1:
        sca=0.001
    else:
        sca=1.0
# flag to print results for elements wihtout contact body, set to 1 in order to print, default 0
    prNone=0
# Open .dat file and print parameters to res-file
#   
    genRem=immRem
    readM=0
    if datF != '' :
       inf=open(datF,'r')
       lin_in=inf.readlines()
       inf.close()
    else :
       readM=1
       print('Read model data from .t16 file')
#       if genRem == 1 :
#           boxTol=-1.0
#           print 'Disable search boxes due to global remeshing'
    inPath=os.path.split(t16F)[0]
#
    comb=0
    fname=t16F
    p = post_open(fname)
    if int(extp) == 0:
        p.extrapolation('translate')
    if int(extp) == 1:
        p.extrapolation('average')
    if int(extp) == 2:
        p.extrapolation('linear')
    if int(extp) == 3:
        p.extrapolation('linear')
        comb=1
    try:
        p.moveto(1)
        print('Results file ',fname,' successfully opened')
    except:
        print('Error opening post file:',fname)
    numET = p.element_tensors()
    numES = p.element_scalars()
#    print numES
    numElem = p.elements()
    numNS = p.node_scalars()
    numNodes = p.nodes()
    ninc = p.increments()
    print('Number of nodes: ',numNodes,', Number of elements: ',numElem )  
    #
    nodelist=[]
    nrbody=0
    foundC=0
    dicx={}
    dicC={}
    dicS={}
    dicEl={}
    dicSetEl={}
    CBNR=[]
    CBNM=[]
    if elset1x != '' :
       elset1=elset1x.split(',')
    else :
       elset1=[]
# tolerance to detect shifting wire numbering, 0.5 means half of wire width angle    
    twistTol=0.5
# tolerance find nodes for curvature calculation, 0.5 means half diameter 
    if DiamTol == '' :    
        DiamTol=0.5
    else :
        DiamTol=float(DiamTol)
# tolerance for bounding box
    absBB=0
    if boxTol == '':
        bbTol=2.0
        proc_out=t16F[:-4] + '_boxsets' + '.proc'
        procOut=open(proc_out,'w')
        procOut.write('| Created by Marc Mentat 2013.0.0 (64bit)\n')
        procOut.write('*prog_option compatibility:prog_version:ment2013\n')
        procOut.write('*prog_analysis_class structural\n')
        procOut.write('*set_update off\n')
    else :
        bbTol=float(boxTol)
        absBB=1
        if bbTol != -1.0 :
            proc_out=t16F[:-4] + '_boxsets' + '.proc'
            procOut=open(proc_out,'w')
            procOut.write('| Created by Marc Mentat 2013.0.0 (64bit)\n')
            procOut.write('*prog_option compatibility:prog_version:ment2013\n')
            procOut.write('*prog_analysis_class structural\n')
            procOut.write('*set_update off\n')
           
# tolerance for finding unique crossing points on edges
    dubtol=1.0e-3*sca
    datVer=0
    readD=0
    fset=0
    elsetListL=[]
    elsetList=[]
    if readM==0:
#
        v13CB={}
        v13NrName={}
        maxX=-1.0e10
        minX=1.0e10
        maxY=-1.0e10
        minY=1.0e10
        maxZ=-1.0e10
        minZ=1.0e10
        elsetListL=[]
        for line in lin_in:
            if line[0:1] == '$':
                continue
            spl=line.split('\n')[0].split()
            if len(spl) < 1:
                spl=[' dummy']
                line=' dummy'
            if spl[0] == 'sizing':
                num_el=spl[2]
                num_nod=spl[3]
            if spl[0] == 'version':
                datVer=int(spl[1])
                if datVer > 12 :
                   readM=1
                   readD=1
                   break
            if line[0:1] != ' ':
                start_print=0
            if start_print==11:
                if line[28:29] == '+' or line[28:29] == '-' :
                    x=line[10:28] + 'E' + line[28:30]
                elif line[27:28] == '+' or line[27:28] == '-':
                    x=line[10:27] + 'E' + line[27:30]
                else:
                    x=line[10:30]
                if line[48:49] == '+' or line[48:49] == '-':
                    y=line[30:48] + 'E' + line[48:50]
                elif line[47:48] == '+' or line[47:48] == '-':
                    y=line[30:47] + 'E' + line[47:50]
                else:
                    y=line[30:50]
                if line[68:69] == '+' or line[68:69] == '-':
                    z=line[50:68] + 'E' + line[68:70]
                elif line[67:68] == '+' or line[67:68] == '-':
                    z=line[50:67] + 'E' + line[67:70]
                else:
                    z=line[50:]
                dicx[int(line[0:10])]=[float(x),float(y),float(z)]
                if float(x) > maxX :
                    maxX=float(x)
                if float(x) < minX :
                    minX=float(x)
                if float(y) > maxY :
                    maxY=float(y)
                if float(y) < minY :
                    minY=float(y)
                if float(z) > maxZ :
                    maxZ=float(z)
                if float(z) < minZ :
                    minZ=float(z)
            if start_print==10:
                start_print=11
            if start_print==2:
                dicEl[int(line[0:10])]=line
                if datVer==13 :
                   v13CB[int(line[0:10])]=currCBi
            if start_print==1:
                if datVer==13 :
                   currCBi=int(spl[7])
                start_print=2
    #
            if start_print==27:
                if len(spl) == 1 :
                    if datVer != 13 :
                        dicC[int(spl[0])]=[cbodynr,cbodyname]
                    start_print=21
                elif spl[1] == 'to':
                    elc=int(spl[0])
                    while elc < (int(spl[2])+1):
                        dicC[elc]=[cbodynr,cbodyname]
                        elc=elc+1
                    start_print=21                
                elif len(spl)== 14 and spl[13]== 'c':
                    for celnr in range(0,len(spl)-1):
                        dicC[int(spl[celnr])]=[cbodynr,cbodyname]
                elif len(spl) < 14 and len(spl) > 1 and int(spl[1]) !=0 :
                    for celnr in range(0,len(spl)):
                        dicC[int(spl[celnr])]=[cbodynr,cbodyname]
                        start_print=21
                if datVer == 13 :
                    start_print=21
            if start_print==26:
                start_print=27
            if start_print==25:
                start_print=26
            if start_print==24:
                start_print=25
            if start_print==235:
                start_print=24
            if start_print==236:
                start_print=235
            if start_print==23:
                if ctype==3:
                    start_print=235
                elif ctype > 3:
                    start_print=236
                else :
                    start_print=24
            if start_print==22:
                if int(spl[1]) == 0:
                    start_print=23
                    cbodynr=int(spl[0])
                    nrbody=nrbody+1
                    cbodyname=line[80:-1]
                    CBNR.append(cbodynr)
                    CBNM.append(cbodyname)
                    if datVer == 13 :
                       v13NrName[cbodynr]=cbodyname
                       
                else:
                    start_print=0
            if start_print==21:
                start_print=22
            if start_print==20:
                num_body=int(spl[0])
                start_print=21
    # contact version 2 or 3
            if start_print==123:
                start_print=22
            if start_print==122:
                start_print=123
            if start_print==121:
                start_print=122
            if start_print==120:
                num_body=int(spl[0])
                start_print=121
            if start_print == 30:
                if len(spl)== 14 and spl[13]== 'c':
                    for setnm in range(0,len(spl)-1):
                        dicS[int(spl[setnm])]=setname
                        elsetListL.append(int(spl[setnm]))
                elif len(spl) < 14 :
                    if len(spl) > 1:
                       if spl[1] == 'to':
                         secnt=0
                         startId=int(spl[0])
                         while secnt < (int(spl[2])-int(spl[0])+1):
                            actId=startId+secnt
                            dicS[actId]=setname
                            elsetListL.append(actId)
                            secnt=secnt+1
                       else:
                         for setnm in range(0,len(spl)):
                            dicS[int(spl[setnm])]=setname
                            elsetListL.append(int(spl[setnm]))
                         start_print=0
                    else:
                      for setnm in range(0,len(spl)):
                         dicS[int(spl[setnm])]=setname
                         elsetListL.append(int(spl[setnm]))
                      start_print=0
                dicSetEl[setname]=elsetListL
            if spl[0] == 'connectivity':
                start_print=1
            if spl[0] == 'coordinates':
                start_print=10
            if spl[0] == 'contact':
                if len(spl) == 1:
                    start_print=20
                    foundC=1
                    ctype=0
                elif spl[1] == '2' or spl[1] == '3' or spl[1] == '4' or spl[1] == '5':
                    start_print=120
                    foundC=1
                    ctype=int(spl[1])
            if spl[0] == 'define' and spl[1] == 'element' and spl[2] == 'set':
                setname=spl[3]
                if setname in elset1:
                    elsetListL=[]
                    start_print=30
        if foundC == 0:
            CBNR.append(0)
            CBNM.append('None')
        if tenCB !='' :
            splTen=tenCB.split(',')
            cntC9=0
            splTeni=[]
            while cntC9 < len(splTen) :
                splTeni.append(int(splTen[cntC9]))
                cntC9=cntC9+1
        if datVer==13 :
            for key in v13CB :
                tmpCB=v13CB[key]
                dicC[key]=[tmpCB,v13NrName[tmpCB]]
    if readM==1 :
        maxX=-1.0e10
        minX=1.0e10
        maxY=-1.0e10
        minY=1.0e10
        maxZ=-1.0e10
        minZ=1.0e10
        ki=0
        while ki < numNodes :
            nod22 = p.node(ki)
            n22x=nod22.x
            n22y=nod22.y
            n22z=nod22.z
            dicx[nod22.id]=[n22x,n22y,n22z]
            if n22x > maxX :
                maxX=n22x
            if n22x < minX :
                minX=n22x
            if n22y > maxY :
                maxY=n22y
            if n22y < minY :
                minY=n22y
            if n22z > maxZ :
                maxZ=n22z
            if n22z < minZ :
                minZ=n22z
            ki=ki+1
        ncbody=p.cbodies()
        ki=0
        while ki < ncbody :
            curCB=p.cbody(ki)
            if curCB.type == 0 :
               foundC=1
               cnid=curCB.id
               cnam=curCB.name
               print(cnid,cnam)
               CBNR.append(cnid)
               CBNM.append(cnam)
               kii=0
               numElC=curCB.nelements
               while kii < numElC :
                   dicC[curCB.elements[kii]]=[cnid,cnam]
                   kii=kii+1
            ki=ki+1
#        if foundC == 0:
#            CBNR.append(0)
#            CBNM.append('None')
        if tenCB !='' :
            splTen=tenCB.split(',')
            cntC9=0
            splTeni=[]
            while cntC9 < len(splTen) :
                splTeni.append(int(splTen[cntC9]))
                cntC9=cntC9+1
        ki=0
        NodeCBid={}
        cbNone=0
        while ki < numElem :
           currEl=p.element_id(ki)
           el22 = p.element(ki)
           numNe=el22.len
           eltype=el22.type
           kii=0
           tmp=str(currEl) + ' ' + str(eltype)
           try :
              currElCB=dicC[currEl][0]
           except :
              if cbNone == 0 :
                 CBNR.append(0)
                 CBNM.append('None')
                 cbNone=1
              dicC[currEl]=[0,'None'] 
              currElCB=dicC[currEl][0]              
           while kii < numNe :
# set the contact body ID for each node.
              NodeCBid[el22.items[kii]]=currElCB
              tmp=tmp  + ' ' +  str(el22.items[kii])
              kii=kii+1
           dicEl[currEl]=tmp
           ki=ki+1
        n = p.sets()
        for i in range(0,n):
            s = p.set(i)
#            print ('Name:',s.name, 'Type: ',s.type)
            if s.type == 'element' and s.name in elset1 :
                for j in range(0,s.len):
                    dicS[s.items[j]]=s.name
                    elsetListL.append(s.items[j])
                    elsetList.append(s.items[j])
                dicSetEl[s.name]=elsetListL
                elsetListL=[]
    if readM == 1 and readD == 1 :
        
        for line in lin_in:
            if line[0:1] == '$':
                continue
            spl=line.split('\n')[0].split()
            if len(spl) < 1:
                spl=[' dummy']
                line=' dummy'
            if (start_print==0 and fset==1 ):
                dicSetEl[setname]=elsetListL
                fset=0
            if line[0:1] != ' ':
                start_print=0        
            if start_print == 30:
                if len(spl)== 14 and spl[13]== 'c':
                    for setnm in range(0,len(spl)-1):
                        dicS[int(spl[setnm])]=setname
                        elsetListL.append(int(spl[setnm]))
                        elsetList.append(int(spl[setnm]))
                elif len(spl) < 14 :
                    if len(spl) > 1:
                       if spl[1] == 'to':
                         secnt=0
                         startId=int(spl[0])
                         while secnt < (int(spl[2])-int(spl[0])+1):
                            actId=startId+secnt
                            dicS[actId]=setname
                            elsetListL.append(actId)
                            elsetList.append(actId)
                            secnt=secnt+1
                       else:
                         for setnm in range(0,len(spl)):
                            dicS[int(spl[setnm])]=setname
                            elsetListL.append(int(spl[setnm]))
                            elsetList.append(int(spl[setnm]))
                         start_print=0
                    else:
                      for setnm in range(0,len(spl)):
                         dicS[int(spl[setnm])]=setname
                         elsetListL.append(int(spl[setnm]))
                         elsetList.append(int(spl[setnm]))
                      start_print=0
            if spl[0] == 'define' and spl[1] == 'element' and spl[2] == 'set':
                setname=spl[3]
                if setname in elset1:
                    elsetListL=[]
                    start_print=30
                    fset=1
#---------------------------------------------------------
# Main loop 
    n1g=[]
    n2g=[]
    n3g=[]
    n4g=[]
    axisym=1
    if nodF !='' :
        inG=open(nodF,'r')
        lin_G=inG.readlines()
        inG.close()
        for linG in lin_G:
            spl=linG.split('\n')[0].split(',')
            curv_calc=0
            if  len(spl) == 3:
                n1g.append(int(spl[0]))
                n2g.append(int(spl[1]))
                n3g.append(int(spl[2]))
                curv_calc=0
                axisym=0
            elif  len(spl) == 4:
                n1g.append(int(spl[0]))
                n2g.append(int(spl[1]))
                n3g.append(int(spl[2]))
                n4g.append(int(spl[3]))
                curv_calc=1
                axisym=0
            elif  len(spl) == 1:
                n1g.append(int(spl[0]))
            else:
                print('Number of nodes in section not equal to 3 or 4 or blank line found in .txt file. Line will be skipped!')
    else:
        n1g.append(int(n1g1))
        if n2g1 != '' :
            n2g.append(int(n2g1))
            n3g.append(int(n3g1))
            axisym=0
        if n4g1 != '' :
            n4g.append(int(n4g1))
            curv_calc=1
        else:
            curv_calc=0
    numSec=len(n1g)
    
    
    if genRem == 1 and curv_calc==1:
       curv_calc=0
       print('Curvature calculation disabled due to global remeshing...')
       
    print('Number of sections to process: ' + str(numSec))
    
    curSec=0
    tot_inc=(ninc-1)*numSec
    work_cnt=0

    stim0=time.clock()
    read_time=stim0-stim_1
    
# branch to differ 3D from axisymmetric case...
#    print "axi: ",axisym
    if axisym == 0 :
    
        boxNodDicS={}   
# variables used for book keeping...
        pn1d={}    
        pn2d={}    
        pn3d={}    
        pn4d={}  
        ten_outd={}
        ten_based={}
        nas_outd={}
        txt_outd={}
        csv_outd={}
        txt_outEd={}
        csv_outEd={}
        csv_outSd={}
        curveNodD={}
        while curSec < numSec:
            print('Processing section number ' + str(curSec+1))
            stim=time.clock()
            time_inc=0.0
            time_disp=0.0
            time_el=0.0
            time_ten=0.0
            pn1=n1g[curSec]
            pn2=n2g[curSec]
            pn3=n3g[curSec]
            pn1d[curSec]=pn1
            pn2d[curSec]=pn2
            pn3d[curSec]=pn3
            boxNodDicS[pn1]=1
            boxNodDicS[pn2]=1
            boxNodDicS[pn3]=1
            pn4=''
            pn4d[curSec]=pn4        
            pn1st=str(pn1)
            oval_calc=0
            diamX=sqrt((dicx[n1g[curSec]][0]-dicx[n2g[curSec]][0])**2+(dicx[n1g[curSec]][1]-dicx[n2g[curSec]][1])**2+(dicx[n1g[curSec]][2]-dicx[n2g[curSec]][2])**2)
# code to try to find point 4 if not given
            if curv_calc==0 :
                try :
                    print('Trying to find point 4..')
                    mid1_2=[(dicx[n1g[curSec]][0]+dicx[n2g[curSec]][0])/2.0,(dicx[n1g[curSec]][1]+dicx[n2g[curSec]][1])/2.0,(dicx[n1g[curSec]][2]+dicx[n2g[curSec]][2])/2.0]
                    p3t=[dicx[n3g[curSec]][0]-mid1_2[0],dicx[n3g[curSec]][1]-mid1_2[1],dicx[n3g[curSec]][2]-mid1_2[2]]
                    pt4dist=[]
                    pt4nod=[]
                    for nodX in dicx:
                        d4=sqrt((dicx[nodX][0]-mid1_2[0]+p3t[0])**2 + (dicx[nodX][1]-mid1_2[1]+p3t[1])**2 + (dicx[nodX][2]-mid1_2[2]+p3t[2])**2)
                        if d4 < 50.0*sca :
                            pt4dist.append(d4)
                            pt4nod.append(nodX)
                    sm1_p3=nsmallest(1,pt4dist)
                    pn4=pt4nod[pt4dist.index(sm1_p3[0])]
                    n4g.append(pn4)
                    pn4d[curSec]=pn4        
                    boxNodDicS[pn4]=1
                    diamY=sqrt((dicx[n3g[curSec]][0]-dicx[n4g[curSec]][0])**2+(dicx[n3g[curSec]][1]-dicx[n4g[curSec]][1])**2+(dicx[n3g[curSec]][2]-dicx[n4g[curSec]][2])**2)
                    print('Point 4 found, ID: ',pn4)
                    oval_calc=1
                except :
                    print('Failed to find point 4..')
                    oval_calc=0
            
# calculate curvature if requested by 4th node..
            if curv_calc == 1:
                pn4=n4g[curSec]
                pn4d[curSec]=pn4        
                boxNodDicS[pn4]=1
                diamY=sqrt((dicx[n3g[curSec]][0]-dicx[n4g[curSec]][0])**2+(dicx[n3g[curSec]][1]-dicx[n4g[curSec]][1])**2+(dicx[n3g[curSec]][2]-dicx[n4g[curSec]][2])**2)
# check contact body ID
                pn1CB=NodeCBid[pn1]
                pn2CB=NodeCBid[pn2]
                pn3CB=NodeCBid[pn3]
                pn4CB=NodeCBid[pn4]
    
# find length and nodes to use for curvature calculation
                pipeLx=maxX-minX
                pipeLy=maxY-minY
                pipeLz=maxZ-minZ
                if pipeLx > 2*pipeLz :
                    localX=0
                    dof1=1
                    dof2=2
                    if curSec==0 :
                        print('Pipe axis is X')
                    plVecStr='1.0,0.0,0.0'
                elif pipeLz > pipeLx :
                    localX=2
                    dof1=1
                    dof2=0
                    if curSec==0 :
                        print('pipe axis is Z')
                    plVecStr='0.0,0.0,1.0'
                elif pipeLy > 2*pipeLz :
                    localX=1
                    dof1=0
                    dof2=2
                    if curSec==0 :
                        print('pipe axis is Y')
                    plVecStr='0.0,1.0,0.0'
                else :
                    print('Structure seems not appropriate for curvature calculation...')
                cDistX=(DiamTol*diamX)+0.001*sca
                yzTol=1.0*sca
                n1Dpyz=[]
                n2Dpyz=[]
                n3Dpyz=[]
                n4Dpyz=[]
                n1Npyz=[]
                n2Npyz=[]
                n3Npyz=[]
                n4Npyz=[]
                n1Dnyz=[]
                n2Dnyz=[]
                n3Dnyz=[]
                n4Dnyz=[]
                n1Nnyz=[]
                n2Nnyz=[]
                n3Nnyz=[]
                n4Nnyz=[]
# find nodes for curvature calculation
                for nodX in dicx:
                    n1yzDist=sqrt((dicx[nodX][dof1]-dicx[n1g[curSec]][dof1])**2+(dicx[nodX][dof2]-dicx[n1g[curSec]][dof2])**2)
                    if (n1yzDist < yzTol and NodeCBid[nodX] == pn1CB):
                        n1Dpyz.append(abs((dicx[nodX][localX]-dicx[n1g[curSec]][localX])-cDistX))
                        n1Npyz.append(nodX)
                        n1Dnyz.append(abs((dicx[nodX][localX]-dicx[n1g[curSec]][localX])+cDistX))
                        n1Nnyz.append(nodX)
                    n2yzDist=sqrt((dicx[nodX][dof1]-dicx[n2g[curSec]][dof1])**2+(dicx[nodX][dof2]-dicx[n2g[curSec]][dof2])**2)
                    if (n2yzDist < yzTol and NodeCBid[nodX] == pn2CB):
                        n2Dpyz.append(abs((dicx[nodX][localX]-dicx[n2g[curSec]][localX])-cDistX))
                        n2Npyz.append(nodX)
                        n2Dnyz.append(abs((dicx[nodX][localX]-dicx[n2g[curSec]][localX])+cDistX))
                        n2Nnyz.append(nodX)
                    n3yzDist=sqrt((dicx[nodX][dof1]-dicx[n3g[curSec]][dof1])**2+(dicx[nodX][dof2]-dicx[n3g[curSec]][dof2])**2)
                    if (n3yzDist < yzTol and NodeCBid[nodX] == pn3CB):
                        n3Dpyz.append(abs((dicx[nodX][localX]-dicx[n3g[curSec]][localX])-cDistX))
                        n3Npyz.append(nodX)
                        n3Dnyz.append(abs((dicx[nodX][localX]-dicx[n3g[curSec]][localX])+cDistX))
                        n3Nnyz.append(nodX)
                    n4yzDist=sqrt((dicx[nodX][dof1]-dicx[n4g[curSec]][dof1])**2+(dicx[nodX][dof2]-dicx[n4g[curSec]][dof2])**2)
                    if (n4yzDist < yzTol and NodeCBid[nodX] == pn4CB):
                        n4Dpyz.append(abs((dicx[nodX][localX]-dicx[n4g[curSec]][localX])-cDistX))
                        n4Npyz.append(nodX)
                        n4Dnyz.append(abs((dicx[nodX][localX]-dicx[n4g[curSec]][localX])+cDistX))
                        n4Nnyz.append(nodX)
                sm1=nsmallest(1,n1Dpyz)
                n1y1p=n1Npyz[n1Dpyz.index(sm1[0])]
                sm1=nsmallest(1,n1Dnyz)
                n1y1n=n1Nnyz[n1Dnyz.index(sm1[0])]
                sm1=nsmallest(1,n2Dpyz)
                n2y2p=n2Npyz[n2Dpyz.index(sm1[0])]
                sm1=nsmallest(1,n2Dnyz)
                n2y2n=n2Nnyz[n2Dnyz.index(sm1[0])]
                sm1=nsmallest(1,n3Dpyz)
                n3z1p=n3Npyz[n3Dpyz.index(sm1[0])]
                sm1=nsmallest(1,n3Dnyz)
                n3z1n=n3Nnyz[n3Dnyz.index(sm1[0])]
                sm1=nsmallest(1,n4Dpyz)
                n4z2p=n4Npyz[n4Dpyz.index(sm1[0])]
                sm1=nsmallest(1,n4Dnyz)
                n4z2n=n4Nnyz[n4Dnyz.index(sm1[0])]
                boxNodDicS[n1y1p]=1
                boxNodDicS[n1y1n]=1
                boxNodDicS[n2y2p]=1
                boxNodDicS[n2y2n]=1
                boxNodDicS[n3z1p]=1
                boxNodDicS[n3z1n]=1
                boxNodDicS[n4z2p]=1
                boxNodDicS[n4z2n]=1
                curveNodD[curSec]=[n1y1p,n1y1n,n2y2p,n2y2n,n3z1p,n3z1n,n4z2p,n4z2n]
    
            if tenCB !='' :
                ten_out=t16F[:-4] + '_wire_' + pn1st + '.bdf'
                ten_outd[curSec]=ten_out
                ten_base=t16F[:-4] + '_' + pn1st  + '_wire_cb'
                ten_based[curSec]=ten_base
                tmp_ten=open(ten_out,'w')
            if nas==1 :
                nas_out=t16F[:-4] + '_' + pn1st + '.bdf'
                nas_outd[curSec]=nas_out
                tmp_nas=open(nas_out,'w')
# write stuff to out file(s)
            if cbpr != 'none':
                txt_out=t16F[:-4] + '_' + pn1st + 'cb.res'
                csv_out=t16F[:-4] + '_' + pn1st + 'cb.csv'
                txt_outd[curSec]=txt_out
                csv_outd[curSec]=csv_out
                tmp_txt=open(txt_out,'w')
                tmp_csv=open(csv_out,'w')
            if elset1x != '':
                txt_outE=t16F[:-4] + '_' + pn1st + 'elset' + '.res'
                csv_outE=t16F[:-4] + '_' + pn1st + 'elset' + '_sum.csv'    
                csv_outS=t16F[:-4] + '_' + pn1st + 'elset' + '_st.csv'    
                txt_outEd[curSec]=txt_outE
                csv_outEd[curSec]=csv_outE
                csv_outSd[curSec]=csv_outS
                tmp_txtE=open(txt_outE,'w')
                tmp_csvE=open(csv_outE,'w')
                tmp_csvS=open(csv_outS,'w')
# write to out file for contact body
            if cbpr != 'none':
                tmp_txt.write('Cross Section Postprocessing tool output file\n')
                tmp_txt.write('---------------------------------------------\n')
                tmp_txt.write('Marc input file used       : ' + datF + '\n')
                tmp_txt.write('Marc results file used     : ' + t16F + '\n')
                tmp_txt.write('Node number 1              : ' + str(pn1) + '\n')
                tmp_txt.write('Node number 2              : ' + str(pn2) + '\n')
                tmp_txt.write('Node number 3              : ' + str(pn3) + '\n')
                tmp_txt.write('Node number 4              : ' + str(pn4) + '\n')
                tmp_txt.write('Diameter Y (Node 1 to 2)   : ' + str(diamX) + '\n')
                if oval_calc==1 :
                    tmp_txt.write('Diameter Z (Node 3 to 4)   : ' + str(diamY) + '\n')
                if curv_calc==1 :
                    tmp_txt.write('Diameter Z (Node 3 to 4)   : ' + str(diamY) + '\n')
                    tmp_txt.write('Tolerance for Axial nodes  : ' + str(DiamTol) + ' (' + str(DiamTol*diamX) + 'mm)\n')
                    tmp_txt.write('Nodes used to calculate K1 : ' + str(n1y1n) + ',' + str(n1y1p) + '\n')
                    tmp_txt.write('Nodes used to calculate K2 : ' + str(n2y2n) + ',' + str(n2y2p) + '\n')
                    tmp_txt.write('Nodes used to calculate K3 : ' + str(n3z1n) + ',' + str(n3z1p) + '\n')
                    tmp_txt.write('Nodes used to calculate K4 : ' + str(n4z2n) + ',' + str(n4z2p) + '\n')                
                if int(si) == 1:
                    tmp_txt.write('Units used                 : SI\n')
                else:
                    tmp_txt.write('Units used                 : mm\n')
                if cbpr == '':
                    tmp_txt.write('Contact bodies in summation: All\n')
                else:
                    tmp_txt.write('Contact bodies in summation: ' + cbpr + '\n')
                if int(sym) == 1:
                    tmp_txt.write('Symmetri option used       : Yes\n')
                else:
                    tmp_txt.write('Symmetri option used       : No\n')
                if int(stressT) == 0:
                    stressName='Cauchy Stress'
                    stressName2='Stress'
                    tmp_txt.write('Stress Tensor used         : Cauchy Stress\n')
                else:
                    stressName='Stress'
                    stressName2='Cauchy Stress'
                    tmp_txt.write('Stress Tensor used         : Stress\n')
                tmp_txt.write('Refinement level used      : ' + str(nas) + '\n')
                if int(calcPl) == 1:
                    tmp_txt.write('Plastic Strain calculation : Yes\n')
                else:
                    tmp_txt.write('Plastic Strain calculation : No\n')
# write to out file for element set
            if elset1x != '':
                tmp_txtE.write('Cross Section Postprocessing tool output file\n')
                tmp_txtE.write('---------------------------------------------\n')
                tmp_txtE.write('Marc input file used       : ' + datF + '\n')
                tmp_txtE.write('Marc results file used     : ' + t16F + '\n')
                tmp_txtE.write('Node number 1              : ' + str(pn1) + '\n')
                tmp_txtE.write('Node number 2              : ' + str(pn2) + '\n')
                tmp_txtE.write('Node number 3              : ' + str(pn3) + '\n')
                tmp_txtE.write('Node number 4              : ' + str(pn4) + '\n')
                tmp_txtE.write('Diameter Y (Node 1 to 2)   : ' + str(diamX) + '\n')
                if oval_calc==1 :
                    tmp_txt.write('Diameter Z (Node 3 to 4)   : ' + str(diamY) + '\n')
                if curv_calc==1 :
                    tmp_txtE.write('Diameter Z (Node 3 to 4)   : ' + str(diamY) + '\n')
                    tmp_txtE.write('Tolerance for Axial nodes  : ' + str(DiamTol) + ' (' + str(DiamTol*diamX) + 'mm)\n')
                    tmp_txtE.write('Nodes used to calculate K1 : ' + str(n1y1n) + ',' + str(n1y1p) + '\n')
                    tmp_txtE.write('Nodes used to calculate K2 : ' + str(n2y2n) + ',' + str(n2y2p) + '\n')
                    tmp_txtE.write('Nodes used to calculate K3 : ' + str(n3z1n) + ',' + str(n3z1p) + '\n')
                    tmp_txtE.write('Nodes used to calculate K4 : ' + str(n4z2n) + ',' + str(n4z2p) + '\n')                
                if int(si) == 1:
                    tmp_txtE.write('Units used                 : SI\n')
                else:
                    tmp_txtE.write('Units used                 : mm\n')
                tmp_txtE.write('Element set in summation   : ' + elset1x + '\n')
                tmp_txtE.write('Refinement level used      : ' + str(refL1) + '\n')
                if int(sym) == 1:
                    tmp_txtE.write('Symmetri option used       : Yes\n')
                else:
                    tmp_txtE.write('Symmetri option used       : No\n')
                if int(stressT) == 0:
                    stressName='Cauchy Stress'
                    stressName2='Stress'
                    tmp_txtE.write('Stress Tensor used         : Cauchy Stress\n')
                else:
                    stressName='Stress'
                    stressName2='Cauchy Stress'
                    tmp_txtE.write('Stress Tensor used         : Stress\n')
                if int(calcPl) == 1:
                    tmp_txt.write('Plastic Strain calculation : Yes\n')
                else:
                    tmp_txt.write('Plastic Strain calculation : No\n')
            if cbpr != 'none':
                if comb==1:
                    tmp_txt.write('Extrapolation method used  : combined\n')
                else:
                    tmp_txt.write('Extrapolation method used  : ' + p.extrapolate + '\n')
                tmp_txt.write('Number of elements         : ' + str(numElem) + '\n')
                tmp_txt.write('Number of nodes            : ' + str(numNodes) + '\n')
                tmp_txt.write('Number of element tensors  : ' + str(numET) + '\n')
                tmp_txt.write('Number of node scalars     : ' + str(numNS) + '\n')
                tmp_txt.write('Number of increments       : ' + str(ninc-1) + '\n')
            if elset1x != '':
                tmp_txtE.write('Extrapolation method used  : ' + p.extrapolate + '\n')
                tmp_txtE.write('Number of elements         : ' + str(numElem) + '\n')
                tmp_txtE.write('Number of nodes            : ' + str(numNodes) + '\n')
                tmp_txtE.write('Number of element tensors  : ' + str(numET) + '\n')
                tmp_txtE.write('Number of node scalars     : ' + str(numNS) + '\n')
                tmp_txtE.write('Number of increments       : ' + str(ninc-1) + '\n')
            curSec=curSec+1
    
# create plane and local coordinate system
        ni=1
        gid=1
        initAngd={}
        boxElDic={}
        boxNodDic={}
        ori0x={}
        ori0y={}
        ori0z={}
        AT0={}
        curSec=0    
#        ninc=4
        for ni in range(1, ninc, nskip):
            dicxd={}
            curSec=0
            stim1=time.clock()
            atime=float(p.time)
            p.moveto(ni)
            print('Scanning increment number ',str(ni-1))
#
            numNodesX = p.nodes()
            numElemX = p.elements()
            irem=0
#            print numNodesX,numElemX
            if (numNodesX != numNodes or numElemX != numElem ) :
                print('Remeshing assumed')
                print('New number of nodes: ',numNodesX,', New number of elements: ',numElemX)
                foundC=0
                dicx={}
                dicC={}
                dicS={}
                dicEl={}
                CBNR=[]
                CBNM=[]
                numNodes=numNodesX
                numElem=numElemX
                irem=1
                ki=0
                while ki < numNodes :
                    nod22 = p.node(ki)
                    n22x=nod22.x
                    n22y=nod22.y
                    n22z=nod22.z
                    dicx[nod22.id]=[n22x,n22y,n22z]
                    if immRem == 0 :
                        if n22x > maxX :
                            maxX=n22x
                        if n22x < minX :
                            minX=n22x
                        if n22y > maxY :
                            maxY=n22y
                        if n22y < minY :
                            minY=n22y
                        if n22z > maxZ :
                            maxZ=n22z
                        if n22z < minZ :
                            minZ=n22z
        
                    ki=ki+1
                ki=0
                ncbody=p.cbodies()
                ki=0
                while ki < ncbody :
                    curCB=p.cbody(ki)
                    if curCB.type == 0 :
                       foundC=1
                       cnid=curCB.id
                       cnam=curCB.name
                       CBNR.append(cnid)
                       CBNM.append(cnam)
                       kii=0
                       numElC=curCB.nelements
                       while kii < numElC :
                           dicC[curCB.elements[kii]]=[cnid,cnam]
                           kii=kii+1
                    ki=ki+1
#                if foundC == 0:
#                    CBNR.append(0)
#                    CBNM.append('None')
                while ki < numElem :
                   currEl=p.element_id(ki)
                   el22 = p.element(ki)
                   numNe=el22.len
                   eltype=el22.type
                   kii=0
                   tmp=str(currEl) + ' ' + str(eltype)
                   try :
                      currElCB=dicC[currEl][0]
                   except :
                      if cbNone == 0 :
                         CBNR.append(0)
                         CBNM.append('None')
                         cbNone=1
                      dicC[currEl]=[0,'None'] 
                      currElCB=dicC[currEl][0]              
                   while kii < numNe :
# set the contact body ID for each node.
                      NodeCBid[el22.items[kii]]=currElCB
                      tmp=tmp  + ' ' +  str(el22.items[kii])
                      kii=kii+1
                   dicEl[currEl]=tmp
                   ki=ki+1
        
#            
            stim2=time.clock()
            time_inc=time_inc+stim2-stim1
            atime=float(p.time)
            while curSec < numSec :
#                print curSec
                pn1=pn1d[curSec]
                pn2=pn2d[curSec]
                pn3=pn3d[curSec]
                pn4=pn4d[curSec]
                tmp_txt=open(txt_outd[curSec],'a')
                tmp_csv=open(csv_outd[curSec],'a')  
                if nas==1 :                
                    tmp_nas=open(nas_outd[curSec],'a')
                    if tenCB !='' :
                        tmp_ten=open(ten_outd[curSec],'a')
                if elset1x != '':
                    tmp_txtE=open(txt_outEd[curSec],'a')
                    tmp_csvE=open(csv_outEd[curSec],'a')
                    tmp_csvS=open(csv_outSd[curSec],'a')
                boxNodDicT={}
                boxElDicT={}
                boxElDicA={}
                if ni == 1 :
                    kk= 0
                    while kk < numNodes:
                        NodId=p.node_id(kk)
                        dicxd[NodId]=[0.0,0.0,0.0]
                        kk=kk+1
                else :
                    if bbTol == -1.0 or irem == 1 :
                        kk= 0
                        while kk < numNodes:
                            NodId=p.node_id(kk)
                            xx1,yy1,zz1=p.node_displacement(kk)
                            dicxd[NodId] = [xx1,yy1,zz1]
                            kk=kk+1
                    else:
                        for NodId in boxNodDic[curSec] :
                            kk=p.node_sequence(NodId)
                            xx1,yy1,zz1=p.node_displacement(kk)
                            dicxd[NodId] = [xx1,yy1,zz1]
                        for NodId in boxNodDicS :
                            kk=p.node_sequence(NodId)
                            xx1,yy1,zz1=p.node_displacement(kk)
                            dicxd[NodId] = [xx1,yy1,zz1]
#           
                stim3=time.clock()
                time_disp=time_disp+stim3-stim2
                p1=[]
                p2=[]
                p3=[]
                Origo=[]
                OrigoY=[]
                AT=[]
                p1.append(dicx[pn1][0]+dicxd[pn1][0])
                p1.append(dicx[pn1][1]+dicxd[pn1][1])
                p1.append(dicx[pn1][2]+dicxd[pn1][2])
                p2.append(dicx[pn2][0]+dicxd[pn2][0])
                p2.append(dicx[pn2][1]+dicxd[pn2][1])
                p2.append(dicx[pn2][2]+dicxd[pn2][2])
                p3.append(dicx[pn3][0]+dicxd[pn3][0])
                p3.append(dicx[pn3][1]+dicxd[pn3][1])
                p3.append(dicx[pn3][2]+dicxd[pn3][2])
                if oval_calc == 1 :
                    p4=[]
                    p4.append(dicx[pn4][0]+dicxd[pn4][0])
                    p4.append(dicx[pn4][1]+dicxd[pn4][1])
                    p4.append(dicx[pn4][2]+dicxd[pn4][2])
                if curv_calc == 1:
# new code to calculate curvature from 3 points, normalized around centre points
# factor to control the contribution from out of bending plane displacements to curvature.
# Earlier versions, 1.0, now in version 20, 0.0, meaning no contribution
                    dir3curv=0.0
                    p4=[]
                    p1_1=[]
                    p1_2=[]
                    p2_1=[]
                    p2_2=[]
                    p3_1=[]
                    p3_2=[]
                    p4_1=[]
                    p4_2=[]
                    n1y1p=curveNodD[curSec][0]
                    n1y1n=curveNodD[curSec][1]
                    n2y2p=curveNodD[curSec][2]
                    n2y2n=curveNodD[curSec][3]
                    n3z1p=curveNodD[curSec][4]
                    n3z1n=curveNodD[curSec][5]
                    n4z2p=curveNodD[curSec][6]
                    n4z2n=curveNodD[curSec][7]
                    p4.append(dicx[pn4][0]+dicxd[pn4][0])
                    p4.append(dicx[pn4][1]+dicxd[pn4][1])
                    p4.append(dicx[pn4][2]+dicxd[pn4][2])
    
                    p1_1.append(-1.0*sqrt((dicx[n1y1n][localX]+dicxd[n1y1n][localX]-p1[localX])**2+dir3curv*(dicx[n1y1n][dof2]+dicxd[n1y1n][dof2]-p1[dof2])**2))
                    p1_1.append(dicx[n1y1n][dof1]+dicxd[n1y1n][dof1]-p1[dof1])
                    p1_2.append(sqrt((dicx[n1y1p][localX]+dicxd[n1y1p][localX]-p1[localX])**2+dir3curv*(dicx[n1y1p][dof2]+dicxd[n1y1p][dof2]-p1[dof2])**2))
                    p1_2.append(dicx[n1y1p][dof1]+dicxd[n1y1p][dof1]-p1[dof1])
                    aK1=sqrt(p1_1[0]**2+p1_1[1]**2)
                    bK1=sqrt(p1_2[0]**2+p1_2[1]**2)
                    cK1=sqrt((p1_2[0]-p1_1[0])**2+(p1_2[1]-p1_1[1])**2)
                    if p1_1[1] < 0.0 :
                        signK1=-1.0
                    else:
                        signK1=1.0
                    
                    K1=signK1*(sqrt(abs((aK1+bK1+cK1)*(bK1+cK1-aK1)*(aK1-bK1+cK1)*(aK1+bK1-cK1)))/(aK1*bK1*cK1))*sca*1000.0
                    
                    
                    p2_1.append(-1.0*sqrt((dicx[n2y2n][localX]+dicxd[n2y2n][localX]-p2[localX])**2+dir3curv*(dicx[n2y2n][dof2]+dicxd[n2y2n][dof2]-p2[dof2])**2))
                    p2_1.append(dicx[n2y2n][dof1]+dicxd[n2y2n][dof1]-p2[dof1])
                    p2_2.append(sqrt((dicx[n2y2p][localX]+dicxd[n2y2p][localX]-p2[localX])**2+dir3curv*(dicx[n2y2p][dof2]+dicxd[n2y2p][dof2]-p2[dof2])**2))
                    p2_2.append(dicx[n2y2p][dof1]+dicxd[n2y2p][dof1]-p2[dof1])
                    aK2=sqrt(p2_1[0]**2+p2_1[1]**2)
                    bK2=sqrt(p2_2[0]**2+p2_2[1]**2)
                    cK2=sqrt((p2_2[0]-p2_1[0])**2+(p2_2[1]-p2_1[1])**2)
                    if p2_1[1] < 0.0 :
                        signK2=-1.0
                    else:
                        signK2=1.0
                    K2=signK2*(sqrt(abs((aK2+bK2+cK2)*(bK2+cK2-aK2)*(aK2-bK2+cK2)*(aK2+bK2-cK2)))/(aK2*bK2*cK2))*sca*1000.0
                    K1_2=(K1+K2)/2.0
                    
                    p3_1.append(-1.0*sqrt((dicx[n3z1n][localX]+dicxd[n3z1n][localX]-p3[localX])**2+dir3curv*(dicx[n3z1n][dof1]+dicxd[n3z1n][dof1]-p3[dof1])**2))
                    p3_1.append(dicx[n3z1n][dof2]+dicxd[n3z1n][dof2]-p3[dof2])
                    p3_2.append(sqrt((dicx[n3z1p][localX]+dicxd[n3z1p][localX]-p3[localX])**2+dir3curv*(dicx[n3z1p][dof1]+dicxd[n3z1p][dof1]-p3[dof1])**2))
                    p3_2.append(dicx[n3z1p][dof2]+dicxd[n3z1p][dof2]-p3[dof2])
                    aK3=sqrt(p3_1[0]**2+p3_1[1]**2)
                    bK3=sqrt(p3_2[0]**2+p3_2[1]**2)
                    cK3=sqrt((p3_2[0]-p3_1[0])**2+(p3_2[1]-p3_1[1])**2)
                    if p3_1[1] < 0.0 :
                        signK3=-1.0
                    else:
                        signK3=1.0
                    K3=signK3*(sqrt(abs((aK3+bK3+cK3)*(bK3+cK3-aK3)*(aK3-bK3+cK3)*(aK3+bK3-cK3)))/(aK3*bK3*cK3))*sca*1000.0
                    
                    p4_1.append(-1.0*sqrt((dicx[n4z2n][localX]+dicxd[n4z2n][localX]-p4[localX])**2+dir3curv*(dicx[n4z2n][dof1]+dicxd[n4z2n][dof1]-p4[dof1])**2))
                    p4_1.append(dicx[n4z2n][dof2]+dicxd[n4z2n][dof2]-p4[dof2])
                    p4_2.append(sqrt((dicx[n4z2p][localX]+dicxd[n4z2p][localX]-p4[localX])**2+dir3curv*(dicx[n4z2p][dof1]+dicxd[n4z2p][dof1]-p4[dof1])**2))
                    p4_2.append(dicx[n4z2p][dof2]+dicxd[n4z2p][dof2]-p4[dof2])
                    aK4=sqrt(p4_1[0]**2+p4_1[1]**2)
                    bK4=sqrt(p4_2[0]**2+p4_2[1]**2)
                    cK4=sqrt((p4_2[0]-p4_1[0])**2+(p4_2[1]-p4_1[1])**2)
                    if p4_1[1] < 0.0 :
                        signK4=-1.0
                    else:
                        signK4=1.0
                    K4=signK4*(sqrt(abs((aK4+bK4+cK4)*(bK4+cK4-aK4)*(aK4-bK4+cK4)*(aK4+bK4-cK4)))/(aK4*bK4*cK4))*sca*1000.0
                    K3_4=(K3+K4)/2.0
# calculate ovality
                    Dia1=sqrt((p1[0]-p2[0])**2+ (p1[1]-p2[1])**2+(p1[2]-p2[2])**2)                  
                    Dia2=sqrt((p3[0]-p4[0])**2+ (p3[1]-p4[1])**2+(p3[2]-p4[2])**2)
                    Oval1=(Dia1-Dia2)/(Dia1+Dia2)                    
                    Oval2=(Dia2-Dia1)/(Dia1+Dia2)                    
# scale curvature with 1000 to get it to 1/m for models in mm            
                else:
                    K1=0.0
                    K2=0.0
                    K3=0.0
                    K4=0.0
                    K1_2=0.0
                    K3_4=0.0
                    Oval1=0.0
                    Oval2=0.0
                if oval_calc==1 :
# calculate ovality
                    Dia1=sqrt((p1[0]-p2[0])**2+ (p1[1]-p2[1])**2+(p1[2]-p2[2])**2)                  
                    Dia2=sqrt((p3[0]-p4[0])**2+ (p3[1]-p4[1])**2+(p3[2]-p4[2])**2)
                    Oval1=(Dia1-Dia2)/(Dia1+Dia2)                    
                    Oval2=(Dia2-Dia1)/(Dia1+Dia2)                                    
                if curv_calc==0 :
                    Np=Xprod(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2],p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2])
                else :
                    Np=Xprod(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2],p3[0]-p4[0],p3[1]-p4[1],p3[2]-p4[2])
                offset=0.1*sca
                uy=0
                while uy < 3:
                    p1[uy]=p1[uy]+offset*Np[uy]
                    p2[uy]=p2[uy]+offset*Np[uy]
                    p3[uy]=p3[uy]+offset*Np[uy]
                    uy=uy+1
                Origo.append((p1[0]+p2[0]+p3[0])/3)
                Origo.append((p1[1]+p2[1]+p3[1])/3)
                Origo.append((p1[2]+p2[2]+p3[2])/3)
                OrigoY.append((p1[0]+p2[0])/2)
                OrigoY.append((p1[1]+p2[1])/2)
                OrigoY.append((p1[2]+p2[2])/2)
                YpLength=sqrt((p2[0]-p1[0])*(p2[0]-p1[0])+(p2[1]-p1[1])*(p2[1]-p1[1])+(p2[2]-p1[2])*(p2[2]-p1[2]))
                Yp=[(p2[0]-p1[0])/YpLength,(p2[1]-p1[1])/YpLength,(p2[2]-p1[2])/YpLength]
                Zp=Xprod(Np[0],Np[1],Np[2],Yp[0],Yp[1],Yp[2])
                if cbpr != 'none':
                    tmp_txt.write('---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                    tmp='-----------------------------------------    CROSS SECTION RESULTS FOR INCREMENT %4i     ---------------------------------------------------------------------------------------------------------------\n\n' %((ni-1))
                    tmp_txt.write(tmp)
                    tmp='Analysis Time = %4.3f\n\n' %(atime)
                    tmp_txt.write(tmp)
                    tmp_txt.write('     Cross Section Geometrical Properties \n')
                    tmp_txt.write('----------------------------------------------\n')
                    tmp='Local X axis: %10.7f %10.7f %10.7f\n' %(Np[0],Np[1],Np[2])
                    tmp_txt.write(tmp)    
                    tmp='Local Y axis: %10.7f %10.7f %10.7f\n' %(Yp[0],Yp[1],Yp[2])
                    tmp_txt.write(tmp)
                    tmp='Local Z axis: %10.7f %10.7f %10.7f\n\n' %(Zp[0],Zp[1],Zp[2])
                    tmp_txt.write(tmp)
                    tmp_txt.write('                            Cross Section Center of Gravity \n')
                    tmp_txt.write('------------------------------------------------------------------------------\n')
                    tmp_txt.write('                Global CS       Rel. Incr 0 Local CS    Mid point of N1 and N2\n')
# write data to txt file for element set
                if elset1x != '':
                    tmp_txtE.write('------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                    tmp='----------------------------------------------------    CROSS SECTION RESULTS FOR INCREMENT %4i     -------------------------------------------------------\n\n' %((ni-1))
                    tmp_txtE.write(tmp)
                    tmp='                                                             Analysis Time = %4.3f\n\n' %(atime)
                    tmp_txtE.write(tmp)
                    tmp_txtE.write('                                                     Global Cross Section Geometrical Properties \n')
                    tmp_txtE.write('                                                ---------------------------------------------------\n')
                    tmp='                                                 Local X axis:  %10.7f  %10.7f  %10.7f\n' %(Np[0],Np[1],Np[2])
                    tmp_txtE.write(tmp)    
                    tmp='                                                 Local Y axis:  %10.7f  %10.7f  %10.7f\n' %(Yp[0],Yp[1],Yp[2])
                    tmp_txtE.write(tmp)
                    tmp='                                                 Local Z axis:  %10.7f  %10.7f  %10.7f\n\n' %(Zp[0],Zp[1],Zp[2])
                    tmp_txtE.write(tmp)
                    tmp_txtE.write('----------------------------------------------------    CROSS SECTION RESULTS FOR ELEMENT SET  -------------------------------------------------------------\n\n')
    
# Create transformation coefficients
#
                AT.append(Np[0])
                AT.append(Np[1])
                AT.append(Np[2])
                AT.append(Yp[0])
                AT.append(Yp[1])
                AT.append(Yp[2])
                AT.append(Zp[0])
                AT.append(Zp[1])
                AT.append(Zp[2])
# create edges
#
                resFM=[]
                resFMa=[]
                resFMe=[]
                resFMeSetEl=[]
                resFMec=[]
                resFMeSetElc=[]
                resFMnp=[]
                resFMns=[]
                resEL=[]
                resCB=[]
                resCBPl=[]
                resCBPl2=[]
                resTVec=[]
                resSTen=[]
                resAngT=[]
                resTenEl=[]
                resTenCB=[]
                resTenCp=[]
                resTenCpD={}
                resTVecD={}
                resSTenD={}
                calcTen=0
                stim4=time.clock()
                if ni == 1 or irem == 1 :
                    cntClEl=0
                    if bbTol != -1.0 :
                        procOut.write('*store_elements\n')
                        procOut.write('crset%i_%i\n' %(n1g[curSec],(ni-1)))
                foundStress=0
                foundPlst=0
                resPlst=[]
                resPlstE=[]
                maxPlst=[]
                maxPlstE=[]
                maxPlstAll=0.0
                maxPlstEl=0
                maxPlstNode=0
                for elem_id in dicEl:
                    spl=dicEl[elem_id].split('\n')[0].split()
                    q=0
                    elem_type=int(spl[1])
                    if ni == 1 or bbTol== -1.0 or irem==1:
                        boxElDicA[elem_id]=1
                    else :
                        boxElDicA=boxElDic[curSec]
                    el_ind=p.element_sequence(elem_id)
                    for ii in range(0,numES):
                        s_sca = p.element_scalar_label(ii)
                        if s_sca=='Total Equivalent Plastic Strain' and calcPl==1 :
                            slist = p.element_scalar(el_ind,ii)
                            for ks in range(0,len(slist)):
                                if slist[ks].value > maxPlstAll :
                                    maxPlstAll=slist[ks].value
                                    maxPlstEl=elem_id
                                    maxPlstNode=slist[ks].id
                            foundPlst=1
                    for i in range(0,numET):
                        s_ten = p.element_tensor_label(i)
                        if s_ten == stressName :
                            if comb==1:
                                p.extrapolation('linear')
                                tlist = p.element_tensor(el_ind,i)
                                p.extrapolation('average')
                                tlist2=p.element_tensor(el_ind,i)
                            else:
                                tlist = p.element_tensor(el_ind,i)
                                tlist2 = tlist
                            foundStress=1
                        if tenCB !='' :
                            if s_ten=='Cauchy Stress in Preferred Sys' or s_ten=='Stress in Preferred Sys':
                                calcTen=1
                                old_ext=p.extrapolate
                                p.extrapolation('linear')
                                tlistT = p.element_tensor(el_ind,i)
                                p.extrapolation(old_ext)
                    if foundStress== 0 :
                        print('Requested Stress Tensor not found in t16 file')
                        print('Trying to find alternative tensor ',stressName2)
                        stressName=stressName2
                        for i in range(0,numET):
                            s_ten = p.element_tensor_label(i)
                            if s_ten == stressName :
                                if comb==1:
                                    p.extrapolation('linear')
                                    tlist = p.element_tensor(el_ind,i)
                                    p.extrapolation('average')
                                    tlist2=p.element_tensor(el_ind,i)
                                else:
                                    tlist = p.element_tensor(el_ind,i)
                                    tlist2 = tlist
                                foundStress=1
                    if foundStress== 0 :
                        print('Requested Stress Tensor not found in t16 file, quiting...')
                        return 1
                    crossp=[]
                    crossS=[]
                    crossPlst=[]
                    crossSa=[]
                    crossT=[]
                    crList=[]
                    pr_tenel=0
#                    print boxElDicA
                    if elem_type==157 or elem_type==127 or elem_type==241:
                        elem_type=134
                    if elem_type==84 or elem_type==120 or elem_type==61 or elem_type==57 or elem_type==35 or elem_type==21:
                        elem_type=7
                    if elem_type==202 :
                        elem_type=136
                    if elem_type==157 or elem_type==127 or elem_type==241 or elem_type==130 or elem_type==184:
                        elem_type=134
                    if (elem_type == 7 or elem_type == 117 or elem_type==136 or elem_type==134) and (boxElDicA[elem_id]==1 or bbTol==-1.0 ):
                        if foundC == 0:
                            cbnr=0
                        else:
                            try:
                                cbnr=dicC[elem_id][0]
                            except:
                                print('element not in any contact body',elem_id)
                                boxElDicT[elem_id]=0
                                continue
                        cbnrP=ni*100+cbnr
                        n1=[int(spl[2]),dicx[int(spl[2])][0]+dicxd[int(spl[2])][0],dicx[int(spl[2])][1]+dicxd[int(spl[2])][1],dicx[int(spl[2])][2]+dicxd[int(spl[2])][2]]
                        n2=[int(spl[3]),dicx[int(spl[3])][0]+dicxd[int(spl[3])][0],dicx[int(spl[3])][1]+dicxd[int(spl[3])][1],dicx[int(spl[3])][2]+dicxd[int(spl[3])][2]]
                        n3=[int(spl[4]),dicx[int(spl[4])][0]+dicxd[int(spl[4])][0],dicx[int(spl[4])][1]+dicxd[int(spl[4])][1],dicx[int(spl[4])][2]+dicxd[int(spl[4])][2]]
                        n4=[int(spl[5]),dicx[int(spl[5])][0]+dicxd[int(spl[5])][0],dicx[int(spl[5])][1]+dicxd[int(spl[5])][1],dicx[int(spl[5])][2]+dicxd[int(spl[5])][2]]
                        if (elem_type == 7 or elem_type == 117 or elem_type==136) :
                            n5=[int(spl[6]),dicx[int(spl[6])][0]+dicxd[int(spl[6])][0],dicx[int(spl[6])][1]+dicxd[int(spl[6])][1],dicx[int(spl[6])][2]+dicxd[int(spl[6])][2]]
                            n6=[int(spl[7]),dicx[int(spl[7])][0]+dicxd[int(spl[7])][0],dicx[int(spl[7])][1]+dicxd[int(spl[7])][1],dicx[int(spl[7])][2]+dicxd[int(spl[7])][2]]
                        if elem_type == 7 or elem_type == 117:
                            n7=[int(spl[8]),dicx[int(spl[8])][0]+dicxd[int(spl[8])][0],dicx[int(spl[8])][1]+dicxd[int(spl[8])][1],dicx[int(spl[8])][2]+dicxd[int(spl[8])][2]]
                            n8=[int(spl[9]),dicx[int(spl[9])][0]+dicxd[int(spl[9])][0],dicx[int(spl[9])][1]+dicxd[int(spl[9])][1],dicx[int(spl[9])][2]+dicxd[int(spl[9])][2]]
                        
                        k1=[n1,n2,(n2[1]-n1[1]),(n2[2]-n1[2]),(n2[3]-n1[3])]
                        k2=[n2,n3,(n3[1]-n2[1]),(n3[2]-n2[2]),(n3[3]-n2[3])]
                        if elem_type == 7 or elem_type == 117:
                            k3=[n3,n4,(n4[1]-n3[1]),(n4[2]-n3[2]),(n4[3]-n3[3])]
                            k4=[n4,n1,(n1[1]-n4[1]),(n1[2]-n4[2]),(n1[3]-n4[3])]
                            k5=[n5,n6,(n6[1]-n5[1]),(n6[2]-n5[2]),(n6[3]-n5[3])]
                            k6=[n6,n7,(n7[1]-n6[1]),(n7[2]-n6[2]),(n7[3]-n6[3])]
                            k7=[n7,n8,(n8[1]-n7[1]),(n8[2]-n7[2]),(n8[3]-n7[3])]
                            k8=[n8,n5,(n5[1]-n8[1]),(n5[2]-n8[2]),(n5[3]-n8[3])]
                            k9=[n1,n5,(n5[1]-n1[1]),(n5[2]-n1[2]),(n5[3]-n1[3])]
                            k10=[n2,n6,(n6[1]-n2[1]),(n6[2]-n2[2]),(n6[3]-n2[3])]
                            k11=[n3,n7,(n7[1]-n3[1]),(n7[2]-n3[2]),(n7[3]-n3[3])]
                            k12=[n4,n8,(n8[1]-n4[1]),(n8[2]-n4[2]),(n8[3]-n4[3])]
                        elif elem_type==136:
                            k3=[n3,n1,(n1[1]-n3[1]),(n1[2]-n3[2]),(n1[3]-n3[3])]
                            k4=[n4,n5,(n5[1]-n4[1]),(n5[2]-n4[2]),(n5[3]-n4[3])]
                            k5=[n5,n6,(n6[1]-n5[1]),(n6[2]-n5[2]),(n6[3]-n5[3])]
                            k6=[n6,n4,(n4[1]-n6[1]),(n4[2]-n6[2]),(n4[3]-n6[3])]
                            k7=[n1,n4,(n4[1]-n1[1]),(n4[2]-n1[2]),(n4[3]-n1[3])]
                            k8=[n2,n5,(n5[1]-n2[1]),(n5[2]-n2[2]),(n5[3]-n2[3])]
                            k9=[n3,n6,(n6[1]-n3[1]),(n6[2]-n3[2]),(n6[3]-n3[3])]
                        elif elem_type==134:
                            k3=[n3,n1,(n1[1]-n3[1]),(n1[2]-n3[2]),(n1[3]-n3[3])]
                            k4=[n1,n4,(n4[1]-n1[1]),(n4[2]-n1[2]),(n4[3]-n1[3])]
                            k5=[n2,n4,(n4[1]-n2[1]),(n4[2]-n2[2]),(n4[3]-n2[3])]
                            k6=[n3,n4,(n4[1]-n3[1]),(n4[2]-n3[2]),(n4[3]-n3[3])]
                        cr1=cross(absBB,Origo,Np,k1)
                        if cr1[0] != 'p' :
                            if cr1[0] != 'n':
                                crossp.append(cr1)
                                crList.append(1.0e-6)
                                crossS.append(get_stress(n1[0],n2[0],cr1[0],tlist))
                                crossSa.append(get_stress(n1[0],n2[0],cr1[0],tlist2))
                                if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                    crossT.append(get_stressT(n1[0],n2[0],cr1[0],tlistT))
                                if foundPlst==1:
                                    crossPlst.append(get_plst(n1[0],n2[0],cr1[0],slist))
                            else:
                                crList.append(cr1[1])
                        cr2=cross(absBB,Origo,Np,k2)
                        dub=0
                        if cr2[0] != 'p' :
                            if cr2[0] != 'n':
                                crList.append(1.0e-6)
                                for i in range (0,len(crossp)):
                                    if abs(cr2[1]-crossp[i][1]) < dubtol and abs(cr2[2]-crossp[i][2]) < dubtol and abs(cr2[3]-crossp[i][3]) < dubtol :
                                        dub=1
                                if dub==0:
                                    crossp.append(cr2)
                                    crossS.append(get_stress(n2[0],n3[0],cr2[0],tlist))
                                    crossSa.append(get_stress(n2[0],n3[0],cr2[0],tlist2))
                                    if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                        crossT.append(get_stressT(n2[0],n3[0],cr2[0],tlistT))
                                        pr_tenel=1
                                    if foundPlst==1:
                                        crossPlst.append(get_plst(n2[0],n3[0],cr2[0],slist))
                            else:
                                crList.append(cr2[1])
                        if elem_type == 7 or elem_type == 117:
                            cr3=cross(absBB,Origo,Np,k3)
                            dub=0
                            if cr3[0] != 'p' :
                                if cr3[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr3[1]-crossp[i][1]) < dubtol and abs(cr3[2]-crossp[i][2]) < dubtol and abs(cr3[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr3)
                                        crossS.append(get_stress(n3[0],n4[0],cr3[0],tlist))
                                        crossSa.append(get_stress(n3[0],n4[0],cr3[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n3[0],n4[0],cr3[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n3[0],n4[0],cr3[0],slist))
                                else:
                                    crList.append(cr3[1])
                            cr4=cross(absBB,Origo,Np,k4)
                            dub=0
                            if cr4[0] != 'p' :
                                if cr4[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr4[1]-crossp[i][1]) < dubtol and abs(cr4[2]-crossp[i][2]) < dubtol and abs(cr4[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr4)
                                        crossS.append(get_stress(n4[0],n1[0],cr4[0],tlist))
                                        crossSa.append(get_stress(n4[0],n1[0],cr4[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n4[0],n1[0],cr4[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n4[0],n1[0],cr4[0],slist))
                                else:
                                    crList.append(cr4[1])
                            cr5=cross(absBB,Origo,Np,k5)
                            dub=0
                            if cr5[0] != 'p' :
                                if cr5[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr5[1]-crossp[i][1]) < dubtol and abs(cr5[2]-crossp[i][2]) < dubtol and abs(cr5[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr5)
                                        crossS.append(get_stress(n5[0],n6[0],cr5[0],tlist))
                                        crossSa.append(get_stress(n5[0],n6[0],cr5[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n5[0],n6[0],cr5[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n5[0],n6[0],cr5[0],slist))
                                else:
                                    crList.append(cr5[1])
                            cr6=cross(absBB,Origo,Np,k6)
                            dub=0
                            if cr6[0] != 'p' :
                                if cr6[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr6[1]-crossp[i][1]) < dubtol and abs(cr6[2]-crossp[i][2]) < dubtol and abs(cr6[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr6)
                                        crossS.append(get_stress(n6[0],n7[0],cr6[0],tlist))
                                        crossSa.append(get_stress(n6[0],n7[0],cr6[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n6[0],n7[0],cr6[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n6[0],n7[0],cr6[0],slist))
                                else:
                                    crList.append(cr6[1])
                            cr7=cross(absBB,Origo,Np,k7)
                            dub=0
                            if cr7[0] != 'p' :
                                if cr7[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr7[1]-crossp[i][1]) < dubtol and abs(cr7[2]-crossp[i][2]) < dubtol and abs(cr7[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr7)
                                        crossS.append(get_stress(n7[0],n8[0],cr7[0],tlist))
                                        crossSa.append(get_stress(n7[0],n8[0],cr7[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n7[0],n8[0],cr7[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n7[0],n8[0],cr7[0],slist))
                                else:
                                    crList.append(cr7[1])
                            cr8=cross(absBB,Origo,Np,k8)
                            dub=0
                            if cr8[0] != 'p' :
                                if cr8[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr8[1]-crossp[i][1]) < dubtol and abs(cr8[2]-crossp[i][2]) < dubtol and abs(cr8[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr8)
                                        crossS.append(get_stress(n8[0],n5[0],cr8[0],tlist))
                                        crossSa.append(get_stress(n8[0],n5[0],cr8[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n8[0],n5[0],cr8[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n8[0],n5[0],cr8[0],slist))
                                else:
                                    crList.append(cr8[1])
                            cr9=cross(absBB,Origo,Np,k9)
                            dub=0
                            if cr9[0] != 'p' :
                                if cr9[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr9[1]-crossp[i][1]) < dubtol and abs(cr9[2]-crossp[i][2]) < dubtol and abs(cr9[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr9)
                                        crossS.append(get_stress(n1[0],n5[0],cr9[0],tlist))
                                        crossSa.append(get_stress(n1[0],n5[0],cr9[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n1[0],n5[0],cr9[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n1[0],n5[0],cr9[0],slist))
                                else:
                                    crList.append(cr9[1])
                            cr10=cross(absBB,Origo,Np,k10)
                            dub=0
                            if cr10[0] != 'p' :
                                if cr10[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr10[1]-crossp[i][1]) < dubtol and abs(cr10[2]-crossp[i][2]) < dubtol and abs(cr10[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr10)
                                        crossS.append(get_stress(n2[0],n6[0],cr10[0],tlist))
                                        crossSa.append(get_stress(n2[0],n6[0],cr10[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n2[0],n6[0],cr10[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n2[0],n6[0],cr10[0],slist))
                                else:
                                    crList.append(cr10[1])
                            cr11=cross(absBB,Origo,Np,k11)
                            dub=0
                            if cr11[0] != 'p' : 
                                if cr11[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr11[1]-crossp[i][1]) < dubtol and abs(cr11[2]-crossp[i][2]) < dubtol and abs(cr11[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr11)
                                        crossS.append(get_stress(n3[0],n7[0],cr11[0],tlist))
                                        crossSa.append(get_stress(n3[0],n7[0],cr11[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n3[0],n7[0],cr11[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n3[0],n7[0],cr11[0],slist))
                                else:
                                    crList.append(cr11[1])
                            cr12=cross(absBB,Origo,Np,k12)
                            dub=0
                            if cr12[0] != 'p' :
                                if cr12[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr12[1]-crossp[i][1]) < dubtol and abs(cr12[2]-crossp[i][2]) < dubtol and abs(cr12[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr12)
                                        crossS.append(get_stress(n4[0],n8[0],cr12[0],tlist))
                                        crossSa.append(get_stress(n4[0],n8[0],cr12[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n4[0],n8[0],cr12[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n4[0],n8[0],cr12[0],slist))
                                else:
                                    crList.append(cr12[1])
                        elif elem_type==136:
                            cr3=cross(absBB,Origo,Np,k3)
                            dub=0
                            if cr3[0] != 'p' : 
                                if cr3[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr3[1]-crossp[i][1]) < dubtol and abs(cr3[2]-crossp[i][2]) < dubtol and abs(cr3[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr3)
                                        crossS.append(get_stress(n3[0],n1[0],cr3[0],tlist))
                                        crossSa.append(get_stress(n3[0],n1[0],cr3[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n3[0],n1[0],cr3[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n3[0],n1[0],cr3[0],slist))
                                else:
                                    crList.append(cr3[1])
                            cr4=cross(absBB,Origo,Np,k4)
                            dub=0
                            if cr4[0] != 'p' : 
                                if cr4[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr4[1]-crossp[i][1]) < dubtol and abs(cr4[2]-crossp[i][2]) < dubtol and abs(cr4[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr4)
                                        crossS.append(get_stress(n4[0],n5[0],cr4[0],tlist))
                                        crossSa.append(get_stress(n4[0],n5[0],cr4[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n4[0],n5[0],cr4[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n4[0],n5[0],cr4[0],slist))
                                else:
                                    crList.append(cr4[1])
                            cr5=cross(absBB,Origo,Np,k5)
                            dub=0
                            if cr5[0] != 'p'  : 
                                if cr5[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr5[1]-crossp[i][1]) < dubtol and abs(cr5[2]-crossp[i][2]) < dubtol and abs(cr5[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr5)
                                        crossS.append(get_stress(n5[0],n6[0],cr5[0],tlist))
                                        crossSa.append(get_stress(n5[0],n6[0],cr5[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n5[0],n6[0],cr5[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n5[0],n6[0],cr5[0],slist))
                                else:
                                    crList.append(cr5[1])
                            cr6=cross(absBB,Origo,Np,k6)
                            dub=0
                            if cr6[0] != 'p'  : 
                                if cr6[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr6[1]-crossp[i][1]) < dubtol and abs(cr6[2]-crossp[i][2]) < dubtol and abs(cr6[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr6)
                                        crossS.append(get_stress(n6[0],n4[0],cr6[0],tlist))
                                        crossSa.append(get_stress(n6[0],n4[0],cr6[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n6[0],n4[0],cr6[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n6[0],n4[0],cr6[0],slist))
                                else:
                                    crList.append(cr6[1])
                            cr7=cross(absBB,Origo,Np,k7)
                            dub=0
                            if cr7[0] != 'p'  : 
                                if cr7[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr7[1]-crossp[i][1]) < dubtol and abs(cr7[2]-crossp[i][2]) < dubtol and abs(cr7[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr7)
                                        crossS.append(get_stress(n1[0],n4[0],cr7[0],tlist))
                                        crossSa.append(get_stress(n1[0],n4[0],cr7[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n1[0],n4[0],cr7[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n1[0],n4[0],cr7[0],slist))
                                else:
                                    crList.append(cr7[1])
                            cr8=cross(absBB,Origo,Np,k8)
                            dub=0
                            if cr8[0] != 'p'  : 
                                if cr8[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr8[1]-crossp[i][1]) < dubtol and abs(cr8[2]-crossp[i][2]) < dubtol and abs(cr8[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr8)
                                        crossS.append(get_stress(n2[0],n5[0],cr8[0],tlist))
                                        crossSa.append(get_stress(n2[0],n5[0],cr8[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n2[0],n5[0],cr8[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n2[0],n5[0],cr8[0],slist))
                                else:
                                    crList.append(cr8[1])
                            cr9=cross(absBB,Origo,Np,k9)
                            dub=0
                            if cr9[0] != 'p'  : 
                                if cr9[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr9[1]-crossp[i][1]) < dubtol and abs(cr9[2]-crossp[i][2]) < dubtol and abs(cr9[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr9)
                                        crossS.append(get_stress(n3[0],n6[0],cr9[0],tlist))
                                        crossSa.append(get_stress(n3[0],n6[0],cr9[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n3[0],n6[0],cr9[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n3[0],n6[0],cr9[0],slist))
                                else:
                                    crList.append(cr9[1])
                        elif elem_type==134:
                            cr3=cross(absBB,Origo,Np,k3)
                            dub=0
                            if cr3[0] != 'p' : 
                                if cr3[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr3[1]-crossp[i][1]) < dubtol and abs(cr3[2]-crossp[i][2]) < dubtol and abs(cr3[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr3)
                                        crossS.append(get_stress(n3[0],n1[0],cr3[0],tlist))
                                        crossSa.append(get_stress(n3[0],n1[0],cr3[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n3[0],n1[0],cr3[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n3[0],n1[0],cr3[0],slist))
                                else:
                                    crList.append(cr3[1])
                            cr4=cross(absBB,Origo,Np,k4)
                            dub=0
                            if cr4[0] != 'p' : 
                                if cr4[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr4[1]-crossp[i][1]) < dubtol and abs(cr4[2]-crossp[i][2]) < dubtol and abs(cr4[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr4)
                                        crossS.append(get_stress(n1[0],n4[0],cr4[0],tlist))
                                        crossSa.append(get_stress(n1[0],n4[0],cr4[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n1[0],n4[0],cr4[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n1[0],n4[0],cr4[0],slist))
                                else:
                                    crList.append(cr4[1])
                            cr5=cross(absBB,Origo,Np,k5)
                            dub=0
                            if cr5[0] != 'p'  : 
                                if cr5[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr5[1]-crossp[i][1]) < dubtol and abs(cr5[2]-crossp[i][2]) < dubtol and abs(cr5[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr5)
                                        crossS.append(get_stress(n2[0],n4[0],cr5[0],tlist))
                                        crossSa.append(get_stress(n2[0],n4[0],cr5[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n2[0],n4[0],cr5[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n2[0],n4[0],cr5[0],slist))
                                else:
                                    crList.append(cr5[1])
                            cr6=cross(absBB,Origo,Np,k6)
                            dub=0
                            if cr6[0] != 'p'  : 
                                if cr6[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr6[1]-crossp[i][1]) < dubtol and abs(cr6[2]-crossp[i][2]) < dubtol and abs(cr6[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr6)
                                        crossS.append(get_stress(n3[0],n4[0],cr6[0],tlist))
                                        crossSa.append(get_stress(n3[0],n4[0],cr6[0],tlist2))
                                        if calcTen==1 and (dicC[elem_id][0] in splTeni) :
                                            crossT.append(get_stressT(n3[0],n4[0],cr6[0],tlistT))
                                            pr_tenel=1
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n3[0],n4[0],cr6[0],slist))
                                else:
                                    crList.append(cr6[1])
#
#
#  if min(t) less than a tolerance, save this element...
                        if (ni == 1 or irem==1) and bbTol != -1.0:
                            if min([sqrt(i**2) for i in crList]) < bbTol :
                                boxElDicT[elem_id]=1
                                boxNodDicT[n1[0]]=1
                                boxNodDicT[n2[0]]=1
                                boxNodDicT[n3[0]]=1
                                boxNodDicT[n4[0]]=1
                                cntClEl=cntClEl+1
                                procOut.write(' %i' %elem_id)
                                if cntClEl == 10 :
                                    procOut.write('\n')
                                    cntClEl=0
                                if (elem_type == 7 or elem_type == 117 or elem_type==136) :
                                    boxNodDicT[n5[0]]=1
                                    boxNodDicT[n6[0]]=1
                                if elem_type == 7 or elem_type == 117:
                                    boxNodDicT[n7[0]]=1
                                    boxNodDicT[n8[0]]=1
                            else :
                                boxElDicT[elem_id]=0
# If number of unique crossing points is more than 2, a surface patch is found
                        if len(crossp) > 2:
                            center=meanP(crossp,crossS)
                            center_a=meanP(crossp,crossSa)
                            if foundPlst==1 :
                               center_pl=meanP_pl(crossPlst)
                               maxPlst.append(max(crossPlst))
                               resCBPl2.append(cbnr)
                            cro0=[0,center[0],center[1],center[2]]
                            cro0s=[center[3],center[4],center[5],center[6],center[7],center[8]]
                            cro0sa=[center_a[3],center_a[4],center_a[5],center_a[6],center_a[7],center_a[8]]
                            v1x=crossp[1][1]-crossp[0][1]
                            v1y=crossp[1][2]-crossp[0][2]
                            v1z=crossp[1][3]-crossp[0][3]
                            v2x=crossp[2][1]-crossp[0][1]
                            v2y=crossp[2][2]-crossp[0][2]
                            v2z=crossp[2][3]-crossp[0][3]
                            alpha1=Sprod(v1x,v1y,v1z,v2x,v2y,v2z)
                            n1=Xprod(v1x,v1y,v1z,v2x,v2y,v2z)
# Number of cross 3
                        if len(crossp) == 3:
                            A1=Xprod(crossp[0][1]-center[0],crossp[0][2]-center[1],crossp[0][3]-center[2],crossp[1][1]-center[0],crossp[1][2]-center[1],crossp[1][3]-center[2])[3]/2
                            A2=Xprod(crossp[1][1]-center[0],crossp[1][2]-center[1],crossp[1][3]-center[2],crossp[2][1]-center[0],crossp[2][2]-center[1],crossp[2][3]-center[2])[3]/2
                            A3=Xprod(crossp[2][1]-center[0],crossp[2][2]-center[1],crossp[2][3]-center[2],crossp[0][1]-center[0],crossp[0][2]-center[1],crossp[0][3]-center[2])[3]/2
                            FM1=getForceMom(cro0,crossp[0],crossp[1],A1,cro0s,crossS[0],crossS[1],AT)
                            FM2=getForceMom(cro0,crossp[1],crossp[2],A2,cro0s,crossS[1],crossS[2],AT)
                            FM3=getForceMom(cro0,crossp[2],crossp[0],A3,cro0s,crossS[2],crossS[0],AT)
                            if foundPlst == 1 :
                                pl1=getPlStA(A1,center_pl,crossPlst[0],crossPlst[1])
                                pl2=getPlStA(A2,center_pl,crossPlst[1],crossPlst[2])
                                pl3=getPlStA(A3,center_pl,crossPlst[2],crossPlst[0])
                                resPlst.append(pl1)
                                resPlst.append(pl2)
                                resPlst.append(pl3)
                                cnt_cbap=0
                                while cnt_cbap < 3 :
                                    resCBPl.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
# tension arm cb
                            if pr_tenel==1 :
                                pVec=[crossp[0],crossp[1],crossp[2]]
                                sTVec=[crossT[0],crossT[1],crossT[2]]
                                angleT=SprodS(Zp[0],Zp[1],Zp[2],(center[0]-OrigoY[0]),(center[1]-OrigoY[1]),(center[2]-OrigoY[2]),Np[0],Np[1],Np[2],0)
                                resTVec.append(pVec)
                                resSTen.append(sTVec)
                                resAngT.append(angleT)
                                resTenEl.append(elem_id)
                                resTenCB.append(cbnr)
                                resTenCp.append(len(crossp))
                                resTenCpD[elem_id]=(len(crossp))
                                resTVecD[elem_id]=pVec
                                resSTenD[elem_id]=sTVec
                            if refLc==1:
                                FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                resFM.extend(FM1s[:])
                                resFM.extend(FM2s[:])
                                resFM.extend(FM3s[:])
                                cnt_cbap=0
                                while cnt_cbap < 12 :
                                    resCB.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
                            else:
                                resFM.append(FM1)
                                resFM.append(FM2)
                                resFM.append(FM3)
                                cnt_cbap=0
                                while cnt_cbap < 3 :
                                    resCB.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
                            if comb==1:
                                FM1a=getForceMom(cro0,crossp[0],crossp[1],A1,cro0sa,cro0sa,cro0sa,AT)
                                FM2a=getForceMom(cro0,crossp[1],crossp[2],A2,cro0sa,cro0sa,cro0sa,AT)
                                FM3a=getForceMom(cro0,crossp[2],crossp[0],A3,cro0sa,cro0sa,cro0sa,AT)
                                if refLc==1:
                                    FM1s=refinePatch(FM1a[13],FM1a[14],FM1a[15],FM1a[16],FM1a[17],FM1a[18],A1)
                                    FM2s=refinePatch(FM2a[13],FM2a[14],FM2a[15],FM2a[16],FM2a[17],FM2a[18],A2)
                                    FM3s=refinePatch(FM3a[13],FM3a[14],FM3a[15],FM3a[16],FM3a[17],FM3a[18],A3)
                                    resFMa.extend(FM1s[:])
                                    resFMa.extend(FM2s[:])
                                    resFMa.extend(FM3s[:])
                                else:
                                    resFMa.append(FM1a)
                                    resFMa.append(FM2a)
                                    resFMa.append(FM3a)
                            if elset1x != '' and elem_id in elsetList:
                                resFMnp.append(FM1[14])
                                resFMnp.append(FM2[14])
                                resFMnp.append(FM3[14])
                                resFMns.append(FM1[17])
                                resFMns.append(FM2[17])
                                resFMns.append(FM3[17])
                                resEL.append([elem_id,1])
                                resEL.append([elem_id,2])
                                resEL.append([elem_id,3])
                                if foundPlst == 1 :
                                   resPlstE.append(pl1)
                                   resPlstE.append(pl2)
                                   resPlstE.append(pl3)
                                   maxPlstE.append(max(crossPlst))
                                if comb==1:
                                    resFMec.append(FM1a)
                                    resFMec.append(FM2a)
                                    resFMec.append(FM3a)
                                    d=0
                                    while d < 3 :
                                       resFMeSetElc.append(dicS[elem_id])
                                       d=d+1
                                if refL1==0:
                                    d=0
                                    while d < 3 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    resFMe.append(FM1)
                                    resFMe.append(FM2)
                                    resFMe.append(FM3)
                                if refL1==1:
                                    d=0
                                    while d < 12 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                    FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                    FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                    resFMe.extend(FM1s[:])
                                    resFMe.extend(FM2s[:])
                                    resFMe.extend(FM3s[:])
                                if refL1==2:                        
                                    d=0
                                    while d < 48 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                    FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                    FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                    FM1s2=refinePatch(FM1s[0][13],FM1s[0][14],FM1s[0][15],FM1s[0][16],FM1s[0][17],FM1s[0][18],A1/4)
                                    FM1s3=refinePatch(FM1s[1][13],FM1s[1][14],FM1s[1][15],FM1s[1][16],FM1s[1][17],FM1s[1][18],A1/4)
                                    FM1s4=refinePatch(FM1s[2][13],FM1s[2][14],FM1s[2][15],FM1s[2][16],FM1s[2][17],FM1s[2][18],A1/4)
                                    FM1s5=refinePatch(FM1s[3][13],FM1s[3][14],FM1s[3][15],FM1s[3][16],FM1s[3][17],FM1s[3][18],A1/4)
                                    FM2s2=refinePatch(FM2s[0][13],FM2s[0][14],FM2s[0][15],FM2s[0][16],FM2s[0][17],FM2s[0][18],A2/4)
                                    FM2s3=refinePatch(FM2s[1][13],FM2s[1][14],FM2s[1][15],FM2s[1][16],FM2s[1][17],FM2s[1][18],A2/4)
                                    FM2s4=refinePatch(FM2s[2][13],FM2s[2][14],FM2s[2][15],FM2s[2][16],FM2s[2][17],FM2s[2][18],A2/4)
                                    FM2s5=refinePatch(FM2s[3][13],FM2s[3][14],FM2s[3][15],FM2s[3][16],FM2s[3][17],FM2s[3][18],A2/4)
                                    FM3s2=refinePatch(FM3s[0][13],FM3s[0][14],FM3s[0][15],FM3s[0][16],FM3s[0][17],FM3s[0][18],A3/4)
                                    FM3s3=refinePatch(FM3s[1][13],FM3s[1][14],FM3s[1][15],FM3s[1][16],FM3s[1][17],FM3s[1][18],A3/4)
                                    FM3s4=refinePatch(FM3s[2][13],FM3s[2][14],FM3s[2][15],FM3s[2][16],FM3s[2][17],FM3s[2][18],A3/4)
                                    FM3s5=refinePatch(FM3s[3][13],FM3s[3][14],FM3s[3][15],FM3s[3][16],FM3s[3][17],FM3s[3][18],A3/4)
                                    resFMe.extend(FM1s2[:])
                                    resFMe.extend(FM1s3[:])
                                    resFMe.extend(FM1s4[:])
                                    resFMe.extend(FM1s5[:])
                                    resFMe.extend(FM2s2[:])
                                    resFMe.extend(FM2s3[:])
                                    resFMe.extend(FM2s4[:])
                                    resFMe.extend(FM2s5[:])
                                    resFMe.extend(FM3s2[:])
                                    resFMe.extend(FM3s3[:])
                                    resFMe.extend(FM3s4[:])
                                    resFMe.extend(FM3s5[:])
# Write to Nastran if nas==1
                            if nas == 1 :
                                pgrid(gid,center[0],center[1],center[2],tmp_nas)
                                pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_nas)
                                pgrid(gid+2,crossp[1][1],crossp[1][2],crossp[1][3],tmp_nas)
                                pgrid(gid+3,crossp[2][1],crossp[2][2],crossp[2][3],tmp_nas)
                                gid=gid+4
                                ptria(gid-2,cbnrP,gid-4,gid-3,gid-2,tmp_nas)
                                ptria(gid-1,cbnrP,gid-4,gid-2,gid-1,tmp_nas)
                                ptria(gid,cbnrP,gid-4,gid-1,gid-3,tmp_nas)
#    
                                pPload(10,FM1,gid-2,tmp_nas)
                                pPload(10,FM2,gid-1,tmp_nas)
                                pPload(10,FM3,gid,tmp_nas)
                                if pr_tenel==1:
                                    pgrid(gid,center[0],center[1],center[2],tmp_ten)
                                    pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_ten)
                                    pgrid(gid+2,crossp[1][1],crossp[1][2],crossp[1][3],tmp_ten)
                                    pgrid(gid+3,crossp[2][1],crossp[2][2],crossp[2][3],tmp_ten)
                                    gid=gid+4
                                    ptria(gid-2,cbnrP,gid-4,gid-3,gid-2,tmp_ten)
                                    ptria(gid-1,cbnrP,gid-4,gid-2,gid-1,tmp_ten)
                                    ptria(gid,cbnrP,gid-4,gid-1,gid-3,tmp_ten)
                                    pTenForce(3,gid,sTVec,plVecStr,tmp_ten)
# end 3 cross points
#                            
# 4 cross points
                        if len(crossp) == 4:
                            v3x=crossp[3][1]-crossp[0][1]
                            v3y=crossp[3][2]-crossp[0][2]
                            v3z=crossp[3][3]-crossp[0][3]
                            alpha2=Sprod(v1x,v1y,v1z,v3x,v3y,v3z)
                            n2=Xprod(v1x,v1y,v1z,v3x,v3y,v3z)
                            order=nodord4(crossp,alpha1,alpha2,n1,n2,Np)
                            A1=Xprod(crossp[0][1]-center[0],crossp[0][2]-center[1],crossp[0][3]-center[2],crossp[order[0]][1]-center[0],crossp[order[0]][2]-center[1],crossp[order[0]][3]-center[2])[3]/2
                            A2=Xprod(crossp[order[0]][1]-center[0],crossp[order[0]][2]-center[1],crossp[order[0]][3]-center[2],crossp[order[1]][1]-center[0],crossp[order[1]][2]-center[1],crossp[order[1]][3]-center[2])[3]/2
                            A3=Xprod(crossp[order[1]][1]-center[0],crossp[order[1]][2]-center[1],crossp[order[1]][3]-center[2],crossp[order[2]][1]-center[0],crossp[order[2]][2]-center[1],crossp[order[2]][3]-center[2])[3]/2
                            A4=Xprod(crossp[order[2]][1]-center[0],crossp[order[2]][2]-center[1],crossp[order[2]][3]-center[2],crossp[0][1]-center[0],crossp[0][2]-center[1],crossp[0][3]-center[2])[3]/2
                            FM1=getForceMom(cro0,crossp[0],crossp[order[0]],A1,cro0s,crossS[0],crossS[order[0]],AT)
                            FM2=getForceMom(cro0,crossp[order[0]],crossp[order[1]],A2,cro0s,crossS[order[0]],crossS[order[1]],AT)
                            FM3=getForceMom(cro0,crossp[order[1]],crossp[order[2]],A3,cro0s,crossS[order[1]],crossS[order[2]],AT)
                            FM4=getForceMom(cro0,crossp[order[2]],crossp[0],A4,cro0s,crossS[order[2]],crossS[0],AT)
                            if foundPlst == 1 :
                                pl1=getPlStA(A1,center_pl,crossPlst[0],crossPlst[order[0]])
                                pl2=getPlStA(A2,center_pl,crossPlst[order[0]],crossPlst[order[1]])
                                pl3=getPlStA(A3,center_pl,crossPlst[order[1]],crossPlst[order[2]])
                                pl4=getPlStA(A4,center_pl,crossPlst[order[2]],crossPlst[0])
                                resPlst.append(pl1)
                                resPlst.append(pl2)
                                resPlst.append(pl3)
                                resPlst.append(pl4)
                                cnt_cbap=0
                                while cnt_cbap < 4 :
                                    resCBPl.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
# tension arm cb
                            if pr_tenel==1 :
                                pVec=[crossp[0],crossp[order[0]],crossp[order[1]],crossp[order[2]]]
                                sTVec=[crossT[0],crossT[order[0]],crossT[order[1]],crossT[order[2]]]
                                angleT=SprodS(Zp[0],Zp[1],Zp[2],(center[0]-OrigoY[0]),(center[1]-OrigoY[1]),(center[2]-OrigoY[2]),Np[0],Np[1],Np[2],0)
                                resTVec.append(pVec)
                                resSTen.append(sTVec)
                                resAngT.append(angleT)
                                resTenEl.append(elem_id)
                                resTenCB.append(cbnr)
                                resTenCp.append(len(crossp))
                                resTenCpD[elem_id]=(len(crossp))
                                resTVecD[elem_id]=pVec
                                resSTenD[elem_id]=sTVec
                            if refLc==1:
                                FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                resFM.extend(FM1s[:])
                                resFM.extend(FM2s[:])
                                resFM.extend(FM3s[:])
                                resFM.extend(FM4s[:])
                                cnt_cbap=0
                                while cnt_cbap < 16 :
                                    resCB.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
                            else:
                                resFM.append(FM1)
                                resFM.append(FM2)
                                resFM.append(FM3)
                                resFM.append(FM4)
                                cnt_cbap=0
                                while cnt_cbap < 4 :
                                    resCB.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
                            if comb==1:
                                FM1a=getForceMom(cro0,crossp[0],crossp[order[0]],A1,cro0sa,cro0sa,cro0sa,AT)
                                FM2a=getForceMom(cro0,crossp[order[0]],crossp[order[1]],A2,cro0sa,cro0sa,cro0sa,AT)
                                FM3a=getForceMom(cro0,crossp[order[1]],crossp[order[2]],A3,cro0sa,cro0sa,cro0sa,AT)
                                FM4a=getForceMom(cro0,crossp[order[2]],crossp[0],A4,cro0sa,cro0sa,cro0sa,AT)
                                if refLc==1:
                                    FM1s=refinePatch(FM1a[13],FM1a[14],FM1a[15],FM1a[16],FM1a[17],FM1a[18],A1)
                                    FM2s=refinePatch(FM2a[13],FM2a[14],FM2a[15],FM2a[16],FM2a[17],FM2a[18],A2)
                                    FM3s=refinePatch(FM3a[13],FM3a[14],FM3a[15],FM3a[16],FM3a[17],FM3a[18],A3)
                                    FM4s=refinePatch(FM4a[13],FM4a[14],FM4a[15],FM4a[16],FM4a[17],FM4a[18],A4)
                                    resFMa.extend(FM1s[:])
                                    resFMa.extend(FM2s[:])
                                    resFMa.extend(FM3s[:])
                                    resFMa.extend(FM4s[:])
                                else:
                                    resFMa.append(FM1a)
                                    resFMa.append(FM2a)
                                    resFMa.append(FM3a)
                                    resFMa.append(FM4a)                            
                            if elset1x != '' and elem_id in elsetList:
                                resFMnp.append(FM1[14])
                                resFMnp.append(FM2[14])
                                resFMnp.append(FM3[14])
                                resFMnp.append(FM4[14])
                                resFMns.append(FM1[17])
                                resFMns.append(FM2[17])
                                resFMns.append(FM3[17])
                                resFMns.append(FM4[17])
                                resEL.append([elem_id,1])
                                resEL.append([elem_id,2])
                                resEL.append([elem_id,3])
                                resEL.append([elem_id,4])
                                if foundPlst == 1 :
                                   resPlstE.append(pl1)
                                   resPlstE.append(pl2)
                                   resPlstE.append(pl3)
                                   resPlstE.append(pl4)
                                   maxPlstE.append(max(crossPlst))
                                if comb==1:
                                    resFMec.append(FM1a)
                                    resFMec.append(FM2a)
                                    resFMec.append(FM3a)
                                    resFMec.append(FM4a)
                                    d=0
                                    while d < 4 :
                                       resFMeSetElc.append(dicS[elem_id])
                                       d=d+1
                                if refL1==0:
                                    d=0
                                    while d < 4 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    resFMe.append(FM1)
                                    resFMe.append(FM2)
                                    resFMe.append(FM3)
                                    resFMe.append(FM4)
                                if refL1==1:
                                    d=0
                                    while d < 16 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                    FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                    FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                    FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                    resFMe.extend(FM1s[:])
                                    resFMe.extend(FM2s[:])
                                    resFMe.extend(FM3s[:])
                                    resFMe.extend(FM4s[:])
                                if refL1==2:                        
                                    d=0
                                    while d < 64 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                    FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                    FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                    FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                    FM1s2=refinePatch(FM1s[0][13],FM1s[0][14],FM1s[0][15],FM1s[0][16],FM1s[0][17],FM1s[0][18],A1/4)
                                    FM1s3=refinePatch(FM1s[1][13],FM1s[1][14],FM1s[1][15],FM1s[1][16],FM1s[1][17],FM1s[1][18],A1/4)
                                    FM1s4=refinePatch(FM1s[2][13],FM1s[2][14],FM1s[2][15],FM1s[2][16],FM1s[2][17],FM1s[2][18],A1/4)
                                    FM1s5=refinePatch(FM1s[3][13],FM1s[3][14],FM1s[3][15],FM1s[3][16],FM1s[3][17],FM1s[3][18],A1/4)
                                    FM2s2=refinePatch(FM2s[0][13],FM2s[0][14],FM2s[0][15],FM2s[0][16],FM2s[0][17],FM2s[0][18],A2/4)
                                    FM2s3=refinePatch(FM2s[1][13],FM2s[1][14],FM2s[1][15],FM2s[1][16],FM2s[1][17],FM2s[1][18],A2/4)
                                    FM2s4=refinePatch(FM2s[2][13],FM2s[2][14],FM2s[2][15],FM2s[2][16],FM2s[2][17],FM2s[2][18],A2/4)
                                    FM2s5=refinePatch(FM2s[3][13],FM2s[3][14],FM2s[3][15],FM2s[3][16],FM2s[3][17],FM2s[3][18],A2/4)
                                    FM3s2=refinePatch(FM3s[0][13],FM3s[0][14],FM3s[0][15],FM3s[0][16],FM3s[0][17],FM3s[0][18],A3/4)
                                    FM3s3=refinePatch(FM3s[1][13],FM3s[1][14],FM3s[1][15],FM3s[1][16],FM3s[1][17],FM3s[1][18],A3/4)
                                    FM3s4=refinePatch(FM3s[2][13],FM3s[2][14],FM3s[2][15],FM3s[2][16],FM3s[2][17],FM3s[2][18],A3/4)
                                    FM3s5=refinePatch(FM3s[3][13],FM3s[3][14],FM3s[3][15],FM3s[3][16],FM3s[3][17],FM3s[3][18],A3/4)
                                    FM4s2=refinePatch(FM4s[0][13],FM4s[0][14],FM4s[0][15],FM4s[0][16],FM4s[0][17],FM4s[0][18],A4/4)
                                    FM4s3=refinePatch(FM4s[1][13],FM4s[1][14],FM4s[1][15],FM4s[1][16],FM4s[1][17],FM4s[1][18],A4/4)
                                    FM4s4=refinePatch(FM4s[2][13],FM4s[2][14],FM4s[2][15],FM4s[2][16],FM4s[2][17],FM4s[2][18],A4/4)
                                    FM4s5=refinePatch(FM4s[3][13],FM4s[3][14],FM4s[3][15],FM4s[3][16],FM4s[3][17],FM4s[3][18],A4/4)
                                    resFMe.extend(FM1s2[:])
                                    resFMe.extend(FM1s3[:])
                                    resFMe.extend(FM1s4[:])
                                    resFMe.extend(FM1s5[:])
                                    resFMe.extend(FM2s2[:])
                                    resFMe.extend(FM2s3[:])
                                    resFMe.extend(FM2s4[:])
                                    resFMe.extend(FM2s5[:])
                                    resFMe.extend(FM3s2[:])
                                    resFMe.extend(FM3s3[:])
                                    resFMe.extend(FM3s4[:])
                                    resFMe.extend(FM3s5[:])
                                    resFMe.extend(FM4s2[:])
                                    resFMe.extend(FM4s3[:])
                                    resFMe.extend(FM4s4[:])
                                    resFMe.extend(FM4s5[:])
# Write to Nastran
                            if nas == 1 :
                                pgrid(gid,center[0],center[1],center[2],tmp_nas)
                                pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_nas)
                                pgrid(gid+2,crossp[order[0]][1],crossp[order[0]][2],crossp[order[0]][3],tmp_nas)
                                pgrid(gid+3,crossp[order[1]][1],crossp[order[1]][2],crossp[order[1]][3],tmp_nas)
                                pgrid(gid+4,crossp[order[2]][1],crossp[order[2]][2],crossp[order[2]][3],tmp_nas)
                                gid=gid+5
                                ptria(gid-3,cbnrP,gid-5,gid-4,gid-3,tmp_nas)
                                ptria(gid-2,cbnrP,gid-5,gid-3,gid-2,tmp_nas)
                                ptria(gid-1,cbnrP,gid-5,gid-2,gid-1,tmp_nas)
                                ptria(gid,cbnrP,gid-5,gid-1,gid-4,tmp_nas)
                                pPload(10,FM1,gid-3,tmp_nas)
                                pPload(10,FM2,gid-2,tmp_nas)
                                pPload(10,FM3,gid-1,tmp_nas)
                                pPload(10,FM4,gid,tmp_nas)
                                if pr_tenel==1:
                                    pgrid(gid,center[0],center[1],center[2],tmp_ten)
                                    pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_ten)
                                    pgrid(gid+2,crossp[order[0]][1],crossp[order[0]][2],crossp[order[0]][3],tmp_ten)
                                    pgrid(gid+3,crossp[order[1]][1],crossp[order[1]][2],crossp[order[1]][3],tmp_ten)
                                    pgrid(gid+4,crossp[order[2]][1],crossp[order[2]][2],crossp[order[2]][3],tmp_ten)
                                    gid=gid+5
                                    ptria(gid-3,cbnrP,gid-5,gid-4,gid-3,tmp_ten)
                                    ptria(gid-2,cbnrP,gid-5,gid-3,gid-2,tmp_ten)
                                    ptria(gid-1,cbnrP,gid-5,gid-2,gid-1,tmp_ten)
                                    ptria(gid,cbnrP,gid-5,gid-1,gid-4,tmp_ten)
                                    pTenForce(4,gid,sTVec,plVecStr,tmp_ten)
# end 4 cross points
#                            
# 5 cross points
                        if len (crossp) ==5:
                            v3x=crossp[3][1]-crossp[0][1]
                            v3y=crossp[3][2]-crossp[0][2]
                            v3z=crossp[3][3]-crossp[0][3]
                            v4x=crossp[4][1]-crossp[0][1]
                            v4y=crossp[4][2]-crossp[0][2]
                            v4z=crossp[4][3]-crossp[0][3]
                            alpha2=Sprod(v1x,v1y,v1z,v3x,v3y,v3z)
                            alpha3=Sprod(v1x,v1y,v1z,v4x,v4y,v4z)
                            n2=Xprod(v1x,v1y,v1z,v3x,v3y,v3z)
                            n3=Xprod(v1x,v1y,v1z,v4x,v4y,v4z)
                            order=nodord5(crossp,alpha1,alpha2,alpha3,n1,n2,n3,Np)
                            A1=Xprod(crossp[0][1]-center[0],crossp[0][2]-center[1],crossp[0][3]-center[2],crossp[order[0]][1]-center[0],crossp[order[0]][2]-center[1],crossp[order[0]][3]-center[2])[3]/2
                            A2=Xprod(crossp[order[0]][1]-center[0],crossp[order[0]][2]-center[1],crossp[order[0]][3]-center[2],crossp[order[1]][1]-center[0],crossp[order[1]][2]-center[1],crossp[order[1]][3]-center[2])[3]/2
                            A3=Xprod(crossp[order[1]][1]-center[0],crossp[order[1]][2]-center[1],crossp[order[1]][3]-center[2],crossp[order[2]][1]-center[0],crossp[order[2]][2]-center[1],crossp[order[2]][3]-center[2])[3]/2
                            A4=Xprod(crossp[order[2]][1]-center[0],crossp[order[2]][2]-center[1],crossp[order[2]][3]-center[2],crossp[order[3]][1]-center[0],crossp[order[3]][2]-center[1],crossp[order[3]][3]-center[2])[3]/2
                            A5=Xprod(crossp[order[3]][1]-center[0],crossp[order[3]][2]-center[1],crossp[order[3]][3]-center[2],crossp[0][1]-center[0],crossp[0][2]-center[1],crossp[0][3]-center[2])[3]/2
                            FM1=getForceMom(cro0,crossp[0],crossp[order[0]],A1,cro0s,crossS[0],crossS[order[0]],AT)
                            FM2=getForceMom(cro0,crossp[order[0]],crossp[order[1]],A2,cro0s,crossS[order[0]],crossS[order[1]],AT)
                            FM3=getForceMom(cro0,crossp[order[1]],crossp[order[2]],A3,cro0s,crossS[order[1]],crossS[order[2]],AT)
                            FM4=getForceMom(cro0,crossp[order[2]],crossp[order[3]],A4,cro0s,crossS[order[2]],crossS[order[3]],AT)
                            FM5=getForceMom(cro0,crossp[order[3]],crossp[0],A5,cro0s,crossS[order[3]],crossS[0],AT)
                            if foundPlst == 1 :
                                pl1=getPlStA(A1,center_pl,crossPlst[0],crossPlst[order[0]])
                                pl2=getPlStA(A2,center_pl,crossPlst[order[0]],crossPlst[order[1]])
                                pl3=getPlStA(A3,center_pl,crossPlst[order[1]],crossPlst[order[2]])
                                pl4=getPlStA(A4,center_pl,crossPlst[order[2]],crossPlst[order[3]])
                                pl5=getPlStA(A5,center_pl,crossPlst[order[3]],crossPlst[0])
                                resPlst.append(pl1)
                                resPlst.append(pl2)
                                resPlst.append(pl3)
                                resPlst.append(pl4)
                                resPlst.append(pl5)
                                cnt_cbap=0
                                while cnt_cbap < 5 :
                                    resCBPl.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
# tension arm cb
                            if pr_tenel==1 :
                                pVec=[crossp[0],crossp[order[0]],crossp[order[1]],crossp[order[2]],crossp[order[3]]]
                                sTVec=[crossT[0],crossT[order[0]],crossT[order[1]],crossT[order[2]],crossT[order[3]]]
                                angleT=SprodS(Zp[0],Zp[1],Zp[2],(center[0]-OrigoY[0]),(center[1]-OrigoY[1]),(center[2]-OrigoY[2]),Np[0],Np[1],Np[2],0)
                                resTVec.append(pVec)
                                resSTen.append(sTVec)
                                resAngT.append(angleT)
                                resTenEl.append(elem_id)
                                resTenCB.append(cbnr)
                                resTenCp.append(len(crossp))
                                resTenCpD[elem_id]=(len(crossp))
                                resTVecD[elem_id]=pVec
                                resSTenD[elem_id]=sTVec
                            if refLc==1:
                                FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                FM5s=refinePatch(FM5[13],FM5[14],FM5[15],FM5[16],FM5[17],FM5[18],A5)
                                resFM.extend(FM1s[:])
                                resFM.extend(FM2s[:])
                                resFM.extend(FM3s[:])
                                resFM.extend(FM4s[:])
                                resFM.extend(FM5s[:])
                                cnt_cbap=0
                                while cnt_cbap < 20 :
                                    resCB.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
                            else:
                                resFM.append(FM1)
                                resFM.append(FM2)
                                resFM.append(FM3)
                                resFM.append(FM4)
                                resFM.append(FM5)
                                cnt_cbap=0
                                while cnt_cbap < 5 :
                                    resCB.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
                            if comb==1:
                                FM1a=getForceMom(cro0,crossp[0],crossp[order[0]],A1,cro0sa,cro0sa,cro0sa,AT)
                                FM2a=getForceMom(cro0,crossp[order[0]],crossp[order[1]],A2,cro0sa,cro0sa,cro0sa,AT)
                                FM3a=getForceMom(cro0,crossp[order[1]],crossp[order[2]],A3,cro0sa,cro0sa,cro0sa,AT)
                                FM4a=getForceMom(cro0,crossp[order[2]],crossp[order[3]],A4,cro0sa,cro0sa,cro0sa,AT)
                                FM5a=getForceMom(cro0,crossp[order[3]],crossp[0],A5,cro0sa,cro0sa,cro0sa,AT)
                                if refLc==1:
                                    FM1s=refinePatch(FM1a[13],FM1a[14],FM1a[15],FM1a[16],FM1a[17],FM1a[18],A1)
                                    FM2s=refinePatch(FM2a[13],FM2a[14],FM2a[15],FM2a[16],FM2a[17],FM2a[18],A2)
                                    FM3s=refinePatch(FM3a[13],FM3a[14],FM3a[15],FM3a[16],FM3a[17],FM3a[18],A3)
                                    FM4s=refinePatch(FM4a[13],FM4a[14],FM4a[15],FM4a[16],FM4a[17],FM4a[18],A4)
                                    FM5s=refinePatch(FM5a[13],FM5a[14],FM5a[15],FM5a[16],FM5a[17],FM5a[18],A5)
                                    resFMa.extend(FM1s[:])
                                    resFMa.extend(FM2s[:])
                                    resFMa.extend(FM3s[:])
                                    resFMa.extend(FM4s[:])
                                    resFMa.extend(FM5s[:])
                                else:
                                    resFMa.append(FM1a)
                                    resFMa.append(FM2a)
                                    resFMa.append(FM3a)
                                    resFMa.append(FM4a)                            
                                    resFMa.append(FM5a)
                            if elset1x != '' and elem_id in elsetList:
                                resFMnp.append(FM1[14])
                                resFMnp.append(FM2[14])
                                resFMnp.append(FM3[14])
                                resFMnp.append(FM4[14])
                                resFMnp.append(FM5[14])
                                resFMns.append(FM1[17])
                                resFMns.append(FM2[17])
                                resFMns.append(FM3[17])
                                resFMns.append(FM4[17])
                                resFMns.append(FM5[17])
                                resEL.append([elem_id,1])
                                resEL.append([elem_id,2])
                                resEL.append([elem_id,3])
                                resEL.append([elem_id,4])
                                resEL.append([elem_id,5])
                                if foundPlst == 1 :
                                   resPlstE.append(pl1)
                                   resPlstE.append(pl2)
                                   resPlstE.append(pl3)
                                   resPlstE.append(pl4)
                                   resPlstE.append(pl5)
                                   maxPlstE.append(max(crossPlst))
                                if comb==1:
                                    resFMec.append(FM1a)
                                    resFMec.append(FM2a)
                                    resFMec.append(FM3a)
                                    resFMec.append(FM4a)
                                    resFMec.append(FM5a)
                                    d=0
                                    while d < 5 :
                                       resFMeSetElc.append(dicS[elem_id])
                                       d=d+1
                                if refL1==0:
                                    d=0
                                    while d < 5 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    resFMe.append(FM1)
                                    resFMe.append(FM2)
                                    resFMe.append(FM3)
                                    resFMe.append(FM4)
                                    resFMe.append(FM5)
                                if refL1==1:
                                    d=0
                                    while d < 20 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                    FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                    FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                    FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                    FM5s=refinePatch(FM5[13],FM5[14],FM5[15],FM5[16],FM5[17],FM5[18],A5)
                                    resFMe.extend(FM1s[:])
                                    resFMe.extend(FM2s[:])
                                    resFMe.extend(FM3s[:])
                                    resFMe.extend(FM4s[:])
                                    resFMe.extend(FM5s[:])
                                if refL1==2:                        
                                    d=0
                                    while d < 80 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                    FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                    FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                    FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                    FM5s=refinePatch(FM5[13],FM5[14],FM5[15],FM5[16],FM5[17],FM5[18],A5)
                                    FM1s2=refinePatch(FM1s[0][13],FM1s[0][14],FM1s[0][15],FM1s[0][16],FM1s[0][17],FM1s[0][18],A1/4)
                                    FM1s3=refinePatch(FM1s[1][13],FM1s[1][14],FM1s[1][15],FM1s[1][16],FM1s[1][17],FM1s[1][18],A1/4)
                                    FM1s4=refinePatch(FM1s[2][13],FM1s[2][14],FM1s[2][15],FM1s[2][16],FM1s[2][17],FM1s[2][18],A1/4)
                                    FM1s5=refinePatch(FM1s[3][13],FM1s[3][14],FM1s[3][15],FM1s[3][16],FM1s[3][17],FM1s[3][18],A1/4)
                                    FM2s2=refinePatch(FM2s[0][13],FM2s[0][14],FM2s[0][15],FM2s[0][16],FM2s[0][17],FM2s[0][18],A2/4)
                                    FM2s3=refinePatch(FM2s[1][13],FM2s[1][14],FM2s[1][15],FM2s[1][16],FM2s[1][17],FM2s[1][18],A2/4)
                                    FM2s4=refinePatch(FM2s[2][13],FM2s[2][14],FM2s[2][15],FM2s[2][16],FM2s[2][17],FM2s[2][18],A2/4)
                                    FM2s5=refinePatch(FM2s[3][13],FM2s[3][14],FM2s[3][15],FM2s[3][16],FM2s[3][17],FM2s[3][18],A2/4)
                                    FM3s2=refinePatch(FM3s[0][13],FM3s[0][14],FM3s[0][15],FM3s[0][16],FM3s[0][17],FM3s[0][18],A3/4)
                                    FM3s3=refinePatch(FM3s[1][13],FM3s[1][14],FM3s[1][15],FM3s[1][16],FM3s[1][17],FM3s[1][18],A3/4)
                                    FM3s4=refinePatch(FM3s[2][13],FM3s[2][14],FM3s[2][15],FM3s[2][16],FM3s[2][17],FM3s[2][18],A3/4)
                                    FM3s5=refinePatch(FM3s[3][13],FM3s[3][14],FM3s[3][15],FM3s[3][16],FM3s[3][17],FM3s[3][18],A3/4)
                                    FM4s2=refinePatch(FM4s[0][13],FM4s[0][14],FM4s[0][15],FM4s[0][16],FM4s[0][17],FM4s[0][18],A4/4)
                                    FM4s3=refinePatch(FM4s[1][13],FM4s[1][14],FM4s[1][15],FM4s[1][16],FM4s[1][17],FM4s[1][18],A4/4)
                                    FM4s4=refinePatch(FM4s[2][13],FM4s[2][14],FM4s[2][15],FM4s[2][16],FM4s[2][17],FM4s[2][18],A4/4)
                                    FM4s5=refinePatch(FM4s[3][13],FM4s[3][14],FM4s[3][15],FM4s[3][16],FM4s[3][17],FM4s[3][18],A4/4)
                                    FM5s2=refinePatch(FM5s[0][13],FM5s[0][14],FM5s[0][15],FM5s[0][16],FM5s[0][17],FM5s[0][18],A5/4)
                                    FM5s3=refinePatch(FM5s[1][13],FM5s[1][14],FM5s[1][15],FM5s[1][16],FM5s[1][17],FM5s[1][18],A5/4)
                                    FM5s4=refinePatch(FM5s[2][13],FM5s[2][14],FM5s[2][15],FM5s[2][16],FM5s[2][17],FM5s[2][18],A5/4)
                                    FM5s5=refinePatch(FM5s[3][13],FM5s[3][14],FM5s[3][15],FM5s[3][16],FM5s[3][17],FM5s[3][18],A5/4)
                                    resFMe.extend(FM1s2[:])
                                    resFMe.extend(FM1s3[:])
                                    resFMe.extend(FM1s4[:])
                                    resFMe.extend(FM1s5[:])
                                    resFMe.extend(FM2s2[:])
                                    resFMe.extend(FM2s3[:])
                                    resFMe.extend(FM2s4[:])
                                    resFMe.extend(FM2s5[:])
                                    resFMe.extend(FM3s2[:])
                                    resFMe.extend(FM3s3[:])
                                    resFMe.extend(FM3s4[:])
                                    resFMe.extend(FM3s5[:])
                                    resFMe.extend(FM4s2[:])
                                    resFMe.extend(FM4s3[:])
                                    resFMe.extend(FM4s4[:])
                                    resFMe.extend(FM4s5[:])
                                    resFMe.extend(FM5s2[:])
                                    resFMe.extend(FM5s3[:])
                                    resFMe.extend(FM5s4[:])
                                    resFMe.extend(FM5s5[:])
# Write to Nastran                               
                            if nas == 1:
                                pgrid(gid,center[0],center[1],center[2],tmp_nas)
                                pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_nas)
                                pgrid(gid+2,crossp[order[0]][1],crossp[order[0]][2],crossp[order[0]][3],tmp_nas)
                                pgrid(gid+3,crossp[order[1]][1],crossp[order[1]][2],crossp[order[1]][3],tmp_nas)
                                pgrid(gid+4,crossp[order[2]][1],crossp[order[2]][2],crossp[order[2]][3],tmp_nas)
                                pgrid(gid+5,crossp[order[3]][1],crossp[order[3]][2],crossp[order[3]][3],tmp_nas)
                                gid=gid+6
                                ptria(gid-4,cbnrP,gid-6,gid-5,gid-4,tmp_nas)
                                ptria(gid-3,cbnrP,gid-6,gid-4,gid-3,tmp_nas)
                                ptria(gid-2,cbnrP,gid-6,gid-3,gid-2,tmp_nas)
                                ptria(gid-1,cbnrP,gid-6,gid-2,gid-1,tmp_nas)
                                ptria(gid,cbnrP,gid-6,gid-1,gid-5,tmp_nas)
                                pPload(10,FM1,gid-4,tmp_nas)
                                pPload(10,FM2,gid-3,tmp_nas)
                                pPload(10,FM3,gid-2,tmp_nas)
                                pPload(10,FM4,gid-1,tmp_nas)
                                pPload(10,FM5,gid,tmp_nas)
                                if pr_tenel==1:
                                    pgrid(gid,center[0],center[1],center[2],tmp_ten)
                                    pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_ten)
                                    pgrid(gid+2,crossp[order[0]][1],crossp[order[0]][2],crossp[order[0]][3],tmp_ten)
                                    pgrid(gid+3,crossp[order[1]][1],crossp[order[1]][2],crossp[order[1]][3],tmp_ten)
                                    pgrid(gid+4,crossp[order[2]][1],crossp[order[2]][2],crossp[order[2]][3],tmp_ten)
                                    pgrid(gid+5,crossp[order[3]][1],crossp[order[3]][2],crossp[order[3]][3],tmp_ten)
                                    gid=gid+6
                                    ptria(gid-4,cbnrP,gid-6,gid-5,gid-4,tmp_ten)
                                    ptria(gid-3,cbnrP,gid-6,gid-4,gid-3,tmp_ten)
                                    ptria(gid-2,cbnrP,gid-6,gid-3,gid-2,tmp_ten)
                                    ptria(gid-1,cbnrP,gid-6,gid-2,gid-1,tmp_ten)
                                    ptria(gid,cbnrP,gid-6,gid-1,gid-5,tmp_ten)
                                    pTenForce(5,gid,sTVec,plVecStr,tmp_ten)
# end 5 cross points
#                            
# 6 cross points
                        if len(crossp) == 6:
                            v3x=crossp[3][1]-crossp[0][1]
                            v3y=crossp[3][2]-crossp[0][2]
                            v3z=crossp[3][3]-crossp[0][3]
                            v4x=crossp[4][1]-crossp[0][1]
                            v4y=crossp[4][2]-crossp[0][2]
                            v4z=crossp[4][3]-crossp[0][3]
                            v5x=crossp[5][1]-crossp[0][1]
                            v5y=crossp[5][2]-crossp[0][2]
                            v5z=crossp[5][3]-crossp[0][3]
                            alpha2=Sprod(v1x,v1y,v1z,v3x,v3y,v3z)
                            alpha3=Sprod(v1x,v1y,v1z,v4x,v4y,v4z)
                            alpha4=Sprod(v1x,v1y,v1z,v5x,v5y,v5z)
                            n2=Xprod(v1x,v1y,v1z,v3x,v3y,v3z)
                            n3=Xprod(v1x,v1y,v1z,v4x,v4y,v4z)
                            n4=Xprod(v1x,v1y,v1z,v5x,v5y,v5z)
                            order=nodord6(crossp,alpha1,alpha2,alpha3,alpha4,n1,n2,n3,n4,Np)
                            A1=Xprod(crossp[0][1]-center[0],crossp[0][2]-center[1],crossp[0][3]-center[2],crossp[order[0]][1]-center[0],crossp[order[0]][2]-center[1],crossp[order[0]][3]-center[2])[3]/2
                            A2=Xprod(crossp[order[0]][1]-center[0],crossp[order[0]][2]-center[1],crossp[order[0]][3]-center[2],crossp[order[1]][1]-center[0],crossp[order[1]][2]-center[1],crossp[order[1]][3]-center[2])[3]/2
                            A3=Xprod(crossp[order[1]][1]-center[0],crossp[order[1]][2]-center[1],crossp[order[1]][3]-center[2],crossp[order[2]][1]-center[0],crossp[order[2]][2]-center[1],crossp[order[2]][3]-center[2])[3]/2
                            A4=Xprod(crossp[order[2]][1]-center[0],crossp[order[2]][2]-center[1],crossp[order[2]][3]-center[2],crossp[order[3]][1]-center[0],crossp[order[3]][2]-center[1],crossp[order[3]][3]-center[2])[3]/2
                            A5=Xprod(crossp[order[3]][1]-center[0],crossp[order[3]][2]-center[1],crossp[order[3]][3]-center[2],crossp[order[4]][1]-center[0],crossp[order[4]][2]-center[1],crossp[order[4]][3]-center[2])[3]/2
                            A6=Xprod(crossp[order[4]][1]-center[0],crossp[order[4]][2]-center[1],crossp[order[4]][3]-center[2],crossp[0][1]-center[0],crossp[0][2]-center[1],crossp[0][3]-center[2])[3]/2
                            FM1=getForceMom(cro0,crossp[0],crossp[order[0]],A1,cro0s,crossS[0],crossS[order[0]],AT)
                            FM2=getForceMom(cro0,crossp[order[0]],crossp[order[1]],A2,cro0s,crossS[order[0]],crossS[order[1]],AT)
                            FM3=getForceMom(cro0,crossp[order[1]],crossp[order[2]],A3,cro0s,crossS[order[1]],crossS[order[2]],AT)
                            FM4=getForceMom(cro0,crossp[order[2]],crossp[order[3]],A4,cro0s,crossS[order[2]],crossS[order[3]],AT)
                            FM5=getForceMom(cro0,crossp[order[3]],crossp[order[4]],A5,cro0s,crossS[order[3]],crossS[order[4]],AT)
                            FM6=getForceMom(cro0,crossp[order[4]],crossp[0],A6,cro0s,crossS[order[4]],crossS[0],AT)
                            if foundPlst == 1 :
                                pl1=getPlStA(A1,center_pl,crossPlst[0],crossPlst[order[0]])
                                pl2=getPlStA(A2,center_pl,crossPlst[order[0]],crossPlst[order[1]])
                                pl3=getPlStA(A3,center_pl,crossPlst[order[1]],crossPlst[order[2]])
                                pl4=getPlStA(A4,center_pl,crossPlst[order[2]],crossPlst[order[3]])
                                pl5=getPlStA(A5,center_pl,crossPlst[order[3]],crossPlst[order[4]])
                                pl6=getPlStA(A6,center_pl,crossPlst[order[4]],crossPlst[0])
                                resPlst.append(pl1)
                                resPlst.append(pl2)
                                resPlst.append(pl3)
                                resPlst.append(pl4)
                                resPlst.append(pl5)
                                resPlst.append(pl6)
                                cnt_cbap=0
                                while cnt_cbap < 6 :
                                    resCBPl.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
# tension arm cb
                            if pr_tenel==1 :
                                pVec=[crossp[0],crossp[order[0]],crossp[order[1]],crossp[order[2]],crossp[order[3]],crossp[order[4]]]
                                sTVec=[crossT[0],crossT[order[0]],crossT[order[1]],crossT[order[2]],crossT[order[3]],crossT[order[4]]]
                                angleT=SprodS(Zp[0],Zp[1],Zp[2],(center[0]-OrigoY[0]),(center[1]-OrigoY[1]),(center[2]-OrigoY[2]),Np[0],Np[1],Np[2],0)
                                resTVec.append(pVec)
                                resSTen.append(sTVec)
                                resAngT.append(angleT)
                                resTenEl.append(elem_id)
                                resTenCB.append(cbnr)
                                resTenCp.append(len(crossp))
                                resTenCpD[elem_id]=(len(crossp))
                                resTVecD[elem_id]=pVec
                                resSTenD[elem_id]=sTVec
                            if refLc==1:
                                FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                FM5s=refinePatch(FM5[13],FM5[14],FM5[15],FM5[16],FM5[17],FM5[18],A5)
                                FM6s=refinePatch(FM6[13],FM6[14],FM6[15],FM6[16],FM6[17],FM6[18],A6)
                                resFM.extend(FM1s[:])
                                resFM.extend(FM2s[:])
                                resFM.extend(FM3s[:])
                                resFM.extend(FM4s[:])
                                resFM.extend(FM5s[:])
                                resFM.extend(FM6s[:])
                                cnt_cbap=0
                                while cnt_cbap < 24 :
                                    resCB.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
                            else:
                                resFM.append(FM1)
                                resFM.append(FM2)
                                resFM.append(FM3)
                                resFM.append(FM4)
                                resFM.append(FM5)
                                resFM.append(FM6)
                                cnt_cbap=0
                                while cnt_cbap < 6 :
                                    resCB.append(cbnr)
                                    cnt_cbap=cnt_cbap+1
                            if comb==1:
                                FM1a=getForceMom(cro0,crossp[0],crossp[order[0]],A1,cro0sa,cro0sa,cro0sa,AT)
                                FM2a=getForceMom(cro0,crossp[order[0]],crossp[order[1]],A2,cro0sa,cro0sa,cro0sa,AT)
                                FM3a=getForceMom(cro0,crossp[order[1]],crossp[order[2]],A3,cro0sa,cro0sa,cro0sa,AT)
                                FM4a=getForceMom(cro0,crossp[order[2]],crossp[order[3]],A4,cro0sa,cro0sa,cro0sa,AT)
                                FM5a=getForceMom(cro0,crossp[order[3]],crossp[order[4]],A5,cro0sa,cro0sa,cro0sa,AT)
                                FM6a=getForceMom(cro0,crossp[order[4]],crossp[0],A6,cro0sa,cro0sa,cro0sa,AT)
                                if refLc==1:
                                    FM1s=refinePatch(FM1a[13],FM1a[14],FM1a[15],FM1a[16],FM1a[17],FM1a[18],A1)
                                    FM2s=refinePatch(FM2a[13],FM2a[14],FM2a[15],FM2a[16],FM2a[17],FM2a[18],A2)
                                    FM3s=refinePatch(FM3a[13],FM3a[14],FM3a[15],FM3a[16],FM3a[17],FM3a[18],A3)
                                    FM4s=refinePatch(FM4a[13],FM4a[14],FM4a[15],FM4a[16],FM4a[17],FM4a[18],A4)
                                    FM5s=refinePatch(FM5a[13],FM5a[14],FM5a[15],FM5a[16],FM5a[17],FM5a[18],A5)
                                    FM6s=refinePatch(FM6a[13],FM6a[14],FM6a[15],FM6a[16],FM6a[17],FM6a[18],A6)
                                    resFMa.extend(FM1s[:])
                                    resFMa.extend(FM2s[:])
                                    resFMa.extend(FM3s[:])
                                    resFMa.extend(FM4s[:])
                                    resFMa.extend(FM5s[:])
                                    resFMa.extend(FM6s[:])
                                else:
                                    resFMa.append(FM1a)
                                    resFMa.append(FM2a)
                                    resFMa.append(FM3a)
                                    resFMa.append(FM4a)                            
                                    resFMa.append(FM5a)
                                    resFMa.append(FM6a)
                            if elset1x != '' and elem_id in elsetList:
                                resFMnp.append(FM1[14])
                                resFMnp.append(FM2[14])
                                resFMnp.append(FM3[14])
                                resFMnp.append(FM4[14])
                                resFMnp.append(FM5[14])
                                resFMnp.append(FM6[14])
                                resFMns.append(FM1[17])
                                resFMns.append(FM2[17])
                                resFMns.append(FM3[17])
                                resFMns.append(FM4[17])
                                resFMns.append(FM5[17])
                                resFMns.append(FM6[17])
                                resEL.append([elem_id,1])
                                resEL.append([elem_id,2])
                                resEL.append([elem_id,3])
                                resEL.append([elem_id,4])
                                resEL.append([elem_id,5])
                                resEL.append([elem_id,6])
                                if foundPlst == 1 :
                                   resPlstE.append(pl1)
                                   resPlstE.append(pl2)
                                   resPlstE.append(pl3)
                                   resPlstE.append(pl4)
                                   resPlstE.append(pl5)
                                   resPlstE.append(pl6)
                                   maxPlstE.append(max(crossPlst))
                                if comb==1:
                                    resFMec.append(FM1a)
                                    resFMec.append(FM2a)
                                    resFMec.append(FM3a)
                                    resFMec.append(FM4a)
                                    resFMec.append(FM5a)
                                    resFMec.append(FM6a)
                                    d=0
                                    while d < 6 :
                                       resFMeSetElc.append(dicS[elem_id])
                                       d=d+1
                                if refL1==0:
                                    d=0
                                    while d < 6 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    resFMe.append(FM1)
                                    resFMe.append(FM2)
                                    resFMe.append(FM3)
                                    resFMe.append(FM4)
                                    resFMe.append(FM5)
                                    resFMe.append(FM6)
                                if refL1==1:
                                    d=0
                                    while d < 24 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                    FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                    FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                    FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                    FM5s=refinePatch(FM5[13],FM5[14],FM5[15],FM5[16],FM5[17],FM5[18],A5)
                                    FM6s=refinePatch(FM6[13],FM6[14],FM6[15],FM6[16],FM6[17],FM6[18],A6)
                                    resFMe.extend(FM1s[:])
                                    resFMe.extend(FM2s[:])
                                    resFMe.extend(FM3s[:])
                                    resFMe.extend(FM4s[:])
                                    resFMe.extend(FM5s[:])
                                    resFMe.extend(FM6s[:])
                                if refL1==2:                        
                                    d=0
                                    while d < 96 :
                                       resFMeSetEl.append(dicS[elem_id])
                                       d=d+1
                                    FM1s=refinePatch(FM1[13],FM1[14],FM1[15],FM1[16],FM1[17],FM1[18],A1)
                                    FM2s=refinePatch(FM2[13],FM2[14],FM2[15],FM2[16],FM2[17],FM2[18],A2)
                                    FM3s=refinePatch(FM3[13],FM3[14],FM3[15],FM3[16],FM3[17],FM3[18],A3)
                                    FM4s=refinePatch(FM4[13],FM4[14],FM4[15],FM4[16],FM4[17],FM4[18],A4)
                                    FM5s=refinePatch(FM5[13],FM5[14],FM5[15],FM5[16],FM5[17],FM5[18],A5)
                                    FM6s=refinePatch(FM6[13],FM6[14],FM6[15],FM6[16],FM6[17],FM6[18],A6)
                                    FM1s2=refinePatch(FM1s[0][13],FM1s[0][14],FM1s[0][15],FM1s[0][16],FM1s[0][17],FM1s[0][18],A1/4)
                                    FM1s3=refinePatch(FM1s[1][13],FM1s[1][14],FM1s[1][15],FM1s[1][16],FM1s[1][17],FM1s[1][18],A1/4)
                                    FM1s4=refinePatch(FM1s[2][13],FM1s[2][14],FM1s[2][15],FM1s[2][16],FM1s[2][17],FM1s[2][18],A1/4)
                                    FM1s5=refinePatch(FM1s[3][13],FM1s[3][14],FM1s[3][15],FM1s[3][16],FM1s[3][17],FM1s[3][18],A1/4)
                                    FM2s2=refinePatch(FM2s[0][13],FM2s[0][14],FM2s[0][15],FM2s[0][16],FM2s[0][17],FM2s[0][18],A2/4)
                                    FM2s3=refinePatch(FM2s[1][13],FM2s[1][14],FM2s[1][15],FM2s[1][16],FM2s[1][17],FM2s[1][18],A2/4)
                                    FM2s4=refinePatch(FM2s[2][13],FM2s[2][14],FM2s[2][15],FM2s[2][16],FM2s[2][17],FM2s[2][18],A2/4)
                                    FM2s5=refinePatch(FM2s[3][13],FM2s[3][14],FM2s[3][15],FM2s[3][16],FM2s[3][17],FM2s[3][18],A2/4)
                                    FM3s2=refinePatch(FM3s[0][13],FM3s[0][14],FM3s[0][15],FM3s[0][16],FM3s[0][17],FM3s[0][18],A3/4)
                                    FM3s3=refinePatch(FM3s[1][13],FM3s[1][14],FM3s[1][15],FM3s[1][16],FM3s[1][17],FM3s[1][18],A3/4)
                                    FM3s4=refinePatch(FM3s[2][13],FM3s[2][14],FM3s[2][15],FM3s[2][16],FM3s[2][17],FM3s[2][18],A3/4)
                                    FM3s5=refinePatch(FM3s[3][13],FM3s[3][14],FM3s[3][15],FM3s[3][16],FM3s[3][17],FM3s[3][18],A3/4)
                                    FM4s2=refinePatch(FM4s[0][13],FM4s[0][14],FM4s[0][15],FM4s[0][16],FM4s[0][17],FM4s[0][18],A4/4)
                                    FM4s3=refinePatch(FM4s[1][13],FM4s[1][14],FM4s[1][15],FM4s[1][16],FM4s[1][17],FM4s[1][18],A4/4)
                                    FM4s4=refinePatch(FM4s[2][13],FM4s[2][14],FM4s[2][15],FM4s[2][16],FM4s[2][17],FM4s[2][18],A4/4)
                                    FM4s5=refinePatch(FM4s[3][13],FM4s[3][14],FM4s[3][15],FM4s[3][16],FM4s[3][17],FM4s[3][18],A4/4)
                                    FM5s2=refinePatch(FM5s[0][13],FM5s[0][14],FM5s[0][15],FM5s[0][16],FM5s[0][17],FM5s[0][18],A5/4)
                                    FM5s3=refinePatch(FM5s[1][13],FM5s[1][14],FM5s[1][15],FM5s[1][16],FM5s[1][17],FM5s[1][18],A5/4)
                                    FM5s4=refinePatch(FM5s[2][13],FM5s[2][14],FM5s[2][15],FM5s[2][16],FM5s[2][17],FM5s[2][18],A5/4)
                                    FM5s5=refinePatch(FM5s[3][13],FM5s[3][14],FM5s[3][15],FM5s[3][16],FM5s[3][17],FM5s[3][18],A5/4)
                                    FM6s2=refinePatch(FM6s[0][13],FM6s[0][14],FM6s[0][15],FM6s[0][16],FM6s[0][17],FM6s[0][18],A6/4)
                                    FM6s3=refinePatch(FM6s[1][13],FM6s[1][14],FM6s[1][15],FM6s[1][16],FM6s[1][17],FM6s[1][18],A6/4)
                                    FM6s4=refinePatch(FM6s[2][13],FM6s[2][14],FM6s[2][15],FM6s[2][16],FM6s[2][17],FM6s[2][18],A6/4)
                                    FM6s5=refinePatch(FM6s[3][13],FM6s[3][14],FM6s[3][15],FM6s[3][16],FM6s[3][17],FM6s[3][18],A6/4)
                                    resFMe.extend(FM1s2[:])
                                    resFMe.extend(FM1s3[:])
                                    resFMe.extend(FM1s4[:])
                                    resFMe.extend(FM1s5[:])
                                    resFMe.extend(FM2s2[:])
                                    resFMe.extend(FM2s3[:])
                                    resFMe.extend(FM2s4[:])
                                    resFMe.extend(FM2s5[:])
                                    resFMe.extend(FM3s2[:])
                                    resFMe.extend(FM3s3[:])
                                    resFMe.extend(FM3s4[:])
                                    resFMe.extend(FM3s5[:])
                                    resFMe.extend(FM4s2[:])
                                    resFMe.extend(FM4s3[:])
                                    resFMe.extend(FM4s4[:])
                                    resFMe.extend(FM4s5[:])
                                    resFMe.extend(FM5s2[:])
                                    resFMe.extend(FM5s3[:])
                                    resFMe.extend(FM5s4[:])
                                    resFMe.extend(FM5s5[:])
                                    resFMe.extend(FM6s2[:])
                                    resFMe.extend(FM6s3[:])
                                    resFMe.extend(FM6s4[:])
                                    resFMe.extend(FM6s5[:])
# Write to Nastran                               
                            if nas == 1:
                                pgrid(gid,center[0],center[1],center[2],tmp_nas)
                                pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_nas)
                                pgrid(gid+2,crossp[order[0]][1],crossp[order[0]][2],crossp[order[0]][3],tmp_nas)
                                pgrid(gid+3,crossp[order[1]][1],crossp[order[1]][2],crossp[order[1]][3],tmp_nas)
                                pgrid(gid+4,crossp[order[2]][1],crossp[order[2]][2],crossp[order[2]][3],tmp_nas)
                                pgrid(gid+5,crossp[order[3]][1],crossp[order[3]][2],crossp[order[3]][3],tmp_nas)
                                pgrid(gid+6,crossp[order[4]][1],crossp[order[4]][2],crossp[order[4]][3],tmp_nas)
                                gid=gid+7
                                ptria(gid-5,cbnrP,gid-7,gid-6,gid-5,tmp_nas)
                                ptria(gid-4,cbnrP,gid-7,gid-5,gid-4,tmp_nas)
                                ptria(gid-3,cbnrP,gid-7,gid-4,gid-3,tmp_nas)
                                ptria(gid-2,cbnrP,gid-7,gid-3,gid-2,tmp_nas)
                                ptria(gid-1,cbnrP,gid-7,gid-1,gid-6,tmp_nas)
                                ptria(gid,cbnrP,gid-7,gid-2,gid-1,tmp_nas)
                                pPload(10,FM1,gid-5,tmp_nas)
                                pPload(10,FM2,gid-4,tmp_nas)
                                pPload(10,FM3,gid-3,tmp_nas)
                                pPload(10,FM4,gid-2,tmp_nas)
                                pPload(10,FM5,gid-1,tmp_nas)
                                pPload(10,FM6,gid,tmp_nas)
                                if pr_tenel==1:
                                    pgrid(gid,center[0],center[1],center[2],tmp_ten)
                                    pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_ten)
                                    pgrid(gid+2,crossp[order[0]][1],crossp[order[0]][2],crossp[order[0]][3],tmp_ten)
                                    pgrid(gid+3,crossp[order[1]][1],crossp[order[1]][2],crossp[order[1]][3],tmp_ten)
                                    pgrid(gid+4,crossp[order[2]][1],crossp[order[2]][2],crossp[order[2]][3],tmp_ten)
                                    pgrid(gid+5,crossp[order[3]][1],crossp[order[3]][2],crossp[order[3]][3],tmp_ten)
                                    pgrid(gid+6,crossp[order[4]][1],crossp[order[4]][2],crossp[order[4]][3],tmp_ten)
                                    gid=gid+7
                                    ptria(gid-5,cbnrP,gid-7,gid-6,gid-5,tmp_ten)
                                    ptria(gid-4,cbnrP,gid-7,gid-5,gid-4,tmp_ten)
                                    ptria(gid-3,cbnrP,gid-7,gid-4,gid-3,tmp_ten)
                                    ptria(gid-2,cbnrP,gid-7,gid-3,gid-2,tmp_ten)
                                    ptria(gid-1,cbnrP,gid-7,gid-1,gid-6,tmp_ten)
                                    ptria(gid,cbnrP,gid-7,gid-2,gid-1,tmp_ten)
                                    pTenForce(6,gid,sTVec,plVecStr,tmp_ten)
# end 6 cross points
#-----------------------------------------------------------------                
                        if len(crossp) > 7:
                            print('Number of edge cuts larger than 6!')
                stim5=time.clock()
                time_el=time_el+stim5-stim4
                if ni == 1 or irem==1:
                    boxNodDic[curSec]=dict(boxNodDicT)
                    boxElDic[curSec]=dict(boxElDicT)
                if nas == 1:
                    for ut in range(0,len(CBNR)):
                        tmp='PSHELL  %8i       1     1.0       1               1\n' %(CBNR[ut]+100*ni)
                        tmp_nas.write(tmp)
                    if pr_tenel == 1:
                        for ut in range(0,len(CBNR)):
                            tmp='PSHELL  %8i       1     1.0       1               1\n' %(CBNR[ut]+100*ni)
                            tmp_ten.write(tmp)
    
#
# complicated part to find wires and stresses that considers the possibility of a 
# general cross section cut...
#
                if tenCB !='' :
                    for cwnr in splTeni :
                        el1=[]
                        el2=[]
                        angT=[]
                        angT2=[]
                        slave=[]
                        cntTen=0
                        single1=0
                        res4str={}
                        res4pos={}
                        resAng=[]
                        w1d={}
                        w1n=0
                        w2d={}
                        while cntTen < len(resTenEl) :
                            if resTenCB[cntTen] == cwnr:
                                cntT2=0
                                single1=0
                                while cntT2 < len(resTenEl) :
                                    if cntTen != cntT2 :
                                        cntT3=0
                                        while cntT3 < resTenCp[cntTen] :
                                            cntT4=0
                                            while cntT4 < resTenCp[cntT2] :
                                                yd=abs(resTVec[cntTen][cntT3][dof1+1]-resTVec[cntT2][cntT4][dof1+1])
                                                if yd < 0.001*sca :
                                                    zd=abs(resTVec[cntTen][cntT3][dof2+1]-resTVec[cntT2][cntT4][dof2+1])
                                                    if zd < 0.001*sca :
                                                        if (resTenEl[cntTen] in el1) and ((resTenEl[cntT2] not in el2) or (resTenEl[cntT2] not in el1)) :
                                                            w1d[resTenEl[cntT2]]=w1d[resTenEl[cntTen]]
                                                        elif (resTenEl[cntT2] in el1) and ((resTenEl[cntTen] not in el2) or (resTenEl[cntTen] not in el1)) :
                                                            w1d[resTenEl[cntTen]]=w1d[resTenEl[cntT2]]
                                                        elif (resTenEl[cntTen] in el2) and ((resTenEl[cntT2] not in el1) or (resTenEl[cntT2] not in el2)) :
                                                            w1d[resTenEl[cntT2]]=w1d[resTenEl[cntTen]]
                                                        elif (resTenEl[cntT2] in el2) and ((resTenEl[cntTen] not in el1) or (resTenEl[cntTen] not in el2)):
                                                            w1d[resTenEl[cntTen]]=w1d[resTenEl[cntT2]]
                                                        elif (((resTenEl[cntTen] in el1) or (resTenEl[cntTen] in el2)) and ((resTenEl[cntT2] in el1) or (resTenEl[cntT2] in el2))) :
                                                            w1n=w1n
                                                        else:
                                                            w1n=w1n+1
                                                            w1d[resTenEl[cntTen]]=w1n
                                                            w1d[resTenEl[cntT2]]=w1n
                                                        el1.append(resTenEl[cntTen])
                                                        el2.append(resTenEl[cntT2])
                                                        angT.append(resAngT[cntTen])
                                                        angT2.append(resAngT[cntT2])
                                                        single1=1
                                                cntT4=cntT4+1
                                            cntT3=cntT3+1
                                    cntT2=cntT2+1
                                if single1==0 :
                                    res4str[resAngT[cntTen]]=resSTen[cntTen]
                                    res4pos[resAngT[cntTen]]=resTVec[cntTen]
                                    resAng.append(resAngT[cntTen])
                            cntTen=cntTen+1
                        cntAng=0
                        for ang in resAng :
                            strTenW=[]
                            angTemp=[]
                            radTemp=[]
                            cnt=0
                            for pnt in res4pos[ang] :
                                pert=1.0e-6*cnt
                                angxx=SprodS(Zp[0],Zp[1],Zp[2],(pnt[1]-OrigoY[0]),(pnt[2]-OrigoY[1]),(pnt[3]-OrigoY[2]),Np[0],Np[1],Np[2],0)
                                angTemp.append((angxx+pert))
                                radxx=sqrt((pnt[1]-OrigoY[0])**2+(pnt[2]-OrigoY[1])**2+(pnt[3]-OrigoY[2])**2)
                                radTemp.append((radxx+(pert*sca)))
                                cnt=cnt+1
                            if (max(angTemp)-min(angTemp)) > 100.0:
                                cntFixA=0
                                while cntFixA < len(angTemp) :
                                    if angTemp[cntFixA] < 0.0 :
                                        angTemp[cntFixA]=angTemp[cntFixA]+360.0
                                    cntFixA=cntFixA+1
                            sm2=nsmallest(2,angTemp)
                            sm2r=nsmallest(2,radTemp)
                            lm2=nlargest(2,angTemp)
                            lm2r=nlargest(2,radTemp)
# find I1
                            foundI1=0
                            ind=0
                            while ind < len(sm2r) :
                                if angTemp.index(sm2[0]) == radTemp.index(sm2r[ind]) :
                                    foundI1=angTemp.index(sm2[0])
                                ind=ind+1
                            ind=0
                            if foundI1==0:
                                while ind < len(sm2r) :
                                    if angTemp.index(sm2[1]) == radTemp.index(sm2r[ind]) :
                                        foundI1=angTemp.index(sm2[1])
                                    ind=ind+1
# find O1
                            foundO1=0
                            ind=0
                            while ind < len(lm2r) :
                                if angTemp.index(sm2[0]) == radTemp.index(lm2r[ind]) :
                                    foundO1=angTemp.index(sm2[0])
                                ind=ind+1
                            ind=0
                            if foundO1==0:
                                while ind < len(lm2r) :
                                    if angTemp.index(sm2[1]) == radTemp.index(lm2r[ind]) :
                                        foundO1=angTemp.index(sm2[1])
                                    ind=ind+1
# find I2
                            foundI2=0
                            ind=0
                            while ind < len(sm2r) :
                                if angTemp.index(lm2[0]) == radTemp.index(sm2r[ind]) :
                                    foundI2=angTemp.index(lm2[0])
                                ind=ind+1
                            ind=0
                            if foundI2==0:
                                while ind < len(sm2r) :
                                    if angTemp.index(lm2[1]) == radTemp.index(sm2r[ind]) :
                                        foundI2=angTemp.index(lm2[1])
                                    ind=ind+1
# find O2
                            foundO2=0
                            ind=0
                            while ind < len(lm2r) :
                                if angTemp.index(lm2[0]) == radTemp.index(lm2r[ind]) :
                                    foundO2=angTemp.index(lm2[0])
                                ind=ind+1
                            ind=0
                            if foundO2==0:
                                while ind < len(lm2r) :
                                    if angTemp.index(lm2[1]) == radTemp.index(lm2r[ind]) :
                                        foundO2=angTemp.index(lm2[1])
                                    ind=ind+1
                            strTW=res4str[ang]
                            strTenW.extend([strTW[foundO1],strTW[foundO2],strTW[foundI1],strTW[foundI2]])
# make sure no negative angles are sorted and written...
                            if ang < 0.0 :
                                ang=360.0+ang
                                resAng[cntAng] = ang
                            res4str[ang]=strTenW
                            cntAng=cntAng+1    
#
# sort out the subdivided wires...
#
                        cnt9=1
                        while cnt9 < (w1n+1) :
                            w2d[cnt9]=[]
                            cnt9=cnt9+1
                        for ele in w1d :
                            w2d[w1d[ele]].append(ele)
                        for wir1 in w2d :
                            pList1=[]
                            STlist1=[]
                            angW=[]
                            radW=[]
                            posW=[]
                            strTW=[]
                            strTenW=[]
                            for elew in w2d[wir1] :
                                cnt10=0
                                while cnt10 < resTenCpD[elew] :
                                    currP=resTVecD[elew][cnt10]
                                    currST=resSTenD[elew][cnt10]
                                    pert=1.0e-12*cnt10*elew
                                    angleT1=SprodS(Zp[0],Zp[1],Zp[2],(currP[1]-OrigoY[0]),(currP[2]-OrigoY[1]),(currP[3]-OrigoY[2]),Np[0],Np[1],Np[2],0)
                                    angW.append(angleT1+pert)
                                    radiusT1=sqrt((currP[1]-OrigoY[0])**2+(currP[2]-OrigoY[1])**2+(currP[3]-OrigoY[2])**2)
                                    radW.append(radiusT1+(pert*sca))
                                    strTW.append(currST)
                                    posW.append(currP)
                                    cnt10=cnt10+1
                            if (max(angW)-min(angW)) > 100.0:
                                cntFixA=0
                                while cntFixA < len(angW) :
                                    if angW[cntFixA] < 0.0 :
                                        angW[cntFixA]=angW[cntFixA]+360.0
                                    cntFixA=cntFixA+1
                            sm2=nsmallest(3,angW)
                            smx=nsmallest((len(radW)/2),radW)
                            lm2=nlargest(3,angW)
                            lmx=nlargest((len(radW)/2),radW)
# find I1
                            foundI1=0
                            ind=0
                            while ind < len(smx) :
                                if angW.index(sm2[0]) == radW.index(smx[ind]) :
                                    foundI1=angW.index(sm2[0])
                                ind=ind+1
                            ind=0
                            if foundI1==0:
                                while ind < len(smx) :
                                    if angW.index(sm2[1]) == radW.index(smx[ind]) :
                                        foundI1=angW.index(sm2[1])
                                    ind=ind+1
                            ind=0
                            if foundI1==0:
                                while ind < len(smx) :
                                    if angW.index(sm2[2]) == radW.index(smx[ind]) :
                                        foundI1=angW.index(sm2[2])
                                    ind=ind+1
# find O1
                            foundO1=0
                            ind=0
                            while ind < len(lmx) :
                                if angW.index(sm2[0]) == radW.index(lmx[ind]) :
                                    foundO1=angW.index(sm2[0])
                                ind=ind+1
                            ind=0
                            if foundO1==0:
                                while ind < len(lmx) :
                                    if angW.index(sm2[1]) == radW.index(lmx[ind]) :
                                        foundO1=angW.index(sm2[1])
                                    ind=ind+1
                            ind=0
                            if foundO1==0:
                                while ind < len(lmx) :
                                    if angW.index(sm2[2]) == radW.index(lmx[ind]) :
                                        foundO1=angW.index(sm2[2])
                                    ind=ind+1
# find I2
                            foundI2=0
                            ind=0
                            while ind < len(smx) :
                                if angW.index(lm2[0]) == radW.index(smx[ind]) :
                                    foundI2=angW.index(lm2[0])
                                ind=ind+1
                            ind=0
                            if foundI2==0:
                                while ind < len(smx) :
                                    if angW.index(lm2[1]) == radW.index(smx[ind]) :
                                        foundI2=angW.index(lm2[1])
                                    ind=ind+1
                            ind=0
                            if foundI2==0:
                                while ind < len(smx) :
                                    if angW.index(lm2[2]) == radW.index(smx[ind]) :
                                        foundI2=angW.index(lm2[2])
                                    ind=ind+1
# find O2
                            foundO2=0
                            ind=0
                            while ind < len(lmx) :
                                if angW.index(lm2[0]) == radW.index(lmx[ind]) :
                                    foundO2=angW.index(lm2[0])
                                ind=ind+1
                            ind=0
                            if foundO2==0:
                                while ind < len(lmx) :
                                    if angW.index(lm2[1]) == radW.index(lmx[ind]) :
                                        foundO2=angW.index(lm2[1])
                                    ind=ind+1
                            ind=0
                            if foundO2==0:
                                while ind < len(lmx) :
                                    if angW.index(lm2[2]) == radW.index(lmx[ind]) :
                                        foundO2=angW.index(lm2[2])
                                    ind=ind+1
                            strTenW.extend([strTW[foundO1],strTW[foundO2],strTW[foundI1],strTW[foundI2]])
                            newCPx=(posW[foundO1][1]+posW[foundO2][1]+posW[foundI1][1]+posW[foundI2][1])/4
                            newCPy=(posW[foundO1][2]+posW[foundO2][2]+posW[foundI1][2]+posW[foundI2][2])/4
                            newCPz=(posW[foundO1][3]+posW[foundO2][3]+posW[foundI1][3]+posW[foundI2][3])/4
                            newAngW=SprodS(Zp[0],Zp[1],Zp[2],(newCPx-OrigoY[0]),(newCPy-OrigoY[1]),(newCPz-OrigoY[2]),Np[0],Np[1],Np[2],1)
                            resAng.append(newAngW)
                            res4str[newAngW]=strTenW
#
                        sortAng=[i[0] for i in sorted(enumerate(resAng), key=lambda x:x[1])]
                        outwnam=ten_based[curSec] + str(cwnr) + '.csv'
                        num_w=len(sortAng)
                        secAng=360.0/num_w
                        if ni == 1 :
                            initAngd[cwnr]=resAng[sortAng[0]]
                            outaa=open(outwnam,'w')
                            outaa.write('inc,time,Curvature K1,Curvature K2,Curvature K3,Curvature K4,Curvature K1_2,Curvature K3_4')
                            for wir in resAng :
                                outaa.write(',Wire number,Angle,Stress O1,Stress O2,Stress I1,Stress I2')
                            outaa.write('\n')
                        else:
                            outaa=open(outwnam,'a')
                        cntWW=1
                        tmp='%i,%6.3f,%10.8e,%10.8e,%10.8e,%10.8e,%10.8e,%10.8e' %((ni-1),atime,K1,K2,K3,K4,K1_2,K3_4)
                        outaa.write(tmp)
                        shiftFlag=0
                        while cntWW < (num_w+1) :
                            ind=sortAng[cntWW-1]
                            if (abs(resAng[ind]-initAngd[cwnr])/secAng > twistTol) and cntWW==1:
                                print('shifted wire')
                                shiftFlag=1
#                        if shiftFlag==1:
                                if resAng[ind] < initAngd[cwnr] :
                                    sortAngNew=sortAng[1:]
                                    sortAngNew.append(sortAng[0])
                                elif resAng[ind] > initAngd[cwnr] :
                                    sortAngNew=[sortAng[-1]]
                                    sortAngNew.extend(sortAng[:-1])
                                sortAng=sortAngNew
                                ind=sortAng[cntWW-1]
                            tmp=',%3i,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f' %(cntWW,resAng[ind],res4str[resAng[ind]][0],res4str[resAng[ind]][1],res4str[resAng[ind]][2],res4str[resAng[ind]][3])
                            outaa.write(tmp)
                            cntWW+=1
                        outaa.write('\n')
                        outaa.close()
                stim6=time.clock()
                time_ten=time_ten+stim6-stim5
    
#-------------------------------------------------------------------------    
# if summation is to be done on contact body levle cbpr != 'none'
                if cbpr != 'none':
                    if txt==1 and ni==1:
                        tmp_csv.write('inc,time,K1,K2,K3,K4,K1_2,K3_4,Ovality Y,Origo X GCS,Origo Y GCS,Origo Z GCS,Origo X LCS 0,Origo Y LCS 0,Origo Z LCS 0,Mid N1&N2 X LCS 0,Mid N1&N2 Y LCS 0,Mid N1&N2 Z LCS 0,Iyy,Izz,Iyz,Area.SUM,Fx.SUM,Fy.SUM,Fz.SUM,Mx.SUM,My.SUM,Mz.SUM,Av Pl Strain.SUM,Max Pl Strain.SUM,Max Pl Strain All.SUM,At Element.SUM,At Node.SUM,')    
                        if cbpr == '':
                            for ub in range(0,len(CBNR)):
                                currCbx=str(CBNM[ub]).replace(' ','')
                                if currCbx == 'None' and prNone == 1:
                                    tmp_csv.write('Area.%s,Fx.%s,Fy.%s,Fz.%s,Mx.%s,My.%s,Mz.%s,Av Pl Strain.%s,Max Pl Strain.%s,' %(currCbx,currCbx,currCbx,currCbx,currCbx,currCbx,currCbx,currCbx,currCbx))
                                elif currCbx != 'None' :
                                    tmp_csv.write('Area.%s,Fx.%s,Fy.%s,Fz.%s,Mx.%s,My.%s,Mz.%s,Av Pl Strain.%s,Max Pl Strain.%s,' %(currCbx,currCbx,currCbx,currCbx,currCbx,currCbx,currCbx,currCbx,currCbx))
                        else:
                            splcb=cbpr.split(',')
                            for dcb in splcb:
                                dcb=int(dcb)
                                indCb=CBNR.index(dcb)
                                print(CBNR[indCb],CBNM[indCb])
                                if CBNM[indCb] == 'None' and prNone == 1:
                                    tmp_csv.write('Area.%s,Fx.%s,Fy.%s,Fz.%s,Mx.%s,My.%s,Mz.%s,Av Pl Strain.%s,Max Pl Strain.%s,' %(CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb]))
                                elif CBNM[indCb] != 'None' :
                                    tmp_csv.write('Area.%s,Fx.%s,Fy.%s,Fz.%s,Mx.%s,My.%s,Mz.%s,Av Pl Strain.%s,Max Pl Strain.%s,' %(CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb],CBNM[indCb]))
                        tmp_csv.write('\n')
                    sumFx=0
                    sumFy=0
                    sumFz=0
                    sumMx=0
                    sumMy=0
                    sumMz=0
                    sumPlA=0.0
                    maxPlstA=0.0
                    Atot=sum([l[9] for l in resFM])
                    if sym == 0:
                        oriX=sum([l[10] for l in resFM])/Atot
                        oriY=sum([l[11] for l in resFM])/Atot
                        oriZ=sum([l[12] for l in resFM])/Atot
                    else:
                        oriX=(p2[0]+p1[0])/2
                        oriY=(p2[1]+p1[1])/2
                        oriZ=(p2[2]+p1[2])/2
                    ori2X=(p2[0]+p1[0])/2
                    ori2Y=(p2[1]+p1[1])/2
                    ori2Z=(p2[2]+p1[2])/2
                    if nas == 1:
                        tmp='CORD2R  %8i       0%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n' %(100*ni,oriX,oriY,oriZ,Zp[0]+oriX,Zp[1]+oriY,Zp[2]+oriZ)
                        tmp_nas.write(tmp)
                        tmp='        %8.2f%8.2f%8.2f\n' %(Np[0]+oriX,Np[1]+oriY,Np[2]+oriZ)
                        tmp_nas.write(tmp)
                    if ni == 1 :
                        ori0x[curSec]=oriX
                        ori0y[curSec]=oriY
                        ori0z[curSec]=oriZ
                        AT0[curSec]=AT
                    Iyy=0
                    Izz=0
                    Iyz=0
                    dicCBF={}
                    dicCBM={}
                    dicCBA={}
                    dicCBPl={}
                    dicCBPl2={}
                    for dcb in range(0,len(CBNR)):
                        dicCBF[CBNR[dcb]]=[0.0,0.0,0.0]
                        dicCBM[CBNR[dcb]]=[0.0,0.0,0.0]
                        dicCBA[CBNR[dcb]]=[0.0]
                        dicCBPl[CBNR[dcb]]=[0.0]
                        dicCBPl2[CBNR[dcb]]=0.0
                    for w in range(0,len(resFM)):
                        L1=getLocalVec(resFM[w][0],resFM[w][1],resFM[w][2],oriX,oriY,oriZ,AT)
                        Iyy=Iyy+L1[2]*L1[2]*resFM[w][9]
                        Izz=Izz+L1[1]*L1[1]*resFM[w][9]
                        Iyz=Iyz+L1[1]*L1[2]*resFM[w][9]
                        sumMx=sumMx+(resFM[w][7]*L1[2])-(resFM[w][8]*L1[1])
                        sumMy=sumMy+resFM[w][6]*L1[2]
                        sumMz=sumMz+resFM[w][6]*L1[1]
                        if comb==1:
                            sumFx=sumFx+resFMa[w][6]
                            sumFy=sumFy+resFMa[w][7]
                            sumFz=sumFz+resFMa[w][8]
                            dicCBF[resCB[w]]=[dicCBF[resCB[w]][0]+resFMa[w][6],dicCBF[resCB[w]][1]+resFMa[w][7],dicCBF[resCB[w]][2]+resFMa[w][8]]
                        else:
                            sumFx=sumFx+resFM[w][6]
                            sumFy=sumFy+resFM[w][7]
                            sumFz=sumFz+resFM[w][8]
                            dicCBF[resCB[w]]=[dicCBF[resCB[w]][0]+resFM[w][6],dicCBF[resCB[w]][1]+resFM[w][7],dicCBF[resCB[w]][2]+resFM[w][8]]
                        dicCBM[resCB[w]]=[dicCBM[resCB[w]][0]+(resFM[w][7]*L1[2])-(resFM[w][8]*L1[1]),dicCBM[resCB[w]][1]+resFM[w][6]*L1[2],dicCBM[resCB[w]][2]+resFM[w][6]*L1[1]]
                        dicCBA[resCB[w]]=[dicCBA[resCB[w]][0]+(resFM[w][9])]
                    if foundPlst == 1 :
                        for w in range(0,len(resPlst)):
                            sumPlA=sumPlA+resPlst[w]
                            dicCBPl[resCBPl[w]]=[dicCBPl[resCBPl[w]][0]+(resPlst[w])]
                        for w in range(0,len(maxPlst)):
                            dicCBPl2[resCBPl2[w]]=max(dicCBPl2[resCBPl2[w]],(maxPlst[w]))
                        maxPlstA=0.0
                        for key in dicCBPl2:
                            if dicCBPl2[key] > maxPlstA :
                                maxPlstA=dicCBPl2[key]
#                    print dicCBPl2,maxPlstA,dicCBPl
                    cogL=getLocalVec(oriX,oriY,oriZ,ori0x[curSec],ori0y[curSec],ori0z[curSec],AT0[curSec])
                    cogL2=getLocalVec(ori2X,ori2Y,ori2Z,ori0x[curSec],ori0y[curSec],ori0z[curSec],AT0[curSec])
                    tmp='Origin X      %10.3f          %10.3f              %10.3f\n' %(oriX,cogL[0],cogL2[0])
                    tmp_txt.write(tmp)
                    tmp='Origin Y      %10.3f          %10.3f              %10.3f\n' %(oriY,cogL[1],cogL2[1])
                    tmp_txt.write(tmp)
                    tmp='Origin Z      %10.3f          %10.3f              %10.3f\n' %(oriZ,cogL[2],cogL2[2])
                    tmp_txt.write(tmp)
                    tmp='Area          %10.3f\n' %(Atot)
                    tmp_txt.write(tmp)
                    tmp='Moment of Inertia Iyy  %14.3f\n' %(Iyy)
                    tmp_txt.write(tmp)
                    tmp='Moment of Inertia Izz  %14.3f\n' %(Izz)
                    tmp_txt.write(tmp)
                    tmp='Moment of Inertia Iyz  %14.3f\n\n' %(Iyz)
                    tmp_txt.write(tmp)
                    if genRem == 0 :
                        tmp_txt.write('                                         Cross Section Curvature \n')
                        tmp_txt.write('---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                        tmp='Pipe curvature at section node 1 (K1)       %10.8e\n' %(K1)
                        tmp_txt.write(tmp)
                        tmp='Pipe curvature at section node 2 (K2)       %10.8e\n' %(K2)
                        tmp_txt.write(tmp)
                        tmp='Pipe curvature at section node 3 (K3)       %10.8e\n' %(K3)
                        tmp_txt.write(tmp)
                        tmp='Pipe curvature at section node 4 (K4)       %10.8e\n' %(K4)
                        tmp_txt.write(tmp)
                        tmp='Pipe curvature average of K1 and K2 (K1_2)  %10.8e\n' %(K1_2)
                        tmp_txt.write(tmp)
                        tmp='Pipe curvature average of K3 and K4 (K3_4)  %10.8e\n' %(K3_4)
                        tmp_txt.write(tmp)
                        tmp='Ovality in local Y axis                     %10.8e\n\n' %(Oval1)
                        tmp_txt.write(tmp)
         #               tmp='Ovality in local Z axis                     %10.8e\n\n' %(Oval2)
         #               tmp_txt.write(tmp)
                    else :
                        tmp_txt.write('                                         Cross Section Curvature (Disabled due to Global Remeshing)\n')
                        tmp_txt.write('---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n')
                    tmp_txt.write('                                         Cross Section Resulting Forces and Moments\n')
                    tmp_txt.write('---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                    tmp_txt.write('Body Number  Body Name                            Area        Force X        Force Y        Force Z       Moment X       Moment Y       Moment Z    Av Pl Strain   Max Pl Strain   Max Pl Str All @Element    @Node\n')
                    if cbpr == '':
                        tmp='Sum of all   Sum of all                     %10.3f %14.2f %14.2f %14.2f %14.2f %14.2f %14.2f  %14.8e  %14.8e   %14.8e %8i %8i\n' %(Atot,sumFx,sumFy,sumFz,sumMx,sumMy,sumMz,sumPlA/Atot,maxPlstA,maxPlstAll,maxPlstEl,maxPlstNode) 
                        tmp_txt.write(tmp)
                        if txt==1:
                            tmpc='%5i,%4.4f,%8.6e,%8.6e,%8.6e,%8.6e,%8.6e,%8.6e,%8.6e,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e,%14.8e,%8i,%8i,' %(ni-1,atime,K1,K2,K3,K4,K1_2,K3_4,Oval1,oriX,oriY,oriZ,cogL[0],cogL[1],cogL[2],cogL2[0],cogL2[1],cogL2[2],Iyy,Izz,Iyz,Atot,sumFx,sumFy,sumFz,sumMx,sumMy,sumMz,sumPlA/Atot,maxPlstA,maxPlstAll,maxPlstEl,maxPlstNode)
                            tmp_csv.write(tmpc)
                        for dcb in range(0,len(CBNR)):
                            if not( CBNR[dcb] == 0 and prNone == 0 ) :
                                try :
                                   plStDivA=dicCBPl[CBNR[dcb]][0]/dicCBA[CBNR[dcb]][0]
                                except :
                                   plStDivA=0.0
                                tmp='%10i   %-30s %10.3f %14.2f %14.2f %14.2f %14.2f %14.2f %14.2f  %14.8e  %14.8e\n' %(CBNR[dcb],CBNM[dcb],dicCBA[CBNR[dcb]][0],dicCBF[CBNR[dcb]][0],dicCBF[CBNR[dcb]][1],dicCBF[CBNR[dcb]][2],dicCBM[CBNR[dcb]][0],dicCBM[CBNR[dcb]][1],dicCBM[CBNR[dcb]][2],plStDivA,dicCBPl2[CBNR[dcb]])
                                tmp_txt.write(tmp)
                                if txt==1:
                                    tmpc='%4.4f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e,' %(dicCBA[CBNR[dcb]][0],dicCBF[CBNR[dcb]][0],dicCBF[CBNR[dcb]][1],dicCBF[CBNR[dcb]][2],dicCBM[CBNR[dcb]][0],dicCBM[CBNR[dcb]][1],dicCBM[CBNR[dcb]][2],plStDivA,dicCBPl2[CBNR[dcb]])
                                    tmp_csv.write(tmpc)
                        tmp_csv.write('\n')
                        if txt==1:
                            tmp_txt.write('\n')
                    else:
                        Atot=0
                        sumFx=0
                        sumFy=0
                        sumFz=0
                        sumMx=0
                        sumMy=0
                        sumMz=0
                        sumPlA=0.0
                        maxPlstA=0.0
                        splcb=cbpr.split(',')
                        for sumCB in splcb:
                            Atot=Atot+dicCBA[int(sumCB)][0]
                            sumFx=sumFx+dicCBF[int(sumCB)][0]
                            sumFy=sumFy+dicCBF[int(sumCB)][1]
                            sumFz=sumFz+dicCBF[int(sumCB)][2]
                            sumMx=sumMx+dicCBM[int(sumCB)][0]
                            sumMy=sumMy+dicCBM[int(sumCB)][1]
                            sumMz=sumMz+dicCBM[int(sumCB)][2]
                            if foundPlst == 1 :
                               sumPlA=sumPlA+dicCBPl[int(sumCB)][0]
                               if dicCBPl2[int(sumCB)] > maxPlstA :
                                      maxPlstA=dicCBPl2[int(sumCB)]
                        tmp='Sum selected Sum selected                   %10.3f %14.2f %14.2f %14.2f %14.2f %14.2f %14.2f  %14.8e  %14.8e\n' %(Atot,sumFx,sumFy,sumFz,sumMx,sumMy,sumMz,sumPlA/Atot,maxPlstA) 
                        tmp_txt.write(tmp)
                        if txt==1:
                            tmpc='%5i,%4.4f,%8.6e,%8.6e,%8.6e,%8.6e,%8.6e,%8.6e,%8.6e,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e,%14.8e,%8i,%8i,' %(ni-1,atime,K1,K2,K3,K4,K1_2,K3_4,Oval1,oriX,oriY,oriZ,cogL[0],cogL[1],cogL[2],cogL2[0],cogL2[1],cogL2[2],Iyy,Izz,Iyz,Atot,sumFx,sumFy,sumFz,sumMx,sumMy,sumMz,sumPlA/Atot,maxPlstA,maxPlstAll,maxPlstEl,maxPlstNode)
                            tmp_csv.write(tmpc)
                        for dcb in splcb:
#                            dcb=int(dcb)-1 # don't understand why I did this..
                            dcb=int(dcb) 
                            indCb=CBNR.index(dcb)
                            dcb=indCb
                            try :
                               plStDivA=dicCBPl[CBNR[dcb]][0]/dicCBA[CBNR[dcb]][0]
                            except :
                               plStDivA=0.0
                            tmp='%10i   %-30s %10.3f %14.2f %14.2f %14.2f %14.2f %14.2f %14.2f  %14.8e  %14.8e \n' %(CBNR[dcb],CBNM[dcb],dicCBA[CBNR[dcb]][0],dicCBF[CBNR[dcb]][0],dicCBF[CBNR[dcb]][1],dicCBF[CBNR[dcb]][2],dicCBM[CBNR[dcb]][0],dicCBM[CBNR[dcb]][1],dicCBM[CBNR[dcb]][2],plStDivA,dicCBPl2[CBNR[dcb]])
                            tmp_txt.write(tmp)
                            if txt==1:
                                tmpc='%4.4f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e,' %(dicCBA[CBNR[dcb]][0],dicCBF[CBNR[dcb]][0],dicCBF[CBNR[dcb]][1],dicCBF[CBNR[dcb]][2],dicCBM[CBNR[dcb]][0],dicCBM[CBNR[dcb]][1],dicCBM[CBNR[dcb]][2],plStDivA,dicCBPl2[CBNR[dcb]])
                                tmp_csv.write(tmpc)
                        tmp_csv.write('\n')
                        if txt==1:
                            tmp_txt.write('\n')
#
#-------------------------------------------------------------------------    
# if summation is to be done on element set level elset1!= ''
                if elset1x != '':
                    sumFxE={}
                    sumFyE={}
                    sumFzE={}
                    sumMxE={}
                    sumMyE={}
                    sumMzE={}
                    sumPlAE={}
                    maxPlstAE={}
                    Atot={}
                    oriX={}
                    oriY={}
                    oriZ={}
                    Iyy={}
                    Izz={}
                    Iyz={}
                    for snam in elset1 :
                        sumFxE[snam]=0.0
                        sumFyE[snam]=0.0
                        sumFzE[snam]=0.0
                        sumMxE[snam]=0.0
                        sumMyE[snam]=0.0
                        sumMzE[snam]=0.0
                        sumPlAE[snam]=0.0
                        maxPlstAE[snam]=0.0
                        Atot[snam]=0.0
                        oriX[snam]=0.0
                        oriY[snam]=0.0
                        oriZ[snam]=0.0
                        Iyy[snam]=0.0
                        Izz[snam]=0.0
                        Iyz[snam]=0.0
                    d=0
                    while d < len(resFMeSetEl) :
                        Atot[resFMeSetEl[d]]=Atot[resFMeSetEl[d]]+resFMe[d][9]
                        d = d + 1
                    d=0
                    while d < len(resFMeSetEl) :
                        oriX[resFMeSetEl[d]]=oriX[resFMeSetEl[d]]+(resFMe[d][10]/Atot[resFMeSetEl[d]])
                        oriY[resFMeSetEl[d]]=oriY[resFMeSetEl[d]]+(resFMe[d][11]/Atot[resFMeSetEl[d]])
                        oriZ[resFMeSetEl[d]]=oriZ[resFMeSetEl[d]]+(resFMe[d][12]/Atot[resFMeSetEl[d]])
                        d = d + 1
#                    print Atot,oriX,oriY,oriZ
                    if ni == 1:
                        ori0xE=oriX
                        ori0yE=oriY
                        ori0zE=oriZ
                        Np0E=Np
                        Yp0E=Yp
                        Zp0E=Zp
                        AT0E=AT
                        angleX=0.0
                        angleY=0.0
                        angleZ=0.0
                    else:
                        angleX=Sprod(Np0E[0],Np0E[1],Np0E[2],Np[0],Np[1],Np[2])
                        angleY=Sprod(Yp0E[0],Yp0E[1],Yp0E[2],Yp[0],Yp[1],Yp[2])
                        angleZ=Sprod(Zp0E[0],Zp0E[1],Zp0E[2],Zp[0],Zp[1],Zp[2])
                    for w in range(0,len(resFMe)):
                        L1=getLocalVec(resFMe[w][0],resFMe[w][1],resFMe[w][2],oriX[resFMeSetEl[w]],oriY[resFMeSetEl[w]],oriZ[resFMeSetEl[w]],AT)
                        Iyy[resFMeSetEl[w]]=Iyy[resFMeSetEl[w]]+L1[2]*L1[2]*resFMe[w][9]
                        Izz[resFMeSetEl[w]]=Izz[resFMeSetEl[w]]+L1[1]*L1[1]*resFMe[w][9]
                        Iyz[resFMeSetEl[w]]=Iyz[resFMeSetEl[w]]+L1[1]*L1[2]*resFMe[w][9]
                        sumMxE[resFMeSetEl[w]]=sumMxE[resFMeSetEl[w]]+(resFMe[w][7]*L1[2])-(resFMe[w][8]*L1[1])
                        sumMyE[resFMeSetEl[w]]=sumMyE[resFMeSetEl[w]]+resFMe[w][6]*L1[2]
                        sumMzE[resFMeSetEl[w]]=sumMzE[resFMeSetEl[w]]+resFMe[w][6]*L1[1]
                        if comb == 0 :
                           sumFxE[resFMeSetEl[w]]=sumFxE[resFMeSetEl[w]]+resFMe[w][6]
                           sumFyE[resFMeSetEl[w]]=sumFyE[resFMeSetEl[w]]+resFMe[w][7]
                           sumFzE[resFMeSetEl[w]]=sumFzE[resFMeSetEl[w]]+resFMe[w][8]
                    if comb==1 :
                        for w in range(0,len(resFMec)):
                           sumFxE[resFMeSetElc[w]]=sumFxE[resFMeSetElc[w]]+resFMec[w][6]
                           sumFyE[resFMeSetElc[w]]=sumFyE[resFMeSetElc[w]]+resFMec[w][7]
                           sumFzE[resFMeSetElc[w]]=sumFzE[resFMeSetElc[w]]+resFMec[w][8]
                    if foundPlst == 1 :
                        for w in range(0,len(resPlstE)):
                            sumPlAE[resFMeSetEl[w]]=sumPlAE[resFMeSetEl[w]]+(resPlstE[w]/Atot[resFMeSetEl[w]])
                        for w in range(0,len(maxPlstE)):
                            maxPlstAE[resFMeSetEl[w]]=max(maxPlstAE[resFMeSetEl[w]],(maxPlstE[w]))
                    if txt==1:
                        if ni==1:
                            tmp_csvE.write('inc,time')
                            for snam in elset1 :
                                tmp_csvE.write(',Element Set,Origo X GCS,Origo Y GCS,Origo Z GCS,Origo X LCS 0,Origo Y LCS 0,Origo Z LCS 0,Angle X rel LCS 0,Angle Y rel LCS 0,Angle Z LCS rel 0,Iyy,Izz,Iyz,Area,Fx,Fy,Fz,Mx,My,Mz,Av Pl Strain,Max Pl Strain')  
                            tmp_csvE.write('\n') 
                        else :
                            tmp_csvE.write('\n') 
                        tmpe='%5i,%4.4f' %(ni-1,atime)
                        tmp_csvE.write(tmpe)                        
                    for snam in elset1 :
                        cogL=getLocalVec(oriX[snam],oriY[snam],oriZ[snam],ori0xE[snam],ori0yE[snam],ori0zE[snam],AT0E)
#                        cogL2=getLocalVec(ori2X,ori2Y,ori2Z,ori0xE,ori0yE,ori0zE,AT0E)
                        tmp='                                                        SET NAME = %s  \n\n' %(snam)
                        tmp_txtE.write(tmp)
                        tmp_txtE.write('                                                        Cross Section Center of Gravity \n')
                        tmp_txtE.write('              ------------------------------------------------------------------------------------------------------------------------\n')
                        tmp_txtE.write('                           Global CS       Rel. Incr 0 Local CS    Angle rel Incr 0 LCS\n')
                        tmp='              Origin X      %10.3f          %10.3f              %10.3f (Angle between local X-axis and Incr 0 local X-axis)\n' %(oriX[snam],cogL[0],angleX)
                        tmp_txtE.write(tmp)
                        tmp='              Origin Y      %10.3f          %10.3f              %10.3f (Angle between local Y-axis and Incr 0 local Y-axis)\n' %(oriY[snam],cogL[1],angleY)
                        tmp_txtE.write(tmp)
                        tmp='              Origin Z      %10.3f          %10.3f              %10.3f (Angle between local Z-axis and Incr 0 local Z-axis)\n\n' %(oriZ[snam],cogL[2],angleZ)
                        tmp_txtE.write(tmp)
                        tmp='              Area          %10.3f\n\n' %(Atot[snam])
                        tmp_txtE.write(tmp)
                        tmp='              Moment of Inertia Iyy  %10.3f\n' %(Iyy[snam])
                        tmp_txtE.write(tmp)
                        tmp='              Moment of Inertia Izz  %10.3f\n' %(Izz[snam])
                        tmp_txtE.write(tmp)
                        tmp='              Moment of Inertia Iyz  %10.3f\n\n' %(Iyz[snam])
                        tmp_txtE.write(tmp)
                        tmp_txtE.write('                                                    Cross Section Resulting Forces and Moments\n')
                        tmp_txtE.write('------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                        tmp_txtE.write('Element Set Name                          Area        Force X        Force Y        Force Z       Moment X       Moment Y       Moment Z   Av Pl Strain  Max Pl Strain \n')
                        tmp='%-30s      %10.3f %14.2f %14.2f %14.2f %14.2f %14.2f %14.2f %14.7e %14.7e\n\n' %(snam,Atot[snam],sumFxE[snam],sumFyE[snam],sumFzE[snam],sumMxE[snam],sumMyE[snam],sumMzE[snam],sumPlAE[snam],maxPlstAE[snam]) 
                        tmp_txtE.write(tmp)
                        tmp_txtE.write('------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n')
                        if txt==1:
                            if ni==1:
                                tmp_csvS.write('inc,time,Element,Cut Point Number,X,Y,Z,Sx,Sy,Sz\n')    
                            tmpe=',%s,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e' %(snam,oriX[snam],oriY[snam],oriZ[snam],cogL[0],cogL[1],cogL[2],angleX,angleY,angleZ,Iyy[snam],Izz[snam],Iyz[snam],Atot[snam],sumFxE[snam],sumFyE[snam],sumFzE[snam],sumMxE[snam],sumMyE[snam],sumMzE[snam],sumPlAE[snam],maxPlstAE[snam])
                            tmp_csvE.write(tmpe)
                            for ww in range(0,len(resEL)):
                                tmps='%5i,%4.4f,%i,%i,%14.6f,%14.6f,%14.6f,%14.6f,%14.6f,%14.6f\n' %(ni-1,atime,resEL[ww][0],resEL[ww][1],resFMnp[ww][0],resFMnp[ww][1],resFMnp[ww][2],resFMns[ww][0],resFMns[ww][1],resFMns[ww][2])
                                tmp_csvS.write(tmps)
                work_cnt=work_cnt+1
                print('Total work done: ',str(round(((work_cnt*100.0)/tot_inc),1)),'%')
                if bbTol != -1.0 and (ni==1 or irem==1):
                    if cntClEl==0:
                        procOut.write('# | End of List\n')
                    else :
                        procOut.write('\n# | End of List\n')
                curSec=curSec+1
            ftim=time.clock()
            etim=ftim-stim
    
            if nas == 1 :
                tmp_nas.write('MAT1           1 210000.             0.3\n')
#                tmp_nas.close()
            if pr_tenel==1:
                tmp_ten.write('MAT1           1 210000.             0.3\n')
#                tmp_ten.close()
#            if cbpr != 'none':
#                tmp_txt.close()
#                tmp_csv.close()
#            if elset1 != '':
#                tmp_txtE.close()
#                tmp_csvE.close()
#                tmp_csvS.close()

            ni=ni+1
        p.close()
        if bbTol != -1.0  :
            procOut.write('*set_update on\n')
            procOut.close()
        ftim2=time.clock()
        etim2=ftim2-stim_1
        print('Total time: ',etim2)
# ----------------------------------------------
# End 3D case, following axisymmetric...
#
    else : # axisymmetric model
        if boxTol == '':
           bbTol=3.0 # more generous tolerance for axi...
        while curSec < numSec:
            boxNodDic={}    
            print('Processing section number ' + str(curSec+1))
            stim=time.clock()
            time_inc=0.0
            time_disp=0.0
            time_el=0.0
            pn1=n1g[curSec]
            boxNodDic[pn1]=1
            pn1st=str(pn1)
            if bbTol != -1.0 and curSec==0:
                procOut.write('*store_elements\n')
                procOut.write('crset%s_0\n' %pn1st)
            comb=0
            if nas==1 :
                nas_out=t16F[:-4] + '_' + pn1st + '.bdf'
                tmp_nas=open(nas_out,'w')
# write stuff to out file(s)
            if cbpr != 'none':
                txt_out=t16F[:-4] + '_' + pn1st + 'cb.res'
                csv_out=t16F[:-4] + '_' + pn1st + 'cb.csv'
                tmp_txt=open(txt_out,'w')
                tmp_csv=open(csv_out,'w')
# write to out file for contact body
            if cbpr != 'none':
                tmp_txt.write('Cross Section Postprocessing tool output file\n')
                tmp_txt.write('---------------------------------------------\n')
                tmp_txt.write('Marc input file used       : ' + datF + '\n')
                tmp_txt.write('Marc results file used     : ' + t16F + '\n')
                tmp_txt.write('Node number 1              : ' + str(pn1) + '\n')
                if int(si) == 1:
                    tmp_txt.write('Units used                 : SI\n')
                else:
                    tmp_txt.write('Units used                 : mm\n')
                if cbpr == '':
                    tmp_txt.write('Contact bodies in summation: All\n')
                else:
                    tmp_txt.write('Contact bodies in summation: ' + cbpr + '\n')
                if int(sym) == 1:
                    tmp_txt.write('Symmetri option used       : Yes\n')
                else:
                    tmp_txt.write('Symmetri option used       : No\n')
                if int(stressT) == 0:
                    stressName='Cauchy Stress'
                    stressName2='Stress'
                    tmp_txt.write('Stress Tensor used         : Cauchy Stress\n')
                else:
                    stressName='Stress'
                    stressName2='Cauchy Stress'
                    tmp_txt.write('Stress Tensor used         : Stress\n')
                tmp_txt.write('Refinement level used      : ' + str(nas) + '\n')
                if comb==1:
                    tmp_txt.write('Extrapolation method used  : combined\n')
                else:
                    tmp_txt.write('Extrapolation method used  : ' + p.extrapolate + '\n')
                if int(calcPl) == 1:
                    tmp_txt.write('Plastic Strain calculation : Yes\n')
                else:
                    tmp_txt.write('Plastic Strain calculation : No\n')
                tmp_txt.write('Number of elements         : ' + str(numElem) + '\n')
                tmp_txt.write('Number of nodes            : ' + str(numNodes) + '\n')
                tmp_txt.write('Number of element tensors  : ' + str(numET) + '\n')
                tmp_txt.write('Number of node scalars     : ' + str(numNS) + '\n')
                tmp_txt.write('Number of increments       : ' + str(ninc-1) + '\n')
            ni=1
            gid=1
            boxElem=[]
            boxElDic={}
            test_list=[]
#            ninc=2
            for ni in range(1, ninc, nskip):
                dicxd={}
                stim1=time.clock()
                atime=float(p.time)
                p.moveto(ni)
                irem=0
                print('Scanning increment number ',str(ni-1))
                numNodesX = p.nodes()
                numElemX = p.elements()
    #            print numNodesX,numElemX
                if (numNodesX != numNodes or numElemX != numElem ) :
                    print('Remeshing assumed')
                    print('New number of nodes: ',numNodesX,', New number of elements: ',numElemX)
                    if cntClEl==0:
                        procOut.write('# | End of List\n')
                    else :
                        procOut.write('\n# | End of List\n')
                    procOut.write('*store_elements\n')
                    procOut.write('crset%s_%i\n' %(pn1st,(ni-1)))
                    foundC=0
                    boxNodDic={}    
                    boxNodDic[pn1]=1
                    boxElem=[]
                    boxElDic={}
                    dicx={}
                    dicC={}
                    dicS={}
                    dicEl={}
                    CBNR=[]
                    CBNM=[]
                    irem=1
                    numNodes=numNodesX
                    numElem=numElemX
                    ki=0
                    while ki < numNodes :
                        nod22 = p.node(ki)
                        n22x=nod22.x
                        n22y=nod22.y
                        n22z=nod22.z
                        dicx[nod22.id]=[n22x,n22y,n22z]
                        if immRem == 0 :
                            if n22x > maxX :
                                maxX=n22x
                            if n22x < minX :
                                minX=n22x
                            if n22y > maxY :
                                maxY=n22y
                            if n22y < minY :
                                minY=n22y
                            if n22z > maxZ :
                                maxZ=n22z
                            if n22z < minZ :
                                minZ=n22z
            
                        ki=ki+1
                    ki=0
                    while ki < numElem :
                       currEl=p.element_id(ki)
                       el22 = p.element(ki)
                       numNe=el22.len
                       eltype=el22.type
                       kii=0
                       tmp=str(currEl) + ' ' + str(eltype)
                       while kii < numNe :
                          tmp=tmp  + ' ' +  str(el22.items[kii])
                          kii=kii+1
                       dicEl[currEl]=tmp
                       ki=ki+1
                    ncbody=p.cbodies()
                    ki=0
                    while ki < ncbody :
                        curCB=p.cbody(ki)
                        if curCB.type == 0 :
                           foundC=1
                           cnid=curCB.id
                           cnam=curCB.name
                           CBNR.append(cnid)
                           CBNM.append(cnam)
                           kii=0
                           numElC=curCB.nelements
                           while kii < numElC :
                               dicC[curCB.elements[kii]]=[cnid,cnam]
                               kii=kii+1
                        ki=ki+1
                    if foundC == 0:
                        CBNR.append(0)
                        CBNM.append('None')
                stim2=time.clock()
                time_inc=time_inc+stim2-stim1
                atime=float(p.time)
                if ni == 1 :
                    kk= 0
                    while kk < numNodes:
                        NodId=p.node_id(kk)
                        dicxd[NodId]=[0.0,0.0,0.0]
                        kk=kk+1
                else :
                    if bbTol == -1.0 or irem == 1:
                        kk= 0
                        while kk < numNodes:
                            NodId=p.node_id(kk)
                            xx1,yy1,zz1=p.node_displacement(kk)
                            dicxd[NodId] = [xx1,yy1,zz1]
                            kk=kk+1
                    else:
                        for NodId in boxNodDic :
                            kk=p.node_sequence(NodId)
                            xx1,yy1,zz1=p.node_displacement(kk)
                            dicxd[NodId] = [xx1,yy1,zz1]
#           
                stim3=time.clock()
                time_disp=time_disp+stim3-stim2
                p1=[]
                p2=[]
                p3=[]
                AT=[]
    
                p1.append(dicx[pn1][0]+dicxd[pn1][0])
                p1.append(dicx[pn1][1]+dicxd[pn1][1])
                p1.append(dicx[pn1][2]+dicxd[pn1][2])
                p2.append(dicx[pn1][0]+dicxd[pn1][0])
                p2.append(dicx[pn1][1]+dicxd[pn1][1]+1.0)
                p2.append(0.0)
                p3.append(dicx[pn1][0]+dicxd[pn1][0])
                p3.append(dicx[pn1][1]+dicxd[pn1][1])
                p3.append(1.0)
                Np=Xprod(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2],p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2])
                offset=0.1*sca
                uy=0
                while uy < 3:
                    p1[uy]=p1[uy]+offset*Np[uy]
                    p2[uy]=p2[uy]+offset*Np[uy]
                    p3[uy]=p3[uy]+offset*Np[uy]
                    uy=uy+1
                Origo=[p1[0],0.0,0.0]
                Yp=[0.0,1.0,0.0]
                Zp=[0.0,0.0,1.0]
                if cbpr != 'none':
                    tmp_txt.write('----------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                    tmp='-----------------------------------------    CROSS SECTION RESULTS FOR INCREMENT %4i     ----------------------------------------------------------------------\n\n' %((ni-1))
                    tmp_txt.write(tmp)
                    tmp='Analysis Time = %4.3f\n\n' %(atime)
                    tmp_txt.write(tmp)
                    tmp_txt.write('     Cross Section Geometrical Properties \n')
                    tmp_txt.write('----------------------------------------------\n')
                    tmp='Local X axis: %10.7f %10.7f %10.7f\n' %(Np[0],Np[1],Np[2])
                    tmp_txt.write(tmp)    
                    tmp='Local Y axis: %10.7f %10.7f %10.7f\n' %(Yp[0],Yp[1],Yp[2])
                    tmp_txt.write(tmp)
                    tmp='Local Z axis: %10.7f %10.7f %10.7f\n\n' %(Zp[0],Zp[1],Zp[2])
                    tmp_txt.write(tmp)
                    tmp_txt.write('                            Cross Section Center of Gravity \n')
                    tmp_txt.write('------------------------------------------------------------------------------\n')
                    tmp_txt.write('                Global CS       Rel. Incr 0 Local CS    Mid point of N1 and N2\n')
                AT.append(Np[0])
                AT.append(Np[1])
                AT.append(Np[2])
                AT.append(Yp[0])
                AT.append(Yp[1])
                AT.append(Yp[2])
                AT.append(Zp[0])
                AT.append(Zp[1])
                AT.append(Zp[2])
                resFM=[]
                resCB=[]
                resCBPl2=[]
                resCBPl=[]
                foundStress=0
                foundPlst=0
                resPlst=[]
                maxPlst=[]
                maxPlstAll=0.0
                maxPlstEl=0
                maxPlstNode=0
                if ni == 1 or irem == 1 :
                    cntClEl=0
                for elem_id in dicEl:
                    spl=dicEl[elem_id].split('\n')[0].split()
                    q=0
                    elem_type=int(spl[1])
                    if ni == 1  or irem == 1 or bbTol== -1.0:
                        boxElDic[elem_id]=1
                    el_ind=p.element_sequence(elem_id)
                    for ii in range(0,numES):
                        s_sca = p.element_scalar_label(ii)
                        if s_sca=='Total Equivalent Plastic Strain' and calcPl==1:
                            slist = p.element_scalar(el_ind,ii)
                            for ks in range(0,len(slist)):
                                if slist[ks].value > maxPlstAll :
                                    maxPlstAll=slist[ks].value
                                    maxPlstEl=elem_id
                                    maxPlstNode=slist[ks].id
                            foundPlst=1
                    for i in range(0,numET):
                        s_ten = p.element_tensor_label(i)
                        if s_ten == stressName :
                            if comb==1:
                                p.extrapolation('linear')
                                tlist = p.element_tensor(el_ind,i)
                                p.extrapolation('average')
                                tlist2=p.element_tensor(el_ind,i)
                            else:
                                tlist = p.element_tensor(el_ind,i)
                                tlist2 = tlist
                            foundStress=1
                    if foundStress== 0 :
                        print('Requested Stress Tensor not found in t16 file')
                        print('Trying to find alternative tensor ',stressName2)
                        stressName=stressName2
                        for i in range(0,numET):
                            s_ten = p.element_tensor_label(i)
                            if s_ten == stressName :
                                if comb==1:
                                    p.extrapolation('linear')
                                    tlist = p.element_tensor(el_ind,i)
                                    p.extrapolation('average')
                                    tlist2=p.element_tensor(el_ind,i)
                                else:
                                    tlist = p.element_tensor(el_ind,i)
                                    tlist2 = tlist
                                foundStress=1
                    if foundStress== 0 :
                        print('Requested Stress Tensor not found in t16 file, quiting...')
                        return 1
                    crossp=[]
                    crossS=[]
                    crossPlst=[]
                    crList=[]
                    if (elem_type == 10 or elem_type==2) and (boxElDic[elem_id]==1 or bbTol == -1.0 ):
                        if foundC == 0:
                            cbnr=0
                        else:
                            try:
                                cbnr=dicC[elem_id][0]
                            except:
                                print('element not in any contact body')
                                continue
                        cbnrP=ni*100+cbnr
                        n1=[int(spl[2]),dicx[int(spl[2])][0]+dicxd[int(spl[2])][0],dicx[int(spl[2])][1]+dicxd[int(spl[2])][1],dicx[int(spl[2])][2]+dicxd[int(spl[2])][2]]
                        n2=[int(spl[3]),dicx[int(spl[3])][0]+dicxd[int(spl[3])][0],dicx[int(spl[3])][1]+dicxd[int(spl[3])][1],dicx[int(spl[3])][2]+dicxd[int(spl[3])][2]]
                        n3=[int(spl[4]),dicx[int(spl[4])][0]+dicxd[int(spl[4])][0],dicx[int(spl[4])][1]+dicxd[int(spl[4])][1],dicx[int(spl[4])][2]+dicxd[int(spl[4])][2]]
                        k1=[n1,n2,(n2[1]-n1[1]),(n2[2]-n1[2]),(n2[3]-n1[3])]
                        k2=[n2,n3,(n3[1]-n2[1]),(n3[2]-n2[2]),(n3[3]-n2[3])]
                        if (elem_type == 10) :
                            n4=[int(spl[5]),dicx[int(spl[5])][0]+dicxd[int(spl[5])][0],dicx[int(spl[5])][1]+dicxd[int(spl[5])][1],dicx[int(spl[5])][2]+dicxd[int(spl[5])][2]]
                            k3=[n3,n4,(n4[1]-n3[1]),(n4[2]-n3[2]),(n4[3]-n3[3])]
                            k4=[n4,n1,(n1[1]-n4[1]),(n1[2]-n4[2]),(n1[3]-n4[3])]
                        elif (elem_type == 2) :
                            k3=[n3,n1,(n1[1]-n3[1]),(n1[2]-n3[2]),(n1[3]-n3[3])]  
                        else :
                            print('unsupported element type...')                        
                        cr1=cross(absBB,Origo,Np,k1)
                        if cr1[0] != 'p' :
                            if cr1[0] != 'n':
                                crossp.append(cr1)
                                crList.append(1.0e-6)
                                crossS.append(get_stress(n1[0],n2[0],cr1[0],tlist))
                                if foundPlst==1:
                                    crossPlst.append(get_plst(n1[0],n2[0],cr1[0],slist))
                            else:
                                crList.append(cr1[1])
                        cr2=cross(absBB,Origo,Np,k2)
                        dub=0
                        if cr2[0] != 'p' :
                            if cr2[0] != 'n':
                                crList.append(1.0e-6)
                                for i in range (0,len(crossp)):
                                    if abs(cr2[1]-crossp[i][1]) < dubtol and abs(cr2[2]-crossp[i][2]) < dubtol and abs(cr2[3]-crossp[i][3]) < dubtol :
                                        dub=1
                                if dub==0:
                                    crossp.append(cr2)
                                    crossS.append(get_stress(n2[0],n3[0],cr2[0],tlist))
                                    if foundPlst==1:
                                        crossPlst.append(get_plst(n2[0],n3[0],cr2[0],slist))
                            else:
                                crList.append(cr2[1])
                        if (elem_type == 10) :
                            cr3=cross(absBB,Origo,Np,k3)
                            dub=0
                            if cr3[0] != 'p' :
                                if cr3[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr3[1]-crossp[i][1]) < dubtol and abs(cr3[2]-crossp[i][2]) < dubtol and abs(cr3[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr3)
                                        crossS.append(get_stress(n3[0],n4[0],cr3[0],tlist))
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n3[0],n4[0],cr3[0],slist))
                                else:
                                    crList.append(cr3[1])
                            cr4=cross(absBB,Origo,Np,k4)
                            dub=0
                            if cr4[0] != 'p' :
                                if cr4[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr4[1]-crossp[i][1]) < dubtol and abs(cr4[2]-crossp[i][2]) < dubtol and abs(cr4[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr4)
                                        crossS.append(get_stress(n4[0],n1[0],cr4[0],tlist))
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n4[0],n1[0],cr4[0],slist))
                                else:
                                    crList.append(cr4[1])
                        elif (elem_type == 2) :
                            cr3=cross(absBB,Origo,Np,k3)
                            dub=0
                            if cr3[0] != 'p' :
                                if cr3[0] != 'n':
                                    crList.append(1.0e-6)
                                    for i in range (0,len(crossp)):
                                        if abs(cr3[1]-crossp[i][1]) < dubtol and abs(cr3[2]-crossp[i][2]) < dubtol and abs(cr3[3]-crossp[i][3]) < dubtol :
                                            dub=1
                                    if dub==0:
                                        crossp.append(cr3)
                                        crossS.append(get_stress(n3[0],n1[0],cr3[0],tlist))
                                        if foundPlst==1:
                                            crossPlst.append(get_plst(n3[0],n1[0],cr3[0],slist))
                                else:
                                    crList.append(cr3[1])                        
#  if min(t) less than a tolerance, save this element...
#                        if elem_id == 410361 :
#                           print crList,bbTol
                        if (ni == 1 or irem==1) and bbTol != -1.0:
                            if min([sqrt(i**2) for i in crList]) < bbTol :
                                boxElDic[elem_id]=1
                                boxNodDic[n1[0]]=1
                                boxNodDic[n2[0]]=1
                                boxNodDic[n3[0]]=1
                                cntClEl=cntClEl+1
                                procOut.write(' %i' %elem_id)
                                if cntClEl == 10 :
                                    procOut.write('\n')
                                    cntClEl=0
                                if elem_type == 10:
                                    boxNodDic[n4[0]]=1
                            else :
                                boxElDic[elem_id]=0
# If number of unique crossing points is  2, a line patch is found
                    if len(crossp) == 2:
                        center=meanP(crossp,crossS)
                        cro0=[0,center[0],center[1],center[2]]
                        cro0s=[center[3],center[4],center[5],center[6],center[7],center[8]]                           
                        FM1=getForceMomAxi(crossp[0][2],crossp[1][2],cro0[2],crossS[0],crossS[1],cro0s)
                        resFM.append(FM1)
                        resCB.append(cbnr)
                        if foundPlst==1:
                           center_pl=meanP_pl(crossPlst)
                           maxPlst.append(max(crossPlst))
                           resCBPl2.append(cbnr)
                           resCBPl.append(cbnr)
                           pl1=getPlStA(FM1[3],center_pl,center_pl,center_pl)
                           resPlst.append(pl1)
                        if nas == 1:
                            pgrid(gid,center[0],center[1],center[2],tmp_nas)
                            pgrid(gid+1,crossp[0][1],crossp[0][2],crossp[0][3],tmp_nas)
                            pgrid(gid+2,crossp[1][1],crossp[1][2],crossp[1][3],tmp_nas)
                            pbeam(gid,cbnrP,gid,gid+1,tmp_nas)
                            pbeam(gid+1,cbnrP,gid,gid+2,tmp_nas)
                            pPloadB(10,FM1,gid,tmp_nas) 
                            gid=gid+3
                if nas == 1:
                    for ut in range(0,len(CBNR)):
                        tmp='PBARL   %8i       1             ROD\n             1.0\n' %(CBNR[ut]+100*ni)
                        tmp_nas.write(tmp)
                if cbpr != 'none':
                    if txt==1 and ni==1:
                        tmp_csv.write('inc,time,Origo X GCS,Origo Y GCS,Origo Z GCS,Iyy,Izz,Iyz,Contact Body,Area,Fx,Fy,Fz,Av Pl Strain,Max Pl Strain,Max Pl Strain All,@Element,@Node,')    
                        if cbpr == '':
                            for ub in range(0,len(CBNR)):
                                tmp_csv.write('Contact Body,Area,Fx,Fy,Fz,Av Pl Strain,Max Pl Strain,')
                        else:
                            splcb=cbpr.split(',')
                            for dcb in splcb:
                                tmp_csv.write('Contact Body,Area,Fx,Fy,Fz,Av Pl Strain,Max Pl Strain,')
                        tmp_csv.write('\n')
                    sumFx=0.
                    sumFy=0.
                    sumFz=0.
                    sumMx=0.
                    sumMy=0.
                    sumMz=0.
                    sumPlA=0.0
                    maxPlstA=0.0
                    Atot=sum([l[3] for l in resFM])
                    oriX=Origo[0]
                    oriY=Origo[1]
                    oriZ=Origo[2]
                    Iyy=0
                    Izz=0
                    Ixx=0
                    dicCBF={}
                    dicCBM={}
                    dicCBA={}
                    dicCBPl={}
                    dicCBPl2={}
                    for dcb in range(0,len(CBNR)):
                        dicCBF[CBNR[dcb]]=[0.0,0.0,0.0]
                        dicCBM[CBNR[dcb]]=[0.0,0.0,0.0]
                        dicCBA[CBNR[dcb]]=[0.0]
                        dicCBPl[CBNR[dcb]]=[0.0]
                        dicCBPl2[CBNR[dcb]]=0.0
                    for w in range(0,len(resFM)):
                        Iyy=Iyy+resFM[w][5]
                        Izz=Izz+resFM[w][6]
                        Ixx=Ixx+resFM[w][4]
                        sumFx=sumFx+resFM[w][0]
                        sumFy=sumFy+resFM[w][1]
                        sumFz=sumFz+resFM[w][2]
                        dicCBF[resCB[w]]=[dicCBF[resCB[w]][0]+resFM[w][0],dicCBF[resCB[w]][1]+resFM[w][1],dicCBF[resCB[w]][2]+resFM[w][2]]
                        dicCBA[resCB[w]]=[dicCBA[resCB[w]][0]+(resFM[w][3])]
                    if foundPlst == 1:
                        for w in range(0,len(resPlst)):
                            sumPlA=sumPlA+resPlst[w]
                            dicCBPl[resCBPl[w]]=[dicCBPl[resCBPl[w]][0]+(resPlst[w])]
                        for w in range(0,len(maxPlst)):
                            dicCBPl2[resCBPl2[w]]=max(dicCBPl2[resCBPl2[w]],(maxPlst[w]))
                        maxPlstA=0.0
#                        print maxPlstAll,maxPlstEl,maxPlstNode
                        for key in dicCBPl2:
                            if dicCBPl2[key] > maxPlstA :
                                maxPlstA=dicCBPl2[key]
                    tmp='Origin X      %10.3f          \n' %(oriX)
                    tmp_txt.write(tmp)
                    tmp='Origin Y      %10.3f          \n' %(oriY)
                    tmp_txt.write(tmp)
                    tmp='Origin Z      %10.3f          \n' %(oriZ)
                    tmp_txt.write(tmp)
                    tmp='Area          %10.3f\n' %(Atot)
                    tmp_txt.write(tmp)
                    tmp='Moment of Inertia Iyy  %14.3f\n' %(Iyy)
                    tmp_txt.write(tmp)
                    tmp='Moment of Inertia Izz  %14.3f\n' %(Izz)
                    tmp_txt.write(tmp)
                    tmp='Moment of Inertia Ixx  %14.3f\n\n' %(Ixx)
                    tmp_txt.write(tmp)
                    tmp_txt.write('                          Cross Section Resulting Forces and Moments (cylindrical directions)\n')
                    tmp_txt.write('----------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
                    tmp_txt.write('Body Number  Body Name                  Area        Force X        Force Y        Force Z        Av Pl Strain   Max Pl Strain   Max Pl Str All @Element    @Node\n')
                    tmp_txt.write('                                                    (axial)       (radial)   (circumferential)   \n')
                    if cbpr == '':
                        tmp='Sum of all   Sum of all           %10.3f %14.2f %14.2f %14.2f      %14.8e  %14.8e   %14.8e %8i %8i\n' %(Atot,sumFx,sumFy,sumFz,sumPlA/Atot,maxPlstA,maxPlstAll,maxPlstEl,maxPlstNode) 
                        tmp_txt.write(tmp)
                        if txt==1:
                            tmpc='%5i,%4.4f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,  SUM,%10.3f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e,%14.8e,%8i,%8i,' %(ni-1,atime,oriX,oriY,oriZ,Iyy,Izz,Ixx,Atot,sumFx,sumFy,sumFz,sumPlA/Atot,maxPlstA,maxPlstAll,maxPlstEl,maxPlstNode)
                            tmp_csv.write(tmpc)
                        for dcb in range(0,len(CBNR)):
                            try :
                               plStDivA=dicCBPl[CBNR[dcb]][0]/dicCBA[CBNR[dcb]][0]
                            except :
                               plStDivA=0.0
                            tmp='%10i   %-20s %10.3f %14.2f %14.2f %14.2f      %14.8e  %14.8e \n' %(CBNR[dcb],CBNM[dcb],dicCBA[CBNR[dcb]][0],dicCBF[CBNR[dcb]][0],dicCBF[CBNR[dcb]][1],dicCBF[CBNR[dcb]][2],plStDivA,dicCBPl2[CBNR[dcb]])
                            tmp_txt.write(tmp)
                            if txt==1:
                                tmpc='%20s,%4.4f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e,' %(CBNM[dcb],dicCBA[CBNR[dcb]][0],dicCBF[CBNR[dcb]][0],dicCBF[CBNR[dcb]][1],dicCBF[CBNR[dcb]][2],plStDivA,dicCBPl2[CBNR[dcb]])
                                tmp_csv.write(tmpc)
                        tmp_csv.write('\n')
                        if txt==1:
                            tmp_txt.write('\n')
                    else:
                        sumFx=0.
                        sumFy=0.
                        sumFz=0.
                        sumMx=0.
                        sumMy=0.
                        sumMz=0.
                        sumPlA=0.0
                        maxPlstA=0.0
                        splcb=cbpr.split(',')
                        for sumCB in splcb:
                            sumFx=sumFx+dicCBF[int(sumCB)][0]
                            sumFy=sumFy+dicCBF[int(sumCB)][1]
                            sumFz=sumFz+dicCBF[int(sumCB)][2]
                            if foundPlst == 1 :
                               sumPlA=sumPlA+dicCBPl[int(sumCB)][0]
                               if dicCBPl2[int(sumCB)] > maxPlstA :
                                      maxPlstA=dicCBPl2[int(sumCB)]
                        tmp='Sum selected Sum selected         %10.3f %14.2f %14.2f %14.2f      %14.8e  %14.8e \n' %(Atot,sumFx,sumFy,sumFz,sumPlA/Atot,maxPlstA) 
                        tmp_txt.write(tmp)
                        if txt==1:
                            tmpc='%5i,%4.4f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,  SUM,%10.3f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e,' %(ni-1,atime,oriX,oriY,oriZ,Iyy,Izz,Ixx,Atot,sumFx,sumFy,sumFz,sumPlA/Atot,maxPlstA)
                            tmp_csv.write(tmpc)
                        for dcb in splcb:
#                            dcb=int(dcb)-1
                            dcb=int(dcb)
                            try :
                               plStDivA=dicCBPl[CBNR[dcb]][0]/dicCBA[CBNR[dcb]][0]
                            except :
                               plStDivA=0.0
                            tmp='%10i   %-20s %10.3f %14.2f %14.2f %14.2f      %14.8e  %14.8e \n' %(CBNR[dcb],CBNM[dcb],dicCBA[CBNR[dcb]][0],dicCBF[CBNR[dcb]][0],dicCBF[CBNR[dcb]][1],dicCBF[CBNR[dcb]][2],plStDivA,dicCBPl2[CBNR[dcb]])
                            tmp_txt.write(tmp)
                            if txt==1:
                                tmpc='%20s,%4.4f,%14.2f,%14.2f,%14.2f,%14.8e,%14.8e,' %(CBNM[dcb],dicCBA[CBNR[dcb]][0],dicCBF[CBNR[dcb]][0],dicCBF[CBNR[dcb]][1],dicCBF[CBNR[dcb]][2],plStDivA,dicCBPl2[CBNR[dcb]])
                                tmp_csv.write(tmpc)
                        tmp_csv.write('\n')
                        if txt==1:
                            tmp_txt.write('\n')

                        
                # ni=ni+1 # TAG commented this line
            if nas == 1 :
                tmp_nas.write('MAT1           1 210000.             0.3\n')
                tmp_nas.close()
            if bbTol != -1.0 :
                if cntClEl==0:
                    procOut.write('# | End of List\n')
                else :
                    procOut.write('\n# | End of List\n')
            curSec=curSec+1    
        p.close()
        if bbTol != -1.0 :
            procOut.write('*set_update on\n')
            procOut.close()
    
    return 0