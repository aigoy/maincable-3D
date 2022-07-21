from math import sqrt
import numpy as np
import pandas as pd
#from scipy.interpolate import lagrange

# 此处输入模型信息文件的路径
nodepath='C:\\Users\\chaizhimin\\Documents\\主缆找形\\maincable-python\\test3D\\NODE.txt' 
elempath='C:\\Users\\chaizhimin\\Documents\\主缆找形\\maincable-python\\test3D\\ELEM.txt'
hangerpath='C:\\Users\\chaizhimin\\Documents\\主缆找形\\maincable-python\\test3D\\HANGER.txt'
infopath='C:\\Users\\chaizhimin\\Documents\\主缆找形\\maincable-python\\test3D\\info.txt'
fixednodepath='C:\\Users\\chaizhimin\\Documents\\主缆找形\\maincable-python\\test3D\\FIXEDNODE.txt'
outputpath='C:\\Users\\chaizhimin\\Documents\\主缆找形\\maincable-python\\test3D\\result.xlsx'

#读取节点坐标信息,fname为坐标信息txt文件，内部格式为3/4列矩阵，依次存放节点编号、x坐标、y坐标（、z坐标），n表示2维/3维
def loadnodes(fname):
    A=np.loadtxt(fname)
    # n=A.shape[1]
    nodeid=A[:,0]
    # if n==3:
    #     B=pd.DataFrame(A,columns=['nodeid','x','y'],index=nodeid)
    # elif n==4:
    B=pd.DataFrame(A,columns=['nodeid','x','y','z'],index=nodeid)
    B['nodeid']=pd.to_numeric(B['nodeid'],downcast='integer')
    B.index=B['nodeid']
    return(B)

#读取吊杆信息,fname为txt文件，内部格式为3/4列矩阵，依次存放节点编号、吊杆y方向分力、y坐标（、z坐标），n表示2维/3维
#此处节点编号应该与nodeid一一对应
def loadhanger(fname):
    A=np.loadtxt(fname)
    # n=A.shape[1]
    hanid=A[:,0]
    # if n==3:
    #     B=pd.DataFrame(A,columns=['hanid','Ty','y'],index=hanid)
    # elif n==4:
    B=pd.DataFrame(A,columns=['hanid','Ty','y','z'],index=hanid)
    return(B)

#读取单元信息，fname为txt文件，内部格式为3列矩阵，存放单元编号与对应节点编号。
def loadelem(fname):
    A=np.loadtxt(fname)
    elemid=A[:,0]
    B=pd.DataFrame(A,columns=['elemid','nodeid1','nodeid2'],index=elemid)
    return(B)

#读取其他信息，包括单位缆索重量q、跨中节点号midnodeid、主缆跨中垂度f
def loadinfo(fname):
    A=np.loadtxt(fname)
    #n=A.shape[0]
    item=['q','midnodeid','fmid']
    B=pd.DataFrame(A,columns=['value'],index=item)
    q=B.loc['q','value']
    midnodeid=B.loc['midnodeid','value']
    fmid=B.loc['fmid','value']
    return(q,midnodeid,fmid)
    
#读取固定节点坐标号
def loadfixednode(fname):
    A=np.loadtxt(fname)
    return(A)

#求解缆索轴向长度，输入x为节点x坐标矩阵
def getL(x):
    max=x.max()
    min=x.min()
    L=max-min
    return(L)

#设定初始H0,输入吊杆y方向分力矩阵Ty，缆索轴向长度L，主缆垂度f
'''
def iniH0(Ty,L,f):
    sty=Ty.sum()
    H0=L*sty/(8*f)
    H1=H0
    H2=1.1*H0
    return(H1,H2)
'''
def iniH0(NODE,HANGER,FIXEDNODE,q,ymid):
    nodeid1=FIXEDNODE[1]
    nodeid2=FIXEDNODE[2]
    x1=NODE.loc[nodeid1,'x']
    x2=NODE.loc[nodeid2,'x']
    L=abs(x2-x1)

    f=(NODE.loc[nodeid1,'y']+NODE.loc[nodeid2,'y'])/2-ymid

    Ty=0
    Ty+=q*2*sqrt(L*L/4+f*f)
    for i in HANGER.loc[:,'hanid'].values:
        if (nodeid1-i)*(nodeid2-i)<=0:
            Ty+=HANGER.loc[i,'Ty']
    
    H0=L*Ty/f/8
    H1=H0
    H2=1.1*H0
    return(H1,H2)

#求解三维K
def get3DKyy(elemid,ELEM,NODE,H0,q):
    nodeid1=ELEM.loc[elemid,'nodeid1']
    nodeid2=ELEM.loc[elemid,'nodeid2']

    x1=NODE.loc[nodeid1,'x']
    x2=NODE.loc[nodeid2,'x']
    y1=NODE.loc[nodeid1,'y']
    y2=NODE.loc[nodeid2,'y']
    z1=NODE.loc[nodeid1,'z']
    z2=NODE.loc[nodeid2,'z']

    dL=abs(x2-x1)
    J=dL/2
    dphi1=-1/dL
    dphi2=1/dL

    ck1=dphi1*z1+dphi2*z2
    ck2=dphi1*y1+dphi2*y2

    K11=H0*dphi1*dphi1*dL+q*ck2/sqrt(1+ck1*ck1+ck2*ck2)*J*dphi1
    K12=H0*dphi1*dphi2*dL+q*ck2/sqrt(1+ck1*ck1+ck2*ck2)*J*dphi2
    K21=H0*dphi2*dphi1*dL+q*ck2/sqrt(1+ck1*ck1+ck2*ck2)*J*dphi1
    K22=H0*dphi2*dphi2*dL+q*ck2/sqrt(1+ck1*ck1+ck2*ck2)*J*dphi2

    return(K11,K12,K21,K22)

def get3DKyz(elemid,ELEM,NODE,H0,q):
    nodeid1=ELEM.loc[elemid,'nodeid1']
    nodeid2=ELEM.loc[elemid,'nodeid2']

    x1=NODE.loc[nodeid1,'x']
    x2=NODE.loc[nodeid2,'x']
    y1=NODE.loc[nodeid1,'y']
    y2=NODE.loc[nodeid2,'y']
    z1=NODE.loc[nodeid1,'z']
    z2=NODE.loc[nodeid2,'z']

    dL=abs(x2-x1)
    J=dL/2
    dphi1=-1/dL
    dphi2=1/dL

    ck1=dphi1*z1+dphi2*z2
    ck2=dphi1*y1+dphi2*y2

    K11=q*ck1/sqrt(1+ck1*ck1+ck2*ck2)*J*dphi1
    K12=q*ck1/sqrt(1+ck1*ck1+ck2*ck2)*J*dphi2
    K21=q*ck1/sqrt(1+ck1*ck1+ck2*ck2)*J*dphi1
    K22=q*ck1/sqrt(1+ck1*ck1+ck2*ck2)*J*dphi2

    return(K11,K12,K21,K22)

def get3DKzy(elemid,ELEM,NODE,H0,q):

    K11=0
    K12=0
    K21=0
    K22=0

    return(K11,K12,K21,K22)

def get3DKzz(elemid,ELEM,NODE,H0,q):
    nodeid1=ELEM.loc[elemid,'nodeid1']
    nodeid2=ELEM.loc[elemid,'nodeid2']

    x1=NODE.loc[nodeid1,'x']
    x2=NODE.loc[nodeid2,'x']
    y1=NODE.loc[nodeid1,'y']
    y2=NODE.loc[nodeid2,'y']
    z1=NODE.loc[nodeid1,'z']
    z2=NODE.loc[nodeid2,'z']

    dL=abs(x2-x1)
    J=dL/2
    dphi1=-1/dL
    dphi2=1/dL

    ck1=dphi1*z1+dphi2*z2
    ck2=dphi1*y1+dphi2*y2

    K11=H0*dphi1*dphi1*dL
    K12=H0*dphi1*dphi2*dL
    K21=H0*dphi2*dphi1*dL
    K22=H0*dphi2*dphi2*dL

    return(K11,K12,K21,K22)


#求解三维fy
def solve3Dfy(elemid,ELEM,NODE,H0,q):
    nodeid1=ELEM.loc[elemid,'nodeid1']
    nodeid2=ELEM.loc[elemid,'nodeid2']

    x1=NODE.loc[nodeid1,'x']
    x2=NODE.loc[nodeid2,'x']
    y1=NODE.loc[nodeid1,'y']
    y2=NODE.loc[nodeid2,'y']
    z1=NODE.loc[nodeid1,'z']
    z2=NODE.loc[nodeid2,'z']

    dL=abs(x2-x1)
    J=dL/2
    dphi1=-1/dL
    dphi2=1/dL

    ck1=dphi1*z1+dphi2*z2
    ck2=dphi1*y1+dphi2*y2

    fy1=H0*ck2*dL*dphi1+q*sqrt(1+ck1*ck1+ck2*ck2)*J
    fy2=H0*ck2*dL*dphi2+q*sqrt(1+ck1*ck1+ck2*ck2)*J

    return(fy1,fy2)

#求解三维fz
def solve3Dfz(elemid,ELEM,NODE,H0):
    nodeid1=ELEM.loc[elemid,'nodeid1']
    nodeid2=ELEM.loc[elemid,'nodeid2']

    x1=NODE.loc[nodeid1,'x']
    x2=NODE.loc[nodeid2,'x']
    y1=NODE.loc[nodeid1,'y']
    y2=NODE.loc[nodeid2,'y']
    z1=NODE.loc[nodeid1,'z']
    z2=NODE.loc[nodeid2,'z']

    dL=abs(x2-x1)
    J=dL/2
    dphi1=-1/dL
    dphi2=1/dL

    ck1=dphi1*z1+dphi2*z2
    ck2=dphi1*y1+dphi2*y2

    fz1=H0*ck1*dL*dphi1
    fz2=H0*ck1*dL*dphi2

    return(fz1,fz2)

'''
    fz=H0*(z1-z2)*dL

    return(fz)
'''

def assemK(ELEM,NODE,HANGER,FIXEDNODE,q,H0):
    nodenum=NODE.shape[0]
    elemnum=ELEM.shape[0]

    Kyy=np.zeros((nodenum,nodenum))
    Kyz=np.zeros((nodenum,nodenum))
    Kzy=np.zeros((nodenum,nodenum))
    Kzz=np.zeros((nodenum,nodenum))

    for i in range(elemnum):
        (Kyy11,Kyy12,Kyy21,Kyy22)=get3DKyy(i+1,ELEM,NODE,H0,q)
        (Kyz11,Kyz12,Kyz21,Kyz22)=get3DKyz(i+1,ELEM,NODE,H0,q)
        (Kzy11,Kzy12,Kzy21,Kzy22)=get3DKzy(i+1,ELEM,NODE,H0,q)
        (Kzz11,Kzz12,Kzz21,Kzz22)=get3DKzz(i+1,ELEM,NODE,H0,q)

        nodeid1=int(ELEM.loc[i+1,'nodeid1'])
        nodeid2=int(ELEM.loc[i+1,'nodeid2'])

        Kyy[nodeid1-1,nodeid1-1]+=Kyy11	
        Kyy[nodeid1-1,nodeid2-1]+=Kyy12
        Kyy[nodeid2-1,nodeid1-1]+=Kyy21
        Kyy[nodeid2-1,nodeid2-1]+=Kyy22

        Kyz[nodeid1-1,nodeid1-1]+=Kyz11	
        Kyz[nodeid1-1,nodeid2-1]+=Kyz12
        Kyz[nodeid2-1,nodeid1-1]+=Kyz21
        Kyz[nodeid2-1,nodeid2-1]+=Kyz22

        Kzy[nodeid1-1,nodeid1-1]+=Kzy11	
        Kzy[nodeid1-1,nodeid2-1]+=Kzy12
        Kzy[nodeid2-1,nodeid1-1]+=Kzy21
        Kzy[nodeid2-1,nodeid2-1]+=Kzy22

        Kzz[nodeid1-1,nodeid1-1]+=Kzz11	
        Kzz[nodeid1-1,nodeid2-1]+=Kzz12
        Kzz[nodeid2-1,nodeid1-1]+=Kzz21
        Kzz[nodeid2-1,nodeid2-1]+=Kzz22

    for j in range(nodenum):
        enodeid=j+1
        if enodeid in HANGER.loc[:,'hanid']:
            Ty=HANGER.loc[enodeid,'Ty']
            y_cable=NODE.loc[enodeid,'y']
            y_deck=HANGER.loc[enodeid,'y']
            z_cable=NODE.loc[enodeid,'z']
            z_deck=HANGER.loc[enodeid,'z']

            Kzy[enodeid,enodeid]-=(Ty*(z_cable-z_deck)/(y_cable-y_deck)**2)
            Kzz[enodeid,enodeid]+=Ty/(y_cable-y_deck)

    for k in FIXEDNODE:
        k=int(k)
        Kyy[:,k-1]=0
        Kyy[k-1,:]=0
        Kyy[k-1,k-1]=1.0

        Kzz[:,k-1]=0
        Kzz[k-1,:]=0
        Kzz[k-1,k-1]=1.0

        Kyz[k-1,:]=0
        Kzy[:,k-1]=0

    K1=np.append(Kyy,Kyz,axis=1)
    K2=np.append(Kzy,Kzz,axis=1)
    K=np.append(K1,K2,axis=0)

    return(K)

def assemf(ELEM,NODE,HANGER,FIXEDNODE,H0,q):
    nodenum=NODE.shape[0]
    elemnum=ELEM.shape[0]

    fy=np.zeros((nodenum,1))
    fz=np.zeros((nodenum,1))

    for i in range(elemnum):
        (fy1,fy2)=solve3Dfy(i+1,ELEM,NODE,H0,q)
        (fz1,fz2)=solve3Dfz(i+1,ELEM,NODE,H0)

        nodeid1=int(ELEM.loc[i+1,'nodeid1'])
        nodeid2=int(ELEM.loc[i+1,'nodeid2'])

        fy[nodeid1-1,0]-=fy1
        fy[nodeid2-1,0]-=fy2

        fz[nodeid1-1,0]-=fz1
        fz[nodeid2-1,0]-=fz2

    for j in range(nodenum):
        enodeid=j+1
        if enodeid in HANGER.loc[:,'hanid']:
            Tyy=HANGER.loc[enodeid,'Ty']
            Tzz=Tyy*(NODE.loc[enodeid,'z']-HANGER.loc[enodeid,'z'])/(NODE.loc[enodeid,'y']-HANGER.loc[enodeid,'y'])
            fy[enodeid-1,0]-=Tyy
            fz[enodeid-1,0]-=Tzz
    
    for k in FIXEDNODE:
        k=int(k)
        fy[k-1]=0
        fz[k-1]=0

    f=np.append(fy,fz,axis=0)
    return(f)

def innerlayer(ELEM,NODE,HANGER,FIXEDNODE,q,H0,midnodeid):
    nodenum=NODE.shape[0]
    elemnum=ELEM.shape[0]

    K=assemK(ELEM,NODE,HANGER,FIXEDNODE,q,H0)
    f=assemf(ELEM,NODE,HANGER,FIXEDNODE,H0,q)
    delta=np.linalg.solve(K,f)

    resnum=int(delta.shape[0])
    aa=int(resnum/2)
    deltay=delta[0:aa,:]
    deltaz=delta[aa:resnum+1,:]

    y0=NODE.loc[:,'y'].values
    z0=NODE.loc[:,'z'].values
    NODE.loc[:,'y']=y0.reshape(nodenum,1)+deltay
    NODE.loc[:,'z']=z0.reshape(nodenum,1)+deltaz

    nn1=1
    #f=NODE.loc[midnodeid,'y']
    #print("内循环第{}次，f={}".format(nn1,f))

    while np.max(delta)>1e-5:

        K=assemK(ELEM,NODE,HANGER,FIXEDNODE,q,H0)
        f=assemf(ELEM,NODE,HANGER,FIXEDNODE,H0,q)
        delta=np.linalg.solve(K,f)

        resnum=int(delta.shape[0])
        aa=int(resnum/2)
        deltay=delta[0:aa,:]
        deltaz=delta[aa:resnum+1,:]

        y0=NODE.loc[:,'y'].values
        z0=NODE.loc[:,'z'].values
        NODE.loc[:,'y']=y0.reshape(nodenum,1)+deltay
        NODE.loc[:,'z']=z0.reshape(nodenum,1)+deltaz

        nn1+=1
        if nn1>50000:
            print('内循环未收敛')
            break
        #f=NODE.loc[midnodeid,'y']
        #print("内循环第{}次，f={}".format(nn1,f))

    f=NODE.loc[midnodeid,'y']
    
    return(NODE,f)

def solveH(H0,H1,f0,f1,fmid):
    g0=f0-fmid
    g1=f1-fmid
    H2=(H0*g1-H1*g0)/(g1-g0)
    return(H2)

def main(nodepath,elempath,hangerpath,infopath,fixednodepath,outputpath):
    #设置输入文件位置
    #iteration=''

    NODE=loadnodes(nodepath)
    ELEM=loadelem(elempath)
    HANGER=loadhanger(hangerpath)
    (q,midnodeid,ymid)=loadinfo(infopath)
    FIXEDNODE=loadfixednode(fixednodepath)

    #初始化H0
    '''
    Ty=ELEM.loc[:,'Ty']
    x=NODE.loc[:,'x']
    L=getL(x)
    H0=iniH0(Ty,L,f)
    '''
    (H0,H1)=iniH0(NODE,HANGER,FIXEDNODE,q,ymid)

    (NODE,f0)=innerlayer(ELEM,NODE,HANGER,FIXEDNODE,q,H0,midnodeid)
    (NODE,f1)=innerlayer(ELEM,NODE,HANGER,FIXEDNODE,q,H1,midnodeid)

    H2=solveH(H0,H1,f0,f1,ymid)
    H0=H1
    H1=H2
    nn2=1
    #print("外循环第{}次".format(nn2))

    while abs(f1-ymid)>1e-5:
    
        (NODE,f0)=innerlayer(ELEM,NODE,HANGER,FIXEDNODE,q,H0,midnodeid)
        (NODE,f1)=innerlayer(ELEM,NODE,HANGER,FIXEDNODE,q,H1,midnodeid)
        
        H2=solveH(H0,H1,f0,f1,ymid)
        H0=H1
        H1=H2
        nn2+=1
        if nn2>10000:
            print('外循环未收敛')
            break
        #print("外循环第{}次".format(nn2))
    NODE.to_excel(outputpath)

main()