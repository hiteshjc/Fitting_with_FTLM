import subprocess as sp
#import pylab
import numpy as np
import scipy.optimize as opt
import sys
#import pylab as P
import scipy
import time
import copy
from operator import itemgetter
import numpy
import os


def get_cv(Ts,fdir,seeds):
	cv=[]
	er=[]
	m_el=[]

	
	for i in range(seeds):
		f=open("./"+str(fdir)+"/out"+str(i)+".out",'r')
		count=0
		for line in f:
			count=count+1
			if count>1492 and count<1543:
				line=line.strip()
				#print (line)
				line=line.split(' ')
				er.append(float(line[2]))
				m_el.append(float(line[6]))
		f.close()		
	for t in Ts :
			beta=1.0/(t*0.086173)	
			numerator1=0.0
			numerator2=0.0
			denominator=0.0
			mag=0.0

			for m in range(len(er)):
				p=np.exp(-beta*er[m])*m_el[m]*m_el[m]
				numerator1=numerator1+(p*er[m])
				numerator2=numerator2+(p*er[m]*er[m])
				denominator=denominator+(p)
			e=numerator1/denominator
			e2=numerator2/denominator

			spheat=(e2-e*e)/(t*t)/16.0*1119.165   # Joules per mole per kelvin
			cv.append(spheat)
	#print("Magnetizations")
	cv=np.array(cv)

	return cv

def get_mag(exph,fdir,beta,seeds):

        mag=[]
        free_en=np.zeros((len(exph),2))
        #idx=0
        #h=0.0
        for idx,h in enumerate(exph):
            for pm in range(2):
                if pm==0:
                    folder='./'+str(fdir)+'/field'+str(idx)+'m'
                if pm==1:
                    folder='./'+str(fdir)+'/field'+str(idx)+'p'
                er=[]
                m_el=[]
                for i in range(seeds):
                        f=open(str(folder)+"/out"+str(i)+".out",'r')
                        count=0
                        for line in f:
                                count=count+1
                                if count>1492 and count<1543:
                                        line=line.strip()
                                        #print (line)
                                        line=line.split(' ')
                                        er.append(float(line[2]))
                                        m_el.append(float(line[6]))
                        f.close()
                denominator=0.0
                for m in range(len(er)):
                    p=np.exp(-beta*er[m])*m_el[m]*m_el[m]
                    denominator=denominator+(p)
                free_en[idx][pm]=1.0/beta*np.log(denominator)*17.277/16.0   # per site
            mag.append((free_en[idx][1]-free_en[idx][0])/0.1)
            #print("Magnetizations")
        return mag
 

        

def create_input_do_lanczos(J1,J2,J3,J4,gxy,gz,hx,hy,hz,seeds,fname,folder):

    for i in range(seeds):
        f1=open(fname,'w')
        f1.write("hamiltonian=ross_sym\n")
        f1.write('\n')
        f1.write("J1="+str(J1)+"\n")
        f1.write("J2="+str(J2)+"\n")
        f1.write("J3="+str(J3)+"\n")
        f1.write("J4="+str(J4)+"\n")
        f1.write("gxy="+str(gxy)+"\n")
        f1.write("gz="+str(gz)+"\n")
        f1.write('\n')
        f1.write("hx="+str(hx)+"\n")
        f1.write("hy="+str(hy)+"\n")
        f1.write("hz="+str(hz)+"\n")
        f1.write('\n')
        f1.write('seed='+str(i))
        f1.write('\n')
        f1.write("nkrylov=50\n")
        f1.write("diagonalize=true\n")
        f1.write("pairs=[[0,1],[0,2],[0,3],[0,7],[0,10],[0,13],[1,2],[1,3],[1,6],[1,11],[1,12],[2,3],[2,5],[2,8],[2,15],[3,4],[3,9],[3,14],[4,5],[4,6],[4,7],[4,9],[4,14],[5,6],[5,7],[5,8],[5,15],[6,7],[6,11],[6,12],[7,10],[7,13],[8,9],[8,10],[8,11],[8,15],[9,10],[9,11],[9,14],[10,11],[10,13],[11,12],[12,13],[12,14],[12,15],[13,14],[13,15],[14,15]]")
        f1.write('\n')
        f1.write('\n')
        f1.write("[Translationsforcomputingmomentumeigenvalues]\n")
        f1.write("t1=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]\n")
        f1.write("t2=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]\n")
        f1.write("t3=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]\n")
        f1.write('\n')
        f1.write("sublattice=[0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]")
        f1.close()
        output=str(folder)+"/out"+str(i)+".out"
        #sp.call(["g++","test.C","-o",output])
        sp.call(["./ed",fname,output])






        
def get_error(parms, Tss, expcvs, cvhs, exph111s, expmag111s, magTs, fdir, printout=0):
        err=0.0

        #J1=parms[0]
        #J2=parms[1]
        #J3=parms[2]
        #J4=parms[3]
        #J2=-J1
        #J3=J1
        #J4=J1
        #gxy=parms[4]
        gxy=0.0
        #gz=parms[5]
        Jx=parms[0]
        Jy=parms[1]
        Jz=parms[2]
        #Jxz=parms[3]
        Jxz=0.0
        #J4=J1
        #gxy=parms[4]
        #gxy=-0.2324
        #gz=2.35
        #gxy=parms[3]
        gz=2.401
        J1=1/6.0*(-2.0*Jx+np.sqrt(2.0)*Jxz+2.0*Jz)
        J2=1/3.0*(-2.0*Jx-2.0*np.sqrt(2.0)*Jxz-Jz)
        J3=1/6.0*(Jx-2.0*np.sqrt(2.0)*Jxz-3.0*Jy+2.0*Jz)
        J4=1/6.0*(-Jx+2.0*np.sqrt(2.0)*Jxz-3.0*Jy-2.0*Jz)

        temp=[0.,0.,0.,0.,0.]
        for num in range(len(cvhs)):
                h=cvhs[num]
                Ts=copy.deepcopy(Tss[num])
                expcv=copy.deepcopy(expcvs[num])

                hx=h/np.sqrt(3)
                hy=hx
                hz=hx
                fname='./'+str(fdir)+'/input.dat'
                folder='./'+str(fdir)
                seeds=25
                create_input_do_lanczos(J1,J2,J3,J4,gxy,gz,hx,hy,hz,seeds,fname,folder)
	
        
                cv=get_cv(Ts,fdir,seeds)
                # Compute error with respect to experiment 
                temperr=0.0
                ctr=0
                for n in range(len(Ts)):
                        temperature=Ts[n]
                        if temperature>1.0 :
                            temperr=temperr+np.power(expcv[n]-cv[n],2.0)
                            ctr=ctr+1
                temperr=temperr/float(ctr)
                temperr=np.sqrt(temperr)
                #temp[num]=temperr
                err=err+temperr

        for num in range(len(magTs)):
                T=magTs[num]
                exph=copy.deepcopy(exph111s[num])
                expmag=copy.deepcopy(expmag111s[num])
                beta=1.0/(T*0.0861738)

                for idx,h in enumerate(exph):
                    for i in range(2):
                        if i==0: 
                            hx=(h-0.05)/np.sqrt(3)
                            folder='./'+str(fdir)+'/field'+str(idx)+'m'
                        if i==1:
                            hx=(h+0.05)/np.sqrt(3)
                            folder='./'+str(fdir)+'/field'+str(idx)+'p'
                        if not os.path.exists(folder): os.mkdir(folder)    
                        hy=hx
                        hz=hx
                        fname=str(folder)+'/input.dat'
                        seeds=5
                        create_input_do_lanczos(J1,J2,J3,J4,gxy,gz,hx,hy,hz,seeds,fname,folder)
                mag111=get_mag(exph,fdir,beta,seeds)
                temperr=0.0
                ctr=0
                for n in range(len(expmag)):
                        temperr=temperr+np.power(expmag[n]-mag111[n],2.0)
                        ctr=ctr+1
                temperr=temperr/float(ctr)
                temperr=np.sqrt(temperr)
                #temp[num]=temperr
                err=err+temperr

        
          


        #err=err+np.power(abs(gxy-4.0),2.0)+np.power(abs(gz-2.0),2.0)
        #err=0.0*temp[0]+0.0*temp[1]+0.0*temp[2]+0.0*temp[3]+1.0*temp[4]
        print("Err = ",err)
        return err
        
        
Tss=[]
expcvs=[]
ndir=sys.argv[1]
it=sys.argv[2]
expmag111s=[]
exph111s=[]
#magTs=[0.5,4.0]
magTs=[0.5]
        
fnames=['./exp_data/Gao_Cv_0field_lesspts.dat','./exp_data/Gao_Cv_magField_1.2T_111_lesspts.dat','./exp_data/Gao_Cv_magField_2T_111_lesspts.dat','./exp_data/Gao_Cv_magField_4T_111_lesspts.dat']#,'./exp_data/Gao_Cv_magField_8T_111.dat','./exp_data/Gao_Cv_magField_14T_111.dat']
cvhs=[0.0,1.2,2.0,4.0]#,8.0,14.0]
#fnames=['./param_try2/Gao_Cv_0field.dat','./param_try2/Gao_Cv_magField_2T_111.dat','./param_try2/Gao_Cv_magField_4T_111.dat']
for fname in fnames:
    f=open(fname,'r')
    T=[]
    expcv=[]
    for line in f:
        line=line.strip()
        line=line.split(" ")
        T.append(float(line[0]))
        expcv.append(float(line[-1]))
    f.close()
    Tss.append(copy.deepcopy(T))
    expcvs.append(copy.deepcopy(expcv))
        
#fnames111=['./exp_data/Gao_magnetization_T=0.5K_lesspts.dat','./exp_data/Gao_magnetization_T=4K_lesspts.dat']
#fnames111=['./exp_data/Gao_magnetization_T=0.5K.dat','./exp_data/Gao_magnetization_T=1.8K.dat','./exp_data/Gao_magnetization_T=4K.dat']
fnames111=['./exp_data/Gao_magnetization_T=0.5K_lesspts.dat']
#fnames111=['./exp_data/Gao_magnetization_T=4K_lesspts.dat']

for fname in fnames111:
    f=open(fname,'r')
    h111=[]
    expmags=[]
    for line in f:
        line=line.strip()
        line=line.split(" ")
        h111.append(float(line[0]))
        expmags.append(float(line[-1]))
    f.close()
    exph111s.append(copy.deepcopy(h111))
    expmag111s.append(copy.deepcopy(expmags))



	
#cv=get_cv(Ts,er,m_el)	

scale1=0.1
scale2=0.1
scale3=0.1
scale4=0.4
ftempname='./'+str(ndir)+'/param_init_final_0gx0jxz_'+str(it)+'it.dat'
f3=open(ftempname,'w')
it=int(it)
for j in range(4*it):       #to generate different random numbers
    rn=np.random.random()

for i in range(it,it+3):
    fdir=str(ndir)+'/trial'+str(i)
    if not os.path.exists(fdir): os.mkdir(fdir)
    Jx=scale1*(2.0*np.random.random()-1)
    Jy=scale1*(2.0*np.random.random()-1)
    Jz=scale1*(2.0*np.random.random()-1)
    #gz=2.4+scale4*(2.0*np.random.random()-1.0)
    gz=2.4
    f3.write(str(round(Jx,3)))
    f3.write('  ')
    f3.write(str(round(Jy,3)))
    f3.write('  ')
    f3.write(str(round(Jz,3)))
    f3.write('  ')
    f3.write(str(round(gz,3)))
    f3.write('  ')
    #x=[scale1*(2.0*np.random.random()-1),scale2*(2.0*np.random.random()-1),scale3*(2.0*np.random.random()-1),2.4+scale4*(2.0*np.random.random()-1.0)]
    x=[Jx,Jy,Jz]#,gz]
    x=np.array(x)


    bnds=((-scale1,scale1),(-scale2,scale2),(-scale3,scale3))#,(2.0,2.8))

    bnds=np.array(bnds)
    (xopt)=opt.minimize(get_error,x,args=(Tss,expcvs,cvhs,exph111s,expmag111s,magTs,fdir),method='SLSQP',bounds=bnds,options={'ftol':0.0000001,'maxiter':10000,'disp':True})


    params=[]


    fname='./'+str(fdir)+'/input.dat'
    f=open(fname,'r')
    count=0
    for line in f:
        count=count+1
        if count > 2 and count < 9 :
            line=line.strip()
            line=line.split('=')
            params.append(float(line[1]))
    f.close()
    J1=params[0]
    J2=params[1]
    J3=params[2]
    J4=params[3]
    gxy=params[4]
    gz=params[5]

    Jx=(-4.0*J1-2.0*J2+J3-J4)/3.0
    Jy=-J3-J4
    Jz=(4.0*J1-J2+2.0*J3-2.0*J4)/3.0
    Jxz=(J1-J2-J3+J4)*np.sqrt(2.0)/3.0

    #f3.write("Final parameters")
    #f3.write('\n')
    #f3.write("Jx="+str(Jx)+"\n")
    #f3.write("Jy="+str(Jy)+"\n")
    #f3.write("Jz="+str(Jz)+"\n")
    #f3.write("Jxz="+str(Jxz)+"\n")
    #f3.write("gxy="+str(gxy)+"\n")
    #f3.write("gz="+str(gz)+"\n")
    f3.write(str(round(Jx,3)))
    f3.write('  ')
    f3.write(str(round(Jy,3)))
    f3.write('  ')
    f3.write(str(round(Jz,3)))
    f3.write('  ')
    f3.write(str(round(gz,3)))
    f3.write('\n')
f3.close()

#pylab.plot(h100,expmags100,color="red",lw=5,label="Expt 100")
#pylab.plot(h110,expmags110,color="blue",lw=5,label="Expt 110")
#pylab.plot(Ts,expcv,color="green",lw=5,label="Experiment at h = "+str(h)+" Tesla")
#pylab.plot(h100,mags100,color="red",markersize=10,marker="o",label="100")
#pylab.plot(h110,mags110,color="blue",markersize=10,marker="o",label="110")
#pylab.plot(Ts,cv,color="red",markersize=10,marker="o",label="best Param")
#pylab.plot(Ts,cv4,color="blue",markersize=10,marker="o",label="Cv exact diagonalization 4 site")

#pylab.legend(numpoints=1,loc="best")
#pylab.xlim([0,7])
#pylab.title("Cv vs T at h = "+str(h)+" Kelvin")
#pylab.ylabel("Cv")
#pylab.xlabel("T (K)")
#pylab.show()

	

