#====================================================================
#
# YuE comment (03/29/19):
#
# Script realize viewer for *.h5 files, which are results of IOTA simulation. 
#
#====================================================================
#
# Import and setup IPython magics:
#
%matplotlib inline
%reload_ext autoreload
%autoreload 2
import sys, os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import h5py

printOutputFlag = 1

#
# Functions 'rewrite...' are used to check rhat manual writing of *.h5 files gives
# the same results as from the corresponding rssynergia methods
#
def rewriteBasicFile(crrntName,charge,mass,max,mean,min,num_particles, \
                     pz,real_num_particles,repetition,s,s_n,std):
    place_h5 = crrntName.find('.h5')
    newNameFile = crrntName[0:place_h5]+'_new.h5'
    print 'newNameFile = ',newNameFile
    f1 = h5py.File(newNameFile,'w')
    f1.create_dataset('charge',shape=(1,),data=charge)
    f1.create_dataset('mass',shape=(1,),data=mass)
    f1.create_dataset('max',shape=max.shape,data=max)
    f1.create_dataset('mean',shape=mean.shape,data=mean)
    f1.create_dataset('min',shape=min.shape,data=min)
    f1.create_dataset('num_particles',shape=num_particles.shape,data=num_particles)
    f1.create_dataset('pz',shape=pz.shape,data=pz)
    f1.create_dataset('real_num_particles',shape=real_num_particles.shape,data=real_num_particles)
    f1.create_dataset('repetition',shape=repetition.shape,data=repetition)
    f1.create_dataset('s',shape=s.shape,data=s)
    f1.create_dataset('s_n',shape=s_n.shape,data=s_n)
    f1.create_dataset('std',shape=std.shape,data=std)
    f1.flush()
    f1.close()
    return

def rewriteParticleFile(crrntName,charge,mass,particles,pz,rep,s_n,tlen):
    place_h5 = crrntName.find('.h5')
    newNameFile = crrntName[0:place_h5]+'_new.h5'
    print 'newNameFile = ',newNameFile
    f1 = h5py.File(newNameFile,'w')
    f1.create_dataset('charge',shape=(1,),data=charge)
    f1.create_dataset('mass',shape=(1,),data=mass)
    f1.create_dataset('particles',shape=particles.shape,data=particles)
    f1.create_dataset('pz',shape=(1,),data=pz)
    f1.create_dataset('rep',shape=(1,),data=rep)
    f1.create_dataset('s_n',shape=(1,),data=s_n)
    f1.create_dataset('tlen',shape=(1,),data=tlen)
    f1.flush()
    f1.close()
    return

def rewriteFull2File(crrntName,charge,corr,emitx,emitxy,emitxyz,emity,emitz,mass,max,mean, \
                     min,mom2,num_particles,pz,real_num_particles,repetition,s,s_n,std):
    place_h5 = crrntName.find('.h5')
    newNameFile = crrntName[0:place_h5]+'_new.h5'
    print 'newNameFile = ',newNameFile
    f1 = h5py.File(newNameFile,'w')
    f1.create_dataset('charge',shape=(1,),data=charge)
    f1.create_dataset('corr',shape=corr.shape,data=corr)
    f1.create_dataset('emitx',shape=emitx.shape,data=emitx)
    f1.create_dataset('emitxy',shape=emitxy.shape,data=emitxy)
    f1.create_dataset('emitxyz',shape=emitxyz.shape,data=emitxyz)
    f1.create_dataset('emity',shape=emity.shape,data=emity)
    f1.create_dataset('emitz',shape=emitz.shape,data=emitz)
    f1.create_dataset('mass',shape=(1,),data=mass)
    f1.create_dataset('max',shape=max.shape,data=max)
    f1.create_dataset('mean',shape=mean.shape,data=mean)
    f1.create_dataset('min',shape=min.shape,data=min)
    f1.create_dataset('mom2',shape=mom2.shape,data=mom2)
    f1.create_dataset('num_particles',shape=num_particles.shape,data=num_particles)
    f1.create_dataset('pz',shape=pz.shape,data=pz)
    f1.create_dataset('real_num_particles',shape=real_num_particles.shape,data=real_num_particles)
    f1.create_dataset('repetition',shape=repetition.shape,data=repetition)
    f1.create_dataset('s',shape=s.shape,data=s)
    f1.create_dataset('s_n',shape=s_n.shape,data=s_n)
    f1.create_dataset('std',shape=std.shape,data=std)
    f1.flush()
    f1.close()
    return

def plotParticles(particles,inputFile,rep,tlen):        
    particlesNumber = len(particles)
    x = np.zeros(particlesNumber)
    px = np.zeros(particlesNumber)
    y = np.zeros(particlesNumber)
    py =np.zeros(particlesNumber)
    z = np.zeros(particlesNumber)
    pz = np.zeros(particlesNumber)
    for n in range(particlesNumber):
        x[n] = particles[n,0]
        px[n] = particles[n,1]
        y[n] = particles[n,2]
        py[n] = particles[n,3]
        z[n] = particles[n,4]
        pz[n] = particles[n,5]    
    txtTitle = 'File "'+inputFile+'": '+str(particlesNumber)+' particles'
    strPrmtrs = str(rep)+' turns, s='+'{: 7.3f}'.format(tlen)+' m' 
    selection = 'plot'
    while selection == 'plot':
        selection = raw_input('\nYour selection: x-px, x-y, y-py, end')
        print 'You select ', selection
        if (selection == 'end'):
            return
        if (selection == 'x-px'):
            minX = 110.*np.min(x)
            maxX = 110.*np.max(x)
            minPx = 1.1e3*np.min(px)
            maxPx = 1.1e3*np.max(px)
#            print 'minX = ',minX,', maxX = ',maxX,'minPx = ',minPx,', maxPx = ',maxPx
            fig = plt.figure()
#            print ('fig = ', fig)
            plt.plot(100.*x,1.e3*px,'.r')
            plt.xlabel('x, cm',fontsize=12,color='m')
            plt.ylabel('px, mrad',fontsize=12,color='m')
            plt.xlim([minX,maxX])
            plt.ylim([minPx,maxPx])
            plt.title(txtTitle,fontsize=14, color='m')
            plt.text(.5*minX,.9*maxPx,strPrmtrs,fontsize=14, color='m')         
            plt.grid(True)
            plt.show()

        if (selection == 'x-y'):
            minX = 110.*np.min(x)
            maxX = 110.*np.max(x)
            minY = 110.*np.min(y)
            maxY = 110.*np.max(y)
#            print 'minX = ',minX,', maxX = ',maxX,'minY = ',minY,', maxY = ',maxY
            fig = plt.figure()
#            print ('fig = ', fig)
            plt.plot(100.*x,100.*y,'.r')
            plt.xlabel('x, cm',fontsize=12,color='m')
            plt.ylabel('y, cm',fontsize=12,color='m')
            plt.xlim([minX,maxX])
            plt.ylim([minY,maxY])
            plt.title(txtTitle,fontsize=14, color='m')
            strText = str(rep)+' turns, s='+'{: 7.3f}'.format(tlen)+' m' 
            plt.text(.5*minX,.9*maxY,strPrmtrs,fontsize=14, color='m')         
            plt.grid(True)
            plt.show()
        if (selection == 'y-py'):
            minY = 110.*np.min(y)
            maxY = 110.*np.max(y)
            minPy = 1.1e3*np.min(py)
            maxPy = 1.1e3*np.max(py)
#            print 'minY = ',minY,', maxY = ',maxY,'minPy = ',minPy,', maxPy = ',maxPy
            fig = plt.figure()
#            print ('fig = ', fig)
            plt.plot(100.*y,1.e3*py,'.r')
            plt.xlabel('y, cm',fontsize=12,color='m')
            plt.ylabel('py, mrad',fontsize=12,color='m')
            plt.xlim([minY,maxY])
            plt.ylim([minPy,maxPy])
            plt.title(txtTitle,fontsize=14, color='m')
            strText = str(rep)+' turns, s='+'{: 7.3f}'.format(tlen)+' m' 
            plt.text(.5*minY,.9*maxPy,strPrmtrs,fontsize=14, color='m')         
            plt.grid(True)
            plt.show()

        selection = 'plot'
    return

workDir = '/home/vagrant/jupyter/eidelyur/iota/example_run/'
workDir = '/home/vagrant/jupyter/eidelyur/iota/checkpoint/' # to check files from this directory

fileLS = os.popen('ls').readlines()
numbFiles = len(fileLS)
# print ('numbFiles = %d, fileLS = %s' % (numbFiles,fileLS))
fileList = [' ']*numbFiles
numb_h5 = 0
numb_basic = 0
numb_particles = 0
numb_full2 = 0
numb_local = 0
k = 0
for line in fileLS:
    if line[0:-1].count('.h5') == 1:
        fileList[k] = line[0:-1]
        if fileList[k].count('basic') == 1:
            numb_basic += 1
        if fileList[k].count('particles') == 1:
            numb_particles += 1
        if fileList[k].count('full2') == 1:
            numb_full2 += 1
        if fileList[k].count('0_') == 1:
            numb_local += 1
        k += 1
        numb_h5 += 1
print 'numb_h5 = ', numb_h5,'; "basic" = ',numb_basic,', "full2" = ', \
numb_full2,', "particles" = ',numb_particles,', "local" = ',numb_local        
# fileList[0:numb_h5]

#
# Input fileName:
#
existFileFlag = 0
print ('Files: %s' % fileList[0:numb_h5])
while existFileFlag == 0:
    crrntName = raw_input("\nSelect a name of file: ")
    print "Your file is: ", crrntName
    for name in fileList[0:numb_h5]:
#        print 'name = ', name, ', check = ',crrntName            
        if crrntName == name:
            existFileFlag = 1
            break
#    print('existFileFlag == ', existFileFlag)            
    if existFileFlag == 0:
        if ((crrntName == 'exit') or (crrntName == 'end') or (crrntName == 'stop')):
            print 'End of work!'
            break
        else:
            print 'File does not exist! Try again!'
    else:
        print 'File is found!'
#
# Type of selected file:
#
        nameFile = workDir + crrntName
        print 'nameFile = ',nameFile
        if (crrntName.count('basic') > 0):
            typeOfFile = 1
        if (crrntName.count('particles') > 0):
            typeOfFile = 2
        if (crrntName.count('full2') > 0):
            typeOfFile = 3
        if (crrntName.count('local') > 0):
            typeOfFile = 4
        print 'typeOfFile = ',typeOfFile
#
# Structure of the selected file:
#
        fp0 = h5py.File(nameFile,'r')
        fp0_len = len(fp0)              # = 7
        print 'fp0_len = ',fp0_len
        namesAttr = [' ']*fp0_len
        indxAttr = [' ']*fp0_len
        fp0_keys = fp0.keys()
        print 'fp0_keys = ',fp0_keys
        k = 0
        for name in fp0_keys:
            shapeName = fp0[name].shape
            namesAttr[k] = str(name)
            stringHelp = str(shapeName)[1:len(str(shapeName))-1] 
            countCommas = stringHelp.count(',')
            print 'For name %s shape = %s => %s; commas = %d' % (name,shapeName,stringHelp,countCommas)
            if countCommas == 0:
                dimnsn = 0
                numbIndx = ['0']
            else:
                numbIndx = stringHelp.split(',')
                dimnsn = len(numbIndx)
                print 'numbIndx = ',numbIndx,', dimnsn = ',dimnsn
                if numbIndx[dimnsn-1] == '':
                    dimnsn = dimnsn-1
                    numbIndx = numbIndx[0:dimnsn]
                    print 'For name ',name,' numbIndx = ',numbIndx,' => dim = ',dimnsn
            indxAttr[k] = numbIndx
            print 'For name ',name,' numbAttr = ',indxAttr[k],' => dim = ',dimnsn, \
                  '; lenIndx = ',len(indxAttr[k])
            k += 1
#
# Values of data:
#
        if (typeOfFile == 1):
#
# Files like 'basic_strct*.h5':
#
            charge = float(fp0[namesAttr[0]].value)
            mass = float(fp0[namesAttr[1]].value)
            max = np.asfarray(fp0[namesAttr[2]].value,float)
            mean = np.asfarray(fp0[namesAttr[3]].value,float)
            min = np.asfarray(fp0[namesAttr[4]].value,float)
            num_particles = fp0[namesAttr[5]].value
            pz = np.asfarray(fp0[namesAttr[6]].value,float)
            real_num_particles = np.asfarray(fp0[namesAttr[7]].value,float)
            repetition = fp0[namesAttr[8]].value
            s = np.asfarray(fp0[namesAttr[9]].value,float)
            s_n = np.asfarray(fp0[namesAttr[10]].value,float)
            std = np.asfarray(fp0[namesAttr[11]].value,float)
            if (printOutputFlag == 1):
                print '    charge = ',charge
                print '    mass = ',mass
                print '    max',indxAttr[2],' = ',max
                print '    min',indxAttr[3],' = ',min
                print '    mean',indxAttr[4],' = ',mean
                print '    num_particles = ',num_particles
                print '    pz = ',pz
                print '    real_num_particles = ',real_num_particles
                print '    repetition = ',repetition
                print '    s = ',s
                print '    s_n = ',s_n
                print '    std',indxAttr[11],' = ',std
            selection = 'loop'
            while selection == 'loop':
                selection = raw_input('\nYour selection: file, end')
                print 'You select ', selection
                if (selection == 'end'):
                    break
                if (selection == 'file'):
                    rewriteBasicFile(crrntName,charge,mass,max,mean,min,num_particles, \
                                     pz,real_num_particles,repetition,s,s_n,std)

        if (typeOfFile == 2):
#
# Files like 'particles_strct*_*.h5':
#
            charge = float(fp0[namesAttr[0]].value)
            mass = float(fp0[namesAttr[1]].value)
            particles = np.asfarray(fp0[namesAttr[2]].value,float)
            pz = float(fp0[namesAttr[3]].value)
            rep = int(fp0[namesAttr[4]].value)
            s_n = float(fp0[namesAttr[5]].value)
            tlen = float(fp0[namesAttr[6]].value)
            if (printOutputFlag == 1):
                print '    charge = ',charge
                print '    mass = ',mass
                print '    particles',indxAttr[2],' = ',particles
                print '    pz = ',pz
                print '    rep = ',rep
                print '    s_n = ',s_n
                print '    tlen = ',tlen
            selection = 'loop'
            while selection == 'loop':
                selection = raw_input('\nYour selection: plot, file, end')
                print 'You select ', selection
                if (selection == 'end'):
                    break
                if (selection == 'plot'):
                    plotParticles(particles,crrntName,rep,tlen)
                if (selection == 'file'):
                    rewriteParticleFile(crrntName,charge,mass,particles,pz,rep,s_n,tlen)

        if (typeOfFile == 3):
#
# Files like 'full2_strct*.h5':
#
            charge = float(fp0[namesAttr[0]].value)
            corr = np.asfarray(fp0[namesAttr[1]].value,float)
            emitx = np.asfarray(fp0[namesAttr[2]].value,float)
            emitxy = np.asfarray(fp0[namesAttr[3]].value,float)
            emitxyz = np.asfarray(fp0[namesAttr[4]].value,float)
            emity = np.asfarray(fp0[namesAttr[5]].value,float)
            emitz = np.asfarray(fp0[namesAttr[6]].value,float)
            mass = float(fp0[namesAttr[7]].value)
            max = np.asfarray(fp0[namesAttr[8]].value,float)
            mean = np.asfarray(fp0[namesAttr[9]].value,float)
            min = np.asfarray(fp0[namesAttr[10]].value,float)
            mom2 = np.asfarray(fp0[namesAttr[11]].value,float)
            num_particles = fp0[namesAttr[12]].value
            pz = np.asfarray(fp0[namesAttr[13]].value,float)
            real_num_particles = np.asfarray(fp0[namesAttr[14]].value,float)
            repetition = fp0[namesAttr[15]].value
            s = np.asfarray(fp0[namesAttr[16]].value,float)
            s_n = np.asfarray(fp0[namesAttr[17]].value,float)
            std = np.asfarray(fp0[namesAttr[18]].value,float)
            if (printOutputFlag == 1):
                print '    charge = ',charge
                print '    corr',indxAttr[1],' = ',corr
                print '    emitx',indxAttr[2],' = ',emitx
                print '    emitxy',indxAttr[3],' = ',emitxy
                print '    emitxyz',indxAttr[4],' = ',emitxyz
                print '    emity',indxAttr[5],' = ',emity
                print '    emitz',indxAttr[6],' = ',emitz
                print '    mass = ',mass
                print '    max',indxAttr[8],' = ',max
                print '    mean',indxAttr[10],' = ',mean
                print '    min',indxAttr[9],' = ',min
                print '    mom2',indxAttr[11],' = ',mom2
                print '    num_particles',indxAttr[12],' = ',num_particles
                print '    pz',indxAttr[13],' = ',pz
                print '    real_num_particles',indxAttr[14],' = ',real_num_particles
                print '    repetition',indxAttr[15],' = ',repetition
                print '    s',indxAttr[16],' = ',s
                print '    s_n',indxAttr[17],' = ',s_n
                print '    std',indxAttr[18],' = ',std
            selection = 'loop'
            while selection == 'loop':
                selection = raw_input('\nYour selection: file, end')
                print 'You select ', selection
                if (selection == 'end'):
                    break
                if (selection == 'file'):
                    rewriteFull2File(crrntName,charge,corr,emitx,emitxy,emitxyz,emity,emitz,mass, \
                                     max,mean,min,mom2,num_particles,pz,real_num_particles, \
                                     repetition,s,s_n,std)

        if (typeOfFile == 4):
            loc_num = fp0[namesAttr[0]].value
            loc_particles = np.asfarray(fp0[namesAttr[1]].value,float)
            if (printOutputFlag == 1):
                print '    loc_num = ',loc_num
                print '    loc_particles',indxAttr[1],' = ',loc_particles
            while selection == 'loop':
                selection = raw_input('\nYour selection: end')
                print 'You select ', selection
                if (selection == 'end'):
                    break    
#
# Files like '0_local_strct*.h5':
#
        fp0.close()
        print ('Files: %s' % fileList[0:numb_h5])
        existFileFlag = 0


''''
#-------------------------------------
# #This is doesn't work(!):
#
# from Tkinter import Label
# widget = Label(None,'trxt='Hello, Yury')
# widget.pack()
# widget.mainloop()
#-------------------------------------

'''
    