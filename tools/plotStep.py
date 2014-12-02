import os
def plotStep(step):
        s1=str(step)
        s2=str(1000000+step)
        xfile1=s2[1:]+'.p3d'
        print xfile1
        jpgfile='jpgs/p'+s2[1:]+'.png'
        str1="sed 's/000000.p3d/"+xfile1+"/g'"+" p.lay.0>pp.lay"
        os.system(str1)
        os.system('tec360 pp.lay -b -p export.mcr -y try.png')
        str1="cp try.png "+jpgfile
        os.system(str1)
        os.system('rm try.png')
