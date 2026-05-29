import os
import numpy as np
import re
import scipy.integrate as integrate
from scipy import optimize
from scipy.interpolate import interp1d

chix = np.array([2.7475,121.58,264.87,299.82,443.12,446.61,450.11,451.25,453.6,523.5])
zy = np.array([2067.6,2049.1,2025.8,2017.4,1908.2,1905.5,1902.7,1901.9,1900.2,1811.9])
mzf = interp1d(chix, zy)

n1 = 228
n2 = 99
n3 = 29
n4 = 17
b1 = 198
b2 = 120
b3 = 20.8
b4 = 9.3
s1 = 37/198
s2 = 16/120
s3 = 4.1/20.8
s4 = 2.3/9.3

def lognorm(xi,si):
    f1 = 1/np.sqrt(2*np.pi*si)
    fin = -1/2*(np.log(xi)/si)**2
    return f1/xi*np.exp(fin)
    
def chi(n,b,s,si):
    def f(xi,n,b,s,si):
        f1 = n*np.log(xi*(s+b))
        f2 = -xi*(s+b)
        f3 = n*np.log(n)-n
        return np.exp(f1+f2-f3)*lognorm(xi,si)
    return integrate.quad(f,0,50,args=(n,b,s,si))[0]
    
chi01 = -2*np.log(chi(n1,b1,0,s1))
chi02 = -2*np.log(chi(n2,b2,0,s2))
chi03 = -2*np.log(chi(n3,b3,0,s3))
chi04 = -2*np.log(chi(n4,b4,0,s4))

resultsfile = open("results.txt", "w")

#mass = np.array([0.001])
#mass = np.logspace(np.log10(1E-3),np.log10(med),10)
#mass = np.logspace((np.log10(200)+3.)/20.*6-3.,(np.log10(200)+3.)/20.*20.-3,15)
mchi = np.linspace(20,400,20)
for mc in mchi:
    if (2*mc-1) < 65:
       Zp = np.linspace(20,2*mc-1,8)
    else: 
       Zp = np.linspace(20,64,10)
       for m in Zp: 
        run = "optCMScont_"+str(mc)+"mzp"+str(m)+"_v1r2"
        print("This is run ",mc, m)
        #CMS
        mza = mzf(mc)
        print(mc,mza)
        '''
        os.chdir("/home/akmal/MadGraph/MG5_aMC_v3_4_1/pp_chichizp_zpll_CMS_light/")
        os.system('cp Cards/param_card_default.dat Cards/param_card.dat')
        lines = open("./Cards/param_card.dat", 'r').readlines()
        lines[18] = "   32 "+str(m)+" # mx \n"
        lines[19] = "   33 "+str(mza)+" # MZA\n"
        lines[20] = "  5000521 "+str(mc)+" # Mchi\n"
        out = open("./Cards/param_card.dat", 'w')
        out.writelines(lines)
        out.close()
        os.system('./bin/madevent generate_events '+run)
        os.system('../Delphes/root2lhco Events/'+run+'/tag_1_delphes_events.root Events/'+run+'/tag_1_delphes_events.lhco')
        os.system('rm ./Events/'+run+'/tag_1_delphes_events.root')
        os.system('rm ./Events/'+run+'/tag_1_pythia8_events.hepmc.gz')
        '''
        os.chdir("/home/akmal/MadGraph/MG5_aMC_v3_4_1/pp_chichizp_ZA_light/")
        '''
        os.system('cp Cards/param_card_default.dat Cards/param_card.dat')
        lines = open("./Cards/param_card.dat", 'r').readlines()
        lines[18] = "   32 "+str(m)+" # mx \n"
        lines[19] = "   33 "+str(mza)+" # MZA\n"
        lines[20] = "  5000521 "+str(mc)+" # Mchi\n"
        out = open("./Cards/param_card.dat", 'w')
        out.writelines(lines)
        out.close()
        os.system('./bin/madevent generate_events '+run)
        '''
        file2 = open('./Events/'+run+'/'+run+'_tag_1_banner.txt')
        banner = file2.readlines()
        xs = float(banner[-4][35:])
        #This only work if only one particle decays to -11 and 11
        for i in range(len(banner)):
            if "11-11" in banner[i].replace(" ", ""):
                epos = i
            if "-1111" in banner[i].replace(" ", ""):
                print(banner[i])
                epos = i
        bft = banner[epos]
        ind = len(bft) - len(bft.lstrip())
        bf = float(bft[ind:ind+12])
        print("Here",xs,bf)
        os.chdir("/home/akmal/Madanalysis/madanalysis5-main/mllcmsdilepton/Build/SampleAnalyzer/User/")
        lines = open("./Analyzer/mllcmsdilepton.cpp", 'r').readlines()
        lines[210] = "  if ((ose == 1 || osmu == 1) && MT2 > 100 && Nb == 0 && Missingenergy > 100 && ll.Pt() > 50 && (abs(mll-("+str(m)+")) < 5))\n"
        out = open("./Analyzer/mllcmsdilepton.cpp", 'w')
        out.writelines(lines)
        out.close()
        os.chdir("/home/akmal/Madanalysis/madanalysis5-main/cmsdilepton/Build/")
        my_file = open("list.txt", "w")
        my_file.write("/home/akmal/MadGraph/MG5_aMC_v3_4_1/pp_chichizp_zpll_CMS_light//Events/"+run+"/tag_1_delphes_events.lhco")
        my_file.close()
        os.system('./MadAnalysis5job list.txt')
        file = open('output.txt')
        content = file.readlines()
        sf01 = int(content[0])
        sf02 = int(content[1])
        sf03 = int(content[2])
        sf04 = int(content[3])
        print(sf01,sf02,sf03,sf04)
        nevent1 = 137*bf*2*xs*1000*sf01/10000
        nevent2 = 137*bf*2*xs*1000*sf02/10000
        nevent3 = 137*bf*2*xs*1000*sf03/10000
        nevent4 = 137*bf*2*xs*1000*sf04/10000
    
        def chi2(g):
            return -2*np.log(chi(n1,b1,(g/0.302822)**2*nevent1,s1))-chi01-2*np.log(chi(n2,b2,(g/0.302822)**2*nevent2,s2))-chi02-2*np.log(chi(n3,b3,(g/0.302822)**2*nevent3,s3))-chi03-2*np.log(chi(n4,b4,(g/0.302822)**2*nevent4,s4))-chi04-3.841
        
        if(chi2(0)*chi2(2)<0):
            gsol = optimize.brentq(chi2,0,2)
        else:
            gsol = 3
    
        resultsfile.write(str(mc)+" "+str(m)+" "+str(gsol)+"\n")
    
resultsfile.close()   

        
        
    
    
    

    
    




