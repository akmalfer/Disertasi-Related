import os
import numpy as np
import re
import scipy.integrate as integrate
from scipy import optimize


def lognorm(xi,si):
    f1 = 1/np.sqrt(2*np.pi)*si
    fin = -1/2*(np.log(xi)/si)**2
    return f1/xi*np.exp(fin)
    
def chi(n,b,s,si):
    def f(xi,n,b,s,si):
        f1 = n*np.log(xi*(s+b))
        f2 = -xi*(s+b)
        f3 = n*np.log(n)-n
        return np.exp(f1+f2-f3)*lognorm(xi,si)
    return integrate.quad(f,0,50,args=(n,b,s,si))[0]
    

resultsfile = open("resultsfixmchi.txt", "w")

#mass = np.array([0.001])
#mass = np.logspace(np.log10(1E-3),np.log10(med),10)
#mass = np.logspace((np.log10(200)+3.)/20.*6-3.,(np.log10(200)+3.)/20.*20.-3,15)
Zp = np.linspace(120,240,12)
for m in Zp: 
    nsim = 0
    j = 0
    while nsim < 9000:
        runback = "ATLAS_mll_check_mzp" + str(m) + "_v15r" + str(j)
        os.chdir("/home/akmal/MadGraph/MG5_aMC_v3_4_1/pp_llvv_ATLAS/")
        os.system('cp Cards/run_card_default.dat Cards/run_card.dat')
        lines = open("./Cards/run_card.dat", 'r').readlines()
        lines[126] = "   "+str(m-10)+"   = mmll    ! min invariant mass of l+l- (same flavour) lepton pair\n"
        lines[127] = "   "+str(m+10)+"  = mmllmax ! max invariant mass of l+l- (same flavour) lepton pair\n"
        out = open("./Cards/run_card.dat", 'w')
        out.writelines(lines)
        out.close()
        os.system('./bin/madevent generate_events ' + runback)
        os.system('../Delphes/root2lhco Events/' + runback + '/tag_1_delphes_events.root Events/' + runback + '/tag_1_delphes_events.lhco')
        os.system('rm ./Events/' + runback + '/tag_1_delphes_events.root')
        os.system('rm ./Events/' + runback + '/tag_1_pythia8_events.hepmc.gz')
        file2 = open('./Events/' + runback + '/' + runback + '_tag_1_banner.txt')
        banner = file2.readlines()
        xs = float(banner[-10][35:])
        print("Here", xs)
        os.chdir("/home/akmal/Madanalysis/madanalysis5-main/mllatlasdilepton/Build/SampleAnalyzer/User/")
        lines = open("./Analyzer/mllatlasdilepton.cpp", 'r').readlines()
        lines[125] = "                                                  if (mee < "+str(m-5)+" || mee > "+str(m+5)+")\n"
        lines[156] = "                                                  if (muu < "+str(m-5)+" || muu > "+str(m+5)+")\n"
        out = open("./Analyzer/mllatlasdilepton.cpp", 'w')
        out.writelines(lines)
        out.close()
        os.chdir("/home/akmal/Madanalysis/madanalysis5-main/mllatlasdilepton/Build/")
        os.system('make')
        my_file = open("list.txt", "w")
        my_file.write("/home/akmal/MadGraph/MG5_aMC_v3_4_1/pp_llvv_ATLAS/Events/" + runback + "/tag_1_delphes_events.lhco")
        my_file.close()
        os.system('./MadAnalysis5job list.txt')
        file = open('output.txt')
        content = file.readlines()
        sf = int(content[0])
        sfj = int(content[1])
        nsim = int(content[2])
        print(nsim, sf, sfj)   
        j = j + 1
        n1 = 1.67116*(134.5/111.6)*139*xs*1000*sf/10000
        n2 = 2.23985*(111.5/70.8)*139*xs*1000*sfj/10000
        b1 = n1
        b2 = n2
        s1 = 12/144
        s2 = 12/124
        chi01 = -2*np.log(chi(n1, b1, 0, s1))
        chi02 = -2*np.log(chi(n2, b2, 0, s2))
        print(n1, n2, chi01, chi02)
    
    run = "testATLAS_"+str(m)+"_v10"
    print("This is run ",m)
    #ATLAS
    os.chdir("/home/akmal/MadGraph/MG5_aMC_v3_4_1/pp_chichizp_zpll_ATLAS/")
    os.system('cp Cards/param_card_default.dat Cards/param_card.dat')
    lines = open("./Cards/param_card.dat", 'r').readlines()
    lines[18] = "   32 "+str(m)+" # mx \n"
    out = open("./Cards/param_card.dat", 'w')
    out.writelines(lines)
    out.close()
    os.system('./bin/madevent generate_events '+run)
    os.system('../Delphes/root2lhco Events/'+run+'/tag_1_delphes_events.root Events/'+run+'/tag_1_delphes_events.lhco')
    os.system('rm ./Events/'+run+'/tag_1_delphes_events.root')
    os.system('rm ./Events/'+run+'/tag_1_pythia8_events.hepmc.gz')
    os.chdir("/home/akmal/MadGraph/MG5_aMC_v3_4_1/pp_chichizp_ZA/")
    os.system('cp Cards/param_card_default.dat Cards/param_card.dat')
    lines = open("./Cards/param_card.dat", 'r').readlines()
    lines[18] = "   32 "+str(m)+" # mx \n"
    out = open("./Cards/param_card.dat", 'w')
    out.writelines(lines)
    out.close()
    os.system('./bin/madevent generate_events '+run)
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
    os.chdir("/home/akmal/Madanalysis/madanalysis5-main/atlasdilepton/Build/")
    my_file = open("list.txt", "w")
    my_file.write("/home/akmal/MadGraph/MG5_aMC_v3_4_1/pp_chichizp_zpll_ATLAS/Events/"+run+"/tag_1_delphes_events.lhco")
    my_file.close()
    os.system('./MadAnalysis5job list.txt')
    file = open('output.txt')
    content = file.readlines()
    sf = int(content[0])
    sfj = int(content[1])
    print(sf,sfj)
    nevent1 = 139*bf*2*xs*1000*sf/10000
    nevent2 = 139*bf*2*xs*1000*sfj/10000
    
    def chi2(g):
        return -2*np.log(chi(n1,b1,(g/0.302822)**2*nevent1,s1))-chi01-2*np.log(chi(n2,b2,(g/0.302822)**2*nevent2,s2))-chi02-3.841
        
    gsol = optimize.brentq(chi2,0,2)
    
    resultsfile.write(str(m)+" "+str(gsol)+"\n")
    
resultsfile.close()    
        
        
    
    
    

    
    




