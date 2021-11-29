#!/usr/bin/env python3
from __future__ import print_function

import os,time,re,math,argparse

import numpy as np
from os import popen#,popen4
import ase.io
import os.path
from copy import deepcopy

from subprocess import Popen,PIPE # as popen # check_output
from sys import stdout,version_info
def Popen4(cmd,timeout=None):
    proc=Popen(cmd,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out,err=proc.communicate(timeout=timeout)
    #out=proc.stdout.readline() #does not work.
    return out,err

def Popen4(cmd):
    proc=Popen(cmd,shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    out,err=proc.communicate()
    #out=proc.stdout.readline() #does not work.

    if version_info[0] < 3: out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\n") if x!=""] #python3 adds \\n to the end, python2 adds \n
    else:out2=[x for x in str(out).replace("b'",'').replace("'",'').split("\\n") if x!=""]

    return out2,err


def grep(key,fname,n=-1):
    #Uses grep to get the nth line with the keyword from a given file. By default the last occurence is returned.
    try:
        #return popen4('grep -m %d "%s" %s '%(n,key,fname),"r")[1].readlines()[-1][0:-1] #don't take \n at the end.  #Fastest one!!!
        return Popen4('grep -m %d "%s" %s '%(n,key,fname))[0][0]
    except:
        return ""

def unique(mylist): #return the unique elements of a list without changing the first apperance order
    #np.unqie a sorted unique list.
    mylist=list(mylist)
    return reduce(lambda l, x: l.append(x) or l if x not in l else l, mylist, [])

def compare(atoms,atoms_ref,aT):#Compatible with ASE str. format 
    #atomic order of two input structure does not matter any longer for comparison.
    flag=0;    cnt=0#atom id
    atomIDs=[];    coords=[];     types={}

    tol=1e-6
    diff=len(atoms_ref)-len(atoms)
    k=0#shift of atom index btw
    apos=atoms.get_scaled_positions()
    arefpos=atoms_ref.get_scaled_positions()

    for i,ref_crd in enumerate(arefpos):
        flag=0
        for j,crd in enumerate(apos):
            if np.dot((ref_crd-crd).T,(ref_crd-crd))<tol:#matched with an atom from atoms.
                #print("Matched: ",ref_crd,crd,atoms_ref[i],atoms[j])
                if atoms[j].symbol != atoms_ref[i].symbol: #For catching the dopants.
                    types[str(i+1)]='D'
                    coords.append(list(ref_crd))
                    atomIDs.append([aT[i],i+1])

                flag=1;break

        if not flag:
            types[str(i+1)]='V'
            k+=1;
            coords.append(list(ref_crd))
            atomIDs.append([aT[i],i+1])

    if k>diff: print("compare: something wrong with comparing two input structures.")
    #print (diff,k,atomIDs,coords)
        
    return atomIDs,coords,types

def compare_old(atoms,atoms_ref,aT):#Compatible with ASE str. format 
    flag=0;    cnt=0#atom id
    atomIDs=[];    coords=[];     types={}

    tol=1e-6
    diff=len(atoms_ref)-len(atoms)
    k=0#shift of atom index btw
    apos=atoms.get_scaled_positions()
    arefpos=atoms_ref.get_scaled_positions()
    for i in range(len(atoms_ref)):
        try:
            crd=apos[i-k]
            ref_crd=arefpos[i]

            if atoms[i-k].symbol != atoms_ref[i].symbol: #For catching the dopants.
                atomIDs.append([aT[i],i+1])
                coords.append(ref_crd)
                types[str(i+1)]='D'
                continue

            for j in range(3): #For catching the missing atoms.
                if abs(crd[j]-ref_crd[j])>tol:
                    atomIDs.append([aT[i],i+1])
                    coords.append(ref_crd)
                    #types.append('V')
                    types[str(i+1)]='V'
                    k+=1;break
        except: #to catch the last one
            atomIDs.append([aT[i],i+1])
            coords.append(ref_crd)
            if atoms[-1].symbol != atoms_ref[-1].symbol:types[str(i+1)]='D'
            else:types[str(i+1)]='V'
            k+=1;

    if k>diff: print("compare: something wrong with comparing two input structures.")
    #print (diff,k,atomIDs,coords)
    return atomIDs,coords,types

def compare_old2(data,ref_data,aT):#Only compatible with VASP5 format (w/ extra line of atomic types)
    flag=0
    cnt=0#atom id
    atomIDs=[]
    coords=[]

    for i in range(len(ref_data)):
        ln=ref_data[i]
        if 'dir' in ln.lower() or 'car' in ln.lower(): flag=1;continue
        #print ln
        if flag: 
            cnt+=1
            if ln not in data:
                #if "SCEL" in ln: continue
                atomIDs.append([aT[cnt-1],cnt])
                coords.append(ln)
    return atomIDs,coords

def sorted_nicely( l ):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def sort_dict_by_val(data,index=0):
    #index: index to use if values of the dictionary is not scalar.
    keys=list(data.keys())
    if len(keys)==1: return keys

    skeys=[keys[0]]#sorted keys
    for i in range(len(keys)):
        key=keys[i]
        val=data[key]
        try: val=val[index]
        except:None
        flag=0
        for j in range(len(skeys)):
            key2=skeys[j]
            val2=data[key2]
            try: val2=val2[index]
            except:None
            if val <= val2 and i!=0: skeys.insert(j,key);flag=1;break
        if not flag and key not in skeys:skeys.append(key)#not to miss the largest element.

    return skeys

def saver(atoms,outf,otype,del_list=None):
#function for (optionally deleting a list of atoms) and saving the structures in a given format.
    #print ("Deleting: ",del_list)
    #print ("Output: ",outf)
    atoms_cp=atoms.copy()
    if del_list != None:
        try:
            #del atoms[[x for x in del_list]]
            del atoms_cp[del_list]
        except: print ("Error in deleting atoms: ",del_list)#;raise Exception 

    if otype=="vasp": ase.io.write(outf,atoms_cp,format=otype,direct=True,vasp5=True)
    else: ase.io.write(outf,atoms_cp,format=otype)

if __name__ == '__main__':
    
    #Common variables.
    outdir="./"
    seed='config'

    parser = argparse.ArgumentParser(description="This script is for analysing the [CASM] configuration enumeration results.")
    parser.add_argument('-t','--type', type=str,choices=['casm','castep','vasp'],required=False, default='casm',help='Input file type. Def: determined automatically from the extension.')
    parser.add_argument('-sug','--suggest', type=int,required=False,default=0, help='Suggest structures for the given vacancy/dopant occupancy using the current data. Def: No suggestion.')
    parser.add_argument('-tol', '--tol',type=float, default=1e6,help="Hull distance tolerance (in eV) for printing/saving the analysed configurations. Def: All configs are reported.")

    parser.add_argument('-ref','--ref_file', type=str,required=False, help='Reference structure for determining the number of dopants/vacancies. No default value. If using -casm the config #0 will be assigned.')

    #parser.add_argument('-s','--save', default=False,action='store_true',  help='For saving the structures in the "%s" folder. Default: .'%outdir)
    parser.add_argument('-s','--save', default=0,type=int, help='For saving the structures in the "%s" folder. Default: No structure saved. e.g. -s 3 will save 3 lowest energy for each vacancy/dopant composition. One can also seperately use -tol for saving multiple structures up to certain hull distance from the ground state at each stoichometry. -s option can be used jointly with suggestion to save the suggested structures.'%outdir)
    parser.add_argument('-o','--outf', type=str,required=False, help='Output file name. Def: Atuomatic naming: e.g. %s_ID.'%seed)
    parser.add_argument('-it','--itype', type=str,required=False, help='Input file type. Def: determined automatically from the extension.')
    parser.add_argument('-ot','--otype', default='res',type=str,required=False, help='Output file type, default: .res (SHELX format)')
    parser.add_argument('-ow','--overwrite', default=False,action='store_true', help='overwrite if output folder exists. Def: No')
    parser.add_argument('-mins','--reportMins', default=False,action='store_true', help='Reports the missing sites in the minimum energy configurations for each stoichiometry in mins.dat. Def:No')
    parser.add_argument('-hc','--hide_common', default=False,action='store_true', help='Hides the common vacancies/dopants in all configurations of a given stoichiometry to ease the analysis.. Def:No')

    args = parser.parse_args()

    if args.outf: seed=args.outf
    #Outfile extension.
    #if args.otype: 
    ext=args.otype.split("-")[-1]
    if args.outf: outf=args.outf+"."+ext

    #Determine the reference structure to find the number of vacancy/dopants.
    if args.ref_file: ref_file=args.ref_file
    else: 
        if args.type=="casm":ref_file='training_data/SCEL1_1_1_1_0_0_0/0/calctype.default/POSCAR'
        elif args.type=="castep" or args.type=="vasp":
            print ("A reference structure file should be given using -ref option for determining the vacancy/dopant composition") 
        #ref_file='training_data/SCEL1_1_1_1_0_0_0/0/POS' 
            exit()
        
    if args.reportMins:
        minsout=open("mins.dat",'w')

    outf=open('skipped_files','w')

    #Read the reference data, common for all types.

    atoms=ase.io.read(ref_file)#,format="vasp")#format determined from extension.
    ase.io.write("./ref_file",atoms,format="vasp",direct=True,vasp5=True)
    ref_data=[x[0:-1] for x in open("./ref_file",'r').readlines()]#don't take the end of line.
    #print (atoms)

    #Create the atom symbols list to be used in compare.
    aTypes=[]
    aT=ref_data[5].split()
    noAtoms=[int(x) for x in ref_data[6].split()]
    cn=0
    for no in noAtoms:
        for i in range(no):   aTypes.append("%s%s"%(aT[cn],i+1))
        cn+=1


    if args.type=="vasp" or args.type=="castep":
        if args.type=="castep":#Get the list of castep files in subfolders for CASTEP.
            fileList=[x for x in Popen4('find ./ -type f -name "*.castep"')[0]]
        
        elif args.type=="vasp": #Get the list of xml files in subfolders for VASP.
            fileList=[x for x in Popen4('find ./ -type f -name "vasprun.xml"')[0]]
        #xmlList=[x[0:-1] for x in popen4('grep -l "" `find ./ -type f -name "vasprun.xml"` | tail -n 1','r')[1].readlines()]

        if len(fileList)==0: 
            if args.type=="vasp":print ("No vasprun.xml was found. Terminating..."); exit()
            else: print ("No *.castep was found. Terminating..."); exit()

        print ('All atomic indices start from 1.')
        fileList=sorted_nicely(fileList)#We don't need to sort the xmlList, as they will later be sorted by energy value.
        #print(fileList)
        
        configs=[];vacs=[];hull_dists=[];energies=[]#most of these not needed to be stored.
        dic={};mins={}
        
        missing={}#This is the suggestion dictionary for the target vacancy composition.
        present={}
        common={} #dict reporting which missing atoms are common througout the set of configs for a given stoichiometry.
        types={}
        for x in fileList:#Loop over each config
            if args.type =='vasp':
                flag=0
                dirname="/".join(x.split("/")[0:-1])
                #check if the calculation has converged.
                if  grep("F=",dirname+"/OSZICAR",1)==" ": print ("%s has not converged."%dirname);continue

                if os.path.exists(dirname+"/OUTCAR"):
                    try:
                        a=Popen4("grep '\-\ Iteration' %s |wc -l"%(dirname+"/OUTCAR"))[0]
                        #print(a)
                        noSCFiter=int(a)

                    except: noSCFiter=0
                    #try: 
                    a=grep("NELM",dirname+"/OUTCAR",-1).split("=")[1].split(";")[0]
                    nomaxSCF=int(a)
                    #except: nomaxSCF=0
                else:
                    noSCFiter=0; nomaxSCF=0
                if noSCFiter==nomaxSCF: #means did not converge after max allowed steps.
                    print ("%s is not SCF-converged and will be excluded from analysis..."%dirname)
                    outf.write("%s\n"%dirname)
                    continue

                for y in x.split("/")[0:-1]:
                    try: int(y); configs.append("config_%s"%y);flag=1;break
                    except:None 
                if not flag: configs.append(dirname)
                #skey="e_fr_energy" #The free energy with the entropy corr.
                skey="e_wo_entrp" #CASM uses this instead.
                a=Popen4('grep  "%s" %s | tail -n 1'%(skey,x))[0][0]
                #print(a)
                energies.append(float(a.split()[-2])) 

                atoms_cur=ase.io.read('%s/POSCAR'%(dirname))

            elif args.type =='castep':
                dirname="/".join(x.split("/")[0:-1])
                #check if the calculation has converged.
                #skey="Final energy, E"
                skey="Final free energy (E-TS)"
                #skey="Total energy corrected for finite basis set"  #If IPRINT=2 in param.
                if 1:
                #try:
                #    try:
                    if 1:
                        a=Popen4('grep  "%s" %s | tail -n 1'%(skey,x))[0][0]
                        print(a)
                        e=float(a.split()[-2]) #If SCF converged then should be there.
                #    except:
                        
                        #e=float(Popen4('grep  "%s" %s | tail -n 1'%(skey,x))[0][0].split()[-2])

                    energies.append(e)
                #except:
                #    continue

                if grep("Geometry optimization completed successfully.",x) == ' ' : None #not geometry converged

                configs.append(x.split('/')[-1].split('.castep')[0])
                print(configs)
    
                if os.path.exists('%s/%s.res'%(dirname,configs[-1])):
                    atoms_cur=ase.io.read('%s/%s.res'%(dirname,configs[-1]))
                elif os.path.exists('%s/%s.cell'%(dirname,configs[-1])):
                    atoms_cur=ase.io.read('%s/%s.cell'%(dirname,configs[-1]))
                elif os.path.exists('%s/%s-out.cell'%(dirname,configs[-1])):
                    atoms_cur=ase.io.read('%s/%s-out.cell'%(dirname,configs[-1]))

            

            #still in file/config loop!!

            atIDs,coords,typ=compare(atoms_cur,atoms,aTypes)
            key=str(len(atIDs))#no of vancies/dopants.

            if key in dic:
                #dic[key][configs[-1]]=hull_dists[-1]
                dic[key][configs[-1]]=energies[-1]
            else:
                #dic[key]={configs[-1]:hull_dists[-1]}
                dic[key]={configs[-1]:energies[-1]}

            if key not in missing:
                missing[key]={}
                present[key]={}
                types[key]={}
                common[key]=[]
                mins[key]=[1e6,'']
                configs=[configs[-1]]#restarting configs for new vacancy.

            missing[key][configs[-1]]=[x[1] for x in atIDs]
            types[key][configs[-1]]=deepcopy(typ)

            if args.hide_common: 
                if len(configs)==1:#First config of the stoi. add all
                    common[key]=deepcopy(missing[key][configs[-1]])
                else:
                    for ms in common[key]: #error when new key is foudn !! CHECK!!!!
                        if ms not in missing[key][configs[-1]]: common[key].remove(ms)
                #print(configs[-1],common[key])
            if energies[-1]<=mins[key][0]:mins[key]=[energies[-1],configs[-1]]

            #Determine present atoms for each stoichometry(key).
            tmp=[ii for ii in range(1,len(aTypes)+1)]
            for ml in missing[key][configs[-1]]:  tmp.remove(ml)
            present[key][configs[-1]]=tmp

 
        ## END OF file/config iteration ##


#Convert absolute energies into hull_distance or relative energies here.
        #minE=min(energies)      
        keys=dic.keys()
        #keys.sort()
        #keys=sorted(keys)
        keys=sorted_nicely(keys)
        for key in keys:#Iterate over vacancy nos or stoichiometries.
            #minE=min(dic[key][:])
            minE=mins[key][0]
            minC=mins[key][1]
            if args.hide_common: 
                """ #Seems to be not needed !!!
                #manually fix the false positive missing atoms from the common list not covered in the first iteration. 
                for ms in missing[key][configs[0]]: 
                    for conf in configs:
                        if ms in common[key] and ms not in missing[key][conf]: common[key].remove(ms);break #break the conf loop.
                """
                print ("\n%d common vacancies/dopants are hidden for clarity: \n%s"%(len(common[key]),", ".join([str(x)+types[key][configs[0]][str(x)] for x in common[key]])))
            
            if args.reportMins:
                minsout.write("%s #no dopants/vacancies: %s\n"%(" ".join([str(x) for x in missing[key][minC]]),key))

            cnt=0 #no of saved structures
            print ("\n   %d dopant/vacancy:"%(float(key)))
            print ("%-25s %-18s %-s  \n%s"%('Config Id','Hull_dist[eV]','Dopant/Vacancy IDs',64*'-'))

            keys2=dic[key].keys()

            if 1:  #the function seems to be working well!!!
        #Sort by hull_dist and apply the hull_dist filter.
                keys2=sort_dict_by_val(dic[key])
            elif 1:
        #sort by config name.
                keys2=sorted_nicely(keys2) 
            else: None  #Random order

            #cind=0#config index
            for conf in keys2:#Iterate over configs.

                #Remove common missing atoms to ease the analysis.
                #if args.hide_common: 
                    #print ('bk',common[key])
                    #for cm in common[key]:  missing[key][conf].remove(cm)

                dic[key][conf]-=minE
                hd=dic[key][conf]
                if hd <= args.tol:
                    #print ("%-25s %-18.8f %s"%(conf,hd, ", ".join(['%5s [%s]'%(str(x[0]),str(x[1])) for x in atIDs])))
                    if args.hide_common:#Remove common missing atoms to ease the analysis.
                        print ("%-25s %-18.8f %s"%(conf,hd, ", ".join([' %s'%str(x)+types[key][conf][str(x)] for x in missing[key][conf] if x not in common[key]])))
                    else:
                        print ("%-25s %-18.8f %s"%(conf,hd, ", ".join([' %s'%str(x)+types[key][conf][str(x)] for x in missing[key][conf]])))
                    #cind+=1


                if args.save != 0:
                    if args.tol==1e6 :
                        if cnt<args.save:
                            fname=outdir+"/%s_%d_%d."%(seed,int(float(key)),int(cnt))+ext
                            saver(atoms,fname,args.otype,del_list=[x-1 for x in missing[key][conf]]) #atIDs not updates, using always the last one.
                            cnt+=1

                    else:#If user given a tolerance.
                        if hd <= args.tol and cnt<args.save:
                            fname=outdir+"/%s_%d_%d."%(seed,int(float(key)),int(cnt))+ext
                            saver(atoms,fname,args.otype,del_list=[x-1 for x in missing[key][conf]])
                            cnt+=1
        #exit()
        if args.reportMins: minsout.close()


    elif args.type=="casm":                    
        data=Popen4("casm query -k 'comp_n(Va)' hull_dist | grep SCEL1")[0]

        configs=[];vacs=[];hull_dists=[]
        dic={}
        for dt in data:
            x=dt.split()
            configs.append(x[0]);vacs.append(x[2]);hull_dists.append(float(x[-1])) #no need to store these.

            #xx={configs[-1]:hull_dists[-1]}
            key=vacs[-1]
            if key in dic:
                dic[key][configs[-1]]=hull_dists[-1]
            else:
                dic[key]={configs[-1]:hull_dists[-1]}
        del data


        missing={}#This is the suggestion dictionary for the target vacancy composition.
        present={}

        keys=dic.keys()
        keys=sorted(keys)
        for key in keys:
            cnt=0 #no of saved structures
            print ("\n   %d dopant/vacancy:"%(float(key)))
            print ("%-25s %-18s %-s  \n%s"%('Config Id','Hull_dist[eV]','Missing atom IDs',61*'-'))
            keys2=dic[key].keys()      

            if 1:
            #Sort by hull_dist and apply the hull_dist filter.
                keys2=sort_dict_by_val(dic[key])
            else:
            #sort by config name.
                keys2=sorted_nicely(keys2) 

            missing[key]=[]
            present[key]=[]
            for conf in keys2:#config names
                hd=dic[key][conf]
                #data=[x[0:-1] for x in open('training_data/%s/POS'%conf).readlines()]                
                #atIDs,coords=compare_old(data,ref_data,aTypes)
                atoms_cur=ase.io.read('training_data/%s/POS'%conf,format="vasp")
                atIDs,coords=compare(atoms_cur,atoms,aTypes)

                missing[key].append([x[1] for x in atIDs])

                #Determine present atoms for each stoichometry(key).
                tmp=range(1,len(aTypes)+1)
                for ml in missing[key][-1]:  tmp.remove(ml)
                present[key].append(tmp)

                if hd <= args.tol:
                    print ("%-25s %-18.8f %s"%(conf,hd, ", ".join(['%5s [%s]'%(str(x[0]),str(x[1])) for x in atIDs])))
                    """
                    if len(present[key])>len(missing[key]): ##BUNDA HATA!!
                        print ("%-25s %-18.8f %s"%(conf,hd, ", ".join(['%5s [%s]'%(str(x[0]),str(x[1])) for x in atIDs])))
                    else:

                        print ("%-25s %-18.8f %s"%(conf,hd, present[key][-1]))
    """

                if args.save != 0:
                    if args.tol==1e6 :
                        if cnt<args.save:
                            fname=outdir+"/%s_%d_%d."%(seed,int(float(key)),int(cnt))+ext
                            #saver(atoms,fname,args.otype,del_list=[x[1]-1 for x in atIDs])
                            saver(atoms,fname,args.otype,del_list=[x-1 for x in missing[key][conf]])
                            cnt+=1

                    else:#If user given a tolerance.
                        if hd <= args.tol and cnt < args.save:
                            fname=outdir+"/%s_%d_%d."%(seed,int(float(key)),int(cnt))+ext
                            #saver(atoms,fname,args.otype,del_list=[x[1]-1 for x in atIDs])
                            saver(atoms,fname,args.otype,del_list=[x-1 for x in missing[key][conf]])
                            cnt+=1



    #START OF COMMON PART!!
    keys=sorted_nicely(missing.keys())
    suggestions={}
    #MAKE COMPATIBLE with the dictionary format of missing and present!!!
    if args.suggest:
        cnt=0
        print ("\nSuggested structures with %d vacancy/dopant based on the available data..."%args.suggest)
        print ("Reference Vacancy  Missing Atom IDs")
        for key in keys:
            ikey=int(float(key))
            if ikey==0:continue
            if ikey>args.suggest: #swith to using present atoms rather than missing ones,
                nAtoms=len(aTypes)
                suggest=present[key]
                ikey=len(aTypes)-ikey
                suggestions[key]=[] #add missing atoms to this to be compatible with part below.
                tmp=[]
                for sg in suggest:
                    tmp=list(set(tmp+sg))
                    #print (len(tmp))
                    if len(tmp)==nAtoms-args.suggest: break
                if len(tmp)!=nAtoms-args.suggest: continue #to skip the upper boundary.
                for ind in range(1,nAtoms+1):
                    if ind not in tmp:
                        suggestions[key].append(ind)

            elif args.suggest%ikey==0: 
                rpt=args.suggest/ikey
                tmp=[]
                for i in range(rpt):
                    tmp.extend(missing[key][i])
                suggestions[key]=tmp
            else: continue
            print (int(float(key)),suggestions[key])
            if args.save != 0:
                #fname=outdir+"/%s_%d_%d."%(seed,int(float(key)),int(cnt))+ext
                #print ("Writing suggested structures in %s_ID files."%seed)
                seed='suggest'
                fname=outdir+"/%s_%d."%(seed,int(cnt))+ext
                saver(atoms,fname,args.otype,del_list=[x-1 for x in suggestions[key]])
                cnt+=1

outf.close()
exit()		    
                        
   # print (suggestions)        

###############
#DELETED PARTS#
###############
    #print dic['1']
    #sorted(dic.items(), key=lambda x: x[1])#sort ascending by values.
    #print dic['1']
    #for key,val in zip(dic.keys(),dic.values()):
    #for key,val in enumerate(dic[key]):

    #ref_data=open('training_data/SCEL1_1_1_1_0_0_0/0/POS','r').readlines()



    #c = list(set(a + b))


        #configs=[int(y)  for y in x.split("/")[0:-1] for x in xmlList]
        #configs=sorted_nicely(configs)
        #print (configs)

        #configs=["/".join(x.split("/")[-3:-1]) for x in xmlList]


        #if args.ref_file: ref_file=args.ref_file
        #else:ref_file='training_data/SCEL1_1_1_1_0_0_0/0/POS'


"""
    #Read the reference data, common for all types.
    ref_data=[x[0:-1] for x in open(ref_file,'r').readlines()]#don't take the end of line.
    atoms=ase.io.read(ref_file,format="vasp")
    #print (atoms)

    #Create the atom symbols list to be used in compare.
    aTypes=[]
    aT=ref_data[5].split()
    noAtoms=[int(x) for x in ref_data[6].split()]
    cn=0
    for no in noAtoms:
        for i in range(no):   aTypes.append("%s%s"%(aT[cn],i+1))
        cn+=1
        """    

"""
        #Determine present atoms for each stoichometry(key).
        present[key]=[]#range(1,len(aTypes)+1)
        for ml in missing[key]:
            tmp=[]
            for ind in range(1,len(aTypes)+1): #range(present[key]):
            #val=present[key][ind]
                if ind not in ml: 
                    tmp.append(ind)
            present[key].append(tmp)
            """
        #print (key,missing[key],present[key])
