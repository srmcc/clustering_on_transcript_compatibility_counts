# Runs the full Zeisel pipeline. 
# Start from a folder of downloaded SRR files. Zeisel's 3005 mouse brain cell dataset can be found at http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361.
# The pipeline runs kallisto to get the TCCs for each cell, ultimately outputting the TCC matrix (3005-by-#Eq classes).
# The pipeline also uses the TCC matrix to compute the pairwise distances between cells, resulting in a 3005-by-3005 distance matrix.

### Modified by Shannon McCurdy to not overwrite files, not subsample, and end analysis early.

import os
import getopt
import sys

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = dir_path +'/'

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:k:t:",["idir=","njobs=","hacked-kallisto-path=","reference-transcriptome="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python Zeisel_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-human-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)
    
SRA_dir=''
num_proc=1
kallipso_path=''
ref_transcriptome=''

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        SRA_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-k","--hacked-kallisto-path"):
        kallipso_path=arg
    elif opt in ("-t","--reference-transcriptome"):
        ref_transcriptome=arg
    
        
if (not SRA_dir) or (not kallipso_path) or (not ref_transcriptome):
    print ('usage is : \n python Zeisel_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-mouse-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)

if not os.path.exists(SRA_dir):
    os.system('python '+dir_path+'get_files.py')

if not os.path.exists(dir_path+'reads_with_UMIs/'):
    print('Extracting reads from SRAs...')
    os.system('mkdir -p ' + dir_path+'reads_with_UMIs/')
    os.system('rm -f ' + dir_path+'reads_with_UMIs/*')
    os.system('python '+dir_path+'process_SRA.py -i '+SRA_dir+' -o ' + dir_path+'reads_with_UMIs/ -n '+str(num_proc))
else:
    print('Reads already extracted from SRAs')
read_dir_base= dir_path+'reads_and_UMI_subsample'


sampling_suffix=['100']
sampling_suffix_ref='100'
RANGE=range(1)
RANGES=range(1,1)


read_dir_to_pass=read_dir_base+sampling_suffix_ref+"/"
if not os.path.exists(read_dir_to_pass):
    print('Separating reads and UMIs...')
    os.system('mkdir -p '+read_dir_to_pass)
    os.system('mkdir -p ' +dir_path+'tmp_dir/')
    os.system('rm -f '+ read_dir_to_pass+'*')
    os.system('rm -f ' +dir_path+'tmp_dir/*')
    os.system('python '+dir_path+'Clean_reads.py -i ' +dir_path+'reads_with_UMIs/ -o '+read_dir_to_pass+' '+
              '-t ' +dir_path+'tmp_dir/ -n '+str(num_proc))
    os.system('rmdir ' +dir_path+'tmp_dir')
else:
    print("Reads and UMIs already separated")


print('Sampling reads')
sampling_rates=['1']
read_dir_to_pass=read_dir_base+sampling_suffix_ref+"/"
for index in RANGES:
    print('Sampling '+sampling_rates[index]+' fraction of reads...')
    out_dir_to_pass=read_dir_base+sampling_suffix[index]+"/"
    os.system('mkdir -p '+out_dir_to_pass)
    os.system('rm -f '+out_dir_to_pass+'*')
    cmd='python '+dir_path+'sample_reads.py -i '+read_dir_to_pass+' -o '+out_dir_to_pass+' -k 100 -r '+sampling_rates[index]+' -n '+str(num_proc)
    os.system(cmd)
    

print('Generating the Kallisto index (with hacked kallisto)...')
index_path=dir_path+'kallisto_index/Zeisel_index.idx'
if not os.path.exists(dir_path+'kallisto_index'):
    os.system('mkdir -p ' +dir_path+'kallisto_index')
    os.system('rm -f ' +dir_path+'kallisto_index/*')
    os.system(kallipso_path+' index -i '+index_path+' '+ref_transcriptome)
    metadata_cmd=kallipso_path+' metadata '+index_path
    os.system(metadata_cmd)
else:
    print("Kallisto already indexed")
num_ec = sum(1 for line in open(dir_path+'kallisto_index/Zeisel_index.idx_ecmap.txt'))
print(num_ec)

print('Generating TCC (with hacked kallisto)...')
TCC_base_dir=dir_path+'transcript_compatibility_counts_subsample'
for index in RANGE:
    print('Running hacked kallisto on '+sampling_rates[index]+' fraction of reads...')
    TCC_dir=TCC_base_dir+sampling_suffix[index]+"/"
    read_dir_to_pass=read_dir_base+sampling_suffix[index]+"/"
    if not os.path.exists(TCC_dir):
        os.system('mkdir -p '+TCC_dir)
        os.system('rm -f '+TCC_dir+'*')
        os.system('python '+dir_path+'get_pseudoalignments.py -i '+read_dir_to_pass+' -o '+TCC_dir+' -k '+kallipso_path+ ' -t '+ index_path +' -n '+ str(num_proc))
    else:
        print("already generated TCC")
    ## making a dictionary of equivlance class index to transcript index  
    if not os.path.exists(dir_path +'eq_dict.dat'):
        eq_dict={}
        flnames=sorted([x for x in os.listdir(TCC_dir) if x.endswith('.class')])
        eq_class_hash=sum(1 for line in open(dir_path+'kallisto_index/Zeisel_index.idx_ecmap.txt'))
        for flname in flnames:
            with open(TCC_dir+flname) as flptr:
                for line in flptr:
                    line = line.strip()
                    vect = line.split()
                    if not vect[0].isdigit():
                        if vect[0] not in eq_dict:
                            eq_dict[vect[0]]=eq_class_hash
                            eq_class_hash+=1

        with open(dir_path +'eq_dict.dat', 'wb') as outfile:
            pickle.dump(eq_dict, outfile, pickle.HIGHEST_PROTOCOL)

print('Generating TCC distribution...')
TCC_dist_base_flname=dir_path+'Zeisel_TCC_distribution_subsample'
TCC_base_flname=dir_path+'Zeisel_TCC_subsample'
for index in RANGE:
    print('Getting the TCC dist for '+sampling_rates[index]+' fraction of reads...')
    TCC_dir=TCC_base_dir+sampling_suffix[index]+"/"
    TCC_dist_flname=TCC_dist_base_flname+sampling_suffix[index]+".dat"
    TCC_flname=TCC_base_flname+sampling_suffix[index]+".dat"
    if not os.path.exists(TCC_dist_flname):
        os.system('python '+dir_path+'get_tcc_dist.py -i '+TCC_dir+' -m '+str(num_ec)+' -t '+TCC_flname+' -d '+ TCC_dist_flname)
    else:
        print('TCC distribution already generated')

print('Generating pairwise distances...')
TCC_distance_base_flname='Zeisel_TCC_pairwise_JS_distance_subsample'
for index in RANGE:
    TCC_dist_flname=TCC_dist_base_flname+sampling_suffix[index]+".dat"
    TCC_distance_flname=TCC_distance_base_flname+sampling_suffix[index]+".dat"
    print('Getting  pairwise distances for '+sampling_rates[index]+' fraction of reads...')
    if not os.path.exists(TCC_distance_flname):
        os.system('python '+dir_path+'get_pairwise_distances.py '+TCC_dist_flname+' '+TCC_distance_flname+' '+str(num_proc))
    else:
        print('TCC distances already generated')

