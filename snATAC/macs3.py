from itertools import repeat
from multiprocessing import Pool
from datetime import datetime
import subprocess
import os 
import logging


def run_hmmratac(ipath, opath):
    # No option to return summits 
    # Cant apply till then for fixed let calls
    # unless we build summits ourself 
    file=os.path.basename(ipath)
    name=file[:-4] 

    logging.info(f"Run hmmratac: {name}")    
    cmd=f"macs3 hmmratac -i {ipath} -f BEDPE -n {name} --outdir {opath}"
    with open(opath + name + '_log.txt', 'w') as log:
        log.write(name)
        log.write("command: " + cmd)
        subprocess.run(cmd, shell=True, stdout=log, stderr=log)



def runp_hmmratac(p, idir, odir):
    files=os.listdir(idir)
    nfile=len(files)
    file_paths=[idir + file for file in files]
    file_outs=list(repeat(odir, nfile))
    # format args as 2d list for starmap
    args=list(zip(file_paths, file_outs))

    with Pool(processes=p) as pool:
        pool.starmap(func=run_hmmratac, iterable=args)


def callpeak_atac(ipath, opath):
    file=os.path.basename(ipath)
    name=file[:-4] 

    logging.info(f"Run callpeak narrow: {name}")    
    cmd=f"macs3 callpeak -t {ipath} -f BED -n {name} -g hs --nomodel --shift -100 --extsize 200 --call-summits --outdir {opath}"
    with open(opath + name + '_log.txt', 'w') as log:
        subprocess.run(cmd, shell=True, stdout=log, stderr=log)


def callpeak_ac(ipath, opath):
    # No shift 
    file=os.path.basename(ipath)
    name=file[:-4] 

    logging.info(f"Run callpeak narrow: {name}")    
    cmd=f"macs3 callpeak -t {ipath} -f BED -n {name} -g hs --nomodel --extsize 200 --call-summits --outdir {opath}"
    with open(opath + name + '_log.txt', 'w') as log:
        subprocess.run(cmd, shell=True, stdout=log, stderr=log)


def callpeak_me3(ipath, opath):
    file=os.path.basename(ipath)
    name=file[:-4] 

    logging.info(f"Run callpeak broad: {name}")    
    cmd=f"macs3 callpeak -t {ipath} -f BED -n {name} -g hs --broad --broad-cutoff 0.1 --nomodel --outdir {opath}"
    with open(opath + name + '_log.txt', 'w') as log:
        subprocess.run(cmd, shell=True, stdout=log, stderr=log)


def run_callpeak(p, idir, odir, mode='atac'):
    files=os.listdir(idir)
    nfile=len(files)
    file_paths=[idir + file for file in files]
    file_outs=list(repeat(odir, nfile))
    # format args as 2d list for starmap
    args=list(zip(file_paths, file_outs))

    if mode == 'atac':
        with Pool(processes=p) as pool:
            pool.starmap(func=callpeak_narrow, iterable=args)
    elif mode == 'ac':
        with Pool(processes=p) as pool:
            pool.starmap(func=callpeak_broad, iterable=args)
    elif mode == 'me3':
        with Pool(processes=p) as pool:
            pool.starmap(func=callpeak_broad, iterable=args)
    else:
        print("Improper mode. Pick one of atac, ac, or me3.") 




if __name__ == '__main__':
    
    logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.INFO)
    logging.info('Start.')
    startTime = datetime.now()

    ### Testing 
    #ipath="/tscc/nfs/home/rlancione/ps-epigen/users/rlan/LFNIH_Pt2/ATAC/Callpeak_Heps/mergefrags_out/17.bed"
    #opath="/tscc/nfs/home/rlancione/ps-epigen/users/rlan/LFNIH_Pt2/ATAC/Callpeak_Heps/callpeak_out/"
    #run_callpeak(ipath=ipath, opath=opath)

    # Parallel version 
    #ipath="/tscc/nfs/home/rlancione/ps-epigen/users/rlan/LFNIH_Pt2/ATAC/Callpeak_Heps/mergefrags_out/"
    #opath="/tscc/nfs/home/rlancione/ps-epigen/users/rlan/LFNIH_Pt2/ATAC/Callpeak_Heps/callpeak_out/"

    #os.makedirs(opath, exist_ok=True)
    #runp_callpeak(p=11, idir=ipath, odir=opath)

    
    ### High mem version
    ipath="/tscc/projects/ps-epigen/users/rlan/LFNIH_Pt2/Pairtag/SplitFrags/me3_mergefrags_out/"
    opath="/tscc/nfs/home/rlancione/ps-epigen/users/rlan/LFNIH_Pt2/Pairtag/Callpeak_me3/me3_callpeak_out/"
    os.makedirs(opath, exist_ok=True)
    
    files = os.listdir(ipath)
    
    for f in files:
        callpeak_me3(ipath=f'{ipath}{f}', opath=opath)


    logging.info((datetime.now() - startTime)) 
    logging.info('Finish.')


