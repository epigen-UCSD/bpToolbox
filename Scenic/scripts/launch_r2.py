


### This script will launch the secondary scripts to setup the scenic run 
# 1. cistopic database merging
# 2. topic modeling and binarization
# 3. scenic outdir and config prep

def cistarget_merge(wdir, project):
    # fasta path
    fasta = f"{wdir}prepare_fasta/consensus.fa"    # Format project path

    # Replace strings in file with extracted idir & pdir
    with open("/tscc/projects/ps-epigen/users/rlan/sandbox/scenic/templates/cistarget_merge.sh", 'rt') as fh:
        shell = fh.read()
        #replace all occurrences of the required string
        shell = shell.replace('[IDIR]', f'{wdir}/cistarget_db/out/partial_db/')
        shell = shell.replace('[PROJECT]', project)
        shell = shell.replace('[OUT_DIR]', f"{wdir}/cistarget_db/out/")


    with open(f"{wdir}cistartet_db/cistarget_merge.sh", "w") as fh:
        fh.write(shell)
