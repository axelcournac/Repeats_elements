# Script to align iteratively reads with the Mirny Lab library - May 2013
# It will create an object of type hdf5 (lib) containing all the information necessary
# To have the information in text format, convert the hdf5 file : python h5dictToTxt.py /home/axel/Bureau/test_mirny/tmp/hic_lib.hdf5 /home/axel/Bureau/test_mirny/tmp/mat

import mapping_MQ
from h5dict import h5dict
import genome
import os
import logging
import sys
import tempfile

print 'gettempdir():', tempfile.gettempdir()
print 'gettempprefix():', tempfile.gettempprefix()

tempfile.tempdir = '/data/temporary'
print 'We change the temporary directory.'
print 'gettempdir():', tempfile.gettempdir()


logging.basicConfig(level=logging.DEBUG)

if len(sys.argv) > 4:
     print(sys.argv[1])
     print(sys.argv[2])
     print(sys.argv[3])
     print(sys.argv[4])
     print(sys.argv[5])
     print(sys.argv[6])
else:
     if len(sys.argv) > 1:
        print( sys.argv[1] )
     else:
        print "No argument entered."

bank = sys.argv[1];
fast1= sys.argv[2];
fast2= sys.argv[3];
path_to_index=sys.argv[4];
path_to_fasta=sys.argv[5];
name_of_enzyme=sys.argv[6];

## A. Map the reads iteratively.
mapping_MQ.iterative_mapping(
    bowtie_path='/home/axel/Bureau/tools/bowtie2-2.2.4/bowtie2',
    bowtie_index_path=path_to_index,
    fastq_path=fast1,
    out_sam_path=bank+'1.bam',
    min_seq_len=20,
    len_step=5,
    nthreads=6,
    max_reads_per_chunk=5000000,   #optional, to split reads into smaller groups
    temp_dir=os.path.join(bank,'tmp'),  # optional, keep temporary files here
    #bowtie_flags='--very-sensitive  --score-min L,-0.6,-0.2')
    bowtie_flags='--very-sensitive')
print "Done!  Alignment 1 mate."
    
mapping_MQ.iterative_mapping(
    bowtie_path='/home/axel/Bureau/tools/bowtie2-2.2.4/bowtie2',
    bowtie_index_path=path_to_index,
    fastq_path=fast2,
    out_sam_path=bank+'/2.bam',
    min_seq_len=20,
    len_step=5,
    nthreads=6,
    max_reads_per_chunk=5000000,
    temp_dir=os.path.join(bank,'tmp'),
    #bowtie_flags='--very-sensitive  --score-min L,-0.6,-0.2')
    bowtie_flags='--very-sensitive')
print "Done!  Alignment 2d mate."
  

## B. Parse the mapped sequences into a Python data structure.
print 'Parse the generated BAMs...'
lib = h5dict(os.path.join(bank,'tmp')+'/hic_lib.hdf5')
genome_db = genome.Genome(path_to_fasta,gapFile='gap.txt',chrmFileTemplate='chr%s.fa', readChrms=['#','M'] )

print 'Done! Affectation lib hdf5.'

mapping_MQ.parse_sam(
    sam_basename1=bank+'1.bam',
    sam_basename2=bank+'2.bam',
    out_dict=lib,  
    genome_db=genome_db,
    keep_ids='True')

print 'Done! Mapping Parse.'

# C. Assign the ultra-sonic fragments to restriction fragments.
mapping_MQ.fill_rsites(
    lib=lib,
    genome_db=genome_db,
    enzyme_name=name_of_enzyme)
    
print 'Done! Assignment to restriction fragments.'

print 'End of the alignment.'

