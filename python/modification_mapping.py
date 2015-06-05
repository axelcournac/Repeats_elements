# we modificate the _filter_unmapped_fastq function in Mirny Library 
# to include a threshold in the mapping quality


def _filter_unmapped_fastq(in_fastq, in_sam, nonunique_fastq):
    '''Read raw sequences from **in_fastq** and alignments from
    **in_sam** and save the non-uniquely aligned and unmapped sequences
    to **unique_sam**.
    '''
    threshold_MQ = 30
    samfile = pysam.Samfile(in_sam)

    nonunique_ids = set()
    for read in samfile:
        tags_dict = dict(read.tags)
        read_id = read.qname
        read_mapq = read.mapq


        # If exists, the option 'XS' contains the score of the second
        # best alignment. Therefore, its presence means a non-unique alignment.
        if 'XS' in tags_dict or read.is_unmapped or read_mapq < threshold_MQ :
            nonunique_ids.add(read_id)
            #print(read_id)

    num_total, num_filtered = _filter_fastq(
        nonunique_ids, in_fastq, nonunique_fastq)
    
    #print(len(nonunique_ids))
    print(num_total, num_filtered)

    return num_total, num_filtered
