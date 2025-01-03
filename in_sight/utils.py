import pysam
from pysam import FastaFile, VariantFile, AlignmentFile

def check_chromosome_prefix(bam_path, chromosome):
    """
    Check if the BAM file uses 'chr' prefix for chromosome names.
    :param bam_path: Path to the BAM file.
    :param chromosome: Chromosome number as a string (e.g., "15").
    :return: The correct chromosome name to use with this BAM file.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        try:
            # Try fetching reads with 'chr' prefix
            if bam_file.count('chr' + chromosome) > 0:
                return 'chr' + chromosome
        except ValueError:
            pass

        try:
            # Try fetching reads without 'chr' prefix
            if bam_file.count(chromosome) > 0:
                return chromosome
        except ValueError:
            pass

    return None  # If neither worked, chromosome name is not found

# Example usage
# bam_path = 'data/bam/SRR19794569_200000000.bam'
# chromosome = '15'
# correct_chromosome_name = check_chromosome_prefix(bam_path, chromosome)
# if correct_chromosome_name:
#     print(f"Correct chromosome name to use: {correct_chromosome_name}")
# else:
#     print("Chromosome name not found in BAM file.")


def adjust_chromosome_name(file_path, chromosome, file_type='fasta'):
    """
    Adjust the chromosome name based on the naming convention used in the file.
    
    :param file_path: Path to the file (BAM, Fasta, or VCF).
    :param chromosome: The chromosome name to adjust.
    :param file_type: The type of file ('bam', 'fasta', or 'vcf').
    :return: Adjusted chromosome name or None if not found.
    """
    try:
        if file_type == 'fasta':
            with FastaFile(file_path) as fasta_file:
                if chromosome in fasta_file.references:
                    return chromosome
                alternative_chrom = ('chr' + chromosome) if not chromosome.startswith('chr') else chromosome[3:]
                if alternative_chrom in fasta_file.references:
                    return alternative_chrom
        elif file_type == 'bam':
            with AlignmentFile(file_path, "rb") as bam_file:
                if chromosome in bam_file.references:
                    return chromosome
                alternative_chrom = ('chr' + chromosome) if not chromosome.startswith('chr') else chromosome[3:]
                if alternative_chrom in bam_file.references:
                    return alternative_chrom
        elif file_type == 'vcf':
            with VariantFile(file_path) as vcf_file:
                chrom_names = [contig for contig in vcf_file.header.contigs]
                if chromosome in chrom_names:
                    return chromosome
                alternative_chrom = ('chr' + chromosome) if not chromosome.startswith('chr') else chromosome[3:]
                if alternative_chrom in chrom_names:
                    return alternative_chrom
    except ValueError as e:
        print(f"Error adjusting chromosome name: {e}")
        return None

    return None  # If the chromosome name doesn't fit any convention

# Example usage
# fasta_path = 'data/ref/hg38.fa.gz'
# bam_path = 'data/bam/SRR19794569_200000000.bam'
# chromosome = '15'
# adjusted_chrom_fasta = adjust_chromosome_name(fasta_path, chromosome, 'fasta')

def get_flanking_sequences(ref_path, motif, flank_length=3, start=0):
    """
    Fetches the left and right flanking sequences of a specified length for a given STR motif using pysam.

    Parameters:
    - ref_path: Path to the reference genome.
    - motif: The motif object containing chromosome, start, and end positions.
    - flank_length: The number of bases to fetch from both sides of the STR.
    - start: The starting number of coordinates

    Returns:
    - A tuple containing the left and right flanking sequences.
    """
    ref = pysam.FastaFile(ref_path)
    adjusted_chrom_name = adjust_chromosome_name(ref_path, motif['chrom'], 'fasta')

    # Adjust start and end positions to include flanking sequences
    # Ensuring the start position does not go below 0
    start_with_flank = max(0, motif['start']-(1-start) - flank_length)
    end_with_flank = motif['end'] + flank_length

    # Fetch the extended sequence
    extended_sequence = ref.fetch(adjusted_chrom_name, start_with_flank, end_with_flank)

    # Calculate the positions of the STR and flanking sequences within the extended sequence
    left_flank_start = start
    left_flank_end = flank_length
    right_flank_start = len(extended_sequence) - flank_length
    right_flank_end = len(extended_sequence)

    # Extract the flanking sequences
    left_flank = extended_sequence[left_flank_start:left_flank_end].upper()
    right_flank = extended_sequence[right_flank_start:right_flank_end].upper()

    ref.close()  # Don't forget to close the reference file after fetching the sequences

    return left_flank, right_flank

def get_head_softclip(read):
    """
    Extracts the soft-clipped sequence and its length at the head (5' end) of a read.
    
    Parameters:
    read (pysam.AlignedSegment): The read from which to extract the soft-clipped sequence.
    
    Returns:
    tuple: A tuple containing the soft-clipped sequence and its length at the head of the read.
           If no soft-clipped sequence exists or the read has no valid CIGAR string, returns an empty string and 0.
    """
    # Extract the CIGAR tuple
    cigar = read.cigartuples
    
    # Check if the CIGAR string is valid
    if cigar is None or len(cigar) == 0:
        return "", 0

    # Check if the first CIGAR operation is a soft clip (5' end)
    if cigar[0][0] == 4:  # 4 corresponds to the S operation in CIGAR
        softclip_length = cigar[0][1]
        # Extract the soft-clipped sequence from the read
        softclip_sequence = read.query_sequence[:softclip_length]
        return softclip_sequence, softclip_length
    return "", 0

def get_head_hardclip(read):
    """
    Extracts the hard-clipped length at the head (5' end) of a read.
    
    Parameters:
    read (pysam.AlignedSegment): The read from which to extract the hard-clipped length.
    
    Returns:
    int: The length of the hard-clipped sequence at the head of the read.
         If no hard-clipped sequence exists or the read has no valid CIGAR string, returns 0.
    """
    # Extract the CIGAR tuple
    cigar = read.cigartuples
    
    # Check if the CIGAR string is valid
    if cigar is None or len(cigar) == 0:
        return 0

    # Check if the first CIGAR operation is a hard clip (5' end)
    if cigar[0][0] == 5:  # 5 corresponds to the H operation in CIGAR
        return cigar[0][1]
    return 0