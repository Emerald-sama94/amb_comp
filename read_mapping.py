import argparse
import os
import shutil
import subprocess
import sys
import tempfile

def parse_arguments():
    """Process inputs into arguments."""
    
    parser = argparse.ArgumentParser(description="Process genome assembly and genomic reads inputs.")
    
    # Genome assembly (required argument)
    parser.add_argument('assembly', type=str, help='Path to the genome assembly.')

    # Short reads (optional, with flags)
    parser.add_argument('-s', '--short1', type=str, help='Path to the first (or only) short read file in FASTQ format.')
    parser.add_argument('-S', '--short2', type=str, help='Path to the second short read file in FASTQ format (if available).')

    # Long reads (optional, with flags)
    parser.add_argument('-pb', '--pacbio', type=str, help='Path to PacBio read file in FASTQ format (if available).')
    parser.add_argument('-on', '--nanopore', type=str, help='Path to Oxford Nanopore file in FASTQ format (if available).')

    args = parser.parse_args()

    # Validate if at least one type of reads is provided:
    if not args.short1 and not args.pacbio and not args.nanopore:
        parser.error("You must provide at least one type of reads: either short reads or long reads!")

    return args

def validate_file_exists(file_path, description):
    """Checks if provided input files do exist."""
    
    if file_path and not os.path.isfile(file_path):
        raise FileNotFoundError(f"The {description} file '{file_path}' does not exist!")
    
def validate_file_content(file_path, marker, description):
    """Checks if input files have the correct format."""

    if file_path:
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
            if not first_line.startswith(marker):
                raise ValueError(f"The {description} file '{file_path}' does not have a valid format, for it does not begin with {marker}!")
    
def validate_files(assembly, short1, short2, long_pb, long_on):
    """Condensates the two validation steps into one single step for the user to call."""

    # Validate files' existence:
    validate_file_exists(assembly, 'genome assembly')
    validate_file_exists(short1, 'first short read')
    validate_file_exists(short2, 'second short read')
    validate_file_exists(long_pb, 'pacbio read')
    validate_file_exists(long_on, 'nanopore read')

    # Validade files' format:
    validate_file_content(assembly, '>', 'genome assembly')
    validate_file_content(short1, '@', 'first short read')
    validate_file_content(short2, '@', 'second short read')
    validate_file_content(long_pb, '@', 'pacbio read')
    validate_file_content(long_on, '@', 'nanopore read')

def bwa_index(assembly):
    """Indexes assemblies for further analysis, such as read mapping."""
    
    print("Indexing assembly...")
    
    # Assembly indexation command:
    index_command = ['bwa', 'index', assembly]

    try:
        # Run indexation command:
        subprocess.run(index_command, check=True)
        print("Assembly has been properly indexed.")

    except subprocess.CalledProcessError as error:
        print("An error occurred while running bwa index command: ", error)

def samtools_sort(mapped_read, sorted_read):
    """Sort mapped reads."""
    
    print("Sorting mapped reads...")

    # Defining sort command:
    sort_command = ['samtools', 'sort', mapped_read, '-o', sorted_read]

    try:
        subprocess.run(sort_command, check=True)

    except subprocess.CalledProcessError as error:
        print('An error occurred while running samtools sort command: ', error)

def samtools_merge(output_bam, sorted_pair_bam, sorted_unp_bam):
    """Merges sorted .bam files (paired and unpaired) into a single one."""

    print("Merging bam files...")
    
    # Defining merge command:
    merge_command = ['samtools', 'merge', output_bam, sorted_pair_bam, sorted_unp_bam]

    try:
        # Run command:
        subprocess.run(merge_command, check=True)

    except subprocess.CalledProcessError as error:
        print('An error occurred while running samtools merge command: ', error)

def bwa_mapping_single(assembly, trimmed_read, out_bam):
    """Calling BWA to map single-end reads to the assembly, returning a sorted BAM file."""
    
    print("Mapping reads to assembly...")

    # Define command for paired-end reads:
    bwa_command = ['bwa', 'mem', assembly, trimmed_read, '-t', '10']
    
    # Define temporary file to store BAM output:
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as bam_tmp:
        try:
            # Running BWA command:
            subprocess.run(bwa_command, stdout=bam_tmp, stderr=subprocess.DEVNULL, check=True)
            
        except subprocess.CalledProcessError as error:
            print('An error occurred while running BWA mem command: ', error) 

    # Sort mapped reads:
    samtools_sort(bam_tmp.name, out_bam) # Samtools sort requires os.Pathlike object, not a file!

    # Clean up temporary file:
    os.remove(bam_tmp)

    print("Reads have been properly mapped to the assembly.")

    return out_bam

def bwa_mapping_paired(assembly, trimmed_read1, trimmed_read2, unpaired_reads, out_bam):
    """Calling BWA to map paired-end reads to the assembly, returning a merged BAM file."""
    
    print("Mapping reads to assembly...")

    # Define command for paired-end reads:
    paired_command = ['bwa', 'mem', assembly, trimmed_read1, trimmed_read2, '-t', '10']
    unpaired_command = ['bwa', 'mem', assembly, unpaired_reads, '-t', '10']
    
    # Define temporary files to store paired and unpaired BAM output:
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as paired_bam_tmp, \
    tempfile.NamedTemporaryFile(mode='w', delete=False) as unpaired_bam_tmp:
        try:
            # Running BWA command for paired reads:
            subprocess.run(paired_command, stdout=paired_bam_tmp, stderr=subprocess.DEVNULL, check=True)
            # Running BWA command for unpaired reads:
            subprocess.run(unpaired_command, stdout=unpaired_bam_tmp, stderr=subprocess.DEVNULL, check=True)

        except subprocess.CalledProcessError as error:
            print('An error occurred while running BWA mem command: ', error)

    # Define temporary files' paths to store sorted mapped reads, as samtools itself will create the files:
    sorted_paired_bam = tempfile.NamedTemporaryFile(delete=False).name
    sorted_unpaired_bam = tempfile.NamedTemporaryFile(delete=False).name
    
    # Running samtools sort for paired BAM:
    samtools_sort(paired_bam_tmp.name, sorted_paired_bam)
    # Running samtools sort for unpaired BAM:
    samtools_sort(unpaired_bam_tmp.name, sorted_unpaired_bam)

    # Merging sorted reads:
    samtools_merge(out_bam, sorted_paired_bam, sorted_unpaired_bam)

    # Clean up temporary files:
    os.remove(paired_bam_tmp)
    os.remove(unpaired_bam_tmp)
    os.remove(sorted_paired_bam)
    os.remove(sorted_unpaired_bam)

    print("Reads have been properly mapped to the assembly.")

    return out_bam

def minimap2_pb(assembly, trimmed_long_pb, out_bam):
    """Calls Minimap2 and map PacBio long reads to the assembly, returning a sorted BAM file."""
    
    print("Mapping PacBio reads to assembly...")

    # Define commands for PacBio and Nanopore reads:
    pacbio_command = ['minimap2', '-ax', 'map-pb', assembly, trimmed_long_pb]
    
    # Define temporary file to store SAM PacBio output:
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as sam_pb_tmp:
        try:
            # Running Minimap2 command:
            subprocess.run(pacbio_command, stdout=sam_pb_tmp, stderr=subprocess.DEVNULL, check=True)
            
        except subprocess.CalledProcessError as error:
            print('An error occurred while running minimap2 command: ', error) 

    # Sort mapped reads:
    samtools_sort(sam_pb_tmp.name, out_bam)

    # Clean up temporary file:
    os.remove(sam_pb_tmp)

    print("Reads have been properly mapped to the assembly.")

    return out_bam

def minimap2_on(assembly, trimmed_long_on, out_bam):
    """Calls Minimap2 and map Oxford Nanopore long reads to the assembly, returning a sorted BAM file."""
    
    print("Mapping Oxford Nanopore reads to assembly...")

    # Define commands for PacBio and Nanopore reads:
    nanopore_command = ['minimap2', '-ax', 'map-on', assembly, trimmed_long_on]
    
    # Define temporary file to store SAM PacBio output:
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as sam_on_tmp:
        try:
            # Running Minimap2 command:
            subprocess.run(nanopore_command, stdout=sam_on_tmp, stderr=subprocess.DEVNULL, check=True)
            
        except subprocess.CalledProcessError as error:
            print('An error occurred while running minimap2 command: ', error) 

    # Sort mapped reads:
    samtools_sort(sam_on_tmp.name, out_bam)

    # Clean up temporary file:
    os.remove(sam_on_tmp)

    print("Reads have been properly mapped to the assembly.")

    return out_bam

def single_end_processing(short1, trimmed1, html_out, json_out):
    """Run fastp command to single-end reads, trimming adapters and bases with low quality score."""

    # Defining fastp command:
    fastp_command = [
            'fastp',
            '-i', short1,
            '-o', trimmed1,
            '-h', html_out,
            '-j', json_out,
            '-5', '-3', '-M', '25'
        ]
    
    # Running command:
    try:
        subprocess.run(fastp_command, check=True)
        print("Fastp command excecuted successfully!")
        
    except subprocess.CalledProcessError as error:
        print("An error occurred while running fastp command: ", error)

def paired_end_processing(short1, short2, trimmed1, trimmed2, trimmed_unp1, trimmed_unp2, html_out, json_out):
    """Run fastp command to paired-end reads, trimming adapters and bases with low quality score."""

    # Defining fastp command:
    fastp_command = [
            'fastp',
            '-i', short1,
            '-I', short2,
            '-o', trimmed1,
            '-O', trimmed2,
            '--unpaired1', trimmed_unp1,
            '--unpaired2', trimmed_unp2,
            '-h', html_out,
            '-j', json_out,
            '-5', '-3', '-M', '25', '-c'
        ]
    
    # Running command:
    try:
        subprocess.run(fastp_command, check=True)
        print("Fastp command excecuted successfully!")
    
    except subprocess.CalledProcessError as error:
        print("An error occurred while running fastp command: ", error)

def long_read_processing(long, trimmed, html_out, json_out):
    """Run fastp command to PacBio or Oxford Nanopore reads, trimming bases with low quality score."""

    # Defining fastp command:
    fastp_command = [
            'fastp',
            '-i', long,
            '-o', trimmed,
            '-h', html_out,
            '-j', json_out,
            '-5', '-3', '-M', '25'
        ]
    
    # Running command:
    try:
        subprocess.run(fastp_command, check=True)
        print("Fastp command excecuted successfully!")
        
    except subprocess.CalledProcessError as error:
        print("An error occurred while running fastp command: ", error)

def main():
    args = parse_arguments()

    # Access the parsed arguments:
    assembly = args.assembly
    short1 = args.short1
    short2 = args.short2
    pacbio = args.pacbio
    nanopore = args.nanopore

    # Validate input files:
    try:
        validate_files(assembly, short1, short2, pacbio, nanopore)
    except (FileNotFoundError, ValueError) as error:
        print(error)
        sys.exit(1)

    # State the validated inputs:
    print(f"Assembly: {assembly}")
    if short1:
        print(f"First short read: {short1}")
    if short2:
        print(f"Second short read: {short2}")
    if pacbio:
        print(f"PacBio long read: {pacbio}")
    if nanopore:
        print(f"Oxford Nanopore long read: {nanopore}")

    # Recover basenames of assembly and reads' files in order to rename downstream files:
    assembly_base_name = os.path.splitext(os.path.basename(assembly))[0]

    if short1:
        read_base_name1 = os.path.splitext(os.path.basename(short1))[0]
        
        if short2:
            read_base_name2 = os.path.splitext(os.path.basename(short2))[0]
            read_base_name = read_base_name1.split('_')[0] # To recover basename without the _1 or _2 from the fastq
    
    if pacbio:
        pacbio_base_name = os.path.splitext(os.path.basename(pacbio))[0]

    if nanopore:
        nanopore_base_name = os.path.splitext(os.path.basename(nanopore))[0]

    # Read processing:
    if short1:
        if short2:
            print("Processing paired-end short reads...")

            # Define paths for the output files:
            trimmed_paired1 = f'{read_base_name1}_trimmed.fastq'
            trimmed_paired2 = f'{read_base_name2}_trimmed.fastq'
            trimmed_unp1 = tempfile.NamedTemporaryFile(delete=False).name
            trimmed_unp2 = tempfile.NamedTemporaryFile(delete=False).name
            paired_html = f'{read_base_name}_qc.html'
            paired_json = f'{read_base_name}_qc.json'
                
            # Calling read processing function:
            paired_end_processing(short1, short2, trimmed_paired1, trimmed_paired2, trimmed_unp1, trimmed_unp2, paired_html, paired_json)

            # Merging unpaired trimmed reads:
            trimmed_unp = f'{read_base_name}_unp_trimmed.fastq'

            with open(trimmed_unp, 'wb') as trimmed_unp_file: # b as it is a binary file
                for file in [trimmed_unp1, trimmed_unp2]:
                    with open(file, 'rb') as fd:
                        shutil.copyfileobj(fd, trimmed_unp_file)
            
            print("Paired-end reads have been correctly processed.")
        
        else:
            print("Processing single-end short reads...")
            
            # Defining paths for output files:
            trimmed_single = f'{read_base_name1}_trimmed.fastq'
            single_html = f'{read_base_name1}_qc.html'
            single_json = f'{read_base_name1}_qc.json'
                
            single_end_processing(short1, trimmed_single, single_html, single_json)
            print("Single-end reads have been correctly processed.")

    if pacbio:
        print("Processing Pacbio reads...")

        # Paths for output files:
        trimmed_long_pb = f'{pacbio_base_name}_trimmed.fastq'
        pacbio_html = f'{pacbio_base_name}_qc.html'
        pacbio_json = f'{pacbio_base_name}_qc.json'
                
        long_read_processing(pacbio, trimmed_long_pb, pacbio_html, pacbio_json)
        print("PacBio reads have been properly processed.")

    if nanopore:
        print("Processing Oxford Nanopore reads...")

        # Paths for output files:
        trimmed_long_on = f'{nanopore_base_name}_trimmed.fastq'
        nanopore_html = f'{nanopore_base_name}_qc.html'
        nanopore_json = f'{nanopore_base_name}_qc.json'
                
        long_read_processing(nanopore, trimmed_long_on, nanopore_html, nanopore_json)
        print("Oxford Nanopore reads have been properly processed.")

    # Read mapping:
    if short1:
        if short2:
            # Assembly indexation:
            bwa_index(assembly)

            # Read mapping:
            paired_bam = f"{read_base_name}_{assembly_base_name}.bam"

            bwa_mapping_paired(assembly, trimmed_paired1, trimmed_paired2, trimmed_unp, paired_bam)
        
        # Assembly indexation:
        bwa_index(assembly)

        # Read mapping:
        single_bam = f"{read_base_name1}_{assembly_base_name}.bam"
        
        bwa_mapping_single(assembly, trimmed_single, single_bam)

    if pacbio:
        # Assembly indexation and read mapping:
        pacbio_bam = f"{pacbio_base_name}_{assembly_base_name}.bam"

        minimap2_pb(assembly, trimmed_long_pb, pacbio_bam)

    if nanopore:
        # Assembly indexation and read mapping:
        nanopore_bam = f"{nanopore_base_name}_{assembly_base_name}.bam"
        
        minimap2_on(assembly, trimmed_long_on, nanopore_bam)

if __name__ == "__main__":
    main()
