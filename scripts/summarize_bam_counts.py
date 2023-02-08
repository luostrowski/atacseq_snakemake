
# Read in read_count.txt files for each bam subset per sample
# assumes directory strucutre of bowtie2 with multiple subset subdirs

# Import libraries
import optparse
import os
import subprocess

# Initialize command line options for input directory and output file
p = optparse.OptionParser()
p.add_option("-d", action = "store", dest = "directory")
p.add_option("-o", action = "store", dest = "outfile")

opts,args=p.parse_args()
in_dir=opts.directory
outfile=opts.outfile

# open outfile for writing
fhw = open(outfile, "w+")
# write colnames to outfile
fhw.write('Sample_Name' + '\t' + 'Sorted_Bam' + '\t' + 'noMT_Bam' + '\t' + 'rmDups_Bam' + '\t' + 'Mapq30_Bam' + '\t' + 'Proper_Pair_Bam' + '\n')

# Loop through each subdirectory in input directory,
# for each subdirectory, create a list to hold the sample names and read numbers
# Within each subdirectory, capture sample name and read number for each file ending in *.read_count.txt
# Store subdirectory results in a master list
master_outs = []
for subdir in sorted(os.listdir(in_dir)):
    newdir = in_dir + subdir
    list_of_sample_metrics = []
    for file in sorted(os.listdir(newdir)):
        filename = os.fsdecode(file)
        if filename.endswith(".read_count.txt"):
            # capture sample name from filename
            sample_name = filename.split('.')[0]

            # capture contents of read_count.txt file
            #   use .strip() to remove newline character
            #   use .decode('ascii') to convert the returned bytes to a string
            read_num = subprocess.check_output("cat " + newdir + '/' + filename, shell = True).strip().decode('ascii')

            # store the sample metrics as a tuple
            sample_metrics = (sample_name, read_num)
            # add this tuple to the list of sample metrics for the subdirectory
            list_of_sample_metrics.append(sample_metrics)
    # add the subdirectory list of tuples to the master list of output from all subdirectories
    master_outs.append(list_of_sample_metrics)

# master_outs is a list of lists of tuples
# master_outs has a length of 5 (5 subdirectories from bowtie2)
# each element of master_outs is a list of length(n), where n is the number of samples
# each tuple has two elements: sample name and read number

# Loop through master outs and export desired information to produce a table

# the master loop needs to be for how many samples there are
for samp in range(len(master_outs[0])):
    line_to_print = []
    for subdir in range(len(master_outs)):
        if subdir == 0:
            sample_name = master_outs[subdir][samp][0]
            bam_reads = master_outs[subdir][samp][1]
            line_to_print.append(sample_name + '\t' + bam_reads + '\t')
        if 0 < subdir < 4:
            bam_reads = master_outs[subdir][samp][1]
            line_to_print.append(bam_reads + '\t')
        if subdir == 4:
            bam_reads = master_outs[subdir][samp][1]
            line_to_print.append(bam_reads + '\n')
    for line in line_to_print:
        fhw.write(line)

fhw.close()
