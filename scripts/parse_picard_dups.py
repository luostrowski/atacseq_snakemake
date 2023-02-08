
# Use dup_metrics.txt files from Picard MarkDuplicates to capture percent duplication

# Import libraries
import optparse
import os

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
fhw.write('Sample_Name' + '\t' + 'Percent_Duplication' + '\n')

# Loop through files in input directory, and work on files ending in *report.txt
#   to grab lines of interest and record the metric to an output file
for file in sorted(os.listdir(in_dir)):
    filename = os.fsdecode(file)
    printList = []
    if filename.endswith(".dup_metrics.txt"):
        # capture sample name from filename
        sample_name = filename.split('.')[0]
        path_file=in_dir + '/' + filename
        theFile = open(path_file,'r')
        FILE = theFile.readlines()
        theFile.close()
        for line in FILE:
            # Obtain sample name from bam file)
            if (line.startswith("Unknown")):
                percent_dup = line.split()[9]
        printList.append(sample_name + '\t' + percent_dup + '\n')
    for item in printList:
        fhw.write(item)

fhw.close()
