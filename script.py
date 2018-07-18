import pysam
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#collect all BAM files from a given directory
#return a list of all BAM files
def scan_files(directory):
	bam_files = []
	for filename in os.listdir(directory):
		if filename.endswith(".bam"):
				bam_files.append(os.path.join(directory, filename))
	return (bam_files)




#Using a list of BAM files, count the number of each base at a specific location in each BAM
#Print counts and percentages of each base at the given position from each BAM file
#Returns relevant information to use in print file
def get_counts(files):

	master_list = []
	
	ref_seq = pysam.FastaFile("/home/upload/pysam_project/hg19.fa")	

	#Ask for chromosome and base position information
	chromosome = raw_input("Chromsome Number: ")
	if chromosome is 'x':
		chromosome = 'X'
	base_position = raw_input("Base Position: ")

	#Quality filter input and input validation
	quality_filter = 0
	#response = raw_input("Use a quality filter? (Y/N): ")	
	#while ((response != 'Y') and (response != 'y') and (response != 'N') and (response != 'n')):
		#response = raw_input("Use a quality filter? (Y/N): ")
	#if (response == 'Y') or (response == 'y'):
		#quality_filter = raw_input("Quality Filter: ")
		#while (int(quality_filter) < 0) or (int(quality_filter) > 60):	
			#quality_filter = raw_input("Quality Filter: ") 
	
	#Search reference genome for the reference base
	ref_base = ref_seq.fetch("chr" + str(chromosome), int(base_position) - 1, int(base_position))
	print("Reference Base: %s" % (ref_base))

	#Instantiate lists that will store the percentage of reads of each base for each BAM file
        a_list = []
        c_list = []
        g_list = []
        t_list = []
	
	#Loop through every BAM file given
	for file in files:
		
		coverage = 0		

		a_count = 0
        	c_count = 0
        	g_count = 0
        	t_count = 0
		
		a_percent = 0
        	c_percent = 0
	    	g_percent = 0
	        t_percent = 0

		is_covered = True

		samfile = pysam.AlignmentFile(file, "rb")

		#For the single base position given, loop through every read
		for pileupcolumn in samfile.pileup(str(chromosome), int(base_position) - 1, int(base_position), truncate = True, max_depth = 50000):
			for pileupread in pileupcolumn.pileups:

				###PROBLEM: This item can't get the index [pileupread.query_position] - NoneType error
                                #print(type(pileupread.alignment))
                                #print(type(pileupread.query_position))
                                #print(pileupread.query_position)
				#print(type(pileupread.alignment.query_qualities))
                                #print(type(pileupread.alignment.query_sequence))
				#print dir(pileupread.alignment)
				#print pileupread.alignment.query_alignment_qualities
                          

                                #if pileupread.alignment.query_qualities[pileupread.query_position] < quality_filter:
                                        #print("Base quality low")
					#continue

				#Sort by base, count each read by base
				if pileupread.alignment.query_sequence[pileupread.query_position] is "C":
					c_count += 1
				elif pileupread.alignment.query_sequence[pileupread.query_position] is "A":
					a_count += 1
				elif pileupread.alignment.query_sequence[pileupread.query_position] is "T":
					t_count += 1
				elif pileupread.alignment.query_sequence[pileupread.query_position] is "G":
					g_count += 1
				else: 
					pass
				#Report the coverage of the position in the given BAM
				coverage = pileupcolumn.n
		
		#If there are NO reads at this location in a given BAM file, do not report 
		if (coverage == 0):
			is_covered = False
		
		#If reads are present, compute the percentage of the total number of reads that is each base	
		else:
			a_percent = float(a_count * 100.0 / coverage)
			c_percent = float(c_count * 100.0 / coverage)
			g_percent = float(g_count * 100.0 / coverage)
			t_percent = float(t_count * 100.0 / coverage)
		
		#Print some relevant information for the console
		print("\nBAM Name: %s" % (file))
		print("Coverage: %d" % (coverage)) 
		print("Number of A bases: %d (%f)" % (a_count, a_percent)) 
		print("Number of C bases: %d (%f)" % (c_count, c_percent)) 
		print("Number of G bases: %d (%f)" % (g_count, g_percent)) 
		print("Number of T bases: %d (%f)" % (t_count, t_percent))

		#For each BAM, add to the base-level lists the percentage of the time that base shows up
		a_list.append(a_percent)
		c_list.append(c_percent)
		g_list.append(g_percent)
		t_list.append(t_percent)	

		samfile.close()
	
	#master_list now contains the percent of times each base appears in each BAM file
	master_list.extend([a_list, c_list, g_list, t_list])

	return(master_list, chromosome, base_position, is_covered, ref_base)


#print_graph uses relevant information from get_counts and uses it to form a histogram
#NO return, but saves PNG files using the chromosome number and base position to label
def print_graph(master_list, chromosome, base_position, is_covered, ref_base):
	
	
	gs = gridspec.GridSpec(3,2)
	fig, axs = plt.subplots(3, 2, sharey = True, tight_layout = True)

	ax1 = plt.subplot(gs[0, :])
	ax2 = plt.subplot(gs[1, 0], sharey = ax1)
	ax3 = plt.subplot(gs[1, 1], sharey = ax1)
	ax4 = plt.subplot(gs[2, 0], sharey = ax1)
	ax5 = plt.subplot(gs[2, 1], sharey = ax1)

	axes = [ax1, ax2, ax3, ax4, ax5]

	bin_array = []
	x = 0
	while x <= 100:
		bin_array.append(x)
		x += 5
	bin_array2 = []
        x = 0.0
        while x <= 5:
                bin_array2.append(x)
                x += 0.05
	
	#Given the reference base, properly set the colors and the legend labels. Reference base will always be GREY
	if ref_base is 'A':
        	ax1.hist([master_list[0], master_list[1], master_list[2], master_list[3]], color = ['gray', 'blue', 'orange', 'red'], bins = bin_array, label = ['Reference = A', 'C', 'G', 'T'])
		ax2.hist(master_list[0], color = 'grey', bins = bin_array2, label = 'A')
                ax3.hist(master_list[1], color = 'blue', bins = bin_array2, label = 'C')
                ax4.hist(master_list[2], color = 'orange', bins = bin_array2, label = 'G')
                ax5.hist(master_list[3], color = 'red', bins = bin_array2, label = 'T')	
	elif ref_base is 'C':
        	ax1.hist([master_list[1], master_list[0], master_list[2], master_list[3]], color = ['gray', 'green', 'orange', 'red'], bins = bin_array, label = ['Reference = C', 'A', 'G', 'T'])
		ax2.hist(master_list[1], color = 'grey', bins = bin_array2, label = 'C')
                ax3.hist(master_list[0], color = 'green', bins = bin_array2, label = 'A')
                ax4.hist(master_list[2], color = 'orange', bins = bin_array2, label = 'G')
                ax5.hist(master_list[3], color = 'red', bins = bin_array2, label = 'T')
	elif ref_base is 'G':
        	ax1.hist([master_list[2], master_list[0], master_list[1], master_list[3]], color = ['gray', 'green', 'blue', 'red'], bins = bin_array, label = ['Reference = G', 'A', 'C', 'T'])
		ax2.hist(master_list[2], color = 'grey', bins = bin_array2, label = 'G')
                ax3.hist(master_list[0], color = 'green', bins = bin_array2, label = 'A')
                ax4.hist(master_list[1], color = 'blue', bins = bin_array2, label = 'C')
                ax5.hist(master_list[3], color = 'red', bins = bin_array2, label = 'T')
	elif ref_base is 'T':
        	ax1.hist([master_list[3], master_list[0], master_list[1], master_list[2]], color = ['gray', 'green', 'blue', 'orange'], bins = bin_array, label = ['Reference = T', 'A', 'C', 'G'])
		ax2.hist(master_list[3], color = 'grey', bins = bin_array2, label = 'T')
                ax3.hist(master_list[0], color = 'green', bins = bin_array2, label = 'A')
                ax4.hist(master_list[1], color = 'blue', bins = bin_array2, label = 'C')
                ax5.hist(master_list[2], color = 'orange', bins = bin_array2, label = 'G')	
	
	for ax in axes:
		ax.legend(loc = 'upper center', fontsize = 'small')
		ax.grid(True)
		ax.minorticks_on()
		
	fig.text(0.5, 0.01, '% of base per read', ha = 'center')
	fig.text(0.0099, 0.5, 'number of bams', va = 'center', rotation = 'vertical')	
	fig.set_size_inches(11, 8.5)	

	#graph title is in the format 'chr#_bp#_dist'
        title = ("chr%s_bp%s_dist" % (chromosome, base_position))
	fig.suptitle(title, y = .995)

	#Location at which to save the histogram, using the title of the graph
	save_loc = ("/home/upload/pysam_project/%s.png" % title)
	
	#Decide whether or not to draw the graph
	if (is_covered):
	        print("Graph saved at %s" % save_loc)
        	plt.savefig(save_loc)
	else:
        	print("No coverage. No graph produced.")
	return()


#Main routine to run given a directory input
def main(directory):
	master_list, chromosome, base_position, is_covered, ref_base = get_counts(scan_files(directory))
	print_graph(master_list, chromosome, base_position, is_covered, ref_base)

main("/home/upload/pysam_project/myeloid_bams")
