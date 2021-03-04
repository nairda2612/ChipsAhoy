# ChIPs Ahoy: An automatic pipeline to process and analyze ChIP-seq data

Authors: 

Andrea Fernández Veloso, Adrián Perera Bonaño and Emma Serrano Pérez
University of Seville

## 1. Summary

## 2. Dependencies

## 3. Input

## 4. Usage

### a. Parameters read

### b. Workspace generation

### c. Index building

### d. Sample load

### e. Sample processing

#### i. Quality control

#### ii. Reads alignment

### f. Peaks calling

### g. Regulome determination and Gene Set Enrichment Analysis

### h. Binding motifs search

## 7. Case study


## 1. Summary

RSUUUMEN AL FINAL DECIR 1 HAY 2 FICEROS DE BASH Y 1 DE R
PERMITE ANALIZAR NUMERO DE MUESTRAS EN PLAN DE TF Q KIERAS PEEERO SIEMPRE Q SE COMPAREN TODOS CON LOS MISMOS (YA SEAN 1 O MAS, CONTROLES

ESTO ES PA FACTTS Y PA ATHALIANA

## 2. Dependencies

MI PARTE PAL FINALLL

## 3. Input

The pipeline chipsahoy requires only one parameter. This parameter, the so-called “parameters file” is a TXT file with the parameters needed by the pipeline. To know more about the file content, go to section 4.a. 

When no parameters are passed to the script, a help message is printed out on the screen to show to the user how to make the pipeline work. Also, if the user tries to pass more than one parameter to the script, a message is printed out on the screen to notify and indicate to the user how to access the usage and run the script properly. 

## 4. Usage

To run the pipeline, write on the command line:

./chipsahoy /path_to_parameters_file/parameters_file.txt

To store the output of the whole script on a file, including errors, write on the command line: 

./chipsahoy /path_to_parameters_file/parameters_file.txt | tee output_file_name.txt

### a. Parameters read

The parameters passed to the function need to be written on a TXT file, the “parameters file”. An example of a parameters file can be found on the directory test. When the pipeline is run with this file, a case study is generated. To see the case study, go to section 7.

The parameters file includes the following information:

-installation_directory: The path to the directory where the ChIPs Ahoy repository shall be installed.

-working_directory: The path to the directory where will be generated the experiment directory.

-experiment_name: The name of the experiment directory, which will be located in the working directory. Every file created by the pipeline will be saved in a subdirectory of the experiment directory.

-number_rep_tf: Number of replicas for every transcription factor. This pipeline is designed to work with the same number of replicas for each transcription factor.

-number_rep_control: Number of replicas for every control. The information given above for the parameter “Number of replicas for every transcription factor” is applicable. 

-path_genome: The path to the file to be used as reference genome, in FASTA format.

-path_annotation: The path to the file to be used as reference genome annotation, in GTF format.

-number_TF: Number of transcription factors. This pipeline can process more than one transcription factors at once.

-number_control: Number of control. This pipeline can use more than one sample as a control. Know that every control sample given to the pipeline will be used for every transcription factor.

-path_chip_X: The path to every sample file, in GZ format. Add to the parameters file, as many lines of this parameters as number chip samples (being X the number of the sample). To use more than one replica for any sample, write the path for every replica in the same line, all of them separated by a space.

-path_control_X: The path to every control sample. The information given above above for the parameter “Path for every chip sample” is applicable.

-up1: Number of base pairs upstream the narrow peaks in which will be defined the promoters and the transcription start site (tss).

-down1: Number of base pairs downstream the narrow peaks in which will be defined the promoters and the transcription start site (tss).

-up2: Number of base pairs upstream the summits in which will be defined the promoters and the transcription start site (tss).

-down2: Number of bas epairs downstram the summits in which will be defined the promoters and the transcription start site (tss).

-pvalueCutoffgo: p-valor used in the Gene Ontology Gene Set Enrichment Analysis.

-pvalueCutoffkegg: p-valor used in the KEGG Pathways Gene Set Enrichment Analysis.

-paired: Write 0 for single-end sequenced samples and 1 for paired-end sequenced samples.

When a parameters file is passed to the pipeline chipsahoy, firstly, the pipeline stores the lines of the parameters file in different variables, and prints out on the screen every the variable for each parameter.

To read de chip and control samples, as their number is changeable, two identical while loops are used. These loops create an array, CHIP or CONTROL, in which the routes for every sample are stored. To break every line read in its different elements (replicas, separated by a space), the array is passed to string, broken, and then every element is stored, separately, in a new array of the same name, CHIP or CONTROL. The lengths of the arrays and string created in every step are printed out on the screen. Thanks to these loops, chipsahoy can process, at once, more than one chip and control samples, each one with more than one replica.

### b. Workspace generation

To create the workspace, the route of the working directory is run. Once there, the experiment directory is created with the name given on the parameter file, then accessed, and there, five directories are created: genome, annotation, samples, results and scripts.

After that, using their paths, written on the parameters file, the genome and annotation files are copied to the directories genome and annotation, under the names “genome.fa” and “annotation.gtf”, respectively. 

### c. Index building

To build the index, the pipeline leaves the annotation directory and accesses the genome directory.

Using the function bowtie2-build, which receives the genome.fa and a word to use as suffix of the index files, the pipeline builds an index of the reference genome, and saves it on the genome directory, under the suffix .index.

### d. Sample load

To generate the samples directories, first, the pipeline leaves the genome directory and accesses the samples directory. There, two subdirectories, chip and control, are created. Since the number of chips and controls (paired-end sequenced or not) and their replicas can change, two identical while loops, one for chip samples and the other for control samples, access those directories and create as many directories in the chip and control directories, as the numbers of each type of samples, with their replicas, are specified in the parameters file. These directories are named chip_sample_X and control_sample_X, being X the number of the sample. In these same loops, the elements from the arrays CHIP and CONTROL, which are the paths to each sample, are copied, in order, to their respectives directories. 

Every sample file will be copied under the name chip_X_Y.gz or control_Xl_Y.gz, being X the number of the sample, and Y the number of the replica. Also, if the samples were paired-end sequenced, they will be saved under the name chip_X_Y_Z.gz or control_X_Y_Z.gz, being Z = 1 and 2, because of the two files obtained from the paired-end sequencing.

### e. Sample processing

To process samples, the pipeline leaves the last control sample directory created, then leaves the control directory, and the samples directory, and accesses the results directory.

Then, two identical loops can be found, one for chip samples, and the other for control samples. These loops go through every replica, of every chip or control sample (also takes into account if the sequencing was paired-end or not) and generates the variables SAMPLENAME1 and SAMPLENAME2, with the names of the corresponding sample files, which include information of the number of sample (X), replica (Y) and sequencing method followed (Z), as previously explaind. Then, these loops access the ChipsAhoy directory from the installation directory, to run the bash script “sample_processing” using the command bash. The arguments passed to this script are the following: path to the pertinent chip or control directory from the working directory, the variables SAMPLENAME1, SAMPLENAME2, PAIRED (depending if the samples are paired-end or not), and chip_X_Y or control_X_Y, which will be the new sample name after the reads are mapped to the reference genome. To know more about the reads alignment, go to section 4.e.ii.

In the sample processing script, a loop to make de quality control and reads alignment can be found. This loop takes into account if the samples are paired-end or not. 

#### i. Quality control

Once the pipeline runs the sample_processing bash script, a quality control is made for the sample passed to the script, using the function fastqc, which receives the name of the sample (SAMPLE1) or samples (if their paired-end, SAMPLE1 and SAMPLE2) . As the last folder accessed was the pertinent sample directory, for each sample, the quality control results will be saved in that same directory.

#### ii. Reads alignment

The next step is to align the reads to the index of the reference genome, previously built. This way, the different zones of the genome where reads accumulate, which probably correspond to binding regions of the transcription factor, can be identified. To make the alignment, the bowtie2 function is used. The parameters passed to this function are the following: 

-x    /path/to/index_of_reference_genome

-U or -1 -2   , depending if the samples are single or paired-end, respectively. This argument receives the GZ file of the sample being processed.

-S   output name of the SAM file generated. This name corresponds to the last parameter passed to the while loop that runs the sample_processing bash script. To remember this, go to section 2e.

During the SAM file generation, the code will be printed on a mapping file on TXT format, and a loading bar will be printed out on the screen to let the user know that the SAM is being created. The loading bar adds a lightning symbol per second that the SAM file is being created.

After the SAM file is created, since it is a very heavy file, it is convert to a BAM file under the same name, given by the variable SAMPLENAME, as explained before. To make this conversion, the function samtools sort is used. The parameters passed to this function are the output BAM file name (SAMPLENAME.bam) and the SAM file. Then, the SAM file is removed.

Also, to ease the access to the mapping information to the user, from the mapping file created previously, the lines that correspond to the number of reads and the overall alignment percentage, are saved on a new TXT file created on the results directory, the mapping values file. This way, when all samples are processed, the mapping values file will have the mapping information for every sample. The mapping file is removed after the relevant information is passed to the mapping values file.

After this, the sample_processing script is done, so the chipsahoy pipeline continues to be run. Only when every sample is processed, the sample processing loop that runs the sample_processing bash script finishes. It is important to remember that, by this point, we’re still in the results directory.

We recommend that the user visualises the alignment made on an interactive tool, like Integrative Genomic Viewer (IGV).

### f. Peaks calling

Even though the reads are aligned on the genome, the user cannot know, only with that information, which are the possible target genes for the transcription factors analyzed. To check if the result of the mapping is significant, a statistical analysis is made. Using a contrast of hypothesis, the error probability can be determined. 

To call peaks, a sliding window goes through the genome, counting reads on the window defined, for chips and control samples. Then, the reads of each pair of samples are compared, and the statistic test previously explained, determines a fold-change and a p-value. This values represent the difference in the number of reads for the chip and control sample, and how significant is that difference, respectively. It is important to know that, since, this process is a multiple testing, a correction of the p-value, called q-value or FDR, needs to be made. Therefore, peaks are regions in which there is a significant difference of reads between the chip and the control samples.

To guarantee that this code is run for every chip and control samples, a double while loop that goes through every chip and every control sample was configured. To call peaks, the function macs2 with the command callpeak is used. The arguments passed to this function are the following:

-t   path to the BAM files for the replicas of the chip sample.

-c   path to the BAM files for the replicas of the control sample.

-f   BAM format, which is the format of the input files.

-n   name of the output file, which is the experiment_name_X_Y, being X and Y the number of chip and control samples, respectively, that are being compared.

-- outdir to specify the output directory, which is still the results directory.

After the peaks calling, the output files can be found on the results directory. Of all the files created, the ones of special interest in this work are the .narrowPeak and .bed files, with the significant narrow peaks and summits, respectively.

We recommend that the user visualises on an interactive tool, like IGV, the result of the alignment, and also the peaks called.

### g. Regulome determination and Gene Set Enrichment Analysis

After the peaks are called, for the next steps the pipeline needs to run an R script for every narrowPeak and bed file. Before that, these files are stored in two arrays, NARROWPEAK and SUMMITSBED, respectively.

This way, a while loop goes through every file with results, narrow peaks and summits, and runs de R script created to determine the regulome and make the Gene Set Enrichment Analysis. The script if run with the command Rscript, and the arguments passed, in addition to the pair of peaks file, the parameters defined in the parameters text up1, down1, up2, down2, pvalueCutoffgo, pvalueCutoffkegg an identification for every sample, which is  experiment_name_X_Y, being X and Y the number of chip and control samples, respectively, that were compared.

The first thing that the R script, called peaks_script.R does, is to store the arguments passed on the bash script on variables. Then, the script loads the necessary packages for the analysis. To know the packages needed, go to section 2. 

Then, the script gets the genetic information from the organism, in this case Arabidopsis thaliana, and gets de universe, in this case, the genes id of the whole organism. After this, the script reads the different files, and, for the narrow peaks, plots the peak coverage.

The next step is to define the promoters for narrow peaks and summits, using the base pair up and downstream given in the parameters file, and to annotate those peaks. A plot of the annotation is generated for each peaks file.

Then, the annotation is converted to data frame. Now, from each data frame, the regulome is determined, simply by extracting the elements whose annotation is “Promoter”. The genes ID of those elements is the regulome, which the pipeline safes on a TXT file.

Once the regulome is determined, 

### h. Binding motifs search

## 7. Case study
