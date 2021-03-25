report: "report/workflow.rst"

rule all:
	input:
		["phy3_test_metadata_plots.pdf",
		"phy3_test.cfml.pdf",
		#"phy2_test_treeWAS_plots.pdf",
		"phy3_test_HGT_candidates_summary.pdf"]

# todo
# BASH problem: strickt bash? add R and conda to PATH!
# 1. replace inputs directory with a config file
# 2. tree of dia hits and correspoding OG!
# (hits of protein x and og of protein x -> tree!)
# 3. cfml for removing recombinant regions from of-msa -> treetime
# 4. ancestoral state prediction (gene gain vs loss!)
# 5. reconstruct of ancestoral seq for each OG and run diamond

# 6. recombinant dna position mapping
# - instead of search of an exact position range
# - check if the protein is in the recombinant region or not!


# OF MSA = get proteins for OGs + concat + additional alignment???
# take species MSA - remove gaps - get proteins?


# OrthoFinder performs light trimming of the MSA to prevent overly long runtimes
# & RAM usage caused by very long, gappy alignemnts. A column is trimmed from
# the alignment if is it greater than 90% gaps and provided two conditions are met.
# 1. The length of the trimmed alignment cannot go below 500 AA
# 2. No more than 25% of non-gap characters can be removed from the alignment.
# If either of these conditions are not met then the threshold for the percentage
# of gaps in removed columns is progressively #increased beyond 90% until both conditions are met.
# The trimming can be turned off using the option "-z".

# Part 0: metadata
#rule prepare_metadata:
#	input:
#		"inputs/{sample}_metadata"
#	output:
#		"inputs/{sample}_metadata/{sample}_cleaned_data.csv"
#	shell:
#		"touch {output}"

rule metadata_vis:
	input:
		"inputs/{sample}_metadata/{sample}_cleaned_data.csv"
	output:
		report("{sample}_metadata_plots.pdf",caption="report/metadata.rst",category="Metadata")
	shell:
		"touch {output}"

# Part 1: DNA data
rule parse_GBs:
	input:
		"inputs/panX/{sample}"
		#"inputs/{sample}_directory_with_proteins"
	output:
		protGB="outputs/GENOME_PROTEIN_DB_{sample}.pkl"
	conda:
		"env/biopython_env.yml"
	notebook:
		"scripts/gb2pkl.ipynb"

rule OrthoFinderAll:
# Initial Orthofinder run with all assemblies and all reference NCBI assemblies
	input:
		protGB="outputs/GENOME_PROTEIN_DB_{sample}.pkl"
	output:
		dir=directory("outputs/{sample}_OrthoFinderAll_output_directory")
	shell:
		"mkdir {output}"

rule RhierBAPS:
	input:
		tree="inputs/panX/{sample}/strain_tree.nwk",
		msa="inputs/panX/{sample}/SNP_whole_matrix_names_corrected.aln"
	output:
		plot=report("{sample}_rhierBAPS_subclusters.pdf",
		caption="report/rhierBAPS.rst", category="Population structure"),
		table="{sample}_rhierBAPS_subclusters.csv"
	shell:
		"touch {output.plot};"
		"cp DEBUG/rso_test_rhierBAPS_subclusters.csv {output.table};"
		"cp DEBUG/rso_test_rhierBAPS_subclusters.pdf {output.plot}"
	#script:
	#	"scripts/rhierBAPS.R"

rule get_subtrees_RhierBAPS:
	input:
		table = rules.RhierBAPS.output.table,
		vis = rules.RhierBAPS.output.plot
	output:
		dir=directory("outputs/{sample}_phylotypes_clusters")
	shell:
		"mkdir {output}"

# Can generate initial tree using OrthoFinderAll single copy orthologs and/or refer to PanX run
# Identify individual phylotype membership for each strain (add to metadata?)
# Repeat OrthoFinder step for individual phylotypes

# Part 2: orthofinder and core genome
rule divide_proteins_by_phylotypes:
	input:
		rules.get_subtrees_RhierBAPS.output.dir,
		tree="inputs/panX/{sample}/strain_tree.nwk",
		all_genomes="inputs/{sample}_directory_with_proteins"
	output:
		dir=directory("outputs/{sample}_directory_with_proteins_ordered_by_phylotype")
	shell:
		"mkdir {output}"

rule OrthoFinderPhylotype:
	input:
		rules.divide_proteins_by_phylotypes.output.dir
	output:
		dir=directory("outputs/{sample}_OrthoFinder_output_directory")
	shell:
		"mkdir {output}"

rule get_Xenologs:
	input:
		rules.OrthoFinderPhylotype.output.dir,
		phy=rules.OrthoFinderAll.output.dir
	output:
		"outputs/{sample}_xenologs_summary.txt"
	shell:
		"touch {output}"

# use OG rthogroups_for_concatenated_alignment.txt to get the nucleotides
# use SpeciesTree/SpeciesTree_rooted.txt as tree
# use MultipleSequenceAlignments/SpeciesTreeAlignment.fa as msa
#rule nucl_for_concatenated_alignment:
#	input:
#		rules.OrthoFinderPhylotype.output.dir
#	output:
#		"outputs/{sample}_nucl_for_pmsa.fasta"
#	shell:
#		"touch {output}"

#rule pal2nal:
#	input:
#		msa=rules.OrthoFinderPhylotype.output.dir,
#		nucl=rules.nucl_for_concatenated_alignment.output
#	output:
#		"outputs/{sample}_nucl.aln"
#	shell:
		#"touch {output}"
		#'og=\"/Users/devseeva/Desktop/work/rso/OrthoFinder/Results_Dec03/Species_Tree/Orthogroups_for_concatenated_alignment.txt\";'
		#"cp $og {output}"

# extract single copy orthologs (prot) alignments
# extract corresponding nuc sequences using ID
# generate nucleotide alignments with pal2nal (with reference to pro alignment and corresponding nuc fasta seqs)
# concatenate all single copy nucleotide alignments (retain positional info?) >> e.g. phylotype1_msa
# remove invariant sites from msa and generate initial tree with raxml >> phylotype1_snps + phylotype1_snp_tree
# run clonalframeml using phylotype1_msa + phylotype1_snp_tree

rule assign_OGs_to_time:
	input:
		ogs=rules.OrthoFinderPhylotype.output.dir,
		tree="outputs/{sample}_treetime"
	output:
		"outputs/{sample}_dated_OGs.csv"
	shell:
		"touch {output}"

# Part 3: clonalFrameML
# recombination prediction from tree and MSA (for now using PanX results)

# the rule will be potentially replaced by the OrthoFinder results
### temporal_panX_names_correction
rule panX:
	input:
		"inputs/panX/{sample}/SNP_whole_matrix.aln"
	output:
		"inputs/panX/{sample}/SNP_whole_matrix_names_corrected.aln"
	shell:
		"sed 's/ <unknown description>//g' {input} > {output}"


## TODO for phy2!!! ERROR: FASTA_to_nucleotide(): unsupported base K in sequence 5 (ASM359060v) position 7737
### orthoF_msa_phy2RSO/Results_Mar10
rule OGs_to_MSA:
	input:
		tree_in=rules.OrthoFinderPhylotype.output.dir,
		core_ogs="/Users/devseeva/Desktop/work/rso/OrthoFinder/orthoF_msa_{sample}/Species_Tree/Orthogroups_for_concatenated_alignment.txt",
		ogs2prot="/Users/devseeva/Desktop/work/rso/OrthoFinder/orthoF_msa_{sample}/Orthogroups/Orthogroups.tsv",
		path2aln="/Users/devseeva/Desktop/work/rso/OrthoFinder/orthoF_msa_{sample}/MultipleSequenceAlignments/",
		prot_db="outputs/GENOME_PROTEIN_DB_{sample}.pkl"
	output:
		msa_fasta="{sample}_concat.fasta",
		msa_phylip="{sample}_concat.phylip",
		msa_xmfa="{sample}_concat.xmfa"
	params:
		path2OF="/Users/devseeva/Desktop/work/rso/OrthoFinder/"
	conda:
		"env/biopython_env.yml"
	notebook:
		"scripts/ogs2msa.ipynb"

#rule Modeltest_NG:
#	input:
#		msa="{sample}_concat.fasta"
#	output:
#		dir=directory("outputs/{sample}_modeltest_1")
#	conda:
#		"env/ml_tree_env.yml"
#	shell:
#		"mkdir {output.dir};"
#		"modeltest-ng  -i {input.msa};"
#		"mv {input.msa}.* {output.dir}"

rule PhyML:
	input:
		msa="{sample}_concat.phylip",
		# TODO model=rules.Modeltest_NG.output.dir
		#modet_test_log="outputs/{sample}_modeltest_1/{sample}_concat.fasta.log"
		#tmp="/Users/devseeva/Desktop/work/sm_workflow/snakefiles/scripts/core_OG_nucl_MSA.phylip_phyml_tree.txt"
	output:
		"{sample}_ML_tree_with_RR.nwk"
	conda:
		"env/ml_tree_env.yml"
	shell:
		##"phyml_command=$(grep -m1 \"> phyml\" {input.modet_test_log}) ; "
		##"phyml_command=${phyml_command#>} ; "
		##"$phyml_command"
		"phyml  -i {input.msa} -m 012345 -f m -v 0 -a e -c 4 -o tlr; "
		"mv {input.msa}_phyml_tree.txt {output}"

# TODO: tree from OF is different for MSA and not MSA methods!
# TODO: model test on MSA, RAxML_or_PhyML, -> tree for CFML
# + one tree after removing recombination
# originally used panX SNPs tree (inputs/panX...)

# TODO recombination length = msa length!!!
# tried kappa = 0.46, 0.93 from Modeltest_NG substituion rates
# (ag+ct)/(ac+at+cg+gt)

# TODO number of chars is different for the alignments!!!
# 6Mb for CFML, 2,7 Mb for original MSA from OF (phy2)
rule clonalFrameML:
	input:
		#tree_old="inputs/mafft/{sample}/mafft_nucl_test_UPGMA_st5.nwk",
		#msa_old="inputs/mafft/{sample}/mafft_nucl_test.fasta",
		msa="{sample}_concat.xmfa",
		#"/Users/devseeva/Desktop/work/sm_workflow/snakefiles/inputs/panX/phy2_test/mauve/phy2_aln.xmfa",
		#"{sample}_concat.fasta",
		tree="{sample}_ML_tree_with_RR.nwk"
	output:
		dir=directory("outputs/{sample}_cfml")
	conda:
		"env/cfml_env.yml"
	shell:
		"mkdir {output};"
		#"cp DEBUG/rso_test_cfml/* {output}"
		"ClonalFrameML {input.tree} {input.msa} {output} -xmfa_file true ; "
		"mv {output}.* {output}/"

# here I use the absolute path to my R instalation
# may cause errors, should be replaced with R conda env
rule cfmlResultR_plot:
	input:
		dir=rules.clonalFrameML.output.dir
	params:
		files_prefix="{sample}_cfml",
		path2R="/Library/Frameworks/R.framework/Versions/4.0/Resources/Rscript"
	output:
		report("{sample}.cfml.pdf",
		caption="report/cfml.rst", category="CFML: recombination prediction")
	shell:
		#"cp DEBUG/rso_test.cfml.pdf {output}"
		"{params.path2R} scripts/cfml_results.R {input.dir}/{params.files_prefix};"
		"mv {input.dir}/{params.files_prefix}*.pdf {output}"

rule remove_recombinations_from_MSA:
	input:
		msa="{sample}_concat.fasta",
		#"inputs/mafft/{sample}/mafft_nucl_test.fasta",
		rec=rules.clonalFrameML.output.dir
	output:
		"{sample}_concat_no_RR.fasta"
	shell:
		"cp {input.msa} {output}"


# TODO replace mafft with og-core msa!
# TODO optimize branch length:
#	TreeAnc.optimize_branch_lengths_joint: THIS TREE HAS LONG BRANCHES.
#      	****TreeTime's JOINT IS NOT DESIGNED TO OPTIMIZE LONG BRANCHES.
#      	****PLEASE OPTIMIZE BRANCHES USING:
#      	****branch_length_mode='input' or 'marginal'

# treetime --aln core_OG_nucl_MSA.fasta --tree SpeciesTree_rooted_node_labels.txt --dates metainfo.csv --date-column collection_date --plot-tree test_rso_treeTime.pdf --outdir OUTDIR
### !!! TreeTime.reroot: rerooting will ignore covariance and shared ancestry. !!!
rule treeTime:
	input:
		debug=rules.OrthoFinderPhylotype.output.dir,
		tree="{sample}_ML_tree_with_RR.nwk",
		#"/Users/devseeva/Desktop/work/sm_workflow/snakefiles/outputs/rso_test_cfml/rso_test_cfml.labelled_tree.newick",
		msa="{sample}_concat_no_RR.fasta",
		dates="inputs/{sample}_metadata/{sample}_cleaned_data.csv"
	params:
		plot_name="{sample}_treetime.pdf"
	output:
		dir=directory("outputs/{sample}_treetime"),
		plot=report("outputs/{sample}_treetime/{sample}_treetime.pdf",
		caption="report/treeTime.rst", category="Population structure")
	conda:
		"env/biopython_env.yml"
	shell:
		#"mkdir {output.dir}; touch {output.plot}"
		"treetime --aln {input.msa} --tree {input.tree} --dates {input.dates} --date-column collection_date --name-column phy_id --plot-tree {params.plot_name} --outdir {output.dir}"

#TODO Fehler in treeWAS(snps = snps, phen = test_phen_CD, tree = tree, n.subs = n.subs,  :
#  Some elements of names(phen)
#                                                 are absent from tree$tip.label.
#Ausführung angehalten
rule treeWAS:
	input:
		cfml=rules.clonalFrameML.output.dir,
		tttttt="outputs/{sample}_treetime/{sample}_treetime.pdf",
		meta="inputs/{sample}_metadata/{sample}_cleaned_data_TW.csv"
	params:
		prefix="outputs/{sample}_cfml/{sample}_cfml"
	output:
		report("{sample}_treeWAS_plots.pdf",caption="report/treeWAS.rst",category="CFML: recombination prediction")
	script:
		"scripts/treeWAS.R"
	#shell:
	#	"cp DEBUG/rso_test_treeWAS_plots.pdf {output}"

# Part 4: diamond
rule nr_db_preparation:
	input:
		"inputs/nr"
	output:
		"outputs/nr.dmnd"
	shell:
		"touch {output}"

rule diamond_search:
	input:
		queries=rules.divide_proteins_by_phylotypes.output.dir,
		db="outputs/nr.dmnd"
	output:
		"outputs/{sample}_nr_aln_original.daa"
	shell:
		"touch {output}"

rule diamond_meganizer:
	input:
		"outputs/{sample}_nr_aln_original.daa"
	output:
		report("outputs/{sample}_nr_aln_meganized.daa",
		caption="report/megan.rst", category="Diamond: NCBI alignment results")
	shell:
		"touch {output}"

rule diamond_filter:
	input:
		"outputs/{sample}_nr_aln_meganized.daa"
	output:
		temp("{sample}.diamond_interesting_hits.csv")
	shell:
		"touch {output}"

rule assign_OGs_to_hits:
	input:
		ogs="outputs/{sample}_dated_OGs.csv",
		hits="{sample}.diamond_interesting_hits.csv"
	output:
		"{sample}.diamond_interesting_hits_with_OGs.csv"
	shell:
		"touch {output}"

rule get_hits_metadata:
	input:
		"{sample}.diamond_interesting_hits.csv"
	output:
		"outputs/{sample}.diamond_interesting_hits_metadata.csv"
	shell:
		"touch {output}"

rule diamond_plot_interesting_hits:
	input:
		hits="{sample}.diamond_interesting_hits.csv",
		meta="outputs/{sample}.diamond_interesting_hits_metadata.csv"
	output:
		report("{sample}.diamond_interesting_hits.pdf",
		caption="report/diamond.rst", category="Diamond: NCBI alignment results")
	shell:
		"touch {output}"

# Summary
rule HGT_candidates_description:
	input:
		dia_og_time="{sample}.diamond_interesting_hits_with_OGs.csv",
		xeno="outputs/{sample}_xenologs_summary.txt",
		recombination=rules.clonalFrameML.output.dir,
		hits_meta_plots="{sample}.diamond_interesting_hits.pdf"
	output:
		report("{sample}_HGT_candidates_summary.pdf",
		caption="report/hgt_info.rst", category="HGT candidates")
	shell:
		"touch {output}"
