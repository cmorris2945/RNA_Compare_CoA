#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# RNA SEQ comparison workflow
# Description: Compares the RNAseq analysis results from DNA Nexus and Cromwell
#
# -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
workflow rnaseq_compare {

    File    workflow_dnanexusbbdukmetrics
    File    workflow_cromwellbbdukmetrics
    String  workflow_sampleId
    File    workflow_dnanexusgenesresults
    File    workflow_cromwellgenesresults
    File    workflow_dnanexusgenomebam
    File    workflow_cromwellgenomebam
    File    workflow_dnanexustranscriptomebam
    File    workflow_cromwelltranscriptomebam
    File    workflow_dnanexusisoformsresults
    File    workflow_cromwellisoformsresults
    File    workflow_dnanexusstarfusion
    File    workflow_cromwellstarfusion
    File    workflow_dnanexusarribafusion
    File    workflow_cromwellarribafusion
    File    workflow_dnanexusmergefusion
    File    workflow_cromwellmergefusion
    File    workflow_dnanexusrnaseqqcmetrics
    File    workflow_cromwellrnaseqqcmetrics

    call trimcompare {
        input:
            sampleId = workflow_sampleId,
            dnanexusmetrics = workflow_dnanexusbbdukmetrics,
            cromwellmetrics = workflow_cromwellbbdukmetrics
    }
    call tpmcompare {
        input:
            sampleId = workflow_sampleId,
            dnanexustpm = workflow_dnanexusgenesresults,
            cromwelltpm = workflow_cromwellgenesresults,
			dnanexusisotpm = workflow_dnanexusisoformsresults,
			cromwellisotpm = workflow_cromwellisoformsresults
    }
    call bamcompare {
    	input:
            sampleId = workflow_sampleId,
            dnanexusgenomebam = workflow_dnanexusgenomebam,
            cromwellgenomebam = workflow_cromwellgenomebam,
            dnanexustranscriptomebam = workflow_dnanexustranscriptomebam,
            cromwelltranscriptomebam = workflow_cromwelltranscriptomebam
    }
    call fusioncompare {
		input:
			sampleId = workflow_sampleId,
			dnanexusstarfusion = workflow_dnanexusstarfusion,
			cromwellstarfusion = workflow_cromwellstarfusion,
			dnanexusarribafusion = workflow_dnanexusarribafusion,
			cromwellarribafusion = workflow_cromwellarribafusion,
			dnanexusmergefusion = workflow_dnanexusmergefusion,
			cromwellmergefusion = workflow_cromwellmergefusion
    }
    call rnaseqqccompare {
		input:
			sampleId = workflow_sampleId,
			dnanexusrnaseqqcmetrics = workflow_dnanexusrnaseqqcmetrics,
			cromwellrnaseqqcmetrics = workflow_cromwellrnaseqqcmetrics
    }
	call eccompare {
        input:
            sampleId = workflow_sampleId,
            dnanexusec = workflow_dnanexusgenesresults,
            cromwellec = workflow_cromwellgenesresults,
			dnanexusisoec = workflow_dnanexusisoformsresults,
			cromwellisoec = workflow_cromwellisoformsresults
    }
}

#-----------------------------------------------------------------------------
# TASK: TRIM METRICS COMPARE
# diff -y --suppress-common-lines dnanexusR1 cromwellR1 > <sampleID>.rnaseq.trimcompare.txt
#   
#-----------------------------------------------------------------------------

task trimcompare {
    String  sampleId
    File    dnanexusmetrics
    File    cromwellmetrics
    String  dockerImage
    String  numCpus
    String  memorySize
    String  diskSize
    String  isPreemptible
    
    command {
    	     	
        (echo "1.RNASEQ TRIMMING COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                                     CROMWELL";\
        echo "--------                                                                     ---------";\
        diff -y --suppress-common-lines <(sed "1d" ${dnanexusmetrics}) <(sed "1d" ${cromwellmetrics})) > ${sampleId}.rnaseq.trimcompare.txt
    }   

	runtime {
        docker: '${dockerImage}'
        cpu: '${numCpus}'
        memory: '${memorySize}'
        disk: '${diskSize}'
        preemptible: '${isPreemptible}'
    }
    
    output {
        File trimcompare = "${sampleId}.rnaseq.trimcompare.txt"
    }
}

#-----------------------------------------------------------------------------
# TASK: TPM COMPARE FOR GENES AND ISOFORMS
# diff -y --suppress-common-lines <(cut -f1,8 ${dnanexustpm}) <(cut -f1,8 ${cromwelltpm})) > ${sampleId}.rnaseq.tpmcompare.txt
#-----------------------------------------------------------------------------
task tpmcompare {
    String  sampleId
    File    dnanexustpm
    File    cromwelltpm
    File    dnanexusisotpm
    File    cromwellisotpm
    String  dockerImage
    String  numCpus
    String  memorySize
    String  diskSize
    String  isPreemptible
    
    	
    command {
        (echo "RNASEQ GENES TPM COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             ---------";\
        echo "gene_symbol  TPM                                                     gene_symbol  TPM";\
        diff -y --suppress-common-lines <(cut -f1,8 ${dnanexustpm}) <(cut -f1,8 ${cromwelltpm});\
		echo "RNASEQ ISOFORMS TPM COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             ---------";\
        echo "trascript_id  TPM                                                     transcript_id  TPM";\
        diff -y --suppress-common-lines <(cut -f1,6 ${dnanexusisotpm}) <(cut -f1,6 ${cromwellisotpm})) > ${sampleId}.rnaseq.tpmcompare.txt
    }
    
	runtime {
        docker: '${dockerImage}'
        cpu: '${numCpus}'
        memory: '${memorySize}'
        disk: '${diskSize}'
        preemptible: '${isPreemptible}'
    }
    
    output {
        File tpmcompare = "${sampleId}.rnaseq.tpmcompare.txt"
    }
}
#-----------------------------------------------------------------------------
# TASK: BAM FILE COMPARE
# diff -y --suppress-common-lines <(samtools flagstat ${dnanexusgenomebam}) <(samtools flagstat ${cromwellgenomebam})) > ${sampleId}.rnaseq.bamcompare.txt
#-----------------------------------------------------------------------------
task bamcompare {
    String  sampleId
    File    dnanexusgenomebam
    File    cromwellgenomebam
    File    dnanexustranscriptomebam
    File    cromwelltranscriptomebam
    String  dockerImage
    String  numCpus
    String  memorySize
    String  diskSize
    String  isPreemptible
    
    
    command{
        (echo "GENOME BAM/CRAM COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             --------";\
        diff -y --suppress-common-lines <(samtools flagstat ${dnanexusgenomebam}) <(samtools flagstat ${cromwellgenomebam});\
        echo;\
        echo "TRANSCRIPTOME BAM COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             --------";\
        diff -y --suppress-common-lines <(samtools flagstat ${dnanexustranscriptomebam}) <(samtools flagstat ${cromwelltranscriptomebam})) > ${sampleId}.rnaseq.bamcompare.txt
    }
    runtime {
        docker: '${dockerImage}'
        cpu: '${numCpus}'
        memory: '${memorySize}'
        disk: '${diskSize}'
        preemptible: '${isPreemptible}'   
    }
    
    output {
        File bamcompare = "${sampleId}.rnaseq.bamcompare.txt"
    }
}

#-----------------------------------------------------------------------------
# TASK: GENE FUSION AND MERGE FUSION COMPARISON
# diff -y --suppress-common-lines <(cut -f1,2,3,6,8 ${dnanexusstarfusion}) <(cut -f1,2,3,6,8 ${cromwellstarfusion}) > ${sampleId}.rnaseq.fusioncompare.txt
#-----------------------------------------------------------------------------
task fusioncompare {
    String  sampleId
    File    dnanexusstarfusion
    File    cromwellstarfusion
    File    dnanexusarribafusion
    File    cromwellarribafusion
    File    dnanexusmergefusion
    File    cromwellmergefusion
    String  dockerImage
    String  numCpus
    String  memorySize
    String  diskSize
    String  isPreemptible
    
    	
    command {
        (echo "RNASEQ STAR FUSION COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             ---------";\
        echo "#FusionName JunctionReadCount SpanningFragCount LeftBreakpoint RightBreakpoint             #FusionName JunctionReadCount SpanningFragCount LeftBreakpoint RightBreakpoint";\
        diff -y --suppress-common-lines <(cut -f1,2,3,6,8 ${dnanexusstarfusion}) <(cut -f1,2,3,6,8 ${cromwellstarfusion});\
		echo;\
		echo "RNASEQ ARRIBA FUSION COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             ---------";\
        echo "#gene1 gene2 breakpoint1 breakpoint2 split_reads1 discordant_mates                          #gene1 gene2 breakpoint1 breakpoint2 split_reads1 discordant_mates";\
        diff -y --suppress-common-lines <(cut -f1,2,5,6,12,14 ${dnanexusarribafusion}) <(cut -f1,2,5,6,12,14 ${cromwellarribafusion});\
		echo;\
		echo "RNASEQ MERGED FUSION COMPARISON REPORT:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             ---------";\
        echo "#FusionName  LeftBreakpoint RightBreakpoint JunctionReadCount	SpanningFragCount           #FusionName LeftBreakpoint RightBreakpoint JunctionReadCount SpanningFragCount";\
        diff -y --suppress-common-lines <(cut -f1,2,3,4,5 ${dnanexusmergefusion}) <(cut -f1,2,3,4,5 ${cromwellmergefusion})) > ${sampleId}.rnaseq.fusioncompare.txt
    }
    
	runtime {
        docker: '${dockerImage}'
        cpu: '${numCpus}'
        memory: '${memorySize}'
        disk: '${diskSize}'
        preemptible: '${isPreemptible}'    
    }
    
    output {
        File fusioncompare = "${sampleId}.rnaseq.fusioncompare.txt"
    } 
}
#-----------------------------------------------------------------------------
# TASK: RNA SEQ QC COMPARISON
# diff -y --suppress-common-lines ${dnanexusrnaseqqcmetrics} ${cromwellrnaseqqcmetrics}) > ${sampleId}.rnaseq.rnaseqqccompare.txt
#-----------------------------------------------------------------------------
task rnaseqqccompare {
    String  sampleId
    File    dnanexusrnaseqqcmetrics
    File    cromwellrnaseqqcmetrics
    String  dockerImage
    String  numCpus
    String  memorySize
    String  diskSize
    String  isPreemptible
    
     		 
    command {
        (echo "RNASEQ QC COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             ---------";\
        diff -y --suppress-common-lines ${dnanexusrnaseqqcmetrics} ${cromwellrnaseqqcmetrics}) > ${sampleId}.rnaseq.rnaseqqccompare.txt
    }

    runtime {
        docker: '${dockerImage}'
        cpu: '${numCpus}'
        memory: '${memorySize}'
        disk: '${diskSize}'
        preemptible: '${isPreemptible}' 
    }
    
    output {
        File rnaseqqccompare = "${sampleId}.rnaseq.rnaseqqccompare.txt"
    } 
}
#-----------------------------------------------------------------------------
# TASK: EXPECTED COUNTS COMPARE FOR GENES AND ISOFORMS
# diff -y --suppress-common-lines <(cut -f1,7 ${dnanexusec}) <(cut -f1,7 ${cromwellec})) > ${sampleId}.rnaseq.expected_countcompare.txt
#-----------------------------------------------------------------------------
task eccompare {
    String  sampleId
    File    dnanexusec
    File    cromwellec
    File    dnanexusisoec
    File    cromwellisoec
    String  dockerImage
    String  numCpus
    String  memorySize
    String  diskSize
    String  isPreemptible
    
    	
    command {
        (echo "RNASEQ GENES EXPECTED COUNTS COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             ---------";\
        echo "gene_symbol  expected_count                                                     gene_symbol  expected_counts";\
        diff -y --suppress-common-lines <(cut -f1,7 ${dnanexusec}) <(cut -f1,7 ${cromwellec});\
		echo "RNASEQ ISOFORMS EXPECTED COUNTS COMPARISON REPORT FOR:${sampleId}";\
        echo "----------------------------------------------------";\
        echo "DNANEXUS                                                             CROMWELL";\
        echo "--------                                                             ---------";\
        echo "trascript_id  expected_count                                                    transcript_id  expected_count";\
        diff -y --suppress-common-lines <(cut -f1,5 ${dnanexusisoec}) <(cut -f1,5 ${cromwellisoec})) > ${sampleId}.rnaseq.expected_countcompare.txt
    }
    
	runtime {
        docker: '${dockerImage}'
        cpu: '${numCpus}'
        memory: '${memorySize}'
        disk: '${diskSize}'
        preemptible: '${isPreemptible}'
    }
    
    output {
        File eccompare = "${sampleId}.rnaseq.expected_countcompare.txt"
    }
}