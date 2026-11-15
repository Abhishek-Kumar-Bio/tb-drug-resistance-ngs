# Learning log

Answer each question in your own words before making the repository public.

## The biology
- [ ] How does rifampicin kill *M. tuberculosis*, and why do changes in the rpoB RRDR cause resistance?
- [ ] katG activates isoniazid; inhA is its target. Explain why katG loss and *fabG1* promoter mutations both cause
      resistance, by different routes.
- [ ] What makes an isolate MDR? pre-XDR? XDR?
- [ ] Why are gyrA E21Q, S95T and G668D present in every isolate?

## The pipeline
- [ ] What does each tool do: fastp, BWA-MEM, samtools, bcftools? What is in a FASTQ, a BAM and a VCF file?
- [ ] Why `--ploidy 1`? What would the default diploid setting do here?
- [ ] What do QUAL and depth mean for a variant call, and why filter on them?
- [ ] How does the code turn a genome position into "rpoB_S450L"? What changes for a gene on the minus strand?
- [ ] Why does the notebook check codons like rpoB 450 = serine before doing anything else?

## The results
- [ ] Explain sensitivity and specificity, and why 5/6 has a 95 % CI of 44–97 %.
- [ ] For each discordant isolate, give two possible explanations and say what extra data would decide between them.
- [ ] Why was PTAN-0359 only 56 % covered, and why are its RIF/INH calls still trustworthy?
- [ ] Why were the rules not changed after seeing the katG in-frame deletion?
- [ ] The katG frameshift in 0095 is a 1-bp deletion in a run of Cs. Why are such calls error-prone, and why is
      this one probably real?

## Next steps
- [ ] How would you subsample reads at random instead of taking the first N?
- [ ] How would you use the full WHO catalogue? What does TB-Profiler add?
