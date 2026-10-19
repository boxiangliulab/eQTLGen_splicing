# bulk sQTL workflow for eQTLGen: 1–4. DataQC → Imputation → HDF5 → per-cohort encoding

We implement the sQTL workflow as a set of sequential per-cohort stages. This is the cookbook for the **four genotype + encoding stages** that follow splicing-phenotype preparation:

- **Stage 1: `1_DataQC`** — genotype QC **and** splicing sample QC + covariate preparation, in one pipeline.
- **Stage 2: `2_Imputation`** — harmonisation → phasing → imputation of the QC'd genotypes to the 1000G-30x GRCh38 reference.
- **Stage 3: `3_ConvertVcf2Hdf5`** — convert the imputed VCFs into the HASE genotype HDF5 format.
- **Stage 4: `4_PerCohortPreparations`** — HASE-encode genotypes + the combined splice matrix + covariates into the privacy-preserving package you upload to central analysis.

This cookbook **assumes you have already run [Stage 0 (`0_SplicingPhenotype`)](#)** and therefore have the combined splice-site matrix `combined_matrix/{cohort_name}_combined_splice_site.txt.gz` sitting in `0_SplicingPhenotype/output/`, together with the shared `tools/` launcher. These four stages mirror the eQTLGen Phase II genotype pipeline almost exactly; the only sQTL-specific change is that the **phenotype (splicing) sample QC lives in Stage 1 alongside genotype QC**, and the combined splice matrix is carried through to Stage 4 so that genotype, splicing phenotype and covariates are HASE-encoded exactly once.

The **final deliverable** of this whole workflow is produced by Stage 4 — `{cohort_name}_IntermediateFilesEncoded_to_upload/` (+ `.md5`) — and that is the package you share with central analysis. Importantly, **Stage 1 is a review gate**: you share the `1_DataQC` report with the central analysts and wait for approval before running Stage 2.

For this cookbook, we follow the style of the [eQTLGen Phase II cookbook](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html) and aim to require minimal additional user actions.



## Eligibility

If you have genotype data on the same cohort(s) for which you produced a splicing phenotype matrix in Stage 0, you will be able to run this analysis. Genotypes may be either:

- **SNP-array genotypes** (PLINK `.bed/.bim/.fam`) — the default, or
- **whole-genome sequencing (WGS)** genotypes (a bgzipped `.vcf.gz` with unique variant IDs).

By default we assume genotypes are in genome build **hg38 (GRCh38)**. If your genotypes are in a different build, set `genome_build` accordingly in Stage 1 (the pipeline lifts over where needed); the bundled imputation reference is GRCh38.



## Prerequisites

If you have run the workflow in the eQTLGen Phase II cookbook (and Stage 0 here), you should already have all the following. The pipelines and configurations are tailored for the **PBS scheduler** by default; if your HPC uses a different scheduler, modify the job scripts (`#PBS` headers and the `-profile`) to request the appropriate resources.

- High-Performance Computing (HPC) environment
- Singularity / Apptainer
- Nextflow
- The completed **Stage 0** output (`0_SplicingPhenotype/output/combined_matrix/{cohort_name}_combined_splice_site.txt.gz`)
- Your cohort **genotype** data (PLINK bfile or WGS VCF)
- A **genotype-to-expression (GTE)** linking file (two tab-separated columns: genotype ID, RNA/splice sample ID)



## Getting Started

### Setup

Work inside the **same** analysis root folder you created for Stage 0 (e.g. `eQTLGen_splicing_phase1`), so that all stages end up as siblings. Later stages read earlier outputs via `../../` relative paths, and the shared `tools/` launcher plus `singularitycache/` / `nextflowcache/` are shared across stages.

```
cd eQTLGen_splicing_phase1
```

There are **two wrapper flavours** — pick one:

- **Default (lean, recommended if your HPC has internet):** small download. The Singularity **containers auto-pull from quay.io**, and you **download the ~31 GB imputation reference once** (Stage 2).
- **Offline (self-contained):** large download that bundles the `.sif` images **and** the `hg38/` reference — for air-gapped HPCs.

#### Option A — default (lean) wrapper

```
cd eQTLGen_splicing_phase1
wget -O eQTLGen_splicing_1-4_per_cohort_prep_wrapper.tar.gz     "https://www.dropbox.com/scl/fi/n62xk43usxbe5tqxjv306/eQTLGen_splicing_1-4_per_cohort_prep_wrapper.tar.gz?rlkey=uvur02agpdn9h532nfjje353w&dl=1"
wget -O eQTLGen_splicing_1-4_per_cohort_prep_wrapper.tar.gz.md5 "https://www.dropbox.com/scl/fi/3n6io7mfr5h5qe0qm4dtx/eQTLGen_splicing_1-4_per_cohort_prep_wrapper.tar.gz.md5?rlkey=7in4q7zdgikbx4hd6lysdt0cq&dl=1"
md5sum -c eQTLGen_splicing_1-4_per_cohort_prep_wrapper.tar.gz.md5     # -> ...: OK
tar -xzvf eQTLGen_splicing_1-4_per_cohort_prep_wrapper.tar.gz
```

Then **download the imputation reference once** (needed by Stage 2, ~31 GB). From inside `2_Imputation/`, either run the helper:

```
cd 2_Imputation
qsub download_reference.pbs   # wget + md5-check + untar -> hg38/, and copies the HASE ref into Stage 4
cd ..
```

or do it manually (same as the eQTLGen cookbook):

```
cd 2_Imputation
wget -O eQTLGenReferenceFiles.tar.gz     "https://www.dropbox.com/s/6g58ygjg9d2fvbi/eQTLGenReferenceFiles.tar.gz?dl=1"
wget -O eQTLGenReferenceFiles.tar.gz.md5 "https://www.dropbox.com/s/ekfciajzevn6o1l/eQTLGenReferenceFiles.tar.gz.md5?dl=1"
md5sum --check eQTLGenReferenceFiles.tar.gz.md5     # -> eQTLGenReferenceFiles.tar.gz: OK
tar -xzf eQTLGenReferenceFiles.tar.gz               # -> hg38/
cp hg38/hase_reference/1000G-30x.* ../4_PerCohortPreparations/PerCohortDataPreparations/bin/hase/data/
cd ..
```

The **Singularity images are pulled automatically** from `quay.io/eqtlgen/…` on first run (cached under `singularitycache/`) — no `--container` is passed. (Nextflow must have outbound network where it runs the pull; if not, use the offline wrapper.)

#### Option B — offline (self-contained) wrapper

For HPC environments whose compute nodes have **no internet access** (so containers cannot be pulled and the reference cannot be `wget`'d), a self-contained offline bundle is available that ships the `.sif` images, the `hg38/` reference, the Stage-4 HASE reference, PLINK/liftOver, and Nextflow inside one archive. Because it is large (~35 GB) it is **not hosted for direct download — please request it by email from e1101920@u.nus.edu.**

You should now see the four stage folders as siblings of `0_SplicingPhenotype/` and `tools/`. In the **default** wrapper each stage holds only its pipeline work folder (code/config/submit); the `singularity_img/` `.sif` and Stage-2 `hg38/` are **absent** (pulled / downloaded as above). In the **offline** wrapper those are bundled. Layout:

```
.
|-- 0_SplicingPhenotype/                      (from Stage 0 — already present)
|-- 1_DataQC/DataQC/                          DataQC.nf, bin, conf, data, deps (plink + chain), submit_dataqc_template.sh
|       `-- [offline only] ../singularity_img/popassign_v0.6.sif
|                                             (splicing PCA on top-N most-variable junctions + chunked; handles large cohorts via OOM-retry, standard node)
|-- 1_DataQC_GRAM_PCA/DataQC/                 large-cohort Stage-1 ALTERNATIVE (full-matrix Gram/dual PCA; needs a large-memory node)
|-- 2_Imputation/
|   |-- eQTLGenImpute/                         eQTLGenImpute.nf, bin, conf, submit_imputation_template.sh
|   |-- download_reference.pbs                 (default) fetches + unpacks the reference -> hg38/
|   `-- hg38/                                  reference: bundled (offline) OR downloaded (default)
|       |-- genome_reference/  harmonizing_reference/  phasing_reference/  imputation_reference/  hase_reference/
|-- 3_ConvertVcf2Hdf5/ConvertVcf2Hdf5/        ConvertVcf2Hdf5.nf, bin, conf, submit_convert_template.sh
|-- 4_PerCohortPreparations/PerCohortDataPreparations/   PerCohortDataPreparations.nf, bin (incl. hase/data + helperscripts), conf, submit
|                                             (batches the splice phenotype into <=1000-junction files -> HASE streams them; handles large cohorts on a standard node)
`-- tools/nextflow                            (shared launcher — from Stage 0)
```
(In the **offline** wrapper each stage also has `singularity_img/<img>.sif`, Stage 2 has `hg38/`, and Stage 4 ships the HASE reference in `bin/hase/data/`.)

Test the nextflow executable (load Java >=11 first, e.g. `module load [Java >=11 module]`). We encountered issues with Java 21 on some nextflow versions; if it fails, try an older Java (Java 17 worked perfectly).

```
cd tools
./nextflow run info
# If this works you will see:  N E X T F L O W  ~  version 21.10.6
cd ..
```

As in Stage 0, inside each pipeline work folder is a `submit_*_template.sh`. You only ever edit the lines inside the **`USER CONFIGURATION REQUIRED`** block, and within it only the values shown as `[ ... ]` placeholders. Everything else (the container — auto-pulled from quay.io in the default wrapper, or `--container …sif` in the offline one — the `reference_path=../hg38`, resources, `-profile`, and the `../../` input paths to the previous stage) is pre-filled to the bundle's relative layout. **Submit each pipeline from inside its own work folder.**



## Configure and run

Run the stages **in order**. Stages 2–4 each read the previous stage's `../output/`, so do not start one until the previous has finished successfully.

### Stage 1 — `1_DataQC` (genotype QC + splicing QC + expression QC + covariates + reports)

Enter the work folder and edit the template:

```
cd 1_DataQC/DataQC
```


In `submit_dataqc_template.sh`, under `USER CONFIGURATION REQUIRED`, provide:

1. `cohort_name` — your cohort name (used to name outputs and to locate the Stage-0 matrix).
2. `gte_path` — full path to your GTE file (col 1 = genotype ID, col 2 = RNA/splice sample ID).
3. **Genotype input** — PLINK bfile is the default: set `bfile_path` to your genotype prefix (no `.bed/.bim/.fam`). For WGS/VCF instead, leave `bfile_path` empty, uncomment `input_vcf` and point it at your bgzipped `.vcf.gz`. (The command auto-selects whichever is set; `input_vcf` takes priority.)
4. `gen_qc_steps` — `Array` for SNP-array genotypes, or `WGS` for whole-genome-sequencing genotypes.
5. `genome_build` — `hg38` by default.
6. Covariate/QC knobs (sensible defaults pre-filled): `num_gen_pcs=10`, `num_splicing_pcs=50`, `num_expression_pcs=100`, `sd_splice=4` (SD threshold for splicing outlier removal), `int_transform=true` (per-junction rank-INT).
7. `-profile` in the Nextflow command (bottom of the script) and the matching config.

Two Stage-0 inputs are pre-filled to their standard locations by `cohort_name` — leave them unless your data live elsewhere:
- `splice_pheno_dir=../../0_SplicingPhenotype/output` — the combined splice matrix (`combined_matrix/${cohort_name}_combined_splice_site.txt.gz`) and the `.site` files used for the read-depth/missingness summary.
- `expression_matrix=../../0_SplicingPhenotype/output/expression_matrix/${cohort_name}_gene_counts.txt.gz` — the raw gene-count matrix, now a **specified input** (feeds expression QC and the Stage-4 permuted / LD path).

To add cohort-specific covariates (e.g. age, reported sex, batch), use the optional `--AdditionalCovariates` argument (tab-delimited file, **with a header row**; first column = genotype ID; genetically-determined sex is already included automatically).

```
qsub submit_dataqc_template.sh   # or equivalent for your scheduler
```

This single pipeline (`DataQC.nf`) runs:

- **Genotype QC** — per-chromosome normalisation/QC, PLINK conversion/merge, genotype PCA + 1000G population assignment, sex check, and ancestry/heterozygosity/relatedness filtering.
- **Splicing QC** (`ProcessSplicing.R`) — iterative PCA/MDS outlier detection + per-junction rank-INT, restricted to the genotyped samples. **This defines the analysed sQTL sample set.**
- **Expression QC** (`ProcessExpression.R`, a faithful eQTLGen port) — CPM-filter + TMM + per-gene INT, **aligned to the splicing-QC'd samples**. It is *diagnostic*: PCA/MDS outliers are flagged, not removed (unless `--expr_remove_outliers true`). The normalised expression matrix — same samples, same order as splicing — feeds the Stage-4 permuted / LD-estimation path. A hard guard stops the run if any splicing-QC'd sample lacks expression data, so the LD reference always matches the sQTL sample set.

It then writes **two covariate files** and renders **two reports** to `../output/`:

- `CovariatePCs.Splicing.txt` — genotype PCs + splicing PCs + genetic sex (+ any AdditionalCovariates) → empirical / meta-analysis (splicing) path.
- `CovariatePCs.Expression.txt` — genotype PCs + expression PCs + genetic sex (+ any AdditionalCovariates) → permuted / LD (expression) path.
- `Report_DataQc_{cohort}.html` — the QC gate: genotype QC + splicing QC + the expression-based **sex check**.
- `Report_ExpressionDiagnostic_{cohort}.html` — diagnostic only (expression PCA before/after, MDS, scree).

> **Review gate.** When Stage 1 finishes, **share `../output/Report_DataQc_{cohort_name}.html` with the central analysts and wait for their review before running Stage 2.** (`Report_ExpressionDiagnostic_{cohort_name}.html` is for your own QC sanity-checking.)

### Stage 2 — `2_Imputation` (harmonise → phasing → imputation)

Run only after the Stage-1 report is approved.

```
cd ../../2_Imputation/eQTLGenImpute
```

In `submit_imputation_template.sh`, under `USER CONFIGURATION REQUIRED`, provide:

1. `cohort_name` — same cohort name.
2. `qc_input_folder` — pre-filled to `../../1_DataQC/output` (the Stage-1 output). No change needed with the recommended layout.
3. `output_path` — pre-filled to `../output`.
4. `genome_build` — `hg38`.
5. `-profile` in the Nextflow command and the matching config.

The bundled GRCh38 reference (`../hg38`) and container are pre-filled.

```
qsub submit_imputation_template.sh   # or equivalent for your scheduler
```

This runs `eQTLGenImpute.nf`: harmonise the QC'd genotypes to the reference panel, Eagle **phasing**, and Minimac **imputation**, per chromosome. This is the **longest stage** — phasing/imputation of all autosomes can take many hours to a day-plus for large cohorts. The driver is submitted to a long-walltime queue; **if the driver is killed before all chromosomes finish, just resubmit the same script — Nextflow's `-resume` reuses completed chromosomes.** Output: `../output/postimpute/*.filtered.vcf.gz` (one per chromosome), the input to Stage 3.

### Stage 3 — `3_ConvertVcf2Hdf5` (imputed VCF → HASE HDF5)

```
cd ../../3_ConvertVcf2Hdf5/ConvertVcf2Hdf5
```

In `submit_convert_template.sh`, under `USER CONFIGURATION REQUIRED`, provide:

1. `cohort_name` — same cohort name.
2. `genopath` — pre-filled to `../../2_Imputation/output/postimpute`.
3. `outputpath` — pre-filled to `../output`.
4. `-profile` and matching config.

```
qsub submit_convert_template.sh   # or equivalent for your scheduler
```

This runs `ConvertVcf2Hdf5.nf`, converting the imputed VCFs into the HASE genotype HDF5 format (chunked genotype `.h5` + individuals + probes + SNP QC metrics) in `../output/`.


### Stage 4 — `4_PerCohortPreparations` (HASE encoding → upload package)

```
cd ../../4_PerCohortPreparations/PerCohortDataPreparations
```

In `submit_per_cohort_prep_template.sh`, under `USER CONFIGURATION REQUIRED`, the inputs are pre-filled:

1. `genotypes_hdf5` — pre-filled to `../../3_ConvertVcf2Hdf5/output` (Stage-3 HDF5).
2. `qc_data_folder` — pre-filled to `../../1_DataQC/output` (Stage-1 splice matrix, expression matrix, and the two covariate files `CovariatePCs.Splicing.txt` / `CovariatePCs.Expression.txt`).
3. `output_path` — pre-filled to `../output`.
4. `-profile` and matching config.

```
qsub submit_per_cohort_prep_template.sh   # or equivalent for your scheduler
```

This runs `PerCohortDataPreparations.nf`: it HASE-encodes the genotypes, the combined (junction + IR) splice matrix, and the covariates into one privacy-preserving package. Because junction and IR were combined into a single matrix upstream, the genotype is encoded **once**; the two modalities stay separable by their `junction_` / `IR_` feature-ID prefixes. The final `*_IntermediateFilesEncoded_to_upload/` package (+ `.md5`) is written to `../output/`.



## Output

Each stage writes to its own `../output/` folder.

**Stage 1 (`1_DataQC/output/`)** — the DataQC gate deliverables (report shared for review; the rest flows to later stages):
```
output/
├── outputfolder_gen/              genotype QC (gen_data_QCd/*_ToImputation.{bed,bim,fam}, GenotypePCs, pop assignment, plots)
├── outputfolder_splice/           splicing QC (splice_data_QCd/splice_data_preprocessed.txt, splice_PCs/splicing_PCs.txt, plots)
├── outputfolder_exp/              expression QC — DIAGNOSTIC (exp_data_QCd/exp_data_preprocessed.txt = Stage-4 permuted/LD input; exp_PCs, plots)
├── CovariatePCs.Splicing.txt      genotype PCs + splicing PCs + genetic sex (+ AdditionalCovariates)   -> empirical/meta-analysis
├── CovariatePCs.Expression.txt    genotype PCs + expression PCs + genetic sex (+ AdditionalCovariates)  -> permuted/LD
├── Report_DataQc_{cohort}.html              <-- SHARE with central analysts; await review BEFORE Stage 2
└── Report_ExpressionDiagnostic_{cohort}.html   diagnostic only (expression PCA/MDS/scree)
```

**Stage 2 (`2_Imputation/output/`)**:
```
output/
└── postimpute/
    └── chr{1..22}.filtered.vcf.gz   <-- imputed, filtered genotypes (input to Stage 3)
```

**Stage 3 (`3_ConvertVcf2Hdf5/output/`)**:
```
output/
├── genotype/        *_{cohort}.h5   (chunked HASE genotype HDF5)
├── individuals/     {cohort}.h5
├── probes/          {cohort}.h5 (+ hash table)
└── SNPQC/
```

**Stage 4 (`4_PerCohortPreparations/output/`)** — the **final deliverable** to upload to central analysis:
```
output/
├── {cohort}_IntermediateFilesEncoded_to_upload/   <-- UPLOAD THIS (empirical/ + permuted/ + mapper/ + SumStats/)
└── {cohort}_OrganizedFiles.md5
```

Only the Stage-4 `*_IntermediateFilesEncoded_to_upload/` package (with its `.md5`) is uploaded to central analysis. The Stage-1 report is shared for the review gate; all other outputs are intermediates consumed by the next stage.



## Troubleshooting

- Check the scheduler stdout/stderr logs (e.g. `1_dataqc.o*`, `2_imputation.o*`, `3_convert_hdf5.o*`, `4_hase_encode.o*` in PBS) and each work folder's `.nextflow.log` / `logs/*_trace.txt`.
- **Resource retries are normal.** The pipelines retry a failed task up to 3 times, doubling its memory and walltime each attempt; a task that is killed (e.g. exit 137 = out-of-memory, exit 143 = walltime) and then succeeds on a retry is not a failure. Only intervene if a task **exhausts all retries** — then raise that process's `memory`/`time` in the stage's `conf/base.config` (and `max_memory`/`max_time` in `nextflow.config` if they would clamp it) and resubmit with `-resume`.
- **Stage 2 is long** — if the driver hits its walltime, resubmit the same script; `-resume` continues from the last completed chromosome.
- Very large cohorts: **both stages scale automatically** — Stage 1 uses top-feature + chunked splicing PCA (OOM-retry to 192 GB on a standard node; use `1_DataQC_GRAM_PCA` only if you want full-matrix PCs), and Stage 4 batches the splice phenotype so HASE streams it (all steps ≤~24 GB even at 673k junctions). No separate large-cohort folder is needed.
- Please also refer to the "Running, monitoring and debugging Nextflow pipelines" section of the [eQTLGen Phase II cookbook](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html).



## Citation

Nextflow:
[Di Tommaso, P., Chatzou, M., Floden, E. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319 (2017). https://doi.org/10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820)

eQTLGen Phase II cookbook:
[https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html](https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html)

PLINK:
[Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience 2015;4:7.](https://doi.org/10.1186/s13742-015-0047-8)

Eagle (phasing):
[Loh PR, Danecek P, Palamara PF, et al. Reference-based phasing using the Haplotype Reference Consortium panel. Nat Genet 2016;48:1443–1448.](https://doi.org/10.1038/ng.3679)

Minimac4 (imputation):
[Das S, Forer L, Schönherr S, et al. Next-generation genotype imputation service and methods. Nat Genet 2016;48:1284–1287.](https://doi.org/10.1038/ng.3656)

HASE (encoding):
[Roshchupkin GV, Adams HHH, Vernooij MW, et al. HASE: Framework for efficient high-dimensional association analyses. Sci Rep 2016;6:36076.](https://doi.org/10.1038/srep36076)


## Contact

Please contact e1101920@u.nus.edu if you have any other questions.
