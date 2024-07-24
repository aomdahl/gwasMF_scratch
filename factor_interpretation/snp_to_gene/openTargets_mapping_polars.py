#!/usr/bin/env python3
#TODO- generalize this script, to specify better output files, liftover settings, etc.
import polars as pl
import pandas as pd
import subprocess
import os
import gc
import psutil
import time
# Load query file
query_file = pl.read_csv("/scratch16/abattle4/ashton/snp_networks/custom_l1_factorization/results/panUKBB_complete_41K_OLD/snpids.txt")['SNP'].to_list()
odir = "/scratch16/abattle4/ashton/snp_networks/scratch/factor_interpretation/snp_to_gene"
#parquet_dir = "/data/abattle4/lab_data/openTargets/v2g/"
parquet_dir = "/data/abattle4/lab_data/openTargets/v2g_scores/"

#Some other code I took that helps track usage::
#from https://stackoverflow.com/questions/938733/total-memory-used-by-python-process
def elapsed_since(start):
    return time.strftime("%H:%M:%S", time.gmtime(time.time() - start))


def get_process_memory():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss


def track(func):
    def wrapper(*args, **kwargs):
        mem_before = get_process_memory()
        start = time.time()
        result = func(*args, **kwargs)
        elapsed_time = elapsed_since(start)
        mem_after = get_process_memory()
        print("{}: memory before: {:,}, after: {:,}, consumed: {:,}; exec time: {}".format(
            func.__name__,
            mem_before, mem_after, mem_after - mem_before,
            elapsed_time))
        return result
    return wrapper


def convertRSIDToHG38(query_file, o_dir):
    # Map to HG37
    conv_file = pl.read_csv('/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz', has_header=True, separator="\t",infer_schema_length=10000,schema_overrides={"chrom": pl.String}) \
        .filter(pl.col('rsid').is_in(query_file)) \
        .select(['chrom', 'pos', 'ref', 'alt', 'rsid']) \
        .sort(by=['chrom', 'pos']) \
        .with_columns([
            ( "chr" + pl.col('chrom').cast(str)).alias('chrom'),
            (pl.col('pos') - 1).alias('start'),
            pl.col('pos').alias('end'),
            pl.col('rsid').alias('dat')
        ]) \
    .select(['chrom', 'start', 'end', 'dat']) 
    # Write to a bed file
    conv_file.to_pandas().to_csv('./local_liftover.tmp', sep='\t', index=False, header=False)
    # Liftover
    command = "~/.bin/liftOver local_liftover.tmp /data/abattle4/aomdahl1/reference_data/liftOver_chains/hg19ToHg38.over.chain.gz " + o_dir + "/local_liftover.hg38.tmp failures.tmp"
    subprocess.run(command, shell=True)
    # Check for failures
    command = "grep -v '#' ./failures.tmp | wc -l"
    missed_snps = subprocess.run(command, shell=True, capture_output=True, text=True).stdout.strip()
    print(f"{missed_snps} SNPs could not be lifted over, sorry.")
    # Clean up
    subprocess.run("rm local_liftover.tmp failures.tmp", shell=True)
    # Read the lifted over file
    ret = pl.read_csv('./local_liftover.hg38.tmp', separator='\t', has_header=False)
    #subprocess.run("rm local_liftover.hg38.tmp", shell=True)
    return ret


def iter_slices(df, batch_size):
    def get_batch(df, offset, batch_size):
        batch = df.slice(offset, batch_size)
        batch = batch.collect(streaming=True)
        return batch
    batch = get_batch(df, 0, batch_size)
    # Yield once even if we got passed an empty LazyFrame
    yield batch
    offset = len(batch)
    if offset:
        while True:
            batch = get_batch(df, offset, batch_size)
            len_ = len(batch)
            if len_:
                offset += len_
                yield batch
            else:
                break

#@profile
@track
def parse_vToG_files(parquet_file_path, chunk_size):
    # Create a lazy frame from the Parquet file
    lazy_df = pl.scan_parquet(parquet_file_path)
    ret_tab = None
    counter=0
    for batch in iter_slices(lazy_df, chunk_size):
        #seaerch for SNPsprint(batch)
        counter += 1
        if ret_tab is None:
            """
            ret_tab = batch.with_columns([
                    (pl.col('chr_id').cast(str) + ":" + pl.col('position').cast(str)).alias('hg38'),
                    pl.col('fpred_labels').list.join(",").alias("labels_list"),
                    pl.col('fpred_scores').cast(pl.List(pl.Utf8)).list.join(",").alias("scores_list")
                ]) \
                .filter(pl.col('hg38').is_in(query_ids)).select(['hg38', 'ref_allele', 'alt_allele', 'labels_list','scores_list'])
            """
            ret_tab = batch.with_columns([
                    (pl.col('chr_id').cast(str) + ":" + pl.col('position').cast(str)).alias('hg38') ]) \
                .filter(pl.col('hg38').is_in(query_ids)).select(['hg38', 'ref_allele', 'alt_allele', 'gene_id','overall_score']).unique()
        else:
            """
            new_add= batch.with_columns([
                    (pl.col('chr_id').cast(str) + ":" + pl.col('position').cast(str)).alias('hg38'),
                    pl.col('fpred_labels').list.join(",").alias("labels_list"),
                    pl.col('fpred_scores').cast(pl.List(pl.Utf8)).list.join(",").alias("scores_list")
                ]) \
                .filter(pl.col('hg38').is_in(query_ids)).select(['hg38', 'ref_allele', 'alt_allele', 'labels_list','scores_list'])
            """
            new_add = batch.with_columns([
                    (pl.col('chr_id').cast(str) + ":" + pl.col('position').cast(str)).alias('hg38') ]) \
                .filter(pl.col('hg38').is_in(query_ids)).select(['hg38', 'ref_allele', 'alt_allele', 'gene_id','overall_score']).unique()
            ret_tab = pl.concat([ret_tab, new_add], rechunk=True)
        if counter % 10 == 0:
            print("Currently on batch " + str(counter))
    return ret_tab

#Generators


lifted_over_snps = convertRSIDToHG38(query_file)
query_ids = lifted_over_snps.select(pl.col('column_1').str.replace("chr", "") + ":" + pl.col('column_3').cast(str)).to_series().to_list()

snp_gene_mapping = pl.DataFrame()
opentargets_files = [f for f in os.listdir(parquet_dir) if f.endswith('.parquet')]
chunk_size =100000
for i,parq in enumerate(opentargets_files):
    print(f"Reading in file {parq}")
    test = parse_vToG_files(parquet_dir + parq, chunk_size)
    print(f"Adding {test.shape[0]} new entries")
    snp_gene_mapping = pl.concat([snp_gene_mapping, test], rechunk=True)
    #print(snp_gene_mapping)
    #input()
    print("Search complete in " + str(i+1) + " of " + str(len(opentargets_files)))
    snp_gene_mapping.write_csv(f"{odir}/41K_openTargets.csv", separator='\t')
    gc.collect()
    
print("Finished parsing all!")

#No gc

#parse_vToG_files: memory before: 215,252,992, after: 226,906,112, consumed: 11,653,120; exec time: 00:00:03
# parse_vToG_files: memory before: 226,906,112, after: 227,946,496, consumed: 1,040,384; exec time: 00:00:03
#parse_vToG_files: memory before: 227,946,496, after: 229,253,120, consumed: 1,306,624; exec time: 00:00:03
#parse_vToG_files: memory before: 229,253,120, after: 232,923,136, consumed: 3,670,016; exec time: 00:00:03
#parse_vToG_files: memory before: 232,923,136, after: 224,911,360, consumed: -8,011,776; exec time: 00:00:03
#With GC:
#parse_vToG_files: memory before: 214,302,720, after: 226,471,936, consumed: 12,169,216; exec time: 00:00:03
#parse_vToG_files: memory before: 226,914,304, after: 221,540,352, consumed: -5,373,952; exec time: 00:00:03
#parse_vToG_files: memory before: 221,544,448, after: 224,931,840, consumed: 3,387,392; exec time: 00:00:03
#parse_vToG_files: memory before: 224,931,840, after: 221,802,496, consumed: -3,129,344; exec time: 00:00:03
#parse_vToG_files: memory before: 221,802,496, after: 232,443,904, consumed: 10,641,408; exec time: 00:00:03
#parse_vToG_files: memory before: 232,443,904, after: 224,862,208, consumed: -7,581,696; exec time: 00:00:04
#huh... not really sure what to make of this.