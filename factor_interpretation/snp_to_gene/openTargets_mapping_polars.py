#!/usr/bin/env python3
"""
July 23, 2024
Author: Ashton Omdahl, with code scraped from various online posting as noted below
Code review, cleanup, and scriptifying with the aid of ChatGPT, 2024.

To call from the command line
python3 your_script.py --query_file path/to/query_file.txt --output_dir path/to/output_dir --parquet_dir path/to/parquet_dir --chunk_size 100000

"""
import polars as pl
import pandas as pd
import subprocess
import os
import gc
import psutil
import time
import argparse



# Load query file
def load_query_file(query_file_path):
    return pl.read_csv(query_file_path)['SNP'].to_list()

# Memory and execution time tracking utilities
# Taken from  https://stackoverflow.com/questions/938733/total-memory-used-by-python-process
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
    """
    If your input file is in RSIDs, this function will convert them to CHR:POSITION in HG38
    """
    conv_file = pl.read_csv('/data/abattle4/lab_data/GWAS_summary_statistics/PanUKBB/high_quality_common_variants_EUR.txt.bgz', has_header=True, separator="\t", infer_schema_length=10000, schema_overrides={"chrom": pl.String}) \
        .filter(pl.col('rsid').is_in(query_file)) \
        .select(['chrom', 'pos', 'ref', 'alt', 'rsid']) \
        .sort(by=['chrom', 'pos']) \
        .with_columns([
            ("chr" + pl.col('chrom').cast(str)).alias('chrom'),
            (pl.col('pos') - 1).alias('start'),
            pl.col('pos').alias('end'),
            pl.col('rsid').alias('dat')
        ]) \
        .select(['chrom', 'start', 'end', 'dat'])
    # Write to a bed file
    conv_file.to_pandas().to_csv('./local_liftover.tmp', sep='\t', index=False, header=False)
    # Liftover
    command = f"~/.bin/liftOver local_liftover.tmp /data/abattle4/aomdahl1/reference_data/liftOver_chains/hg19ToHg38.over.chain.gz {o_dir}/local_liftover.hg38.tmp failures.tmp"
    subprocess.run(command, shell=True)
    # Check for failures
    command = "grep -v '#' ./failures.tmp | wc -l"
    missed_snps = subprocess.run(command, shell=True, capture_output=True, text=True).stdout.strip()
    print(f"{missed_snps} SNPs could not be lifted over, sorry.")
    # Clean up
    subprocess.run("rm local_liftover.tmp failures.tmp", shell=True)
    # Read the lifted over file
    ret = pl.read_csv('./local_liftover.hg38.tmp', separator='\t', has_header=False)
    subprocess.run("rm local_liftover.hg38.tmp", shell=True)
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

@track
def parse_vToG_files(parquet_file_path, chunk_size, query_ids):
    """
    This lazy frame iterator is a handy tool, based on code from :
    https://github.com/pola-rs/polars/issues/10683
    """
    # Create a lazy frame from the Parquet file
    lazy_df = pl.scan_parquet(parquet_file_path)
    ret_tab = None
    counter = 0
    for batch in iter_slices(lazy_df, chunk_size):
        counter += 1
        if ret_tab is None:
            ret_tab = batch.with_columns([
                (pl.col('chr_id').cast(str) + ":" + pl.col('position').cast(str)).alias('hg38')
            ]) \
            .filter(pl.col('hg38').is_in(query_ids)).select(['hg38', 'ref_allele', 'alt_allele', 'gene_id', 'overall_score']).unique()
        else:
            new_add = batch.with_columns([
                (pl.col('chr_id').cast(str) + ":" + pl.col('position').cast(str)).alias('hg38')
            ]) \
            .filter(pl.col('hg38').is_in(query_ids)).select(['hg38', 'ref_allele', 'alt_allele', 'gene_id', 'overall_score']).unique()
            ret_tab = pl.concat([ret_tab, new_add], rechunk=True)
        if counter % 10 == 0:
            print("Currently on batch " + str(counter))
    return ret_tab

def main(is_rsids, query_file_path, output_dir, parquet_dir, chunk_size):
    query_file = load_query_file(query_file_path)
    if is_rsids:
        snp_coords = convertRSIDToHG38(query_file, output_dir)
    else:
        snp_coords = query_file
    query_ids = snp_coords.select(pl.col('column_1').str.replace("chr", "") + ":" + pl.col('column_3').cast(str)).to_series().to_list()
    snp_gene_mapping = pl.DataFrame()
    opentargets_files = [f for f in os.listdir(parquet_dir) if f.endswith('.parquet')]
    
    for i, parq in enumerate(opentargets_files):
        print(f"Reading in file {parq}")
        snp_dat = parse_vToG_files(os.path.join(parquet_dir, parq), chunk_size, query_ids)
        print(f"Adding {snp_dat.shape[0]} new entries")
        snp_gene_mapping = pl.concat([snp_gene_mapping, snp_dat], rechunk=True)
        print("Search complete in " + str(i + 1) + " of " + str(len(opentargets_files)))
        snp_gene_mapping.write_csv(f"{output_dir}/41K_openTargets.csv", separator='\t')
        gc.collect()

    print("Finished parsing all!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process large Parquet files in chunks.")
    parser.add_argument('--query_file', type=str, required=True, help="Path to the query file. Must contain column with header 'SNP'. This should be either RSIDs or in hg38")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory to store output files.")
    parser.add_argument('--parquet_dir', type=str, default="/data/abattle4/lab_data/openTargets/v2g_scores/", help="Directory containing Parquet files to parse.")
    parser.add_argument('--chunk_size', type=int, default=100000, help="Number of rows per chunk. Default is 100000.")
    parser.add_argument('--rsids', type=bool, default=False, action="store_true", help="Specify this if you are passing in RSIDs. Otherwise, we expect hg38 CHR:POS format")
    args = parser.parse_args()
    
    main(args.rsids, args.query_file, args.output_dir, args.parquet_dir, args.chunk_size)
