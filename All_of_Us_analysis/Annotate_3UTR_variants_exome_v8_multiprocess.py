import os, pysam, subprocess, pybedtools, gzip, math
from time import time
from datetime import datetime
from multiprocessing import Pool, cpu_count
def timer(t_start):
    return f"{int(time() - t_start)}s"

def pwrite(f, text, timestamp=None):
    """Write printed output to a file and print it to console."""
    log_text = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {text}" if timestamp else text
    if f:
        f.write(log_text + '\n')
        f.flush()
    print(log_text)

def process_batch(batch_indices, ref_bed_file, gcs_dir, local_dir):
    # Processes a batch of VCF files (called by each worker)
    batch_results = []
    
    for idx in batch_indices:
        # Download interval list
        interval_list_file = f"{gcs_dir}{idx}.interval_list"
        local_interval_list = os.path.join(local_dir, f"{idx}.interval_list")
        subprocess.run(f"gsutil -u $GOOGLE_PROJECT cp {interval_list_file} {local_interval_list}", shell=True, check=True)

        # Convert interval_list to a BED file
        interval_bed_file = os.path.join(local_dir, f"{idx}.bed")
        with open(local_interval_list, 'r') as infile, open(interval_bed_file, 'w') as outfile:
            for line in infile:
                if not line.startswith("@"):
                    chrom, start, end = line.strip().split("\t")[:3]
                    outfile.write(f"{chrom}\t{int(start)-1}\t{end}\n")

        # Intersect with reference BED
        ref_bed = pybedtools.BedTool(ref_bed_file)
        intervals = pybedtools.BedTool(interval_bed_file)
        filtered_bed = ref_bed.intersect(intervals, u=True)

        if len(filtered_bed) > 0:
            # Download VCF and index
            vcf_file = f"{gcs_dir}{idx}.vcf.bgz"
            local_vcf = os.path.join(local_dir, f"{idx}.vcf.bgz")
            subprocess.run(f"gsutil -u $GOOGLE_PROJECT cp {vcf_file} {local_vcf}", shell=True, check=True)
            subprocess.run(f"gsutil -u $GOOGLE_PROJECT cp {vcf_file}.tbi {local_vcf}.tbi", shell=True, check=True)

            # Process variants
            vcf_in = pysam.VariantFile(local_vcf)
            for bed_entry in filtered_bed:
                for rec in vcf_in.fetch(bed_entry.chrom, int(bed_entry.start), int(bed_entry.stop)):
                    if 'AN' in rec.info.keys():
                        v_an = rec.info['AN']
                    else:
                        v_an = 'NA'
                    
                    if rec.id is None:
                        v_id = 'NA'
                    else:
                        v_id = rec.id
                    
                    if 'variant_type' in rec.info.keys():
                        v_type = rec.info['variant_type']
                    else:
                        v_type = 'NA'
                    
                    # each allele may have multiple alternative alleles
                    for j in range(len(rec.alts)):
                        # check the filter (sometimes the filter has only one entry even if there are multiple alt alleles)
                        if len(rec.filter) == 0:
                            v_filter = 'PASS'
                        else:
                            v_filter = ','.join([x for x in rec.filter])

                        v_alt = rec.alts[j]
                        if 'AF' in rec.info.keys():
                            v_af = rec.info['AF'][j]
                        else:
                            v_af = 'NA'
                        if 'AC' in rec.info.keys():
                            v_ac = rec.info['AC'][j]
                        else:
                            v_ac = 'NA'
                                
                        row = [
                            rec.chrom, rec.pos, rec.ref, 
                            v_alt,
                            v_af,
                            v_ac,
                            v_an,
                            v_id,
                            v_type,
                            v_filter,
                            bed_entry.name.split('__')[0],  # pa_site
                            bed_entry.name.split('__')[1]   # transcript_id
                        ]
                        batch_results.append("\t".join(map(str, row)))

            # Cleanup
            os.remove(local_vcf)
            os.remove(f"{local_vcf}.tbi")

        # Cleanup
        os.remove(local_interval_list)
        os.remove(interval_bed_file)
    
    return batch_results


t_start = time()
local_dir = "Temp_files/"
gcs_dir = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/vcf/"
f_log = open("Annotate_3UTR_variants_exome_v8_log.txt", 'w')

ref_bed_file = 'gencode.v25.annotation_annotation_fixed_UTR3_polished_by_PA_site_PolyA_DB.bed'
if not os.path.exists(ref_bed_file):
    subprocess.run(f"gsutil cp {os.getenv('WORKSPACE_BUCKET')}/{ref_bed_file} .")
    
output_file = "All_of_us_exome_v8_variant_call_stats.txt.gz"

# List all VCF files
list_cmd = f"gsutil -u $GOOGLE_PROJECT ls {gcs_dir}*.vcf.bgz"
vcf_files = subprocess.check_output(list_cmd, shell=True).decode('utf-8').splitlines()
idx_lst = [os.path.basename(f).split('.')[0] for f in vcf_files]
total_files = len(idx_lst)
pwrite(f_log, f"Total number of files to process: {total_files}")


# Configurable settings
NUM_CORES = cpu_count()          # Use all available cores
FILES_PER_CORE_PER_BATCH = 10    # Process this many files per core at a time
BATCH_SIZE = NUM_CORES * FILES_PER_CORE_PER_BATCH  # Total files per batch
pwrite(f_log, f"Total number of cores: {NUM_CORES}")

# Split into batches
num_batches = math.ceil(total_files / BATCH_SIZE)
processed_files = 0

# Write output header
with gzip.open(output_file, 'wt') as f:
    f.write('\t'.join([
        'chromosome', 'position', 'ref', 'alt', 'af', 'ac', 'an', 
        'id', 'type', 'filter', 'pa_site', 'transcript_id'
    ]) + '\n')

# Process batches sequentially (with parallel processing within each batch)
for batch_num in range(num_batches):
    batch_start = batch_num * BATCH_SIZE
    batch_end = min((batch_num + 1) * BATCH_SIZE, total_files)
    batch_indices = idx_lst[batch_start:batch_end]

    pwrite(f_log, f"\nProcessing batch {batch_num + 1}/{num_batches} (files {batch_start + 1}-{batch_end})...")

    # Split batch into sub-batches for each core
    sub_batches = [batch_indices[i::NUM_CORES] for i in range(NUM_CORES)]

    with Pool(NUM_CORES) as pool:
        batch_results = pool.starmap(
            process_batch,
            [(sub_batch, ref_bed_file, gcs_dir, local_dir) for sub_batch in sub_batches]
        )

    # Write results for this batch
    with gzip.open(output_file, 'at') as f:  # 'at' = append in text mode
        for sub_batch_result in batch_results:
            for row in sub_batch_result:
                f.write(row + '\n')

    processed_files += len(batch_indices)
    pwrite(f_log, f"Completed: {processed_files}/{total_files} files | Time elapsed: {timer(t_start)}")

# Upload final output to GCS
subprocess.run(f"gsutil -u $GOOGLE_PROJECT cp {output_file} {os.getenv('WORKSPACE_BUCKET')}/Output/", shell=True)
pwrite(f_log, f"\nDone! Total time: {timer(t_start)}")

