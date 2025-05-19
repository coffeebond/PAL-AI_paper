import os, pysam, subprocess, pybedtools, gzip, math
from time import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from multiprocessing import Manager, cpu_count
from collections import deque

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
                    
                    for j, alt in enumerate(rec.alts):
                        # check the filter (sometimes the filter has only one entry even if there are multiple alt alleles)
                        if len(rec.filter) == 0:
                            v_filter = 'PASS'
                        else:
                            v_filter = ','.join([x for x in rec.filter])
                                
                        row = [
                            rec.chrom, rec.pos, rec.ref, alt,
                            rec.info.get('AF', ['NA'])[j],
                            rec.info.get('AC', ['NA'])[j],
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
os.makedirs(local_dir, exist_ok=True)
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


NUM_CORES = cpu_count()          # Use all available cores
pwrite(f_log, f"Total number of cores: {NUM_CORES}")

FILES_PER_TASK = 10
MAX_PENDING = NUM_CORES * 2
LOG_INTERVAL = 1000

counter = 0
idx_iter = iter(idx_lst)
futures = []

with gzip.open(output_file, 'wt') as f_out, ProcessPoolExecutor(NUM_CORES) as executor:
    f_out.write('\t'.join([
        'chromosome', 'position', 'ref', 'alt', 'af', 'ac', 'an',
        'id', 'type', 'filter', 'pa_site', 'transcript_id'
    ]) + '\n')

    while True:
        try:
            batch = [next(idx_iter) for _ in range(FILES_PER_TASK)]
        except StopIteration:
            break

        futures.append(executor.submit(process_batch, batch, ref_bed_file, gcs_dir, local_dir))

        if len(futures) >= MAX_PENDING:
            done, _ = wait(futures, return_when=FIRST_COMPLETED)
            for f in done:
                results = f.result()
                for row in results:
                    f_out.write(row + '\n')
                counter += FILES_PER_TASK
                if counter % LOG_INTERVAL == 0:
                    pwrite(f_log, f"Processed: {counter}/{total_files} files | Time elapsed: {timer(t_start)}")
                futures.remove(f)

    # Finish any remaining tasks
    done, _ = wait(futures)
    for f in done:
        results = f.result()
        for row in results:
            f_out.write(row + '\n')
        counter += FILES_PER_TASK
        if counter % LOG_INTERVAL == 0 or counter >= total_files:
            pwrite(f_log, f"Processed: {counter}/{total_files} files | Time elapsed: {timer(t_start)}")

subprocess.run(f"gsutil -u $GOOGLE_PROJECT cp {output_file} {os.getenv('WORKSPACE_BUCKET')}/Output/", shell=True)
pwrite(f_log, f"\nDone! Total time: {timer(t_start)}")
f_log.close()