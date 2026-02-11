import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from tempfile import TemporaryDirectory

def count_sequences(fasta_path):
    """Count number of sequences in FASTA (fast grep)."""
    cmd = ["grep", "-c", "^>", fasta_path]
    try:
        return int(subprocess.check_output(cmd).decode().strip())
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to count sequences in {fasta_path}: {e}")

def split_fasta(input_fasta, num_chunks=20, output_dir='nt_chunks'):
    """Split FASTA into N chunks by sequence count using awk (fast for huge files)."""
    os.makedirs(output_dir, exist_ok=True)

    total_seqs = count_sequences(input_fasta)
    print(f"Total sequences: {total_seqs}")

    seqs_per_chunk = total_seqs // num_chunks
    remainder = total_seqs % num_chunks

    chunk_files = []
    current_seq = 1
    for i in range(1, num_chunks + 1):
        chunk_size = seqs_per_chunk + (1 if i <= remainder else 0)
        chunk_file = os.path.join(output_dir, f"nt_chunk_{i:03d}.fasta")

        # Use awk to extract exact range of sequences (lines: 2*(seq-1) +1 to 2*seq_size)
        start_line = (current_seq - 1) * 2 + 1
        end_line = (current_seq + chunk_size - 1) * 2

        awk_cmd = f"awk 'NR >= {start_line} && NR <= {end_line}' {input_fasta} > {chunk_file}"
        subprocess.run(awk_cmd, shell=True, check=True)

        # Verify no data loss (count sequences in chunk)
        chunk_seqs = count_sequences(chunk_file)
        if chunk_seqs != chunk_size:
            raise RuntimeError(f"Chunk {chunk_file} has {chunk_seqs} seqs, expected {chunk_size} — data loss!")

        chunk_files.append(chunk_file)
        current_seq += chunk_size

    print(f"Split into {len(chunk_files)} chunks in {output_dir}")
    return chunk_files

def build_mmi(chunk_file):
    """Build .mmi index for a FASTA chunk (minimap2 -d)."""
    mmi_file = chunk_file + ".mmi"
    cmd = ["minimap2", "-d", mmi_file, chunk_file]
    try:
        subprocess.run(cmd, check=True)
        print(f"Built {mmi_file}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to build mmi for {chunk_file}: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 split_build_nt.py /path/to/nt.fa num_chunks [output_dir=nt_chunks] [max_workers=20]")
        sys.exit(1)

    nt_fasta = sys.argv[1]
    num_chunks = int(sys.argv[2])
    output_dir = sys.argv[3] if len(sys.argv) > 3 else 'nt_chunks'
    max_workers = int(sys.argv[4]) if len(sys.argv) > 4 else 20  # Parallel build threads — adjust for your machine

    # Split in temp dir to avoid partial writes on failure
    with TemporaryDirectory() as tmpdir:
        chunks = split_fasta(nt_fasta, num_chunks, tmpdir)

        # Build in parallel
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(build_mmi, chunk) for chunk in chunks]
            for future in as_completed(futures):
                future.result()  # Raise on error — ensures correctness

        # Move completed chunks + .mmi to final output_dir
        os.makedirs(output_dir, exist_ok=True)
        for chunk in chunks:
            os.rename(chunk, os.path.join(output_dir, os.path.basename(chunk)))
            os.rename(chunk + ".mmi", os.path.join(output_dir, os.path.basename(chunk) + ".mmi"))

    print(f"All done! Chunked FASTA + .mmi files in {output_dir}")