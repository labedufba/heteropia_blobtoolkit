python3 -c '
import json

fasta_path = "/home/beagle/HD1/Heteropia_genome/ASSEMBLY_FAROL/assembly/heteropia_sp.target.fasta"

print("Processando scaffolds do FASTA...")
scaffolds = []
gc_total = 0
at_total = 0
n_total = 0

current_seq = []
with open(fasta_path, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                seq = "".join(current_seq).upper()
                l = len(seq)
                scaffolds.append(seq)
                gc_total += seq.count("G") + seq.count("C")
                at_total += seq.count("A") + seq.count("T")
                n_total += seq.count("N")
                current_seq = []
        else:
            current_seq.append(line)
    if current_seq:
        seq = "".join(current_seq).upper()
        l = len(seq)
        scaffolds.append(seq)
        gc_total += seq.count("G") + seq.count("C")
        at_total += seq.count("A") + seq.count("T")
        n_total += seq.count("N")

scaffolds.sort(key=len, reverse=True)
scaffold_lengths = [len(s) for s in scaffolds]
total_len = sum(scaffold_lengths)

# Ajusta o tamanho do bin para gerar cerca de 1000 a 1500 bins no total (arquivo leve de ~80KB)
target_bins = 1000
bin_size = max(1000, total_len // target_bins)

binned_lengths = []
binned_gc = []
binned_at = []
binned_n = []

for seq in scaffolds:
    l = len(seq)
    for i in range(0, l, bin_size):
        sub = seq[i:i+bin_size]
        sub_len = len(sub)
        binned_lengths.append(sub_len)
        gc = sub.count("G") + sub.count("C")
        at = sub.count("A") + sub.count("T")
        n = sub.count("N")
        binned_gc.append(round(gc / sub_len, 4) if sub_len else 0)
        binned_at.append(round(at / sub_len, 4) if sub_len else 0)
        binned_n.append(round(n / sub_len, 4) if sub_len else 0)

def get_nx(target_pct):
    target = total_len * target_pct
    acc = 0
    for l in scaffold_lengths:
        acc += l
        if acc >= target:
            return l
    return 0

data = {
    "assembly": "Heteropia sp.",
    "total_length": total_len,
    "scaffold_count": len(scaffold_lengths),
    "n50": get_nx(0.50),
    "n90": get_nx(0.90),
    "scaffold_lengths": scaffold_lengths,
    "binned_scaffold_lengths": binned_lengths,
    "binned_gc_content": binned_gc,
    "binned_at_content": binned_at,
    "binned_n_content": binned_n
}

# Salva com o nome do exemplo e também como assembly_stats.json
for filename in ["output.assembly-stats_SSA.json", "assembly_stats.json"]:
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)

print(f"SUCESSO! Arquivo leve gerado com tamanho total de {total_len} bp e {len(binned_lengths)} bins.")
'
