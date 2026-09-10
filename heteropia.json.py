python3 -c '
import json

fasta_path = "/home/beagle/HD1/Heteropia_genome/ASSEMBLY_FAROL/assembly/heteropia_sp.target.fasta"

lengths = []
gc_count = 0
at_count = 0
n_count = 0
current_len = 0

print("Lendo o arquivo FASTA da Heteropia sp...")
with open(fasta_path, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_len > 0:
                lengths.append(current_len)
                current_len = 0
        else:
            seq = line.upper()
            current_len += len(seq)
            gc_count += seq.count("G") + seq.count("C")
            at_count += seq.count("A") + seq.count("T")
            n_count += seq.count("N")
    if current_len > 0:
        lengths.append(current_len)

lengths.sort(reverse=True)
total_len = sum(lengths)

def get_nx(target_pct):
    target = total_len * target_pct
    acc = 0
    for l in lengths:
        acc += l
        if acc >= target:
            return l
    return 0

stats = {
    "assembly": "Heteropia sp.",
    "length": total_len,
    "scaffolds": len(lengths),
    "n50": get_nx(0.50),
    "n90": get_nx(0.90),
    "gc": round((gc_count / total_len) * 100, 2) if total_len else 0,
    "at": round((at_count / total_len) * 100, 2) if total_len else 0,
    "n": round((n_count / total_len) * 100, 2) if total_len else 0
}

with open("assembly_stats.json", "w") as f:
    json.dump(stats, f, indent=2)

print("SUCESSO! Dados reais de Heteropia sp. gravados em assembly_stats.json:")
print(json.dumps(stats, indent=2))
'
