python3 -c '
import json

fasta_path = "/home/beagle/HD1/Heteropia_genome/ASSEMBLY_FAROL/assembly/heteropia_sp.target.fasta"

scaffolds = []
gc_count = 0
at_count = 0
n_count = 0
current_name = ""
current_seq = []

print("Lendo scaffolds do arquivo FASTA...")
with open(fasta_path, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                seq_str = "".join(current_seq).upper()
                l = len(seq_str)
                scaffolds.append({"name": current_name, "length": l})
                gc_count += seq_str.count("G") + seq_str.count("C")
                at_count += seq_str.count("A") + seq_str.count("T")
                n_count += seq_str.count("N")
            current_name = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_seq:
        seq_str = "".join(current_seq).upper()
        l = len(seq_str)
        scaffolds.append({"name": current_name, "length": l})
        gc_count += seq_str.count("G") + seq_str.count("C")
        at_count += seq_str.count("A") + seq_str.count("T")
        n_count += seq_str.count("N")

# Ordena do maior para o menor
scaffolds.sort(key=lambda x: x["length"], reverse=True)
lengths = [s["length"] for s in scaffolds]
total_len = sum(lengths)

def get_nx(target_pct):
    target = total_len * target_pct
    acc = 0
    for l in lengths:
        acc += l
        if acc >= target:
            return l
    return 0

# Estrutura completa esperada pelo SnailPlot em D3
full_data = {
    "assembly": "Heteropia sp.",
    "checkm": {},
    "busco": {},
    "general": {
        "total_length": total_len,
        "scaffold_count": len(scaffolds),
        "n50": get_nx(0.50),
        "n90": get_nx(0.90),
        "gc_content": round((gc_count / total_len) * 100, 2) if total_len else 0,
        "at_content": round((at_count / total_len) * 100, 2) if total_len else 0,
        "n_content": round((n_count / total_len) * 100, 2) if total_len else 0
    },
    "scaffolds": scaffolds,
    "scaffold_lengths": lengths
}

with open("assembly_stats.json", "w") as f:
    json.dump(full_data, f, indent=2)

print(f"SUCESSO! assembly_stats.json gerado com {len(scaffolds)} scaffolds e tamanho total de {total_len} bp.")
'
