python3 -c '
import json

fasta_path = "/home/beagle/HD1/Heteropia_genome/ASSEMBLY_FAROL/assembly/heteropia_sp.target.fasta"

bin_size = 1000  # Tamanho do bin para subdividir a espiral
scaffold_lengths = []
binned_lengths = []
binned_gc = []
binned_at = []
binned_n = []

print("Lendo e fatiando o genoma para o formato binned do Snail Plot...")
with open(fasta_path, "r") as f:
    current_seq = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                seq = "".join(current_seq).upper()
                l = len(seq)
                scaffold_lengths.append(l)
                # Subdivide o scaffold em bins
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
                current_seq = []
        else:
            current_seq.append(line)
            
    if current_seq:
        seq = "".join(current_seq).upper()
        l = len(seq)
        scaffold_lengths.append(l)
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

scaffold_lengths.sort(reverse=True)
total_len = sum(scaffold_lengths)

def get_nx(target_pct):
    target = total_len * target_pct
    acc = 0
    for l in scaffold_lengths:
        acc += l
        if acc >= target:
            return l
    return 0

# Estrutura completa de 6.000+ linhas esperada pelo index.html
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

with open("assembly_stats.json", "w") as f:
    json.dump(data, f, indent=2)

print(f"SUCESSO! assembly_stats.json gerado com {len(binned_lengths)} bins no formato binned!")
'
