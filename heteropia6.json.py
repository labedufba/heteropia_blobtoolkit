python3 -c '
import json

fasta_path = "/home/beagle/HD1/Heteropia_genome/ASSEMBLY_FAROL/assembly/heteropia_sp.target.fasta"

print("Lendo e calculando métricas reais de Heteropia sp...")
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
                scaffolds.append(seq)
                gc_total += seq.count("G") + seq.count("C")
                at_total += seq.count("A") + seq.count("T")
                n_total += seq.count("N")
                current_seq = []
        else:
            current_seq.append(line)
    if current_seq:
        seq = "".join(current_seq).upper()
        scaffolds.append(seq)
        gc_total += seq.count("G") + seq.count("C")
        at_total += seq.count("A") + seq.count("T")
        n_total += seq.count("N")

scaffolds.sort(key=len, reverse=True)
scaffold_lengths = [len(s) for s in scaffolds]
total_len = sum(scaffold_lengths)

bin_size = max(1000, total_len // 1000)

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

# Função para calcular NX (Length) e o número acumulado de sequências (Count)
def get_nx_stats(target_pct):
    target = total_len * target_pct
    acc = 0
    for idx, l in enumerate(scaffold_lengths, 1):
        acc += l
        if acc >= target:
            return l, idx
    return 0, 0

scaffold_n50_len, scaffold_n50_count = get_nx_stats(0.50)
scaffold_n90_len, scaffold_n90_count = get_nx_stats(0.90)

gc_pct = round((gc_total / total_len) * 100, 2) if total_len else 0
at_pct = round((at_total / total_len) * 100, 2) if total_len else 0
n_pct = round((n_total / total_len) * 100, 2) if total_len else 0

# Estrutura completa com todas as equivalências de chaves que o JS busca na tabela
data = {
    "assembly": "Heteropia sp.",
    "span": total_len,
    "total_length": total_len,
    
    # Porcentagens de base (elimina o NaN)
    "gc": gc_pct,
    "at": at_pct,
    "n": n_pct,
    "gc_percent": gc_pct,
    "at_percent": at_pct,
    "n_percent": n_pct,
    "gc_content": gc_pct,
    "at_content": at_pct,
    "n_content": n_pct,
    
    # Scaffolds
    "scaffold_count": len(scaffold_lengths),
    "scaffolds": len(scaffold_lengths),
    "longest_scaffold": scaffold_lengths[0],
    "scaffold_n50": scaffold_n50_len,
    "scaffold_n50_length": scaffold_n50_len,
    "scaffold_n50_count": scaffold_n50_count,
    "scaffold_n90": scaffold_n90_len,
    "scaffold_n90_length": scaffold_n90_len,
    "scaffold_n90_count": scaffold_n90_count,
    "n50": scaffold_n50_len,
    "n90": scaffold_n90_len,
    
    # Contigs (como são idênticos aos scaffolds no seu FASTA)
    "contig_count": len(scaffold_lengths),
    "contigs": len(scaffold_lengths),
    "longest_contig": scaffold_lengths[0],
    "contig_n50": scaffold_n50_len,
    "contig_n50_length": scaffold_n50_len,
    "contig_n50_count": scaffold_n50_count,
    "contig_n90": scaffold_n90_len,
    "contig_n90_length": scaffold_n90_len,
    "contig_n90_count": scaffold_n90_count,
    
    # Arrays de distribuição e binning
    "scaffold_lengths": scaffold_lengths,
    "contig_lengths": scaffold_lengths,
    "binned_scaffold_lengths": binned_lengths,
    "binned_gc_content": binned_gc,
    "binned_at_content": binned_at,
    "binned_n_content": binned_n
}

# Salva nos nomes padronizados que o site usa
with open("assembly_stats.json", "w") as f:
    json.dump(data, f, indent=2)

print("SUCESSO! assembly_stats.json regerado com todas as métricas de Scaffolds e Contigs salvas.")
'
