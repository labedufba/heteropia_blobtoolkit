python3 -c '
import json, os

fasta_path = "/home/beagle/HD1/Heteropia_genome/ASSEMBLY_FAROL/assembly/heteropia_sp.target.fasta"

# 1. Encontra o arquivo JSON de exemplo que funciona
sample_file = "output.assembly-stats_SSA.json"
if not os.path.exists(sample_file):
    # Procura qualquer arquivo json na pasta
    jsons = [f for f in os.listdir(".") if f.endswith(".json")]
    sample_file = jsons[0] if jsons else None

if not sample_file:
    print("ERRO: Nenhum arquivo JSON de exemplo encontrado na pasta.")
    exit(1)

print(f"Usando {sample_file} como modelo de estrutura...")
with open(sample_file, "r") as f:
    template = json.load(f)

# 2. Calcula as estatisticas reais do FASTA da Heteropia sp.
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

# 3. Substitui os valores no molde sem quebrar a estrutura D3
if "assembly" in template: template["assembly"] = "Heteropia sp."
if "total_length" in template: template["total_length"] = total_len
if "scaffold_count" in template: template["scaffold_count"] = len(scaffold_lengths)
if "n50" in template: template["n50"] = get_nx(0.50)
if "n90" in template: template["n90"] = get_nx(0.90)
if "scaffold_lengths" in template: template["scaffold_lengths"] = scaffold_lengths
if "binned_scaffold_lengths" in template: template["binned_scaffold_lengths"] = binned_lengths
if "binned_gc_content" in template: template["binned_gc_content"] = binned_gc
if "binned_at_content" in template: template["binned_at_content"] = binned_at
if "binned_n_content" in template: template["binned_n_content"] = binned_n

# Grava com todos os nomes que o index.html possa estar buscando
out_names = ["assembly_stats.json", "output.assembly-stats_SSA.json", "heteropia_sp.assembly-stats.json"]
for name in out_names:
    with open(name, "w") as f:
        json.dump(template, f, indent=2)

print(f"SUCESSO! Moldes atualizados com {total_len} bp e {len(scaffolds)} scaffolds para Heteropia sp.")
'
