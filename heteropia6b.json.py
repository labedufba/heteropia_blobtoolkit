cd /home/bruno_cajado/FAROL_genome/blobtoolkit/site_blobtoolkit

python3 -c '
import json

# Carrega o JSON existente
with open("assembly_stats.json", "r") as f:
    data = json.load(f)

# Insere o bloco do BUSCO fornecido
busco_stats = {
    "C": 80.1,
    "D": 3.0,
    "F": 4.9,
    "M": 15.0,
    "S": round(80.1 - 3.0, 1),  # Single-copy (C - D) = 77.1
    "n": 672,
    "count": 672,
    "complete": 80.1,
    "duplicated": 3.0,
    "fragmented": 4.9,
    "missing": 15.0,
    "single": round(80.1 - 3.0, 1)
}

data["busco"] = busco_stats
data["busco_stats"] = busco_stats

# Salva o JSON atualizado
with open("assembly_stats.json", "w") as f:
    json.dump(data, f, indent=2)

print("SUCESSO! Dados do BUSCO integrados com sucesso ao assembly_stats.json:")
print(json.dumps(busco_stats, indent=2))
'
