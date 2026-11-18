from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

enzymes = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT"}
dna = input("DNA sequence: ").upper().strip()

if not dna or any(b not in "ATGC" for b in dna):
    print("Invalid DNA sequence")
else:
    print("Length:", len(dna), "GC%:", round(gc_fraction(Seq(dna))*100, 2))
    print("Reverse complement:", Seq(dna).reverse_complement())
    enzyme = input("Enzyme (EcoRI/BamHI/HindIII): ").strip()

    if enzyme not in enzymes:
        print("Invalid enzyme")
    else:
        site = enzymes[enzyme]
        positions = [i for i in range(len(dna)-len(site)+1) if dna[i:i+len(site)] == site]
        print("Site positions:", [p+1 for p in positions])

        if not positions:
            print("Restriction site not found")
        else:
            pos = positions[0]
            for i in range(pos, pos+len(site)):
                for base in "ATGC":
                    if base != dna[i]:
                        new = dna[:i] + base + dna[i+1:]
                        if new[pos:pos+len(site)] != site:
                            print("Mutation:", dna[i], "->", base, "at", i+1)
                            print("Mutated DNA:", new)
                            print("Selected site disrupted")
                            break
                else:
                    continue
                break
