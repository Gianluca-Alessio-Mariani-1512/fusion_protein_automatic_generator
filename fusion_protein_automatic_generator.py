import requests
from Bio.Seq import Seq
import argparse

def get_input_from_vcf(data_file):
    found_row = None
    genome_reference_string = None

    with open(data_file, 'r') as file:
        for line in file:
            if line.startswith('##FusionSampleOverallCall=POSITIVE'):
                found_row = line.strip()  # Rimuove eventuali spazi bianchi
            if line.startswith('##reference='):
                genome_reference_string = line.strip()  # Rimuove eventuali spazi bianchi
            
            # Esci dal ciclo se entrambe le righe sono state trovate
            if found_row and genome_reference_string:
                break

    start_index = genome_reference_string.find("reference=") + len("reference=")
    if start_index != -1:
        genome_reference = genome_reference_string[start_index:].strip()  # Rimuove eventuali spazi bianchi
        print("Reference Genome:", genome_reference)
    else:
        print("No Reference Genome Found")

    # Trova l'indice di inizio e fine della sottostringa
    start_index = found_row.find("IsoformsDetected=") + len("IsoformsDetected=")
    end_index = found_row.find("]", start_index)
    # Estrai la sottostringa
    if start_index != -1 and end_index != -1:
        variant_found = found_row[start_index:end_index].strip()  # Rimuove eventuali spazi bianchi
        print("Found Variant:", variant_found)
    else:
        print("No Variant Found")

    found_lines = []

    # Apri il file e cerca la stringa
    with open(data_file, 'r') as file:
        for line in file:
            if variant_found in line:
                found_lines.append(line.strip())  # Aggiungi la riga alla lista, rimuovendo spazi bianchi

    # Ignora la prima riga trovata, se esiste
    if found_lines:
        found_lines = found_lines[1:]

    couple_list = []

    for i in found_lines:
        # Suddividi la stringa in base ai caratteri di tabulazione
        elements = i.split('\t')
        # Ottieni i primi due elementi
        chrnum_str = elements[0]
        chrnum = ''.join(char for char in chrnum_str if char.isdigit())
        break_point = elements[1]

        eighth_element = elements[7].split(';')

        if eighth_element[3].startswith('GENE_NAME='):
        # Estrai il valore dopo "="
            gene_name = eighth_element[3].split('=')[1]

        # Stampa i risultati
        couple_list.append((gene_name, chrnum, break_point))

    return genome_reference, couple_list

def get_exons_from_ucsc(genename, chrnum, start, end, assembly):
    url = f"http://api.genome.ucsc.edu/getData/track?track=knownGene;genome={assembly};chrom=chr{chrnum};start={start};end={end}"
    response = requests.get(url)
    data = response.json()
    genes = data["knownGene"]
    for i in genes:
        tags = i["tag"].split(',')
        if 'MANE_Select' in tags and i["geneName"] == genename:
            chosen_gene = i
            break

    return chosen_gene

def check_for_MANE_zhan(assembly, genename, gene_id):
    url = f"https://github.com/zhanyinx/variantalker/tree/main/resources/mane.transcript.gencode.v43.{assembly}.txt"
    response = requests.get(url)
    MANE_text = response.text
    if gene_id in MANE_text:
        print(f"The gene {genename} identified with the UCSC Genome code {gene_id} is present in the MANE list in the file {url}")
    else:
        print(f"The gene {genename} identified with the UCSC Genome code {gene_id} is NOT present in the MANE list in the file {url}")

def get_seq2check(gene_UCSC_entry):
    accession_number = gene_UCSC_entry['geneName2']
    # Costruisci l'URL per ottenere le informazioni sulla proteina
    url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"
    
    # Effettua la richiesta GET
    response = requests.get(url)
    
    if response.ok:
        # Estrai la sequenza dal contenuto della risposta
        fasta_content = response.text
        # La sequenza è contenuta nelle righe successive all'intestazione
        sequence_lines = fasta_content.splitlines()[1:]  # Salta la prima riga (intestazione)
        protein_sequence = ''.join(sequence_lines)  # Unisci le righe della sequenza
        
        seq2check = protein_sequence[-10:]
        return seq2check
    else:
        print("Errore nella richiesta:", response.text)
        return None

def get_exon_sequence(gene_UCSC_entry, break_end, term_type, strand):
    # Estrai e converti i dati
    chrnum=gene_UCSC_entry["chrom"]
    block_sizes = list(map(int, gene_UCSC_entry['blockSizes'].strip(',').split(',')))
    chrom_starts = list(map(int, gene_UCSC_entry['chromStarts'].strip(',').split(',')))

    # Calcola le posizioni di inizio e fine degli esoni
    exons = []
    for start, size in zip(chrom_starts, block_sizes):
        exon_start = start + gene_UCSC_entry['chromStart']
        exon_end = exon_start + size
        if term_type == "N" and strand == "+":
            if exon_end <= break_end:  # Solo esoni che terminano prima del breakpoint
                exons.append((exon_start, exon_end))
        if term_type == "C" and strand == "+":
            if exon_start >= break_end:  # Solo esoni che terminano prima del breakpoint
                exons.append((exon_start, exon_end))
        if term_type == "N" and strand == "-":
            if exon_start >= break_end:  # Solo esoni che terminano prima del breakpoint
                exons.append((exon_start, exon_end))
        if term_type == "C" and strand == "-":
            if exon_end <= break_end:  # Solo esoni che terminano prima del breakpoint
                exons.append((exon_start, exon_end))

    codingseq_fragments = []

    if term_type == "N" and strand == "+":
        first_exon = exons.pop(0)
        first_exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=first_exon[0], end=first_exon[1])
        trimmed_start_exon = trim_start_cod(dna_seq=first_exon_seq)
        codingseq_fragments.append(trimmed_start_exon)
        for i in range(len(exons)):
            exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=exons[i][0], end=exons[i][1])
            codingseq_fragments.append(exon_seq)
        coding_seq = "".join(codingseq_fragments)
        return coding_seq
    
    if term_type == "C" and strand == "+":
        last_exon = exons.pop(-1)
        for i in range(len(exons)):
            exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=exons[i][0], end=exons[i][1])
            codingseq_fragments.append(exon_seq)
        last_exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=last_exon[0], end=last_exon[1])
        trimmed_end_exon = trim_stop_cod(dna_seq=last_exon_seq, seq2check=get_seq2check(gene_UCSC_entry=gene_UCSC_entry))
        codingseq_fragments.append(trimmed_end_exon)
        coding_seq = "".join(codingseq_fragments)
        return coding_seq
    
    if term_type == "N" and strand == "-":
        exons = exons[::-1]
        first_exon = exons.pop(0)
        first_exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=first_exon[0], end=first_exon[1])
        dna_seq = Seq(first_exon_seq)
        complement_seq = dna_seq.complement()
        reverse_complement_seq = complement_seq[::-1]
        trimmed_start_exon = trim_start_cod(dna_seq=reverse_complement_seq)
        codingseq_fragments.append(trimmed_start_exon)
        for i in range(len(exons)):
            exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=exons[i][0], end=exons[i][1])
            dna_seq = Seq(exon_seq)
            complement_seq = dna_seq.complement()
            reverse_complement_seq = complement_seq[::-1]
            codingseq_fragments.append(reverse_complement_seq)
        coding_seq = "".join(codingseq_fragments)
        return coding_seq
    
    if term_type == "C" and strand == "-":
        exons = exons[::-1]
        last_exon = exons.pop(-1)
        for i in range(len(exons)):
            exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=exons[i][0], end=exons[i][1])
            dna_seq = Seq(exon_seq)
            complement_seq = dna_seq.complement()
            reverse_complement_seq = complement_seq[::-1]
            codingseq_fragments.append(str(reverse_complement_seq))
        last_exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=last_exon[0], end=last_exon[1])
        dna_seq = Seq(last_exon_seq)
        complement_seq = dna_seq.complement()
        reverse_complement_seq = complement_seq[::-1]
        trimmed_end_exon = trim_stop_cod(dna_seq=reverse_complement_seq, seq2check=get_seq2check(gene_UCSC_entry=gene_UCSC_entry))
        codingseq_fragments.append(str(trimmed_end_exon))
        coding_seq = "".join(codingseq_fragments)
        return coding_seq

def get_sequence(chrnum, assembly, start, end):
    url = f"http://api.genome.ucsc.edu/getData/sequence?track=knownGene;genome={assembly};chrom={chrnum};start={start};end={end}"
    response = requests.get(url)
    data_ex = response.json()

    dna_seq = (data_ex["dna"])
    return dna_seq

def trim_start_cod(dna_seq):
    start_codon = "ATG"

    # Trovare l'indice del codone di inizio
    start_codon_index = dna_seq.find(start_codon)

    # Rimuovere la sequenza di basi prima dell'AUG
    if start_codon_index != -1:
        trimmed_start_exon_seq = dna_seq[start_codon_index:]
    else:
        trimmed_start_exon_seq = None  # Se il codone non viene trovato
    return trimmed_start_exon_seq

def trim_stop_cod(dna_seq, seq2check):
    last_exon_frames = [dna_seq, f"A{dna_seq}", f"AA{dna_seq}"]
    last_exon_frames_prot = []

    couples = []

    for index, i in enumerate(last_exon_frames):
        
        final_dnaseq = Seq(i)

        protein = final_dnaseq.translate()

        # Trova l'indice del primo asterisco
        stop_cod = protein.find('*')

        # Se l'asterisco è presente, fai slicing
        if stop_cod != -1:
            trimmed_end_exon = protein[:stop_cod]  # Parte prima dell'asterisco
        else:
            trimmed_end_exon = protein  # Se non c'è asterisco, prendi la stringa originale
        last_exon_frames_prot.append(trimmed_end_exon)
        couples.append((index, trimmed_end_exon))

    valid_exon_frame = [s for s in last_exon_frames_prot if s.endswith(seq2check)]

    first_element = None
    for first, second in couples:
        if second == str(valid_exon_frame[0]):
            first_element = first
            break

    coding_len = ((len(str(valid_exon_frame[0]))*3) - first_element)

    trimmed_end_exon_seq = dna_seq[:coding_len]
    return trimmed_end_exon_seq

def check_correct_frame(N_term_codingseq, C_term_codingseq, seq2check):
    codingseq_frames = [C_term_codingseq, f"A{C_term_codingseq}", f"AA{C_term_codingseq}"]
    codingseq_frames_prot = []
    couples = []

    for index, i in enumerate(codingseq_frames):
        final_dnaseq = Seq(f"{N_term_codingseq}{i}")
        protein = final_dnaseq.translate()
        codingseq_frames_prot.append(protein)
        couples.append((index, protein))

    valid_exon_frame = [s for s in codingseq_frames_prot if s.endswith(seq2check)]

    first_element = None
    for first, second in couples:
        if second == str(valid_exon_frame[0]):
            first_element = first
            break

    return codingseq_frames[first_element], first_element

def main():
    parser = argparse.ArgumentParser(description="Path of the vcf input file")
    parser.add_argument('-p', '--path', type=str, required=False, 
                        help="Path of the vcf input file")
    
    args = parser.parse_args()
    
    # data_file = "data/22-C-004462_v1_22-C-004462_RNA_v1_Non-Filtered_2025-02-20_03_33_25.vcf"
    data_file = args.path

    vcf_output = get_input_from_vcf(data_file=data_file)
    assembly = str(vcf_output[0])

    genename1 = vcf_output[1][0][0]
    chrnum1=int(vcf_output[1][0][1])
    break_end1 = int(vcf_output[1][0][2])

    genename2 = vcf_output[1][1][0]
    chrnum2=int(vcf_output[1][1][1])
    break_end2 = int(vcf_output[1][1][2])

    genes1 = get_exons_from_ucsc(chrnum=chrnum1, start=break_end1 - 1, end=break_end1 + 1, assembly=assembly, genename=genename1)["name"]
    genes1_entry = get_exons_from_ucsc(chrnum=chrnum1, start=break_end1 - 1, end=break_end1 + 1, assembly=assembly, genename=genename1)

    genes2 = get_exons_from_ucsc(chrnum=chrnum2, start=break_end2 - 1, end=break_end2 + 1, assembly=assembly, genename=genename2)["name"]
    genes2_entry = get_exons_from_ucsc(chrnum=chrnum2, start=break_end2 - 1, end=break_end2 + 1, assembly=assembly, genename=genename2)

    gene_id_base1 = genes1.split('.')[0]
    gene_id_base2 = genes2.split('.')[0]

    check_for_MANE_zhan(assembly=assembly, gene_id=gene_id_base1, genename=genename1)
    check_for_MANE_zhan(assembly=assembly, gene_id=gene_id_base2, genename=genename2)

    codingseq1 = get_exon_sequence(gene_UCSC_entry=genes1_entry, break_end=break_end1, term_type="N", strand=genes1_entry['strand'])
    codingseq2 = get_exon_sequence(gene_UCSC_entry=genes2_entry, break_end=break_end2, term_type="C", strand=genes2_entry['strand'])

    codingseq2, first_element = check_correct_frame(N_term_codingseq=codingseq1, C_term_codingseq=codingseq2, seq2check=get_seq2check(gene_UCSC_entry=genes2_entry))

    rem1 = len(codingseq1) % 3
    rem2 = len(codingseq2) % 3

    codingseq1_adj = codingseq1[:rem1]
    final_dnaseq1 = Seq(codingseq1)
    codingseq2_adj = f"{codingseq1[-rem1:]}{codingseq2}"
    final_dnaseq2 = Seq(codingseq2_adj)
    final_dnaseq_tot = Seq(f"{codingseq1}{codingseq2}")

    protein1 = final_dnaseq1.translate()
    protein2 = final_dnaseq2.translate()
    protein_tot = final_dnaseq_tot.translate()

    print(f"\nCoding sequence for protein 1: \n{codingseq1}")
    print(f"\nCoding sequence for protein 2: \n{codingseq2}")
    print(f"\nAdditional bases at the end of the 1st coding sequence: {rem1}")
    print(f"Additional bases at the beginning of the 2nd coding sequence: {rem2}")
    print(f"Necessary frameshift to maintain correct frame for 2nd protein: {first_element}")
    print(f"\nTranslated sequence for protein 1: \n{protein1}")
    print(f"\nTranslated sequence for protein 2: \n{protein2}")
    print(f"\nFinal fusion protein sequence: \n{protein_tot}")

if __name__ == "__main__":
    main()
























# gene1_UCSC_entry = get_gene_from_ucsc(genename=genename1, chrnum=chrnum1, start=(break_end1-1), end=(break_end1+1), assembly=assembly)
# check_for_MANE_zhan(assembly=assembly, genename=genename1, gene_id=gene1_UCSC_entry["name"].split('.')[0])
# exons1 = get_exon_list(break_end=break_end1, gene_UCSC_entry=gene1_UCSC_entry)
# # print(exons1)
# codingseq1 = get_seq_from_exons_N_term(exons=exons1, chrnum=chrnum1, assembly=assembly)


# gene2_UCSC_entry = get_gene_from_ucsc(genename=genename2, chrnum=chrnum2, start=(break_end2-1), end=(break_end2+1), assembly=assembly)
# check_for_MANE_zhan(assembly=assembly, genename=genename2, gene_id=gene2_UCSC_entry["name"].split('.')[0])
# exons2 = get_exon_list(break_end=break_end2, gene_UCSC_entry=gene2_UCSC_entry)
# # print(exons2)
# seq2check2 = get_lastresid2check(exons=exons2, target_name=gene2_UCSC_entry["name"], assembly=assembly)
# codingseq2 = get_seq_from_exons_C_term(exons=exons2, chrnum=chrnum2, seq2check=seq2check2, assembly=assembly)
# codingseq2, first_element = check_correct_frame(N_term_codingseq=codingseq1, C_term_codingseq=codingseq2, seq2check=seq2check2)

# rem1 = len(codingseq1) % 3
# rem2 = len(codingseq2) % 3

# codingseq1_adj = codingseq1[:rem1]
# final_dnaseq1 = Seq(codingseq1)
# codingseq2_adj = f"{codingseq1[-rem1:]}{codingseq2}"
# final_dnaseq2 = Seq(codingseq2_adj)
# final_dnaseq_tot = Seq(f"{codingseq1}{codingseq2}")

# protein1 = final_dnaseq1.translate()
# protein2 = final_dnaseq2.translate()
# protein_tot = final_dnaseq_tot.translate()

# print(f"\nCoding sequence for protein 1: \n{codingseq1}")
# print(f"\nCoding sequence for protein 2: \n{codingseq2}")
# print(f"\nAdditional bases at the end of the 1st coding sequence: {rem1}")
# print(f"Additional bases at the beginning of the 2nd coding sequence: {rem2}")
# print(f"Necessary frameshift to maintain correct frame for 2nd protein: {first_element}")
# print(f"\nTranslated sequence for protein 1: \n{protein1}")
# print(f"\nTranslated sequence for protein 2: \n{protein2}")
# print(f"\nFinal fusion protein sequence: \n{protein_tot}")














