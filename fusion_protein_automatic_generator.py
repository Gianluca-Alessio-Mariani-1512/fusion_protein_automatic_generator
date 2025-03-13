import requests
from Bio.Seq import Seq
import argparse
import sys

def get_input_from_vcf(data_file):
    """
    Get input data for fusion protein reconstruction from vcf file
    
    Args:
        data_file: .vcf file path
        
    Returns:
        assembly: reference genome (e.g. "hg19")
        breakpoint_list: tuple made of (gene name e.g. CLTC or ALK, chromosome number e.g. 17 or 2, exact base position of the breakpoint in the chromosome)
    """

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
        assembly = genome_reference_string[start_index:].strip()  # Rimuove eventuali spazi bianchi
        print("Reference Genome:", assembly)
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

    breakpoint_list = []

    for i in found_lines:
        # print(i)
        # Suddividi la stringa in base ai caratteri di tabulazione
        elements = i.split('\t')
        # Ottieni i primi due elementi
        chrnum_str = elements[0]
        chrnum = ''.join(char for char in chrnum_str if char.isdigit())
        break_point = elements[1]

        eighth_element = elements[7].split(';')

        for i in eighth_element:
            if i.startswith('GENE_NAME='):
            # Estrai il valore dopo "="
                gene_name = i.split('=')[1]
                print("Found Gene:", gene_name)
        if gene_name == -1:
            print("No Gene Found")

        # Stampa i risultati
        breakpoint_list.append((gene_name, chrnum, break_point))

    return assembly, breakpoint_list, variant_found

def get_UCSC_entry(genename, chrnum, start, end, assembly):
    """
    Get the UCSC entry of the gene by using its common name, genome assembly and position on the genome; only MANE select genes will be selected
    
    Args:
        genename: common gene name e.g. ALK
        chrnum: chromosome number
        start - end: base start and end position for query; every gene included in this region will be listed in the query output
        assembly: reference genome (e.g. "hg19")
    Returns:
        chosen_gene: UCSC entry for the chosen gene (e.g. "hg19")
    """

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
    """
    check if the selected gene is actually MANE select by comparing with MANE select list
    
    Args:
        genename: common gene name e.g. ALK
        gene_id: UCSC Genome gene code e.g. ENST00000017003
        assembly: reference genome (e.g. "hg19")
    Returns:
        Print 
    """

    url = f"https://github.com/zhanyinx/variantalker/tree/main/resources/mane.transcript.gencode.v43.{assembly}.txt"
    response = requests.get(url)
    MANE_text = response.text
    if gene_id in MANE_text:
        print(f"The gene {genename} identified with the UCSC Genome code {gene_id} is present in the MANE list in the file {url}")
    else:
        print(f"The gene {genename} identified with the UCSC Genome code {gene_id} is NOT present in the MANE list in the file {url}")

def get_prot_seq2check(gene_UCSC_entry):
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

def get_coding_sequence(gene_UCSC_entry, break_end, term_type, exons_Nterm=None, exons_Cterm=None):

    # Estrai e converti i dati
    chrnum = gene_UCSC_entry["chrom"]
    strand = gene_UCSC_entry["strand"]
    block_sizes = list(map(int, gene_UCSC_entry['blockSizes'].strip(',').split(',')))
    chrom_starts = list(map(int, gene_UCSC_entry['chromStarts'].strip(',').split(',')))

    exons_total = []
    for start, size in zip(chrom_starts, block_sizes):
        exon_start = start + gene_UCSC_entry['chromStart']
        exon_end = exon_start + size
        exons_total.append((exon_start, exon_end))

    if strand == "-":
        exons_total = exons_total[::-1]

    exons_tot_dict = {exons_total[i]:i+1 for i in range(len(exons_total))}
    # print(exons_tot_dict)

    if exons_Nterm is None and term_type == "N":
        exons = []
        for start, size in zip(chrom_starts, block_sizes):
            exon_start = start + gene_UCSC_entry['chromStart']
            exon_end = exon_start + size
            if term_type == "N" and strand == "+":
                if exon_end <= break_end +1:  # Solo esoni che terminano prima del breakpoint
                    exons.append((exon_start, exon_end))
            if term_type == "N" and strand == "-":
                if exon_start >= break_end -1:  # Solo esoni che terminano prima del breakpoint
                    exons.append((exon_start, exon_end))

        new_exon_dict = {}

        for i in exons:
            if i in exons_tot_dict:
                new_exon_dict[i] = exons_tot_dict[i]

    elif exons_Cterm is None and term_type == "C":
        exons = []
        for start, size in zip(chrom_starts, block_sizes):
            exon_start = start + gene_UCSC_entry['chromStart']
            exon_end = exon_start + size
            if term_type == "C" and strand == "+":
                if exon_start >= break_end -1:  # Solo esoni che terminano prima del breakpoint
                    exons.append((exon_start, exon_end))
            if term_type == "C" and strand == "-":
                if exon_end <= break_end +1:  # Solo esoni che terminano prima del breakpoint
                    exons.append((exon_start, exon_end))

        new_exon_dict = {}

        for i in exons:
            if i in exons_tot_dict:
                new_exon_dict[i] = exons_tot_dict[i]
    
    else:
        new_exon_dict = exons_tot_dict

    # print(new_exon_dict)
    new_exon_dict_inv = {value: key for key, value in new_exon_dict.items()}
    # print("inverted:", new_exon_dict_inv)

    if exons_Nterm is not None:
        if term_type == "N":
            nuovo_dizionario = {}
            for chiave, valore in new_exon_dict_inv.items():
                if chiave <= exons_Nterm:
                    nuovo_dizionario[chiave] = valore
            # print("no_exons_dict:", nuovo_dizionario)
            new_exon_dict_inv = nuovo_dizionario
    if exons_Cterm is not None:
        if term_type == "C":
            nuovo_dizionario = {}
            for chiave, valore in new_exon_dict_inv.items():
                if chiave >= exons_Cterm:
                    nuovo_dizionario[chiave] = valore
            # print("no_exons_dict:", nuovo_dizionario)
            new_exon_dict_inv = nuovo_dizionario

    max_exon = max(new_exon_dict_inv, key=new_exon_dict_inv.get)
    min_exon = min(new_exon_dict_inv, key=new_exon_dict_inv.get)
    range_string = f"{max_exon}-{min_exon}"

    ordered_exons = [new_exon_dict_inv[key] for key in sorted(new_exon_dict_inv.keys())]

    # print(exons)
    codingseq_fragments = []

    if term_type == "N" and strand == "+":
        first_exon = ordered_exons.pop(0)
        first_exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=first_exon[0], end=first_exon[1])
        trimmed_start_exon = trim_start_cod(dna_seq=first_exon_seq)
        codingseq_fragments.append(str(trimmed_start_exon))
        for i in range(len(ordered_exons)):
            exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=ordered_exons[i][0], end=ordered_exons[i][1])
            codingseq_fragments.append(str(exon_seq))
        coding_seq = "".join(codingseq_fragments)
        return coding_seq, range_string
    
    if term_type == "C" and strand == "+":
        last_exon = ordered_exons.pop(-1)
        for i in range(len(ordered_exons)):
            exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=ordered_exons[i][0], end=ordered_exons[i][1])
            codingseq_fragments.append(str(exon_seq))
        last_exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=last_exon[0], end=last_exon[1])
        trimmed_end_exon = trim_stop_cod(dna_seq=last_exon_seq, seq2check=get_prot_seq2check(gene_UCSC_entry=gene_UCSC_entry))
        codingseq_fragments.append(str(trimmed_end_exon))
        coding_seq = "".join(codingseq_fragments)
        return coding_seq, range_string
    
    if term_type == "N" and strand == "-":
        first_exon = ordered_exons.pop(0)
        first_exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=first_exon[0], end=first_exon[1])
        dna_seq = Seq(first_exon_seq)
        complement_seq = dna_seq.complement()
        reverse_complement_seq = complement_seq[::-1]
        trimmed_start_exon = trim_start_cod(dna_seq=reverse_complement_seq)
        codingseq_fragments.append(str(trimmed_start_exon))
        for i in range(len(ordered_exons)):
            exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=ordered_exons[i][0], end=ordered_exons[i][1])
            dna_seq = Seq(exon_seq)
            complement_seq = dna_seq.complement()
            reverse_complement_seq = complement_seq[::-1]
            codingseq_fragments.append(str(reverse_complement_seq))
        coding_seq = "".join(codingseq_fragments)
        return coding_seq, range_string
    
    if term_type == "C" and strand == "-":
        last_exon = ordered_exons.pop(-1)
        for i in range(len(ordered_exons)):
            exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=ordered_exons[i][0], end=ordered_exons[i][1])
            dna_seq = Seq(exon_seq)
            complement_seq = dna_seq.complement()
            reverse_complement_seq = complement_seq[::-1]
            codingseq_fragments.append(str(reverse_complement_seq))
        last_exon_seq = get_sequence(chrnum=chrnum, assembly="hg19", start=last_exon[0], end=last_exon[1])
        dna_seq = Seq(last_exon_seq)
        complement_seq = dna_seq.complement()
        reverse_complement_seq = complement_seq[::-1]
        trimmed_end_exon = trim_stop_cod(dna_seq=reverse_complement_seq, seq2check=get_prot_seq2check(gene_UCSC_entry=gene_UCSC_entry))
        codingseq_fragments.append(str(trimmed_end_exon))
        coding_seq = "".join(codingseq_fragments)
        return coding_seq, range_string

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

def check_correct_frame(N_term_dna_seq, C_term_dna_seq, seq2check):
    codingseq_frames = [C_term_dna_seq, f"A{C_term_dna_seq}", f"AA{C_term_dna_seq}"]
    codingseq_frames_prot = []
    couples = []

    for index, i in enumerate(codingseq_frames):
        final_dnaseq = Seq(f"{N_term_dna_seq}{i}")
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
    parser = argparse.ArgumentParser(description="Automatic script to get a fusion protein sequence from a vcf file")
    parser.add_argument('-p', '--path', type=str, required=False, 
                        help="Path of the vcf input file")
    
    parser.add_argument('-bN', '--break_Nterm', type=int, required=False, 
                        help="Manually override the Nterminal protein breakpoint, Nterminal protein ends in specified base, specified base included")
    
    parser.add_argument('-bC', '--break_Cterm', type=int, required=False, 
                        help="Manually override the Cterminal protein breakpoint, Cterminal protein starts in specified base, specified base included")
    
    parser.add_argument('-e', '--exon_count', action='store_true', required=False, 
                        help="Enable to only get the two proteins exon count and nothing else")
    
    parser.add_argument('-eN', '--exons_Nterm', type=int, required=False, 
                        help="Manually override the Nterminal exon count, take all exons from first to the specified exon, specified exon included")
    
    parser.add_argument('-eC', '--exons_Cterm', type=int, required=False, 
                        help="Manually override the Cterminal exon count, take all exons from specified to the last exon, specified exon included")
    
    parser.add_argument('-o', '--output', type=str, required=False, 
                        help="Manually override the output name")
    
    args = parser.parse_args()
    
    # data_file = "data/22-C-004462_v1_22-C-004462_RNA_v1_Non-Filtered_2025-02-20_03_33_25.vcf"
    data_file = args.path

    if args.exons_Nterm:
        namestring_exons_Nterm = f"_Nto_Exon{str(args.exons_Nterm)}"
    else:
        namestring_exons_Nterm = ""

    if args.exons_Cterm:
        namestring_exons_Cterm = f"_Cfrom_Exon{str(args.exons_Cterm)}"
    else:
        namestring_exons_Cterm = ""

    vcf_output = get_input_from_vcf(data_file=data_file)
    assembly = str(vcf_output[0])
    variant_found = str(vcf_output[2])

    genenameN = vcf_output[1][0][0]
    chrnumN=int(vcf_output[1][0][1])

    if args.break_Nterm:
        break_endN = args.break_Nterm
        namestring_break_endN = f"_breakN{str(args.break_Nterm)}"
    else:
        break_endN = int(vcf_output[1][0][2])
        namestring_break_endN = ""

    genenameC = vcf_output[1][1][0]
    chrnumC=int(vcf_output[1][1][1])

    if args.break_Cterm:
        break_endC = args.break_Cterm
        namestring_break_endC = f"_breakC{str(args.break_Cterm)}"
    else:
        break_endC = int(vcf_output[1][1][2])
        namestring_break_endC = ""

    genesN_entry = get_UCSC_entry(chrnum=chrnumN, start=break_endN - 1, end=break_endN + 1, assembly=assembly, genename=genenameN)
    # print(genesN_entry)
    genesC_entry = get_UCSC_entry(chrnum=chrnumC, start=break_endC - 1, end=break_endC + 1, assembly=assembly, genename=genenameC)
    # print(genesC_entry)

    gene_id_baseN = genesN_entry['name'].split('.')[0]
    gene_id_baseC = genesC_entry['name'].split('.')[0]

    check_for_MANE_zhan(assembly=assembly, gene_id=gene_id_baseN, genename=genenameN)
    check_for_MANE_zhan(assembly=assembly, gene_id=gene_id_baseC, genename=genenameC)

    codingseq1, rangeN = get_coding_sequence(gene_UCSC_entry=genesN_entry, break_end=break_endN, term_type="N", exons_Nterm=args.exons_Nterm)
    codingseq2, rangeC = get_coding_sequence(gene_UCSC_entry=genesC_entry, break_end=break_endC, term_type="C", exons_Cterm=args.exons_Cterm)

    if args.exon_count:
        print(f"The {genenameN} gene exon count is:", genesN_entry["blockCount"])
        print(f"The {genenameC} gene exon count is:", genesC_entry["blockCount"])
        print(f"With the provided vcf input file the selected {genenameN} exon in the final fusion protein will be:", rangeN)
        print(f"With the provided vcf input file the selected {genenameC} exon in the final fusion protein will be:", rangeC)
        sys.exit(0)

    rem1 = len(codingseq1) % 3
    rem2 = len(codingseq2) % 3

    codingseq2, first_element = check_correct_frame(N_term_dna_seq=codingseq1, C_term_dna_seq=codingseq2, seq2check=get_prot_seq2check(gene_UCSC_entry=genesC_entry))

    codingseq1_adj = codingseq1[:-rem1] if rem1 != 0 else codingseq1
    final_dnaseq1 = Seq(codingseq1_adj)
    codingseq2_adj = f"{codingseq1[-rem1:]}{codingseq2}" if rem1 != 0 else codingseq2
    final_dnaseq2 = Seq(codingseq2_adj)
    final_dnaseq_tot = Seq(f"{codingseq1}{codingseq2}")

    protein1 = final_dnaseq1.translate()
    protein2 = final_dnaseq2.translate()
    protein_tot = final_dnaseq_tot.translate()

    if args.output:
        outputname = f"{args.output}.txt"
    else:
        outputname = f"{variant_found}{namestring_exons_Nterm}{namestring_exons_Cterm}{namestring_break_endN}{namestring_break_endC}___FusionProteinSequence.txt"

    with open(outputname, 'w') as file:
        file.write("Output file initialized\n")

    with open(outputname, 'a') as file:
        file.write(f"\nCoding sequence for protein 1: \n")
        file.write(f"{codingseq1}\n")

        file.write(f"\nCoding sequence for protein 2: \n")
        file.write(f"{codingseq2}\n")

        file.write(f"\nAdditional bases at the end of the 1st coding sequence: {rem1}")
        file.write(f"\nAdditional bases at the beginning of the 2nd coding sequence: {rem2}")
        file.write(f"\nNecessary frameshift to maintain correct frame for 2nd protein: {first_element}\n")

        file.write(f"\nTranslated sequence for protein 1: \n")
        file.write(f"{str(protein1)}\n")
        file.write(f"\nTranslated sequence for protein 2: \n")
        file.write(f"{str(protein2)}\n")
        file.write(f"\nFinal fusion protein sequence: \n")
        file.write(f"{str(protein_tot)}\n")

    print(f"\nCoding sequence for protein 1: \n{codingseq1}\nCoding sequence 1 length: {len(codingseq1)}")
    print(f"\nCoding sequence for protein 2: \n{codingseq2}\nCoding sequence 2 length: {len(codingseq2)}")
    print(f"\nAdditional bases at the end of the 1st coding sequence: {rem1}")
    print(f"Additional bases at the beginning of the 2nd coding sequence: {rem2}")
    print(f"Necessary frameshift to maintain correct frame for 2nd protein: {first_element}")
    print(f"\nTranslated sequence for protein 1: \n{protein1}")
    print(f"\nTranslated sequence for protein 2: \n{protein2}")
    print(f"\nFinal fusion protein sequence: \n{protein_tot}")

    print(f"\n#####\n\nOutput data saved in {variant_found}___FusionProteinSequence.txt\n\n#####")

if __name__ == "__main__":
    main()