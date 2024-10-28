from Bio.Seq import Seq
import csv
import subprocess
from glob import glob
import os

# Globals
anarciexe = 'ANARCI'
outpath = './03-AnnotatedResults/'

os.system('mkdir -p ' + outpath)

def is_valid_sequence(seq):
    """Check if the sequence contains only valid nucleotide characters."""
    return all(base in 'ATCG' for base in seq.upper())

def process_file(chain):
    """Process ASV records through ANARCI software tool and write new table to outpath."""
    tsvs = glob('./02-Results/*_' + chain + '.tsv')
    assert (len(tsvs) > 0), "NO TSV file found, check that R scripts ran."
    assert (len(tsvs) == 1), "More than one TSV file was found for " + './02-Results/*_' + chain + '.tsv'
    inf = tsvs[0]
    infDR = csv.DictReader(open(inf, 'r'), delimiter='\t')
    outf = open(os.path.join(outpath, os.path.basename(inf)), 'w')
    outfDW = csv.DictWriter(outf, delimiter='\t', restval='-', fieldnames=infDR.fieldnames + 
                                ['chain_type', 'e-value', 'score', 'seqstart_index', 'seqend_index', 'scheme', 'frame', 'AA', 'numbering', 'domain'])
    outfDW.writeheader()

    for record in infDR:
        seq = Seq(record['ASV'].strip())
        if not is_valid_sequence(str(seq)):
            print(f"Skipping invalid sequence: {seq}")
            continue

        # Build translations of each frame until one has a predicted domain
        for i in range(3):
            try:
                aa = (seq[i:] + ("N" * (3 - len(seq[i:]) % 3))).translate()
            except CodonTable.TranslationError as e:
                print(f"Translation error for sequence {seq}: {e}")
                continue
            
            anarci = run_anarci(aa, i)
            if anarci is not None:
                record.update(anarci)
                break
        outfDW.writerow(record)
    outf.close()



def run_anarci(aa):
    """Run ANARCI program and parse results into a dictionary. Return FIRST hit."""
    result = None
    subAAs = str(aa).strip('X').split('*')  # Parse amino acid sequence, removing terminal 'X's and splitting by '*'

    for i, s in enumerate(subAAs):
        if len(s) > 0:
            # Define command and run ANARCI
            cmd = "anarci --scheme imgt -i " + s
            test = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            output = test.stdout.splitlines()
            
            # Check and parse output
            if len(output) > 5:
                try:
                    mcols = output[4].split('|')
                    mdata = output[5].split('|')
                    result = {
                        'chain_type': mdata[mcols.index('chain_type')],
                        'e-value': mdata[mcols.index('e-value')],
                        'score': mdata[mcols.index('score')],
                        'seqstart_index': int(mdata[mcols.index('seqstart_index')]),
                        'seqend_index': int(mdata[mcols.index('seqend_index')]),
                        'scheme': output[6].split('=')[1].strip(),
                        'AA': '*'.join(subAAs),
                        'numbering': ','.join([l.split()[1] for l in output[7:] if l != '//']),
                        'domain': ','.join([l.split()[-1] for l in output[7:] if l != '//'])
                    }

                    # Highlight the annotated sequence portion
                    splitAA = s[0:result['seqstart_index']]
                    splitAA += '`' + s[result['seqstart_index']:result['seqend_index']] + '`'
                    splitAA += s[result['seqend_index']:]
                    subAAs[i] = splitAA
                    result['AA'] = '*'.join(subAAs)

                except (IndexError, ValueError) as e:
                    print(f"Parsing error: {e} in ANARCI output for sequence: {s}")
                break
            else:
                print(f"Unexpected ANARCI output format for sequence: {s}")
                print("ANARCI Output:", output)
                print("ANARCI Errors:", test.stderr)

    return result



# Run for each chain type
for chain in ['HeavyChain', 'LightChain']:
    process_file(chain)
