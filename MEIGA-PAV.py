## DEPENDENCIES ##
# External
import os
import sys
import argparse
import copy
import re

# Internal
import formats
import alignment
import unix
import sequences
import gRanges
import annotation
import bamtools
import retrotransposons

###############
## Functions ##
###############


def search4polyA(sequence):
    '''
    Search for poly(A) at sequence end 
    '''
    ## Configuration for monomere search:
    windowSize = 8
    maxWindowDist = 2
    minMonomerSize = 10
    minPurity = 80  

    ## Seach poly(A) monomers
    targetMonomer = 'A'
    monomersA = sequences.find_monomers(sequence, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    if not monomersA:
        return False, None
        
    ## Select monomer closest to sequence end
    candidate = monomersA[-1]

    ## Filter out monomer if more than Xbp from end
    seqLen = len(sequence)
    dist2end = seqLen - candidate.end

    if dist2end <= 30:
        polyA = True

    else:
        polyA = False 

    return polyA, candidate

def search4polyT(sequence):
    '''
    Search for poly(T) at sequence begin 
    '''
    ## Configuration for monomere search:
    windowSize = 8
    maxWindowDist = 2
    minMonomerSize = 10
    minPurity = 80  

    ## Seach poly(T) monomers
    targetMonomer = 'T'
    monomersT = sequences.find_monomers(sequence, targetMonomer, windowSize, maxWindowDist, minMonomerSize, minPurity)

    if not monomersT:
        return False, None

    ## Select monomer closest to sequence beg        
    candidate = monomersT[0]

    ## Filter out monomer if more than Xbp from beg
    dist2beg = candidate.beg

    if dist2beg <= 30:
        polyT = True

    else:
        polyT = False 

    return polyT, candidate

def call_MEI_candidate(VCF):
    '''
    '''
    ## Create VCF with candidate MEI calls
    candidateVCF = formats.VCF()
    candidateVCF.header = VCF.header

    ## Create VCF with non-candidate MEI
    filteredVCF = formats.VCF()
    filteredVCF.header = VCF.header

    ## For each variant
    for variant in VCF.variants:
        
        ## Filter calls distinct from INS
        if variant.info['SVTYPE'] != 'INS':
            continue

        ## Search for poly(A) tail at inserted sequence
        polyA, monomerA = search4polyA(variant.alt)

        ## Search for poly(T) tail at inserted sequence
        polyT, monomerT = search4polyT(variant.alt)

        ## a) Filter calls if poly(A) nor poly(T) found
        if not polyA and not polyT:
            filteredVCF.add(variant)
        
        ## b) Add variant passing all the filters
        else:
            candidateVCF.add(variant)

    return candidateVCF, filteredVCF

def call_MEI(vcf, consensus, reference, sourceDb, outDir):
    '''
    '''    
    ## 0. Create temporary folder
    tmpDir = outDir + '/tmp'
    unix.mkdir(tmpDir)

    ## 1. Write inserted sequences into fasta file
    fastaPath = tmpDir + '/MEI_candidate.fa'
    fasta = ins2fasta(vcf, tmpDir)
    fasta.write(fastaPath)

    ## 2. Create index for consensus sequences
    fileName = 'consensus'  
    consensusIndex = alignment.index_minimap2(consensus, fileName, tmpDir)

    ## 3. Align inserted sequences against consensus:
    PAF_path = alignment.alignment_minimap2(fastaPath, consensusIndex, 'hits2consensus', 1, tmpDir)
    PAF_consensus = formats.PAF()
    PAF_consensus.read(PAF_path)

    ## Temporary
    index="/Users/brodriguez/Research/References/Annotations/H.sapiens/hg38/Repetitive_dna/smallRNAs.mmi"
    PAF_path = alignment.alignment_minimap2(fastaPath, index, 'hits2small_MEI', 1, tmpDir)

    ## Align inserted sequences against the reference genome
    #SAM_path = alignment.alignment_bwa(fastaPath, reference, 'hits2genome', 1, tmpDir)
    #PAF_path = alignment.sam2paf(SAM_path, 'hits2genome', tmpDir)
    #PAF_genome = formats.PAF()
    #PAF_genome.read(PAF_path)

    ## 4. Generate single PAF objects per inserted sequence:
    PAFs_consensus = group_alignments(PAF_consensus)
    #PAFs_genome = group_alignments(PAF_genome)

    ## 5. Resolve structure for each insertion with matches on retrotransposon consensus sequences
    structures = {}

    for insId in PAFs_consensus:
        structures[insId] = MEI_structure(PAFs_consensus[insId], fasta.seqDict[insId])
        seqBeg, seqEnd = structures[insId]['CHAIN'].interval()

    ## 6. Resolve 3' partnered transductions
    structures = resolve_partnered_3prime(structures, fasta, reference, sourceDb, tmpDir)

    ## 6. Search for 5' partnered transductions
    structures = search4partnered_5prime(structures, fasta, reference, tmpDir)

    ## 7. Search for orphan transductions
    ## Remove resolved insertions
    #for insId in structures:
    #    if structures[insId]['PASS']:
    #        del PAFs_genome[insId]

    ## Do orphan transduction search
    #search4orphan(PAFs_genome, sourceDb, fasta) # TO FINISH LATER (Only two L1 orphan transductions so far..)

    ## 8. Generate output VCF containing MEI calls
    ## Create header for output dictionary
    outVCF = formats.VCF()
    outVCF.header = vcf.header

    ## Add MEI specific fields to the VCF header
    info2add = {'ITYPE': ['.', 'String', 'Type of insertion (solo, partnered,  orphan or NUMT)'], \
                '3PRIME': ['0', 'Flag', 'Partnered 3-prime transduction'], \
                '5PRIME': ['0', 'Flag', 'Partnered 5-prime transduction'], \
                'FAM': ['.', 'String', 'Repeat family'], \
                'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
                'RETRO_LEN': ['1', 'Integer', 'Inserted retrotransposon length'], \
                'TRUNCATION_5_LEN': ['1', 'Integer', 'Size of 5prime truncation'], \
                'TRUNCATION_3_LEN': ['1', 'Integer', 'Size of 3prime truncation'], \
                'INVERSION_LEN': ['1', 'Integer', '5-inversion length'], \
                'RETRO_COORD': ['.', 'String', 'Coordinates for inserted retrotransposon piece of sequence'], \
                'IS_FULL': ['0', 'Flag', 'Full length mobile element'], \
                'ORF1': ['0', 'Flag', 'ORF1 identified'], \
                'ORF2': ['0', 'Flag', 'ORF2 identified'], \
                'COMPETENT': ['0', 'Flag', 'Potential competent full L1 with intact ORFs'], \
                'TDCOORD_5PRIME': ['1', 'Integer', '5-prime transduced sequence coordinates'], \
                'TDCOORD_3PRIME': ['1', 'Integer', '3-prime transduced sequence coordinates'], \
                'TDLEN_5PRIME': ['1', 'Integer', '5-prime transduction length'], \
                'TDLEN_3PRIME': ['1', 'Integer', '3-prime transduction length'], \
                'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], \
                'MT_COORD': ['.', 'String', 'Coordinates for the piece of MT genome integrated']                    
                }

    outVCF.header.info.update(info2add)

    ## Select INS corresponding to MEI calls and add update info field with MEI features
    for variant in vcf.variants:
        insId = variant.chrom + ':' + str(variant.pos) 
                            
        # Discard unresolved inserted sequences
        if (insId not in structures) or ((insId in structures) and (structures[insId]['PASS'] is False)):
            continue

        variant2add = copy.deepcopy(variant)
        variant2add.info.update(structures[insId])
        outVCF.add(variant2add)

    ## 9. Do cleanup
    #unix.rm([tmpDir])

    return outVCF

def search4orphan(hits, sourceDb, fasta):
    '''
    TO FINISH LATER
    '''
    # For each partnered event
    for insId in hits:
    
        # For each hit
        for hit in hits[insId].alignments:

            hit.tName = 'chr' + hit.tName

            ## Filter out hits 
            if (hit.alignmentPerc() < 75) or (hit.MAPQ < 30) or (hit.tName not in sourceDb):
                continue            

def search4partnered_5prime(structures, fasta, reference, outDir):
    '''
    '''
    ## 1. Create Fasta with sequences to realign
    seq2realign = formats.FASTA()

    for insId in structures:

        # Discard if strand not determined
        if structures[insId]['STRAND'] is None:
            continue

        ## Extract unresolved 5' sequence if any
        qBeg, qEnd = structures[insId]['CHAIN'].interval()

        if structures[insId]['STRAND'] == '+':
            seq2realign.seqDict[insId] = fasta.seqDict[insId][:qBeg]

        else:
            seq2realign.seqDict[insId] = fasta.seqDict[insId][qEnd:]

    fastaPath = outDir + '/seq2realign.5prime.fasta'
    seq2realign.write(fastaPath)

    ## 2. Realign sequences on the reference with BWA-mem
    SAM_path = alignment.alignment_bwa(fastaPath, reference, 'hits2genome.5prime', 1, outDir)
    PAF_path = alignment.sam2paf(SAM_path, 'hits2genome.5prime', outDir)
    
    PAF = formats.PAF()
    PAF.read(PAF_path)

    ## 3. Make 5' transduction calls
    # For each hit
    for hit in PAF.alignments:
        
        hit.tName = 'chr' + hit.tName
        iRef, coord = hit.qName.split(':')
        iBeg = int(coord) - 500
        iEnd = int(coord) + 500

        ## Filter out hits 
        if (hit.alignmentPerc() < 75) or (hit.MAPQ < 30) or (iRef == hit.tName and gRanges.overlap(iBeg, iEnd, hit.tBeg, hit.tEnd)[0]):
            continue

        ## Make call
        structures[hit.qName]['ITYPE'] = 'partnered'
        structures[hit.qName]['5PRIME'] = True
        structures[hit.qName]['TDCOORD_5PRIME'] = hit.tName + ':' + str(hit.tBeg) + '-' +  str(hit.tEnd)
        structures[hit.qName]['TDLEN_5PRIME'] = hit.tEnd - hit.tBeg
                
    return structures
        

def resolve_partnered_3prime(structures, fasta, reference, sourceDb, outDir):
    '''
    '''
    ## 1. Create Fasta with sequences to realign
    seq2realign = formats.FASTA()
    pattern = re.compile("Partnered_[0-9]+")
    partneredDict = {}

    for insId in structures:

        # Discard solo
        if structures[insId]['ITYPE'] != 'partnered':
            continue

        ## Initialize partnered dict for the ins
        partneredDict[insId] = {}
        partneredDict[insId]['NB_PARTNERED'] = 0
        partneredDict[insId]['NB_RESOLVED'] = 0

        # For each hit
        for hit in structures[insId]['CHAIN'].alignments:

            # Discard if not partnered
            if not pattern.match(hit.tName):
                continue

            # Add candidate partnered sequence to the fasta
            seqId = insId + '|' + hit.tName
            seq = fasta.seqDict[insId][hit.qBeg:hit.qEnd]
            seq2realign.seqDict[seqId] = seq

            ## Update partnered dictionary
            partneredDict[insId]['NB_PARTNERED'] += 1       
            partneredDict[insId][hit.tName] = hit


    fastaPath = outDir + '/seq2realign.3prime.fasta'
    seq2realign.write(fastaPath)

    ## 2. Realign sequences on the reference with BWA-mem
    SAM_path = alignment.alignment_bwa(fastaPath, reference, 'hits2genome.3prime', 1, outDir)
    PAF_path = alignment.sam2paf(SAM_path, 'hits2genome.3prime', outDir)
    PAF = formats.PAF()
    PAF.read(PAF_path)

    ## 3. Add hit information to partnered transduction candidates
    hits = PAF.hits2dict()

    # For each partnered event
    for ID in hits:
    
        insId, tdId = ID.split('|')
        partneredDict[insId]['CYTOID'] = None

        # For each hit
        for hit in hits[ID]:

            hit.tName = 'chr' + hit.tName

            ## Check if it´s a partnered transduction from a known source element
            if (hit.tName in sourceDb) and sourceDb[hit.tName].collect_interval(hit.tBeg - 200, hit.tEnd + 200, 'ALL'):
                source = sourceDb[hit.tName].collect_interval(hit.tBeg - 200, hit.tEnd + 200, 'ALL')[0][0]
                partneredDict[insId]['CYTOID'] = source.optional['cytobandId']

            ## Filter out hits 
            if (hit.alignmentPerc() < 75 or hit.MAPQ < 30) and (partneredDict[insId]['CYTOID'] is None):
                continue
            
            ## Add hit information  
            partneredDict[insId]['NB_RESOLVED'] += 1
            partneredDict[insId][tdId].tName = hit.tName + ':' + str(hit.tBeg) + '-' + str(hit.tEnd)
            
    ## 4. Add transduction information
    for insId in partneredDict:

        # a) Make transduction call
        if (partneredDict[insId]['NB_PARTNERED'] > 0) and (partneredDict[insId]['NB_PARTNERED'] == partneredDict[insId]['NB_RESOLVED']):

            tdIds = [key for key in partneredDict[insId].keys() if key not in ['NB_PARTNERED', 'NB_RESOLVED', 'CYTOID']]
            structures[insId]['CYTOID'] = partneredDict[insId]['CYTOID']
            structures[insId]['TDCOORD_3PRIME'] = ','.join([partneredDict[insId][tdId].tName for tdId in tdIds])
            structures[insId]['TDLEN_3PRIME'] = ','.join([str(partneredDict[insId][tdId].qEnd - partneredDict[insId][tdId].qBeg) for tdId in tdIds])

        # b) Make solo call
        else:
            structures[insId]['ITYPE'] = 'solo'
    
    return structures


def MEI_structure(PAF, insertSeq):
    '''
    '''
    structure = {}
    structure['LEN'] = len(insertSeq)

    ## 1. Chain alignments
    structure['CHAIN'] = PAF.chain(20, 50)

    ## 2. Determine insertion family
    families = list(set([hit.tName.split('|')[1] for hit in structure['CHAIN'].alignments]))
    subfamilies = list(set([hit.tName.split('|')[2] for hit in structure['CHAIN'].alignments]))

    # a) L1 insertion
    if 'L1' in families:
        structure['FAM'] = 'L1'
        structure['SUBFAM'] = ','.join(subfamilies) 

    # b) SVA insertion
    elif 'SVA' in families:
        structure['FAM'] = 'SVA'
        structure['SUBFAM'] = ','.join(subfamilies) 

    # c) Alu insertion
    elif 'Alu' in families:
        structure['FAM'] = 'Alu'
        structure['SUBFAM'] = ','.join(subfamilies) 

    elif ('L1' in families) and ('Alu' in families):
        print('FUSION: ', insertSeq)

    ## 3. Search for polyA/T tails at unresolved insert ends and determine insertion type
    structure['ITYPE'] = 'solo'
    rtBeg, rtEnd = structure['CHAIN'].interval()

    ## Set parameters for monomer search
    windowSize = 8
    maxWindowDist = 2
    minMonomerSize = 10
    minPurity = 80  

    ### 3.1 PolyA search
    ## Search for monomers
    targetSeq = insertSeq[rtEnd:]
    monomersA = sequences.find_monomers(targetSeq, 'A', windowSize, maxWindowDist, minMonomerSize, minPurity)

    ## Map to insert sequence coordinates
    for monomer in monomersA:
        monomer.beg = monomer.beg + rtEnd
        monomer.end = monomer.end + rtEnd

    ## Make polyA calls
    structure['POLYA'] = 0
    structure['STRAND'] = None

    # a) Single polyA
    if (len(monomersA) == 1):
        dist2rt = monomersA[0].beg - rtEnd
        dist2end = structure['LEN'] - monomersA[0].end

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'solo'
            structure['POLYA'] = 1
            structure['STRAND'] = '+'

    # b) Multiple polyA (Transduction candidate)
    # AAAAAAAAAAAAAA-------------AAAAAAAAAAAAAA
    elif (len(monomersA) > 1):
        dist2rt = monomersA[0].beg - rtEnd
        dist2end = structure['LEN'] - monomersA[-1].end

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'partnered'
            structure['3PRIME'] = True
            structure['POLYA'] = len(monomersA)
            structure['STRAND'] = '+'

    ## 3.2 PolyT search 
    ## Search for monomers
    targetSeq = insertSeq[:rtBeg]
    monomersT = sequences.find_monomers(targetSeq, 'T', windowSize, maxWindowDist, minMonomerSize, minPurity)

    ## Make polyT calls
    structure['POLYT'] = 0

    # a) Single polyT
    if (len(monomersT) == 1):
        dist2end = monomersT[0].beg
        dist2rt = rtBeg - monomersT[0].end 

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'solo'
            structure['POLYT'] = 1
            structure['STRAND'] = '-'

    # b) Multiple polyT (Partnered transduction candidate)
    # TTTTTTTTTTTTT-------------TTTTTTTTTTTTT
    elif (len(monomersT) == 2):
        dist2end = monomersT[0].beg
        dist2rt = rtBeg - monomersT[-1].end 

        # Apply filter
        if (dist2rt <= 30) and (dist2end <= 30):
            structure['ITYPE'] = 'partnered'
            structure['3PRIME'] = True
            structure['POLYT'] = len(monomersT)
            structure['STRAND'] = '-'
        
    ## 3.3 Determine if polyA or polyT found
    # a) PolyA found
    if (structure['POLYA'] != 0) and (structure['POLYT'] == 0):

        tail = 'polyA'
        structure['NBPOLY'] = structure['POLYA']
        monomers = monomersA

    # b) PolyT found
    elif (structure['POLYT'] != 0) and (structure['POLYA'] == 0):
        tail = 'polyT'
        structure['NBPOLY'] = structure['POLYT']
        monomers = monomersT

    # c) No tail or ambiguous
    else:
        tail = None
        structure['NBPOLY'] = 0

    ## 3.4 Determine candidate insertion type based on the number of polyA/T tails found
    hits2add = []

    # a) Solo
    if structure['NBPOLY'] == 1:
        fields = [structure['CHAIN'].alignments[0].qName, structure['LEN'], monomers[0].beg, monomers[0].end, None, 'PolyA/T', 0, 0, 0, 0, 0, 0]
        hit = formats.PAF_alignment(fields)
        hits2add.append(hit)

    # b) Partnered
    elif structure['NBPOLY'] > 1:

        ## First polyA/T
        fields = [structure['CHAIN'].alignments[0].qName, structure['LEN'], monomers[0].beg, monomers[0].end, None, 'PolyA/T', 0, 0, 0, 0, 0, 0]
        hit = formats.PAF_alignment(fields)
        hits2add.append(hit)

        ## Add partnered region/s plus polyAT/s
        counter = 1

        for monomer1, monomer2 in zip(monomers, monomers[1:]):

            # Partnered 
            fields = [structure['CHAIN'].alignments[0].qName, structure['LEN'], monomer1.end, monomer2.beg, None, 'Partnered_' + str(counter), 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            hits2add.append(hit)

            # Next polyA/T
            fields = [structure['CHAIN'].alignments[0].qName, structure['LEN'], monomer2.beg, monomer2.end, None, 'PolyA/T', 0, 0, 0, 0, 0, 0]
            hit = formats.PAF_alignment(fields)
            hits2add.append(hit)

            counter += 1        
      
    ## 3.5 Add polyA/T plus transduced annotation to the chain
    if tail == 'polyA':
        structure['CHAIN'].alignments = structure['CHAIN'].alignments + hits2add

    elif tail == 'polyT':
        structure['CHAIN'].alignments = hits2add + structure['CHAIN'].alignments 

    ## 4. Infer inserted sequence length
    lengths = retrotransposons.infer_lengths('solo', structure['CHAIN'], structure['STRAND'])
    structure.update(lengths)

    ## 5. Assess ORFs status for L1 insertions
    if structure['FAM'] == 'L1':
        orf1, orf2 = retrotransposons.find_orf(insertSeq)[0:2]
        structure['ORF1'] = True if orf1 is not None else False
        structure['ORF2'] = True if orf2 is not None else False
        structure['COMPETENT'] = True if (structure['ORF1'] and structure['ORF2'] and structure['IS_FULL']) else False

    ## 6. Apply filters 
    failed = []

    # 6.1 Percentage resolved filter
    # Compute % of insertion resolved
    structure['PERC-RESOLVED'] = structure['CHAIN'].perc_query_covered()

    if structure['PERC-RESOLVED'] < 60:
        failed.append('PERC-RESOLVED')

    # 6.2 Length filtering for solo insertions
    if structure['ITYPE'] == 'solo':

        if (structure['FAM'] == 'L1') and (structure['LEN'] > 6500):
            failed.append('LEN')
        
        elif (structure['FAM'] == 'Alu') and (structure['LEN'] > 400): 
            failed.append('LEN')

        elif (structure['FAM'] == 'SVA') and (structure['LEN'] > 5000): 
            failed.append('LEN')

    structure['FAILED'] = failed

    # a) Insertions passes all the filters
    if not failed:
        structure['PASS'] = True
    
    # b) At least one failed filter
    else:
        structure['PASS'] = False

    return structure


def ins2fasta(vcf, outDir):
    '''
    Write inserted sequences into a fasta

    Input:
        1. vcf: Path to VCF file
        2. outDir: Output directory

    Output: 
        1. fastaPath: path to fasta object inserted sequences
    '''     
    ## 1. Initialize fasta object
    fasta = formats.FASTA() 

    ## 2. Collect inserted sequences
    for variant in vcf.variants:
 
        insId = variant.chrom + ':' + str(variant.pos) 
        seq = variant.alt

        fasta.seqDict[insId] = seq

    return fasta

def group_alignments(paf):
    '''

    '''     
    pafDict = {}

    ## For each hit
    for hit in paf.alignments:

        # Initialize paf object for this inserted sequence
        if hit.qName not in pafDict:
            pafDict[hit.qName] = formats.PAF()
    
        # Add hit to the corresponding paf
        pafDict[hit.qName].alignments.append(hit)

    return pafDict


def call_NUMT(vcf, mtGenome, outDir):
    '''
    '''    
    ## 0. Create temporary folder
    tmpDir = outDir + '/tmp'
    unix.mkdir(tmpDir)

    ## 1. Write inserted sequences into fasta file
    fastaPath = tmpDir + '/insertions.fa'
    fasta = ins2fasta(vcf, tmpDir)
    fasta.write(fastaPath)

    ## 2. Create index for the mitochondrial genome
    fileName = 'mtGenome'  
    mtIndex = alignment.index_minimap2(mtGenome, fileName, tmpDir)

    ## 3. Align inserted sequences against the mitochondrial genome
    PAF_path = alignment.alignment_minimap2(fastaPath, mtIndex, 'hits2mt', 1, tmpDir)
    PAF_mt = formats.PAF()
    PAF_mt.read(PAF_path)

    ## 4. Generate single PAF objects per inserted sequence:
    PAFs_mt = group_alignments(PAF_mt)

    ## 5. Make NUMTs calls
    NUMTs = {}

    for insId in PAFs_mt:
        chain = PAFs_mt[insId].chain(20, 50)

        # Make NUMT call if enough % of sequence resolved
        if chain.perc_query_covered() >= 60:

            coords = chain.interval_template() 

            NUMT = {}
            NUMT['ITYPE'] = 'NUMT'
            NUMT['MT_COORD'] = str(coords[0]) + '-' + str(coords[1])
            NUMTs[insId] = NUMT
    
    ## 6. Generate output VCF containing NUMT calls
    ## Create header for output dictionary
    outVCF = formats.VCF()
    outVCF.header = vcf.header

    ## Add MEI specific fields to the VCF header
    info2add = {'ITYPE': ['.', 'String', 'Type of insertion (solo, partnered or orphan)'], \
                '3PRIME': ['0', 'Flag', 'Partnered 3-prime transduction'], \
                '5PRIME': ['0', 'Flag', 'Partnered 5-prime transduction'], \
                'FAM': ['.', 'String', 'Repeat family'], \
                'CYTOID': ['.', 'String', 'Source element cytoband identifier'], \
                'RETRO_LEN': ['1', 'Integer', 'Inserted retrotransposon length'], \
                'TRUNCATION_5_LEN': ['1', 'Integer', 'Size of 5prime truncation'], \
                'TRUNCATION_3_LEN': ['1', 'Integer', 'Size of 3prime truncation'], \
                'INVERSION_LEN': ['1', 'Integer', '5-inversion length'], \
                'RETRO_COORD': ['.', 'String', 'Coordinates for inserted retrotransposon piece of sequence'], \
                'IS_FULL': ['0', 'Flag', 'Full length mobile element'], \
                'ORF1': ['0', 'Flag', 'ORF1 identified'], \
                'ORF2': ['0', 'Flag', 'ORF2 identified'], \
                'COMPETENT': ['0', 'Flag', 'Potential competent full L1 with intact ORFs'], \
                'TDCOORD_5PRIME': ['1', 'Integer', '5-prime transduced sequence coordinates'], \
                'TDCOORD_3PRIME': ['1', 'Integer', '3-prime transduced sequence coordinates'], \
                'TDLEN_5PRIME': ['1', 'Integer', '5-prime transduction length'], \
                'TDLEN_3PRIME': ['1', 'Integer', '3-prime transduction length'], \
                'STRAND': ['.', 'String', 'Insertion DNA strand (+ or -)'], 
                }
    outVCF.header.info.update(info2add)

    ## Select INS corresponding to MEI calls and add update info field with MEI features
    for variant in vcf.variants:
        insId = variant.chrom + ':' + str(variant.pos) 

        # Discard unresolved inserted sequences
        if (insId not in NUMTs):
            continue


        variant2add = copy.deepcopy(variant)
        variant2add.info.update(NUMTs[insId])
        outVCF.add(variant2add)

    ## 9. Do cleanup
    #unix.rm([tmpDir])

    return outVCF


######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('vcf', help='Path to VCF file with assembly-based SV calls')
parser.add_argument('consensus', help='Path to FASTA file containing consensus retrotransposon sequences')
parser.add_argument('reference', help='Path to FASTA file containing the reference genome')
parser.add_argument('mtGenome', help='Path to FASTA file containing the mitochondrial genome')
parser.add_argument('fileName', help='Output file name')
parser.add_argument('-o', '--outDir', default=os.getcwd(), dest='outDir', help='output directory. Default: current working directory' )

## 2. Parse user input ##
args = parser.parse_args()
vcf = args.vcf
consensus = args.consensus
reference = args.reference
mtGenome = args.mtGenome
fileName = args.fileName
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]


print()
print('***** ', scriptName, 'configuration *****')
print('vcf: ', vcf)
print('consensus: ', consensus)
print('reference: ', reference)
print('mtGenome: ', mtGenome)
print('fileName: ', fileName)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########
## Note: NUMT detection is disabled

## 1. Read VCF
VCF = formats.VCF()
VCF.read(vcf)

## 2. Load source elements database

annotDir = '/Users/brodriguez/Research/Projects/HGSVC2/Analysis/Source_L1/V2/data/'
annotations = annotation.load_annotations(['TRANSDUCTIONS'], VCF.header.refLengths, annotDir, None, 1, outDir)

## 2. Filter VCF by selecting retrotransposition insertion candidates 
# (inserted sequences with polyA/T tails at their ends)
candidateVCF, filteredVCF = call_MEI_candidate(VCF)

## 3. Search for NUMTs
#NUMT_VCF = call_NUMT(filteredVCF, mtGenome, outDir)

## 4. Do MEI calling for candidate insertions
MEI_VCF = call_MEI(candidateVCF, consensus, reference, annotations['TRANSDUCTIONS'], outDir)

## 5. Merge NUMT and MEI VCFs
outVCF = MEI_VCF
#outVCF.variants = MEI_VCF.variants + NUMT_VCF.variants

## 6. Write VCF containing MEI calls
infoIds = ['VARTYPE', 'SVTYPE', 'SVLEN', 'ID', 'LEAD_SAMPLE', 'TIG_REGION', 'TIG_STRAND', 'ITYPE', '3PRIME', '5PRIME', 'FAM', 'CYTOID', 'RETRO_LEN', 'TRUNCATION_5_LEN', 'TRUNCATION_3_LEN', 'INVERSION_LEN', 'RETRO_COORD', 'IS_FULL', 'ORF1', 'ORF2', 'COMPETENT', 'TDCOORD_5PRIME', 'TDCOORD_3PRIME', 'TDLEN_5PRIME', 'TDLEN_3PRIME', 'STRAND']
	
#infoIds = ['VARTYPE', 'SVTYPE', 'SVLEN', 'ID', 'LEAD_SAMPLE', 'TIG_REGION', 'TIG_STRAND', 'ITYPE', '3PRIME', '5PRIME', 'FAM', 'CYTOID', 'RETRO_LEN', 'TRUNCATION_5_LEN', 'TRUNCATION_3_LEN', 'INVERSION_LEN', 'RETRO_COORD', 'IS_FULL', 'ORF1', 'ORF2', 'COMPETENT', 'TDCOORD', 'TDLEN', 'STRAND', 'MT_COORD']
formatIds = ['GT']

outVCF.write(infoIds, formatIds, fileName, outDir)
