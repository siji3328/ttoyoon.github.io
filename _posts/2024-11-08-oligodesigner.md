---
layout: post
title: "Oligo designer"
date: 2024-11-08
categories: [Bioinformatics, Tool]
tags: [Oliognucleotide, Oiligo designer]
---

# PCA Primer Design Script

This Python script designs primers based on a provided DNA sequence. It uses the Wallace formula to calculate the melting temperature (Tm) of DNA fragments and filters fragments based on specified length and Tm criteria.

## Code

```python
def calculate_tm(sequence):
    # Wallace 공식: Tm = 2(A+T) + 4(G+C)
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    tm = 2 * (a_count + t_count) + 4 * (g_count + c_count)
    return tm

def design_pca_primers(sequence, min_fragment_length=75, max_fragment_length=95, min_tm=50, max_tm=55):
    primers = []
    fragments = []

    i = 0
    while i < len(sequence):
        selected_fragment = None
        overlap_tm = 0
        for fragment_length in range(max_fragment_length, min_fragment_length - 1, -1):
            if i + fragment_length > len(sequence):
                fragment_length = len(sequence) - i

            fragment = sequence[i:i + fragment_length]

            if i + fragment_length < len(sequence):
                overlap_sequence = sequence[i + fragment_length - (fragment_length - min_fragment_length):i + fragment_length]
                overlap_tm = calculate_tm(overlap_sequence)

                if min_tm <= overlap_tm <= max_tm:
                    selected_fragment = fragment
                    i += fragment_length - (fragment_length - min_fragment_length)
                    break
            else:
                selected_fragment = fragment
                i += fragment_length

        if selected_fragment:
            fragments.append((selected_fragment, overlap_tm))

    # 짝수 개의 조각을 유지하기 위해 조정
    if len(fragments) % 2 != 0:
        fragments.pop()

    # 프라이머 설계 (Forward 및 Reverse)
    for i, (fragment, overlap_tm) in enumerate(fragments):
        forward_primer = fragment
        if i < len(fragments) - 1:
            reverse_primer = fragments[i + 1][0][:min_fragment_length]  # 다음 조각의 오버랩 서열을 사용
        else:
            reverse_primer = ""  # 마지막 조각의 Reverse Primer는 필요하지 않음

        primers.append({
            "fragment_number": i + 1,
            "fragment_sequence": fragment,
            "fragment_length": len(fragment),
            "overlap_tm": overlap_tm,
            "forward_primer": forward_primer,
            "reverse_primer": reverse_primer
        })

    return primers

# 제공된 전체 DNA 서열
sequence = "ATGAGCAAAGGTGAAGAACTGTTTACCGGCGTTGTGCCGATTCTGGTGGAACTGGATGGCGATGTGAACGGTCACAAATTCAGCGTGCGTGGTGAAGGTGAAGGCGATGCCACGATTGGCAAACTGACGCTGAAATTTATCTGCACCACCGGCAAACTGCCGGTGCCGTGGCCGACGCTGGTGACCACCCTGACCTATGGCGTTCAGTGTTTTAGTCGCTATCCGGATCACATGAAACGTCACGATTTCTTTAAATCTGCAATGCCGGAAGGCTATGTGCAGGAACGTACGATTAGCTTTAAAGATGATGGCAAATATAAAACGCGCGCCGTTGTGAAATTTGAAGGCGATACCCTGGTGAACCGCATTGAACTGAAAGGCACGGATTTTAAAGAAGATGGCAATATCCTGGGCCATAAACTGGAATACAACTTTAATAGCCATAATGTTTATATTACGGCGGATAAACAGAAAAATGGCATCAAAGCGAATTTTACCGTTCGCCATAACGTTGAAGATGGCAGTGTGCAGCTGGCAGATCATTATCAGCAGAATACCCCGATTGGTGATGGTCCGGTGCTGCTGCCGGATAATCATTATCTGAGCACGCAGACCGTTCTGTCTAAAGATCCGAACGAAAAACGGGACCACATGGTTCTGCACGAATATGTGAATGCGGCAGGTATTACGTGGAGCCATCCGCAGTTCGAAAAATAA"

# 프라이머 설계
primers = design_pca_primers(sequence, min_fragment_length=75, max_fragment_length=95, min_tm=50, max_tm=54)

# 결과 출력
for primer in primers:
    print(f"Fragment {primer['fragment_number']}:")
    print(f"  Sequence: {primer['fragment_sequence']}")
    print(f"  Fragment Length: {primer['fragment_length']} bp")
    print(f"  Overlap Tm: {primer['overlap_tm']} °C")
    print(f"  Forward Primer: {primer['forward_primer']}")
    print(f"  Reverse Primer: {primer['reverse_primer']}")
    print("\n")
```

## Description
This script calculates the melting temperature (Tm) of DNA fragments using the Wallace formula and designs primers for PCR based on specified criteria. It iterates through a given DNA sequence to find suitable fragments that meet length and Tm requirements.
