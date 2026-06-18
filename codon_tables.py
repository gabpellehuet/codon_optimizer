"""Bundled Kazusa Codon Usage Database tables for common expression-host organisms.

Source: https://www.kazusa.or.jp/codon/ -- raw codon usage tables, parsed and
converted to DNA triplets / three-letter amino acid codes to match the
CODON_TABLE_TEXT format used by reverse_translation.parse_codon_table.
Bundled locally so the GUI works fully offline -- no network dependency.
"""

from reverse_translation import CODON_TABLE_TEXT as _E_COLI_TABLE

S_CEREVISIAE = """
AmAcid  Codon      Number    /1000     Fraction

Phe     TTT       170666.00       26.1      0.59
Ser     TCT       153557.00       23.5      0.26
Tyr     TAT       122728.00       18.8      0.56
Cys     TGT        52903.00        8.1      0.63
Phe     TTC       120510.00       18.4      0.41
Ser     TCC        92923.00       14.2      0.16
Tyr     TAC        96596.00       14.8      0.44
Cys     TGC        31095.00        4.8      0.37
Leu     TTA       170884.00       26.2      0.28
Ser     TCA       122028.00       18.7      0.21
End     TAA         6913.00        1.1      0.47
End     TGA         4447.00        0.7      0.30
Leu     TTG       177573.00       27.2      0.29
Ser     TCG        55951.00        8.6      0.10
End     TAG         3312.00        0.5      0.23
Trp     TGG        67789.00       10.4      1.00
Leu     CTT        80076.00       12.3      0.13
Pro     CCT        88263.00       13.5      0.31
His     CAT        89007.00       13.6      0.64
Arg     CGT        41791.00        6.4      0.14
Leu     CTC        35545.00        5.4      0.06
Pro     CCC        44309.00        6.8      0.15
His     CAC        50785.00        7.8      0.36
Arg     CGC        16993.00        2.6      0.06
Leu     CTA        87619.00       13.4      0.14
Pro     CCA       119641.00       18.3      0.42
Gln     CAA       178251.00       27.3      0.69
Arg     CGA        19562.00        3.0      0.07
Leu     CTG        68494.00       10.5      0.11
Pro     CCG        34597.00        5.3      0.12
Gln     CAG        79121.00       12.1      0.31
Arg     CGG        11351.00        1.7      0.04
Ile     ATT       196893.00       30.1      0.46
Thr     ACT       132522.00       20.3      0.35
Asn     AAT       233124.00       35.7      0.59
Ser     AGT        92466.00       14.2      0.16
Ile     ATC       112176.00       17.2      0.26
Thr     ACC        83207.00       12.7      0.22
Asn     AAC       162199.00       24.8      0.41
Ser     AGC        63726.00        9.8      0.11
Ile     ATA       116254.00       17.8      0.27
Thr     ACA       116084.00       17.8      0.30
Lys     AAA       273618.00       41.9      0.58
Arg     AGA       139081.00       21.3      0.48
Met     ATG       136805.00       20.9      1.00
Thr     ACG        52045.00        8.0      0.14
Lys     AAG       201361.00       30.8      0.42
Arg     AGG        60289.00        9.2      0.21
Val     GTT       144243.00       22.1      0.39
Ala     GCT       138358.00       21.2      0.38
Asp     GAT       245641.00       37.6      0.65
Gly     GGT       156109.00       23.9      0.47
Val     GTC        76947.00       11.8      0.21
Ala     GCC        82357.00       12.6      0.22
Asp     GAC       132048.00       20.2      0.35
Gly     GGC        63903.00        9.8      0.19
Val     GTA        76927.00       11.8      0.21
Ala     GCA       105910.00       16.2      0.29
Glu     GAA       297944.00       45.6      0.70
Gly     GGA        71216.00       10.9      0.22
Val     GTG        70337.00       10.8      0.19
Ala     GCG        40358.00        6.2      0.11
Glu     GAG       125717.00       19.2      0.30
Gly     GGG        39359.00        6.0      0.12
"""

H_SAPIENS = """
AmAcid  Codon      Number    /1000     Fraction

Phe     TTT       714298.00       17.6      0.46
Ser     TCT       618711.00       15.2      0.19
Tyr     TAT       495699.00       12.2      0.44
Cys     TGT       430311.00       10.6      0.46
Phe     TTC       824692.00       20.3      0.54
Ser     TCC       718892.00       17.7      0.22
Tyr     TAC       622407.00       15.3      0.56
Cys     TGC       513028.00       12.6      0.54
Leu     TTA       311881.00        7.7      0.08
Ser     TCA       496448.00       12.2      0.15
End     TAA        40285.00        1.0      0.30
End     TGA        63237.00        1.6      0.47
Leu     TTG       525688.00       12.9      0.13
Ser     TCG       179419.00        4.4      0.05
End     TAG        32109.00        0.8      0.24
Trp     TGG       535595.00       13.2      1.00
Leu     CTT       536515.00       13.2      0.13
Pro     CCT       713233.00       17.5      0.29
His     CAT       441711.00       10.9      0.42
Arg     CGT       184609.00        4.5      0.08
Leu     CTC       796638.00       19.6      0.20
Pro     CCC       804620.00       19.8      0.32
His     CAC       613713.00       15.1      0.58
Arg     CGC       423516.00       10.4      0.18
Leu     CTA       290751.00        7.2      0.07
Pro     CCA       688038.00       16.9      0.28
Gln     CAA       501911.00       12.3      0.27
Arg     CGA       250760.00        6.2      0.11
Leu     CTG      1611801.00       39.6      0.40
Pro     CCG       281570.00        6.9      0.11
Gln     CAG      1391973.00       34.2      0.73
Arg     CGG       464485.00       11.4      0.20
Ile     ATT       650473.00       16.0      0.36
Thr     ACT       533609.00       13.1      0.25
Asn     AAT       689701.00       17.0      0.47
Ser     AGT       493429.00       12.1      0.15
Ile     ATC       846466.00       20.8      0.47
Thr     ACC       768147.00       18.9      0.36
Asn     AAC       776603.00       19.1      0.53
Ser     AGC       791383.00       19.5      0.24
Ile     ATA       304565.00        7.5      0.17
Thr     ACA       614523.00       15.1      0.28
Lys     AAA       993621.00       24.4      0.43
Arg     AGA       494682.00       12.2      0.21
Met     ATG       896005.00       22.0      1.00
Thr     ACG       246105.00        6.1      0.11
Lys     AAG      1295568.00       31.9      0.57
Arg     AGG       486463.00       12.0      0.21
Val     GTT       448607.00       11.0      0.18
Ala     GCT       750096.00       18.4      0.27
Asp     GAT       885429.00       21.8      0.46
Gly     GGT       437126.00       10.8      0.16
Val     GTC       588138.00       14.5      0.24
Ala     GCC      1127679.00       27.7      0.40
Asp     GAC      1020595.00       25.1      0.54
Gly     GGC       903565.00       22.2      0.34
Val     GTA       287712.00        7.1      0.12
Ala     GCA       643471.00       15.8      0.23
Glu     GAA      1177632.00       29.0      0.42
Gly     GGA       669873.00       16.5      0.25
Val     GTG      1143534.00       28.1      0.46
Ala     GCG       299495.00        7.4      0.11
Glu     GAG      1609975.00       39.6      0.58
Gly     GGG       669768.00       16.5      0.25
"""

M_MUSCULUS = """
AmAcid  Codon      Number    /1000     Fraction

Phe     TTT       422153.00       17.2      0.44
Ser     TCT       398250.00       16.2      0.20
Tyr     TAT       298518.00       12.2      0.43
Cys     TGT       279729.00       11.4      0.48
Phe     TTC       535439.00       21.8      0.56
Ser     TCC       444041.00       18.1      0.22
Tyr     TAC       394074.00       16.1      0.57
Cys     TGC       301384.00       12.3      0.52
Leu     TTA       165150.00        6.7      0.07
Ser     TCA       289799.00       11.8      0.14
End     TAA        23403.00        1.0      0.28
End     TGA        40148.00        1.6      0.49
Leu     TTG       329668.00       13.4      0.13
Ser     TCG       103815.00        4.2      0.05
End     TAG        19126.00        0.8      0.23
Trp     TGG       306619.00       12.5      1.00
Leu     CTT       329757.00       13.4      0.13
Pro     CCT       450637.00       18.4      0.31
His     CAT       260637.00       10.6      0.41
Arg     CGT       114854.00        4.7      0.08
Leu     CTC       495018.00       20.2      0.20
Pro     CCC       446868.00       18.2      0.30
His     CAC       375626.00       15.3      0.59
Arg     CGC       229758.00        9.4      0.17
Leu     CTA       198032.00        8.1      0.08
Pro     CCA       423707.00       17.3      0.29
Gln     CAA       293318.00       12.0      0.26
Arg     CGA       161412.00        6.6      0.12
Leu     CTG       969515.00       39.5      0.39
Pro     CCG       151521.00        6.2      0.10
Gln     CAG       836320.00       34.1      0.74
Arg     CGG       250836.00       10.2      0.19
Ile     ATT       377698.00       15.4      0.34
Thr     ACT       335039.00       13.7      0.25
Asn     AAT       382284.00       15.6      0.43
Ser     AGT       311331.00       12.7      0.15
Ile     ATC       552184.00       22.5      0.50
Thr     ACC       465115.00       19.0      0.35
Asn     AAC       499149.00       20.3      0.57
Ser     AGC       483013.00       19.7      0.24
Ile     ATA       180467.00        7.4      0.16
Thr     ACA       391437.00       16.0      0.29
Lys     AAA       537723.00       21.9      0.39
Arg     AGA       297135.00       12.1      0.22
Met     ATG       559953.00       22.8      1.00
Thr     ACG       138180.00        5.6      0.10
Lys     AAG       825270.00       33.6      0.61
Arg     AGG       299472.00       12.2      0.22
Val     GTT       262535.00       10.7      0.17
Ala     GCT       491093.00       20.0      0.29
Asp     GAT       515049.00       21.0      0.45
Gly     GGT       280522.00       11.4      0.18
Val     GTC       377902.00       15.4      0.25
Ala     GCC       637878.00       26.0      0.38
Asp     GAC       638504.00       26.0      0.55
Gly     GGC       520069.00       21.2      0.33
Val     GTA       182733.00        7.4      0.12
Ala     GCA       388723.00       15.8      0.23
Glu     GAA       661498.00       27.0      0.41
Gly     GGA       411344.00       16.8      0.26
Val     GTG       696158.00       28.4      0.46
Ala     GCG       157124.00        6.4      0.09
Glu     GAG       965963.00       39.4      0.59
Gly     GGG       372099.00       15.2      0.23
"""

P_PASTORIS = """
AmAcid  Codon      Number    /1000     Fraction

Phe     TTT         1963.00       24.1      0.54
Ser     TCT         1983.00       24.4      0.29
Tyr     TAT         1300.00       16.0      0.47
Cys     TGT          626.00        7.7      0.64
Phe     TTC         1675.00       20.6      0.46
Ser     TCC         1344.00       16.5      0.20
Tyr     TAC         1473.00       18.1      0.53
Cys     TGC          356.00        4.4      0.36
Leu     TTA         1265.00       15.6      0.16
Ser     TCA         1234.00       15.2      0.18
End     TAA           69.00        0.8      0.51
End     TGA           27.00        0.3      0.20
Leu     TTG         2562.00       31.5      0.33
Ser     TCG          598.00        7.4      0.09
End     TAG           40.00        0.5      0.29
Trp     TGG          834.00       10.3      1.00
Leu     CTT         1289.00       15.9      0.16
Pro     CCT         1282.00       15.8      0.35
His     CAT          960.00       11.8      0.57
Arg     CGT          564.00        6.9      0.17
Leu     CTC          620.00        7.6      0.08
Pro     CCC          553.00        6.8      0.15
His     CAC          737.00        9.1      0.43
Arg     CGC          175.00        2.2      0.05
Leu     CTA          873.00       10.7      0.11
Pro     CCA         1540.00       18.9      0.42
Gln     CAA         2069.00       25.4      0.61
Arg     CGA          340.00        4.2      0.10
Leu     CTG         1215.00       14.9      0.16
Pro     CCG          320.00        3.9      0.09
Gln     CAG         1323.00       16.3      0.39
Arg     CGG          158.00        1.9      0.05
Ile     ATT         2532.00       31.1      0.50
Thr     ACT         1820.00       22.4      0.40
Asn     AAT         2038.00       25.1      0.48
Ser     AGT         1020.00       12.5      0.15
Ile     ATC         1580.00       19.4      0.31
Thr     ACC         1175.00       14.5      0.26
Asn     AAC         2168.00       26.7      0.52
Ser     AGC          621.00        7.6      0.09
Ile     ATA          906.00       11.1      0.18
Thr     ACA         1118.00       13.8      0.24
Lys     AAA         2433.00       29.9      0.47
Arg     AGA         1634.00       20.1      0.48
Met     ATG         1517.00       18.7      1.00
Thr     ACG          491.00        6.0      0.11
Lys     AAG         2748.00       33.8      0.53
Arg     AGG          539.00        6.6      0.16
Val     GTT         2188.00       26.9      0.42
Ala     GCT         2351.00       28.9      0.45
Asp     GAT         2899.00       35.7      0.58
Gly     GGT         2075.00       25.5      0.44
Val     GTC         1210.00       14.9      0.23
Ala     GCC         1348.00       16.6      0.26
Asp     GAC         2103.00       25.9      0.42
Gly     GGC          655.00        8.1      0.14
Val     GTA          804.00        9.9      0.15
Ala     GCA         1228.00       15.1      0.23
Glu     GAA         3043.00       37.4      0.56
Gly     GGA         1550.00       19.1      0.33
Val     GTG          998.00       12.3      0.19
Ala     GCG          314.00        3.9      0.06
Glu     GAG         2360.00       29.0      0.44
Gly     GGG          468.00        5.8      0.10
"""

B_SUBTILIS = """
AmAcid  Codon      Number    /1000     Fraction

Phe     TTT        24450.00       30.0      0.68
Ser     TCT        10320.00       12.7      0.20
Tyr     TAT        18967.00       23.3      0.65
Cys     TGT         2939.00        3.6      0.46
Phe     TTC        11677.00       14.3      0.32
Ser     TCC         6766.00        8.3      0.13
Tyr     TAC        10290.00       12.6      0.35
Cys     TGC         3472.00        4.3      0.54
Leu     TTA        16167.00       19.8      0.21
Ser     TCA        11874.00       14.6      0.23
End     TAA         1562.00        1.9      0.61
End     TGA          621.00        0.8      0.24
Leu     TTG        12914.00       15.8      0.16
Ser     TCG         5266.00        6.5      0.10
End     TAG          381.00        0.5      0.15
Trp     TGG         8765.00       10.7      1.00
Leu     CTT        17772.00       21.8      0.23
Pro     CCT         8640.00       10.6      0.28
His     CAT        12832.00       15.7      0.68
Arg     CGT         5903.00        7.2      0.18
Leu     CTC         8697.00       10.7      0.11
Pro     CCC         2850.00        3.5      0.09
His     CAC         6141.00        7.5      0.32
Arg     CGC         6720.00        8.2      0.20
Leu     CTA         3975.00        4.9      0.05
Pro     CCA         5796.00        7.1      0.19
Gln     CAA        16620.00       20.4      0.52
Arg     CGA         3537.00        4.3      0.10
Leu     CTG        18757.00       23.0      0.24
Pro     CCG        13316.00       16.3      0.44
Gln     CAG        15089.00       18.5      0.48
Arg     CGG         5664.00        6.9      0.17
Ile     ATT        29487.00       36.2      0.49
Thr     ACT         7082.00        8.7      0.16
Asn     AAT        18702.00       22.9      0.56
Ser     AGT         5577.00        6.8      0.11
Ile     ATC        22167.00       27.2      0.37
Thr     ACC         7342.00        9.0      0.17
Asn     AAC        14522.00       17.8      0.44
Ser     AGC        11766.00       14.4      0.23
Ile     ATA         7960.00        9.8      0.13
Thr     ACA        17642.00       21.6      0.40
Lys     AAA        39449.00       48.4      0.70
Arg     AGA         8561.00       10.5      0.25
Met     ATG        21424.00       26.3      1.00
Thr     ACG        12110.00       14.9      0.27
Lys     AAG        16997.00       20.8      0.30
Arg     AGG         3341.00        4.1      0.10
Val     GTT        15193.00       18.6      0.28
Ala     GCT        15150.00       18.6      0.24
Asp     GAT        27108.00       33.2      0.64
Gly     GGT        10566.00       13.0      0.19
Val     GTC        14067.00       17.3      0.26
Ala     GCC        13433.00       16.5      0.22
Asp     GAC        15485.00       19.0      0.36
Gly     GGC        18967.00       23.3      0.34
Val     GTA        10638.00       13.0      0.20
Ala     GCA        17205.00       21.1      0.28
Glu     GAA        39217.00       48.1      0.68
Gly     GGA        17743.00       21.8      0.31
Val     GTG        14105.00       17.3      0.26
Ala     GCG        16176.00       19.8      0.26
Glu     GAG        18429.00       22.6      0.32
Gly     GGG         9094.00       11.2      0.16
"""

CODON_TABLES = {
    "Escherichia coli K12": _E_COLI_TABLE,
    "Saccharomyces cerevisiae": S_CEREVISIAE,
    "Homo sapiens": H_SAPIENS,
    "Mus musculus": M_MUSCULUS,
    "Pichia pastoris": P_PASTORIS,
    "Bacillus subtilis": B_SUBTILIS,
}
