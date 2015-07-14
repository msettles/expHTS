import re
import os


def flashParser(f, h, d):
    if not os.path.isfile(f):
        return

    f = open(f, "r")
    strSearch = ["[FLASH]     Total pairs:", "[FLASH]     Discarded pairs:", "[FLASH]     Combined pairs:", "[FLASH]         Innie pairs:", "[FLASH]         Outie pairs:", "[FLASH]     Uncombined pairs:", "[FLASH]     Percent combined:"]

    header = ["Flash_Pairs", "Flash_Discarded", "Combined_Pairs", "Innie_Combined", "Innie_Combined_Percent",  "Outie_Combined", "Outie_Combined_Percent", "Uncombined", "Percent_Combined"]

    data = []
    lines = f.readlines()
    if lines == []:
        return

    index = 0
    for e in lines:
        if strSearch[index] in e:
            if strSearch[index] == "[FLASH]     Percent combined:":
                data.append(float(re.findall("\d+.\d+", e)[-1]))
            else:
                data.append(int(re.search(r'\d+', e).group()))

            if "[FLASH]         Innie pairs:" == strSearch[index] or "[FLASH]         Outie pairs:" == strSearch[index]:
                data.append(float(re.findall("\d+.\d+", e)[-1]))

            index += 1
            if index >= len(strSearch):
                break
    h += header
    d += data


def sickleParser(f, h, d):
    if not os.path.isfile(f):
        return

    f = open(f, "r")

    strSearch = ["FastQ paired records kept:", "FastQ single records kept:", "FastQ paired records discarded:", "FastQ single records discarded:", "PE1 Base pairs left removed:",
        "PE1 Base pairs right removed:", "PE2 Base pairs left removed", "PE2 Base pairs right removed:"]

    header = ["Sickle_Pairs_Kept", "Sickle_Single_Kept", "Sickle_Pairs_Discarded", "Sickle_Single_Discarded", "Sickle_R1_BP_Left_Removed", "Sickle_R1_BP_Right_Removed", "Sickle_R2_BP_Left_Removed", "Sickle_R2_BP_Right_Removed"]
    headerPolyA = ["Sickle_R1_PolyAT_Trim", "Sickle_R1_PolyAT_Trim"]

    data = []
    dataPolyA = []

    lines = f.readlines()
    index = 0
    if lines == []:
        return

    for e in lines:
        if strSearch[index] in e:
            e = e.split(":")[1]
            data.append(int(re.search(r'\d+', e).group()))
            index += 1
            if index >= len(strSearch):
                break
        if "Poly AT tails" in e:
            e = e.split(":")[1]
            dataPolyA.append(int(re.search(r'\d+', e).group()))

    h += header
    d += data

    if len(dataPolyA) != 0:
        h += headerPolyA
        d += dataPolyA


def scythParser(f, h, d):
    if not os.path.isfile(f):
        return
    f = open(f, "r")

    strSearch = ["contaminated: "]
    header = ["Scythe_contaminated", "Scythe_no_contamination", "Scythe_total"]
    data = []

    lines = f.readlines()
    if lines == []:
        return

    index = 0

    for e in lines:
        if strSearch[index] in e:
            info = e.split(',')
            data.append(int(re.search(r'\d+', info[0]).group()))
            data.append(int(re.search(r'\d+', info[1]).group()))
            data.append(int(re.search(r'\d+', info[2]).group()))
            break

    h += header
    d += data


def filterParser(f, h, d):
    if not os.path.isfile(f):
        return
    f = open(f, "r")

    strSearch = ["PE_written:", "SE_written:", "Discarded:"]

    header = ["Screening_PE_written", "Screening_SE_written", "Screening_discarded"]
    data = []

    lines = f.readlines()
    index = 0
    if lines == []:
        return

    lines = lines[-1]
    lines = lines.split("|")

    for e in lines:
        if strSearch[index] in e:
            data.append(int(re.search(r'\d+', e).group()))
            index += 1
            if index >= len(strSearch):
                break

    h += header
    d += data


def deduperParser(f, h, d):
    #Final:| reads: 4998890 | duplicates: 393582 | percent: 7.87 | discarded: 1110 | total_seconds: 63.79 | reads/sec: 783
    if not os.path.isfile(f):
        return
    f = open(f, "r")  # open the file

    strSearch = ["reads:", "duplicates:", "percent:", "discarded:"]

    header = ["Deduper_Reads", "Deduper_Duplicates", "Deduper_Percent_Duplicates", "Deduper_Discarded"]
    data = []

    lines = f.readlines()  # read in each line

    if lines == []:
        return

    lines = lines[-1]  # remove the \n at the end of line
    lines = lines.split("|")[1:5]  # split by pipe char
    index = 0

    for e in lines:  # for each element in line
        e = e.strip()
        data.append(e.split(' ')[1])

    h += header
    d += data

def finalParser(f, h, d):
    if not os.path.isfile(f):
        return
    f = open(f, "r")

    lines = f.readlines()

    if lines == []:
        return
    header = lines[0].split("\t")
    data = lines[1].split("\t")
    h += header
    d += data


def htseqCountParse(f, h, d):
    if not os.path.isfile(f):
        return
    f = open(f, 'r')
    lines = f.readlines()
    if lines == []:
        return

    header = ["features", "zero_count_features", "in_feature", "no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique", "percent_in_feature"]
    data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    index = 3
    for e in lines:
        e = e.split("\t")
        if '__' in e[0]:
            data[index] = int(e[1])
            index += 1
        elif e[1] == 0:
            data[0] += 1
            data[1] += 1
        else:
            data[0] += 1
            data[2] + int(e[1])

    data[9] = float(data[2]/sum(data[2:9]))

    h += header
    d += data


def flagstatsParse(f, h, d):
    if not os.path.isfile(f):
        return
    f = open(f, "r")
    lines = f.readlines()
    if lines == []:
        return

    header = ["Total", "Secondary", "Supplementary", "Duplicates", "Mapped", "Mapped_Percent", "Paired","Read1", "Read2", "Properly_Paired", "Properly_Paired_Percent", "With_itself_and_Mate_mapped", "Singletons", "Singletons_Percent", "With_Mate_Mapped_To_A_Different_chr", "With_Mate_mapped_to_a_different_chr_Mapq>=5"]
    data = []

    for e in lines:  # for each element in line
                data.append(float(re.findall("\d+", e)[0]))
                tmp = re.findall("\d+.\d+%", e)
                if tmp != []:
                    data.append(tmp[0])

    h += header
    d += data


def printToFile(out, header, data):
    f = open(out, "w")

    f.write('\t'.join(map(str, header)))
    f.write('\t'.join(map(str, data)))

    f.close()


def parseOut(base, sample):
    import os
    data = []
    header = []

    header.append("Sample")
    data.append(sample)

    filterParser(os.path.join(base, "PE_filter_info.log"), header, data)
    deduperParser(os.path.join(base, "PE_deduper_info.log"), header, data)
    sickleParser(os.path.join(base, "PE_sickle_info.log"), header, data)
    flashParser(os.path.join(base, "flash_info.log"), header, data)
    finalParser(os.path.join(base, "finalCleanup.log"), header, data)

    printToFile(os.path.join(base, sample + "_SummaryStats.log"), header, data)

    return os.path.join(base, sample + "_SummaryStats.log")


def parseOutMapping(base, sample):
    import os
    data = []
    header = []

    header.append("Sample")
    data.append(sample)

    flagstatsParse(os.path.join(base, sample + ".flagstats"), header, data)

    printToFile(os.path.join(base, sample + "_MappingSummaryStats.log"), header, data)

    return os.path.join(base, sample + "_MappingSummaryStats.log")


def parseOutHTseq(base, sample):
    import os
    data = []
    header = []

    header.append("Sample")
    data.append(sample)

    htseqCountParse(os.path.join(base, sample + ".counts"), header, data)

    printToFile(os.path.join(base, sample + "_CountingSummaryStats.log"), header, data)

    return os.path.join(base, sample + "_CountingSummaryStats.log")


def bringTogether(listFiles, out):
    first = 0
    outFile = open(out, "w")

    for e in listFiles:
        f = open(e, "r")
        lines = f.readlines()
        if first == 0:
            outFile.write(lines[0])
            if lines[0][-1] == '\n':
                outFile.write('\n')

            outFile.write(lines[1])
            if lines[1][-1] == '\n':
                outFile.write('\n')
        else:
            outFile.write(lines[1])
            if lines[1][-1] == '\n':
                outFile.write('\n')
        first += 1
        f.close()

    outFile.close()
