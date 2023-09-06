# Writing a GFF file parser

Let's start by looking at the GFF file format.
Because these are meant to be a flexible file format, it contains some required information and lots of optional information. The rows include:

1. **seqid**  The name of the sequence where the feature is located
2. **source** The procedure or algorithm used to generate the feature
3. **type** The type of feature
4. **start** The starting base position (1-based)
5. **end** The ending base position (1-based, inclusive)
6. **score** Numeric value generally indicating the confidence of the source in this feature (. for null value)
7. **strand** Single character indicating the strand of the feature: +, -, . for undertermined, ? for unknown
8. **phase** The phase of CDS features, 0, 1, or 2 (for CDS features) or . for everything else
9. **attributes** A list of tag-value pairs separated by semicolons with additional information

---

## And example entry

Here is the first line of one of the homework files:

    NC_060925.1
    Curated Genomic
    pseudogene
    144134
    146717
    .
    -
    .
    ID=gene-SEPTIN14P14;Dbxref=GeneID:107105354,HGNC:HGNC:51701;Name=SEPTIN14P14;description=septin 14 pseudogene 14;gbkey=Gene;gene=SEPTIN14P14;gene_biotype=pseudogene;gene_synonym=SEPT14P14;pseudo=true

---

## Which fields are relevant?

In order to focus on the information we are interested in, we first need to filter out lines that aren't of the correct type of feature.

Where will we find that information?

---

## Which fields are relevant?

In order to focus on the information we are interested in, we first need to filter out lines that aren't of the correct type of feature.

Where will we find that information?

    !python
    for line in open(fname):
        line = line.rstrip().split("\t")
        if line[2] != "pseudogene":
            continue

---

## How long?

We also want to know how long is each pseudogene. How can we get that information? 


---

## How long?

We also want to know how long is each pseudogene. How can we get that information? 


    !python
    for line in open(fname):
        line = line.rstrip().split("\t")
        if line[2] != "pseudogene":
            continue
        length = int(line[6]) - int(line[5])

---

## How can we break up the attributes field?

Finally, we need to separate all of the attributes. Since we know that they should be in tag-value pairs, this seems like an ideal case for a dictionary, especially since there is no gaurentee that each line will have all the same attributes.

To begin with, we need a break apart the attributes field. We know that each pair is separated by a semicolon so let's make a list using that as the delimiter.

    !python
    fields = line[8].split(";")

Now if we look at fields[0], we see it is equal to "ID=gene-SEPTIN14P14"

---

## Going from a list to a dictionary

How do we handle working on one list item at a time? A *for* loop! And inside the loop, we need to break up the tag-value pairs into their two parts.

    !python
    info = {}
    for pair in fields:
        key, value = pair.split("=")
        info[key] = value

---

## Keeping track of one entry per gene

The gene name is one of the attributes under the tag "gene". But before we add it to our master dictionary, we need to check if there has already been an entry with the gene name (there are multiple occurances of some pseudogenes). Also, if it is a duplicate, we want to keep track of how many copies and the largest size.

To do this, we can use a conditional statement. Assuming that we have a dictionary named *genes* to keep track of gene entries...

    !python
    name = info['gene']
    if name not in genes:               # We've never seen this gene
        genes[name] = info
        genes[name]['count'] = 1        # We need to be able to keep count
        genes[name]['length'] = length  # We need to keep track of lengths
    else:
        genes[name]['count'] += 1       # We add another sighting
        # Only keep the biggest length
        genes[name]['length'] = max(genes[name]['length'], length)

---

## Putting it all together

    !python
    genes = {}
    for line in open(fname):
        line = line.rstrip().split("\t")
        if line[2] != "pseudogene":
            continue
        length = int(line[6]) - int(line[5])
        info = {}
        for pair in fields:
            key, value = pair.split("=")
            info[key] = value
        name = info['gene']
        if name not in genes:
            genes[name] = info
            genes[name]['count'] = 1
            genes[name]['length'] = length
        else:
            genes[name]['count'] += 1
            genes[name]['length'] = max(genes[name]['length'], length)
