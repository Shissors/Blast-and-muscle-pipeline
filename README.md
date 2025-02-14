# BLAST-and-Muscle-pipeline
## Description
Takes a protein from a sequence performs blastp against another sequence from the database.
Aligns top hits with each other using muscle.
Makes a phylogenetic tree of those hits.
Puts it all into a nice report. (i should add cute characters to the report).


## Requirements
###Packages
[Python]([url](https://www.python.org/)) and something to run on (like VScode), you could also run this through a command prompt I think.
Packages:
[Biopython](https://biopython.org/)
[Matplotlib](https://matplotlib.org/)
[blast.exe](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
[muscle.exe](https://drive5.com/muscle/downloads_v3.htm)

### Data
A huge-ass fasta file that contains sequences that will be used for the blast search
Your query sequence, also a fasta (the sequence can also be within a huge database). Note: when you are uploading your own query make sure they follow the standard naming format like this: sp|A0A075B706|TRDJ1_HUMAN T

###This image(just keep it)
![cat](https://th.bing.com/th/id/OIP.bGYHV9PXpVn4CFAPFnzmbgHaD2?rs=1&pid=ImgDetMain)


## Running the program
use this
```
protname [query] -blastcom [blastp/n] -musclecom muscle.exe -fastafile [database containing query] -databasefasta [blast search database]
```
you can add multiple queries, separating each by a comma.


## Additional files
I have the example input files in the repository.
The reports are the outputs(PDF files).

### Things i dont know
If this works with blastn 
Why brown cows don't produce chocolate milk
![real](https://i0.wp.com/onpasture.com/wp-content/uploads/2018/01/Cow-on-bike.png?fit=375%2C293&ssl=1)



