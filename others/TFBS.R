library(data.table)
tfbs = fread("~/biodata/ENCODE/ENCODE.tf.bound.union.bed", header = FALSE)
tfbs

## Bed format <https://bedtools.readthedocs.io/en/latest/content/general-usage.html>
## For example, start=9, end=20 is interpreted to span bases 10 through 20,inclusive.
##
## the midpoint of each TFBS was determined by averaging the start
## and end positions of the binding site and, if necessary, rounding up to the nearest
## integer
tfbs$midpoint = round((tfbs$V2 + 1 + tfbs$V3) / 2)
tfbs

tfbs_midpoint = copy(tfbs)
tfbs_midpoint$V2 = tfbs_midpoint$midpoint - 1
tfbs_midpoint$V3 = tfbs_midpoint$midpoint

tfbs_midpoint$midpoint = NULL
tfbs_midpoint

fwrite(tfbs_midpoint, file = "~/biodata/ENCODE/tfbs_midpoint.union.bed.gz", sep = "\t", col.names = FALSE)
