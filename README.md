# mark-nonconverted-reads
This program examines a BAM or SAM stream for reads that contain multiple Cs in non-CpG context. It markes these reads with a tag (XX:Z:UC) and sets the vendor failed bit to prevent these reads from being used in downstream methylation calling. It also prints summary information to stderr.
