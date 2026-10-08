Recommended read lengths and site offsets for translation initiation or
termination peaks, derived from this library's own metagene profiles.

For each read length the report gives the distance from the mapped read end to
the first nucleotide of the P-site codon for RIBO/TIS start peaks, or the A-site
codon for TTS stop peaks. TTS advice also includes derived P-site distances.
The adjacent orfbounder/ directory contains a JSON manifest and, for each end
with usable advice, sample-keyed read-length and offset JSONs. Five-prime
exports negate the reported distance; three-prime exports keep it positive.
RIBO/TIS exports use the directly estimated P-site, and TTS exports use the
directly estimated A-site. Each end retains its own recommendation. These
files prepare inputs only; HRIBO does not execute ORFBounder. The report shows
peak enrichment, reading frames, and periodicity. A recommendation is only
made when the corresponding boundary signal supports one; otherwise the
report explains why no setup is suggested.
