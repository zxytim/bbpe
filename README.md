# bbpe: Byte BPE

This project is created to study the byte-pair-encoding starting from byte to
support multiple modalities.

After later found that BPE is actually just a pretty weak lossless compressor
(even inferior than zip), while Transformer + SGD is also a lossless compressor,
but is much much better, I deprecated the idea of "vocabulary" and turn into a
"nocabularist" .

Nevertheless, this is a quite efficient  O(|corpus| + |vocab|)  and accuate implementation of BPE (without
approximations used by sentencepiece). 

Please inspect `./bbpe_v7_efficient.cpp` for details.
