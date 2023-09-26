### MY Graph500 Generator for 32bits Vertex ID

> This is a modified version of the Graph500 generator. It is modified to generate graphs with 32bits vertex ID. 
> Also called Kronecker generator.

To use this generator, 
1. cd `text2bin` and `make`
2. cd `src` and `make generator32`
3. `./src/generator32 N edgefactor [outputfile, e.g. Kron25_16.txt]` to generate .txt file
4. `.text2bin/text2bin [inputfile] [outputfile]` to convert .txt file to .bin file
