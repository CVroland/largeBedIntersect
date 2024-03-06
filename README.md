# largeBedIntersect : a simple tool to intersect large bed files

largeBedIntersect is a simple tool that processes the intersection of two large BED files. It can be used for BEDPE files as well. Return the regions of the first file that intersect with at least one region of the second file. It is the equivalent of bedtools intersect -wa. The main interest of this tool is that it can handle very large files and BEDPE files. For human genome it requires about 3.6Go of ram.  
Because it is possible to save an index of the second file, the intersection can be done in a very short time. The tool is written in python and only need the pandas and numpy libraries.

## Usage

```text
python largeBedIntersect.py intersect -a <bedA> -b <bedB> [-o <output>] [--indexed] [--strand] [--type first_first|second_second|first_second|second_first] [--aType bed|bedpe] [--bType bed|bedpe] 
```

If you have to intersect many times the same file B with different file A, you can save the index of the file B and use it for the next intersections.

```text
python largeBedIntersect.py index <bedB> [-o <index>] [--strand] [--bType bed|bedpe] [--btype first|second|both] 
```

This will save the index of the file B in a numpy compressed npz file. You can then use this index for the next intersections by using the --indexed option.

## Evolutions

- [ ] add the option tp compute the coverage of the intersection
- [ ] simulate the equivalent of bedtools intersect without the -wa option