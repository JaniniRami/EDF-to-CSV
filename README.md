# EDF-to-CSV

Convert your EDF (European Data Format) files into CSV (comma-separated value) files with ease!

## Why?

Ever tried reading an EDF file with 19 channels using pyedflib, only to be thwarted by German characters? That's exactly what happened to me! So I created my own EDF reader.

## Usage
```
python edf_to_csv.py <EDF File Path > <Is Digital Output?>
```

Example:
```
python edf_to_csv.py "/Users/ramijanini/Desktop/CVP-OSA-Project/Embla/2007/49010001 (1)/edf_file.edf" False
```

## References

[Dive into the world of EDF+](https://www.edfplus.info/)
