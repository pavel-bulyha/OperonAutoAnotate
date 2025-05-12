# Operon Annotation Tool

A Python utility to annotate operon regions in GenBank files using OperonDB data. It uses Biopython for parsing and Playwright + BeautifulSoup for web scraping.

## Files

- **annotate_operons.py**: Reads a GenBank file, maps genes, computes operon boundaries, and adds annotations.
- **operondbSearch.py**: Scrapes OperonDB and returns a DataFrame with Operon ID, Species, Name, Genes, and Definition.

## Installation & Run

Requires Python 3.7+ and the packages: Biopython, Pandas, BeautifulSoup4, and Playwright.

```bash
pip install biopython pandas beautifulsoup4 playwright
playwright install
```
## Usage

Run the script from the command line. You must specify the path to the input GenBank file. Optionally, provide the output file name and modify the behavior with the --include_all or --no_include_all flag.

```bash
python annotate_operons.py Haemophilus_influenzae_Rd_KW20.gb -o annotated.gb --include_all
```