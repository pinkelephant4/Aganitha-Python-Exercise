# 🧬 PubMed Paper Fetcher

[![TestPyPI](https://img.shields.io/badge/TestPyPI-pubmed--paper--fetcher--ananya-purple)](https://test.pypi.org/project/pubmed-paper-fetcher-ananya)

A Python command-line tool to fetch research papers from **PubMed** based on a user-defined query, filtering for papers that include **non-academic authors affiliated with pharmaceutical or biotech companies**.


## 🚀 Features

- Supports full PubMed query syntax
- Identifies non-academic (industry) authors using heuristics
- Filters and returns company-affiliated papers
- Outputs to CSV or nicely formatted console table
- Supports up to 10,000 articles per query


## 📁 Project Structure
```
paper_fetcher/
│
├── src/
│ └── paper_fetcher/
│ ├── core.py # Core logic to search, fetch, and extract data from PubMed
│ └── cli.py # Command-line interface using argparse
│
├── pyproject.toml # Poetry dependency and packaging config
├── README.md # This file
```

## ⚙️ Installation


### ▶️ From TestPyPI:

```bash
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple \
            pubmed-paper-fetcher-ananya
```     

> 📦 TestPyPI Page: [https://test.pypi.org/project/pubmed-paper-fetcher-ananya](https://test.pypi.org/project/pubmed-paper-fetcher-ananya)

#### 🧪 Example Usage

```bash
get-papers-list "mRNA vaccine" --limit 50 --debug
```

### ▶️ From Github Repo (Development)

#### ✅ Step 1: Clone the Repository

```bash
git clone https://github.com/pinkelephant4/paper-fetcher.git
cd paper-fetcher
```

#### ✅ Step 2: Install with [Poetry](https://python-poetry.org/)

If you don't have Poetry installed:

```bash
pipx install poetry
```

Then run:

```bash
poetry install
```
Run locally with:

```bash
poetry run get-papers-list "cancer therapy" --limit 100
```

## 🚀 Usage

Run the program using Poetry's CLI wrapper:

```bash
poetry run get-papers-list "your pubmed query here"
```

### 🔧 Command-Line Options

| Option          | Description                                                |
| --------------- | ---------------------------------------------------------- |
| `query`         | **(Required)** Your PubMed query string (quoted)           |
| `-f`, `--file`  | Output CSV file (if omitted, prints to console)            |
| `-l`, `--limit` | Max number of articles to fetch (default: 1000, max: 9999) |
| `-d`, `--debug` | Print debug logs to stderr                                 |
| `-h`, `--help`  | Show help message                                          |


### 📌 Example

```bash
poetry run get-papers-list "covid 19" -f output.csv --limit 3000 --debug
```


## 🧰 Tools & Libraries Used

* **[Biopython](https://biopython.org/)** (`Bio.Entrez`) – For accessing the PubMed API
* **[Poetry](https://python-poetry.org/)** – For dependency management and packaging
* **[argparse](https://docs.python.org/3/library/argparse.html)** – For CLI
* **GitHub Copilot & ChatGPT** – Used for auto-completion, suggestions, and documentation generation


## 📜 License

MIT License — free to use, modify, and distribute.


