# 🧬 PubMed Paper Fetcher

A Python command-line tool to fetch research papers from **PubMed** based on a user-defined query, filtering for papers that include **non-academic authors affiliated with pharmaceutical or biotech companies**.


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


### ✅ Step 1: Clone the Repository

```bash
git clone https://github.com/yourusername/paper-fetcher.git
cd paper-fetcher
```

### ✅ Step 2: Install with [Poetry](https://python-poetry.org/)

If you don't have Poetry installed:

```bash
pipx install poetry
```

Then run:

```bash
poetry install
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


